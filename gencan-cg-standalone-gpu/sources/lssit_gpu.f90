! Hessiana da face na GPU (fase 6). O estado mora em lssit%gpu (lssit.f90, sob
! #ifdef GENCAN_GPU), assim como a Hessiana da CPU mora em lssit%hrow/hcol/hval.
!
! Chamado por preparehv (so com truehp e compressed), NO LUGAR da compactacao do
! COO feita na CPU. Parte do COO do evalh (indices 1..n, um triangulo) e de frvar:
!   - padrao novo (estrutura do COO ou face mudou): renumera e filtra pela face
!     (mesma regra de preparehv) e monta o CSR com o espelho das entradas fora da
!     diagonal (a API generica do cuSPARSE nao faz SpMV "simetrica, so um
!     triangulo"). A conversao segue o coo2csr do orientador (conta por linha,
!     acumula o inicio de cada linha), mas coloca as entradas fora do lugar, com
!     um cursor por linha, em vez da permutacao in-place por ciclos (3x mais
!     lenta aqui; versao com o coo2csr no branch fase6-coo2csr). A ordem dentro
!     da linha nao importa. Guarda pos1/pos2 = posicao no CSR de cada entrada do
!     COO original; envia row_ptr/col/pos1/pos2 e refaz descritor, buffer e preprocess;
!   - todo passo: envia os hne valores do COO e espalha no CSR por kernel.
! Toda a memoria (host e device) e alocada uma vez, na primeira chamada, com a
! capacidade maxima (n, hnnzmax); o COO do lssit nao e alterado.
module lssit_gpu_mod
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipsparse_generic
  implicit none
  private
  public :: lssit_gpu_type, lssit_gpu_prepare, lssit_gpu_spmv, lssit_gpu_free

  real(c_double), target, save :: one = 1.0d0, zero = 0.0d0

  interface
    ! cgm_kernels.cpp: val(p1(e)) = v(e) e, se p2(e) > 0, val(p2(e)) = v(e)
    subroutine launch_scatter(hne, v, p1, p2, val) bind(C, name="launch_scatter")
      import :: c_int, c_ptr
      integer(c_int), value :: hne
      type(c_ptr), value :: v, p1, p2, val
    end subroutine
  end interface

  type :: lssit_gpu_type
    logical :: initialized = .false.
    logical :: active = .false.        ! a Hessiana da face atual esta na GPU
    logical :: changed = .false.       ! o ultimo prepare refez o padrao
    logical :: has_pattern = .false.
    integer :: n = 0, hnnzmax = 0
    integer :: nfr = 0, hne = 0        ! padrao atual (hne: entradas do COO original)
    integer :: nnz = 0                 ! entradas da matriz cheia (com o espelho)
    logical :: empty = .true., spmv_ready = .false.
    real(8) :: t_prep_ms = 0.0d0       ! tempo do SpMV_preprocess no ultimo padrao novo
    real(8) :: t_build_ms = 0.0d0      ! tempo de build_pattern (host) no ultimo padrao novo

    ! host: estrutura do ultimo padrao (para detectar mudanca)
    integer, allocatable :: srow(:), scol(:)
    logical, allocatable :: sfrvar(:)
    ! host: montagem do padrao
    integer, allocatable :: renum(:), alen(:)          ! alen: contagem, depois cursor por linha
    integer, allocatable :: col(:)                     ! 2*hnnzmax (integer = int32)
    integer, allocatable :: row_ptr(:)                 ! n+1 (integer = int32)
    integer, allocatable :: pos1(:), pos2(:)           ! hnnzmax

    ! device. d_x/d_y: vetores do produto y = H x, capacidade n
    type(c_ptr) :: handle = c_null_ptr
    type(c_ptr) :: d_row_ptr = c_null_ptr, d_col = c_null_ptr, d_val = c_null_ptr
    type(c_ptr) :: d_hcoo = c_null_ptr, d_pos1 = c_null_ptr, d_pos2 = c_null_ptr
    type(c_ptr) :: d_x = c_null_ptr, d_y = c_null_ptr
    type(c_ptr) :: mat = c_null_ptr, vx = c_null_ptr, vy = c_null_ptr, d_buf = c_null_ptr

    ! device: vetores de trabalho do CG (cgm_gpu.f90), capacidade n; p = d_x, hp = d_y.
    ! Alocados aqui, uma vez, porque o cgm recebe o lssit como intent(in).
    type(c_ptr) :: d_g = c_null_ptr, d_d = c_null_ptr, d_r = c_null_ptr
    type(c_ptr) :: d_xf = c_null_ptr, d_lbf = c_null_ptr, d_ubf = c_null_ptr
    type(c_ptr) :: d_lindf = c_null_ptr, d_uindf = c_null_ptr
    type(c_ptr) :: d_dg = c_null_ptr, d_s = c_null_ptr, d_u = c_null_ptr
    type(c_ptr) :: d_part = c_null_ptr, d_cnt = c_null_ptr, d_red = c_null_ptr

    ! contadores (medicao): chamadas de prepare, padroes novos, tempo total de prepare
    integer :: nprep = 0, nnew = 0
    real(8) :: t_prep_tot_ms = 0.0d0
  end type lssit_gpu_type

  ! blocos maximos das reducoes do CG (tem que bater com cgm_kernels.cpp)
  integer, parameter, public :: cg_maxb = 1024

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine lssit_gpu_init(gpu, n, hnnzmax)
    type(lssit_gpu_type), intent(inout) :: gpu
    integer, intent(in) :: n, hnnzmax

    gpu%n = n;  gpu%hnnzmax = hnnzmax
    allocate(gpu%srow(hnnzmax), gpu%scol(hnnzmax), gpu%sfrvar(n), gpu%renum(n), gpu%alen(n), &
             gpu%col(2*hnnzmax), gpu%row_ptr(n+1), &
             gpu%pos1(hnnzmax), gpu%pos2(hnnzmax))
    call hipCheck(hipMalloc(gpu%d_row_ptr, (n+1) * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_col, 2 * hnnzmax * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_val, 2 * hnnzmax * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_hcoo, hnnzmax * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_pos1, hnnzmax * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_pos2, hnnzmax * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_x, n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_y, n * 8_c_size_t))
    call hipsparseCheck(hipsparseCreate(gpu%handle))
    ! CG
    call hipCheck(hipMalloc(gpu%d_g,  n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_d,  n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_r,  n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_xf,  n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_lbf, n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_ubf, n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_lindf, n * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_uindf, n * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_dg, n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_s,  n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_u,  n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_part, cg_maxb * 7 * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_cnt, 4_c_size_t))
    call hipCheck(hipMemset(gpu%d_cnt, 0, 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_red, 16 * 8_c_size_t))
    call hipCheck(hipMemset(gpu%d_red, 0, 16 * 8_c_size_t))
    gpu%has_pattern = .false.
    gpu%nprep = 0;  gpu%nnew = 0;  gpu%t_prep_tot_ms = 0.0d0
    gpu%initialized = .true.
  end subroutine lssit_gpu_init

  ! *****************************************************************
  ! *****************************************************************

  ! hrow/hcol/hval(1:hne): COO do evalh (indices 1..n); frvar: face livre.
  subroutine lssit_gpu_prepare(gpu, n, hnnzmax, frvar, hne, hrow, hcol, hval)
    type(lssit_gpu_type), intent(inout), target :: gpu
    integer, intent(in) :: n, hnnzmax, hne
    logical, intent(in) :: frvar(n)
    integer, intent(in) :: hrow(hnnzmax), hcol(hnnzmax)
    real(kind=8), intent(in), target :: hval(hnnzmax)

    logical :: same
    integer :: e
    integer(8) :: t0, t1, rate

    if (.not. gpu%initialized) call lssit_gpu_init(gpu, n, hnnzmax)
    call system_clock(t0, rate)

    ! mesmo padrao? (estrutura do COO e face; vale para qualquer evalh, nao so CUTEst)
    same = gpu%has_pattern .and. hne == gpu%hne
    if (same) same = all(frvar(1:n) .eqv. gpu%sfrvar(1:n))
    if (same) then                     ! um so laco, para na primeira diferenca
      do e = 1, hne
        if (hrow(e) /= gpu%srow(e) .or. hcol(e) /= gpu%scol(e)) then
          same = .false.
          exit
        end if
      end do
    end if
    gpu%changed = .not. same

    if (.not. same) then
      call teardown_spmv(gpu)
      call build_pattern(gpu, n, frvar, hne, hrow, hcol)
      if (.not. gpu%empty) then
        call hipCheck(hipMemcpy(gpu%d_row_ptr, c_loc(gpu%row_ptr), (gpu%nfr+1) * 4_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(gpu%d_col, c_loc(gpu%col), gpu%nnz * 4_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(gpu%d_pos1, c_loc(gpu%pos1), hne * 4_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(gpu%d_pos2, c_loc(gpu%pos2), hne * 4_c_size_t, hipMemcpyHostToDevice))
      end if
    end if

    ! valores: COO -> GPU, espalhados no CSR por kernel
    if (.not. gpu%empty) then
      call hipCheck(hipMemcpy(gpu%d_hcoo, c_loc(hval), hne * 8_c_size_t, hipMemcpyHostToDevice))
      call launch_scatter(int(hne, c_int), gpu%d_hcoo, gpu%d_pos1, gpu%d_pos2, gpu%d_val)
    end if

    if (.not. same .and. .not. gpu%empty) call setup_spmv(gpu)
    gpu%active = .true.

    ! contadores (o scatter e assincrono: o tempo nao espera o kernel terminar)
    call system_clock(t1)
    gpu%nprep = gpu%nprep + 1
    if (.not. same) gpu%nnew = gpu%nnew + 1
    gpu%t_prep_tot_ms = gpu%t_prep_tot_ms + real(t1 - t0, 8) / real(rate, 8) * 1d3
  end subroutine lssit_gpu_prepare

  ! *****************************************************************
  ! *****************************************************************

  ! Monta o CSR nfr x nfr da face: row_ptr, col e pos1/pos2.
  subroutine build_pattern(gpu, n, frvar, hne, hrow, hcol)
    type(lssit_gpu_type), intent(inout) :: gpu
    integer, intent(in) :: n, hne
    logical, intent(in) :: frvar(n)
    integer, intent(in) :: hrow(hne), hcol(hne)
    integer :: i, j, e, r, c, k, nfr
    integer(8) :: t0, t1, rate

    call system_clock(t0, rate)
    ! renumeracao identica a de preparehv: i-esimo .true. de frvar -> i
    k = 0
    do i = 1, n
      if (frvar(i)) then
        k = k + 1
        gpu%renum(i) = k
      else
        gpu%renum(i) = 0
      end if
    end do
    nfr = k
    gpu%nfr = nfr

    ! como no coo2csr: conta as entradas de cada linha (com o espelho) ...
    gpu%alen(1:nfr) = 0
    do e = 1, hne
      r = gpu%renum(hrow(e));  c = gpu%renum(hcol(e))
      if (r == 0 .or. c == 0) cycle
      gpu%alen(r) = gpu%alen(r) + 1
      if (r /= c) gpu%alen(c) = gpu%alen(c) + 1
    end do

    ! ... e acumula o inicio de cada linha
    gpu%row_ptr(1) = 1
    do i = 1, nfr
      gpu%row_ptr(i+1) = gpu%row_ptr(i) + gpu%alen(i)
    end do
    gpu%nnz = gpu%row_ptr(nfr+1) - 1
    gpu%empty = (nfr == 0 .or. gpu%nnz == 0)

    ! colocacao fora do lugar: alen vira o cursor de cada linha; pos1/pos2 guardam
    ! onde cada entrada do COO original (e seu espelho) caiu no CSR
    gpu%alen(1:nfr) = gpu%row_ptr(1:nfr)
    do e = 1, hne
      gpu%pos1(e) = 0;  gpu%pos2(e) = 0
      r = gpu%renum(hrow(e));  c = gpu%renum(hcol(e))
      if (r == 0 .or. c == 0) cycle
      j = gpu%alen(r);  gpu%col(j) = c;  gpu%pos1(e) = j;  gpu%alen(r) = j + 1
      if (r /= c) then
        j = gpu%alen(c);  gpu%col(j) = r;  gpu%pos2(e) = j;  gpu%alen(c) = j + 1
      end if
    end do

    gpu%hne = hne
    gpu%srow(1:hne) = hrow(1:hne);  gpu%scol(1:hne) = hcol(1:hne)
    gpu%sfrvar(1:n) = frvar(1:n)
    gpu%has_pattern = .true.
    call system_clock(t1)
    gpu%t_build_ms = real(t1 - t0, 8) / real(rate, 8) * 1d3
  end subroutine build_pattern

  ! *****************************************************************
  ! *****************************************************************

  ! y = H_ff x, com x em gpu%d_x e y em gpu%d_y (tamanho nfr)
  subroutine lssit_gpu_spmv(gpu)
    type(lssit_gpu_type), intent(in) :: gpu
    if (gpu%nfr <= 0) return
    if (gpu%empty) then
      call hipCheck(hipMemset(gpu%d_y, 0, gpu%nfr * 8_c_size_t))
    else
      call hs_check(hs_spmv(gpu%handle, HS_OP_N, c_loc(one), gpu%mat, gpu%vx, c_loc(zero), &
                    gpu%vy, HS_R_64F, HS_SPMV_ALG_DEFAULT, gpu%d_buf), "SpMV")
    end if
  end subroutine lssit_gpu_spmv

  ! *****************************************************************
  ! *****************************************************************

  subroutine lssit_gpu_free(gpu)
    type(lssit_gpu_type), intent(inout) :: gpu
    if (.not. gpu%initialized) return
    write(*,'(A,I0,A,I0,A,F12.3)') ' (GPU) preparehv: chamadas = ', gpu%nprep, &
         ', padroes novos = ', gpu%nnew, ', tempo total (ms) = ', gpu%t_prep_tot_ms
    call teardown_spmv(gpu)
    call hipCheck(hipFree(gpu%d_row_ptr))
    call hipCheck(hipFree(gpu%d_col))
    call hipCheck(hipFree(gpu%d_val))
    call hipCheck(hipFree(gpu%d_hcoo))
    call hipCheck(hipFree(gpu%d_pos1))
    call hipCheck(hipFree(gpu%d_pos2))
    call hipCheck(hipFree(gpu%d_x))
    call hipCheck(hipFree(gpu%d_y))
    call hipsparseCheck(hipsparseDestroy(gpu%handle))
    call hipCheck(hipFree(gpu%d_g));  call hipCheck(hipFree(gpu%d_d));  call hipCheck(hipFree(gpu%d_r))
    call hipCheck(hipFree(gpu%d_xf)); call hipCheck(hipFree(gpu%d_lbf)); call hipCheck(hipFree(gpu%d_ubf))
    call hipCheck(hipFree(gpu%d_lindf)); call hipCheck(hipFree(gpu%d_uindf))
    call hipCheck(hipFree(gpu%d_dg)); call hipCheck(hipFree(gpu%d_s)); call hipCheck(hipFree(gpu%d_u))
    call hipCheck(hipFree(gpu%d_part)); call hipCheck(hipFree(gpu%d_cnt)); call hipCheck(hipFree(gpu%d_red))
    deallocate(gpu%srow, gpu%scol, gpu%sfrvar, gpu%renum, gpu%alen, gpu%col, &
               gpu%row_ptr, gpu%pos1, gpu%pos2)
    gpu%initialized = .false.;  gpu%active = .false.;  gpu%has_pattern = .false.
  end subroutine lssit_gpu_free

  ! *****************************************************************
  ! *****************************************************************

  ! descritores, buffer e preprocess do SpMV para o padrao atual
  subroutine setup_spmv(gpu)
    type(lssit_gpu_type), intent(inout) :: gpu
    integer(c_size_t) :: bsz
    integer(c_int) :: st
    integer(8) :: t0, t1, rate

    gpu%t_prep_ms = 0.0d0
    call hs_check(hs_create_csr(gpu%mat, int(gpu%nfr, c_int64_t), int(gpu%nfr, c_int64_t), &
                  int(gpu%nnz, c_int64_t), gpu%d_row_ptr, gpu%d_col, gpu%d_val, &
                  HS_INDEX_32I, HS_INDEX_32I, HS_BASE_ONE, HS_R_64F), "CreateCsr")
    call hs_check(hs_create_dnvec(gpu%vx, int(gpu%nfr, c_int64_t), gpu%d_x, HS_R_64F), "CreateDnVec x")
    call hs_check(hs_create_dnvec(gpu%vy, int(gpu%nfr, c_int64_t), gpu%d_y, HS_R_64F), "CreateDnVec y")
    call hs_check(hs_spmv_buffer_size(gpu%handle, HS_OP_N, c_loc(one), gpu%mat, gpu%vx, &
                  c_loc(zero), gpu%vy, HS_R_64F, HS_SPMV_ALG_DEFAULT, bsz), "SpMV_bufferSize")
    call hs_alloc_buffer(gpu%d_buf, bsz)          ! nulo se bsz == 0

    call system_clock(count_rate=rate)
    call system_clock(t0)
    st = hs_spmv_preprocess(gpu%handle, HS_OP_N, c_loc(one), gpu%mat, gpu%vx, &
                c_loc(zero), gpu%vy, HS_R_64F, HS_SPMV_ALG_DEFAULT, gpu%d_buf)
    call hipCheck(hipDeviceSynchronize())
    call system_clock(t1)
    gpu%t_prep_ms = real(t1 - t0, 8) / real(rate, 8) * 1000.0d0
    if (st /= 0) then
      write(*,'(A,I0,A,I0)') "AVISO: SpMV_preprocess falhou (codigo ", st, "; bufferSize=", bsz
      gpu%t_prep_ms = -1.0d0
    end if
    gpu%spmv_ready = .true.
  end subroutine setup_spmv

  ! *****************************************************************
  ! *****************************************************************

  subroutine teardown_spmv(gpu)
    type(lssit_gpu_type), intent(inout) :: gpu
    if (.not. gpu%spmv_ready) return
    call hs_check(hs_destroy_spmat(gpu%mat), "DestroySpMat")
    call hs_check(hs_destroy_dnvec(gpu%vx), "DestroyDnVec x")
    call hs_check(hs_destroy_dnvec(gpu%vy), "DestroyDnVec y")
    call hs_free_buffer(gpu%d_buf)
    gpu%spmv_ready = .false.
  end subroutine teardown_spmv

end module lssit_gpu_mod
