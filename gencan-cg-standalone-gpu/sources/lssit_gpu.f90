! Hessiana da face na GPU (fase 6). O estado mora em lssit%gpu (lssit.f90, sob
! #ifdef GENCAN_GPU), assim como a Hessiana da CPU mora em lssit%hrow/hcol/hval.
!
! Chamado por preparehv (so com truehp e compressed), NO LUGAR da compactacao do
! COO feita na CPU. Parte do COO do evalh (indices 1..n, um triangulo) e de frvar:
!   - padrao novo (estrutura do COO ou face mudou): renumera e filtra pela face
!     (mesma regra de preparehv), espelha as entradas fora da diagonal (a API
!     generica do cuSPARSE nao faz SpMV "simetrica, so um triangulo"; a ordem
!     dentro da linha do CSR nao importa, entao o espelho vai junto da entrada),
!     converte para CSR com o coo2csr do orientador carregando a origem de cada
!     entrada, e guarda pos1/pos2 = posicao no CSR de cada entrada do COO original;
!     envia row_ptr/col/pos1/pos2 e refaz descritor, buffer e preprocess;
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
    real(8) :: t_conv_ms = 0.0d0       ! parte de t_build_ms gasta no coo2csr

    ! host: estrutura do ultimo padrao (para detectar mudanca)
    integer, allocatable :: srow(:), scol(:)
    logical, allocatable :: sfrvar(:)
    ! host: montagem do padrao
    integer, allocatable :: renum(:), alen(:)
    integer, allocatable :: wrow(:), wcol(:), wid(:)   ! 2*hnnzmax: COO da face com espelho
    integer, allocatable :: row_ptr(:)                 ! n+1 (integer = int32)
    integer, allocatable :: pos1(:), pos2(:)           ! hnnzmax

    ! device. d_x/d_y: vetores do produto y = H x, capacidade n
    type(c_ptr) :: handle = c_null_ptr
    type(c_ptr) :: d_row_ptr = c_null_ptr, d_col = c_null_ptr, d_val = c_null_ptr
    type(c_ptr) :: d_hcoo = c_null_ptr, d_pos1 = c_null_ptr, d_pos2 = c_null_ptr
    type(c_ptr) :: d_x = c_null_ptr, d_y = c_null_ptr
    type(c_ptr) :: mat = c_null_ptr, vx = c_null_ptr, vy = c_null_ptr, d_buf = c_null_ptr
  end type lssit_gpu_type

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine lssit_gpu_init(gpu, n, hnnzmax)
    type(lssit_gpu_type), intent(inout) :: gpu
    integer, intent(in) :: n, hnnzmax

    gpu%n = n;  gpu%hnnzmax = hnnzmax
    allocate(gpu%srow(hnnzmax), gpu%scol(hnnzmax), gpu%sfrvar(n), gpu%renum(n), gpu%alen(n), &
             gpu%wrow(2*hnnzmax), gpu%wcol(2*hnnzmax), gpu%wid(2*hnnzmax), gpu%row_ptr(n+1), &
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
    gpu%has_pattern = .false.
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

    if (.not. gpu%initialized) call lssit_gpu_init(gpu, n, hnnzmax)

    ! mesmo padrao? (estrutura do COO e face; vale para qualquer evalh, nao so CUTEst)
    same = gpu%has_pattern .and. hne == gpu%hne
    if (same) same = all(frvar(1:n) .eqv. gpu%sfrvar(1:n))
    if (same) same = all(hrow(1:hne) == gpu%srow(1:hne)) .and. all(hcol(1:hne) == gpu%scol(1:hne))
    gpu%changed = .not. same

    if (.not. same) then
      call teardown_spmv(gpu)
      call build_pattern(gpu, n, frvar, hne, hrow, hcol)
      if (.not. gpu%empty) then
        call hipCheck(hipMemcpy(gpu%d_row_ptr, c_loc(gpu%row_ptr), (gpu%nfr+1) * 4_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(gpu%d_col, c_loc(gpu%wcol), gpu%nnz * 4_c_size_t, hipMemcpyHostToDevice))
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
  end subroutine lssit_gpu_prepare

  ! *****************************************************************
  ! *****************************************************************

  ! Monta o CSR nfr x nfr da face: row_ptr, col (em wcol) e pos1/pos2.
  subroutine build_pattern(gpu, n, frvar, hne, hrow, hcol)
    type(lssit_gpu_type), intent(inout) :: gpu
    integer, intent(in) :: n, hne
    logical, intent(in) :: frvar(n)
    integer, intent(in) :: hrow(hne), hcol(hne)
    integer :: i, j, e, r, c, k, nnz
    integer(8) :: t0, t1, t2, t3, rate

    call system_clock(t0, rate)
    t1 = t0;  t2 = t0
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
    gpu%nfr = k

    ! COO da face com espelho; wid guarda a origem (+e direta, -e espelhada)
    nnz = 0
    gpu%pos1(1:hne) = 0;  gpu%pos2(1:hne) = 0
    do e = 1, hne
      r = gpu%renum(hrow(e));  c = gpu%renum(hcol(e))
      if (r == 0 .or. c == 0) cycle
      nnz = nnz + 1
      gpu%wrow(nnz) = r;  gpu%wcol(nnz) = c;  gpu%wid(nnz) = e
      if (r /= c) then
        nnz = nnz + 1
        gpu%wrow(nnz) = c;  gpu%wcol(nnz) = r;  gpu%wid(nnz) = -e
      end if
    end do
    gpu%nnz = nnz
    gpu%empty = (gpu%nfr == 0 .or. nnz == 0)

    if (.not. gpu%empty) then
      call system_clock(t1)
      call coo2csr(gpu%nfr, nnz, gpu%wrow, gpu%wcol, gpu%wid, gpu%alen, gpu%row_ptr)
      call system_clock(t2)
      gpu%row_ptr(gpu%nfr+1) = nnz + 1
      do j = 1, nnz
        e = gpu%wid(j)
        if (e > 0) then
          gpu%pos1(e) = j
        else
          gpu%pos2(-e) = j
        end if
      end do
    end if

    gpu%hne = hne
    gpu%srow(1:hne) = hrow(1:hne);  gpu%scol(1:hne) = hcol(1:hne)
    gpu%sfrvar(1:n) = frvar(1:n)
    gpu%has_pattern = .true.
    call system_clock(t3)
    gpu%t_conv_ms  = real(t2 - t1, 8) / real(rate, 8) * 1d3
    gpu%t_build_ms = real(t3 - t0, 8) / real(rate, 8) * 1d3
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
    deallocate(gpu%srow, gpu%scol, gpu%sfrvar, gpu%renum, gpu%alen, gpu%wrow, gpu%wcol, &
               gpu%wid, gpu%row_ptr, gpu%pos1, gpu%pos2)
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

  ! *****************************************************************
  ! *****************************************************************

  ! Rotina do orientador (docs/coo2csr_orientador.f90), in-place, O(nnz).
  ! Adaptacao: o vetor de valores (aval, real) virou aid, inteiro, que carrega
  ! a origem de cada entrada; o resto do codigo executavel e o original.
  ! Na saida asta(i) = inicio da linha i (base 1), alen(i) = tamanho da linha i;
  ! a ordem das colunas dentro de cada linha nao e garantida (nao precisa).
  subroutine coo2csr(m,nnz,arow,acol,aid,alen,asta)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,nnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: acol(nnz),alen(m),arow(nnz),asta(m)
    integer, intent(inout) :: aid(nnz)

    ! LOCAL SCALARS
    integer :: i,j,col,coltmp,row,rowtmp
    integer :: id,idtmp

    alen(1:m) = 0

    do i = 1,nnz
       row = arow(i)
       alen(row) = alen(row) + 1
    end do

    asta(1) = 1
    do i = 2,m
       asta(i) = asta(i-1) + alen(i-1)
    end do

    do i = 1,nnz

       id  = aid(i)
       col = acol(i)
       row = arow(i)

       arow(i) = - 1

       do while ( row .ge. 0 )
          j = asta(row)
          asta(row) = j + 1

          idtmp  = aid(j)
          coltmp = acol(j)
          rowtmp = arow(j)

          aid(j)  = id
          acol(j) = col
          arow(j) = - 1

          id  = idtmp
          col = coltmp
          row = rowtmp
       end do

    end do

    asta(1:m) = asta(1:m) - alen(1:m)

  end subroutine coo2csr

end module lssit_gpu_mod
