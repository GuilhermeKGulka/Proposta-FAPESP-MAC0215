! Hessiana da face na GPU (fase 6). O estado mora em lssit%gpu (lssit.f90, sob
! #ifdef GENCAN_GPU), assim como a Hessiana da CPU mora em lssit%hrow/hcol/hval.
!
! Fluxo (chamado no fim de preparehv, so com truehp e compressed):
!   1. preparehv original deixa o COO da face compactado e renumerado (1..nfr),
!      um triangulo por entrada (como em computetruehv);
!   2. espelha (c,r) de cada entrada fora da diagonal no fim dos proprios vetores
!      do lssit (capacidade 2*hnnzmax): a API generica do cuSPARSE nao faz SpMV
!      "simetrica, so um triangulo", entao a matriz vai cheia para os dois backends;
!   3. coo2csr in-place (rotina do orientador, docs/coo2csr_orientador.f90);
!   4. envia o CSR. Se o COO compactado e igual ao do ultimo padrao, row_ptr/col,
!      o descritor, o buffer e o preprocess continuam valendo e so os valores sobem.
!
! Depois disso lssit%hrow/hcol/hval NAO sao mais um COO valido: o produto H v da
! face so pode ser feito pela GPU (lssit_gpu_spmv) ate o proximo evalh.
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

  type :: lssit_gpu_type
    logical :: initialized = .false.
    logical :: active = .false.        ! a Hessiana da face atual esta na GPU
    logical :: changed = .false.       ! o ultimo prepare refez o padrao
    integer :: n = 0, cap = 0          ! cap: capacidade de lssit%hrow/hcol/hval
    integer :: nfr = -1, hne = -1      ! padrao atual (hne: triangulo compactado)
    integer :: nnz = 0                 ! entradas da matriz cheia (com o espelho)
    logical :: empty = .true., spmv_ready = .false.
    real(8) :: t_prep_ms = 0.0d0       ! tempo do SpMV_preprocess no ultimo padrao novo

    ! host
    integer, allocatable :: prow(:), pcol(:)   ! triangulo compactado do padrao atual
    integer, allocatable :: alen(:)            ! alen do coo2csr
    integer, allocatable :: row_ptr(:)         ! asta do coo2csr + fim (integer = int32)

    ! device. d_x/d_y: vetores do produto y = H x, capacidade n
    type(c_ptr) :: handle = c_null_ptr
    type(c_ptr) :: d_row_ptr = c_null_ptr, d_col = c_null_ptr, d_val = c_null_ptr
    type(c_ptr) :: d_x = c_null_ptr, d_y = c_null_ptr
    type(c_ptr) :: mat = c_null_ptr, vx = c_null_ptr, vy = c_null_ptr, d_buf = c_null_ptr
  end type lssit_gpu_type

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine lssit_gpu_init(gpu, n, cap)
    type(lssit_gpu_type), intent(inout) :: gpu
    integer, intent(in) :: n, cap

    gpu%n = n;  gpu%cap = cap
    allocate(gpu%prow(cap), gpu%pcol(cap), gpu%alen(n), gpu%row_ptr(n+1))
    call hipCheck(hipMalloc(gpu%d_row_ptr, (n+1) * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_col, cap * 4_c_size_t))
    call hipCheck(hipMalloc(gpu%d_val, cap * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_x, n * 8_c_size_t))
    call hipCheck(hipMalloc(gpu%d_y, n * 8_c_size_t))
    call hipsparseCheck(hipsparseCreate(gpu%handle))
    gpu%nfr = -1;  gpu%hne = -1
    gpu%initialized = .true.
  end subroutine lssit_gpu_init

  ! *****************************************************************
  ! *****************************************************************

  ! hrow/hcol/hval(1:hne): COO da face ja compactado por preparehv (indices 1..nfr).
  ! Na saida, hrow/hcol/hval foram sobrescritos (espelho + coo2csr).
  subroutine lssit_gpu_prepare(gpu, n, cap, nfr, hne, hrow, hcol, hval)
    type(lssit_gpu_type), intent(inout), target :: gpu
    integer, intent(in) :: n, cap, nfr, hne
    integer, intent(inout), target :: hrow(cap), hcol(cap)
    real(kind=8), intent(inout), target :: hval(cap)

    integer :: e, nnz
    logical :: same

    if (.not. gpu%initialized) call lssit_gpu_init(gpu, n, cap)

    ! mesmo padrao? (compara o COO compactado; vale para qualquer evalh, nao so CUTEst)
    same = (nfr == gpu%nfr .and. hne == gpu%hne)
    if (same) same = all(hrow(1:hne) == gpu%prow(1:hne)) .and. all(hcol(1:hne) == gpu%pcol(1:hne))
    gpu%changed = .not. same
    if (.not. same) then
      call teardown_spmv(gpu)
      gpu%nfr = nfr;  gpu%hne = hne
      gpu%prow(1:hne) = hrow(1:hne);  gpu%pcol(1:hne) = hcol(1:hne)
    end if

    ! espelho: (c,r) de cada entrada fora da diagonal, no fim dos vetores
    nnz = hne
    do e = 1, hne
      if (hrow(e) /= hcol(e)) then
        nnz = nnz + 1
        if (nnz > cap) stop 'lssit_gpu_prepare: capacidade insuficiente (esperado 2*hnnzmax)'
        hrow(nnz) = hcol(e);  hcol(nnz) = hrow(e);  hval(nnz) = hval(e)
      end if
    end do
    gpu%nnz = nnz
    gpu%empty = (nfr == 0 .or. nnz == 0)
    gpu%active = .true.
    if (gpu%empty) return

    ! COO -> CSR in-place; row_ptr(1:nfr) recebe o inicio de cada linha
    call coo2csr(nfr, nnz, hrow, hcol, hval, gpu%alen, gpu%row_ptr)

    if (.not. same) then
      gpu%row_ptr(nfr+1) = nnz + 1
      call hipCheck(hipMemcpy(gpu%d_row_ptr, c_loc(gpu%row_ptr), (nfr+1) * 4_c_size_t, hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(gpu%d_col, c_loc(hcol), nnz * 4_c_size_t, hipMemcpyHostToDevice))
    end if
    call hipCheck(hipMemcpy(gpu%d_val, c_loc(hval), nnz * 8_c_size_t, hipMemcpyHostToDevice))

    if (.not. same) call setup_spmv(gpu)
  end subroutine lssit_gpu_prepare

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
    call hipCheck(hipFree(gpu%d_x))
    call hipCheck(hipFree(gpu%d_y))
    call hipsparseCheck(hipsparseDestroy(gpu%handle))
    deallocate(gpu%prow, gpu%pcol, gpu%alen, gpu%row_ptr)
    gpu%initialized = .false.;  gpu%active = .false.
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

  ! Rotina do orientador (docs/coo2csr_orientador.f90); codigo executavel sem
  ! alteracoes (so os comentarios internos foram resumidos aqui).
  ! In-place, O(nnz). Na saida asta(i) = inicio da linha i (base 1), alen(i) =
  ! tamanho da linha i; a ordem das colunas dentro de cada linha nao e garantida.
  ! Adaptacao: virou procedimento do modulo; quem chama ja filtrou a face,
  ! renumerou (preparehv) e espelhou o triangulo (lssit_gpu_prepare).
  subroutine coo2csr(m,nnz,arow,acol,aval,alen,asta)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,nnz

    ! ARRAY ARGUMENTS
    integer,      intent(inout) :: acol(nnz),alen(m),arow(nnz),asta(m)
    real(kind=8), intent(inout) :: aval(nnz)

    ! LOCAL SCALARS
    integer      :: i,j,col,coltmp,row,rowtmp
    real(kind=8) :: val,valtmp

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

       val = aval(i)
       col = acol(i)
       row = arow(i)

       arow(i) = - 1

       do while ( row .ge. 0 )
          j = asta(row)
          asta(row) = j + 1

          valtmp = aval(j)
          coltmp = acol(j)
          rowtmp = arow(j)

          aval(j) = val
          acol(j) = col
          arow(j) = - 1

          val = valtmp
          col = coltmp
          row = rowtmp
       end do

    end do

    asta(1:m) = asta(1:m) - alen(1:m)

  end subroutine coo2csr

end module lssit_gpu_mod
