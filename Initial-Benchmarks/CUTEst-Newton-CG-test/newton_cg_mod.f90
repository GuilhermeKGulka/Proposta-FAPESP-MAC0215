module newton_cg_mod
  use iso_c_binding
  use rocm_context
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipfort_hipBLAS
  use benchmark_cg_mod
  use cutest_iface
  implicit none

  private
  public :: newton_stats, newton_cg_gpu, newton_cg_cpu

  type :: newton_stats
    integer :: newton_iter      = 0
    integer :: fcnt             = 0
    integer :: gcnt             = 0
    integer :: hcnt             = 0
    integer :: istop            = -1
    integer :: cg_fails         = 0
    integer :: first_cg_fail    = -1
    integer :: last_cg_iter     = 0
    real(8) :: cg_time_ms       = 0.0d0
    real(8) :: total_time_ms    = 0.0d0
    real(8) :: avg_iter_time_ms = 0.0d0
  end type newton_stats

contains

  ! -----------------------------------------------------------------
  ! Constrói o PADRAO CSR (hrow_ptr, hcol).  Deve ser chamado UMA VEZ,
  ! porque a esparsidade da Hessiana nao muda entre iteracoes.
  ! Temporarios COO sao passados pelo chamador (alocados uma unica vez).
  ! -----------------------------------------------------------------
  subroutine hess_csr_pattern(x, n, hnnzmax, hrow_ptr, hcol, nnz_out, status, &
                              hval_coo, hrow_coo, hcol_coo, &
                              row_counts, cur_ptr)
    implicit none
    integer,        intent(in)  :: n, hnnzmax
    real(8),        intent(in)  :: x(n)
    integer(c_int), intent(out) :: hrow_ptr(n+1), hcol(:)
    integer,        intent(out) :: nnz_out, status
    real(8),        intent(out) :: hval_coo(:)
    integer,        intent(out) :: hrow_coo(:), hcol_coo(:)
    integer,        intent(out) :: row_counts(:), cur_ptr(:)

    integer :: hnnz_upper, k, r, c, i, idx

    call cutest_ush(status, n, x, hnnz_upper, hnnzmax, &
                    hval_coo, hrow_coo, hcol_coo)
    if (status /= 0) then
      nnz_out = 0
      return
    end if

    nnz_out = hnnz_upper
    do k = 1, hnnz_upper
      if (hrow_coo(k) /= hcol_coo(k)) nnz_out = nnz_out + 1
    end do

    row_counts = 0
    do k = 1, hnnz_upper
      row_counts(hrow_coo(k)) = row_counts(hrow_coo(k)) + 1
      if (hrow_coo(k) /= hcol_coo(k)) &
        row_counts(hcol_coo(k)) = row_counts(hcol_coo(k)) + 1
    end do

    hrow_ptr(1) = 1
    do i = 1, n
      hrow_ptr(i+1) = hrow_ptr(i) + row_counts(i)
    end do

    cur_ptr = hrow_ptr(1:n)
    do k = 1, hnnz_upper
      r = hrow_coo(k); c = hcol_coo(k)
      idx = cur_ptr(r); hcol(idx) = int(c, c_int)
      cur_ptr(r) = cur_ptr(r) + 1
      if (r /= c) then
        idx = cur_ptr(c); hcol(idx) = int(r, c_int)
        cur_ptr(c) = cur_ptr(c) + 1
      end if
    end do
  end subroutine hess_csr_pattern

  ! -----------------------------------------------------------------
  ! Recalcula APENAS hval(x) reusando o padrao ja construido.
  ! Chama cutest_ush, e redistribui os valores pela estrutura CSR.
  ! Requer que hrow_coo/hcol_coo tenham a MESMA ordem da chamada anterior
  ! (verdadeiro para CUTEst, cuja estrutura e determinada pelo problema).
  ! -----------------------------------------------------------------
  subroutine hess_csr_values(x, n, hnnzmax, hrow_ptr, hcol, hval, status, &
                             hval_coo, hrow_coo, hcol_coo, &
                             row_counts, cur_ptr)
    implicit none
    integer,        intent(in)  :: n, hnnzmax
    real(8),        intent(in)  :: x(n)
    integer(c_int), intent(in)  :: hrow_ptr(n+1), hcol(:)
    real(8),        intent(out) :: hval(:)
    integer,        intent(out) :: status
    real(8),        intent(out) :: hval_coo(:)
    integer,        intent(out) :: hrow_coo(:), hcol_coo(:)
    integer,        intent(out) :: row_counts(:), cur_ptr(:)

    integer :: hnnz_upper, k, r, c, idx

    call cutest_ush(status, n, x, hnnz_upper, hnnzmax, &
                    hval_coo, hrow_coo, hcol_coo)
    if (status /= 0) return

    cur_ptr = hrow_ptr(1:n)
    do k = 1, hnnz_upper
      r = hrow_coo(k); c = hcol_coo(k)
      idx = cur_ptr(r); hval(idx) = hval_coo(k)
      cur_ptr(r) = cur_ptr(r) + 1
      if (r /= c) then
        idx = cur_ptr(c); hval(idx) = hval_coo(k)
        cur_ptr(c) = cur_ptr(c) + 1
      end if
    end do
  end subroutine hess_csr_values

  ! -----------------------------------------------------------------
  ! Newton-CG truncado, CG na CPU
  ! -----------------------------------------------------------------
  subroutine newton_cg_cpu(x, n, hnnzmax, eps, eps_cg, cgmaxit, maxit, stats)
    implicit none
    integer,      intent(in)    :: n, hnnzmax, cgmaxit, maxit
    real(8),      intent(in)    :: eps, eps_cg
    real(8),      intent(inout) :: x(n)
    type(newton_stats), intent(out) :: stats

    integer :: status, nnz, iter, cg_iter, istop_cg
    real(8) :: f, gnorm
    real(8), allocatable :: g(:), d(:), hval(:)
    integer(c_int), allocatable :: hrow_ptr(:), hcol(:)
    integer(8) :: t0, t1, trate, tc0, tc1
    logical :: pattern_built

    real(8), allocatable :: hval_coo(:)
    integer, allocatable :: hrow_coo(:), hcol_coo(:)
    integer, allocatable :: row_counts(:), cur_ptr(:)

    allocate(g(n), d(n), hval(2*hnnzmax), hrow_ptr(n+1), hcol(2*hnnzmax))
    allocate(hval_coo(hnnzmax), hrow_coo(hnnzmax), hcol_coo(hnnzmax))
    allocate(row_counts(n), cur_ptr(n))

    call system_clock(count_rate=trate)
    call system_clock(t0)

    pattern_built = .false.

    do iter = 0, maxit
      call cutest_ufn(status, n, x, f)
      stats%fcnt = stats%fcnt + 1
      if (status /= 0) then; stats%istop = -1; exit; end if

      call cutest_ugr(status, n, x, g)
      stats%gcnt = stats%gcnt + 1
      if (status /= 0) then; stats%istop = -1; exit; end if

      gnorm = norm2(g(1:n))
      if (gnorm <= eps) then; stats%istop = 0; exit; end if
      if (iter >= maxit) then; stats%istop = 1; exit; end if

      ! --- Padrao CSR: construido UMA VEZ (nao conta como avaliacao) ---
      if (.not. pattern_built) then
        call hess_csr_pattern(x, n, hnnzmax, hrow_ptr, hcol, nnz, status, &
                              hval_coo, hrow_coo, hcol_coo, &
                              row_counts, cur_ptr)
        if (status /= 0) then; stats%istop = -1; exit; end if
        pattern_built = .true.
      end if

      ! --- Valores hval(x): recalculados a cada iteracao (conta) ---
      call hess_csr_values(x, n, hnnzmax, hrow_ptr, hcol, hval, status, &
                           hval_coo, hrow_coo, hcol_coo, &
                           row_counts, cur_ptr)
      stats%hcnt = stats%hcnt + 1
      if (status /= 0) then; stats%istop = -1; exit; end if

      d(1:n) = 0.0d0

      call system_clock(tc0)
      call cpu_conjugated_gradient(cg_iter, istop_cg, eps_cg, g, cgmaxit, &
                                   n, nnz, d, hval, hcol, hrow_ptr)
      call system_clock(tc1)
      stats%cg_time_ms = stats%cg_time_ms + &
                         real(tc1 - tc0, 8) / real(trate, 8) * 1.0d3

      stats%last_cg_iter = cg_iter
      if (istop_cg == 1) then
        stats%cg_fails = stats%cg_fails + 1
        if (stats%first_cg_fail < 0) stats%first_cg_fail = iter + 1
      end if

      x(1:n) = x(1:n) + d(1:n)
      stats%newton_iter = stats%newton_iter + 1
    end do

    call system_clock(t1)
    stats%total_time_ms = real(t1 - t0, 8) / real(trate, 8) * 1000.0d0
    if (stats%newton_iter > 0) &
      stats%avg_iter_time_ms = stats%total_time_ms / real(stats%newton_iter, 8)

    deallocate(g, d, hval, hrow_ptr, hcol)
    deallocate(hval_coo, hrow_coo, hcol_coo, row_counts, cur_ptr)
  end subroutine newton_cg_cpu

  ! -----------------------------------------------------------------
  ! Newton-CG truncado, CG na GPU
  ! -----------------------------------------------------------------
  subroutine newton_cg_gpu(x, n, hnnzmax, eps, eps_cg, cgmaxit, maxit, &
                           stats, ctx)
    implicit none
    integer,      intent(in)    :: n, hnnzmax, cgmaxit, maxit
    real(8),      intent(in)    :: eps, eps_cg
    real(8),      intent(inout) :: x(n)
    type(newton_stats), intent(out) :: stats
    type(rocm_cg_context), intent(inout) :: ctx

    integer :: status, nnz, iter, cg_iter, istop_cg
    real(8) :: f, gnorm
    real(8), allocatable, target :: g(:), d(:), hval(:)
    integer(c_int), allocatable, target :: hrow_ptr(:), hcol(:)
    integer(8) :: t0, t1, trate, tc0, tc1
    logical :: pattern_built

    real(8), allocatable :: hval_coo(:)
    integer, allocatable :: hrow_coo(:), hcol_coo(:)
    integer, allocatable :: row_counts(:), cur_ptr(:)

    allocate(g(n), d(n), hval(2*hnnzmax), hrow_ptr(n+1), hcol(2*hnnzmax))
    allocate(hval_coo(hnnzmax), hrow_coo(hnnzmax), hcol_coo(hnnzmax))
    allocate(row_counts(n), cur_ptr(n))

    call system_clock(count_rate=trate)
    call system_clock(t0)

    pattern_built = .false.

    do iter = 0, maxit
      call cutest_ufn(status, n, x, f)
      stats%fcnt = stats%fcnt + 1
      if (status /= 0) then; stats%istop = -1; exit; end if

      call cutest_ugr(status, n, x, g)
      stats%gcnt = stats%gcnt + 1
      if (status /= 0) then; stats%istop = -1; exit; end if

      gnorm = norm2(g(1:n))
      if (gnorm <= eps) then; stats%istop = 0; exit; end if
      if (iter >= maxit) then; stats%istop = 1; exit; end if

      ! --- Padrao CSR: construido e copiado UMA VEZ (nao conta) ---
      if (.not. pattern_built) then
        call hess_csr_pattern(x, n, hnnzmax, hrow_ptr, hcol, nnz, status, &
                              hval_coo, hrow_coo, hcol_coo, &
                              row_counts, cur_ptr)
        if (status /= 0) then; stats%istop = -1; exit; end if

        call hipCheck(hipMemcpy(ctx%d_hrow_ptr, c_loc(hrow_ptr), &
                                (n+1)*4_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(ctx%d_hcol, c_loc(hcol), &
                                nnz*4_c_size_t, hipMemcpyHostToDevice))
        pattern_built = .true.
      end if

      ! --- Valores hval(x): recalculados a cada iteracao (conta) ---
      call hess_csr_values(x, n, hnnzmax, hrow_ptr, hcol, hval, status, &
                           hval_coo, hrow_coo, hcol_coo, &
                           row_counts, cur_ptr)
      stats%hcnt = stats%hcnt + 1
      if (status /= 0) then; stats%istop = -1; exit; end if

      ! --- Copia SOMENTE hval e o RHS; zera d ---
      call hipCheck(hipMemcpy(ctx%d_hval, c_loc(hval), &
                              nnz*8_c_size_t, hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(ctx%d_r, c_loc(g), &
                              n*8_c_size_t, hipMemcpyHostToDevice))
      call hipCheck(hipMemset(ctx%d_d, 0, n*8_c_size_t))
      call hipCheck(hipDeviceSynchronize())

      call system_clock(tc0)
      call gpu_conjugated_gradient(cg_iter, istop_cg, eps_cg, gnorm, cgmaxit, &
                                   n, nnz, ctx)
      call hipCheck(hipDeviceSynchronize())
      call system_clock(tc1)
      stats%cg_time_ms = stats%cg_time_ms + &
                         real(tc1 - tc0, 8) / real(trate, 8) * 1.0d3

      stats%last_cg_iter = cg_iter
      if (istop_cg == 1) then
        stats%cg_fails = stats%cg_fails + 1
        if (stats%first_cg_fail < 0) stats%first_cg_fail = iter + 1
      end if

      call hipCheck(hipMemcpy(c_loc(d), ctx%d_d, n*8_c_size_t, &
                              hipMemcpyDeviceToHost))
      x(1:n) = x(1:n) + d(1:n)
      stats%newton_iter = stats%newton_iter + 1
    end do

    call system_clock(t1)
    stats%total_time_ms = real(t1 - t0, 8) / real(trate, 8) * 1000.0d0
    if (stats%newton_iter > 0) &
      stats%avg_iter_time_ms = stats%total_time_ms / real(stats%newton_iter, 8)

    deallocate(g, d, hval, hrow_ptr, hcol)
    deallocate(hval_coo, hrow_coo, hcol_coo, row_counts, cur_ptr)
  end subroutine newton_cg_gpu

end module newton_cg_mod