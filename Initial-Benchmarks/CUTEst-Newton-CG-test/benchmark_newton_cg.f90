program benchmark_newton_cg
  use iso_c_binding
  use rocm_context
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipfort_hipBLAS  
  use cutest_iface
  use newton_cg_mod
  implicit none

  character(len=10) :: pname
  integer :: n, hnnzmax
  integer :: status, allocerr, i
  integer :: cgmaxit, maxit, nnz0
  real(8), parameter :: eps    = 1.0d-08
  real(8)            :: eps_cg

  real(8), allocatable :: x_init(:), x_gpu(:), x_cpu(:), bl(:), bu(:)
  real(8) :: density, speedup, speedup_cg, sol_diff
  type(rocm_cg_context) :: ctx
  type(newton_stats) :: stats_cpu, stats_gpu

  ! -----------------------------------------------------------------
  ! 1) Le problema CUTEst
  ! -----------------------------------------------------------------
  open(10, file='OUTSDIF.d', form='formatted', status='old'); rewind 10

  call cutest_udimen(status, 10, n)
  if (status /= 0) stop 'cutest_udimen falhou'

  allocate(x_init(n), bl(n), bu(n), stat=allocerr)
  if (allocerr /= 0) stop 'falha de alocacao'

  call cutest_usetup(status, 10, 6, 4, n, x_init, bl, bu)
  if (status /= 0) stop 'cutest_usetup falhou'

  call cutest_pname(status, 10, pname)
  call cutest_udimsh(status, hnnzmax)
  close(10)

  if (hnnzmax <= 0) then
    write(*,*) 'Problema ', trim(pname), ' sem Hessiana. Pulando.'
    stop
  end if
  hnnzmax = hnnzmax + n   ! margem de seguranca

  ! -----------------------------------------------------------------
  ! 2) Estimativa do nnz do CSR simetrico (na Hessiana em x0)
  ! -----------------------------------------------------------------
  block
    real(8), allocatable :: hval_coo(:)
    integer, allocatable :: hrow_coo(:), hcol_coo(:)
    integer :: hnnz_upper, k
    allocate(hval_coo(hnnzmax), hrow_coo(hnnzmax), hcol_coo(hnnzmax))
    call cutest_ush(status, n, x_init, hnnz_upper, hnnzmax, &
                    hval_coo, hrow_coo, hcol_coo)
    if (status /= 0) stop 'cutest_ush falhou'
    nnz0 = hnnz_upper
    do k = 1, hnnz_upper
      if (hrow_coo(k) /= hcol_coo(k)) nnz0 = nnz0 + 1
    end do
    deallocate(hval_coo, hrow_coo, hcol_coo)
  end block

  ! -----------------------------------------------------------------
  ! 3) Configuracoes do Newton-CG truncado
  ! -----------------------------------------------------------------
  allocate(x_cpu(n), x_gpu(n), stat=allocerr)
  x_cpu = x_init
  x_gpu = x_init

  maxit   = 100
  cgmaxit = n
  eps_cg  = eps !min(1.0d-01, sqrt(eps))    ! forcing term simples

  ! -----------------------------------------------------------------
  ! 4) Roda Newton-CG na GPU
  ! -----------------------------------------------------------------
  call init_rocm_context(ctx, n, 2*hnnzmax)
  call hipblasCheck(hipblasSetPointerMode(ctx%hblas_handle, 1))

  call newton_cg_gpu(x_gpu, n, hnnzmax, eps, eps_cg, cgmaxit, maxit, &
                     stats_gpu, ctx)

  call finalize_rocm_context(ctx)

  ! -----------------------------------------------------------------
  ! 5) Roda Newton-CG na CPU
  ! -----------------------------------------------------------------
  call newton_cg_cpu(x_cpu, n, hnnzmax, eps, eps_cg, cgmaxit, maxit, &
                     stats_cpu)

  ! -----------------------------------------------------------------
  ! 6) Metricas
  ! -----------------------------------------------------------------
  density   = real(nnz0, 8) / (real(n, 8) * real(n, 8))
  speedup    = stats_cpu%total_time_ms / max(stats_gpu%total_time_ms, 1.0d-30)
  speedup_cg = stats_cpu%cg_time_ms    / max(stats_gpu%cg_time_ms,    1.0d-30)

  !write(*,*)' x_cpu(1:min(10,n)) =', x_cpu(1:min(10,n)), &
  !          ' x_gpu(1:min(10,n)) =', x_gpu(1:min(10,n))

  sol_diff  = maxval(abs(x_cpu(1:n) - x_gpu(1:n)))

  write(*,'(/,A)') '=================================================='
  write(*,'(A,A)') 'Problem: ', trim(pname)
  write(*,'(A,I0,A,I0)') '  n = ', n, '   nnz(CSR) = ', nnz0
  write(*,'(A,ES12.4)') '  density     = ', density
  write(*,'(A,ES12.4)') '  max|x_cpu-x_gpu| = ', sol_diff
  write(*,'(/,A)') '  --- CPU ---'
  write(*,'(A,I0)')    '    newton_iter   = ', stats_cpu%newton_iter
  write(*,'(A,I0,A,I0,A,I0)') '    fcnt/gcnt/hcnt = ', &
       stats_cpu%fcnt, '/', stats_cpu%gcnt, '/', stats_cpu%hcnt
  write(*,'(A,I0)')    '    istop         = ', stats_cpu%istop
  write(*,'(A,I0)')    '    cg_fails      = ', stats_cpu%cg_fails
  write(*,'(A,I0)')    '    first_cg_fail = ', stats_cpu%first_cg_fail
  write(*,'(A,I0)')    '    last_cg_iter  = ', stats_cpu%last_cg_iter
  write(*,'(A,ES12.4)') '    total_ms      = ', stats_cpu%total_time_ms
  write(*,'(A,ES12.4)') '    avg_iter_ms   = ', stats_cpu%avg_iter_time_ms
  write(*,'(/,A)') '  --- GPU ---'
  write(*,'(A,I0)')    '    newton_iter   = ', stats_gpu%newton_iter
  write(*,'(A,I0,A,I0,A,I0)') '    fcnt/gcnt/hcnt = ', &
       stats_gpu%fcnt, '/', stats_gpu%gcnt, '/', stats_gpu%hcnt
  write(*,'(A,I0)')    '    istop         = ', stats_gpu%istop
  write(*,'(A,I0)')    '    cg_fails      = ', stats_gpu%cg_fails
  write(*,'(A,I0)')    '    first_cg_fail = ', stats_gpu%first_cg_fail
  write(*,'(A,I0)')    '    last_cg_iter  = ', stats_gpu%last_cg_iter
  write(*,'(A,ES12.4)') '    total_ms      = ', stats_gpu%total_time_ms
  write(*,'(A,ES12.4)') '    avg_iter_ms   = ', stats_gpu%avg_iter_time_ms
  write(*,'(A,ES12.4)') '  speedup        = ', speedup
  write(*,'(A)') '=================================================='

  ! -----------------------------------------------------------------
  ! 7) CSV
  ! -----------------------------------------------------------------
  open(unit=20, file='results_newton_cg.csv', status='unknown', &
       position='append', action='write')
  write(20, '(A, ",", I0, ",", I0, ",", ES12.5, ",", &
              &I0, ",", I0, ",", &
              &I0, ",", I0, ",", I0, ",", I0, ",", I0, ",", I0, ",", &
              &I0, ",", I0, ",", &
              &I0, ",", I0, ",", I0, ",", I0, ",", &
              &I0, ",", I0, ",", &
              &ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", &
              &ES12.5, ",", ES12.5, ",", &
              &ES12.5, ",", ES12.5, ",", ES12.5)') &
    trim(pname), n, nnz0, density, &
    stats_cpu%newton_iter, stats_gpu%newton_iter, &
    stats_cpu%fcnt, stats_gpu%fcnt, &
    stats_cpu%gcnt, stats_gpu%gcnt, &
    stats_cpu%hcnt, stats_gpu%hcnt, &
    stats_cpu%istop, stats_gpu%istop, &
    stats_cpu%cg_fails, stats_gpu%cg_fails, &
    stats_cpu%first_cg_fail, stats_gpu%first_cg_fail, &
    stats_cpu%last_cg_iter, stats_gpu%last_cg_iter, &
    stats_cpu%total_time_ms, stats_gpu%total_time_ms, &
    stats_cpu%avg_iter_time_ms, stats_gpu%avg_iter_time_ms, &
    stats_cpu%cg_time_ms, stats_gpu%cg_time_ms, &
    speedup, speedup_cg, sol_diff
  close(20)

  deallocate(x_init, x_cpu, x_gpu, bl, bu)

end program benchmark_newton_cg