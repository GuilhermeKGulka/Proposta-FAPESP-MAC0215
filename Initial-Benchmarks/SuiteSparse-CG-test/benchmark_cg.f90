module benchmark_cg_mod
  use iso_c_binding
  use rocm_context
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipfort_hipBLAS
  implicit none

  private
  public :: gpu_conjugated_gradient, cpu_conjugated_gradient

  interface
    subroutine launch_compute_alpha(r2, dot_pHp, alpha) bind(C, name="launch_compute_alpha")
      use iso_c_binding
      type(c_ptr), value :: r2, dot_pHp, alpha
    end subroutine launch_compute_alpha

    subroutine launch_compute_beta(r2, lastr2, beta) bind(C, name="launch_compute_beta")
      use iso_c_binding
      type(c_ptr), value :: r2, lastr2, beta
    end subroutine launch_compute_beta

    subroutine launch_update_p_vector(p_n, p, r, beta) bind(C, name="launch_update_p_vector")
      use iso_c_binding
      integer(c_int), value :: p_n
      type(c_ptr), value :: p, r, beta
    end subroutine launch_update_p_vector
  end interface

contains

! *****************************************************************
! *****************************************************************

  subroutine gpu_conjugated_gradient(cg_iter, istop, eps, gnorm, cgmaxit, n, hnnz, ctx)
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(inout) :: cg_iter, istop
    real(kind=8), intent(in) :: eps, gnorm
    integer, intent(in) :: cgmaxit, n, hnnz

    ! LOCAL SCALARS
    real(kind=8), target :: rnorm

    ! HIP context
    type(rocm_cg_context), intent(in) :: ctx
    
    !Numerical Optimization, Nocedal & Wright, pg.111, Algorithm 5.2

    cg_iter = 0

    ! r = Hd + g !
    ! spmv: r = 1*Hd + 1*r : pois iniciamos r como copia de g
    !
    call hipsparseCheck(hipsparseDcsrmv(ctx%hsparse_handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, &
                                    n, n, hnnz, &
                                    1.0d0, ctx%descrH, ctx%d_hval, ctx%d_hrow_ptr, ctx%d_hcol, &
                                    ctx%d_d, &
                                    1.0d0, ctx%d_r))
    !
    !!!!!!!!!!!!!!

    ! p = -r !
    !
    call hipCheck(hipMemset(ctx%d_p, 0, n * 8_c_size_t))
    !
    ! p = -r + beta*p (mas p = 0)
    call launch_update_p_vector(int(n, c_int), ctx%d_p, ctx%d_r, ctx%d_beta)
    !
    !!!!!!!!!!
            
    ! r2 = r*r !
    call hipblasCheck(hipblasDotEx(ctx%hblas_handle, n, ctx%d_r, HIP_R_64F, 1, &
                                   ctx%d_r, HIP_R_64F, 1, ctx%d_r2, HIP_R_64F, HIP_R_64F))
    !!!!!!!!!!!!

    ! lastr2 = r2 !
    call hipCheck(hipMemcpy(ctx%d_lastr2, ctx%d_r2, 8_c_size_t, hipMemcpyDeviceToDevice))
    !!!!!!!!!!!!!!!!

  10 continue

    ! rnorm = norm2( r ) !!!!
    call hipblasCheck(hipblasNrm2Ex(ctx%hblas_handle, n, ctx%d_r, HIP_R_64F, 1, &
                                    ctx%d_rnorm, HIP_R_64F, HIP_R_64F))
    !!!!!!!!!!!!!!!!!!!!!!!!!
    
    call hipCheck(hipMemcpy(c_loc(rnorm), ctx%d_rnorm, 8_c_size_t, hipMemcpyDeviceToHost))

    !write(*,*) 'cg iter = ',cg_iter,' rnorm = ',rnorm
    
    if ( rnorm .le. eps*max(1.0d0,gnorm) ) then
      !write(*,*) 'Residual norm smaller than the required tolerance.', cg_iter,':',rnorm,':',eps*max(1.0d0,gnorm)
      istop = 0
      return
    end if

    if ( cg_iter .ge. cgmaxit ) then
      !write(*,*) 'Maximum of cg iter reached.', cg_iter,':',rnorm,':',eps*max(1.0d0,gnorm)
      istop = 1
      return
    end if

    cg_iter = cg_iter + 1

    ! Hp = H*p !
    ! spmv: Hp = 1*H*p + 0*Hp
    !
    call hipsparseCheck(hipsparseDcsrmv(ctx%hsparse_handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, &
                                        n, n, hnnz, &
                                        1.0d0, ctx%descrH, ctx%d_hval, ctx%d_hrow_ptr, ctx%d_hcol, &
                                        ctx%d_p, &
                                        0.0d0, ctx%d_Hp))
    !
    !!!!!!!!!!!!
                             
    ! dot_pHp = p*Hp !
    call hipblasCheck(hipblasDotEx(ctx%hblas_handle, n, ctx%d_p, HIP_R_64F, 1, &
                                   ctx%d_Hp, HIP_R_64F, 1, ctx%d_dot_pHp, HIP_R_64F, HIP_R_64F))
    !!!!!!!!!!!!!!!!!!

    ! alpha = lastr2 / dot_pHp !
    call launch_compute_alpha(ctx%d_lastr2, ctx%d_dot_pHp, ctx%d_alpha)

    ! d = alpha*p + d !
    call hipblasCheck(hipblasAxpyEx(ctx%hblas_handle, n, ctx%d_alpha, HIP_R_64F, &
                                    ctx%d_p, HIP_R_64F, 1, ctx%d_d, HIP_R_64F, 1, HIP_R_64F))
    !!!!!!!!!!!!!!!!!!!
                                    
    ! r = alpha*Hp + r !
    call hipblasCheck(hipblasAxpyEx(ctx%hblas_handle, n, ctx%d_alpha, HIP_R_64F, &
                                    ctx%d_Hp, HIP_R_64F, 1, ctx%d_r, HIP_R_64F, 1, HIP_R_64F))
    !!!!!!!!!!!!!!!!!!!!

    ! r2 = r*r !
    call hipblasCheck(hipblasDotEx(ctx%hblas_handle, n, ctx%d_r, HIP_R_64F, 1, &
                                   ctx%d_r, HIP_R_64F, 1, ctx%d_r2, HIP_R_64F, HIP_R_64F))
    !!!!!!!!!!!!

    !beta = r2 / lastr2     E     lastr2 = r2
    call launch_compute_beta(ctx%d_r2, ctx%d_lastr2, ctx%d_beta)  

    ! p = -r + beta*p (CHAMADA DO KERNEL FUNDIDO)
    call launch_update_p_vector(int(n, c_int), ctx%d_p, ctx%d_r, ctx%d_beta)

    go to 10

  end subroutine gpu_conjugated_gradient

  subroutine cpu_conjugated_gradient(cg_iter, istop, eps, g, cgmaxit, n, hnnz, d, hval, hcol, hrow_ptr)
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(inout) :: cg_iter, istop
    real(kind=8), intent(in) :: eps
    integer, intent(in) :: cgmaxit, n, hnnz

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: d(n)
    real(kind=8), intent(in) :: g(n), hval(:)
    integer(c_int), intent(in) :: hrow_ptr(n+1), hcol(:)

    ! LOCAL SCALARS
    real(kind=8) :: gnorm, rnorm, r2, lastr2, p_dot_Hp, alpha, beta, val
    integer :: row, k
    
    ! LOCAL ARRAYS
    real(kind=8), allocatable :: p(:), r(:), Hd(:), Hp(:)

    allocate(p(n), r(n), Hd(n), Hp(n))
    
    cg_iter = 0

    gnorm = norm2(g)

    ! Calculo Hd
    Hd(1:n) = 0.0d0

    r(1:n) = g(1:n) + Hd(1:n)
    p(1:n) = -r(1:n)
    r2 = dot_product(r, r)
    lastr2 = r2

    do
      rnorm = norm2(r)
      
      if (rnorm .le. eps * max(1.0d0, gnorm)) then
        !write(*,*) 'Residual norm smaller than the required tolerance.', cg_iter,':',rnorm,':',eps*max(1.0d0,gnorm)
        istop = 0
        exit
      end if
      
      if (cg_iter .ge. cgmaxit) then
        !write(*,*) 'Maximum of cg iter reached.', cg_iter,':',rnorm,':',eps*max(1.0d0,gnorm)
        istop = 1
        exit
      end if
      
      cg_iter = cg_iter + 1

      Hp(1:n) = 0.0d0
      do row = 1, n
        val = 0.0d0
        do k = hrow_ptr(row), hrow_ptr(row+1) - 1
          val = val + hval(k) * p(hcol(k))
        end do
        Hp(row) = val
      end do

      p_dot_Hp = dot_product(p, Hp)
      alpha = lastr2 / p_dot_Hp

      d(1:n) = d(1:n) + alpha * p(1:n)
      r(1:n) = r(1:n) + alpha * Hp(1:n)
      r2 = dot_product(r, r)
      beta = r2 / lastr2
      lastr2 = r2
      p(1:n) = -r(1:n) + beta * p(1:n)
    end do

    deallocate(p, r, Hd, Hp)

  end subroutine cpu_conjugated_gradient

end module benchmark_cg_mod

program benchmarkma
  use benchmark_cg_mod
  use iso_c_binding
  use rocm_context
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipfort_hipBLAS
  use mmio
  implicit none

  character(len=256) :: filename, problem_name
  character(len=32) :: rep, field, symmetry
  integer :: unit = 10, io_status
  integer :: n, m, nnz_file, nnz_total, i, j, k, r, c, idx
  real(8) :: val, val_im

  ! Arrays COO temporarios
  integer, allocatable :: coo_row(:), coo_col(:), row_counts(:), cur_ptr(:)
  real(8), allocatable :: coo_val(:)

  ! Arrays CSR
  integer(c_int), allocatable, target :: hrow_ptr(:), hcol(:)
  real(8), allocatable, target :: hval(:), d(:), g(:)

  ! Vetores GPU
  real(8), allocatable :: d_gpu(:)

  ! Contexto GPU
  type(rocm_cg_context) :: ctx

  ! Vetores CPU
  real(8), allocatable :: d_cpu(:), p_cpu(:), r_cpu(:), Hp_cpu(:)

  ! Variaveis de controle do CG
  real(8), parameter :: eps = 1.0d-08
  real(8) :: gnorm
  integer :: cg_iter_gpu, cg_iter_cpu, cgmaxit, istop_gpu, istop_cpu

  ! Variaveis de benchmark
  integer(8) :: t_start, t_end, rate
  real(8) :: t_gpu_copy, t_gpu_cg, t_cpu_cg
  real(8) :: speedup, speedup_with_copy, density, iter_diff_rel

  call system_clock(count_rate=rate)
  
  ! LER ARGUMENTO DA LINHA DE COMANDO
  if (command_argument_count() < 1) then
    write(*,*) "Uso: ./benchmark <arquivo.mtx>"
    stop
  end if
  call get_command_argument(1, filename)
  
  ! Extrair nome do problema
  i = index(filename, '/', .true.)
  j = index(filename, '.mtx', .true.)
  if (i == 0) i = 0
  if (j == 0) j = len_trim(filename) + 1
  problem_name = filename(i+1:j-1)
  
  ! LER ARQUIVO MATRIX MARKET VIA MODULO NIST MMIO
  open(unit=unit, file=trim(filename), status='old', action='read', iostat=io_status)
  if (io_status /= 0) then
    write(*,*) "Erro ao abrir arquivo: ", trim(filename)
    stop
  end if

  ! Ler o banner do MatrixMarket
  call mm_read_banner(unit, rep, field, symmetry, io_status)
  if (io_status /= 0) then
    write(*,*) "Erro ao ler o cabeçalho MatrixMarket de: ", trim(filename)
    close(unit)
    stop
  end if

  ! Ler as dimensoes (linhas, colunas, nnz)
  call mm_read_mtx_crd_size(unit, m, n, nnz_file, io_status)
  if (io_status /= 0) then
    write(*,*) "Erro ao ler dimensoes da matriz: ", trim(filename)
    close(unit)
    stop
  end if

  ! Aloca o dobro de espaco para prevencao caso precise espelhar simetria
  allocate(coo_row(nnz_file*2), coo_col(nnz_file*2), coo_val(nnz_file*2))
  
  nnz_total = 0
  do k = 1, nnz_file
    if (trim(field) == 'pattern') then
      read(unit, *, iostat=io_status) r, c
      val = 1.0d0
    else if (trim(field) == 'complex') then
      read(unit, *, iostat=io_status) r, c, val, val_im
    else
      read(unit, *, iostat=io_status) r, c, val
    end if
    
    if (io_status /= 0) exit
    
    nnz_total = nnz_total + 1
    coo_row(nnz_total) = r
    coo_col(nnz_total) = c
    coo_val(nnz_total) = val
    
    ! Espelha para matriz completa apenas se for simetrica/skew-symmetric
    if ((trim(symmetry) == 'symmetric' .or. trim(symmetry) == 'skew-symmetric' .or. &
         trim(symmetry) == 'hermitian') .and. r /= c) then
      nnz_total = nnz_total + 1
      coo_row(nnz_total) = c
      coo_col(nnz_total) = r
      if (trim(symmetry) == 'skew-symmetric') then
        coo_val(nnz_total) = -val
      else
        coo_val(nnz_total) = val
      end if
    end if
  end do
  close(unit)
  
  ! CONVERTER COO PARA CSR (Via Counting Sort O(NNZ))
  allocate(hrow_ptr(n+1), hcol(nnz_total), hval(nnz_total))
  allocate(row_counts(n), cur_ptr(n))
  
  row_counts = 0
  do k = 1, nnz_total
    row_counts(coo_row(k)) = row_counts(coo_row(k)) + 1
  end do
  
  hrow_ptr(1) = 1
  do i = 1, n
    hrow_ptr(i+1) = hrow_ptr(i) + row_counts(i)
  end do
  
  cur_ptr = hrow_ptr(1:n)
  do k = 1, nnz_total
    r = coo_row(k)
    idx = cur_ptr(r)
    hcol(idx) = int(coo_col(k), c_int)
    hval(idx) = coo_val(k)
    cur_ptr(r) = cur_ptr(r) + 1
  end do
  
  deallocate(coo_row, coo_col, coo_val, row_counts, cur_ptr)

  ! CONFIGURARACOES INICIAIS DO METODO
  allocate(d(n), g(n))
  g(1:n) = -1.0d0
  d(1:n) = 0.0d0
  gnorm = maxval(abs(g(1:n)))
  cgmaxit = 2 * n

  ! INICIALIZAR ROCM CONTEXT E COPIA DE MEMORIA
  call init_rocm_context(ctx, n, nnz_total)
  call hipblasCheck(hipblasSetPointerMode(ctx%hblas_handle, 1))

  call system_clock(t_start)
  call hipCheck(hipMemcpy(ctx%d_hrow_ptr, c_loc(hrow_ptr), (n+1)*4_c_size_t, hipMemcpyHostToDevice))
  call hipCheck(hipMemcpy(ctx%d_hcol, c_loc(hcol), nnz_total*4_c_size_t, hipMemcpyHostToDevice))
  call hipCheck(hipMemcpy(ctx%d_hval, c_loc(hval), nnz_total*8_c_size_t, hipMemcpyHostToDevice))
  call hipCheck(hipMemcpy(ctx%d_r, c_loc(g), n*8_c_size_t, hipMemcpyHostToDevice))
  call hipCheck(hipMemcpy(ctx%d_d, c_loc(d), n*8_c_size_t, hipMemcpyHostToDevice))
  call hipCheck(hipDeviceSynchronize())
  call system_clock(t_end)
  t_gpu_copy = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0

  call system_clock(t_start)
  call gpu_conjugated_gradient(cg_iter_gpu, istop_gpu, eps, gnorm, cgmaxit, n, nnz_total, ctx)
  call system_clock(t_end)
  t_gpu_cg = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0

  call hipCheck(hipMemcpy(ctx%d_d, c_loc(d), n*8_c_size_t, hipMemcpyHostToDevice))

  ! BENCHMARK CPU CG
  allocate(d_cpu(n), p_cpu(n), r_cpu(n), Hp_cpu(n))
  d_cpu(1:n) = d(1:n)
  
  call system_clock(t_start)
  call cpu_conjugated_gradient(cg_iter_cpu, istop_cpu, eps, g, cgmaxit, n, nnz_total, d_cpu, hval, hcol, hrow_ptr)
  call system_clock(t_end)
  t_cpu_cg = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0

  ! METRICAS E ESCRITA
  speedup = t_cpu_cg / t_gpu_cg
  speedup_with_copy = t_cpu_cg / (t_gpu_copy + t_gpu_cg)
  density = real(nnz_total, 8) / (real(n, 8) * real(n, 8))
  iter_diff_rel = real(abs(cg_iter_gpu - cg_iter_cpu),8) / (real(max(1, cg_iter_cpu),8))

  open(unit=20, file='results_cg.csv', status='unknown', position='append', action='write')
  write(20, '(A, ",", I0, ",", I0, ",", ES12.5, ",", I0, ",", I0, ",", I0, ",", I0, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5)') &
        trim(problem_name), n, nnz_total, density, cg_iter_gpu, cg_iter_cpu, istop_gpu, istop_cpu, iter_diff_rel, &
        t_gpu_copy, t_gpu_cg, t_cpu_cg, speedup, speedup_with_copy
  close(20)

  call finalize_rocm_context(ctx)

  deallocate(hrow_ptr, hcol, hval, d, g, d_cpu, p_cpu, r_cpu, Hp_cpu)
  
end program benchmarkma