program test_gpu_performance
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  implicit none
  
  integer, parameter :: n_sizes(6) = [100, 1000, 10000, 100000, 1000000, 10000000]
  integer :: test
  
  do test = 1, size(n_sizes)
    call run_test(n_sizes(test))
    print *, "----------------------------------------"
  end do
  
end program test_gpu_performance

subroutine run_test(n)
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  implicit none
  
  integer, intent(in) :: n
  integer :: i, j, k, iter, idx, niter = 100
  integer(8) :: t_start, t_end, rate
  real(8) :: t_cpu, t_kernel, wall_time
  integer(c_int) :: hnnz
  
  ! arrays host (usando c_int para compatibilidade com HIP)
  integer(c_int), allocatable, target :: hrow_ptr(:), hcol(:)
  real(8), allocatable, target :: hval(:), h_x(:), h_y(:), h_y_ref(:)
  
  ! pointers device
  type(c_ptr) :: handle, descr
  type(c_ptr) :: d_row_ptr = c_null_ptr
  type(c_ptr) :: d_col     = c_null_ptr
  type(c_ptr) :: d_val     = c_null_ptr
  type(c_ptr) :: d_x       = c_null_ptr
  type(c_ptr) :: d_y       = c_null_ptr
  
  integer(c_size_t) :: bytes_row, bytes_col, bytes_val, bytes_vec
  
  ! Inicializa a taxa do relógio do sistema
  call system_clock(count_rate=rate)
  
  hnnz = n!3*n - 2
  allocate(hrow_ptr(n+1), hcol(hnnz), hval(hnnz), h_x(n), h_y(n), h_y_ref(n))
  
  h_x(1:n) = 1.0d0
  write(*,*) "Iniciando teste com n = ", n
  
  ! 1. inicializacao dados host
  idx = 0
  do i = 1, n
    hrow_ptr(i) = idx + 1
    
    ! Subdiagonal
    !if (i > 1) then
    !    idx = idx + 1
    !    hcol(idx) = i - 1
    !    hval(idx) = -400.0d0 * h_x(i-1)
    !end if
    ! Diagonal
    idx = idx + 1
    hcol(idx) = i
    if (i == 1) then
        hval(idx) = 1200.0d0 * h_x(1)**2 - 400.0d0 * h_x(2) + 2.0d0
    else if (i == n) then
        hval(idx) = 200.0d0
    else
        hval(idx) = 1200.0d0 * h_x(i)**2 - 400.0d0 * h_x(i+1) + 202.0d0
    end if
    ! Superdiagonal
    !if (i < n) then
    !    idx = idx + 1
    !    hcol(idx) = i + 1
    !    hval(idx) = -400.0d0 * h_x(i)
    !end if
  end do
  hrow_ptr(n+1) = idx + 1
  
  ! 2. benchmark CPU
  t_cpu = 10.0d30

  do iter = 1, niter
    call system_clock(t_start)

    do i = 1, n
      h_y(i) = 0.0d0
      do k = hrow_ptr(i), hrow_ptr(i+1) - 1
        j = hcol(k)
        h_y(i) = h_y(i) + hval(k) * h_x(j)
      end do
    end do

    call system_clock(t_end)
    wall_time = real(t_end - t_start, 8) / real(rate, 8)
    if ( wall_time .le. t_cpu ) then 
      t_cpu = wall_time 
    endif
  end do
  write(*,*) "  CPU min Time: ", t_cpu * 1000.0d0, " ms"
  
  ! 3. tamanhos em bytes
  bytes_row = int(n + 1, c_size_t) * 4_c_size_t
  bytes_col = int(hnnz, c_size_t) * 4_c_size_t
  bytes_val = int(hnnz, c_size_t) * 8_c_size_t
  bytes_vec = int(n, c_size_t) * 8_c_size_t
  
  ! 4. malloc GPU
  call hipcheck(hipMalloc(d_row_ptr, bytes_row))
  call hipcheck(hipMalloc(d_col, bytes_col))
  call hipcheck(hipMalloc(d_val, bytes_val))
  call hipcheck(hipMalloc(d_x, bytes_vec))
  call hipcheck(hipMalloc(d_y, bytes_vec))
  
  ! 5. copia H->D
  call hipcheck(hipMemcpy(d_row_ptr, c_loc(hrow_ptr), bytes_row, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_col, c_loc(hcol), bytes_col, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_val, c_loc(hval), bytes_val, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_x, c_loc(h_x), bytes_vec, hipMemcpyHostToDevice))
  
  ! 6. configurar hipSPARSE
  call hipsparsecheck(hipsparseCreate(handle))
  call hipsparsecheck(hipsparseCreateMatDescr(descr))
  call hipsparsecheck(hipsparseSetMatIndexBase(descr, HIPSPARSE_INDEX_BASE_ONE))
  
  ! 7. execucao do kernel (GPU Benchmark)
  
  ! Warm-up: Garante que o JIT da GPU compile os kernels necessários antes de cronometrar
  call hipsparsecheck(hipsparseDcsrmv(handle, &
                          HIPSPARSE_OPERATION_NON_TRANSPOSE, &
                          n, n, hnnz, &
                          1.0d0, descr, &
                          d_val, d_row_ptr, d_col, &
                          d_x, &
                          0.0d0, d_y))
  call hipcheck(hipDeviceSynchronize())
                          
  t_kernel = 10.0d30
  
  do iter = 1, niter
    call system_clock(t_start)
    
    call hipsparsecheck(hipsparseDcsrmv(handle, &
                            HIPSPARSE_OPERATION_NON_TRANSPOSE, &
                            n, n, hnnz, &
                            1.0d0, descr, &
                            d_val, d_row_ptr, d_col, &
                            d_x, &
                            0.0d0, d_y))
                            
    call hipcheck(hipDeviceSynchronize())
    call system_clock(t_end)
    
    wall_time = real(t_end - t_start, 8) / real(rate, 8)
    if ( wall_time .le. t_kernel ) then
      t_kernel = wall_time
    endif
  end do
  write(*,*) "  GPU min Time: ", t_kernel * 1000.0d0, " ms"
  
  ! 8. limpeza
  call hipcheck(hipFree(d_row_ptr))
  call hipcheck(hipFree(d_col))
  call hipcheck(hipFree(d_val))
  call hipcheck(hipFree(d_x))
  call hipcheck(hipFree(d_y))
  call hipsparsecheck(hipsparseDestroy(handle))
  call hipsparsecheck(hipsparseDestroyMatDescr(descr))
  deallocate(hrow_ptr, hcol, hval, h_x, h_y, h_y_ref)
  
end subroutine run_test