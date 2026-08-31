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
  integer, parameter :: nnz_per_row = 1 ! Exatamente 3 elementos por linha
  integer :: i, j, k, iter, idx, tmp, niter = 100
  integer(8) :: t_start, t_end, rate
  real(8) :: t_cpu, t_kernel, wall_time, u(3), val_rand
  integer(4) :: cols(3)
  integer(c_int) :: hnnz
  
  ! Arrays host
  integer(c_int), allocatable, target :: hrow_ptr(:), hcol(:)
  real(8), allocatable, target :: hval(:), h_x(:), h_y(:)
  
  ! Pointers device (GPU)
  type(c_ptr) :: handle, descr
  type(c_ptr) :: d_row_ptr = c_null_ptr
  type(c_ptr) :: d_col     = c_null_ptr
  type(c_ptr) :: d_val     = c_null_ptr
  type(c_ptr) :: d_x       = c_null_ptr
  type(c_ptr) :: d_y       = c_null_ptr
  
  integer(c_size_t) :: bytes_row, bytes_col, bytes_val, bytes_vec
  
  call system_clock(count_rate=rate)
  
  ! 1. Tamanho fixo e exato
  hnnz = int(n, c_int) * nnz_per_row
  allocate(hrow_ptr(n+1), hcol(hnnz), hval(hnnz), h_x(n), h_y(n))
  h_x(1:n) = 1.0d0
  
  write(*,*) "Iniciando teste com n = ", n, " | Total NNZ = ", hnnz
  
  ! ------------------------------------------------------------------
  ! 2. GERACAO DA MATRIZ (CSR Direto com 3 colunas aleatórias/linha)
  ! ------------------------------------------------------------------
  idx = 0
  do i = 1, n
    hrow_ptr(i) = idx + 1
    
    ! Sorteia 3 colunas aleatórias para a linha i
    call random_number(u)
    cols(1) = min(n, max(1, int(u(1) * real(n, 8)) + 1))
    cols(2) = min(n, max(1, int(u(2) * real(n, 8)) + 1))
    cols(3) = min(n, max(1, int(u(3) * real(n, 8)) + 1))
    
    ! Ordena os 3 indices para manter padrão CSR válido
    if (cols(1) > cols(2)) then; tmp = cols(1); cols(1) = cols(2); cols(2) = tmp; end if
    if (cols(2) > cols(3)) then; tmp = cols(2); cols(2) = cols(3); cols(3) = tmp; end if
    if (cols(1) > cols(2)) then; tmp = cols(1); cols(1) = cols(2); cols(2) = tmp; end if
    
    do k = 1, nnz_per_row
      idx = idx + 1
      hcol(idx) = int(cols(k), c_int)
      call random_number(val_rand)
      hval(idx) = val_rand
    end do
  end do
  hrow_ptr(n+1) = idx + 1
  
  ! ------------------------------------------------------------------
  ! 3. BENCHMARK CPU (SpMV CSR)
  ! ------------------------------------------------------------------
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
    if (wall_time .le. t_cpu) t_cpu = wall_time 
  end do
  write(*,*) "  CPU min Time: ", t_cpu * 1000.0d0, " ms"
  
  ! ------------------------------------------------------------------
  ! 4. SETUP E BENCHMARK GPU (hipSPARSE)
  ! ------------------------------------------------------------------
  bytes_row = int(n + 1, c_size_t) * 4_c_size_t
  bytes_col = int(hnnz, c_size_t) * 4_c_size_t
  bytes_val = int(hnnz, c_size_t) * 8_c_size_t
  bytes_vec = int(n, c_size_t) * 8_c_size_t
  
  call hipcheck(hipMalloc(d_row_ptr, bytes_row))
  call hipcheck(hipMalloc(d_col, bytes_col))
  call hipcheck(hipMalloc(d_val, bytes_val))
  call hipcheck(hipMalloc(d_x, bytes_vec))
  call hipcheck(hipMalloc(d_y, bytes_vec))
  
  call hipcheck(hipMemcpy(d_row_ptr, c_loc(hrow_ptr), bytes_row, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_col, c_loc(hcol), bytes_col, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_val, c_loc(hval), bytes_val, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_x, c_loc(h_x), bytes_vec, hipMemcpyHostToDevice))
  
  call hipsparsecheck(hipsparseCreate(handle))
  call hipsparsecheck(hipsparseCreateMatDescr(descr))
  call hipsparsecheck(hipsparseSetMatIndexBase(descr, HIPSPARSE_INDEX_BASE_ONE))
  
  ! Warm-up (Força compilação do JIT e ativação da GPU)
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
    if (wall_time .le. t_kernel) t_kernel = wall_time
  end do
  write(*,*) "  GPU min Time: ", t_kernel * 1000.0d0, " ms"
  
  ! ------------------------------------------------------------------
  ! 5. LIMPEZA
  ! ------------------------------------------------------------------
  call hipcheck(hipFree(d_row_ptr))
  call hipcheck(hipFree(d_col))
  call hipcheck(hipFree(d_val))
  call hipcheck(hipFree(d_x))
  call hipcheck(hipFree(d_y))
  call hipsparsecheck(hipsparseDestroy(handle))
  call hipsparsecheck(hipsparseDestroyMatDescr(descr))
  deallocate(hrow_ptr, hcol, hval, h_x, h_y)
  
end subroutine run_test