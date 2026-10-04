program benchmark_suitesparse
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use mmio
  use hipsparse_generic
  implicit none

  
  character(len=256) :: filename, problem_name
  character(len=32) :: rep, field, symmetry
  integer :: unit = 10, io_status
  integer :: n, m, nnz_file, nnz_total, i, j, k, r, c, idx, niter
  real(8) :: val, val_im
  
  ! Arrays temporarios COO
  integer, allocatable :: coo_row(:), coo_col(:), row_counts(:), cur_ptr(:)
  real(8), allocatable :: coo_val(:)
  
  ! Arrays CSR (Host)
  integer(c_int), allocatable, target :: hrow_ptr(:), hcol(:)
  real(8), allocatable, target :: hval(:), h_x(:), h_y(:), h_yg(:)
  
  ! Arrays Device (GPU)
  type(c_ptr) :: handle, mat_a, vec_x, vec_y
  type(c_ptr) :: d_buf = c_null_ptr
  integer(c_size_t) :: buf_size
  real(c_double), target :: alpha_s, beta_s
  logical :: do_prep
  integer(c_int) :: st
  character(len=16) :: arg2
  type(c_ptr) :: d_row_ptr = c_null_ptr, d_col = c_null_ptr, d_val = c_null_ptr
  type(c_ptr) :: d_x = c_null_ptr, d_y = c_null_ptr
  integer(c_size_t) :: bytes_row, bytes_col, bytes_val, bytes_vec
  
  ! Timers e Metricas
  integer(8) :: t_start, t_end, rate
  real(8) :: t_gpu_copy, t_gpu_spmv, t_cpu_spmv, speedup, density
  real(8) :: t_gpu_avg, t_cpu_avg, speedup_with_copy, t_gpu_prep, max_err
  
  call system_clock(count_rate=rate)
  
  ! LER ARGUMENTO DA LINHA DE COMANDO
  if (command_argument_count() < 1) then
    write(*,*) "Uso: ./benchmark_suitesparse <arquivo.mtx> [nopre]"
    stop
  end if
  call get_command_argument(1, filename)
  ! 2o argumento opcional "nopre": pula hipsparseSpMV_preprocess (para comparar)
  do_prep = .true.
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg2)
    if (trim(arg2) == 'nopre') do_prep = .false.
  end if
  
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
  allocate(row_counts(n), cur_ptr(n), h_x(n), h_y(n), h_yg(n))
  
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
  
  ! Vetores de entrada/saida
  h_x = 1.0d0
  h_y = 0.0d0
  
  ! SETUP HIPSPARSE E COPIA DE MEMORIA
  niter = min(10000, n)
  
  bytes_row = int(n + 1, c_size_t) * 4_c_size_t
  bytes_col = int(nnz_total, c_size_t) * 4_c_size_t
  bytes_val = int(nnz_total, c_size_t) * 8_c_size_t
  bytes_vec = int(n, c_size_t) * 8_c_size_t
  
  call hipcheck(hipMalloc(d_row_ptr, bytes_row))
  call hipcheck(hipMalloc(d_col, bytes_col))
  call hipcheck(hipMalloc(d_val, bytes_val))
  call hipcheck(hipMalloc(d_x, bytes_vec))
  call hipcheck(hipMalloc(d_y, bytes_vec))
  
  call hipsparsecheck(hipsparseCreate(handle))
  alpha_s = 1.0d0
  beta_s  = 0.0d0
  
  ! Medir tempo de copia H2D
  call system_clock(t_start)
  call hipcheck(hipMemcpy(d_row_ptr, c_loc(hrow_ptr), bytes_row, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_col, c_loc(hcol), bytes_col, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_val, c_loc(hval), bytes_val, hipMemcpyHostToDevice))
  call hipcheck(hipMemcpy(d_x, c_loc(h_x), bytes_vec, hipMemcpyHostToDevice))
  call hipcheck(hipDeviceSynchronize())
  call system_clock(t_end)
  t_gpu_copy = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0 ! Em ms
  
  ! DESCRITORES GENERICOS (CSR 1-based, indices int32, valores double)
  call hs_check(hs_create_csr(mat_a, int(n, c_int64_t), int(n, c_int64_t), &
                int(nnz_total, c_int64_t), d_row_ptr, d_col, d_val, &
                HS_INDEX_32I, HS_INDEX_32I, HS_BASE_ONE, HS_R_64F), "CreateCsr")
  call hs_check(hs_create_dnvec(vec_x, int(n, c_int64_t), d_x, HS_R_64F), "CreateDnVec x")
  call hs_check(hs_create_dnvec(vec_y, int(n, c_int64_t), d_y, HS_R_64F), "CreateDnVec y")

  ! BUFFER EXTERNO
  call hs_check(hs_spmv_buffer_size(handle, HS_OP_N, c_loc(alpha_s), mat_a, vec_x, &
                c_loc(beta_s), vec_y, HS_R_64F, HS_SPMV_ALG_DEFAULT, buf_size), "SpMV_bufferSize")
  call hipcheck(hipMalloc(d_buf, max(buf_size, 1_c_size_t)))

  ! PREPROCESS (analise da matriz, uma vez; cronometrado a parte)
  t_gpu_prep = 0.0d0
  if (do_prep) then
    call system_clock(t_start)
    st = hs_spmv_preprocess(handle, HS_OP_N, c_loc(alpha_s), mat_a, vec_x, &
                  c_loc(beta_s), vec_y, HS_R_64F, HS_SPMV_ALG_DEFAULT, d_buf)
    call hipcheck(hipDeviceSynchronize())
    call system_clock(t_end)
    t_gpu_prep = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0
    if (st /= 0) then   ! preprocess e opcional: avisa e segue sem ele
      write(*,'(A,I0,A)') "AVISO: SpMV_preprocess falhou (codigo ", st, "); seguindo sem ele"
      t_gpu_prep = -1.0d0
    end if
  end if

  ! BENCHMARK GPU
  ! Warmup
  call hs_check(hs_spmv(handle, HS_OP_N, c_loc(alpha_s), mat_a, vec_x, c_loc(beta_s), vec_y, &
                HS_R_64F, HS_SPMV_ALG_DEFAULT, d_buf), "SpMV (warmup)")
  call hipcheck(hipDeviceSynchronize())

  ! Loop Principal GPU
  call system_clock(t_start)
  do i = 1, niter
    call hs_check(hs_spmv(handle, HS_OP_N, c_loc(alpha_s), mat_a, vec_x, c_loc(beta_s), vec_y, &
                  HS_R_64F, HS_SPMV_ALG_DEFAULT, d_buf), "SpMV")
  end do
  call hipcheck(hipDeviceSynchronize())
  call system_clock(t_end)
  t_gpu_spmv = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0

  ! Resultado da GPU de volta ao host (fora dos timers) para checagem
  call hipcheck(hipMemcpy(c_loc(h_yg), d_y, bytes_vec, hipMemcpyDeviceToHost))
  
  ! BENCHMARK CPU
  call system_clock(t_start)
  do i = 1, niter
    do r = 1, n
      val = 0.0d0
      do k = hrow_ptr(r), hrow_ptr(r+1) - 1
        val = val + hval(k) * h_x(hcol(k))
      end do
      h_y(r) = val
    end do
  end do
  call system_clock(t_end)
  t_cpu_spmv = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0
  
  ! CHECAGEM DE CORRETUDE (GPU vs CPU), so aviso na tela
  max_err = maxval(abs(h_y - h_yg)) / max(maxval(abs(h_y)), tiny(1.0d0))
  if (max_err > 1.0d-8) write(*,'(A,A,A,ES10.3)') "AVISO: ", trim(problem_name), &
        " erro relativo GPU vs CPU = ", max_err

  ! CALCULO DAS NOVAS METRICAS
  t_gpu_avg = t_gpu_spmv / real(niter, 8)
  t_cpu_avg = t_cpu_spmv / real(niter, 8)
  speedup = t_cpu_spmv / t_gpu_spmv
  speedup_with_copy = t_cpu_spmv / (t_gpu_copy + t_gpu_spmv)
  density = real(nnz_total, 8) / (real(n, 8) * real(n, 8))
  
  ! ESCRITA NA TABELA CSV
  open(unit=20, file='results.csv', status='unknown', position='append', action='write')
  write(20, '(A, ",", I0, ",", I0, ",", ES12.5, ",", I0, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5)') &
        trim(problem_name), n, nnz_total, density, niter, t_gpu_copy, t_gpu_spmv, t_gpu_avg, t_cpu_spmv, t_cpu_avg, speedup, speedup_with_copy, t_gpu_prep
  close(20)
  
  ! Limpeza
  call hipcheck(hipFree(d_row_ptr))
  call hipcheck(hipFree(d_col))
  call hipcheck(hipFree(d_val))
  call hipcheck(hipFree(d_x))
  call hipcheck(hipFree(d_y))
  call hipcheck(hipFree(d_buf))
  call hs_check(hs_destroy_spmat(mat_a), "DestroySpMat")
  call hs_check(hs_destroy_dnvec(vec_x), "DestroyDnVec x")
  call hs_check(hs_destroy_dnvec(vec_y), "DestroyDnVec y")
  call hipsparsecheck(hipsparseDestroy(handle))
  deallocate(hrow_ptr, hcol, hval, h_x, h_y, h_yg)

end program benchmark_suitesparse