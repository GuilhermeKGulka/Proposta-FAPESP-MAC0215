program benchmark_cg
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

  ! Arrays CSR (Host)
  integer(c_int), allocatable, target :: hrow_ptr(:), hcol(:)
  real(8), allocatable, target :: hval(:), d(:), g(:), d_gpu(:), d_cpu(:)

  ! Contexto GPU
  type(rocm_cg_context) :: ctx

  ! Variaveis de controle do CG
  real(8), parameter :: eps = 1.0d-08
  real(8) :: gnorm
  integer :: cg_iter_gpu, cg_iter_cpu, cgmaxit, istop_gpu, istop_cpu
  logical :: do_prep
  character(len=16) :: arg2

  ! Variaveis de benchmark
  integer(8) :: t_start, t_end, rate
  real(8) :: t_gpu_copy, t_gpu_cg, t_cpu_cg, t_gpu_prep
  real(8) :: speedup, speedup_with_copy, density, iter_diff_rel, res_gpu, res_cpu

  call system_clock(count_rate=rate)

  ! LER ARGUMENTO DA LINHA DE COMANDO
  if (command_argument_count() < 1) then
    write(*,*) "Uso: ./benchmark_cg <arquivo.mtx> [nopre]"
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

  ! VALIDACAO DA ENTRADA: o codigo assume matriz quadrada com indices em [1, n].
  ! Fora disso o counting sort escreveria alem dos vetores (CSR corrompido).
  if (nnz_total > 0) then
    if (m /= n .or. minval(coo_row(1:nnz_total)) < 1 .or. maxval(coo_row(1:nnz_total)) > n &
        .or. minval(coo_col(1:nnz_total)) < 1 .or. maxval(coo_col(1:nnz_total)) > n) then
      write(*,'(A,A,A,I0,A,I0,A)') "IGNORADA: ", trim(problem_name), &
            ": nao quadrada ou indices fora de [1,n] (m=", m, ", n=", n, ")"
      stop 2
    end if
  end if

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

  ! CONFIGURACOES INICIAIS DO METODO
  allocate(d(n), g(n), d_gpu(n), d_cpu(n))
  g(1:n) = -1.0d0
  d(1:n) = 0.0d0
  gnorm = norm2(g(1:n))
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

  ! DESCRITORES GENERICOS, BUFFERS E PREPROCESS (setup unico, cronometrado a parte)
  call setup_spmv(ctx, n, nnz_total, do_prep, t_gpu_prep)

  ! BENCHMARK GPU CG
  call system_clock(t_start)
  call gpu_conjugated_gradient(cg_iter_gpu, istop_gpu, eps, gnorm, cgmaxit, n, nnz_total, ctx)
  call hipCheck(hipDeviceSynchronize())
  call system_clock(t_end)
  t_gpu_cg = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0

  ! Solucao da GPU de volta ao host (fora dos timers)
  call hipCheck(hipMemcpy(c_loc(d_gpu), ctx%d_d, n*8_c_size_t, hipMemcpyDeviceToHost))

  ! BENCHMARK CPU CG
  d_cpu(1:n) = d(1:n)

  call system_clock(t_start)
  call cpu_conjugated_gradient(cg_iter_cpu, istop_cpu, eps, g, cgmaxit, n, nnz_total, d_cpu, hval, hcol, hrow_ptr)
  call system_clock(t_end)
  t_cpu_cg = real(t_end - t_start, 8) / real(rate, 8) * 1000.0d0

  ! RESIDUO VERDADEIRO (||g + H d|| / max(1,||g||)) das duas solucoes, fora dos timers
  res_gpu = rel_residual(d_gpu)
  res_cpu = rel_residual(d_cpu)

  ! METRICAS E ESCRITA
  speedup = t_cpu_cg / t_gpu_cg
  speedup_with_copy = t_cpu_cg / (t_gpu_copy + t_gpu_cg)
  density = real(nnz_total, 8) / (real(n, 8) * real(n, 8))
  iter_diff_rel = real(abs(cg_iter_gpu - cg_iter_cpu),8) / (real(max(1, cg_iter_cpu),8))

  open(unit=20, file='results_cg.csv', status='unknown', position='append', action='write')
  write(20, '(A, ",", I0, ",", I0, ",", ES12.5, ",", I0, ",", I0, ",", I0, ",", I0, ",", ES12.5, ",", &
              &ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", ES12.5, ",", &
              &ES12.5, ",", ES12.5)') &
        trim(problem_name), n, nnz_total, density, cg_iter_gpu, cg_iter_cpu, istop_gpu, istop_cpu, &
        iter_diff_rel, t_gpu_copy, t_gpu_cg, t_cpu_cg, speedup, speedup_with_copy, &
        t_gpu_prep, res_gpu, res_cpu
  close(20)

  call finalize_rocm_context(ctx)

  deallocate(hrow_ptr, hcol, hval, d, g, d_gpu, d_cpu)

contains

  ! Residuo relativo verdadeiro de H*x = -g (g e o termo independente do CG: H d = -g)
  function rel_residual(x) result(res)
    real(8), intent(in) :: x(:)
    real(8) :: res, s
    real(8), allocatable :: rv(:)
    integer :: row, kk

    allocate(rv(n))
    do row = 1, n
      s = 0.0d0
      do kk = hrow_ptr(row), hrow_ptr(row+1) - 1
        s = s + hval(kk) * x(hcol(kk))
      end do
      rv(row) = g(row) + s
    end do
    res = norm2(rv) / max(1.0d0, norm2(g))
    deallocate(rv)
  end function rel_residual

end program benchmark_cg
