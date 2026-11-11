module gencan_gpu_context
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipfort_hipBLAS
  use hipsparse_generic
  implicit none

  private
  public :: gencan_cg_context, init_gencan_gpu, finalize_gencan_gpu, copy_hessian_to_gpu
  
  real(c_double), target, save :: alpha_one = 1.0d0, alpha_zero = 0.0d0

  type :: gencan_cg_context
    type(c_ptr) :: hsparse_handle = c_null_ptr
    type(c_ptr) :: hblas_handle   = c_null_ptr
    type(c_ptr) :: mat_h          = c_null_ptr
    type(c_ptr) :: d_buf          = c_null_ptr
    
    type(c_ptr) :: d_x   = c_null_ptr
    type(c_ptr) :: d_g   = c_null_ptr
    type(c_ptr) :: d_d   = c_null_ptr
    type(c_ptr) :: d_p   = c_null_ptr
    type(c_ptr) :: d_hp  = c_null_ptr
    type(c_ptr) :: d_r   = c_null_ptr
    type(c_ptr) :: d_z   = c_null_ptr
    
    type(c_ptr) :: d_x_lb = c_null_ptr
    type(c_ptr) :: d_x_ub = c_null_ptr
    type(c_ptr) :: d_lind = c_null_ptr
    type(c_ptr) :: d_uind = c_null_ptr
    type(c_ptr) :: d_ifr  = c_null_ptr

    type(c_ptr) :: d_hrow_ptr = c_null_ptr
    type(c_ptr) :: d_hcol     = c_null_ptr
    type(c_ptr) :: d_hval     = c_null_ptr
    
    type(c_ptr) :: d_alpha = c_null_ptr
    type(c_ptr) :: d_beta  = c_null_ptr
    type(c_ptr) :: d_dot   = c_null_ptr
    type(c_ptr) :: d_amax  = c_null_ptr
    type(c_ptr) :: d_abox  = c_null_ptr
    
    integer :: n = 0
    integer :: nnz = 0
    logical :: initialized = .false.
    logical :: spmv_ready = .false.
  end type

contains

  subroutine init_gencan_gpu(ctx, n, nnz_max)
    type(gencan_cg_context), intent(inout) :: ctx
    integer, intent(in) :: n, nnz_max

    if (ctx%initialized) return
    ctx%n = n

    call hipCheck(hipMalloc(ctx%d_x,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_g,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_x_lb, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_x_ub, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_lind, n * 4_c_size_t)) 
    call hipCheck(hipMalloc(ctx%d_uind, n * 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_ifr,  n * 4_c_size_t))

    call hipCheck(hipMalloc(ctx%d_d,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_p,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_hp, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_r,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_z,  n * 8_c_size_t))

    call hipCheck(hipMalloc(ctx%d_alpha, 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_beta,  8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_dot,   8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_amax,  8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_abox,  8_c_size_t))

    ! CORRIGIDO: Alocando o dobro do espaço para a Hessiana espelhada (simétrica completa)
    call hipCheck(hipMalloc(ctx%d_hrow_ptr, (n+1) * 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_hcol, 2 * nnz_max * 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_hval, 2 * nnz_max * 8_c_size_t))

    call hipsparseCheck(hipsparseCreate(ctx%hsparse_handle))
    call hipblasCheck(hipblasCreate(ctx%hblas_handle))

    ctx%initialized = .true.
  end subroutine

  subroutine copy_hessian_to_gpu(ctx, n, hnnz_upper, hval_coo, hrow_coo, hcol_coo)
    type(gencan_cg_context), intent(inout) :: ctx
    integer, intent(in) :: n, hnnz_upper
    real(8), intent(in) :: hval_coo(hnnz_upper)
    integer, intent(in) :: hrow_coo(hnnz_upper), hcol_coo(hnnz_upper)

    integer :: i, k, nnz_total, r, c, idx
    integer(c_size_t) :: bsz
    integer, allocatable, target :: row_counts(:), cur_ptr(:), hcol_csr(:), hrow_ptr(:)
    real(8), allocatable, target :: hval_csr(:)

    ! 1. Contagem
    nnz_total = hnnz_upper
    do k = 1, hnnz_upper
      if (hrow_coo(k) /= hcol_coo(k)) nnz_total = nnz_total + 1
    end do
    ctx%nnz = nnz_total

    allocate(row_counts(n), cur_ptr(n), hrow_ptr(n+1), hcol_csr(nnz_total), hval_csr(nnz_total))
    row_counts = 0
    do k = 1, hnnz_upper
      row_counts(hrow_coo(k)) = row_counts(hrow_coo(k)) + 1
      if (hrow_coo(k) /= hcol_coo(k)) &
        row_counts(hcol_coo(k)) = row_counts(hcol_coo(k)) + 1
    end do

    ! 2. Ponteiros
    hrow_ptr(1) = 1
    do i = 1, n
      hrow_ptr(i+1) = hrow_ptr(i) + row_counts(i)
    end do

    ! 3. Preenchimento
    cur_ptr = hrow_ptr(1:n)
    do k = 1, hnnz_upper
      r = hrow_coo(k); c = hcol_coo(k)
      idx = cur_ptr(r); hcol_csr(idx) = c; hval_csr(idx) = hval_coo(k)
      cur_ptr(r) = cur_ptr(r) + 1
      if (r /= c) then
        idx = cur_ptr(c); hcol_csr(idx) = r; hval_csr(idx) = hval_coo(k)
        cur_ptr(c) = cur_ptr(c) + 1
      end if
    end do

    ! 4. Enviar para GPU
    call hipCheck(hipMemcpy(ctx%d_hrow_ptr, c_loc(hrow_ptr), (n+1)*4_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(ctx%d_hcol, c_loc(hcol_csr), nnz_total*4_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(ctx%d_hval, c_loc(hval_csr), nnz_total*8_c_size_t, hipMemcpyHostToDevice))

    ! 5. Recriar Descritor hipSPARSE
    if (ctx%spmv_ready) then
      call hs_check(hs_destroy_spmat(ctx%mat_h), "DestroySpMat")
      call hs_free_buffer(ctx%d_buf)
    end if
    
    call hs_check(hs_create_csr(ctx%mat_h, int(n, c_int64_t), int(n, c_int64_t), &
                  int(nnz_total, c_int64_t), ctx%d_hrow_ptr, ctx%d_hcol, ctx%d_hval, &
                  HS_INDEX_32I, HS_INDEX_32I, HS_BASE_ONE, HS_R_64F), "CreateCsr")
    
    call hs_check(hs_spmv_buffer_size(ctx%hsparse_handle, HS_OP_N, c_loc(alpha_one), ctx%mat_h, &
                  ctx%d_p, c_loc(alpha_zero), ctx%d_hp, HS_R_64F, HS_SPMV_ALG_DEFAULT, bsz), "BufferSize")
    
    call hs_alloc_buffer(ctx%d_buf, bsz)
    
    ctx%spmv_ready = .true.

    deallocate(row_counts, cur_ptr, hrow_ptr, hcol_csr, hval_csr)
  end subroutine

  subroutine finalize_gencan_gpu(ctx)
    type(gencan_cg_context), intent(inout) :: ctx
    if (.not. ctx%initialized) return
    
    call hipCheck(hipFree(ctx%d_x)); call hipCheck(hipFree(ctx%d_g)); call hipCheck(hipFree(ctx%d_d))
    call hipCheck(hipFree(ctx%d_p)); call hipCheck(hipFree(ctx%d_hp)); call hipCheck(hipFree(ctx%d_r))
    call hipCheck(hipFree(ctx%d_z)); call hipCheck(hipFree(ctx%d_x_lb)); call hipCheck(hipFree(ctx%d_x_ub))
    call hipCheck(hipFree(ctx%d_lind)); call hipCheck(hipFree(ctx%d_uind)); call hipCheck(hipFree(ctx%d_ifr))
    call hipCheck(hipFree(ctx%d_hrow_ptr)); call hipCheck(hipFree(ctx%d_hcol)); call hipCheck(hipFree(ctx%d_hval))
    call hipCheck(hipFree(ctx%d_alpha)); call hipCheck(hipFree(ctx%d_beta)); call hipCheck(hipFree(ctx%d_dot))
    call hipCheck(hipFree(ctx%d_amax)); call hipCheck(hipFree(ctx%d_abox))
    
    if (ctx%spmv_ready) then
      call hs_check(hs_destroy_spmat(ctx%mat_h), "DestroySpMat")
      call hs_free_buffer(ctx%d_buf)
    end if
    
    call hipsparseCheck(hipsparseDestroy(ctx%hsparse_handle))
    call hipblasCheck(hipblasDestroy(ctx%hblas_handle))
    
    ctx%initialized = .false.
    ctx%spmv_ready = .false.
  end subroutine

end module gencan_gpu_context