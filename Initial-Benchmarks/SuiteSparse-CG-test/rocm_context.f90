module rocm_context
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipfort_hipBLAS
  implicit none

  private
  public :: rocm_cg_context, init_rocm_context, finalize_rocm_context

  type :: rocm_cg_context
    ! Handles
    type(c_ptr) :: hsparse_handle, hblas_handle

    ! Descritores
    type(c_ptr) :: descrH, descr_vec_d, descr_vec_r, descr_vec_p, descr_vec_Hp
    
    ! Device pointers
    type(c_ptr) :: d_hrow_ptr = c_null_ptr
    type(c_ptr) :: d_hcol     = c_null_ptr
    type(c_ptr) :: d_hval     = c_null_ptr
    type(c_ptr) :: d_d        = c_null_ptr
    type(c_ptr) :: d_r        = c_null_ptr
    type(c_ptr) :: d_p        = c_null_ptr
    type(c_ptr) :: d_Hp       = c_null_ptr
    type(c_ptr) :: d_buffer

    type(c_ptr) :: d_alpha    = c_null_ptr
    type(c_ptr) :: d_beta     = c_null_ptr
    type(c_ptr) :: d_dot_pHp  = c_null_ptr
    type(c_ptr) :: d_rnorm    = c_null_ptr
    type(c_ptr) :: d_r2       = c_null_ptr
    type(c_ptr) :: d_lastr2   = c_null_ptr
    
    ! Status
    logical :: initialized = .false.
  end type rocm_cg_context
  
contains

  subroutine init_rocm_context(ctx, n, hnnz)
    implicit none
    
    type(rocm_cg_context), intent(inout) :: ctx
    integer, intent(in) :: n, hnnz
    
    if (ctx%initialized) return

    ! Alocar matriz
    call hipCheck(hipMalloc(ctx%d_hrow_ptr, (n+1) * 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_hcol, hnnz * 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_hval, hnnz * 8_c_size_t))
    
    ! Alocar vetores
    call hipCheck(hipMalloc(ctx%d_d, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_r, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_p, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_Hp, n * 8_c_size_t))

    ! Alocar escalares
    call hipCheck(hipMalloc(ctx%d_alpha, 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_beta, 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_dot_pHp, 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_rnorm, 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_r2, 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_lastr2, 8_c_size_t))
    
    ! Criar handles
    call hipsparseCheck(hipsparseCreate(ctx%hsparse_handle))
    call hipblasCheck(hipblasCreate(ctx%hblas_handle))
    
    ! Criar descritores da matriz
    call hipsparseCheck(hipsparseCreateMatDescr(ctx%descrH))
    call hipsparseCheck(hipsparseSetMatIndexBase(ctx%descrH, HIPSPARSE_INDEX_BASE_ONE))
    call hipsparseCheck(hipsparseSetMatType(ctx%descrH, HIPSPARSE_MATRIX_TYPE_GENERAL))

    
    ctx%initialized = .true.
  
  end subroutine init_rocm_context

  subroutine finalize_rocm_context(ctx)
    use hipfort
    use hipfort_check
    use hipfort_hipsparse
    implicit none
    
    type(rocm_cg_context), intent(inout) :: ctx
    
    if (.not. ctx%initialized) return
    
    ! Destruir descritores
    call hipsparseCheck(hipsparseDestroyMatDescr(ctx%descrH)) 

    ! Liberar memória
    call hipCheck(hipFree(ctx%d_hval))
    call hipCheck(hipFree(ctx%d_hrow_ptr))
    call hipCheck(hipFree(ctx%d_hcol))
    call hipCheck(hipFree(ctx%d_d))
    call hipCheck(hipFree(ctx%d_r))
    call hipCheck(hipFree(ctx%d_p))
    call hipCheck(hipFree(ctx%d_Hp))

    call hipCheck(hipFree(ctx%d_alpha))
    call hipCheck(hipFree(ctx%d_beta))
    call hipCheck(hipFree(ctx%d_dot_pHp))
    call hipCheck(hipFree(ctx%d_rnorm))
    call hipCheck(hipFree(ctx%d_r2))
    call hipCheck(hipFree(ctx%d_lastr2))

    ! Destruir handles
    call hipsparseCheck(hipsparseDestroy(ctx%hsparse_handle))
    call hipblasCheck(hipblasDestroy(ctx%hblas_handle))

    ctx%initialized = .false.
  end subroutine

end module rocm_context