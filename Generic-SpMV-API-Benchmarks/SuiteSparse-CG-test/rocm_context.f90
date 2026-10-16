module rocm_context
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipfort_hipBLAS
  use hipsparse_generic
  implicit none

  private
  public :: rocm_cg_context, init_rocm_context, setup_spmv, finalize_rocm_context
  public :: ctx_spmv_residual, ctx_spmv_hp

  ! Escalares do SpMV (ponteiros de host: modo de ponteiro padrao do hipSPARSE)
  real(c_double), target, save :: hs_alpha = 1.0d0, hs_beta0 = 0.0d0, hs_beta1 = 1.0d0

  type :: rocm_cg_context
    ! Handles
    type(c_ptr) :: hsparse_handle, hblas_handle

    ! Descritores genericos do hipSPARSE (matriz CSR e vetores)
    type(c_ptr) :: mat_h  = c_null_ptr
    type(c_ptr) :: vec_d  = c_null_ptr
    type(c_ptr) :: vec_r  = c_null_ptr
    type(c_ptr) :: vec_p  = c_null_ptr
    type(c_ptr) :: vec_Hp = c_null_ptr

    ! Buffer externo do SpMV: UM por matriz, usado nas duas chamadas ((r = H d + r) e
    ! (Hp = H p)). Nulo quando bufferSize devolve 0 (regra do hipSPARSE).
    type(c_ptr) :: d_buf = c_null_ptr

    ! Device pointers
    type(c_ptr) :: d_hrow_ptr = c_null_ptr
    type(c_ptr) :: d_hcol     = c_null_ptr
    type(c_ptr) :: d_hval     = c_null_ptr
    type(c_ptr) :: d_d        = c_null_ptr
    type(c_ptr) :: d_r        = c_null_ptr
    type(c_ptr) :: d_p        = c_null_ptr
    type(c_ptr) :: d_Hp       = c_null_ptr

    type(c_ptr) :: d_alpha    = c_null_ptr
    type(c_ptr) :: d_beta     = c_null_ptr
    type(c_ptr) :: d_dot_pHp  = c_null_ptr
    type(c_ptr) :: d_rnorm    = c_null_ptr
    type(c_ptr) :: d_r2       = c_null_ptr
    type(c_ptr) :: d_lastr2   = c_null_ptr

    ! Status
    logical :: initialized = .false.
    logical :: spmv_ready  = .false.
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

    ctx%initialized = .true.
  end subroutine init_rocm_context

  ! Cria descritores genericos, buffers e (opcionalmente) faz o preprocess.
  ! Deve ser chamada DEPOIS de copiar row_ptr/col/val para a GPU.
  ! t_prep_ms: tempo do preprocess (0 se do_prep = .false.; -1 se falhou).
  subroutine setup_spmv(ctx, n, nnz, do_prep, t_prep_ms)
    implicit none

    type(rocm_cg_context), intent(inout) :: ctx
    integer, intent(in) :: n, nnz
    logical, intent(in) :: do_prep
    real(8), intent(out) :: t_prep_ms

    integer(c_size_t) :: bsz_d, bsz_p, bsz
    integer(c_int) :: st
    integer(8) :: t0, t1, rate

    call hs_check(hs_create_csr(ctx%mat_h, int(n, c_int64_t), int(n, c_int64_t), &
                  int(nnz, c_int64_t), ctx%d_hrow_ptr, ctx%d_hcol, ctx%d_hval, &
                  HS_INDEX_32I, HS_INDEX_32I, HS_BASE_ONE, HS_R_64F), "CreateCsr")
    call hs_check(hs_create_dnvec(ctx%vec_d,  int(n, c_int64_t), ctx%d_d,  HS_R_64F), "CreateDnVec d")
    call hs_check(hs_create_dnvec(ctx%vec_r,  int(n, c_int64_t), ctx%d_r,  HS_R_64F), "CreateDnVec r")
    call hs_check(hs_create_dnvec(ctx%vec_p,  int(n, c_int64_t), ctx%d_p,  HS_R_64F), "CreateDnVec p")
    call hs_check(hs_create_dnvec(ctx%vec_Hp, int(n, c_int64_t), ctx%d_Hp, HS_R_64F), "CreateDnVec Hp")

    ! A analise do preprocess pertence a MATRIZ, nao aos vetores: usa-se um unico buffer
    ! (tamanho = maximo dos dois) e um unico preprocess, reaproveitados nas duas chamadas.
    ! chamada 1: r = 1*H*d + 1*r   (x = d, y = r, beta = 1)
    call hs_check(hs_spmv_buffer_size(ctx%hsparse_handle, HS_OP_N, c_loc(hs_alpha), ctx%mat_h, &
                  ctx%vec_d, c_loc(hs_beta1), ctx%vec_r, HS_R_64F, HS_SPMV_ALG_DEFAULT, bsz_d), &
                  "SpMV_bufferSize (d)")
    ! chamada 2: Hp = 1*H*p + 0*Hp  (x = p, y = Hp, beta = 0)
    call hs_check(hs_spmv_buffer_size(ctx%hsparse_handle, HS_OP_N, c_loc(hs_alpha), ctx%mat_h, &
                  ctx%vec_p, c_loc(hs_beta0), ctx%vec_Hp, HS_R_64F, HS_SPMV_ALG_DEFAULT, bsz_p), &
                  "SpMV_bufferSize (p)")
    bsz = max(bsz_d, bsz_p)
    call hs_alloc_buffer(ctx%d_buf, bsz)      ! nulo se bsz == 0

    t_prep_ms = 0.0d0
    if (do_prep) then
      call system_clock(count_rate=rate)
      call system_clock(t0)
      st = hs_spmv_preprocess(ctx%hsparse_handle, HS_OP_N, c_loc(hs_alpha), ctx%mat_h, &
                  ctx%vec_d, c_loc(hs_beta1), ctx%vec_r, HS_R_64F, HS_SPMV_ALG_DEFAULT, ctx%d_buf)
      call hipCheck(hipDeviceSynchronize())
      call system_clock(t1)
      t_prep_ms = real(t1 - t0, 8) / real(rate, 8) * 1000.0d0
      if (st /= 0) then   ! preprocess e opcional: avisa e segue sem ele
        write(*,'(A,I0,A,I0,A,I0,A)') "AVISO: SpMV_preprocess falhou (codigo ", st, &
                                       "; bufferSize d=", bsz_d, " p=", bsz_p, "); seguindo sem ele"
        t_prep_ms = -1.0d0
      end if
    end if

    ctx%spmv_ready = .true.
  end subroutine setup_spmv

  ! r = H*d + r
  subroutine ctx_spmv_residual(ctx)
    implicit none
    type(rocm_cg_context), intent(in) :: ctx
    call hs_check(hs_spmv(ctx%hsparse_handle, HS_OP_N, c_loc(hs_alpha), ctx%mat_h, ctx%vec_d, &
                  c_loc(hs_beta1), ctx%vec_r, HS_R_64F, HS_SPMV_ALG_DEFAULT, ctx%d_buf), &
                  "SpMV (r = Hd + r)")
  end subroutine ctx_spmv_residual

  ! Hp = H*p
  subroutine ctx_spmv_hp(ctx)
    implicit none
    type(rocm_cg_context), intent(in) :: ctx
    call hs_check(hs_spmv(ctx%hsparse_handle, HS_OP_N, c_loc(hs_alpha), ctx%mat_h, ctx%vec_p, &
                  c_loc(hs_beta0), ctx%vec_Hp, HS_R_64F, HS_SPMV_ALG_DEFAULT, ctx%d_buf), &
                  "SpMV (Hp = Hp)")
  end subroutine ctx_spmv_hp

  subroutine finalize_rocm_context(ctx)
    implicit none

    type(rocm_cg_context), intent(inout) :: ctx

    if (.not. ctx%initialized) return

    ! Destruir descritores genericos e liberar buffers
    if (ctx%spmv_ready) then
      call hs_check(hs_destroy_spmat(ctx%mat_h), "DestroySpMat")
      call hs_check(hs_destroy_dnvec(ctx%vec_d), "DestroyDnVec d")
      call hs_check(hs_destroy_dnvec(ctx%vec_r), "DestroyDnVec r")
      call hs_check(hs_destroy_dnvec(ctx%vec_p), "DestroyDnVec p")
      call hs_check(hs_destroy_dnvec(ctx%vec_Hp), "DestroyDnVec Hp")
      call hs_free_buffer(ctx%d_buf)
      ctx%spmv_ready = .false.
    end if

    ! Liberar memoria
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
  end subroutine finalize_rocm_context

end module rocm_context