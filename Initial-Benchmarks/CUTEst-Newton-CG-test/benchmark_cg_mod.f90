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
    
    ! rnorm0 = norm2(g) capturado uma vez, antes do loop
    if ( rnorm .le. eps * max(1.0d0, gnorm) ) then
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

! *****************************************************************
! *****************************************************************

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

! *****************************************************************
! *****************************************************************

end module benchmark_cg_mod