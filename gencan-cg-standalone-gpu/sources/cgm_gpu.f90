module modcgm_gpu
    use lssit
    use gencan_gpu_context
    use hipfort
    use hipfort_check
    use hipfort_hipsparse
    use hipfort_hipBLAS
    use hipsparse_generic
    use iso_c_binding
    implicit none

    private
    public :: cgm_gpu

    real(c_double), target, save :: alpha_one = 1.0d0, alpha_zero = 0.0d0

    interface
       subroutine launch_compute_amax_abox(nfr, delta, p, d, x, lbnd, ubnd, lind, uind, ifr, amax, abox) bind(C, name="launch_compute_amax_abox")
         import :: c_ptr, c_int
         integer(c_int), value :: nfr
         type(c_ptr), value :: delta, p, d, x, lbnd, ubnd, lind, uind, ifr, amax, abox
       end subroutine
       subroutine launch_update_d_r(nfr, alpha, p, hp, d, r) bind(C, name="launch_update_d_r")
         import :: c_ptr, c_int
         integer(c_int), value :: nfr
         type(c_ptr), value :: alpha, p, hp, d, r
       end subroutine
       subroutine launch_update_p_vector(nfr, p, z, bbeta) bind(C, name="launch_update_p_vector")
         import :: c_ptr, c_int, c_double
         integer(c_int), value :: nfr
         type(c_ptr), value :: p, z
         real(c_double), value :: bbeta
       end subroutine
    end interface

  contains
    subroutine cgm_gpu (n,nfr,ifr,x,lind,lbnd,uind,ubnd,g,delta,eps,maxit,&
            precond,d,iter,iprint,ierr,cginfo,lssit_obj, ctx)
        
        implicit none

        integer, intent(in) :: n,nfr,maxit,iprint
        integer, intent(out):: iter,ierr,cginfo
        real(kind=8), intent(in) :: delta,eps
        type(lssit_type), intent(in) :: lssit_obj
        type(gencan_cg_context), intent(inout) :: ctx
        logical, intent(in) :: precond

        logical, intent(in) :: lind(n),uind(n)
        
        ! CORRIGIDO: Adicionado 'target' no ifr
        integer, intent(in), target :: ifr(nfr)
        
        real(kind=8), intent(in), target :: x(n),g(nfr),lbnd(n),ubnd(n)
        real(kind=8), intent(out) :: d(n)
        
        real(kind=8), target :: alpha, amax, abox, ztr, rnorm2, gnorm2, ptr, pthp, dnorm2, gtd
        real(kind=8), target :: q, qprev, bestprog, currprog, ztrprev, znorm2, bbeta
        real(kind=8), parameter :: macheps = 1.0d-16, epsnqmp = 1.0d-08
        integer, parameter :: maxcgitnp = 5
        real(kind=8), target :: delta_target
        integer :: itnqmp
        logical :: samep, restarted

        integer(c_int), allocatable, target :: lind_int(:), uind_int(:)
        real(kind=8), allocatable, target :: d_gpu(:)
        
        ierr = 0; iter = 0; cginfo = 0; itnqmp = 0
        restarted = .false.
        qprev = huge(1.0d0); bestprog = 0.0d0
        delta_target = delta

        call copy_hessian_to_gpu(ctx, n, lssit_obj%hne, lssit_obj%hval, lssit_obj%hrow, lssit_obj%hcol)

        allocate(lind_int(n), uind_int(n))
        where (lind)
           lind_int = 1
        elsewhere
           lind_int = 0
        end where
        
        where (uind)
           uind_int = 1
        elsewhere
           uind_int = 0
        end where
        
        call hipCheck(hipMemcpy(ctx%d_x, c_loc(x), n*8_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(ctx%d_g, c_loc(g), nfr*8_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(ctx%d_x_lb, c_loc(lbnd), n*8_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(ctx%d_x_ub, c_loc(ubnd), n*8_c_size_t, hipMemcpyHostToDevice))
        
        call hipCheck(hipMemcpy(ctx%d_lind, c_loc(lind_int), n*4_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(ctx%d_uind, c_loc(uind_int), n*4_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemcpy(ctx%d_ifr, c_loc(ifr), nfr*4_c_size_t, hipMemcpyHostToDevice))
        
        call hipCheck(hipMemcpy(ctx%d_r, c_loc(g), nfr*8_c_size_t, hipMemcpyHostToDevice))
        call hipCheck(hipMemset(ctx%d_d, 0, nfr*8_c_size_t))
        
        call hipblasCheck(hipblasDotEx(ctx%hblas_handle, nfr, ctx%d_r, HIP_R_64F, 1, &
                                       ctx%d_r, HIP_R_64F, 1, ctx%d_dot, HIP_R_64F, HIP_R_64F))
        call hipCheck(hipMemcpy(c_loc(rnorm2), ctx%d_dot, 8_c_size_t, hipMemcpyDeviceToHost))
        
        call hipblasCheck(hipblasDotEx(ctx%hblas_handle, nfr, ctx%d_g, HIP_R_64F, 1, &
                                       ctx%d_g, HIP_R_64F, 1, ctx%d_dot, HIP_R_64F, HIP_R_64F))
        call hipCheck(hipMemcpy(c_loc(gnorm2), ctx%d_dot, 8_c_size_t, hipMemcpyDeviceToHost))
        
        ztr = rnorm2
        ztrprev = ztr

100     continue
        if ( iter .ne. 0 .and. &
            ( ( rnorm2 .le. eps * eps * gnorm2 .and. iter .ge. 4 ) .or. &
            ( rnorm2 .le. macheps ) ) ) then
            cginfo = 0; return
        end if

        if ( iter .ge. max(4, maxit) ) then
            cginfo = 8; return
        end if

        call hipCheck(hipMemcpy(ctx%d_z, ctx%d_r, nfr*8_c_size_t, hipMemcpyDeviceToDevice))
        
        if ( iter .eq. 0 ) then
            call hipCheck(hipMemcpy(ctx%d_p, ctx%d_r, nfr*8_c_size_t, hipMemcpyDeviceToDevice))
            call launch_update_p_vector(int(nfr, c_int), ctx%d_p, ctx%d_z, 0.0_c_double)
        else
            bbeta = ztr / ztrprev
            call launch_update_p_vector(int(nfr, c_int), ctx%d_p, ctx%d_z, real(bbeta, c_double))
        end if

        call hipblasCheck(hipblasDotEx(ctx%hblas_handle, nfr, ctx%d_p, HIP_R_64F, 1, &
                                       ctx%d_r, HIP_R_64F, 1, ctx%d_dot, HIP_R_64F, HIP_R_64F))
        call hipCheck(hipMemcpy(c_loc(ptr), ctx%d_dot, 8_c_size_t, hipMemcpyDeviceToHost))
        
        if ( ptr .gt. 0.0d0 ) then
            call launch_update_p_vector(int(nfr, c_int), ctx%d_p, ctx%d_z, 0.0_c_double)
            ptr = -ptr
        end if

        call hs_check(hs_spmv(ctx%hsparse_handle, HS_OP_N, c_loc(alpha_one), ctx%mat_h, ctx%d_p, &
                      c_loc(alpha_zero), ctx%d_hp, HS_R_64F, HS_SPMV_ALG_DEFAULT, ctx%d_buf), "SpMV")
        
        call hipblasCheck(hipblasDotEx(ctx%hblas_handle, nfr, ctx%d_p, HIP_R_64F, 1, &
                                       ctx%d_hp, HIP_R_64F, 1, ctx%d_dot, HIP_R_64F, HIP_R_64F))
        call hipCheck(hipMemcpy(c_loc(pthp), ctx%d_dot, 8_c_size_t, hipMemcpyDeviceToHost))

        call launch_compute_amax_abox(int(nfr, c_int), c_loc(delta_target), ctx%d_p, ctx%d_d, &
                                       ctx%d_x, ctx%d_x_lb, ctx%d_x_ub, ctx%d_lind, ctx%d_uind, &
                                       ctx%d_ifr, ctx%d_amax, ctx%d_abox)
        call hipCheck(hipMemcpy(c_loc(amax), ctx%d_amax, 8_c_size_t, hipMemcpyDeviceToHost))
        call hipCheck(hipMemcpy(c_loc(abox), ctx%d_abox, 8_c_size_t, hipMemcpyDeviceToHost))

        if ( pthp .gt. 0.0d0 ) then
            alpha = min( amax, ztr / pthp )
        else if ( iter .eq. 0 ) then
            alpha = amax
        else
            cginfo = 7; return
        end if

        call launch_update_d_r(int(nfr, c_int), c_loc(alpha), ctx%d_p, ctx%d_hp, ctx%d_d, ctx%d_r)
        
        call hipblasCheck(hipblasDotEx(ctx%hblas_handle, nfr, ctx%d_r, HIP_R_64F, 1, &
                                       ctx%d_r, HIP_R_64F, 1, ctx%d_dot, HIP_R_64F, HIP_R_64F))
        call hipCheck(hipMemcpy(c_loc(rnorm2), ctx%d_dot, 8_c_size_t, hipMemcpyDeviceToHost))

        if ( alpha .eq. amax ) then; cginfo = 1; return; end if
        if ( abox .le. 0.1d0 ) then; cginfo = 5; return; end if
        
        samep = .true.
        if ( abs(alpha * maxval(abs(d))) .gt. macheps * max(1.0d0, maxval(abs(d))) ) samep = .false.
        if ( samep ) then; cginfo = 6; return; end if

        ztrprev = ztr
        ztr = rnorm2
        iter = iter + 1
        goto 100

999     continue
        allocate(d_gpu(n))
        call hipCheck(hipMemcpy(c_loc(d_gpu), ctx%d_d, nfr*8_c_size_t, hipMemcpyDeviceToHost))
        d(1:nfr) = d_gpu(1:nfr)
        deallocate(d_gpu, lind_int, uind_int)
        
        ierr = 0
    end subroutine cgm_gpu
end module