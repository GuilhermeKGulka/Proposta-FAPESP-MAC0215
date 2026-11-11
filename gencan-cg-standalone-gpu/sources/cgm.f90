module modcgm  

    use lssit
    use iso_c_binding, only: c_ptr

    implicit none

    private
    real(kind=8), parameter :: macheps = 1.0d-16,macheps12 = sqrt(macheps), &
                                 lrgstreal = huge( 1.0d0 ), epsnqmp = 1.0d-08, theta = 1.0d-06
    integer, parameter :: maxcgitnp  =  5
    
    public cgm  

  contains
    subroutine cgm (n,nfr,ifr,x,lind,lbnd,uind,ubnd,g,delta,eps,maxit,&
            precond,d,iter,iprint,ierr,cginfo,lssit)
        
        implicit none

        !SCALARS ARGUMENTS
        integer, intent(in) :: n,nfr,maxit,iprint
        integer, intent(out):: iter,ierr,cginfo
        real(kind=8), intent(in) :: delta,eps
        type(lssit_type), intent(in) :: lssit

        !ARRAY ARGUMENTS
        logical, intent(in) :: lind(n),uind(n)
        integer, intent(in) :: ifr(nfr)
        real(kind=8), intent(in)  :: x(n),g(nfr),lbnd(n),ubnd(n)
        real(kind=8), intent(out) :: d(n)
    
        !LOCAL SCALARS
        logical          :: goth,gotp,negcur,samep,precond,restarted
        integer          :: i,itnqmp,itertmp
        real(kind=8)     :: alpha,abox,amax,bestprog,bbeta,currprog,   &
                            dnorm2,gnorm2,gtd,gtp,pnorm2,ptd,pthp,ptr, &
                            q,qprev,rnorm2,ztrprev,ztr,znorm2,xsupn

        !LOCAL ARRAYS
        real(kind=8) :: p(nfr),hp(nfr),r(nfr),z(nfr),xaux(n),rfull(n),zfull(n)
        
        ! ==================================================================
        ! Initialization
        ! ==================================================================
        if ( iprint .ge. 1 ) then 
            if ( precond ) then 
                write(*,*) 'Preconditioner CG start:'
            else
                write(*,*) 'CG start:'
            end if
        end if

        ierr = 0
       
        itertmp = 0
        restarted = .false.

001 continue
        goth   = .false.
        gotp   = .false.

        xaux(1:n)  = x(1:n)
        xsupn      = maxval( abs ( x(ifr(1:nfr)) ) )
        gnorm2     = dot_product( g(1:nfr), g(1:nfr) )

        iter       =      0
        itnqmp     =      0
        qprev      = lrgstreal
        bestprog   =  0.0d0

        d(1:nfr)   = 0.0d0 
        r(1:nfr)   = g(1:nfr)
        rfull(1:n) = 0.0d0
        p(1:nfr)   = 0.0d0
        zfull(1:n) = 0.0d0
        q          =  0.0d0
        gtd        =  0.0d0
        dnorm2     =  0.0d0
        rnorm2     = gnorm2
        
        ztr        =  0.0d0

        if ( iprint .ge. 2 ) then
            write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q
        end if
    100 continue 
        ! ==================================================================
        ! Test stopping criteria
        ! ==================================================================

        ! if ||r||_2 = ||H d + g||_2 <= eps * ||g||_2 then stop

        if ( iter .ne. 0 .and. &
            ( ( rnorm2 .le. eps * eps * gnorm2 .and. iter .ge. 4 ) .or. &
            ( rnorm2 .le. macheps ) ) ) then
            cginfo = 0
            if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo
            iter = iter + itertmp
            return
        end if

        ! if the maximum number of iterations was achieved then stop

        if ( iter .ge. max(4, maxit) ) then
            cginfo = 8
            if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo
            iter = iter + itertmp
            return
        end if

        ! ==================================================================
        ! Preconditioner
        ! ==================================================================


        if ( .not. precond .or. restarted ) then
            z(1:nfr) = r(1:nfr)
            ztrprev = ztr
            ztr     = rnorm2
            znorm2  = rnorm2
        else
            rfull(ifr(1:nfr)) = r(1:nfr)
            call applyhpre(n,rfull(1:n),zfull(1:n),lssit)
            z(1:nfr) = zfull(ifr(1:nfr))

            ztrprev = ztr 
            ztr = dot_product( z(1:nfr), r(1:nfr) )
            znorm2 = dot_product( z(1:nfr), z(1:nfr) )
        end if

        ! ==================================================================
        ! Compute direction
        ! ==================================================================

        if ( iter .eq. 0 ) then
            p(1:nfr) = - z(1:nfr)
            ptr    = - ztr
            pnorm2 =   znorm2
        else
            bbeta = ztr / ztrprev
            p(1:nfr) = -z(1:nfr) + bbeta * p(1:nfr)

            if ( .not. precond .and. restarted) then 
                pnorm2 = rnorm2 - 2.0d0 * bbeta * ( ptr + alpha * pthp ) &
                    + bbeta * bbeta * pnorm2
                ptr = - rnorm2 + bbeta * ( ptr + alpha * pthp )
            else
                ptr = dot_product( p(1:nfr), r(1:nfr) )
                pnorm2 = dot_product( p(1:nfr), p(1:nfr) )
            end if 
        end if

        ! Force p to be a descent direction of q(d), i.e.,
        ! <\nabla q(d), p> = <H d + g, p> = <r, p> \le 0.

        if ( ptr .gt. 0.0d0 ) then
            p(1:nfr) = - p(1:nfr)
            ptr = - ptr
        end if
        ! ==================================================================
        ! Compute p^T H p
        ! ==================================================================
        call computehv(n,nfr,ifr(1:nfr),x(1:n),g(1:nfr),p(1:nfr),hp(1:nfr),ierr,lssit)
        if ( ierr .ne. 0 ) return 
        ! Compute p^T hp
        pthp = dot_product(p(1:nfr), hp(1:nfr)) 

        ! ==================================================================
        ! Compute maximum step
        ! ==================================================================
        amax = lrgstreal
        do i = 1,nfr
            if ( p(i) .gt. 0.0d0 ) then
                amax = min( amax,  (   delta - d(i) ) / p(i) )
            else if ( p(i) .lt. 0.0d0 ) then
                amax = min( amax,  ( - delta - d(i) ) / p(i) )
            end if
        end do
 
        ! ==================================================================
        ! Compute the step
        ! ==================================================================
        negcur = .false.
        ! If p^T H p > 0 then take the conjugate gradients step
        if ( pthp .gt. 0.0d0 ) then
            alpha = min( amax, ztr / pthp )
        ! Else, if we are at iteration zero then take the maximum
        ! positive step in the minus gradient direction
        else if ( iter .eq. 0 ) then
            alpha = amax
            negcur = .true.
        ! Otherwise, stop at the current iterate
        else
            cginfo = 7
            if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo
            iter = iter + itertmp
            return
        end if
        

        ! ==================================================================
        ! Test the angle condition
        ! ==================================================================

        ptd = dot_product(p(1:nfr), d(1:nfr))
        gtp = dot_product(g(1:nfr), p(1:nfr))

        ! These are gtd and dnorm2 for the new direction d which was not
        ! computed yet.

        gtd = gtd + alpha * gtp
        dnorm2 = dnorm2 + alpha * alpha * pnorm2 + 2.0d0 * alpha * ptd


        if ( gtd .gt. 0.0d0 .or. &
        ( gtd * gtd .lt. 1.0d-12 * gnorm2 * dnorm2 ) ) then

            if ( precond .and. .not. restarted .and. iter .eq. 0 ) then
                if ( iprint .ge. 1 ) write(*,*) 'CG will be restarted without preconditioner)'
                restarted  = .true.
                itertmp    = iter
                go to 001
            end if
     
          cginfo = 3 
          if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo 
          return 
        end if
       
        ! ==================================================================
        ! Compute the quadratic model functional value at the new point
        ! ==================================================================
        qprev = q
        q = q + 0.5d0 * alpha * alpha * pthp + alpha * ptr

        ! ==================================================================
        ! Compute new d
        ! ==================================================================
        d(1:nfr) = d(1:nfr) + alpha * p(1:nfr)

        ! ==================================================================
        ! Compute the residual r = H d + g
        ! ==================================================================
        r(1:nfr) = r(1:nfr) + alpha * hp(1:nfr)
        rnorm2 = dot_product( r(1:nfr), r(1:nfr) )

        ! ==================================================================
        ! Increment number of iterations
        ! ==================================================================

        iter = iter + 1

        ! ==================================================================
        ! Print information of this iteration
        ! ==================================================================

        if ( iprint .ge. 2 ) then
            write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q
        end if

        ! ==================================================================
        ! Test other stopping criteria
        ! ==================================================================

        ! Boundary of the "trust region"

        if ( alpha .eq. amax ) then
            cginfo = 1
            if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo
            iter = iter + itertmp
            return
        end if

        ! Small useful proportion

        abox = lrgstreal
        do i = 1,nfr
            if ( d(i) .gt. 0.0d0 .and. uind(ifr(i)) ) then
                abox = min( abox, ( ubnd(ifr(i)) - x(ifr(i)) ) / d(i) )
            else if ( d(i) .lt. 0.0d0 .and. lind(ifr(i)) ) then
                abox = min( abox, ( lbnd(ifr(i)) - x(ifr(i)) ) / d(i) )
            end if
        end do
        if ( abox .le. 0.1d0 ) then
            cginfo = 5
            if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo
            iter = iter + itertmp
            return
        end if

        ! Two consecutive iterates are too much close

        samep = .true.
        do i = 1,nfr
            if ( abs( alpha * p(i) ) .gt. &
                    macheps * max( 1.0d0, abs( d(i) ) ) ) then
                samep = .false.
            end if
        end do
    
        if ( samep ) then
            cginfo = 6
            if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo
            iter = iter + itertmp
            return
        end if

        ! Many iterations without good progress of the quadratic model

        currprog = qprev - q
        bestprog = max( currprog, bestprog )

        if ( currprog .le. epsnqmp * bestprog ) then
            itnqmp = itnqmp + 1
            if ( itnqmp .ge. maxcgitnp ) then
                cginfo = 4
                if ( iprint .ge. 1 ) write(*,*) 'Exit cgm with cginfo = ', cginfo
                iter = iter + itertmp
                return
            endif
        else
            itnqmp = 0
        endif

        ! ==================================================================
        ! Iterate
        ! ==================================================================

        go to 100

        ! ==================================================================
        ! End of main loop
        ! ==================================================================

984  format(  5X,'CG iter = ',I7,' rnorm = ',1P,D10.4,' dnorm = ', &
    1P,D10.4,' q = ',1P,D11.4)
    end subroutine cgm

  end module modcgm
  