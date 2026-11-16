module modgencan

    use lssit 
    use modcgm, only: cgm
    use iso_c_binding, only: c_ptr
    
    implicit none
  
    interface
       subroutine fsub(n,x,f,ierr,pdata)
         use iso_c_binding, only: c_ptr
         implicit none
         integer, intent(in) :: n
         real(kind=8), intent(in) :: x(n)
         real(kind=8), intent(out) :: f 
         integer, intent(inout) :: ierr
         type(c_ptr), optional, intent(in) :: pdata
       end subroutine fsub
  
       subroutine gsub(n,x,g,ierr,pdata)
         use iso_c_binding, only: c_ptr
         implicit none
         integer, intent(in) :: n
         real(kind=8), intent(in) :: x(n)
         real(kind=8), intent(out) :: g(n)
         integer, intent(inout) :: ierr
         type(c_ptr), optional, intent(in) :: pdata
       end subroutine gsub
  
       subroutine hsub(n,x,lim,hnnz,hrow,hcol,hval,ierr,pdata)
         use iso_c_binding, only: c_ptr
         implicit none
         integer, intent(in) :: n,lim
         real(kind=8), intent(in) :: x(n)
         integer, intent(out) :: hnnz
         integer, intent(out) :: hrow(lim),hcol(lim)
         real(kind=8), intent(out) :: hval(lim)
         integer, intent(inout) :: ierr
         type(c_ptr), optional, intent(in) :: pdata
       end subroutine hsub
  
       subroutine stpsub(n,x,gsupn,iter,inhdefstp,stp,ierr,pdata)
         use iso_c_binding, only: c_ptr
         implicit none
         integer, intent(in) :: n,iter
         real(kind=8), intent(in) :: gsupn
         real(kind=8), intent(in) :: x(n)
         logical, intent(out) :: inhdefstp,stp
         integer, intent(inout) :: ierr
         type(c_ptr), optional, intent(in) :: pdata
       end subroutine stpsub
  
       subroutine tabsub(n,f,gsupn,iter,pdata)
          use iso_c_binding, only: c_ptr
          implicit none
          integer, intent(in) :: iter,n
          real(kind=8), intent(in) :: f,gsupn
          type(c_ptr), optional, intent(in) :: pdata
        end subroutine tabsub

    end interface
  
    private
    
    ! ------------------------------------------------------------------
    ! Screen output. Edit the two values below and recompile; they change
    ! only what is printed, never the computation.
    !
    !   iprint      0  nothing from the solver
    !               2  one iteration line every iprintstep iterations, plus
    !                  the stopping reason (recommended for batteries)
    !               3  one iteration line per iteration
    !              10  everything, including the line search, the face
    !                  decisions and every CG iteration
    !   iprintstep  sampling step used at level 2
    ! ------------------------------------------------------------------
    integer, parameter, public :: iprint = 2, iprintstep = 100
    integer, parameter :: iprintevery = merge( 1, iprintstep, iprint .ge. 3 )
    integer :: nprinted = 0

    ! SCALAR PARAMETERS
    integer, parameter :: texmax = 20, iprintcg = min( 2, max( 0, iprint - 2 ) )
    real(kind=8), parameter :: lrgstreal = huge( 1.0d0 ), macheps = epsilon( 1.0d0 ), &
         macheps12 = sqrt( macheps ), macheps14 = sqrt( macheps12 ), macheps23 = macheps ** ( 2.0d0 / 3.0d0 ), beta = 0.5d0, &
         gamma = 1.0d-04, eta = 1.0d+04, lspgmin = 1.0d-16, lspgmax = 1.0d+16, r = 0.1d0, &
         tau1 = 0.1d0, tau2 = 0.9d0, sgmsmall = 1.0d-08, sgmmin = 1.0d-08, sgmmax = 1.0d+16, &
         hmin = 1.0d-08, hmax = 1.0d+08, delmin = 1.0d+04,a1 = 1.0d+8,a2 = 1.0d-16
  
    public genunc, genuncb, gencan, gencanb, project
  
  contains
    
    ! *****************************************************************
    ! *****************************************************************
  
    subroutine genunc(evalf,evalg,evalh,truehp,hnnzmax,n,x,f,g,gsupn, &
         ftarget,eps,maxit,extallowed,precond,iter,ierr,istop,evalstp,tabline,pdata)
  
      implicit none
  
      ! PROCEDURE ARGUMENTS
      procedure(fsub) :: evalf
      procedure(gsub) :: evalg
      procedure(hsub) :: evalh
      procedure(stpsub), optional :: evalstp
      procedure(tabsub), optional :: tabline
  
      ! SCALAR ARGUMENTS
      logical, intent(in) :: extallowed,precond,truehp
      integer, intent(in) :: hnnzmax,maxit,n
      integer, intent(out) :: ierr,istop,iter
      real(kind=8), intent(in) :: ftarget,eps
      real(kind=8), intent(out) :: f,gsupn
      type(c_ptr), optional, intent(in) :: pdata
  
      ! ARRAY ARGUMENTS
      real(kind=8), intent(inout) :: x(n)
      real(kind=8), intent(out) :: g(n)
  
      ! LOCAL SCALARS
      type(lssit_type) :: lssit

      ierr = 0

      call lssitini(evalg,truehp,hnnzmax,0,n,0,0,lssit)
      lssit%pdata = pdata

      call genuncb(evalf,evalg,evalh,truehp,n,x,f,g,gsupn,lssit,ftarget,eps, &
           maxit,extallowed,precond,iter,ierr,istop,evalstp,tabline,pdata)

      call lssitend(lssit) 

    end subroutine genunc
  
    ! *****************************************************************
    ! *****************************************************************
  
    subroutine genuncb(evalf,evalg,evalh,truehp,n,x,f,g,gsupn,lssit,ftarget, &
         eps,maxit,extallowed,precond,iter,ierr,istop,evalstp,tabline,pdata)
  
      implicit none
  
      ! PROCEDURE ARGUMENTS
      procedure(fsub) :: evalf
      procedure(gsub) :: evalg
      procedure(hsub) :: evalh
      procedure(stpsub), optional :: evalstp
      procedure(tabsub), optional :: tabline
  
      ! SCALAR ARGUMENTS
      logical, intent(in) :: extallowed,precond,truehp
      integer, intent(in) :: maxit,n
      integer, intent(inout) :: ierr
      integer, intent(out) :: istop,iter
      real(kind=8), intent(in) :: ftarget,eps
      real(kind=8), intent(out) :: f,gsupn
      type(c_ptr), optional, intent(in) :: pdata
  
      ! ARRAY ARGUMENTS
      real(kind=8), intent(inout) :: x(n)
      real(kind=8), intent(out) :: g(n)
      type(lssit_type), intent(inout) :: lssit
  
      ! istop =  0 means gsupn <= eps
      ! istop =  6 means f <= ftarget
      ! istop =  8 means lack of progress (f did not improve in 3 consecuitove iterations)
      ! istop = 10 means iter >= maxit
      ! istop = 20 means user-provided stopping criterion satisfied
      
      ! LOCAL SCALARS
      logical :: extinhibited,extrecom,inhdefstp,ustp,sameface
      integer :: extfails,i,navoidcy,nfnoprogr,tex, &
      cgiter,cginfo,maxitcg,inform
      real(kind=8) :: eps12,eps14,eps18,fbest,fbtrial,fext,fref,ftrial, &
           gtd,gsupnb,t,ttmp,acgeps,bcgeps,sts,sty, &
           cgeps,xsupn,ssupn,cgdel,tsmall,ysupn,normd1,beta1,beta2
    
      ! LOCAL ARRAYS
      logical :: lind(n),uind(n),frvar(n)
      integer :: countg(4),ifr(n)
      real(kind=8) :: d(n),gtrial(n),xbtrial(n), &
           xext(n),xref(n),xtrial(n),lbnd(n),ubnd(n),s(n),y(n),d2(n)
  
      !PARAMETERS
      real(kind=8), parameter :: cgepsi     =   1.0d-01
      real(kind=8), parameter :: cgepsf     =   1.0d-08
    
      ! ==========================================================================
      ! Initialization
      ! ==========================================================================
      write(*,1000) n,eps

      ! Every call starts its own table, so the header comes out again.
      nprinted = 0

      iter = 0
      lind(1:n) = .false.
      uind(1:n) = .false.
      ifr(1:n)  = (/ (i, i=1,n) /)
      y(1:n)    = 0.0d0 
  
      extfails = 0
      extinhibited = .false.
  
      navoidcy = 0
      nfnoprogr = 0
      fbest = lrgstreal
      gsupnb = lrgstreal
  
      countg(1:4) = 0
      
      eps12 = sqrt( eps )
      eps14 = sqrt( eps12 )
      eps18 = sqrt( eps14 )
  
      maxitcg = n


      call evalf(n,x,f,ierr,pdata)
      if ( ierr .ne. 0 ) return
  
      call evalg(n,x,g,ierr,pdata)
      if ( ierr .ne. 0 ) return
     
      gsupn = maxval( abs( g(1:n) ) )
      xsupn = maxval( abs( x(1:n) ) )
      if ( precond ) then 
         if ( gsupn .ne. 0.0d0 ) then
            tsmall = 1.0d-8 * max( 1.0d0, xsupn / gsupn )
         else
            tsmall = 0.0d0
         end if
         s(1:n) = x(1:n) - tsmall * g(1:n)
         call evalg(n,s,y,ierr,pdata)
         s(1:n) = s(1:n) - x(1:n)
         y(1:n) = y(1:n) - g(1:n)
         sts    = dot_product( s(1:n), s(1:n) )
         sty    = dot_product( s(1:n), y(1:n) )
         ssupn  = maxval( abs ( s(1:n) ) )
         ysupn  = maxval( abs ( y(1:n) ) )
      end if

      if ( gsupn .ne. 0.0d0 ) then
        acgeps = log10( cgepsf / cgepsi ) / log10( sqrt(cgepsf) / gsupn )
        bcgeps = log10( cgepsi ) - acgeps * log10( gsupn )
      else
        acgeps = 0.0d0
        bcgeps = cgepsf
      end if

      frvar(1:n) = .true.
     
      sameface = .false.
      ! ==========================================================================
      ! Main loop
      ! ==========================================================================
  100 continue
    
      ! Check stopping criterion
      if (iter .gt. 0) then 
        sameface = .true.
        gsupn = maxval( abs( g(1:n) ) )
      end if

      ! Verbosity control (see the block at the top of this module).
      ! Nothing below level 2; at level 2 one iteration out of iprintstep.
      ! The header follows the number of lines actually printed.
      if ( iprint .ge. 2 .and. mod( iter, iprintevery ) .eq. 0 ) then
         if ( mod(nprinted,50) .eq. 0 ) write(*,1040)
         nprinted = nprinted + 1
         if ( iter .eq. 0 ) then
            write(*,1050) iter,f,gsupn
         else
            write(*,1060) iter,f,gsupn,xsupn,ssupn,n
         end if
      end if
      if ( iprint .ge. 10 ) then
         write(*,*)
         write(*,*) 'iter = ',iter,' f = ',f,' gsupn = ',gsupn
      end if
      if ( present( tabline ) ) then
          call tabline(n,f,gsupn,iter,pdata)
          if ( ierr .ne. 0 ) return
       end if
  
      inhdefstp = .false.
      
      if ( present( evalstp ) ) then
         call evalstp(n,x,gsupn,iter,inhdefstp,ustp,ierr,pdata)
         if ( ierr .ne. 0 ) return
  
         if ( ustp ) then
            if ( iprint .ge. 2 ) write(*,*) 'User-defined stopping criterion satisfied'
            istop = 20
            return
         end if
      end if
      
      if( .not. inhdefstp )  then
         if ( gsupn .le. eps  ) then
            if ( iprint .ge. 2 ) write(*,*) 'Stopping criterion gsupn <= eps satisfied'
            istop = 0
            return
         end if
    
         if ( f .le. ftarget ) then
            if ( iprint .ge. 2 ) write(*,*) 'Stopping criterion f <= ftarget satisfied'
            istop = 6
            return
         end if
  
         if ( f .lt. fbest ) then
            nfnoprogr = 0
         else
            nfnoprogr = nfnoprogr + 1
            if ( nfnoprogr .gt. 3 ) then
               if ( iprint .ge. 2 ) write(*,*) 'Lack of progress.'
               istop = 8
               return
            end if
          end if
         
         if ( f .lt. fbest .or. ( f .eq. fbest .and. gsupn .lt. gsupnb ) ) then
            fbest = f
            gsupnb = gsupn
         end if
  
         if ( iter .ge. maxit ) then
            if ( iprint .ge. 2 ) write(*,*) 'Maximum of iterations reached.'
            istop = 10
            return
         end if
      end if
  
      ! iterate
  
      iter = iter + 1
  
      extrecom = .false.
      
      fbtrial = lrgstreal  

      if ( iter .eq. 1 ) then
         cgdel = max( 100.0d0, 100.0d0 * xsupn )
      else
         cgdel = max( delmin, 10.0d0 * ssupn )
      end if
  
      cgeps   = 10.0d0 ** ( acgeps * log10( gsupn ) + bcgeps )
      cgeps   = max( cgepsf, min( cgepsi, cgeps ) )

      if ( iprint .ge. 10 ) then
         write (*,*)
         write (*,*) "cgeps   = ", cgeps
         write (*,*) "maxitcg = ", maxitcg
      end if

      if ( truehp ) then 
         lssit%compressed = .true.
         call evalh(n,x,lssit%hnnzmax,lssit%hne,lssit%hrow,lssit%hcol,lssit%hval,ierr,pdata)
         if ( ierr .ne. 0 ) return
      else
         call comphapp(n,s(1:n),y(1:n),sty,sts,lspgmin,lspgmax,sameface,lssit)
      end if

      if ( precond ) then 
        call comphpre(n,s(1:n),y(1:n),sty,sts,lspgmin,lspgmax,sameface,inform,lssit)
      end if

      if ( truehp ) then
         call preparehv(n,n,frvar,lssit)
      end if

      if ( precond .and. inform .eq. 0 ) then 
         call cgm (n,n,ifr(1:n),x(1:n),lind,lbnd,uind,ubnd,g(1:n),cgdel,cgeps,maxitcg,&
               precond,d(1:n),cgiter,iprintcg,ierr,cginfo,lssit)
      else
         call cgm (n,n,ifr(1:n),x(1:n),lind,lbnd,uind,ubnd,g(1:n),cgdel,cgeps,maxitcg,&
               .false.,d(1:n),cgiter,iprintcg,ierr,cginfo,lssit)
      end if
      normd1 = norm2( d(1:n) )
      if ( normd1 .gt. a1*gsupn ) then
         beta1 = ( a1 * gsupn ) / normd1
         d2(1:n) = beta1 * d(1:n)
      else
         d2(1:n) = d(1:n)
      end if
      gtd = dot_product( g(1:n), d2(1:n) )
      if ( gtd .gt. -a2 * ( gsupn**2 ) ) then 
         beta2 = ( 1 - a2 ) / ( 1 + ( gtd / ( gsupn**2 ) ) )
         d(1:n) = beta2 * d2(1:n) + ( 1.0d0 - beta2 ) * ( -g(1:n) )
      else
         d(1:n) = d2(1:n)
      end if
      gtd = dot_product( g(1:n), d(1:n) )

      if ( iprint .ge. 10 ) write(*,*) 'Newton direction dsupn = ',maxval( abs( d(1:n) ) ),' deucn = ',norm2( d(1:n) )

      go to 300
     
      ! Backtracking
  300 continue

      if ( iprint .ge. 10 ) write(*,*) 'gtd = ',gtd
        
      t = 1.0d0
  
  310 continue
        
      xtrial(1:n) = x(1:n) + t * d(1:n)

      call evalf(n,xtrial,ftrial,ierr,pdata)
      if ( ierr .ne. 0 ) return
      
      if ( iprint .ge. 10 ) write(*,*) 'in backtracking, t = ',t,' ftrial = ',ftrial
  
      if ( ftrial .lt. fbtrial ) then
         fbtrial = ftrial
         xbtrial(1:n) = xtrial(1:n)
      end if
        
      if ( t .eq. 0.0d0 ) then
         if ( iprint .ge. 10 ) write(*,*) 'null step in backtracking'
         xtrial(1:n) = x(1:n)
         ftrial = f
         go to 400
         
      else if ( ftrial .le. ftarget ) then
         if ( iprint .ge. 10 ) write(*,*) 'ftrial <= ftarget'
         go to 400
         
      else if ( ftrial .le. f + t * gamma * gtd ) then
         if ( iprint .ge. 10 ) write(*,*) 'Armijo holds'
         if ( t .eq. 1.0d0 ) then
            extrecom = .true.
         end if
         go to 400
  
      else if ( t .eq. 1.0d0 .and. &
           ftrial .le. f + macheps14 * max( 1.0d0, abs( f ) ) .and. &
           norm2( d(1:n) ) .le. macheps14 * norm2( x(1:n) ) ) then
         if ( navoidcy .lt. 2 ) then
            navoidcy = navoidcy + 1
            if ( iprint .ge. 10 ) write(*,*) 'The Newton step did NOT give sufficient decrease.'
            if ( iprint .ge. 10 ) write(*,*) 'However, the Euclidian norm of the step is small and the change in f is negligible.'
            if ( iprint .ge. 10 ) write(*,*) 'Therefore, we accept it anyway.'
            go to 400
         end if
      end if
  
      ttmp = ( - gtd * t ** 2 ) / ( 2.0d0 * ( ftrial - f - t * gtd ) )
        
      if ( .not. ( tau1 * t .le. ttmp .and. ttmp .le. tau2 * t ) ) then
         t = 0.5d0 * t
      else
         t = ttmp
      end if
      
      go to 310
  
  400 continue
      
      if ( fbtrial .lt. ftrial ) then
         ftrial = fbtrial
         xtrial(1:n) = xbtrial(1:n)
      end if
      
      call evalg(n,xtrial,gtrial,ierr,pdata)
      if ( ierr .ne. 0 ) return
  
      ! Extrapolate if recommended
      
      if ( extallowed .and. .not. extinhibited .and. extrecom ) then
         if ( .not. maxval( abs( gtrial(1:n) ) ) .le. eps ) then
            if ( dot_product( gtrial(1:n), d(1:n) ) .le. beta * gtd ) then     
               fref = ftrial
               xref(1:n) = xtrial(1:n)
      
               tex = 1
410 continue
               if ( tex .le. texmax ) then
                  if ( .not. fref .le. ftarget ) then
                     xext(1:n) = x(1:n) + 2.0d0 ** tex * d(1:n)

                     call evalf(n,xext,fext,ierr,pdata)
                     if ( ierr .ne. 0 ) return
                     
                     if ( fext .lt. fref ) then
                        fref = fext
                        xref(1:n) = xext(1:n)
                        
                        if ( iprint .ge. 10 ) write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' gave descent fext = ',fext
                        
                        tex = tex + 1
                        go to 410
                     else
                        if ( iprint .ge. 10 ) write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' failed. fext = ', fext
                        if ( tex .gt. 1 )  then
                           extfails = 0
                        else
                           extfails = extfails + 1
                           if ( extfails .eq. 3 ) then
                              extinhibited = .true.
                           end if
                        end if
                     end if
                  else
                     if ( iprint .ge. 10 ) write(*,*) 'fref <= ftarget'
                  end if
               else
                  if ( iprint .ge. 10 ) write(*,*) 'tex > texmax (tex =',tex,', texmax = ',texmax,')'
               end if
  
               if ( fref .lt. ftrial ) then
                  ftrial = fref
                  xtrial(1:n) = xref(1:n)
                  
                  call evalg(n,xtrial,gtrial,ierr,pdata)
                  if ( ierr .ne. 0 ) return
               end if
            end if
         end if
      end if
   
      if ( iprint .ge. 10 ) write(*,*) 'f = ', f, 'ftrial = ', ftrial
      f = ftrial
      s(1:n) = xtrial(1:n) - x(1:n)
      y(1:n) = gtrial(1:n) - g(1:n)
      x(1:n) = xtrial(1:n)
      g(1:n) = gtrial(1:n)
      
      xsupn = maxval( abs( x(1:n) ) )
      ssupn = maxval( abs( s(1:n) ) )
      ysupn = maxval( abs( y(1:n) ) )
      sts   = dot_product( s(1:n), s(1:n) )
      sty   = dot_product( s(1:n), y(1:n) )
  
      go to 100  

   1040 format(/,5X,'It',9X,'function',2X,' gsupn',3X,'xsupn',1X,'stesupn',7X,'n')
   1050 format(     I7,1X,1P,D16.8,1X,1P,D7.1)
   1060 format(     I7,1X,1P,D16.8,3(1X,1P,D7.1),1X,I7)

   1000 format(/,5X,'Entry to gencan.',         &
   /,5X,'Number of variables: ',I7, &
   /,5X,'Required optimality tolerance: ',1P,D7.1)
    end subroutine genuncb
    
    ! *****************************************************************
    ! *****************************************************************
      
    subroutine gencan(evalf,evalg,evalh,truehp,hnnzmax,n,x,lind, &
      lbnd,uind,ubnd,f,g,gpsupn,ftarget,eps,maxit,extallowed,precond,iter, &
      ierr,istop,evalstp,tabline,pdata)
  
   implicit none
  
   ! PROCEDURE ARGUMENTS
   procedure(fsub) :: evalf
   procedure(gsub) :: evalg
   procedure(hsub) :: evalh
   procedure(stpsub), optional :: evalstp
   procedure(tabsub), optional :: tabline
  
   ! SCALAR ARGUMENTS
   logical, intent(in) :: extallowed,precond,truehp
   integer, intent(in) :: hnnzmax,maxit,n
   integer, intent(out) :: ierr,istop,iter
   real(kind=8), intent(in) :: ftarget,eps
   real(kind=8), intent(out) :: f,gpsupn
   type(c_ptr), optional, intent(in) :: pdata
  
   ! ARRAY ARGUMENTS
   logical, intent(in) :: lind(n),uind(n)
   real(kind=8), intent(in) :: lbnd(n),ubnd(n)
   real(kind=8), intent(inout) :: x(n)
   real(kind=8), intent(out) :: g(n)
  
   ! LOCAL SCALARS
   type(lssit_type) :: lssit
  
    ierr = 0

      call lssitini(evalg,truehp,hnnzmax,0,n,0,0,lssit)
    lssit%pdata = pdata

    call gencanb(evalf,evalg,evalh,truehp,n,x,lind,lbnd,uind,ubnd, &
            f,g,gpsupn,lssit,ftarget,eps,maxit,extallowed,precond,iter,ierr,istop, &
            evalstp,tabline,pdata)
    
    call lssitend(lssit)

  end subroutine gencan
  
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine gencanb(evalf,evalg,evalh,truehp,n,x,lind,lbnd,uind,ubnd, &
      f,g,gpsupn,lssit,ftarget,eps,maxit,extallowed,precond,iter,ierr,istop,evalstp, &
      tabline,pdata)
  
   implicit none
  
   ! PROCEDURE ARGUMENTS
   procedure(fsub) :: evalf
   procedure(gsub) :: evalg
   procedure(hsub) :: evalh
   procedure(stpsub), optional :: evalstp
   procedure(tabsub), optional :: tabline
  
   ! SCALAR ARGUMENTS
   logical, intent(in) :: extallowed,precond,truehp
   integer, intent(in) :: maxit,n
   integer, intent(inout) :: ierr
   integer, intent(out) :: istop,iter
   real(kind=8), intent(in) :: ftarget,eps
   real(kind=8), intent(out) :: f,gpsupn
   type(c_ptr), optional, intent(in) :: pdata
  
   ! ARRAY ARGUMENTS
   logical, intent(in) :: lind(n),uind(n)
   real(kind=8), intent(in) :: lbnd(n),ubnd(n)
   real(kind=8), intent(inout) :: x(n)
   real(kind=8), intent(out) :: g(n)
   type(lssit_type), intent(inout) :: lssit
  
   ! Solves the box-constrained minimization problem
   !
   !        Minimize f(x) subject to l <= x <= u
   !
   ! istop =  0 means gsupn <= eps
   ! istop =  6 means f <= ftarget
   ! istop =  8 means lack of progress (f did not improve in 3 consecuitove iterations)
   ! istop = 10 means iter >= maxit
   ! istop = 20 means user-provided stopping criterion satisfied
   
   ! LOCAL SCALARS
   logical :: extinhibited,extrecom,infrontier,inhdefstp,newface,ustp,samep
   integer :: extfails,i,navoidcy,nfnoprogr,nfr,tex
   real(kind=8) :: eps12,eps14,eps18,fbest,fbtrial,fext,fref,ftrial,gisupn,gtd, &
        gpsupnb,lspg,sts,sty,t,ttmp,acgeps,bcgeps,cgeps,&
        xsupn,gsupn,tref,ssupn,cgdel,tsmall,ysupn
  
   ! LOCAL ARRAYS
   logical :: frvar(n),fvpre(n)
   integer :: countg(4),ifr(n)
   real(kind=8) :: d(n),gp(n),gprev(n),gptrial(n),gtrial(n),s(n),y(n),xbtrial(n), &
        xext(n),xprev(n),xref(n),xtrial(n)

   !PARAMETERS
   real(kind=8), parameter :: cgepsi     =   1.0d-01
   real(kind=8), parameter :: cgepsf     =   1.0d-08
  
   ! ==========================================================================
   ! Initialization
   ! ==========================================================================
   write(*,1000) n,eps

   ! Every call starts its own table, so the header comes out again.
   nprinted = 0

   iter = 0
  
   extfails = 0
   extinhibited = .false.
   
   navoidcy = 0
   nfnoprogr = 0
   fbest = lrgstreal
   gpsupnb = lrgstreal
  
   countg(1:4) = 0
   
   eps12 = sqrt( eps )
   eps14 = sqrt( eps12 )
   eps18 = sqrt( eps14 )

   call project(n,lind,lbnd,uind,ubnd,x)
  
   call evalf(n,x,f,ierr,pdata)
   if ( ierr .ne. 0 ) return
  
   call evalg(n,x,g,ierr,pdata)
   if ( ierr .ne. 0 ) return

   gp(1:n) = x(1:n) - g(1:n)

   call project(n,lind,lbnd,uind,ubnd,gp)

   gp(1:n) = gp(1:n) - x(1:n)

   gpsupn = maxval( abs( gp(1:n) ) )   
   xsupn  = maxval( abs( x(1:n) ) )
   if ( precond ) then 
      if ( gpsupn .ne. 0.0d0 ) then
         tsmall = macheps12 * max( 1.0d0, xsupn / gpsupn )
      else
         tsmall = 0.0d0
      end if
      s(1:n) = x(1:n) + tsmall * gp(1:n)
         call evalg(n,s,y,ierr,pdata)
      s(1:n) = s(1:n) - x(1:n)
      y(1:n) = y(1:n) - g(1:n)
      sts    = dot_product( s(1:n), s(1:n) )
      sty    = dot_product( s(1:n), y(1:n) )
      ssupn  = maxval( abs ( s(1:n) ) )
      ysupn  = maxval( abs ( y(1:n) ) )
   end if

   if ( gpsupn .ne. 0.0d0 ) then
      acgeps = log10( cgepsf / cgepsi ) / log10( sqrt(cgepsf) / gpsupn )
      bcgeps = log10( cgepsi ) - acgeps * log10( gpsupn )
   else
      acgeps = 0.0d0
      bcgeps = cgepsf
   end if
  
   ! ==========================================================================
   ! Main loop
   ! ==========================================================================

  100 continue
  
   ! Check stopping criterion

    if(iter .gt. 0) then 
        gp(1:n) = x(1:n) - g(1:n)
        
        call project(n,lind,lbnd,uind,ubnd,gp)

        gp(1:n) = gp(1:n) - x(1:n)

        gpsupn = maxval( abs( gp(1:n) ) )
        xsupn = maxval( abs( x(1:n) ) )
        ssupn = maxval( abs( s(1:n) ) )
        ysupn  = maxval( abs ( y(1:n) ) )
    end if
    gsupn = maxval( abs( g(1:n) ) )

   ! Verbosity control (see the block at the top of this module).
   if ( iprint .ge. 2 .and. mod( iter, iprintevery ) .eq. 0 ) then
      if ( mod(nprinted,50) .eq. 0 ) write(*,1040)
      nprinted = nprinted + 1
      if ( iter .eq. 0 ) then
         write(*,1050) iter,f,gsupn
      else
         write(*,1060) iter,f,gpsupn,xsupn,ssupn,nfr
      end if
   end if
   if ( iprint .ge. 10 ) then
      write(*,*)
      write(*,*) 'iter = ',iter,' f = ',f,' gpsupn = ',gpsupn
   end if
   if ( present( tabline ) ) then
      call tabline(n,f,gpsupn,iter,pdata)
      if ( ierr .ne. 0 ) return
   end if
  
   inhdefstp = .false.
   
   if ( present( evalstp ) ) then
      call evalstp(n,x,gpsupn,iter,inhdefstp,ustp,ierr,pdata)
      if ( ierr .ne. 0 ) return
  
      if ( ustp ) then
         if ( iprint .ge. 2 ) write(*,*) 'User-defined stopping criterion satisfied'
         istop = 20
         return
      end if
   end if
  
   if ( .not. inhdefstp ) then
      if ( gpsupn .le. eps ) then
         if ( iprint .ge. 2 ) write(*,*) 'Stopping criterion gpsupn <= eps satisfied'
         istop = 0
         return
      end if
      
      if ( f .le. ftarget ) then
         if ( iprint .ge. 2 ) write(*,*) 'Stopping criterion f <= ftarget satisfied'
         istop = 6
         return
      end if
  
      if ( f .lt. fbest ) then
         nfnoprogr = 0
      else
         nfnoprogr = nfnoprogr + 1
         if ( nfnoprogr .gt. 3 ) then
            if ( iprint .ge. 2 ) write(*,*) 'Lack of progress.'
            istop = 8
            return
         end if
      end if
      
      if ( f .lt. fbest .or. ( f .eq. fbest .and. gpsupn .lt. gpsupnb ) ) then
         fbest = f
         gpsupnb = gpsupn
      end if
  
      if ( iter .ge. maxit ) then
         if ( iprint .ge. 2 ) write(*,*) 'Maximum of iterations reached.'
         istop = 10
         return
      end if
   end if
   
   ! opt by an inner-to-the-face or a leaving-face iteration
  
   frvar(1:n) = ( .not. lind(1:n) .or. ( x(1:n) - lbnd(1:n) ) / max( 1.0d0, lbnd(1:n) ) .gt. 1.0d-12 ) .and. &
                ( .not. uind(1:n) .or. ( ubnd(1:n) - x(1:n) ) / max( 1.0d0, ubnd(1:n) ) .gt. 1.0d-12 )
  
   newface = .false.
   if ( iter .eq. 0 .or. any( frvar(1:n) .neqv. fvpre(1:n) .or. &
        ( .not. frvar(1:n) .and. x(1:n) .ne. xprev(1:n) ) ) ) then
      newface = .true.
   end if

   ! iterate
  
   iter = iter + 1
  
   extrecom = .false.
   
   xprev(1:n) = x(1:n)
   gprev(1:n) = g(1:n)
   fvpre(1:n) = frvar(1:n)
  
   gisupn = max( 0.0d0, maxval( abs( gp(1:n) ), frvar(1:n) ) )
  
   if ( iprint .ge. 10 ) write(*,*) 'gisupn = ',gisupn
  
   nfr = count( frvar(1:n) )
  
   if ( iprint .ge. 10 ) write(*,*) 'Number of free variables: ',nfr
  
   if ( gisupn .ge. r * gpsupn .and. norm2( pack(gp, frvar) ) .ge. r * norm2(gp) ) then
  
      ifr(1:nfr) = pack( (/ (i, i=1,n) /), frvar(1:n) )
  
      if ( iprint .ge. 10 ) write(*,*) 'An inner-to-the-face iteration will be done'
  
      call newtonls()
      if ( ierr .ne. 0 ) return
   
      if ( fbtrial .lt. ftrial ) then
         ftrial = fbtrial
         xtrial(1:n) = xbtrial(1:n)
      end if
      
   else
  
      if ( iprint .ge. 10 ) write(*,*) 'A leaving-face iteration will be done'
  
      if ( iter .eq. 1 .or. sty .le. 0.0d0 ) then
         lspg = max( 1.0d0, xsupn ) / gsupn
      else
         lspg = sts / sty
      end if
  
      lspg = max( lspgmin, min( lspg, lspgmax ) )
  
      t = 1.0d0
  
      xtrial(1:n) = x(1:n) - lspg * g(1:n)
      call project(n,lind,lbnd,uind,ubnd,xtrial)
  
      d(1:n) = xtrial(1:n) - x(1:n)
  
      gtd = dot_product( g(1:n), d(1:n) )
  
  200    continue
    
      call evalf(n,xtrial,ftrial,ierr,pdata)
      if ( ierr .ne. 0 ) return
  
      if ( t .eq. 0.0d0 ) then
         if ( iprint .ge. 10 ) write(*,*) 'null step in backtracking'
         xtrial(1:n) = x(1:n)
         ftrial = f
         
      else if ( ftrial .le. ftarget ) then
         if ( iprint .ge. 10 ) write(*,*) 'ftrial <= ftarget'
  
      else if ( .not. ( ftrial .le. f + t * gamma * gtd ) ) then
         ttmp = ( - gtd * t ** 2 ) / ( 2.0d0 * ( ftrial - f - t * gtd ) )
         
         if ( .not. ( tau1 * t .le. ttmp .and. ttmp .le. tau2 * t ) ) then
            t = 0.5d0 * t
         else
            t = ttmp
         end if
  
         xtrial(1:n) = x(1:n) + t * d(1:n)
  
         go to 200
      end if
  
   end if
  
   call evalg(n,xtrial,gtrial,ierr,pdata)
   if ( ierr .ne. 0 ) return
  
   if ( extallowed .and. .not. extinhibited .and. extrecom ) then
      gptrial(1:n) = xtrial(1:n) - gtrial(1:n)
      call project(n,lind,lbnd,uind,ubnd,gptrial)
      gptrial(1:n) = gptrial(1:n) - xtrial(1:n)
      if ( .not. maxval( abs( gptrial(1:n) ) ) .le. eps ) then
         if ( infrontier .or. dot_product( gtrial(ifr(1:nfr)), d(1:nfr) ) .le. beta * gtd ) then          
            fref = ftrial
            xref(1:n) = xtrial(1:n)
            xext(1:n) = x(1:n)
   
            tex = 1
320 continue
            if ( tex .le. texmax ) then
               if ( .not. fref .le. ftarget ) then
                  xext(ifr(1:nfr)) = x(ifr(1:nfr)) + ( 2.0d0 ** tex )* tref * d(1:nfr)
                  call project(n,lind,lbnd,uind,ubnd,xext)

                  samep = all( xext(ifr(1:nfr)) .le. xref(ifr(1:nfr)) + macheps23 * max( 1.0d0, abs( xref(ifr(1:nfr)) ) ) .and. &
                               xext(ifr(1:nfr)) .ge. xref(ifr(1:nfr)) - macheps23 * max( 1.0d0, abs( xref(ifr(1:nfr)) ) ) )

                  if ( .not. samep ) then 
                    call evalf(n,xext,fext,ierr,pdata)
                    if ( ierr .ne. 0 ) return
    
                    if ( fext .lt. fref  ) then
                        fref = fext
                        xref(1:n) = xext(1:n)
                        
                        if ( iprint .ge. 10 ) write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' gave descent fext = ',fext
                        
                        tex = tex + 1
                        go to 320
                    else
                        if ( iprint .ge. 10 ) write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' failed. fext = ', fext
                        if ( tex .gt. 1 )  then
                            extfails = 0
                        else
                            extfails = extfails + 1
                            if ( extfails .eq. 3 ) then
                            extinhibited = .true.
                            end if
                        end if
                    end if
                  else
                    if ( iprint .ge. 10 ) write(*,*) 'Flag of Extrapolation: Very similar consecutive projected points.'
                  end if
               else
                  if ( iprint .ge. 10 ) write(*,*) 'fref <= ftarget'
               end if
            else
               if ( iprint .ge. 10 ) write(*,*) 'tex > texmax (tex =',tex,', texmax = ',texmax,')'
            end if
  
            if ( fref .lt. ftrial ) then
               ftrial = fref
               xtrial(1:n) = xref(1:n)
               
               call evalg(n,xtrial,gtrial,ierr,pdata)
               if ( ierr .ne. 0 ) return
            end if
         end if
      else
            if ( iprint .ge. 10 ) write(*,*) 'maxval( abs( gptrial(1:n) ) ) <= eps then no extrapolation is required.'
      end if
   end if
  
   f = ftrial
   x(1:n) = xtrial(1:n)
   g(1:n) = gtrial(1:n)
   
   s(1:n) = x(1:n) - xprev(1:n)
   y(1:n) = g(1:n) - gprev(1:n)
   sty = dot_product( s(1:n), y(1:n) )
   sts = sum( s(1:n) ** 2 )
  
   go to 100

   1040 format(/,5X,'It',9X,'function',2X,'gpsupn',3X,'xsupn',1X,'stesupn',3X,'nfree')
   1050 format(     I7,1X,1P,D16.8,1X,1P,D7.1)
   1060 format(     I7,1X,1P,D16.8,3(1X,1P,D7.1),1X,I7)

   1000 format(/,5X,'Entry to gencan.',         &
   /,5X,'Number of variables: ',I7, &
   /,5X,'Required optimality tolerance: ',1P,D7.1)
  
  contains
   
   ! *****************************************************************
   ! *****************************************************************
   
   subroutine newtonls()
  
     implicit none
     
     ! LOCAL SCALARS
     integer :: i,rbdnnz,cgiter,maxitcg,cginfo,inform
     real(kind=8) :: amax,normd1,beta1,beta2

     ! LOCAL ARRAYS
     integer :: rbdind(nfr)
     real(kind=8) :: amaxv(nfr),d2(n)
     character :: amaxt(nfr),rbdtype(nfr)

     fbtrial = lrgstreal
     maxitcg = nfr

              
     if ( iter .eq. 1 ) then
        cgdel = max( 100.0d0, 100.0d0 * xsupn )
     else
        cgdel = max( delmin, 10.0d0 * ssupn )
     end if

     cgeps   = 10.0d0 ** ( acgeps * log10( gpsupn ) + bcgeps )
     cgeps   = max( cgepsf, min( cgepsi, cgeps ) )

     if ( truehp ) then 
         lssit%compressed = .true.
         call evalh(n,x,lssit%hnnzmax,lssit%hne,lssit%hrow,lssit%hcol,lssit%hval,ierr,pdata)
         if ( ierr .ne. 0 ) return
     else
         call comphapp(n,s(1:n),y(1:n),sty,sts,lspgmin,lspgmax,.not. newface,lssit)
     end if

     if ( precond ) then 
        call comphpre(n,s(1:n),y(1:n),sty,sts,lspgmin,lspgmax,.not. newface,inform,lssit)
     end if

     if ( truehp ) then 
       call preparehv(n,nfr,frvar,lssit)
     end if

     if ( precond .and. inform .eq. 0 ) then 
         call cgm (n,nfr,ifr(1:nfr),x(1:n),lind,lbnd,uind,ubnd,g(ifr(1:nfr)),cgdel,cgeps,maxitcg,&
             precond,d(1:n),cgiter,iprintcg,ierr,cginfo,lssit)
     else
         call cgm (n,nfr,ifr(1:nfr),x(1:n),lind,lbnd,uind,ubnd,g(ifr(1:nfr)),cgdel,cgeps,maxitcg,&
              .false.,d(1:n),cgiter,iprintcg,ierr,cginfo,lssit)
     end if

      normd1 = norm2( d(1:nfr) ) 
      if ( normd1 .gt. a1*gpsupn ) then
         beta1 = ( a1 * gpsupn ) / normd1
         d2(1:nfr) = beta1 * d(1:nfr)
      else
         d2(1:nfr) = d(1:nfr)
      end if
      gtd = dot_product( g(ifr(1:nfr)), d2(1:nfr) )

      if ( gtd .gt. -a2 * ( gpsupn**2 ) ) then 
         beta2 = ( 1 - a2 ) / (1 + ( gtd / ( gpsupn**2 ) ) )
         d(1:nfr) = beta2 * d2(1:nfr) + (1.0d0 - beta2) * ( -g(ifr(1:nfr)) )
      else
         d(1:nfr) = d2(1:nfr)
      end if
     gtd = dot_product( g(ifr(1:nfr)), d(1:nfr) )

     if ( iprint .ge. 10 ) write(*,*) 'dsupn = ',maxval( abs( d(1:nfr) ) ),' deucn = ',norm2( d(1:nfr) )

     ! Compute maximum feasible step
  
     amaxv(1:nfr) = lrgstreal
     
     where ( d(1:nfr) .lt. 0.0d0 .and. lind(ifr(1:nfr)) )
        amaxt(1:nfr) = 'L'
        amaxv(1:nfr) = ( lbnd(ifr(1:nfr)) - x(ifr(1:nfr)) ) / d(1:nfr)
     end where
  
     where ( d(1:nfr) .gt. 0.0d0 .and. uind(ifr(1:nfr)) )
        amaxt(1:nfr) = 'U'
        amaxv(1:nfr) = ( ubnd(ifr(1:nfr)) - x(ifr(1:nfr)) ) / d(1:nfr)
     end where
  
     amax = minval( amaxv(1:nfr) )
  
     if ( iprint .ge. 10 ) write(*,*) 'amax = ',amax
  
     if ( amax .lt. lrgstreal ) then
        rbdnnz = count( amaxv(1:nfr) .eq. amax )
        rbdind(1:rbdnnz) = pack( (/ (i,i=1,nfr) /), amaxv(1:nfr) .eq. amax )
        rbdtype(1:rbdnnz) = pack( amaxt(1:nfr), amaxv(1:nfr) .eq. amax )
        if ( iprint .ge. 10 ) write(*,*) 'rbdnnz = ',rbdnnz!,' rbdind = ',rbdind(1:rbdnnz),' rbdtype = ',rbdtype(1:rbdnnz)
     else
        rbdnnz = 0
        if ( iprint .ge. 10 ) write(*,*) 'rbdnnz = ',rbdnnz
     end if
     
     ! Test projected unitary step
     
     if ( amax .lt. 1.0d0 ) then
        if ( iprint .ge. 10 ) write(*,*) 'Unitary step is infeasible. Trying projected unitary step.'
        
        xtrial(1:n) = x(1:n)
        xtrial(ifr(1:nfr)) = xtrial(ifr(1:nfr)) + d(1:nfr)
  
        call projinface(n,lind,lbnd,uind,ubnd,xtrial,nfr,ifr)
       
        call evalf(n,xtrial,ftrial,ierr,pdata)
        if ( ierr .ne. 0 ) return
   
        if ( iprint .ge. 10 ) write(*,*) 'ftrial = ',ftrial
  
        if ( ftrial .lt. fbtrial ) then
           fbtrial = ftrial
           xbtrial(1:n) = xtrial(1:n)
        end if
        
        if ( ftrial .le. ftarget ) then
           if ( iprint .ge. 10 ) write(*,*) 'ftrial <= ftarget'
           return
  
        else if ( ftrial .le. f ) then
           if ( iprint .ge. 10 ) write(*,*) 'Simple decrease was obtained and since we are in the frontier, this is enough.'
           extrecom   = .true.
           infrontier = .true.
           tref = 1.0d0
           return

        else if ( ftrial .le. f + macheps14 * max( 1.0d0, abs( f ) ) .and. &
            norm2( xtrial(ifr(1:nfr)) - x(ifr(1:nfr)) ) .le. macheps14 * norm2( x(ifr(1:nfr)) ) ) then
          if ( navoidcy .lt. 2 ) then
             navoidcy = navoidcy + 1
             if ( iprint .ge. 10 ) write(*,*) 'The unitary step is infeasible and its projection did NOT give simple decrease.'
             if ( iprint .ge. 10 ) write(*,*) 'However, the Euclidian norm of the step is small and the change in f is negligible.'
             if ( iprint .ge. 10 ) write(*,*) 'Therefore, we accept it anyway.'
             return
          end if
       end if
     end if
  
     ! Backtracking  
     if ( iprint .ge. 10 ) write(*,*) 'gtd = ',gtd
     
     t = min( 1.0d0, amax )
  
  310   continue
     
     xtrial(1:n) = x(1:n)
     xtrial(ifr(1:nfr)) = xtrial(ifr(1:nfr)) + t * d(1:nfr)
  
     if ( t .eq. amax ) then
        where ( rbdtype(1:rbdnnz) .eq. 'L' )
           xtrial(ifr(rbdind(1:rbdnnz))) = lbnd(ifr(rbdind(1:rbdnnz)))
        elsewhere
           xtrial(ifr(rbdind(1:rbdnnz))) = ubnd(ifr(rbdind(1:rbdnnz)))
        end where
     end if
     
     call evalf(n,xtrial,ftrial,ierr,pdata)
     if ( ierr .ne. 0 ) return
   
     if ( iprint .ge. 10 ) write(*,*) 'in backtracking, t = ',t,' ftrial = ',ftrial
  
     if ( ftrial .lt. fbtrial ) then
        fbtrial = ftrial
        xbtrial(1:n) = xtrial(1:n)
     end if
     
      if ( t .eq. 0.0d0 ) then
         if ( iprint .ge. 10 ) write(*,*) 'null step in backtracking'
         xtrial(1:n) = x(1:n)
         ftrial = f
         return
         
     else if ( ftrial .le. ftarget ) then
        if ( iprint .ge. 10 ) write(*,*) 'ftrial <= ftarget'
        return
  
     else if ( t .eq. amax .and. ftrial .le. f ) then ! Step 2.3
        if ( iprint .ge. 10 ) write(*,*) 'Simple decrease was obtained and since we are in the frontier, this is enough.'
        extrecom   = .true.
        infrontier = .true.
        tref = amax
        return
  
     else if ( ftrial .le. f + t * gamma * gtd ) then
        if ( iprint .ge. 10 ) write(*,*) 'Armijo holds'
        if ( t .eq. min( 1.0d0, amax ) ) then
           extrecom   = .true.
           infrontier = .false.
           tref = t ! = min(1.0d0, amax)
        end if
        return

     else if ( t .eq. min( 1.0d0, amax ) .and. &
        ftrial .le. f + macheps14 * max( 1.0d0, abs( f ) ) .and. &
        norm2( xtrial(ifr(1:nfr)) - x(ifr(1:nfr)) ) .le. macheps14 * norm2( x(ifr(1:nfr)) ) ) then
      if ( navoidcy .lt. 2 ) then
            navoidcy = navoidcy + 1
            if ( iprint .ge. 10 ) write(*,*) 'The unitary or the maximum step did NOT give sufficient decrease.'
            if ( iprint .ge. 10 ) write(*,*) 'However, the Euclidian norm of the step is small and the change in f is negligible.'
            if ( iprint .ge. 10 ) write(*,*) 'Therefore, we accept it anyway.'
            return
        end if
     end if
  
     ttmp = ( - gtd * t ** 2 ) / ( 2.0d0 * ( ftrial - f - t * gtd ) )
     
     if ( .not. ( tau1 * t .le. ttmp .and. ttmp .le. tau2 * t ) ) then
        t = 0.5d0 * t
     else
        t = ttmp
     end if
  
     go to 310
  
   end subroutine newtonls
  
  end subroutine gencanb
  
  
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine projinface(n,lind,lbnd,uind,ubnd,x,nind,ind)
  
   implicit none
  
   ! SCALAR ARGUMENTS
   integer, intent(in) :: n,nind
  
   ! ARRAY ARGUMENTS
   logical, intent(in) :: lind(n),uind(n)
   integer, intent(in) :: ind(nind)
   real(kind=8), intent(in) :: lbnd(n),ubnd(n)
   real(kind=8), intent(inout) :: x(n)
  
   where ( lind(ind(1:nind)) )
      x(ind(1:nind)) = max( lbnd(ind(1:nind)), x(ind(1:nind)) )
   end where
  
   where ( uind(ind(1:nind)) )
      x(ind(1:nind)) = min( ubnd(ind(1:nind)), x(ind(1:nind)) )
   end where
  
  end subroutine projinface
  
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine project(n,lind,lbnd,uind,ubnd,x)
  
   implicit none
  
   ! SCALAR ARGUMENTS
   integer, intent(in) :: n
  
   ! ARRAY ARGUMENTS
   logical, intent(in) :: lind(n),uind(n)
   real(kind=8), intent(in) :: lbnd(n),ubnd(n)
   real(kind=8), intent(inout) :: x(n)
  
   where ( lind(1:n) )
      x(1:n) = max( lbnd(1:n), x(1:n) )
   end where
  
   where ( uind(1:n) )
      x(1:n) = min( ubnd(1:n), x(1:n) )
   end where
  
  end subroutine project


  end module modgencan
  