module lssit
  use iso_c_binding, only: c_ptr
  
  implicit none

  interface
      subroutine gsub(n,x,g,ierr,pdata)
      use iso_c_binding, only: c_ptr
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: g(n)
      integer, intent(inout) :: ierr
      type(c_ptr), optional, intent(in) :: pdata
    end subroutine gsub
  end interface

  type lssit_type
     procedure(gsub), pointer, nopass  :: evalg

     type(c_ptr) :: pdata

     integer, allocatable :: hrow(:), hcol(:),jsta(:),jlen(:),jvar(:)
     real(kind=8), allocatable :: hval(:),jval(:),s(:),y(:),pdiag(:),psmdy(:),hds(:)
     logical, allocatable :: constrin(:)
     integer :: n,m,p,hnnzmax,hne,jnnzmax,jne,infohapp
     logical :: compressed,sameface,truehp
     real(kind=8) :: rho,sty,seucn,yeucn,plspg,psmdyty,hstds,hlspg
  end type lssit_type

  private
  real(kind=8), parameter ::  macheps = epsilon( 1.0d0 ), macheps12 = sqrt( macheps )
  
  public :: lssit_type, lssitini, lssitend, preparehv, computehv, comphpre,applyhpre, &
            comphapp
  
contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine lssitini(evalg,truehp,hnnzmax,jnnzmax,n,m,p,lssit)

    implicit none

    ! PROCEDURE ARGUMENTS
    procedure(gsub) :: evalg

    ! SCALAR ARGUMENTS
    integer, intent(in) :: hnnzmax,jnnzmax,n,m,p
    logical, intent(in) :: truehp
    type(lssit_type), intent(out) :: lssit

    ! LOCAL SCALARS
    integer :: allocerr

    lssit%n  = n
    lssit%m  = m
    lssit%p  = p
    lssit%truehp = truehp

    lssit%hnnzmax  = hnnzmax
    lssit%jnnzmax  = jnnzmax

    lssit%evalg => evalg 

    allocate(lssit%hrow(hnnzmax),lssit%hcol(hnnzmax),lssit%hval(hnnzmax),&
             lssit%jsta(m+p),lssit%jlen(m+p),lssit%jvar(jnnzmax),&
             lssit%jval(jnnzmax),lssit%constrin(m+p),lssit%s(n),&
             lssit%y(n),lssit%pdiag(n),lssit%psmdy(n),lssit%hds(n),stat=allocerr)

    lssit%constrin(1:m+p) = .false.
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In lssitini, allocation error.'
       stop
    end if

  end subroutine lssitini

  ! *****************************************************************
  ! *****************************************************************

  subroutine lssitend(lssit)

    implicit none

    ! SCALAR ARGUMENTS
    type(lssit_type), intent(inout) :: lssit

    ! LOCAL SCALARS
    integer :: allocerr

    deallocate(lssit%hrow,lssit%hcol,lssit%hval,&
             lssit%jsta,lssit%jlen,lssit%jvar,&
             lssit%jval,lssit%constrin,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In lssitend, allocation error.'
       stop
    end if

  end subroutine lssitend

  ! *****************************************************************
  ! *****************************************************************

  subroutine preparehv (n,nfr,frvar,lssit)
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,nfr
    type(lssit_type), intent(inout) :: lssit

    ! ARRAY ARGUMENTS
    logical, intent(in) :: frvar(n)

    ! LOCAL SCALARS
    integer :: i,j,nkeepel,oldsta,a,b

    ! LOCAL ARRAYS
    integer :: renumbered(n)

    if ( n .ne. nfr) then
      renumbered(1:n) = unpack( (/ (i,i=1,nfr) /), frvar(1:n), 0 )

      where ( renumbered(lssit%hrow(1:lssit%hne)) .ne. 0 .and. &
              renumbered(lssit%hcol(1:lssit%hne)) .ne. 0 )
          lssit%hrow(1:lssit%hne) = renumbered(lssit%hrow(1:lssit%hne))
          lssit%hcol(1:lssit%hne) = renumbered(lssit%hcol(1:lssit%hne))
      elsewhere
          lssit%hrow(1:lssit%hne) = 0
          lssit%hcol(1:lssit%hne) = 0
      end where

      nkeepel = count( lssit%hrow(1:lssit%hne) .ne. 0 )

      lssit%hval(1:nkeepel) = pack( lssit%hval(1:lssit%hne), lssit%hrow(1:lssit%hne) .ne. 0 )
      
      lssit%hrow(1:nkeepel) = pack( lssit%hrow(1:lssit%hne), lssit%hrow(1:lssit%hne) .ne. 0 )
      lssit%hcol(1:nkeepel) = pack( lssit%hcol(1:lssit%hne), lssit%hcol(1:lssit%hne) .ne. 0 )

      lssit%hne = nkeepel

      if ( .not. lssit%compressed ) then

          lssit%jne = 0
          oldsta = lssit%jsta(1)
          lssit%jsta(1) = 1
          do j = 1,lssit%m + lssit%p
            a = oldsta
            b = a + lssit%jlen(j) - 1
            if ( lssit%constrin(j) ) then
                lssit%jvar(a:b) = renumbered( lssit%jvar(a:b) )
                nkeepel = count( lssit%jvar(a:b) .ne. 0 )
                lssit%jlen(j) = nkeepel
                lssit%jval(lssit%jsta(j):lssit%jsta(j)+nkeepel-1) = pack( lssit%jval(a:b), lssit%jvar(a:b) .ne. 0 )
                lssit%jvar(lssit%jsta(j):lssit%jsta(j)+nkeepel-1) = pack( lssit%jvar(a:b), lssit%jvar(a:b) .ne. 0 )
            else
                lssit%jvar(a:b) = 0
                lssit%jlen(j) = 0
            end if
            lssit%jne = lssit%jne + lssit%jlen(j)
            if ( j .lt. lssit%m + lssit%p ) then
                oldsta = lssit%jsta(j+1)
                lssit%jsta(j+1) = lssit%jsta(j) + lssit%jlen(j)
            end if
          end do
      end if
    end if

  end subroutine preparehv 

  ! *****************************************************************
  ! *****************************************************************

  subroutine computehv (n,nfr,ifr,x,g,v,hv,ierr,lssit)
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,nfr
    integer, intent(out) :: ierr
    type(lssit_type), intent(in) :: lssit

    ! ARRAY ARGUMENTS
    real(kind=8),intent(in) :: v(nfr),x(n),g(nfr)
    real(kind=8),intent(out) :: hv(nfr)
    integer, intent(in) :: ifr(nfr)

    ierr = 0
    if ( lssit%truehp ) then
      call computetruehv(nfr,v,hv,lssit)
    else
      call computeapprohv(n,nfr,ifr,x,g,v,hv,ierr,lssit)
    end if
    if ( ierr .ne. 0 ) return 

  end subroutine computehv 

  ! *****************************************************************
  ! *****************************************************************

  subroutine computetruehv(nfr,v,hv,lssit)
    implicit none

    ! SCALAR ARGUMENTS
    integer,intent(in) :: nfr
    type(lssit_type), intent(in) :: lssit

    ! ARRAY ARGUMENTS
    real(kind=8),intent(in) :: v(nfr)
    real(kind=8),intent(out) :: hv(nfr)

    ! LOCAL SCALARS
    real(kind=8) :: atv
    integer :: i,j,k,m,p,a,b

    m = lssit%m 
    p = lssit%p

    hv(1:nfr) = 0.0d0
    do k = 1,lssit%hne
      i = lssit%hrow(k)
      j = lssit%hcol(k)

      hv(i) = hv(i) + lssit%hval(k) * v(j)
      if (j .ne. i) then
            hv(j) = hv(j) + lssit%hval(k) * v(i)
      end if 
    end do 

    if ( .not. lssit%compressed ) then 
      if ( any( lssit%constrin(1:m+p) ) ) then
        do j = 1,m + p
          if ( lssit%constrin(j) ) then
            a = lssit%jsta(j)
            b = a + lssit%jlen(j) - 1

            atv = dot_product( lssit%jval(a:b), v(lssit%jvar(a:b)) )
            hv(lssit%jvar(a:b)) = hv(lssit%jvar(a:b)) + ( atv * lssit%rho ) * lssit%jval(a:b)
          end if
        end do
      end if
    end if

  end subroutine computetruehv 

  ! *****************************************************************
  ! *****************************************************************

  subroutine computeapprohv(n,nfr,ifr,x,g,v,hv,ierr,lssit)
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,nfr
    integer, intent(out) :: ierr
    type(lssit_type), intent(in) :: lssit

    ! ARRAY ARGUMENTS
    real(kind=8),intent(in) :: v(nfr),x(n),g(nfr)
    real(kind=8),intent(out) :: hv(nfr)
    integer, intent(in) :: ifr(nfr)

    ! LOCAL SCALARS
    real(kind=8) :: t

    !LOCAL ARRAYS
    real(kind=8) :: vfull(n),hvfull(n),xaux(n),xsmall(nfr),gaux(n)

    ierr = 0
    if ( any ( lssit%constrin(1:lssit%m+lssit%p) ) ) then 
        vfull(1:n) = 0.0d0
        vfull( ifr(1:nfr) ) = v(1:nfr)
        call applyhapp(n,vfull(1:n),hvfull(1:n),lssit)
        hv(1:nfr) = hvfull(ifr(1:nfr))
    else
      xaux(1:n) = x(1:n)
      xsmall(1:nfr) = x(ifr(1:nfr))
      t = norm2( v(1:nfr) )
      if ( t .eq. 0.0d0 ) then
          write(*,*) "Failure to evaluate H*v by forward Difference"
          ierr = -11
          return
      end if
      t = macheps12 * max( maxval( abs ( xsmall(1:nfr) ) ) / maxval( abs( v(1:nfr) ) ), 1.0d0 )

      xaux(ifr(1:nfr)) = xsmall(1:nfr) + t*v(1:nfr)    
      call lssit%evalg(n,xaux,gaux,ierr,lssit%pdata)
      if ( ierr .ne. 0 ) return 
      hv(1:nfr) = (1.0d0 / t) * ( gaux(ifr(1:nfr)) - g(1:nfr) )
    end if

  end subroutine computeapprohv 

  ! *****************************************************************
  ! *****************************************************************

  subroutine comphpre(n,s,y,sty,sts,lspgmin,lspgmax,sameface,inform,lssit)
    use iso_c_binding, only: c_ptr
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)      :: n
    integer, intent(out)     :: inform 
    logical, intent(in)      :: sameface
    real(kind=8), intent(in) :: sty,sts,lspgmin,lspgmax
    type(lssit_type), intent(inout) :: lssit

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: s(n),y(n)

    ! Consider the preconditioner
    !
    !     P = Q + E + diag(rho A^t A)
    !
    ! for matrix
    !
    !     H = B + S + rho A^t A,
    !
    ! where E is the spectral correction of diag(rho A^t A) and Q is
    ! the BFGS correction of (E + diag(rho A^t A)), while S and B are
    ! the spectral and BFGS corrections of matrix (rho A^t A),
    ! respectively.
    !
    ! This subroutine computes:
    !
    ! (a) pdiag = diag(rho A^t A),
    !
    ! (b) plspg such that E = plspg I,
    !
    ! (c) psmdy = s - D^-1 y, where D = E + diag(rho A^t A), and
    !
    ! (d) the inner product psmdty = <psmdy,y>.
    !
    ! These quantities will be used later, in subroutine applyp, to
    ! compute z = P^{-1} r.

    ! LOCAL SCALARS
    integer      :: j,a,b,m,p
    real(kind=8) :: sttmp

    inform = 0

    m = lssit%m
    p = lssit%p
    lssit%s(1:n)   = s(1:n)
    lssit%y(1:n)   = y(1:n)
    lssit%sty      = sty
    lssit%seucn    = norm2( s(1:n) )
    lssit%yeucn    = norm2( y(1:n) )
    lssit%sameface = sameface
     

    ! ------------------------------------------------------------------
    ! Compute diag( rho A^t A )
    ! ------------------------------------------------------------------
    lssit%pdiag(1:n) = 0.0d0
    do j = 1,m+p
        if ( ( lssit%constrin(j) ) ) then
            a = lssit%jsta(j)
            b = a + lssit%jlen(j) - 1
            lssit%pdiag( lssit%jvar(a:b) ) = lssit%pdiag( lssit%jvar(a:b) ) + lssit%rho * lssit%jval(a:b) ** 2
        end if
    end do

    ! ------------------------------------------------------------------
    ! Compute plspg = s^t (y - diag( rho A^t A ) s) / (s^t s)
    ! ------------------------------------------------------------------

    sttmp = dot_product( lssit%pdiag(1:n), s(1:n) ** 2 )
    if ( sty - sttmp .le. 0.0d0 ) then
        lssit%plspg = 1.0d-10
    else
        lssit%plspg = max( 1.0d-10, min( ( sty - sttmp ) / sts, lspgmax ) )
    end if

    ! ------------------------------------------------------------------
    ! Compute the BFGS correction Q of ( E + diag( rho A^t A ) )
    !
    ! Q = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
    !     [ <s - D^-1 y, y> s s^t ] / (s^t y)^2,
    !
    ! where D = ( E + diag( rho A^t A ) )
    ! ------------------------------------------------------------------

    if ( sameface .and. sty .gt. macheps12 * lssit%seucn * lssit%yeucn ) then
        lssit%psmdy(1:n) = s(1:n) - y(1:n) / ( lssit%plspg  + lssit%pdiag(1:n) )
        lssit%psmdyty = dot_product( lssit%psmdy(1:n), y(1:n) ) 
    end if


  end subroutine comphpre


  ! *****************************************************************
  ! *****************************************************************

  subroutine applyhpre(n,r,z,lssit)
    use iso_c_binding, only: c_ptr
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)      :: n
    type(lssit_type), intent(in) :: lssit

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: r(n)
    real(kind=8), intent(out) :: z(n)

    ! This subroutine computes the product of the inverse of the matrix
    ! computed by subroutine comphpre times vector r, i.e., z = P^{-1} r.

    ! LOCAL SCALARS
    integer      :: m,p
    real(kind=8) :: c1,c2,psmdytr,str

    m = lssit%m
    p = lssit%p

    ! ------------------------------------------------------------------
    ! Compute ( E + diag( rho A^T A ) )^{-1} r
    ! ------------------------------------------------------------------

    z(1:n) = r(1:n) / ( lssit%plspg + lssit%pdiag(1:n) )

    ! ------------------------------------------------------------------
    ! Add Q^{-1} r, where
    !
    ! Q^{-1} = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
    !          [ <s - D^-1 y, y> s s^t ] / (s^t y)^2
    !
    ! and D = ( E + diag( rho A^T A ) )
    ! ------------------------------------------------------------------

    if ( lssit%sameface .and. lssit%sty .gt. macheps12 * lssit%seucn * lssit%yeucn ) then
        str = dot_product ( lssit%s(1:n), r(1:n) )
        psmdytr = dot_product( lssit%psmdy(1:n), r(1:n) )
        c1 = str / lssit%sty
        c2 = psmdytr / lssit%sty - lssit%psmdyty * str / lssit%sty ** 2

        z(1:n) = z(1:n) + c1 * lssit%psmdy(1:n) + c2 * lssit%s(1:n)

    end if


  end subroutine applyhpre

  ! *****************************************************************
  ! *****************************************************************

  subroutine comphapp(n,s,y,sty,sts,lspgmin,lspgmax,sameface,lssit)
    use iso_c_binding, only: c_ptr
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)      :: n
    logical, intent(in)      :: sameface
    real(kind=8), intent(in) :: sty,sts,lspgmin,lspgmax
    type(lssit_type), intent(inout) :: lssit

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: s(n),y(n)

    ! This subroutine computes an approximation H of the Hessian of the
    ! Augmented Lagrangian following a very simple idea: "discard the
    ! second order terms and then correct the remaining matrix in order
    ! to satisfy a secant equation".
    !
    ! Hence, H takes the form
    !
    ! H = B + S + rho A^t A,
    !
    ! where S is the spectral correction of (rho A^t A) and B is
    ! the BFGS correction of (S + rho A^t A). More specifically,
    !
    ! S = hlspg I,
    !
    ! where
    !
    ! hlspg = max(lspgmin, min(lspgmax, s^t (y - rho A^t A s) / s^t s)),
    !
    ! D = S + rho A^t A,
    !
    ! and
    !
    ! B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ].
    !
    ! Note that this subroutine does not compute matrix H explicitly,
    ! but computes some quantities that will be used later to compute
    ! the product of H by a vector p.
    !
    ! The quantities computed by this subroutine are:
    !
    ! (a) hlspg = s^t (y - rho A^t A s) / (s^t s)
    !
    ! (b) hds = D s = ( hlspg I + rho A^t A ) s, and
    !
    ! (c) hstds = <s,hds>.

    ! LOCAL SCALARS
    integer      :: j,a,b,m,p
    real(kind=8) :: ats

    lssit%infohapp = 0

    m = lssit%m
    p = lssit%p
    lssit%s(1:n)   = s(1:n)
    lssit%y(1:n)   = y(1:n)
    lssit%sty      = sty
    lssit%seucn    = norm2( s(1:n) )
    lssit%yeucn    = norm2( y(1:n) )
    lssit%sameface = sameface


    ! ------------------------------------------------------------------
    ! Compute hds = rho A^t A s
    ! ------------------------------------------------------------------

    lssit%hds(1:n) = 0.0d0
    do j = 1,m+p
      if ( lssit%constrin(j) ) then
        ! COMPUTE THE INNER PRODUCT <a,s>
        a = lssit%jsta(j)
        b = a + lssit%jlen(j) - 1
        ats = dot_product( lssit%jval(a:b), s( lssit%jvar(a:b) ) )
        ats = ats *  lssit%rho
        lssit%hds( lssit%jvar(a:b) ) = lssit%hds( lssit%jvar(a:b) ) + ats * lssit%jval(a:b)
      end if
    end do

    lssit%hstds = dot_product ( s(1:n), lssit%hds(1:n) )

    ! ------------------------------------------------------------------
    ! Compute hlspg = s^t (y - rho A^t A s) / (s^t s)
    ! ------------------------------------------------------------------

    if ( sty - lssit%hstds .le. 0.0d0 ) then
        lssit%hlspg = 1.0d-10
    else
        lssit%hlspg = max( 1.0d-10, min( (sty - lssit%hstds) / sts, lspgmax ) )
    end if

    lssit%hds(1:n) = lssit%hds(1:n) + lssit%hlspg * s(1:n)
    lssit%hstds = lssit%hstds + lssit%hlspg * sts

  end subroutine comphapp

  ! *****************************************************************
  ! *****************************************************************

  subroutine applyhapp(n,v,hv,lssit)
    use iso_c_binding, only: c_ptr
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)      :: n
    type(lssit_type), intent(in) :: lssit

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: v(n)
    real(kind=8), intent(out) :: hv(n)

    ! This subroutine computes the product of the matrix computed by
    ! subroutine comphapp times vector v.

    ! LOCAL SCALARS
    integer      :: a,b,m,p,j
    real(kind=8) :: atv,c1,c2,vtds,vty

    m = lssit%m
    p = lssit%p

    ! ------------------------------------------------------------------
    ! Compute ( hlspg I ) v
    ! ------------------------------------------------------------------
    hv(1:n) = lssit%hlspg * v(1:n)

    ! ------------------------------------------------------------------
    ! Add ( rho A^T A ) v
    ! ------------------------------------------------------------------
    do j = 1,m+p
        if ( lssit%constrin(j) ) then
            ! COMPUTE THE INNER PRODUCT <a,v>
            a = lssit%jsta(j)
            b = a + lssit%jlen(j) - 1
            atv = dot_product( lssit%jval(a:b), v( lssit%jvar(a:b) ) )
            atv = atv * lssit%rho

            ! ADD rho * atp * a
            hv( lssit%jvar(a:b) ) = hv( lssit%jvar(a:b) ) + atv * lssit%jval(a:b)
        end if
    end do

    ! ------------------------------------------------------------------
    ! Add B v,
    ! where B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ]
    ! ------------------------------------------------------------------

    if ( lssit%sameface .and. lssit%sty .gt. macheps12 * lssit%seucn * lssit%yeucn ) then
        vty  = dot_product( v(1:n), lssit%y(1:n) )
        vtds = dot_product( v(1:n), lssit%hds(1:n) )
        c1 = vty / lssit%sty
        c2 = vtds / lssit%hstds

        hv(1:n) = hv(1:n) + c1 * lssit%y(1:n) - c2 * lssit%hds(1:n)
    end if

  end subroutine applyhapp

end module lssit

