! *****************************************************************
! *****************************************************************

! This is a main program that calls Gencan to solve a simple
! problem. It is intended to be used as an example of usage. Gencan
! applies to problems of the form
!
! Minimize f(x)
!
! subject to
!
!   l <= x <= u   l and u \in R^n are the bound constraints,
!
! and, through the entry point genunc, to unconstrained problems.
!
! The user provides three subroutines:
!
!   evalf  computes f(x)
!   evalg  computes the gradient of f at x
!   evalh  computes the lower triangle of the Hessian of f at x, in
!          sparse coordinate format (row, column, value).
!
! evalh is only called when truehp = .true. (Newton-CG with the true
! Hessian). With truehp = .false. Gencan uses a secant approximation
! of the Hessian-vector product and evalh is never called.
!
! The example below is the generalized Rosenbrock function
!
!   f(x) = sum_{i=1}^{n-1} [ 100 (x_{i+1} - x_i^2)^2 + (1 - x_i)^2 ],
!
! whose global minimizer is x* = (1,...,1) with f(x*) = 0. It is
! solved twice: first without bounds (genunc), then with the box
! -2 <= x_i <= 0.5, which makes the solution x* = (0.5,...) lie on
! the boundary. Change n below to test larger instances.

! *****************************************************************
! *****************************************************************

program toyprob

  use modgencan, only: gencan, genunc
  use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

  implicit none

  ! Re-define this type (pdata_type) anyway you want. Gencan receives
  ! a pointer to a 'structure of this type'. Gencan has no access to
  ! the structure. It simply passes the pointer back to the user
  ! defined subroutines evalf, evalg, and evalh. So, this is a safe
  ! way (no common blocks) of passing to the user-provided routines
  ! any information related to the problem. In this example, it is
  ! only being used for the user to count by itself the number of
  ! calls to each routine.
  type :: pdata_type
     integer :: counters(3) = 0
  end type pdata_type

  ! PROBLEM SIZE (edit here)
  integer, parameter :: n = 100

  ! LOCAL SCALARS
  logical :: extallowed,precond,truehp
  integer :: hnnzmax,ierr,istop,iter,maxit
  real(kind=8) :: eps,f,ftarget,gpsupn

  ! LOCAL ARRAYS
  logical :: lind(n),uind(n)
  real(kind=8) :: g(n),lbnd(n),ubnd(n),x(n)
  type(pdata_type), target :: pdata

  ! Maximum number of nonzeros in the lower triangle of the Hessian.
  ! For the Rosenbrock function the Hessian is tridiagonal: n
  ! diagonal entries plus n-1 subdiagonal entries.
  hnnzmax = 2 * n - 1

  ! Parameters of the solver
  eps        = 1.0d-08   ! stop when || projected gradient ||_inf <= eps
  ftarget    = -1.0d+12  ! stop if f(x) <= ftarget (unbounded problem)
  maxit      = huge(1)   ! maximum number of iterations
  extallowed = .true.    ! allow extrapolation steps
  precond    = .true.    ! use the preconditioner in the CG
  truehp     = .true.    ! true Hessian (evalh) in the Newton-CG

  ! -----------------------------------------------------------------
  ! 1) Unconstrained problem: genunc
  ! -----------------------------------------------------------------

  x(1:n) = -1.2d0
  pdata%counters(1:3) = 0

  write(*,*)
  write(*,*) '===================================================='
  write(*,*) ' Unconstrained Rosenbrock, n = ', n
  write(*,*) '===================================================='

  call genunc(evalf,evalg,evalh,truehp,hnnzmax,n,x,f,g,gpsupn, &
       ftarget,eps,maxit,extallowed,precond,iter,ierr,istop,pdata=c_loc(pdata))

  call report(n,x,f,gpsupn,iter,ierr,istop,pdata%counters)

  ! -----------------------------------------------------------------
  ! 2) Bound constrained problem: gencan
  ! -----------------------------------------------------------------

  ! lind(i) / uind(i) say whether variable i has a lower / upper bound;
  ! lbnd(i) / ubnd(i) hold the bound values (ignored when the
  ! corresponding indicator is .false.).
  lind(1:n) = .true.
  lbnd(1:n) = -2.0d0
  uind(1:n) = .true.
  ubnd(1:n) =  0.5d0

  x(1:n) = -1.2d0
  pdata%counters(1:3) = 0

  write(*,*)
  write(*,*) '===================================================='
  write(*,*) ' Bound constrained Rosenbrock, -2 <= x <= 0.5, n = ', n
  write(*,*) '===================================================='

  call gencan(evalf,evalg,evalh,truehp,hnnzmax,n,x,lind,lbnd,uind,ubnd, &
       f,g,gpsupn,ftarget,eps,maxit,extallowed,precond,iter,ierr,istop, &
       pdata=c_loc(pdata))

  call report(n,x,f,gpsupn,iter,ierr,istop,pdata%counters)

  stop

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine report(n,x,f,gpsupn,iter,ierr,istop,counters)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,iter,ierr,istop
    real(kind=8), intent(in) :: f,gpsupn

    ! ARRAY ARGUMENTS
    integer, intent(in) :: counters(3)
    real(kind=8), intent(in) :: x(n)

    write(*,*)
    write(*,*) 'Solution (first entries): ', x(1:min(n,5))
    write(*,*) 'f(x)                    = ', f
    write(*,*) 'sup-norm of proj. grad. = ', gpsupn
    write(*,*) 'Iterations              = ', iter
    write(*,*) 'Calls to evalf/g/h      = ', counters(1:3)
    write(*,*) 'ierr (0 = ok)           = ', ierr
    write(*,*) 'istop                   = ', istop

  end subroutine report

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalf(n,x,f,ierr,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    real(kind=8), intent(out) :: f
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: i
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(1) = pdata%counters(1) + 1

    ierr = 0

    f = 0.0d0
    do i = 1,n-1
       f = f + 100.0d0 * ( x(i+1) - x(i)**2 )**2 + ( 1.0d0 - x(i) )**2
    end do

  end subroutine evalf

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalg(n,x,g,ierr,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(2) = pdata%counters(2) + 1

    ierr = 0

    g(1)     = 0.0d0
    g(2:n)   = 200.0d0 * ( x(2:n) - x(1:n-1)**2 )
    g(1:n-1) = g(1:n-1) - 400.0d0 * x(1:n-1) * ( x(2:n) - x(1:n-1)**2 ) - 2.0d0 * ( 1.0d0 - x(1:n-1) )

  end subroutine evalg

  ! *****************************************************************
  ! *****************************************************************

  ! The Hessian is returned in sparse coordinate format: hnnz entries,
  ! entry k being the value hval(k) at row hrow(k) and column hcol(k).
  ! Only the lower triangle (hrow >= hcol) is given. lim is the size
  ! of the arrays; if the Hessian has more than lim nonzeros the
  ! routine must set ierr to a nonzero value and return.

  subroutine evalh(n,x,lim,hnnz,hrow,hcol,hval,ierr,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,lim
    integer, intent(out) :: hnnz
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    integer, intent(out) :: hrow(lim),hcol(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    integer :: i
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(3) = pdata%counters(3) + 1

    ierr = 0

    if ( 2 * n - 1 .gt. lim ) then
       ierr = -1
       return
    end if

    ! Diagonal: entry i gets 1200 x_i^2 - 400 x_{i+1} + 2 from the i-th term
    ! of the sum and 200 from the (i-1)-th term.
    hrow(1) = 1
    hcol(1) = 1
    hval(1) = 1200.0d0 * x(1)**2 - 400.0d0 * x(2) + 2.0d0
    do i = 2,n-1
       hrow(i) = i
       hcol(i) = i
       hval(i) = 1200.0d0 * x(i)**2 - 400.0d0 * x(i+1) + 202.0d0
    end do
    hrow(n) = n
    hcol(n) = n
    hval(n) = 200.0d0

    ! Subdiagonal: entry (i+1,i) is -400 x_i
    do i = 1,n-1
       hrow(n+i) = i + 1
       hcol(n+i) = i
       hval(n+i) = -400.0d0 * x(i)
    end do

    hnnz = 2 * n - 1

  end subroutine evalh

end program toyprob
