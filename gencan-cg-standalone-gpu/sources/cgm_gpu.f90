module cgm_gpu_mod
  use iso_c_binding
  use hipfort
  use hipfort_check
  use lssit
  use hess_csr_mod
  implicit none
  private
  public :: cgm_gpu_ctx, cgm_gpu_init, cgm_gpu_free, cgm_gpu

  integer, parameter :: MAXB = 1024
  real(kind=8), parameter :: macheps = 1.0d-16, epsnqmp = 1.0d-08, lrgstreal = huge(1.0d0)
  integer, parameter :: maxcgitnp = 5

  interface
    subroutine launch_cgm_update_p(n, coef, r, p) bind(C, name="launch_cgm_update_p")
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: n
      real(c_double), value :: coef
      type(c_ptr), value :: r, p
    end subroutine
    subroutine launch_cgm_update_p_pre(n, coef, c1, c2, r, D, u, s, p) &
        bind(C, name="launch_cgm_update_p_pre")
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: n
      real(c_double), value :: coef, c1, c2
      type(c_ptr), value :: r, D, u, s, p
    end subroutine
    subroutine launch_cgm_dots0(n, r, D, s, u, part, cnt, out) bind(C, name="launch_cgm_dots0")
      import :: c_int, c_ptr
      integer(c_int), value :: n
      type(c_ptr), value :: r, D, s, u, part, cnt, out
    end subroutine
    subroutine launch_cgm_dots1(n, delta, p, hp, d, g, r, part, cnt, out) bind(C, name="launch_cgm_dots1")
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: n
      real(c_double), value :: delta
      type(c_ptr), value :: p, hp, d, g, r, part, cnt, out
    end subroutine
    subroutine launch_cgm_update(n, coef, macheps, p, hp, d, r, xf, lbf, ubf, lindf, uindf, &
                                 dg, sv, uv, part, cnt, out) bind(C, name="launch_cgm_update")
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: n
      real(c_double), value :: coef, macheps
      type(c_ptr), value :: p, hp, d, r, xf, lbf, ubf, lindf, uindf, dg, sv, uv, part, cnt, out
    end subroutine
  end interface

  type :: cgm_gpu_ctx
    integer :: n = 0
    type(hess_csr_type) :: h
    type(c_ptr) :: d_g = c_null_ptr, d_d = c_null_ptr, d_p = c_null_ptr
    type(c_ptr) :: d_hp = c_null_ptr, d_r = c_null_ptr
    type(c_ptr) :: d_xf = c_null_ptr, d_lbf = c_null_ptr, d_ubf = c_null_ptr
    type(c_ptr) :: d_lindf = c_null_ptr, d_uindf = c_null_ptr
    type(c_ptr) :: d_dg = c_null_ptr, d_s = c_null_ptr, d_u = c_null_ptr
    type(c_ptr) :: d_part = c_null_ptr, d_cnt = c_null_ptr, d_red = c_null_ptr
    real(c_double), allocatable :: xf(:), lbf(:), ubf(:), Df(:), sf(:), uf(:)
    integer(c_int), allocatable :: lindf(:), uindf(:)
  end type cgm_gpu_ctx

contains

  subroutine cgm_gpu_init(ctx, n, hnnzmax, do_prep)
    type(cgm_gpu_ctx), intent(inout) :: ctx
    integer, intent(in) :: n, hnnzmax
    logical, intent(in) :: do_prep
    ctx%n = n
    call hipCheck(hipMalloc(ctx%d_g,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_d,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_p,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_hp, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_r,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_xf,  n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_lbf, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_ubf, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_lindf, n * 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_uindf, n * 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_dg, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_s, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_u, n * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_part, MAXB * 7 * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_cnt, 4_c_size_t))
    call hipCheck(hipMemset(ctx%d_cnt, 0, 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_red, 8 * 8_c_size_t))
    allocate(ctx%xf(n), ctx%lbf(n), ctx%ubf(n), ctx%lindf(n), ctx%uindf(n))
    allocate(ctx%Df(n), ctx%sf(n), ctx%uf(n))
    call hess_csr_init(ctx%h, n, hnnzmax, ctx%d_p, ctx%d_hp, do_prep)
  end subroutine cgm_gpu_init

  subroutine cgm_gpu_free(ctx)
    type(cgm_gpu_ctx), intent(inout) :: ctx
    call hess_csr_free(ctx%h)
    call hipCheck(hipFree(ctx%d_g));  call hipCheck(hipFree(ctx%d_d));  call hipCheck(hipFree(ctx%d_p))
    call hipCheck(hipFree(ctx%d_hp)); call hipCheck(hipFree(ctx%d_r))
    call hipCheck(hipFree(ctx%d_xf)); call hipCheck(hipFree(ctx%d_lbf)); call hipCheck(hipFree(ctx%d_ubf))
    call hipCheck(hipFree(ctx%d_lindf)); call hipCheck(hipFree(ctx%d_uindf))
    call hipCheck(hipFree(ctx%d_dg)); call hipCheck(hipFree(ctx%d_s)); call hipCheck(hipFree(ctx%d_u))
    call hipCheck(hipFree(ctx%d_part)); call hipCheck(hipFree(ctx%d_cnt)); call hipCheck(hipFree(ctx%d_red))
    deallocate(ctx%xf, ctx%lbf, ctx%ubf, ctx%lindf, ctx%uindf, ctx%Df, ctx%sf, ctx%uf)
  end subroutine cgm_gpu_free

  ! Mesma interface do cgm original (lssit no fim). Antes: hess_csr_update(ctx%h, ...) para a face atual;
  ! se precond, tambem comphpre(...) na CPU (preenche lssit%plspg, s, psmdy, ...).
  subroutine cgm_gpu(ctx, n, nfr, ifr, x, lind, lbnd, uind, ubnd, g, delta, eps, maxit, &
                     precond, d, iter, iprint, ierr, cginfo, lssit_obj)
    type(cgm_gpu_ctx), intent(inout), target :: ctx
    integer, intent(in) :: n, nfr, maxit, iprint
    integer, intent(out) :: iter, ierr, cginfo
    real(kind=8), intent(in) :: delta, eps
    logical, intent(in) :: precond
    logical, intent(in) :: lind(n), uind(n)
    integer, intent(in) :: ifr(nfr)
    real(kind=8), intent(in), target :: x(n), g(nfr), lbnd(n), ubnd(n)
    real(kind=8), intent(out), target :: d(n)
    type(lssit_type), intent(in) :: lssit_obj

    logical :: restarted, use_pre, corr
    integer :: i, j, itnqmp, itertmp, sgn, sgnprev
    real(kind=8) :: alpha, amax, abox, bbeta, coef, currprog, bestprog, dnorm2, gnorm2, gtd, gtp
    real(kind=8) :: pnorm2, ptd, pthp, ptr, ptr_raw, q, qprev, rnorm2, ztr, ztrprev
    real(kind=8) :: c1, c2, rDr, sr, ur, sty
    real(c_double), target :: red(8)
    type(c_ptr) :: pD, ps, pu

    ierr = 0
    if (nfr /= ctx%h%nfr) stop 'cgm_gpu: nfr difere da face da Hessiana (faltou hess_csr_update?)'

    do i = 1, nfr
      j = ifr(i)
      ctx%xf(i) = x(j);  ctx%lbf(i) = lbnd(j);  ctx%ubf(i) = ubnd(j)
      ctx%lindf(i) = merge(1, 0, lind(j));  ctx%uindf(i) = merge(1, 0, uind(j))
    end do
    call hipCheck(hipMemcpy(ctx%d_xf,  c_loc(ctx%xf),  nfr * 8_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(ctx%d_lbf, c_loc(ctx%lbf), nfr * 8_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(ctx%d_ubf, c_loc(ctx%ubf), nfr * 8_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(ctx%d_lindf, c_loc(ctx%lindf), nfr * 4_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(ctx%d_uindf, c_loc(ctx%uindf), nfr * 4_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(ctx%d_g, c_loc(g), nfr * 8_c_size_t, hipMemcpyHostToDevice))

    ! ---- precondicionador: dados compactados por ifr ----
    corr = .false.
    sty = 1.0d0
    if (precond) then
      if (lssit_obj%m + lssit_obj%p /= 0) stop 'cgm_gpu: precond so com m = p = 0'
      sty  = lssit_obj%sty
      corr = lssit_obj%sameface .and. &
             sty > sqrt(epsilon(1.0d0)) * lssit_obj%seucn * lssit_obj%yeucn
      do i = 1, nfr
        j = ifr(i)
        ctx%Df(i) = lssit_obj%plspg + lssit_obj%pdiag(j)
        if (corr) then
          ctx%sf(i) = lssit_obj%s(j);  ctx%uf(i) = lssit_obj%psmdy(j)
        else
          ctx%sf(i) = 0.0d0;  ctx%uf(i) = 0.0d0       ! evita 0*NaN
        end if
      end do
      call hipCheck(hipMemcpy(ctx%d_dg, c_loc(ctx%Df), nfr * 8_c_size_t, hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(ctx%d_s, c_loc(ctx%sf), nfr * 8_c_size_t, hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(ctx%d_u, c_loc(ctx%uf), nfr * 8_c_size_t, hipMemcpyHostToDevice))
    end if

    gnorm2 = dot_product(g(1:nfr), g(1:nfr))
    restarted = .false.
    itertmp = 0

    outer: do
      use_pre = precond .and. .not. restarted
      iter = 0;  itnqmp = 0
      qprev = lrgstreal;  bestprog = 0.0d0
      q = 0.0d0;  gtd = 0.0d0;  dnorm2 = 0.0d0
      rnorm2 = gnorm2;  ztr = 0.0d0;  ztrprev = 0.0d0
      c1 = 0.0d0;  c2 = 0.0d0
      sgnprev = 1
      call hipCheck(hipMemcpy(ctx%d_r, ctx%d_g, nfr * 8_c_size_t, hipMemcpyDeviceToDevice))
      call hipCheck(hipMemset(ctx%d_d, 0, nfr * 8_c_size_t))
      call hipCheck(hipMemset(ctx%d_p, 0, nfr * 8_c_size_t))

      if (use_pre) then
        pD = ctx%d_dg;  ps = ctx%d_s;  pu = ctx%d_u
        call launch_cgm_dots0(int(nfr, c_int), ctx%d_r, ctx%d_dg, ctx%d_s, ctx%d_u, &
                              ctx%d_part, ctx%d_cnt, ctx%d_red)
        call hipCheck(hipMemcpy(c_loc(red), ctx%d_red, 3 * 8_c_size_t, hipMemcpyDeviceToHost))
        rDr = red(1);  sr = red(2);  ur = red(3)
      else
        pD = c_null_ptr;  ps = c_null_ptr;  pu = c_null_ptr
      end if

      if (iprint >= 2) write(*, 984) iter, sqrt(rnorm2), sqrt(dnorm2), q

      inner: do
        if (iter /= 0 .and. ((rnorm2 <= eps * eps * gnorm2 .and. iter >= 4) .or. rnorm2 <= macheps)) then
          cginfo = 0;  exit inner
        end if
        if (iter >= max(4, maxit)) then
          cginfo = 8;  exit inner
        end if

        ! ---- z = M^{-1} r  (so os escalares; z em si so aparece fundido em p) ----
        ztrprev = ztr
        if (use_pre) then
          if (corr) then
            c1 = sr / sty
            c2 = ur / sty - lssit_obj%psmdyty * sr / sty**2
          end if
          ztr = rDr + c1 * ur + c2 * sr
        else
          ztr = rnorm2
        end if

        ! ---- direcao ----
        if (iter == 0) then
          coef = 0.0d0
        else
          bbeta = ztr / ztrprev
          coef = bbeta * real(sgnprev, 8)
        end if
        if (use_pre) then
          call launch_cgm_update_p_pre(int(nfr, c_int), coef, c1, c2, ctx%d_r, ctx%d_dg, ctx%d_u, &
                                       ctx%d_s, ctx%d_p)
        else
          call launch_cgm_update_p(int(nfr, c_int), coef, ctx%d_r, ctx%d_p)
        end if

        call hess_csr_spmv(ctx%h)
        call launch_cgm_dots1(int(nfr, c_int), delta, ctx%d_p, ctx%d_hp, ctx%d_d, ctx%d_g, ctx%d_r, &
                              ctx%d_part, ctx%d_cnt, ctx%d_red)
        call hipCheck(hipMemcpy(c_loc(red), ctx%d_red, 7 * 8_c_size_t, hipMemcpyDeviceToHost))

        ptr_raw = red(4)
        pnorm2  = red(5)
        if (ptr_raw > 0.0d0) then
          sgn = -1
        else
          sgn = 1
        end if
        ptr  = -abs(ptr_raw)
        pthp = red(1)
        ptd  = real(sgn, 8) * red(2)
        gtp  = real(sgn, 8) * red(3)
        if (sgn > 0) then
          amax = red(6)
        else
          amax = red(7)
        end if

        if (pthp > 0.0d0) then
          alpha = min(amax, ztr / pthp)
        else if (iter == 0) then
          alpha = amax
        else
          cginfo = 7;  exit inner
        end if

        gtd = gtd + alpha * gtp
        dnorm2 = dnorm2 + alpha * alpha * pnorm2 + 2.0d0 * alpha * ptd

        if (gtd > 0.0d0 .or. (gtd * gtd < 1.0d-12 * gnorm2 * dnorm2)) then
          if (precond .and. .not. restarted .and. iter == 0) then
            restarted = .true.
            itertmp = iter
            cycle outer
          end if
          cginfo = 3;  exit inner
        end if

        qprev = q
        q = q + 0.5d0 * alpha * alpha * pthp + alpha * ptr

        coef = real(sgn, 8) * alpha
        call launch_cgm_update(int(nfr, c_int), coef, macheps, ctx%d_p, ctx%d_hp, ctx%d_d, ctx%d_r, &
                               ctx%d_xf, ctx%d_lbf, ctx%d_ubf, ctx%d_lindf, ctx%d_uindf, &
                               pD, ps, pu, ctx%d_part, ctx%d_cnt, ctx%d_red)
        call hipCheck(hipMemcpy(c_loc(red), ctx%d_red, 6 * 8_c_size_t, hipMemcpyDeviceToHost))
        rnorm2 = red(1)
        if (use_pre) then
          rDr = red(3);  sr = red(4);  ur = red(5)
        end if
        abox = red(6)

        iter = iter + 1
        if (iprint >= 2) write(*, 984) iter, sqrt(rnorm2), sqrt(dnorm2), q

        if (alpha == amax) then
          cginfo = 1;  exit inner
        end if
        if (abox <= 0.1d0) then
          cginfo = 5;  exit inner
        end if
        if (red(2) == 0.0d0) then
          cginfo = 6;  exit inner
        end if

        currprog = qprev - q
        bestprog = max(currprog, bestprog)
        if (currprog <= epsnqmp * bestprog) then
          itnqmp = itnqmp + 1
          if (itnqmp >= maxcgitnp) then
            cginfo = 4;  exit inner
          end if
        else
          itnqmp = 0
        end if

        sgnprev = sgn
      end do inner
      exit outer
    end do outer

    iter = iter + itertmp
    if (iprint >= 1) write(*,*) 'Exit cgm_gpu with cginfo = ', cginfo

    call hipCheck(hipMemcpy(c_loc(d), ctx%d_d, nfr * 8_c_size_t, hipMemcpyDeviceToHost))

984 format(5X, 'CG iter = ', I7, ' rnorm = ', 1P, D10.4, ' dnorm = ', 1P, D10.4, ' q = ', 1P, D11.4)
  end subroutine cgm_gpu

end module cgm_gpu_mod