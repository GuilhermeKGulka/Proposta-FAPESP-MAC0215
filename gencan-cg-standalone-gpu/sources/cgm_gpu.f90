module cgm_gpu_mod
  use iso_c_binding
  use hipfort
  use hipfort_check
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
    subroutine launch_cgm_dots1(n, delta, p, hp, d, g, r, part, cnt, out) bind(C, name="launch_cgm_dots1")
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: n
      real(c_double), value :: delta
      type(c_ptr), value :: p, hp, d, g, r, part, cnt, out
    end subroutine
    subroutine launch_cgm_update(n, coef, macheps, p, hp, d, r, xf, lbf, ubf, lindf, uindf, &
                                 part, cnt, out) bind(C, name="launch_cgm_update")
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: n
      real(c_double), value :: coef, macheps
      type(c_ptr), value :: p, hp, d, r, xf, lbf, ubf, lindf, uindf, part, cnt, out
    end subroutine
  end interface

  type :: cgm_gpu_ctx
    integer :: n = 0
    type(hess_csr_type) :: h                 ! Hessiana CSR (nfr x nfr) + SpMV
    type(c_ptr) :: d_g = c_null_ptr, d_d = c_null_ptr, d_p = c_null_ptr
    type(c_ptr) :: d_hp = c_null_ptr, d_r = c_null_ptr
    type(c_ptr) :: d_xf = c_null_ptr, d_lbf = c_null_ptr, d_ubf = c_null_ptr
    type(c_ptr) :: d_lindf = c_null_ptr, d_uindf = c_null_ptr
    type(c_ptr) :: d_part = c_null_ptr, d_cnt = c_null_ptr, d_red = c_null_ptr
    real(c_double), allocatable :: xf(:), lbf(:), ubf(:)
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
    call hipCheck(hipMalloc(ctx%d_part, MAXB * 7 * 8_c_size_t))
    call hipCheck(hipMalloc(ctx%d_cnt, 4_c_size_t))
    call hipCheck(hipMemset(ctx%d_cnt, 0, 4_c_size_t))
    call hipCheck(hipMalloc(ctx%d_red, 8 * 8_c_size_t))
    allocate(ctx%xf(n), ctx%lbf(n), ctx%ubf(n), ctx%lindf(n), ctx%uindf(n))
    ! SpMV: x = p, y = Hp
    call hess_csr_init(ctx%h, n, hnnzmax, ctx%d_p, ctx%d_hp, do_prep)
  end subroutine cgm_gpu_init

  subroutine cgm_gpu_free(ctx)
    type(cgm_gpu_ctx), intent(inout) :: ctx
    call hess_csr_free(ctx%h)
    call hipCheck(hipFree(ctx%d_g));  call hipCheck(hipFree(ctx%d_d));  call hipCheck(hipFree(ctx%d_p))
    call hipCheck(hipFree(ctx%d_hp)); call hipCheck(hipFree(ctx%d_r))
    call hipCheck(hipFree(ctx%d_xf)); call hipCheck(hipFree(ctx%d_lbf)); call hipCheck(hipFree(ctx%d_ubf))
    call hipCheck(hipFree(ctx%d_lindf)); call hipCheck(hipFree(ctx%d_uindf))
    call hipCheck(hipFree(ctx%d_part)); call hipCheck(hipFree(ctx%d_cnt)); call hipCheck(hipFree(ctx%d_red))
    deallocate(ctx%xf, ctx%lbf, ctx%ubf, ctx%lindf, ctx%uindf)
  end subroutine cgm_gpu_free

  ! Mesma interface do cgm original, trocando lssit por ctx. Antes de chamar e preciso ter
  ! feito hess_csr_update(ctx%h, ...) para a face e os valores atuais (substitui preparehv).
  subroutine cgm_gpu(ctx, n, nfr, ifr, x, lind, lbnd, uind, ubnd, g, delta, eps, maxit, &
                     precond, d, iter, iprint, ierr, cginfo)
    type(cgm_gpu_ctx), intent(inout), target :: ctx
    integer, intent(in) :: n, nfr, maxit, iprint
    integer, intent(out) :: iter, ierr, cginfo
    real(kind=8), intent(in) :: delta, eps
    logical, intent(in) :: precond
    logical, intent(in) :: lind(n), uind(n)
    integer, intent(in) :: ifr(nfr)
    real(kind=8), intent(in), target :: x(n), g(nfr), lbnd(n), ubnd(n)
    real(kind=8), intent(out), target :: d(n)

    logical :: restarted
    integer :: i, j, itnqmp, itertmp, sgn, sgnprev
    real(kind=8) :: alpha, amax, abox, bbeta, coef, currprog, bestprog, dnorm2, gnorm2, gtd, gtp
    real(kind=8) :: pnorm2, ptd, pthp, ptr, ptr_raw, q, qprev, rnorm2, ztr, ztrprev
    real(c_double), target :: red(8)

    ierr = 0
    if (precond) stop 'cgm_gpu: precondicionador ainda nao implementado (etapa 4)'
    if (nfr /= ctx%h%nfr) stop 'cgm_gpu: nfr difere da face da Hessiana (faltou hess_csr_update?)'

    ! --- dados por-chamada: compactados por ifr (host) e enviados uma vez ---
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

    gnorm2 = dot_product(g(1:nfr), g(1:nfr))

    restarted = .false.
    itertmp = 0

    outer: do        ! "001 continue" do original: reinicio sem precondicionador
      iter = 0;  itnqmp = 0
      qprev = lrgstreal;  bestprog = 0.0d0
      q = 0.0d0;  gtd = 0.0d0;  dnorm2 = 0.0d0
      rnorm2 = gnorm2;  ztr = 0.0d0;  ztrprev = 0.0d0
      sgnprev = 1
      call hipCheck(hipMemcpy(ctx%d_r, ctx%d_g, nfr * 8_c_size_t, hipMemcpyDeviceToDevice))
      call hipCheck(hipMemset(ctx%d_d, 0, nfr * 8_c_size_t))
      call hipCheck(hipMemset(ctx%d_p, 0, nfr * 8_c_size_t))

      if (iprint >= 2) write(*, 984) iter, sqrt(rnorm2), sqrt(dnorm2), q

      inner: do
        ! ---------------- criterios de parada ----------------
        if (iter /= 0 .and. ((rnorm2 <= eps * eps * gnorm2 .and. iter >= 4) .or. rnorm2 <= macheps)) then
          cginfo = 0;  exit inner
        end if
        if (iter >= max(4, maxit)) then
          cginfo = 8;  exit inner
        end if

        ! ---------------- z (precondicionador) ----------------
        if (.not. precond .or. restarted) then
          ztrprev = ztr
          ztr = rnorm2
        else
          stop 'cgm_gpu: ramo com precondicionador (etapa 4)'
        end if

        ! ---------------- direcao: p = -z + beta p ----------------
        if (iter == 0) then
          coef = 0.0d0
        else
          bbeta = ztr / ztrprev
          coef = bbeta * real(sgnprev, 8)      ! p anterior foi (talvez) invertido
        end if
        call launch_cgm_update_p(int(nfr, c_int), coef, ctx%d_r, ctx%d_p)

        ! ---------------- hp = H p  e reducoes ----------------
        call hess_csr_spmv(ctx%h)
        call launch_cgm_dots1(int(nfr, c_int), delta, ctx%d_p, ctx%d_hp, ctx%d_d, ctx%d_g, ctx%d_r, &
                              ctx%d_part, ctx%d_cnt, ctx%d_red)
        call hipCheck(hipMemcpy(c_loc(red), ctx%d_red, 7 * 8_c_size_t, hipMemcpyDeviceToHost))

        ptr_raw = red(4)
        pnorm2  = red(5)
        if (ptr_raw > 0.0d0) then          ! original: p = -p, ptr = -ptr
          sgn = -1
        else
          sgn = 1
        end if
        ptr  = -abs(ptr_raw)
        pthp = red(1)                      ! invariante ao sinal
        ptd  = real(sgn, 8) * red(2)
        gtp  = real(sgn, 8) * red(3)
        if (sgn > 0) then
          amax = red(6)
        else
          amax = red(7)
        end if

        ! ---------------- passo ----------------
        if (pthp > 0.0d0) then
          alpha = min(amax, ztr / pthp)
        else if (iter == 0) then
          alpha = amax
        else
          cginfo = 7;  exit inner
        end if

        ! ---------------- condicao de angulo ----------------
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

        ! ---------------- modelo quadratico ----------------
        qprev = q
        q = q + 0.5d0 * alpha * alpha * pthp + alpha * ptr

        ! ---------------- d, r, e testes que dependem do novo d ----------------
        coef = real(sgn, 8) * alpha
        call launch_cgm_update(int(nfr, c_int), coef, macheps, ctx%d_p, ctx%d_hp, ctx%d_d, ctx%d_r, &
                               ctx%d_xf, ctx%d_lbf, ctx%d_ubf, ctx%d_lindf, ctx%d_uindf, &
                               ctx%d_part, ctx%d_cnt, ctx%d_red)
        call hipCheck(hipMemcpy(c_loc(red), ctx%d_red, 3 * 8_c_size_t, hipMemcpyDeviceToHost))
        rnorm2 = red(1)
        abox   = red(3)

        iter = iter + 1

        if (iprint >= 2) write(*, 984) iter, sqrt(rnorm2), sqrt(dnorm2), q

        if (alpha == amax) then                 ! fronteira da regiao de confianca
          cginfo = 1;  exit inner
        end if
        if (abox <= 0.1d0) then                 ! pouca proporcao util
          cginfo = 5;  exit inner
        end if
        if (red(2) == 0.0d0) then               ! samep: iterados muito proximos
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