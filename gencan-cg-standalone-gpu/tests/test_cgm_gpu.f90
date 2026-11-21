program test_cgm_gpu
  use iso_c_binding
  use hipfort
  use hipfort_check
  use lssit
  use modcgm, only: cgm
  use hess_csr_mod
  use cgm_gpu_mod
  implicit none

  integer, parameter :: n = 3000, noff = 9000, hnnzmax = n + noff, nscen = 6
  character(len=34), parameter :: nome(nscen) = [character(len=34) :: &
       'SPD, sem bounds, delta grande', 'SPD, delta pequeno', 'SPD, bounds distantes', &
       'SPD, bounds proximos', 'indefinida', 'SPD, maxit = 5']
  logical, parameter :: spd(nscen)   = [.true., .true., .true., .true., .false., .true.]
  logical, parameter :: hasb(nscen)  = [.false., .false., .true., .true., .false., .false.]
  real(8), parameter :: delta(nscen) = [1d8, 1d-2, 1d8, 1d8, 10d0, 1d8]
  real(8), parameter :: dist(nscen)  = [0d0, 0d0, 1d3, 5d-2, 0d0, 0d0]
  real(8), parameter :: frac(nscen)  = [1d0, 1d0, 0.8d0, 0.8d0, 1d0, 1d0]
  integer, parameter :: maxit(nscen) = [n, n, n, n, n, 5]

  integer :: hrow(hnnzmax), hcol(hnnzmax), hne, i, e, s, nfr, itc, itg, infoc, infog, ierr
  real(8) :: hval(hnnzmax), rowsum(n), u, derr, t0c, t1c, t0g, t1g
  integer(8) :: c0, c1, rate
  logical :: frvar(n), lind(n), uind(n), changed
  real(8), allocatable :: x(:), lbnd(:), ubnd(:), gfull(:), gc(:), dcpu(:), dgpu(:)
  integer, allocatable :: ifr(:)
  type(lssit_type) :: ls
  type(cgm_gpu_ctx) :: ctx

  call system_clock(count_rate=rate)
  call cgm_gpu_init(ctx, n, hnnzmax, .true.)
  allocate(x(n), lbnd(n), ubnd(n), gfull(n), dcpu(n), dgpu(n))

  do s = 1, nscen
    ! ---------- H: diagonal + off-diagonais aleatorias (COO "como vem do CUTEst") ----------
    hne = 0
    rowsum = 0.0d0
    do e = 1, noff
      call random_number(u);  hrow(n + e) = 1 + int(u * n)
      call random_number(u);  hcol(n + e) = 1 + int(u * n)
      if (hrow(n + e) == hcol(n + e)) hcol(n + e) = mod(hcol(n + e), n) + 1
      call random_number(u);  hval(n + e) = 2.0d0 * u - 1.0d0
      rowsum(hrow(n + e)) = rowsum(hrow(n + e)) + abs(hval(n + e))
      rowsum(hcol(n + e)) = rowsum(hcol(n + e)) + abs(hval(n + e))
    end do
    do i = 1, n
      hrow(i) = i;  hcol(i) = i
      if (spd(s)) then
        hval(i) = 1.0d0 + rowsum(i)
      else
        call random_number(u);  hval(i) = 4.0d0 * u - 2.0d0
      end if
    end do
    hne = hnnzmax

    ! ---------- ponto, bounds, face, gradiente ----------
    do i = 1, n
      call random_number(u);  x(i) = 2.0d0 * u - 1.0d0
      call random_number(u);  frvar(i) = (u < frac(s))
      call random_number(u);  gfull(i) = 2.0d0 * u - 1.0d0
      lind(i) = hasb(s);  uind(i) = hasb(s)
      call random_number(u);  lbnd(i) = x(i) - dist(s) * (1.0d0 + u)
      call random_number(u);  ubnd(i) = x(i) + dist(s) * (1.0d0 + u)
    end do
    nfr = count(frvar)
    allocate(ifr(nfr), gc(nfr))
    ifr = pack([(i, i = 1, n)], frvar)
    gc = gfull(ifr)

    ! ---------- CPU: cgm ORIGINAL ----------
    call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, ls)
    ls%compressed = .true.
    ls%hne = hne
    ls%hrow(1:hne) = hrow(1:hne);  ls%hcol(1:hne) = hcol(1:hne);  ls%hval(1:hne) = hval(1:hne)
    call preparehv(n, nfr, frvar, ls)
    dcpu = 0.0d0
    call system_clock(c0)
    call cgm(n, nfr, ifr, x, lind, lbnd, uind, ubnd, gc, delta(s), 1.0d-8, maxit(s), .false., &
             dcpu, itc, 0, ierr, infoc, ls)
    call system_clock(c1);  t0c = real(c1 - c0, 8) / real(rate, 8) * 1d3

    ! ---------- GPU ----------
    dgpu = 0.0d0
    call hess_csr_invalidate(ctx%h)
    call hess_csr_update(ctx%h, hne, hrow, hcol, hval, frvar, changed)
    call system_clock(c0)
    call cgm_gpu(ctx, n, nfr, ifr, x, lind, lbnd, uind, ubnd, gc, delta(s), 1.0d-8, maxit(s), .false., &
                 dgpu, itg, 0, ierr, infog, ls)
    call system_clock(c1);  t0g = real(c1 - c0, 8) / real(rate, 8) * 1d3

    derr = norm2(dgpu(1:nfr) - dcpu(1:nfr)) / max(norm2(dcpu(1:nfr)), tiny(1.0d0))
    write(*,'(A,A,A,I5)') '== ', nome(s), '  nfr=', nfr
    write(*,'(A,I6,A,I2,A,F10.3,A)') '   CPU: iter=', itc, ' cginfo=', infoc, '  t(ms)=', t0c
    write(*,'(A,I6,A,I2,A,F10.3,A,ES9.2)') '   GPU: iter=', itg, ' cginfo=', infog, '  t(ms)=', t0g, &
          '   ||d_gpu-d_cpu||/||d_cpu|| = ', derr

    call lssitend(ls)
    deallocate(ifr, gc)
  end do

  call cgm_gpu_free(ctx)

contains
  subroutine dummy_g(n, x, g, ierr, pdata)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: pdata
    g = 0.0d0
  end subroutine dummy_g
end program test_cgm_gpu