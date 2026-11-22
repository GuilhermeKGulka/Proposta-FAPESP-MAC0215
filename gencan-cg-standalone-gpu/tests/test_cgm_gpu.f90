! Fase 6, passo 4: o cgm com a Hessiana da face na GPU (desvio para cgm_gpu)
! contra o cgm original, em cenarios sinteticos. Os dois lados chamam a MESMA
! rotina cgm do Gencan; o que muda e o lssit: lsc com usegpu = .false. (caminho
! original da CPU) e lsg com usegpu = .true. (preparehv manda a Hessiana para a
! GPU e o cgm desvia para cgm_gpu). Criterio: mesmo iter, mesmo cginfo e
! ||d_gpu - d_cpu|| / ||d_cpu|| <= 1e-12.
program test_cgm_gpu
  use iso_c_binding
  use lssit
  use modcgm, only: cgm
  implicit none

  integer, parameter :: n = 3000, noff = 9000, hnnzmax = n + noff, nscen = 9
  character(len=34), parameter :: nome(nscen) = [character(len=34) :: &
       'SPD, sem bounds, delta grande', 'SPD, delta pequeno', 'SPD, bounds distantes', &
       'SPD, bounds proximos', 'indefinida', 'SPD, maxit = 5', &
       'SPD, precond (com correcao)', 'SPD, bounds prox., precond', 'SPD, precond (sem correcao)']
  logical, parameter :: spd(nscen)   = [.true., .true., .true., .true., .false., .true., .true., .true., .true.]
  logical, parameter :: hasb(nscen)  = [.false., .false., .true., .true., .false., .false., .false., .true., .false.]
  real(8), parameter :: delta(nscen) = [1d8, 1d-2, 1d8, 1d8, 10d0, 1d8, 1d8, 1d8, 1d8]
  real(8), parameter :: dist(nscen)  = [0d0, 0d0, 1d3, 5d-2, 0d0, 0d0, 0d0, 5d-2, 0d0]
  real(8), parameter :: frac(nscen)  = [1d0, 1d0, 0.8d0, 0.8d0, 1d0, 1d0, 0.9d0, 0.8d0, 0.9d0]
  integer, parameter :: maxit(nscen) = [n, n, n, n, n, 5, n, n, n]
  ! cenarios 7-9: precondicionador; sface = argumento sameface do comphpre (liga a correcao BFGS)
  logical, parameter :: prec(nscen)  = [.false., .false., .false., .false., .false., .false., .true., .true., .true.]
  logical, parameter :: sface(nscen) = [.false., .false., .false., .false., .false., .false., .true., .true., .false.]

  integer :: hrow(hnnzmax), hcol(hnnzmax), hne, i, e, s, nfr, itc, itg, infoc, infog, ierr
  real(8) :: hval(hnnzmax), rowsum(n), u, derr, t0c, t0g
  real(8) :: sv(n), yv(n), sty, sts
  integer :: inform, nfalhas
  logical :: ok
  integer(8) :: c0, c1, rate
  logical :: frvar(n), lind(n), uind(n)
  real(8), allocatable :: x(:), lbnd(:), ubnd(:), gfull(:), gc(:), dcpu(:), dgpu(:)
  integer, allocatable :: ifr(:)
  type(lssit_type) :: lsc, lsg

  call system_clock(count_rate=rate)
  nfalhas = 0
  allocate(x(n), lbnd(n), ubnd(n), gfull(n), dcpu(n), dgpu(n))
  call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, lsg)   ! GPU: um objeto so, como no Gencan

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

    if (prec(s)) then                       ! dados do comphpre: y = H s (H simetrica do triangulo)
      do i = 1, n
        call random_number(u);  sv(i) = 2.0d0 * u - 1.0d0
      end do
      yv = 0.0d0
      do e = 1, hne
        yv(hrow(e)) = yv(hrow(e)) + hval(e) * sv(hcol(e))
        if (hrow(e) /= hcol(e)) yv(hcol(e)) = yv(hcol(e)) + hval(e) * sv(hrow(e))
      end do
      sty = dot_product(sv, yv);  sts = dot_product(sv, sv)
    end if

    ! ---------- CPU: cgm original (usegpu = .false.) ----------
    call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, lsc)
    lsc%usegpu = .false.
    call carrega(lsc)
    dcpu = 0.0d0
    call system_clock(c0)
    call cgm(n, nfr, ifr, x, lind, lbnd, uind, ubnd, gc, delta(s), 1.0d-8, maxit(s), prec(s), &
             dcpu, itc, 0, ierr, infoc, lsc)
    call system_clock(c1);  t0c = real(c1 - c0, 8) / real(rate, 8) * 1d3
    if (lsc%gpu%active) stop 'test_cgm_gpu: lsc nao deveria usar a GPU'
    call lssitend(lsc)

    ! ---------- GPU: o mesmo cgm, com a Hessiana da face na GPU ----------
    call carrega(lsg)
    dgpu = 0.0d0
    call system_clock(c0)
    call cgm(n, nfr, ifr, x, lind, lbnd, uind, ubnd, gc, delta(s), 1.0d-8, maxit(s), prec(s), &
             dgpu, itg, 0, ierr, infog, lsg)
    call system_clock(c1);  t0g = real(c1 - c0, 8) / real(rate, 8) * 1d3
    if (.not. lsg%gpu%active) stop 'test_cgm_gpu: lsg deveria usar a GPU'

    derr = norm2(dgpu(1:nfr) - dcpu(1:nfr)) / max(norm2(dcpu(1:nfr)), tiny(1.0d0))
    write(*,'(A,A,A,I5)') '== ', nome(s), '  nfr=', nfr
    write(*,'(A,I6,A,I2,A,F10.3,A)') '   CPU: iter=', itc, ' cginfo=', infoc, '  t(ms)=', t0c
    write(*,'(A,I6,A,I2,A,F10.3,A,ES9.2)') '   GPU: iter=', itg, ' cginfo=', infog, '  t(ms)=', t0g, &
          '   ||d_gpu-d_cpu||/||d_cpu|| = ', derr
    ok = (itc == itg) .and. (infoc == infog) .and. (derr <= 1.0d-12)
    if (.not. ok) nfalhas = nfalhas + 1
    write(*,'(A,A)') '   ', merge('OK     ', 'FALHOU ', ok)

    deallocate(ifr, gc)
  end do

  call lssitend(lsg)
  write(*,'(A,I0,A,I0)') 'test_cgm_gpu: falhas = ', nfalhas, ' de ', nscen

contains
  ! como no gencan: evalh (copia do COO) -> comphpre -> preparehv
  subroutine carrega(ls)
    type(lssit_type), intent(inout) :: ls
    ls%compressed = .true.
    ls%hne = hne
    ls%hrow(1:hne) = hrow(1:hne);  ls%hcol(1:hne) = hcol(1:hne);  ls%hval(1:hne) = hval(1:hne)
    if (prec(s)) call comphpre(n, sv, yv, sty, sts, 1.0d-16, 1.0d16, sface(s), inform, ls)
    call preparehv(n, nfr, frvar, ls)
  end subroutine carrega

  subroutine dummy_g(n, x, g, ierr, pdata)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: pdata
    g = 0.0d0
  end subroutine dummy_g
end program test_cgm_gpu
