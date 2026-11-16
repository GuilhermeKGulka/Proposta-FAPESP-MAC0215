program test_hess_csr
  use iso_c_binding
  use hipfort
  use hipfort_check
  use lssit
  use hess_csr_mod
  implicit none

  integer, parameter :: n = 200, hnnzmax = 2000, noff = 600, ntrial = 5
  integer :: hrow(hnnzmax), hcol(hnnzmax), hne, i, e, t, nfr, ierr
  real(8) :: hval(hnnzmax), u, err, ref
  logical :: frvar(n), changed
  integer, allocatable :: ifr(:)
  real(8), allocatable, target :: v(:), hv_gpu(:)
  real(8), allocatable :: hv_cpu(:), x(:), g(:)
  type(c_ptr) :: d_x, d_y
  type(hess_csr_type) :: h
  type(lssit_type) :: ls
  real(8) :: frac(ntrial) = [1.0d0, 0.7d0, 0.7d0, 0.3d0, 0.0d0]   ! trial 3 repete a face do 2

  call hipCheck(hipMalloc(d_x, n * 8_c_size_t))
  call hipCheck(hipMalloc(d_y, n * 8_c_size_t))
  call hess_csr_init(h, n, hnnzmax, d_x, d_y, .true.)

  ! COO "triangulo superior" com diagonal completa (como o CUTEst devolve)
  hne = 0
  do i = 1, n
    hne = hne + 1;  hrow(hne) = i;  hcol(hne) = i
  end do
  do e = 1, noff
    hne = hne + 1
    call random_number(u);  hrow(hne) = 1 + int(u * n)
    call random_number(u);  hcol(hne) = 1 + int(u * n)
  end do

  do t = 1, ntrial
    ! face: trial 3 reaproveita a do trial 2 (so valores novos)
    if (t /= 3) then
      do i = 1, n
        call random_number(u);  frvar(i) = (u < frac(t))
      end do
    end if
    nfr = count(frvar)
    allocate(ifr(max(nfr,1)))
    if (nfr > 0) ifr(1:nfr) = pack([(i, i=1,n)], frvar)

    do e = 1, hne
      call random_number(u);  hval(e) = 2.0d0 * u - 1.0d0
    end do

    allocate(v(n), hv_gpu(n), hv_cpu(n), x(n), g(n))
    v = 0.0d0;  x = 0.0d0;  g = 0.0d0
    do i = 1, nfr
      call random_number(u);  v(i) = 2.0d0 * u - 1.0d0
    end do

    ! --- GPU ---
    call hess_csr_invalidate(h)
    call hess_csr_update(h, hne, hrow, hcol, hval, frvar, changed)
    call hipCheck(hipMemcpy(d_x, c_loc(v), n * 8_c_size_t, hipMemcpyHostToDevice))
    call hess_csr_spmv(h)
    call hipCheck(hipDeviceSynchronize())
    hv_gpu = 0.0d0
    call hipCheck(hipMemcpy(c_loc(hv_gpu), d_y, n * 8_c_size_t, hipMemcpyDeviceToHost))

    ! --- CPU: codigo real do lssit ---
    hv_cpu = 0.0d0
    call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, ls)
    ls%compressed = .true.
    ls%hne = hne
    ls%hrow(1:hne) = hrow(1:hne);  ls%hcol(1:hne) = hcol(1:hne);  ls%hval(1:hne) = hval(1:hne)
    call preparehv(n, nfr, frvar, ls)
    if (nfr > 0) call computehv(n, nfr, ifr, x, g(1:nfr), v(1:nfr), hv_cpu(1:nfr), ierr, ls)

    err = 0.0d0;  ref = 1.0d0
    if (nfr > 0) then
      err = maxval(abs(hv_gpu(1:nfr) - hv_cpu(1:nfr)))
      ref = max(1.0d0, maxval(abs(hv_cpu(1:nfr))))
    end if
    write(*,'(A,I0,A,I5,A,I7,A,L1,A,F9.3,A,ES10.2)') 'trial ', t, '  nfr=', nfr, '  nnz=', h%nnz, &
          '  padrao refeito=', changed, '  prep(ms)=', h%t_prep_ms, '  erro rel=', err / ref

    call lssitend(ls)
    deallocate(ifr, v, hv_gpu, hv_cpu, x, g)
  end do

  call hess_csr_free(h)
  call hipCheck(hipFree(d_x));  call hipCheck(hipFree(d_y))

contains

  subroutine dummy_g(n, x, g, ierr, pdata)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: pdata
    g = 0.0d0
  end subroutine dummy_g

end program test_hess_csr