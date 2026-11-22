! Fase 6, passo 3: Hessiana da face na GPU via lssit (preparehv com usegpu).
!
! Parte A (correcao): para cada trial, H_ff v pela GPU (preparehv + lssit_gpu_spmv,
!   um so objeto lssit ao longo dos trials, como no Gencan) contra preparehv +
!   computehv da CPU (outro objeto, usegpu = .false.). Confere tambem quando o
!   padrao CSR e refeito.
! Parte B (medicao, D4 do docs/FASE6.md): n grande, mesma estrutura e valores
!   novos a cada "passo de Newton". Compara o tempo de preparehv com GPU (fase 6:
!   compacta + espelha + coo2csr + envio) com preparehv so CPU e com o
!   hess_csr_update da fase 5 (envio do COO + scatter em kernel), e o SpMV dos dois.
program test_hess_csr
  use iso_c_binding
  use hipfort
  use hipfort_check
  use lssit
  use lssit_gpu_mod
  use hess_csr_mod
  implicit none

  integer :: nfalhas = 0

  call parte_a()
  call parte_b()
  write(*,'(A,I0,A)') 'test_hess_csr: falhas = ', nfalhas, ' de 6'

contains

  ! ===================================================================
  subroutine parte_a()
    integer, parameter :: n = 200, hnnzmax = 2000, noff = 600, ntrial = 6
    integer :: hrow(hnnzmax), hcol(hnnzmax), hne, i, e, t, nfr, ierr
    real(8) :: hval(hnnzmax), u, err, ref
    logical :: frvar(n), frvar4(n), ok
    integer, allocatable :: ifr(:)
    real(8), allocatable, target :: v(:), hv_gpu(:)
    real(8), allocatable :: hv_cpu(:), x(:), g(:)
    type(lssit_type) :: lsg, lsc
    ! trial 3 repete a face do 2 (mesmo padrao, valores novos); trial 5: face vazia;
    ! trial 6 repete a face do 4, mas com outra estrutura de COO (padrao tem que mudar)
    real(8), parameter :: frac(ntrial) = [1.0d0, 0.7d0, 0.7d0, 0.3d0, 0.0d0, 0.3d0]
    logical, parameter :: chg_esp(ntrial) = [.true., .true., .false., .true., .true., .true.]

    write(*,'(A)') '######## parte A: preparehv GPU x CPU ########'
    call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, lsg)     ! usegpu = .true. (padrao)
    call gera_estrutura(n, noff, hne, hrow, hcol)

    do t = 1, ntrial
      if (t == 6) then
        frvar = frvar4
        call gera_estrutura(n, noff, hne, hrow, hcol)
      else if (t /= 3) then
        do i = 1, n
          call random_number(u);  frvar(i) = (u < frac(t))
        end do
      end if
      if (t == 4) frvar4 = frvar
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

      ! --- GPU: o mesmo lsg em todos os trials (o padrao pode ser reaproveitado) ---
      lsg%compressed = .true.
      lsg%hne = hne
      lsg%hrow(1:hne) = hrow(1:hne);  lsg%hcol(1:hne) = hcol(1:hne);  lsg%hval(1:hne) = hval(1:hne)
      call preparehv(n, nfr, frvar, lsg)
      hv_gpu = 0.0d0
      if (nfr > 0) then
        call hipCheck(hipMemcpy(lsg%gpu%d_x, c_loc(v), nfr * 8_c_size_t, hipMemcpyHostToDevice))
        call lssit_gpu_spmv(lsg%gpu)
        call hipCheck(hipDeviceSynchronize())
        call hipCheck(hipMemcpy(c_loc(hv_gpu), lsg%gpu%d_y, nfr * 8_c_size_t, hipMemcpyDeviceToHost))
      end if

      ! --- CPU: caminho original do lssit ---
      call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, lsc)
      lsc%usegpu = .false.
      lsc%compressed = .true.
      lsc%hne = hne
      lsc%hrow(1:hne) = hrow(1:hne);  lsc%hcol(1:hne) = hcol(1:hne);  lsc%hval(1:hne) = hval(1:hne)
      call preparehv(n, nfr, frvar, lsc)
      hv_cpu = 0.0d0
      if (nfr > 0) call computehv(n, nfr, ifr, x, g(1:nfr), v(1:nfr), hv_cpu(1:nfr), ierr, lsc)
      call lssitend(lsc)

      err = 0.0d0;  ref = 1.0d0
      if (nfr > 0) then
        err = maxval(abs(hv_gpu(1:nfr) - hv_cpu(1:nfr)))
        ref = max(1.0d0, maxval(abs(hv_cpu(1:nfr))))
      end if
      ok = (lsg%gpu%changed .eqv. chg_esp(t)) .and. lsg%gpu%active .and. err / ref <= 1.0d-14
      if (.not. ok) nfalhas = nfalhas + 1
      write(*,'(A,I0,A,I5,A,I7,A,L1,A,L1,A,ES10.2,A)') 'trial ', t, '  nfr=', nfr, '  nnz=', lsg%gpu%nnz, &
            '  padrao refeito=', lsg%gpu%changed, ' (esperado ', chg_esp(t), ')  erro rel=', err / ref, &
            merge('  OK    ', '  FALHOU', ok)

      deallocate(ifr, v, hv_gpu, hv_cpu, x, g)
    end do
    call lssitend(lsg)
  end subroutine parte_a

  ! COO "triangulo" com diagonal completa e off-diagonais aleatorias (como do CUTEst)
  subroutine gera_estrutura(n, noff, hne, hrow, hcol)
    integer, intent(in) :: n, noff
    integer, intent(out) :: hne, hrow(:), hcol(:)
    integer :: i, e
    real(8) :: u
    hne = 0
    do i = 1, n
      hne = hne + 1;  hrow(hne) = i;  hcol(hne) = i
    end do
    do e = 1, noff
      hne = hne + 1
      call random_number(u);  hrow(hne) = 1 + int(u * n)
      call random_number(u);  hcol(hne) = 1 + int(u * n)
    end do
  end subroutine gera_estrutura

  ! ===================================================================
  subroutine parte_b()
    integer, parameter :: n = 300000, noff = 900000, hnnzmax = n + noff, nrep = 20, nspmv = 50
    integer, allocatable :: hrow(:), hcol(:)
    real(8), allocatable, target :: hval(:), v(:)
    logical, allocatable :: frvar(:)
    real(8) :: u, t_new(nrep), t_cpu(nrep), t_old(nrep), s_new, s_old, d
    integer :: i, k, nfr, hne
    logical :: changed
    type(lssit_type) :: lsg, lsc
    type(hess_csr_type) :: h
    type(c_ptr) :: d_x, d_y
    integer(8) :: c0, c1, rate

    write(*,'(A)') '######## parte B: medicao (D4), n = 300000, noff = 900000, face 80% ########'
    call system_clock(count_rate=rate)
    allocate(hrow(hnnzmax), hcol(hnnzmax), hval(hnnzmax), frvar(n), v(n))
    call gera_estrutura(n, noff, hne, hrow, hcol)
    do i = 1, n
      call random_number(u);  frvar(i) = (u < 0.8d0)
    end do
    nfr = count(frvar)
    call random_number(v)

    call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, lsg)
    call lssitini(dummy_g, .true., hnnzmax, 0, n, 0, 0, lsc)
    lsc%usegpu = .false.
    call hipCheck(hipMalloc(d_x, n * 8_c_size_t))
    call hipCheck(hipMalloc(d_y, n * 8_c_size_t))
    call hess_csr_init(h, n, hnnzmax, d_x, d_y, .true.)

    do k = 1, nrep
      call random_number(hval(1:hne))          ! "evalh": mesma estrutura, valores novos

      ! fase 6: preparehv com GPU (a copia para o lssit faz o papel do evalh, fora do tempo)
      lsg%compressed = .true.;  lsg%hne = hne
      lsg%hrow(1:hne) = hrow(1:hne);  lsg%hcol(1:hne) = hcol(1:hne);  lsg%hval(1:hne) = hval(1:hne)
      call system_clock(c0)
      call preparehv(n, nfr, frvar, lsg)
      call hipCheck(hipDeviceSynchronize())
      call system_clock(c1);  t_new(k) = real(c1 - c0, 8) / real(rate, 8) * 1d3

      ! referencia: preparehv so CPU (o que o Gencan original ja paga)
      lsc%compressed = .true.;  lsc%hne = hne
      lsc%hrow(1:hne) = hrow(1:hne);  lsc%hcol(1:hne) = hcol(1:hne);  lsc%hval(1:hne) = hval(1:hne)
      call system_clock(c0)
      call preparehv(n, nfr, frvar, lsc)
      call system_clock(c1);  t_cpu(k) = real(c1 - c0, 8) / real(rate, 8) * 1d3

      ! fase 5: hess_csr_update (COO original + scatter em kernel)
      call system_clock(c0)
      call hess_csr_update(h, hne, hrow, hcol, hval, frvar, changed)
      call hipCheck(hipDeviceSynchronize())
      call system_clock(c1);  t_old(k) = real(c1 - c0, 8) / real(rate, 8) * 1d3
    end do

    ! SpMV: mesmo vetor nos dois CSR (ordens de coluna diferentes)
    call hipCheck(hipMemcpy(lsg%gpu%d_x, c_loc(v), nfr * 8_c_size_t, hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(d_x, c_loc(v), nfr * 8_c_size_t, hipMemcpyHostToDevice))
    call lssit_gpu_spmv(lsg%gpu);  call hess_csr_spmv(h);  call hipCheck(hipDeviceSynchronize())
    call system_clock(c0)
    do k = 1, nspmv
      call lssit_gpu_spmv(lsg%gpu)
    end do
    call hipCheck(hipDeviceSynchronize())
    call system_clock(c1);  s_new = real(c1 - c0, 8) / real(rate, 8) * 1d3 / nspmv
    call system_clock(c0)
    do k = 1, nspmv
      call hess_csr_spmv(h)
    end do
    call hipCheck(hipDeviceSynchronize())
    call system_clock(c1);  s_old = real(c1 - c0, 8) / real(rate, 8) * 1d3 / nspmv
    d = diff_y(lsg%gpu%d_y, d_y, nfr)

    write(*,'(A,I0,A,I0,A,I0)') 'nfr = ', nfr, '   nnz fase 6 = ', lsg%gpu%nnz, '   nnz fase 5 = ', h%nnz
    write(*,'(A)') '                               1o passo (padrao novo)   media passos 2..20 (ms)'
    write(*,'(A,F12.3,F14.3)') '  fase 6: preparehv + GPU   ', t_new(1), sum(t_new(2:)) / (nrep - 1)
    write(*,'(A,F12.3,F14.3)') '  preparehv so CPU          ', t_cpu(1), sum(t_cpu(2:)) / (nrep - 1)
    write(*,'(A,F12.3,F14.3)') '  fase 5: hess_csr_update   ', t_old(1), sum(t_old(2:)) / (nrep - 1)
    write(*,'(A,F10.4,A,F10.4,A,ES9.2)') '  SpMV (ms): fase 6 ', s_new, '   fase 5 ', s_old, &
          '   max|y6 - y5|/max|y5| = ', d

    call hess_csr_free(h)
    call hipCheck(hipFree(d_x));  call hipCheck(hipFree(d_y))
    call lssitend(lsg);  call lssitend(lsc)
  end subroutine parte_b

  real(8) function diff_y(da, db, m)
    type(c_ptr), intent(in) :: da, db
    integer, intent(in) :: m
    real(8), allocatable, target :: ya(:), yb(:)
    allocate(ya(m), yb(m))
    call hipCheck(hipMemcpy(c_loc(ya), da, m * 8_c_size_t, hipMemcpyDeviceToHost))
    call hipCheck(hipMemcpy(c_loc(yb), db, m * 8_c_size_t, hipMemcpyDeviceToHost))
    diff_y = maxval(abs(ya - yb)) / max(1.0d0, maxval(abs(yb)))
  end function diff_y

  subroutine dummy_g(n, x, g, ierr, pdata)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: pdata
    g = 0.0d0
  end subroutine dummy_g

end program test_hess_csr
