module hess_csr_mod
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipsparse_generic
  implicit none
  private
  public :: hess_csr_type, hess_csr_init, hess_csr_update, hess_csr_spmv, hess_csr_free, hess_csr_invalidate

  real(c_double), target, save :: one = 1.0d0, zero = 0.0d0

  interface
    subroutine launch_scatter(hne, v, p1, p2, val) bind(C, name="launch_scatter")
      import :: c_int, c_ptr
      integer(c_int), value :: hne
      type(c_ptr), value :: v, p1, p2, val
    end subroutine
  end interface

  type :: hess_csr_type
    integer :: n = 0, hnnzmax = 0, nfr = 0, nnz = 0, hne = 0
    logical :: has_face = .false., empty = .true., spmv_ready = .false., do_prep = .true.
    real(8) :: t_prep_ms = 0.0d0

    ! host
    logical, allocatable :: frvar_prev(:)
    integer, allocatable :: renum(:), rcount(:), cur(:)
    integer, allocatable :: pos1(:), pos2(:)   ! target vem do dummy h em hess_csr_update
    integer(c_int), allocatable :: row_ptr(:), col(:)
    real(c_double), allocatable :: val(:)

    ! device (d_x / d_y pertencem ao chamador; capacidade n)
    type(c_ptr) :: handle   = c_null_ptr
    type(c_ptr) :: d_row_ptr = c_null_ptr, d_col = c_null_ptr, d_val = c_null_ptr
    type(c_ptr) :: d_x = c_null_ptr, d_y = c_null_ptr
    type(c_ptr) :: d_hcoo = c_null_ptr, d_pos1 = c_null_ptr, d_pos2 = c_null_ptr
    type(c_ptr) :: mat = c_null_ptr, vx = c_null_ptr, vy = c_null_ptr, d_buf = c_null_ptr
  end type hess_csr_type

contains

  subroutine hess_csr_init(h, n, hnnzmax, d_x, d_y, do_prep)
    type(hess_csr_type), intent(inout) :: h
    integer, intent(in) :: n, hnnzmax
    type(c_ptr), intent(in) :: d_x, d_y
    logical, intent(in) :: do_prep

    h%n = n;  h%hnnzmax = hnnzmax;  h%do_prep = do_prep
    h%d_x = d_x;  h%d_y = d_y
    allocate(h%frvar_prev(n), h%renum(n), h%rcount(n), h%cur(n))
    allocate(h%pos1(hnnzmax), h%pos2(hnnzmax))
    allocate(h%row_ptr(n+1), h%col(2*hnnzmax), h%val(2*hnnzmax))
    h%frvar_prev = .false.

    call hipCheck(hipMalloc(h%d_row_ptr, (n+1) * 4_c_size_t))
    call hipCheck(hipMalloc(h%d_col, 2 * hnnzmax * 4_c_size_t))
    call hipCheck(hipMalloc(h%d_val, 2 * hnnzmax * 8_c_size_t))
    call hipCheck(hipMalloc(h%d_hcoo, hnnzmax * 8_c_size_t))
    call hipCheck(hipMalloc(h%d_pos1, hnnzmax * 4_c_size_t))
    call hipCheck(hipMalloc(h%d_pos2, hnnzmax * 4_c_size_t))
    call hipsparseCheck(hipsparseCreate(h%handle))
  end subroutine hess_csr_init

  ! COO (indices em 1..n, como vem de evalh) + frvar  ->  CSR nfr x nfr na GPU.
  ! changed = .true. se o padrao foi refeito (face nova).
  subroutine hess_csr_update(h, hne, hrow, hcol, hval, frvar, changed)
    type(hess_csr_type), intent(inout), target :: h
    integer, intent(in) :: hne
    integer, intent(in) :: hrow(:), hcol(:)
    real(8), intent(in), target, contiguous :: hval(:)
    logical, intent(in) :: frvar(:)
    logical, intent(out) :: changed

    if (hne > h%hnnzmax) stop 'hess_csr_update: hne > hnnzmax'

    changed = .not. h%has_face .or. hne /= h%hne .or. any(frvar(1:h%n) .neqv. h%frvar_prev(1:h%n))

    if (changed) then
      call teardown_spmv(h)
      call build_pattern(h, hne, hrow, hcol, frvar)
    end if

    if (changed .and. h%nfr > 0) then
      call hipCheck(hipMemcpy(h%d_row_ptr, c_loc(h%row_ptr), (h%nfr+1) * 4_c_size_t, hipMemcpyHostToDevice))
      if (h%nnz > 0) &
        call hipCheck(hipMemcpy(h%d_col, c_loc(h%col), h%nnz * 4_c_size_t, hipMemcpyHostToDevice))
    end if
    if (changed .and. h%nfr > 0 .and. hne > 0) then      ! posicoes COO -> CSR (so quando o padrao muda)
      call hipCheck(hipMemcpy(h%d_pos1, c_loc(h%pos1), hne * 4_c_size_t, hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(h%d_pos2, c_loc(h%pos2), hne * 4_c_size_t, hipMemcpyHostToDevice))
    end if
    if (h%nnz > 0 .and. hne > 0) then                    ! valores: COO -> GPU, scatter em kernel
      call hipCheck(hipMemcpy(h%d_hcoo, c_loc(hval), hne * 8_c_size_t, hipMemcpyHostToDevice))
      call launch_scatter(int(hne, c_int), h%d_hcoo, h%d_pos1, h%d_pos2, h%d_val)
    end if

    if (changed) call setup_spmv(h)
  end subroutine hess_csr_update

  subroutine build_pattern(h, hne, hrow, hcol, frvar)
    type(hess_csr_type), intent(inout) :: h
    integer, intent(in) :: hne
    integer, intent(in) :: hrow(:), hcol(:)
    logical, intent(in) :: frvar(:)
    integer :: i, e, r, c, k

    ! renumeracao identica a de preparehv: i-esimo .true. de frvar -> i
    k = 0
    do i = 1, h%n
      if (frvar(i)) then
        k = k + 1
        h%renum(i) = k
      else
        h%renum(i) = 0
      end if
    end do
    h%nfr = k

    h%rcount(1:h%n) = 0
    do e = 1, hne
      r = h%renum(hrow(e));  c = h%renum(hcol(e))
      if (r == 0 .or. c == 0) cycle
      h%rcount(r) = h%rcount(r) + 1
      if (r /= c) h%rcount(c) = h%rcount(c) + 1
    end do

    h%row_ptr(1) = 1
    do i = 1, h%nfr
      h%row_ptr(i+1) = h%row_ptr(i) + h%rcount(i)
    end do
    h%nnz = h%row_ptr(h%nfr+1) - 1

    h%cur(1:h%nfr) = h%row_ptr(1:h%nfr)
    do e = 1, hne
      r = h%renum(hrow(e));  c = h%renum(hcol(e))
      if (r == 0 .or. c == 0) then
        h%pos1(e) = 0;  h%pos2(e) = 0
        cycle
      end if
      h%pos1(e) = h%cur(r);  h%col(h%cur(r)) = c;  h%cur(r) = h%cur(r) + 1
      if (r /= c) then
        h%pos2(e) = h%cur(c);  h%col(h%cur(c)) = r;  h%cur(c) = h%cur(c) + 1
      else
        h%pos2(e) = 0
      end if
    end do

    h%hne = hne
    h%frvar_prev(1:h%n) = frvar(1:h%n)
    h%has_face = .true.
    h%empty = (h%nfr == 0 .or. h%nnz == 0)
  end subroutine build_pattern

  subroutine setup_spmv(h)
    type(hess_csr_type), intent(inout) :: h
    integer(c_size_t) :: bsz
    integer(c_int) :: st
    integer(8) :: t0, t1, rate

    h%t_prep_ms = 0.0d0
    if (h%empty) return

    call hs_check(hs_create_csr(h%mat, int(h%nfr, c_int64_t), int(h%nfr, c_int64_t), &
                  int(h%nnz, c_int64_t), h%d_row_ptr, h%d_col, h%d_val, &
                  HS_INDEX_32I, HS_INDEX_32I, HS_BASE_ONE, HS_R_64F), "CreateCsr")
    call hs_check(hs_create_dnvec(h%vx, int(h%nfr, c_int64_t), h%d_x, HS_R_64F), "CreateDnVec x")
    call hs_check(hs_create_dnvec(h%vy, int(h%nfr, c_int64_t), h%d_y, HS_R_64F), "CreateDnVec y")

    call hs_check(hs_spmv_buffer_size(h%handle, HS_OP_N, c_loc(one), h%mat, h%vx, &
                  c_loc(zero), h%vy, HS_R_64F, HS_SPMV_ALG_DEFAULT, bsz), "SpMV_bufferSize")
    call hs_alloc_buffer(h%d_buf, bsz)          ! nulo se bsz == 0

    if (h%do_prep) then
      call system_clock(count_rate=rate)
      call system_clock(t0)
      st = hs_spmv_preprocess(h%handle, HS_OP_N, c_loc(one), h%mat, h%vx, &
                  c_loc(zero), h%vy, HS_R_64F, HS_SPMV_ALG_DEFAULT, h%d_buf)
      call hipCheck(hipDeviceSynchronize())
      call system_clock(t1)
      h%t_prep_ms = real(t1 - t0, 8) / real(rate, 8) * 1000.0d0
      if (st /= 0) then
        write(*,'(A,I0,A,I0)') "AVISO: SpMV_preprocess falhou (codigo ", st, "; bufferSize=", bsz
        h%t_prep_ms = -1.0d0
      end if
    end if
    h%spmv_ready = .true.
  end subroutine setup_spmv

  subroutine teardown_spmv(h)
    type(hess_csr_type), intent(inout) :: h
    if (.not. h%spmv_ready) return
    call hs_check(hs_destroy_spmat(h%mat), "DestroySpMat")
    call hs_check(hs_destroy_dnvec(h%vx), "DestroyDnVec x")
    call hs_check(hs_destroy_dnvec(h%vy), "DestroyDnVec y")
    call hs_free_buffer(h%d_buf)
    h%spmv_ready = .false.
  end subroutine teardown_spmv

  ! y = H_ff * x   (vetores de tamanho nfr, em d_x / d_y)
  subroutine hess_csr_spmv(h)
    type(hess_csr_type), intent(in) :: h
    if (h%nfr == 0) return
    if (h%empty) then
      call hipCheck(hipMemset(h%d_y, 0, h%nfr * 8_c_size_t))
    else
      call hs_check(hs_spmv(h%handle, HS_OP_N, c_loc(one), h%mat, h%vx, c_loc(zero), &
                    h%vy, HS_R_64F, HS_SPMV_ALG_DEFAULT, h%d_buf), "SpMV")
    end if
  end subroutine hess_csr_spmv

  subroutine hess_csr_free(h)
    type(hess_csr_type), intent(inout) :: h
    call teardown_spmv(h)
    call hipCheck(hipFree(h%d_row_ptr))
    call hipCheck(hipFree(h%d_col))
    call hipCheck(hipFree(h%d_val))
    call hipCheck(hipFree(h%d_hcoo))
    call hipCheck(hipFree(h%d_pos1))
    call hipCheck(hipFree(h%d_pos2))
    call hipsparseCheck(hipsparseDestroy(h%handle))
  end subroutine hess_csr_free

  subroutine hess_csr_invalidate(h)
    type(hess_csr_type), intent(inout) :: h
    h%has_face = .false.
  end subroutine hess_csr_invalidate

end module hess_csr_mod