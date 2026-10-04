! Varredura: y = A*x (A 3x4, exemplo da doc do hipSPARSE), double, API generica.
! Para cada base (0/1) e algoritmo (0..4) imprime o status de bufferSize,
! preprocess e SpMV (sem abortar). Esperado em y: 6 9 21.
program test_generic
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipsparse
  use hipsparse_generic
  implicit none

  write(*,'(A,5I6)') "valores usados (op, idx32, base0, base1, R64F): ", &
        HS_OP_N, HS_INDEX_32I, HS_BASE_ZERO, HS_BASE_ONE, HS_R_64F
  write(*,'(A)') "status: 0 = ok, 3 = INVALID_VALUE, -1 = nao executado"

  call scan(HS_BASE_ZERO)
  call scan(HS_BASE_ONE)

contains

  subroutine scan(base)
    integer(c_int), intent(in) :: base
    integer(c_int), target :: rp(4), ci(8)
    real(c_double), target :: va(8), x(4), y(3), alpha, beta
    type(c_ptr) :: handle, mat, vx, vy, d_rp, d_ci, d_va, d_x, d_y, buf
    integer(c_size_t) :: bsz
    integer(c_int) :: alg, st_buf, st_pre, st_mv
    integer :: i

    rp = [0, 3, 5, 8] + base
    ci = [0, 1, 3, 1, 2, 0, 2, 3] + base
    va = [(real(i, c_double), i = 1, 8)]
    x = 1.0d0
    alpha = 1.0d0
    beta = 0.0d0

    call hipcheck(hipMalloc(d_rp, 16_c_size_t))
    call hipcheck(hipMalloc(d_ci, 32_c_size_t))
    call hipcheck(hipMalloc(d_va, 64_c_size_t))
    call hipcheck(hipMalloc(d_x, 32_c_size_t))
    call hipcheck(hipMalloc(d_y, 24_c_size_t))
    call hipcheck(hipMemcpy(d_rp, c_loc(rp), 16_c_size_t, hipMemcpyHostToDevice))
    call hipcheck(hipMemcpy(d_ci, c_loc(ci), 32_c_size_t, hipMemcpyHostToDevice))
    call hipcheck(hipMemcpy(d_va, c_loc(va), 64_c_size_t, hipMemcpyHostToDevice))
    call hipcheck(hipMemcpy(d_x, c_loc(x), 32_c_size_t, hipMemcpyHostToDevice))

    call hipsparsecheck(hipsparseCreate(handle))
    call hs_check(hs_create_csr(mat, 3_c_int64_t, 4_c_int64_t, 8_c_int64_t, d_rp, d_ci, d_va, &
                                HS_INDEX_32I, HS_INDEX_32I, base, HS_R_64F), "CreateCsr")
    call hs_check(hs_create_dnvec(vx, 4_c_int64_t, d_x, HS_R_64F), "CreateDnVec x")
    call hs_check(hs_create_dnvec(vy, 3_c_int64_t, d_y, HS_R_64F), "CreateDnVec y")

    do alg = 0, 4
      y = 0.0d0
      call hipcheck(hipMemcpy(d_y, c_loc(y), 24_c_size_t, hipMemcpyHostToDevice))
      st_pre = -1
      st_mv = -1
      bsz = 0
      st_buf = hs_spmv_buffer_size(handle, HS_OP_N, c_loc(alpha), mat, vx, c_loc(beta), vy, &
                                   HS_R_64F, alg, bsz)
      if (st_buf == 0) then
        ! buffer folgado (>= 1 MiB) para descartar problema de tamanho
        call hipcheck(hipMalloc(buf, max(bsz, 1048576_c_size_t)))
        st_pre = hs_spmv_preprocess(handle, HS_OP_N, c_loc(alpha), mat, vx, c_loc(beta), vy, &
                                    HS_R_64F, alg, buf)
        st_mv = hs_spmv(handle, HS_OP_N, c_loc(alpha), mat, vx, c_loc(beta), vy, &
                        HS_R_64F, alg, buf)
        call hipcheck(hipMemcpy(c_loc(y), d_y, 24_c_size_t, hipMemcpyDeviceToHost))
        call hipcheck(hipFree(buf))
      end if
      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,3F7.1)') "base ", base, " alg ", alg, &
            ": bufsize=", bsz, " | bufferSize=", st_buf, " preprocess=", st_pre, &
            " spmv=", st_mv, " | y=", y
    end do

    call hs_check(hs_destroy_spmat(mat), "DestroySpMat")
    call hs_check(hs_destroy_dnvec(vx), "DestroyDnVec x")
    call hs_check(hs_destroy_dnvec(vy), "DestroyDnVec y")
    call hipsparsecheck(hipsparseDestroy(handle))
    call hipcheck(hipFree(d_rp)); call hipcheck(hipFree(d_ci)); call hipcheck(hipFree(d_va))
    call hipcheck(hipFree(d_x));  call hipcheck(hipFree(d_y))
  end subroutine scan

end program test_generic