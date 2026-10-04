! Bindings diretos (bind(c)) para a API generica do hipSPARSE (SpMV), pois o
! hipfort instalado nao tem interfaces utilizaveis para elas. Os valores dos
! enums vem do proprio hipfort (os mesmos que ja funcionam nas chamadas legadas),
! ja que a numeracao varia entre versoes do hipSPARSE (ex.: NON_TRANSPOSE e 0
! aqui, nao 111).
module hipsparse_generic
  use iso_c_binding
  use hipfort
  use hipfort_hipsparse
  implicit none
  private
  public :: HS_OP_N, HS_INDEX_32I, HS_BASE_ZERO, HS_BASE_ONE, HS_R_32F, HS_R_64F
  public :: HS_SPMV_ALG_DEFAULT, hs_check
  public :: hs_create_csr, hs_create_dnvec, hs_destroy_spmat, hs_destroy_dnvec
  public :: hs_spmv_buffer_size, hs_spmv_preprocess, hs_spmv

  integer(c_int), parameter :: HS_OP_N = int(HIPSPARSE_OPERATION_NON_TRANSPOSE, c_int)
  integer(c_int), parameter :: HS_INDEX_32I = int(HIPSPARSE_INDEX_32I, c_int)
  integer(c_int), parameter :: HS_BASE_ZERO = int(HIPSPARSE_INDEX_BASE_ZERO, c_int)
  integer(c_int), parameter :: HS_BASE_ONE = int(HIPSPARSE_INDEX_BASE_ONE, c_int)
  integer(c_int), parameter :: HS_R_32F = int(HIP_R_32F, c_int)
  integer(c_int), parameter :: HS_R_64F = int(HIP_R_64F, c_int)
  integer(c_int), parameter :: HS_SPMV_ALG_DEFAULT = int(HIPSPARSE_MV_ALG_DEFAULT, c_int)

  interface
    function hs_create_csr(spmat, rows, cols, nnz, row_ptr, col_ind, val, &
                           row_type, col_type, base, val_type) &
             bind(c, name="hipsparseCreateCsr") result(stat)
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr) :: spmat
      integer(c_int64_t), value :: rows, cols, nnz
      type(c_ptr), value :: row_ptr, col_ind, val
      integer(c_int), value :: row_type, col_type, base, val_type
      integer(c_int) :: stat
    end function

    function hs_create_dnvec(vec, n, values, val_type) &
             bind(c, name="hipsparseCreateDnVec") result(stat)
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr) :: vec
      integer(c_int64_t), value :: n
      type(c_ptr), value :: values
      integer(c_int), value :: val_type
      integer(c_int) :: stat
    end function

    function hs_destroy_spmat(mat) bind(c, name="hipsparseDestroySpMat") result(stat)
      import :: c_ptr, c_int
      type(c_ptr), value :: mat
      integer(c_int) :: stat
    end function

    function hs_destroy_dnvec(vec) bind(c, name="hipsparseDestroyDnVec") result(stat)
      import :: c_ptr, c_int
      type(c_ptr), value :: vec
      integer(c_int) :: stat
    end function

    function hs_spmv_buffer_size(handle, op, alpha, mat, x, beta, y, &
                                 compute_type, alg, buffer_size) &
             bind(c, name="hipsparseSpMV_bufferSize") result(stat)
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: handle, alpha, mat, x, beta, y
      integer(c_int), value :: op, compute_type, alg
      integer(c_size_t) :: buffer_size
      integer(c_int) :: stat
    end function

    function hs_spmv_preprocess(handle, op, alpha, mat, x, beta, y, &
                                compute_type, alg, buffer) &
             bind(c, name="hipsparseSpMV_preprocess") result(stat)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle, alpha, mat, x, beta, y, buffer
      integer(c_int), value :: op, compute_type, alg
      integer(c_int) :: stat
    end function

    function hs_spmv(handle, op, alpha, mat, x, beta, y, compute_type, alg, buffer) &
             bind(c, name="hipsparseSpMV") result(stat)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle, alpha, mat, x, beta, y, buffer
      integer(c_int), value :: op, compute_type, alg
      integer(c_int) :: stat
    end function
  end interface

contains

  ! Aborta dizendo QUAL chamada falhou
  subroutine hs_check(stat, label)
    integer(c_int), intent(in) :: stat
    character(len=*), intent(in) :: label
    if (stat /= 0) then
      write(*,'(A,A,A,I0)') "HIPSPARSE ERROR em ", label, ": codigo ", stat
      stop 1
    end if
  end subroutine hs_check

end module hipsparse_generic