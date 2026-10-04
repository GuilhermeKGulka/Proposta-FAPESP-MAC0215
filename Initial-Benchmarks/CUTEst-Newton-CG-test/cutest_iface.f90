module cutest_iface
  ! ====================================================================
  ! Interfaces Fortran para a biblioteca CUTEst (build double precision).
  !
  ! Correspondem as rotinas classicas geradas pelo sifdecoder e usadas
  ! em sources/cutestma.f90 do Algencan 3.1.1. Todas retornam `status`
  ! (integer): 0 = sucesso; != 0 indica erro.
  ! ====================================================================
  implicit none

  interface

    subroutine cutest_udimen( status, input, n )
      implicit none
      integer, intent(out) :: status, n
      integer, intent(in)  :: input
    end subroutine cutest_udimen

    subroutine cutest_usetup( status, input, out, io_buffer, n, X, X_l, X_u )
      implicit none
      integer,   intent(out) :: status, n
      integer,   intent(in)  :: input, out, io_buffer
      real(8),   intent(out) :: X(*), X_l(*), X_u(*)
    end subroutine cutest_usetup

    subroutine cutest_ugr( status, n, X, G )
      implicit none
      integer, intent(out) :: status
      integer, intent(in)  :: n
      real(8), intent(in)  :: X(*)
      real(8), intent(out) :: G(*)
    end subroutine cutest_ugr

    subroutine cutest_ush( status, n, X, hnnz, lim, H_val, H_row, H_col )
      implicit none
      integer, intent(out) :: status, hnnz
      integer, intent(in)  :: n, lim
      real(8), intent(in)  :: X(*)
      real(8), intent(out) :: H_val(*)
      integer, intent(out) :: H_row(*), H_col(*)
    end subroutine cutest_ush

    subroutine cutest_udimsh( status, hnnz )
      implicit none
      integer, intent(out) :: status, hnnz
    end subroutine cutest_udimsh

    subroutine cutest_pname( status, input, pname )
      implicit none
      integer,          intent(out) :: status
      integer,          intent(in)  :: input
      character(len=*), intent(out) :: pname
    end subroutine cutest_pname

  end interface

end module cutest_iface