module mmio
  implicit none

  ! ==============================================================================
  ! Modulo NIST Matrix Market I/O em Fortran Moderno
  ! Baseado nas especificacoes do NIST (National Institute of Standards)
  ! ==============================================================================

  public :: mm_read_banner, mm_read_mtx_crd_size

contains

  subroutine mm_read_banner(unit, rep, field, symmetry, status)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: rep, field, symmetry
    integer, intent(out) :: status
    character(len=1024) :: line

    status = 0
    read(unit, '(A)', iostat=status) line
    if (status /= 0) return

    call to_lower(line)

    if (index(line, '%%matrixmarket') == 0) then
      status = -1 ! Formato invalido
      return
    end if

    ! Representacao
    if (index(line, 'coordinate') /= 0) then
      rep = 'coordinate'
    else if (index(line, 'array') /= 0) then
      rep = 'array'
    else
      rep = 'unknown'
    end if

    ! Tipo de Dado (field)
    if (index(line, 'pattern') /= 0) then
      field = 'pattern'
    else if (index(line, 'complex') /= 0) then
      field = 'complex'
    else if (index(line, 'integer') /= 0) then
      field = 'integer'
    else
      field = 'real'
    end if

    ! Simetria
    if (index(line, 'skew-symmetric') /= 0) then
      symmetry = 'skew-symmetric'
    else if (index(line, 'symmetric') /= 0) then
      symmetry = 'symmetric'
    else if (index(line, 'hermitian') /= 0) then
      symmetry = 'hermitian'
    else
      symmetry = 'general'
    end if

  end subroutine mm_read_banner

  subroutine mm_read_mtx_crd_size(unit, m, n, nz, status)
    integer, intent(in) :: unit
    integer, intent(out) :: m, n, nz, status
    character(len=1024) :: line

    status = 0
    do
      read(unit, '(A)', iostat=status) line
      if (status /= 0) return
      ! Pula linhas de comentario
      if (line(1:1) /= '%') exit
    end do

    read(line, *, iostat=status) m, n, nz
  end subroutine mm_read_mtx_crd_size

  subroutine to_lower(str)
    character(len=*), intent(inout) :: str
    integer :: idx, code
    do idx = 1, len_trim(str)
      code = ichar(str(idx:idx))
      if (code >= 65 .and. code <= 90) str(idx:idx) = char(code + 32)
    end do
  end subroutine to_lower

end module mmio