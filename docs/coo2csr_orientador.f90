! Rotina coo2csr fornecida pelo orientador (21/11/2026), copiada sem alteracoes.
! E a base pedida para a conversao COO -> CSR na fase 6 (ver docs/PLANO.md, secao 7.5).
! Arquivo de referencia: NAO faz parte do build.
!
! Observacoes para adaptar:
! - a conversao e feita in-place (sobrescreve arow/acol/aval), com custo O(nnz)
!   e sem memoria extra alem de alen/asta;
! - no fim, asta(i) e o inicio da linha i (base 1) e alen(i) e o numero de
!   elementos dela; a ordem das colunas dentro de cada linha nao e garantida;
! - nao filtra pela face (frvar), nao renumera e nao espelha o triangulo:
!   isso precisa ser feito antes da chamada.

! ******************************************************************
! ******************************************************************

subroutine coo2csr(m,nnz,arow,acol,aval,alen,asta)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: m,nnz

  ! ARRAY ARGUMENTS
  integer,      intent(inout) :: acol(nnz),alen(m),arow(nnz),asta(m)
  real(kind=8), intent(inout) :: aval(nnz)

  ! This subroutines converts a matrix from coordinate format to
  ! compressed sparse row format.

  ! LOCAL SCALARS
  integer      :: i,j,col,coltmp,row,rowtmp
  real(kind=8) :: val,valtmp

  ! The if below is wrong, since arrays alen and asta remain un-initialized
  ! if ( nnz .eq. 0 ) return

  alen(1:m) = 0

  do i = 1,nnz
     row = arow(i)
     alen(row) = alen(row) + 1
  end do

  asta(1) = 1
  do i = 2,m
     asta(i) = asta(i-1) + alen(i-1)
  end do

  do i = 1,nnz

     val = aval(i)
     col = acol(i)
     row = arow(i)

     arow(i) = - 1

     do while ( row .ge. 0 )
        j = asta(row)
        asta(row) = j + 1

        valtmp = aval(j)
        coltmp = acol(j)
        rowtmp = arow(j)

        aval(j) = val
        acol(j) = col
        arow(j) = - 1

        val = valtmp
        col = coltmp
        row = rowtmp
     end do

  end do

  asta(1:m) = asta(1:m) - alen(1:m)

end subroutine coo2csr
