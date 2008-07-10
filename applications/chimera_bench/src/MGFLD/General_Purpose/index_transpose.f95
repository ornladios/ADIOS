SUBROUTINE index_transpose( id, it, iy, ids, its, iys, id_trans, it_trans, &
& iy_trans )
!-----------------------------------------------------------------------
!
!    File:         index_transpose
!    Module:       index_transpose
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/05/08
!
!    Purpose:
!        To compute the independent variable grid indices and
!         interpolation coefficients for subroutine eqstt_x.
!
!    Input arguments:
!
!  id             : densiity grid index
!  it             : temperature grid index
!  iy             : electron fraction grid index
!  ids            : prior densiity grid index
!  its            : prior temperature grid index
!  iys            : prior electron fraction grid index
!
!    Output arguments:
!
!  id_trans       : old density index to be transposes
!  it_trans       : old temperature index to be transposes
!  iy_trans       : old electron fraction index to be transposes
!
!    Subprograms called:
!        none
!
!    Include files:
!  edit_module, parallel_module  
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)    :: id           ! density cube index
INTEGER, INTENT(in)    :: it           ! temperature cube index
INTEGER, INTENT(in)    :: iy           ! electron fraction cube index

INTEGER, INTENT(in)    :: ids          ! prior density cube index
INTEGER, INTENT(in)    :: its          ! prior temperature cube index
INTEGER, INTENT(in)    :: iys          ! prior electron fraction cube index

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(2,2,2)  :: id_trans  ! old density index to be transposes
INTEGER, INTENT(out), DIMENSION(2,2,2)  :: it_trans  ! old temperature index to be transposes
INTEGER, INTENT(out), DIMENSION(2,2,2)  :: iy_trans  ! old electron fraction index to be transposes

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                :: ida           ! density do index
INTEGER                :: ita           ! temperature do index
INTEGER                :: iya           ! electron fraction do index

!-----------------------------------------------------------------------
!        Formate
!-----------------------------------------------------------------------

 1001 FORMAT (' indexing problem in subroutine index_transpose, id =',i6, &
& ' it=',i6, ' iy=',i6, ' ids=',i6, ' its=',i6,' iys=',i6)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

id_trans           = 0
it_trans           = 0
iy_trans           = 0

!-----------------------------------------------------------------------
!
!             \\\\\ COMPUTE THE TRANSPOSE ALGORITHMS /////!             
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Density increases above table
!-----------------------------------------------------------------------

IF ( id - ids == 1 ) THEN
  DO ita = 1, 2
    DO iya = 1, 2
      id_trans(1,ita,iya) = 2
      it_trans(1,ita,iya) = ita
      iy_trans(1,ita,iya) = iya
    END DO ! iya = 1, 2
  END DO ! ita = 1, 2
  RETURN
END IF ! id - ids == -1

!-----------------------------------------------------------------------
!  Density decreases below table
!-----------------------------------------------------------------------

IF ( id - ids == -1 ) THEN
  DO ita = 1, 2
    DO iya = 1, 2
      id_trans(2,ita,iya) = 1
      it_trans(2,ita,iya) = ita
      iy_trans(2,ita,iya) = iya
    END DO ! iya = 1, 2
  END DO ! ita = 1, 2
  RETURN
END IF ! id - ids == 1

!-----------------------------------------------------------------------
!  Temperature increases above table
!-----------------------------------------------------------------------

IF ( it - its == 1 ) THEN
  DO ida = 1, 2
    DO iya = 1, 2
      id_trans(ida,1,iya) = ida
      it_trans(ida,1,iya) = 2
      iy_trans(ida,1,iya) = iya
    END DO ! ita = 1, 2
  END DO ! ida = 1, 2
  RETURN
END IF ! it - its == -1

!-----------------------------------------------------------------------
!  Temperature decrease below table
!-----------------------------------------------------------------------

IF ( it - its == -1 ) THEN
  DO ida = 1, 2
    DO iya = 1, 2
      id_trans(ida,2,iya) = ida
      it_trans(ida,2,iya) = 1
      iy_trans(ida,2,iya) = iya
    END DO ! iya = 1, 2
  END DO ! ida = 1, 2
  RETURN
END IF ! it - its == 1

!-----------------------------------------------------------------------
!  Electron fraction increases above table
!-----------------------------------------------------------------------

IF ( iy - iys == 1 ) THEN
  DO ida = 1, 2
    DO ita = 1, 2
      id_trans(ida,ita,1) = ida
      it_trans(ida,ita,1) = ita
      iy_trans(ida,ita,1) = 2
    END DO ! ita = 1, 2
  END DO ! ida = 1, 2
  RETURN
END IF ! iy - iys == -1

!-----------------------------------------------------------------------
!  Electron fraction decreases below table
!-----------------------------------------------------------------------

IF ( iy - iys == -1 ) THEN
  DO ida = 1, 2
    DO ita = 1, 2
      id_trans(ida,ita,2) = ida
      it_trans(ida,ita,2) = ita
      iy_trans(ida,ita,2) = 1
    END DO ! ita = 1, 2
  END DO ! ida = 1, 2
  RETURN
END IF ! iy - iys == 1

!-----------------------------------------------------------------------
!  Indexing problem if none of the above choices selected
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nprint,1001) id, it, iy, ids, its, iys
WRITE (nlog,1001) id, it, iy, ids, its, iys
STOP

END SUBROUTINE index_transpose
