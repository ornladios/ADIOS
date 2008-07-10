SUBROUTINE rnuplot
!-----------------------------------------------------------------------
!
!    File:         rnuplot
!    Module:       rnuplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To creates files for radii at selected densities vs time
!       plots.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, edit_module, mdl_cnfg_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nrnuplt, nrnu, irnuplt, ncyrnu, intrnu, &
& dtimeplot, rhobar, data_path
USE mdl_cnfg_module, ONLY : jr_max, r
USE t_cntrl_module, ONLY : dtnph, time

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

LOGICAL                           :: first = .true.
LOGICAL                           :: lprint

INTEGER                           :: i            ! fixed density index
INTEGER                           :: istat        ! open-close file error flag
INTEGER                           :: itime        ! used to evaluate time criterion for editing
INTEGER                           :: itimeprev    ! used to evaluate time criterion for editing
INTEGER                           :: j            ! radial zone index
INTEGER                           :: jd           ! radial zone index just outward of rhornu(i)
INTEGER                           :: j_s          ! saved radial zone index
INTEGER                           :: nrnutest     ! used for cycle criterion for editing
INTEGER, PARAMETER                :: nzero = 0    ! 0

REAL(KIND=double)                 :: tmult        ! used to evaluate time criterion for editing

REAL(KIND=double), DIMENSION(6)   :: rhornu       ! selected densities
REAL(KIND=double), DIMENSION(6)   :: rnu          ! radius interpolated to rhornu

!-----------------------------------------------------------------------
!        Formats.
!-----------------------------------------------------------------------

  201 FORMAT (es15.8,6(es11.3))
 1001 FORMAT (' jd cannot be found in subroutine rnuplot for rhornu(', &
& i1,')=',es10.3)
 8001 FORMAT (' File nrnuplt cannot be opened in subroutime rnuplot')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if irnuplt = 0.
!-----------------------------------------------------------------------

IF ( irnuplt == 0 ) RETURN

!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------

IF ( first ) THEN
 
  rhornu(1)        = 1.d+13
  rhornu(2)        = 1.d+12
  rhornu(3)        = 1.d+11
  rhornu(4)        = 1.d+10
  rhornu(5)        = 1.d+9
  rhornu(6)        = 1.d+8
  first            = .false.

END IF

!-----------------------------------------------------------------------
!
!           \\\\\ CRITERIA FOR WRITING TO RNUPLOT FILE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cycle criterion.
!-----------------------------------------------------------------------

lprint             = .false.

IF ( ncycle > ncyrnu(2) ) THEN
  nrnutest         = intrnu(3)
ELSE IF ( ncycle > ncyrnu(1) ) THEN
  nrnutest         = intrnu(2)
ELSE
  nrnutest         = intrnu(1)
END IF ! ncycle

nrnu               = nrnu + 1

IF ( nrnu >= nrnutest ) THEN
  lprint           = .true.
  nrnu             = nzero
END IF ! nrnu > nrnutest
!-----------------------------------------------------------------------
!  Time criterion.
!-----------------------------------------------------------------------

tmult              = 1.d+3/dtimeplot
itime              = INT( time * tmult )
itimeprev          = INT( ( time - dtnph ) * tmult )
IF ( itime > itimeprev ) lprint = .true.

!-----------------------------------------------------------------------
!  Return if lprint = false.
!-----------------------------------------------------------------------

IF ( .not. lprint ) RETURN

!-----------------------------------------------------------------------
!  Start cycle here
!-----------------------------------------------------------------------

DO i = 1,6

!-----------------------------------------------------------------------
!  Find jd such that
!     rhobar(jd) < rhornu(i)
!  and
!     rhobar(jd-1) > rhornu(i)
!-----------------------------------------------------------------------

  IF ( rhobar(2) < rhornu(i) ) THEN
    rnu(i)       = zero
    CYCLE ! i = 1,6 loop
  END IF ! rhobar(2) < rhornu(i)

  j_s            = 0
  DO j = 2,jr_max
    IF ( rhobar(j) <= rhornu(i) ) THEN
      jd        = j
      EXIT
    END IF ! rhobar(j) < rhornu(i)
    j_s         = j
  END DO ! j = 2,jr_max
  
  IF ( j_s == jr_max ) THEN
    WRITE (nprint,1001) i,rhornu(i)
    CYCLE ! i = 1,6 loop
  END IF

!-----------------------------------------------------------------------
!  Compute the radius rnu(i)
!-----------------------------------------------------------------------

   rnu(i)        = rinterp( r(jd), r(jd-1), rhobar(jd), rhornu(i), rhobar(jd-1) )

END DO ! i = 1,6

!-----------------------------------------------------------------------
!  Open rnu.d file.
!-----------------------------------------------------------------------

OPEN (UNIT=nrnuplt,FILE=TRIM(data_path)//'/Plot_Files/rnu.d', &
& STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nrnuplt,FILE=TRIM(data_path)//'/Plot_Files/rnu.d', &
& STATUS='old', POSITION='append',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,8001)
  STOP
END IF

!-----------------------------------------------------------------------
!  Write to rnu.d file.
!-----------------------------------------------------------------------

WRITE (nrnuplt,201) time,(rnu(i),i=1,6)

!-----------------------------------------------------------------------
!  Close rnu.d file.
!-----------------------------------------------------------------------

CLOSE (UNIT=nrnuplt, STATUS='keep')

!-----------------------------------------------------------------------
!  Return.
!-----------------------------------------------------------------------

RETURN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS
REAL (KIND=double) FUNCTION rinterp(a,b,x,y,z)

REAL (KIND=double) :: a
REAL (KIND=double) :: b
REAL (KIND=double) :: x
REAL (KIND=double) :: y
REAL (KIND=double) :: z

rinterp      = b + ( a - b ) * ( y - z )/( x - z )

END FUNCTION rinterp

END
