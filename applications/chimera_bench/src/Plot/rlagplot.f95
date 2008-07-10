SUBROUTINE rlagplot
!-----------------------------------------------------------------------
!
!    File:         rlagplot
!    Module:       rlagplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/24/99
!
!    Purpose:
!      To create files of the radii of selected Lagrangian mass
!       points at selected times for post processing.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!      none
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!  dtimeplot : time (milliseconds) between successive data dumps.
!  dmlag    : enclosed rest mass between selected Lagrangian mass points
!  irlagplt : time when the last dump was printed
!  irlagplt : 0, omit creating rmlag files
!             1, create rmlag files
!  nrlagplt : unit number for the rmlag files
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, mdl_cnfg_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero
USE physcnst_module, ONLY : msolar

USE edit_module, ONLY : nprint, nrlagplt, irlagplt, dtimeplot, head, &
& nlog, dmlag
USE mdl_cnfg_module, ONLY : jr_max, r, rstmss
USE t_cntrl_module, ONLY : time, dtnph

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=28)                :: rlagfile

LOGICAL                           :: lprint
LOGICAL                           :: first = .true.

INTEGER                           :: i             ! Lagrangian mass index
INTEGER                           :: j             ! radial zone index
INTEGER                           :: jmlag         ! rstmss(jmlag) > masslag(i) > rstmss(jd-1)
INTEGER                           :: istat         ! open and close file flag
INTEGER, PARAMETER                :: nmd = 1000    ! upper bound of Lagrangian mass number
INTEGER                           :: nm            ! number of Lagrangian masses
INTEGER                           :: l             ! setof 15 Lagrangian mass points index
INTEGER                           :: lcnt          ! number of sets of 15 Lagranian mass sehlls to edit
INTEGER                           :: lmin          ! lower index of set of 15 Lagrangian mass points
INTEGER                           :: lmax          ! upper index of set of 15 Lagrangian mass points


INTEGER                           :: itime         ! parameter for criterion for writing to nuradplot files
INTEGER                           :: itimeprev     ! parameter for criterion for writing to nuradplot files

REAL(KIND=double)                 :: tmult         ! multiplier for criterion for writing to nuradplot files
REAL(KIND=double), PARAMETER      :: km_cm = 1.0d-05

REAL(KIND=double), DIMENSION(nmd) :: masslag       ! Lagrangian mass shells
REAL(KIND=double), DIMENSION(nmd) :: rmlag         ! radii of the Lagrangian mass shells
!     *        j,jmlag,l,lcnt,lmax,lmin,nm
!      double precision km_cm,masslag,rmlag
!                                                                      c
!                                                                      c

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    3 format (1x,a128)
    5 format (1x,128('*')/)
  101 format (' jmlag cannot be found in subroutine rlagplot, i=',i3, &
& ' masslag(i)=',es10.3,' rstmss(jr_max)=',es10.3)
  221 format ('     time  ',15f7.3)
  223 format ('     sec'/)
  225 format (1x,es11.3,15f7.3)
  501 format (' Unable to cover rstmss in subroutine rlagplot, masslag(1000)=', &
& es10.3,' dmlag=',es10.3,' rstmss(jr_max)/msolar=',es10.3)
 2001 format (' nm too large, nm=',i5,' lmin=',i5,' lmax=',i5)
 8001 format (' File ',a18,' cannot be opened in subroutime rlagplot')
 9001 format (' File ',a18,' cannot be closed in subroutime rlagplot')

!-----------------------------------------------------------------------
!  Return IF irlagplt = 0.                                       c
!-----------------------------------------------------------------------

iF ( irlagplt == 0 ) RETURN

!-----------------------------------------------------------------------
!  Initialize.                                                   c
!                                                                      c
!  Compute Lagrangian masses.                                c
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  masslag(1)       = zero
  DO i = 2,nmd
    masslag(i)     = masslag(i-1) + dmlag
    IF ( masslag(i) > rstmss(jr_max)/msolar ) THEN
      nm           = i - 1
      masslag(1)   = 0.005d0
      GO TO 600
    END IF ! masslag(i) > rstmss(jr_max)/msolar
  END DO
  WRITE (nprint,501) masslag(1000), dmlag, rstmss(jr_max)/msolar
  WRITE (nlog,501) masslag(1000), dmlag, rstmss(jr_max)/msolar
  STOP
END IF ! first = true

!-----------------------------------------------------------------------
!  Time criterion.
!-----------------------------------------------------------------------

  600 CONTINUE
lprint             = .false.
tmult              = 1.d+3/dtimeplot
itime              = INT( time * tmult )
itimeprev          = INT( ( time - dtnph ) * tmult )
IF ( itime > itimeprev ) lprint = .true.

IF ( .not. lprint ) RETURN

!-----------------------------------------------------------------------
!  Find the radii, rmlag(i), corresponding to the lagrangian masses,
!   masslag(i).
!-----------------------------------------------------------------------

DO i = 1,nm

!-----------------------------------------------------------------------
!  (1) Find jmlag such that rstmss(jmlag) > masslag(i) > rstmss(jd-1).
!-----------------------------------------------------------------------

  DO j = 2,jr_max
    IF ( rstmss(j)/msolar >= masslag(i) ) THEN
      jmlag        = j
      go to 110
    END IF ! rstmss(j) > masslag(i)
  END DO ! j = 2,jr_max
  WRITE (nprint,101) i, masslag(i), rstmss(jr_max)
  WRITE (nlog,101) i, masslag(i), rstmss(jr_max)
  RETURN

!-----------------------------------------------------------------------
!  (2) Interpolate rmlag(i) in mass to masslag(i).
!-----------------------------------------------------------------------

  110  CONTINUE
  rmlag(i)         = rinterp( r(jmlag), r(jmlag-1),                     &
&                    rstmss(jmlag  )/msolar, masslag(i),                &
&                    rstmss(jmlag-1)/msolar ) * km_cm

END DO ! i = 1,nm

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Divide output into sets of 15
!-----------------------------------------------------------------------

DO lcnt = 1,100
  lmax             = MIN(   lcnt       * 15     , nm )
  lmin             = MAX( ( lcnt - 1 ) * 15 + 1 ,  1 )

!-----------------------------------------------------------------------
!  Write to rlagfile ("variable lagplot files").
!-----------------------------------------------------------------------

  WRITE (rlagfile,'(a24,i2.2,a2)') 'Data3/Plot_Files/rlagplot',lcnt,'.d'

!-----------------------------------------------------------------------
!  Open lagplot1.d file
!-----------------------------------------------------------------------

  OPEN (UNIT=nrlagplt,FILE=rlagfile,STATUS='new',IOSTAT=istat)
  IF ( istat == 0 ) THEN
    WRITE (nrlagplt,3) head
    WRITE (nrlagplt,5)
    WRITE (nrlagplt,221) (masslag(l),l=lmin,lmax)
    WRITE (nrlagplt,223)
  END IF ! istat == 0
       
  IF ( istat /= 0 ) OPEN (UNIT=nrlagplt, FILE=rlagfile, STATUS='old', &
&  POSITION='append', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8001) rlagfile
    WRITE (nlog,8001) rlagfile
    STOP
  END IF ! istat /= 0

  WRITE (nrlagplt,225) time, (dlog10(rmlag(l)),l=lmin,lmax)

!-----------------------------------------------------------------------
!  Close lagplot1.d file
!-----------------------------------------------------------------------

  CLOSE (UNIT=nrlagplt, STATUS='keep', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,9001) rlagfile
    WRITE (nlog,9001) rlagfile
    STOP ! istat /= 0
  END IF

  IF ( lmax >= nm ) RETURN

END DO ! lcnt = 1,100

WRITE (nprint,2001) nm, lmin, lmax
WRITE (nlog,2001) nm, lmin, lmax
RETURN

CONTAINS
  REAL (KIND=double) FUNCTION rinterp(a,b,x,y,z)

  REAL (KIND=double) :: a
  REAL (KIND=double) :: b
  REAL (KIND=double) :: x
  REAL (KIND=double) :: y
  REAL (KIND=double) :: z

  rinterp        = b + ( a - b ) * ( y - z )/( x - z )
  END FUNCTION rinterp

END
