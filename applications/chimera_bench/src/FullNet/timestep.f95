SUBROUTINE timestep( kstep, dtime )
!===============================================================================
!  This routine calculates the trial timestep, based on the prior timestep.  
!  For tdel >0, this calculation is based on the relative changes of the 
!               abundances during the previous step in the evolution.  
!           =0, the timestep is calculated using the time derivatives of 
!               the abundances based on reaction rates.  
!           <0, the timestep is held constant, tdel=-tdelo.
!  The timestep is further constrained to be no more than a factor tdelmm 
!  larger than the timestep in the previous step.  There is also the 
!  provision for limiting the timestep in the event that the thermodynamic 
!  conditions are changing too rapidly.
!===============================================================================

USE kind_module, ONLY : double

USE abundances, ONLY : y, yo, ydot
USE conditions, ONLY : t, tt, t9t, tdel, rhot
USE controls, ONLY : tdelmm, ytime, changemx, idiag
USE edit_module, ONLY : nlog
USE nuclear_data, ONLY : nname
USE nuc_number, ONLY : ny
USE thermo_data, ONLY : tstop

IMPLICIT none
SAVE

INTEGER           :: ints(1),intso(1)   ! nucleus governing timestep
INTEGER           :: i,kstep

REAL(KIND=double) :: ydotoy(ny)
REAL(KIND=double) :: changeth, changest, tdelo, t9old, rhold, dt, dtherm
REAL(KIND=double) :: tdels, tdelm, tdeln, dtime

!-----------------------------------------------------------------------
!  Retain old values of timestep and thermo and calculate remaining time  
!-----------------------------------------------------------------------

changeth             = .1d0
changest             = .1d0
tdelo                = tdel
t9old                = t9t
rhold                = rhot
tdels                = tstop - t
IF ( tdels <= 0.0d0 ) tdels = dtime
dt                   = 1.0d20

!-----------------------------------------------------------------------
!  If this is not the initial timestep, calculate timestep from changes 
!  in last timestep.
!-----------------------------------------------------------------------

IF ( tdelo > 0.0d0 ) THEN
  tdelm              = tdelmm * tdelo
  WHERE ( y > ytime )
    ydotoy           = DABS((y-yo)/y)
  ELSEWHERE
    ydotoy           = 0.0d0
  END WHERE
  ints               = maxloc(ydotoy)
  tdeln              = changemx * tdelo/ydotoy(ints(1))
  tdel               = MIN( tdeln, tdels, tdelm )

!-----------------------------------------------------------------------
!  If this is an initial timestep, yo does not exist, so calculate 
!  timestep from derivatives.
!-----------------------------------------------------------------------

ELSE IF ( tdelo == 0.0d0 ) THEN
  intso(1)           = 0
!       If(t9t.le.1.0) THEN
!         changest=.0001*changemx
!       ELSE IF(t9t.le.3.0) THEN
!         changest=.001*changemx
!       Else
!         changest=.01*changemx
!       Endif
  CALL cross_sect
  CALL yderiv

  WHERE ( y > ytime )
    ydotoy           = DABS(ydot/y)
  ELSEWHERE
    ydotoy           = 0.0d0
  END WHERE
  ints               = maxloc(ydotoy)
  tdeln              = changest * changemx/ydotoy(ints(1))
  tdel               = MIN( tdels, tdeln )

!-----------------------------------------------------------------------
!  Keep timestep constant
!-----------------------------------------------------------------------

ELSE
  tdel               = -tdelo
END IF

!-----------------------------------------------------------------------
!  Diagnostic Output
!-----------------------------------------------------------------------

IF ( idiag >= 1 ) WRITE (50,"(a4,i5,4es12.4)") &
&      'tdel',kstep,tdel,tdeln,tdelo,tdels
!     If(idiag>=2) THEN
!       Write(50,"(a5,i4,2es12.4)") 
!    &    (nname(k),k,y(k),ydotoy(k),k=1,ny)
!     END IF

!-----------------------------------------------------------------------
!  Retain the index of the species setting the timestep 
!-----------------------------------------------------------------------

IF ( ints(1) /= intso(1) ) THEN
  IF ( idiag >= 1 ) WRITE (50,*) 'ITC ',nname(ints(1)),t,tdel
  intso              = ints
END IF

!-----------------------------------------------------------------------
!  Limit timestep if Thermodynamic variation is too large
!-----------------------------------------------------------------------

DO i=1,10
  tt                 = t + tdel
  CALL t9rhofind(kstep)
  IF ( t9old > 0 ) THEN
    dtherm           = DABS( t9t - t9old )/t9old + DABS( rhot - rhold )/rhold
  END IF
  IF ( dtherm < changeth ) EXIT
  tdel           = .5d0 * tdel
  IF (i==10) WRITE (6,*) 'Error in Thermo variations after ',i, &
&        'reductions',tdel,t9t,rhot
END DO

!     If(idiag>=1) Write(50,"(a5,i5,2es12.4)") 
!    &  'T9del',kstep,tdel,dtherm

RETURN
END SUBROUTINE timestep                                                                 
