SUBROUTINE nuclear_time_step( jr_min, jr_max, t, idim, dtnph )
!-----------------------------------------------------------------------
!
!    File:         nuclear_time_step
!    Module:       nuclear_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/04
!
!    Purpose:
!      To select the time step restricted by nuclear burn.

!      Temperature change due to nuclear burning is stored in dT_burn(j).
!       t_cntl_burn(1) is the maximum permitted abs( dT_burn(j)/t(j) ).
!       Temperature change time step control bypassed if t_cntl_burn(1) <= 0.
!
!      Maximum relative composition change due to nuclear burning is given by dynmax.
!       t_cntl_burn(2) is the maximum permitted abs( dyn/( yn + ynmin ) )
!       Composition change time step control is bypassed if t_cntl_burn(2) <= 0.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min   : minimum radial zone index
!  jr_max   : maximum radial zone index
!  t        : temperature (K)
!  idim     : array dimension
!  dtnph    : current time step
!
!    Output arguments:
!        none
!
!    Include files:
!      kind_module, numerical_module
!      nucbrn_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero

USE nucbrn_module, ONLY: jdynmax, dynmax, jdTmax, dTmax, dtime_burn, dT_burn, t_cntl_burn

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min            ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max            ! maximum radial zone index
INTEGER, INTENT(in)              :: idim              ! array dimension

REAL(KIND=double), INTENT(in)    :: dtnph             ! current time step
REAL(KIND=double), INTENT(in), DIMENSION(idim+1) :: t ! temperature (K)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j                 ! radial zone index

REAL(KIND=double)                :: d                 ! working variable

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ TEMPERATURE CHANGE TIME STEP DUE /////
!            \\\\\ TO NUCLEAR BURNING DTIME_BURN(1) /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Temperature change time step restriction due to nuclear
!         burning is computed if t_cntl_burn(1) > 0
!-----------------------------------------------------------------------

IF ( t_cntl_burn(1) > zero ) THEN

  dTmax            = zero

!........Get maximum relative change of temperature

  DO j = jr_min,jr_max
    d              = DABS( dT_burn(j)/t(j) )
    IF ( d > dTmax ) THEN
      dTmax        = d
      jdTmax       = j
    END IF
  END DO

!-----------------------------------------------------------------------
!        Temperature change time step control due to nuclear burning
!         only if
!
!             dmax = max( abs( dtmpmn(j,4,i_ray) )/t(j) ) > 0.
!-----------------------------------------------------------------------

  IF ( dTmax > zero ) THEN
    dtime_burn(1)  = t_cntl_burn(1) * dtnph/dTmax
  END IF
  
END IF ! t_cntl_burn(1) > zero

!-----------------------------------------------------------------------
!
!            \\\\\ COMPOSITION CHANGE TIME STEP DUE /////
!            \\\\\ TO NUCLEAR BURNING DTIME_BURN(2) /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Composition change time step control is used if
!         t_cntl_burn(2) > 0 and dynmax > zero
!-----------------------------------------------------------------------

IF ( t_cntl_burn(2) > zero  .and.  dynmax > zero ) THEN
  dtime_burn(2)    = t_cntl_burn(2) * dtnph/dynmax
END IF

RETURN
END SUBROUTINE nuclear_time_step
