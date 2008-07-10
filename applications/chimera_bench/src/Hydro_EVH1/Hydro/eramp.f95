!***********************************************************************
!
!  Energy Ramp:  Increases the energy in the zones between the gain 
!  radius and the shock in order to promote an explosion.
!
!  Author:  Ken DeNisco
!           10/23/03   FAU
!
!***********************************************************************
SUBROUTINE eramp( e_bomb, bomb_time, dtime, jexpl_min, jexpl_max,  &
& ij_ray, ik_ray )

USE kind_module
USE numerical_module, ONLY : zero

USE edit_module,      ONLY : nlog
USE mdl_cnfg_module,  ONLY : rho, t, ye, dmrst
USE cycle_module,     ONLY : ncycle
USE incrmnt_module,   ONLY : dtmpmn
USE mdl_cnfg_module,  ONLY : jr_max

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jexpl_min   ! gain radius zone index
INTEGER, INTENT(in)              :: jexpl_max   ! shock zone index
INTEGER, INTENT(in)              :: ij_ray      ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray      ! index denoting the k-index of a specific radial ray

REAL(KIND=double), INTENT(in)    :: e_bomb      ! energy increment
REAL(KIND=double), INTENT(in)    :: bomb_time   ! total time interval
REAL(KIND=double), INTENT(in)    :: dtime       ! time increment

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j           ! zone index
INTEGER                          :: ivar = 2    ! EOS energy index
INTEGER                          :: jjshockmx   ! maximum radial zone index
INTEGER                          :: jjshockmn   ! minimum radial zone index
INTEGER                          :: jbomb_max   ! out edge of thermal bomb

REAL(KIND=double)                :: e_incrmnt   ! energy increment for this cycle
REAL(KIND=double)                :: mass        ! total mass between gain radius and shock
REAL(KIND=double)                :: t_out       ! temperature after increase
REAL(KIND=double)                :: e_intl      ! internal energy before increment
REAL(KIND=double)                :: e_finl      ! internal energy after increment
REAL(KIND=double)                :: dedd        ! derivative of the energy wrt density
REAL(KIND=double)                :: dedt        ! derivative of the temperature wrt density
REAL(KIND=double)                :: dedy        ! derivative of the electron fraction wrt density

REAL(KIND=double)                :: pqmin = 0.5d0 ! minimum value of pseudoviscosity denoting shock

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  If jexpl_max < 0, set jexpl_max at inner edge of shock
!-----------------------------------------------------------------------

IF ( jexpl_max < 0 ) THEN
  CALL findshock_min( jexpl_min, jr_max, ij_ray, ik_ray, pqmin, jjshockmn, &
&  jjshockmx )
  IF ( jjshockmn == jjshockmx ) THEN
    jbomb_max        = ABS(jexpl_max)
  ELSE
    jbomb_max        = jjshockmn
  END IF ! jjshockmn = jjshockmx
ELSE
  jbomb_max          = jexpl_max
END IF ! jexpl_max < 0
WRITE (nlog,3001) jbomb_max
 3001 FORMAT (' jbomb_max=',i3)

!-----------------------------------------------------------------------
!  Compute e_increment, the energy per unit mass added to the model
!   between jexpl_min and jexpl_max for this cycle
!-----------------------------------------------------------------------

mass                 = zero
DO j = jexpl_min, jbomb_max
  mass               = mass + dmrst(j)
END DO

e_incrmnt            = ( e_bomb / mass ) * (dtime / bomb_time)
      
!-----------------------------------------------------------------------
!  Add e_increment to the internal energy and compute the resulting
!   change in temperature
!-----------------------------------------------------------------------

DO j = jexpl_min, jbomb_max
  CALL eqstt_x( ivar, j, ij_ray, ik_ray, rho(j), t(j), ye(j), e_intl, &
&  dedd, dedt, dedy )
  e_finl             = e_intl + e_incrmnt
  CALL tgvndeye_x( j, ij_ray,ik_ray, rho(j), e_finl, ye(j), t(j), t_out )
  dtmpmn(j,10,ij_ray,ik_ray) = t_out - t(j)
END DO
                
RETURN
END SUBROUTINE eramp
