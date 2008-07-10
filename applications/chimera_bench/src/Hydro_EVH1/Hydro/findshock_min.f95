SUBROUTINE findshock_min( jr_min, jr_max, ij_ray, ik_ray, pqmin, jjshockmn, &
& jjshockmx )
!-----------------------------------------------------------------------
!
!    File:         findshock_min
!    Module:       findshock_min
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/09/03
!
!    Purpose:
!      To find the inner and outer zones of the minimum radius shock.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min      : inner radial zone to end search 
!  jr_max      : outer radial zone to start the search
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!  pqmin       : minimum ratio of pseudoviscosity to pressure used
!                 to infer the presence of a shock
!
!    Output arguments:
!  jjshockmx   : outer radial zone of shock
!  jjshockmn   : inner radial zone of shock
!
!    Include files:
!  kind_module
!  eos_snc_x_module, shock_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE eos_snc_x_module, ONLY : aesv
USE shock_module, ONLY : pq_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)    :: pqmin         ! minimum ratio of pseudoviscosity to pressure used

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out)             :: jjshockmx     ! maximum radial zone index
INTEGER, INTENT(out)             :: jjshockmn     ! minimum radial zone index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find jjshockmn
!-----------------------------------------------------------------------

DO j = jr_min+1,jr_max
  jjshockmn         = j - 1
  IF ( pq_x(j,ij_ray,ik_ray)/aesv(j,1,ij_ray,ik_ray) >= pqmin ) EXIT
END DO

IF ( jjshockmn == jr_max - 1 ) THEN
  jjshockmx         = jjshockmn
  RETURN
END IF ! jjshockmn = jr_max - 1

!-----------------------------------------------------------------------
!  Find jjshockmx
!-----------------------------------------------------------------------

DO j = jjshockmn+1,jr_max
  jjshockmx         = j + 1
  IF ( pq_x(j,ij_ray,ik_ray)/aesv(j,1,ij_ray,ik_ray) < pqmin ) RETURN
END DO

RETURN
END SUBROUTINE findshock_min
