SUBROUTINE j_rho(jmn,jmx,rho_t)
!-----------------------------------------------------------------------
!
!    File:         j_rho
!    Module:       j_rho
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/07/02
!
!    Purpose:
!      To determine jmn and jmx, the mass zone index on either side of
!       a given density, rho_t.
!
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  rho_t     : given value of the mass density (g/cm3)
!
!    Output arguments:
!  jmn       : zone index immediately interior to rho_t
!  jmx       : zone index immediately exterior to rho_t
!
!    Input arguments (common):
!        none
!
!    Output arguments (common):
!        none
!
!    Include files:
!  kind_module
!  mdl_cnfg_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE mdl_cnfg_module, ONLY : jr_max, rho

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: rho_t         ! given value of the mass density (g/cm3)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out)             :: jmn           ! radial zone index just interior to rho_t
INTEGER, INTENT(out)             :: jmx           ! radial zone index just exterior to rho_t

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine j such that rho(j) < rho_t and rho(j-1) > rho_t
!-----------------------------------------------------------------------

DO j = 2,jr_max
  IF ( rho(j) <= rho_t  .and.  rho(j-1) >= rho_t ) THEN
    jmn            = j - 1
    jmx            = j
    RETURN
  END IF
END DO

!-----------------------------------------------------------------------
!  If none found, set jmn = jmx = 0
!-----------------------------------------------------------------------

jmn                = 0
jmx                = 0

RETURN
END SUBROUTINE j_rho
