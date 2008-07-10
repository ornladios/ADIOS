SUBROUTINE nu_U( jr_min, jr_max, rho, nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         nu_U
!    Module:       nu_U
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/16/01
!
!    Purpose:
!      To calculate at the current time step the neutrino energy per 
!       volume and per unit mass, and the isotropic pressure and flux.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  jr_min      : inner zone for which calculation w is to be made
!  jr_max      : outer zone for which calculation w is to be made
!  nx          : x-array extent
!  nnu         : neutrino flavor aray extent
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!  nnugp(n)    : number of neutrino energy groups for n-neutrinos
!  ecoefa(j,k) : energy sum neutrino states per unit volume at radial zone j-1/2 energy unu(j,k)
!  psi0(j,k,n) : zeroth angular moment of the neutrino distribution function
!
!    Output arguments (common):
!  aunu(j)     :  neutrino energy/mass at time n (ergs/g)
!  apnu(j)     :  neutrino pressure (dynes/cm2)
!  fnu(j)      :  neutrino energy flux at time t (ergs/cm2*sec)
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, third, half

USE e_advct_module, ONLY : Ef
USE nu_dist_module, ONLY : ecoefa, stwt, psi0, aunu, apnu, fnu, fluxnu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor aray extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)  :: rho   ! density (g cm^{-3})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  aunu(j)          = zero
  apnu(j)          = zero
  fnu(j)           = zero
END DO

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!  Neutrino energy/volume (ergs/cm3)
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max
      aunu(j)      = aunu(j) + ecoefa(j,k) * stwt(n) * psi0(j,k,n)
      apnu(j)      = apnu(j) + ecoefa(j,k) * stwt(n) * psi0(j,k,n) * Ef(j,k,n)
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Neutrino energy/mass (ergs/gm)
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  aunu(j)          = aunu(j)/rho(j)
END DO

!-----------------------------------------------------------------------
!  Neutrino energy flux at current time step (ergs/cm2*sec)
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  fnu(j)           = half * ( SUM(fluxnu(j,:)) + SUM(fluxnu(j-1,:)) )
END DO

RETURN
END SUBROUTINE nu_U
