!-----------------------------------------------------------------------
!    Module:       mdl_cnfg_y_module
!    Author:       S. W. Bruenn
!    Date:         3/20/05
!-----------------------------------------------------------------------

MODULE mdl_cnfg_y_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Radial mass zoning
!-----------------------------------------------------------------------
!
!  ja_min : inner angular zone index.
!  ja_max : outer angular zone index.
!-----------------------------------------------------------------------

INTEGER                                        :: ja_min
INTEGER                                        :: ja_max


!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------
!
!  rho(j) : the density of angular zone j at timestep m (current time)
!   (g/cm**3).
!  t(j) : the temperature of angular zone j at timestep m (K).
!  ye(j) : the electron fraction of angular zone j at timestep m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: ye


END module mdl_cnfg_y_module
