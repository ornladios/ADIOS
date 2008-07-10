!-----------------------------------------------------------------------
!    Module:       mdl_cnfg_z_module
!    Author:       S. W. Bruenn
!    Date:         3/20/05
!-----------------------------------------------------------------------

MODULE mdl_cnfg_z_module

USE kind_module

SAVE



!-----------------------------------------------------------------------
!  Azimuthal zone indices
!-----------------------------------------------------------------------
!
!  jz_min : inner azimuthal zone index.
!  jz_max : outer azimuthal zone index.
!-----------------------------------------------------------------------

INTEGER                                        :: jz_min
INTEGER                                        :: jz_max


!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------
!
!  rho(j) : the density of azimuthal zone j at timestep m (current time) (g/cm**3).
!  t(j)   : the temperature of azimuthal zone j at timestep m (K).
!  ye(j)  : the electron fraction of azimuthal zone j at timestep m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: ye


END module mdl_cnfg_z_module
