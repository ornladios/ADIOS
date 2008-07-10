SUBROUTINE dimension_prb_cntl_arrays( nnu )
!-----------------------------------------------------------------------
!
!    File:         dimension_prb_cntl_arrays
!    Module:       dimension_prb_cntl_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/04
!
!    Purpose:
!      To allocate the dimensions of the prb_cntl arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nnu       : neutrino flavor array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module, prb_cntl_module
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE prb_cntl_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Radhyd problem control arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_prb_cntl_arrays')

!-----------------------------------------------------------------------
!        Allocate prb_cntl_module arrays
!-----------------------------------------------------------------------

ALLOCATE (itnu(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itnu      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE CONTROLS /////
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!  Hydrodynamics controls
!-----------------------------------------------------------------------

irelhy                    = 0
ihydro                    = 0
lagr                      = 'ye'
rezn                      = 'no'

!-----------------------------------------------------------------------
!  Transport controls
!-----------------------------------------------------------------------

inutrn                    = 0
ireltrns                  = 0
idiff                     = 0
jnumin                    = 0
jnumax                    = 0
itnu                      = 1
jexect                    = 0
iyenu                     = 1
jexecy                    = 0

!-----------------------------------------------------------------------
!  Absorption and emission controls
!-----------------------------------------------------------------------

iaefnp                    = 0
i_aeps                    = 0
iaence                    = 0
iaenca                    = 0
iaencnu                   = 0
iaenct                    = 0

rhoaefnp                  = 1.d+15
edmpe                     = zero
edmpa                     = zero
roaencnu                  = 1.d+15
roaenct                   = 1.d+13

!-----------------------------------------------------------------------
!  Charged current weak magnetism switch.
!-----------------------------------------------------------------------

icc_wm                    = 1

!-----------------------------------------------------------------------
!  Master scattering switch
!-----------------------------------------------------------------------

iscat                     = 1

!-----------------------------------------------------------------------
!  Neutrino nucleon scattering controls
!-----------------------------------------------------------------------

in                        = 0
ip                        = 0
ietann                    = 0

!-----------------------------------------------------------------------
!  Neutrino nuclei scattering controls
!-----------------------------------------------------------------------

ihe                       = 0
iheavy                    = 0
iicor                     = 0

!-----------------------------------------------------------------------
!  Neutral current weak magnetism switch.
!-----------------------------------------------------------------------

inc_wm                    = 1

!-----------------------------------------------------------------------
!  Neutrino-electron scattering controls
!-----------------------------------------------------------------------

nes                       = 0
rhonesmn                  = 1.d+7
rhonesmx                  = 1.d+15

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering controls
!-----------------------------------------------------------------------

nncs                      = 0
rhonncsmn                 = 1.d+7
rhonncsmx                 = 1.d+15

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering controls
!-----------------------------------------------------------------------

isctn                     = 0
rhosctnemn                = 1.d+7
rhosctnemx                = 1.d+15
rhosctntmn                = 1.d+7
rhosctntmx                = 1.d+15

!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic scattering controls
!-----------------------------------------------------------------------

isctnn                    = 0
rhosctnnemn               = 1.d+7
rhosctnnemx               = 1.d+15
rhosctnntmn               = 1.d+7
rhosctnntmx               = 1.d+15

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation controls
!-----------------------------------------------------------------------

ipair                     = 0
rhopairemn                = 1.d+7
rhopairemx                = 1.d+15
rhopairtmn                = 1.d+7
rhopairtmx                = 1.d+15

!-----------------------------------------------------------------------
!  Nuclear excitation pair annihilation controls
!-----------------------------------------------------------------------

ipairA                    = 0
rhopairAemn               = 1.d+7
rhopairAemx               = 1.d+15
rhopairAtmn               = 1.d+7
rhopairAtmx               = 1.d+15

!-----------------------------------------------------------------------
!  Bremsstrahlung pair annihilation controls
!-----------------------------------------------------------------------

ibrem                     = 0
rhobrememn                = 1.d+7
rhobrememx                = 1.d+15
rhobremtmn                = 1.d+7
rhobremtmx                = 1.d+15

!-----------------------------------------------------------------------
!  Convection controlsls
!-----------------------------------------------------------------------

ilcnvct                   = 0
iemix                     = 0
iyemix                    = 0
ipsimix                   = 0
ipcnvct                   = 0
iscnvct                   = 0

!-----------------------------------------------------------------------
!  NSE flashing and deflashing controls
!-----------------------------------------------------------------------

tnse                      = 1.d+15
tdnse                     = 1.d+00

!-----------------------------------------------------------------------
!  Explosion simulation control
!-----------------------------------------------------------------------

i_bomb                    = 0

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_prb_cntl_arrays
