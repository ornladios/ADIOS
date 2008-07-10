SUBROUTINE dimension_mdl_cnfg_arrays( nx )
!-----------------------------------------------------------------------
!
!    File:         dimension_mdl_cnfg_arrays
!    Module:       dimension_mdl_cnfg_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the model configuration arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x (radial) array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, mdl_cnfg_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero, one

USE edit_module, ONLY : nlog
USE mdl_cnfg_module
USE parallel_module, ONLY : myid

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER, INTENT(in)              :: nx            ! radial array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Model configuration arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_mdl_cnfg_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE MDL_CNFG_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------

ALLOCATE (rho(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rhor(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhor      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tr        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (yer(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yer       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rho0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (p0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p0        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Mechanical state
!-----------------------------------------------------------------------

ALLOCATE (u(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Zone masses
!-----------------------------------------------------------------------

ALLOCATE (rstmss(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rstmss    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmrst(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmrst     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (grvmss(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grvmss    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmgrv(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmgrv     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  General relativistic parameters
!-----------------------------------------------------------------------

ALLOCATE (gamgr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gamgr     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrr      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrh      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agra(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agra      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrah(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrah     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wgrr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wgrr      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wgr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wgr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gamgrr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gamgrr    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gamgra(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gamgra    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (grvmssa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grvmssa   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmgrva(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmgrva    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gamgr_i(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gamgr_i   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agr_i(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_i     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wgr_i(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wgr_i     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (grvmss_i(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grvmss_i  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmgrv_i(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmgrv_i   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Specific energy density arrays for rezoning
!-----------------------------------------------------------------------

ALLOCATE (uad(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uad       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (utoti(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'utoti     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (utotf(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'utotf     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE MDLCNFG_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------

rho                       = zero
t                         = zero
ye                        = zero

rhor                      = zero
tr                        = zero
yer                       = zero

rho0                      = zero
ye0                       = zero
p0                        = zero

!-----------------------------------------------------------------------
!  Mechanical state
!-----------------------------------------------------------------------

u                         = zero
v                         = zero
w                         = zero
r                         = zero
dr                        = zero

!-----------------------------------------------------------------------
!  Zone masses
!-----------------------------------------------------------------------

rstmss                    = zero
dmrst                     = zero
grvmss                    = zero
dmgrv                     = zero

!-----------------------------------------------------------------------
!  General relativistic parameters
!-----------------------------------------------------------------------

gamgr                     = one
agrr                      = one
agr                       = one
agrh                      = one
wgrr                      = zero
wgr                       = zero
gamgrr                    = zero
gamgra                    = zero
grvmssa                   = zero
dmgrva                    = zero
gamgr_i                   = zero
agr_i                     = zero
wgr_i                     = zero
grvmss_i                  = zero
dmgrv_i                   = zero

!-----------------------------------------------------------------------
!  Specific energy density arrays for rezoning
!-----------------------------------------------------------------------

uad                       = zero
utoti                     = zero
utotf                     = zero

WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_mdl_cnfg_arrays
