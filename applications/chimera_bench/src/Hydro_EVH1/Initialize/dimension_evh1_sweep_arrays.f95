SUBROUTINE dimension_evh1_sweep_arrays( max_12 )
!-----------------------------------------------------------------------
!
!    File:         dimension_evh1_sweep_arrays
!    Module:       dimension_evh1_sweep_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate the dimensions certain of the evh1 arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  max_12    : max(nx,ny,nz)+12
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, evh1_sweep, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero, one

USE edit_module, ONLY : nlog
USE evh1_sweep
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: max_12        ! max(nx,ny,nz)+12

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

INTEGER                           :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' EVH1 sweep arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_evh1_sweep_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!           \\\\\ ALLOCATE EVH!_SWEEP_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Logical variables indicating the presence of a shock and the
!  background density relative to a present value
!-----------------------------------------------------------------------

ALLOCATE (l_shock(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'l_shock   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (l_rho(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'l_rho     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
! State variables
!-----------------------------------------------------------------------

ALLOCATE (r(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (p_mat(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_mat     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (p_nu(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_nu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (p(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xa(max_12+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xa        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xa0(max_12+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xa0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dx(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dx        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dx0(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dx0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol0(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
! EVH1 additions
!-----------------------------------------------------------------------

ALLOCATE (temp(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'temp      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (entrop(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'entrop    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ei(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_nu(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_v(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_v       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ge(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ge        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gc(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gc        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xn(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xp(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xp        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flat_s(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_s    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
! MGFLD additions
!-----------------------------------------------------------------------

ALLOCATE (rhobar(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhobar    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (egrav(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'egrav     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (degrav(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'degrav    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (egrav0(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'egrav0    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ekin(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ekin      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_strc(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_strc   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_stre(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_stre   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (g_force_c(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'g_force_c '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (g_pot_c(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'g_pot_c   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (g_force_e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'g_force_e '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (g_pot_e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'g_pot_e   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_nu_c(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu_c    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f_nu_e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f_nu_e    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (lapse_c(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'lapse_c   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (lapse_e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'lapse_e   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u_edge(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_edge    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!           \\\\\ INITIALIZE EVH!_SWEEP_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! State variables
!-----------------------------------------------------------------------

r                         = zero
p_mat                     = zero
p_nu                      = zero
p                         = zero
e                         = zero
u                         = zero
v                         = zero
w                         = zero
xa                        = zero
xa0                       = zero
dx                        = zero
dx0                       = zero
dvol                      = zero
dvol0                     = zero

!-----------------------------------------------------------------------
! EVH1 additions
!-----------------------------------------------------------------------

temp                      = zero
entrop                    = zero
ei                        = zero
e_nu                      = zero
e_v                       = zero
ge                        = zero
gc                        = zero
ye                        = zero
xn                        = zero
xp                        = zero
flat_s                    = zero

!-----------------------------------------------------------------------
! MGFLD additions
!-----------------------------------------------------------------------

rhobar                    = zero
egrav                     = zero
degrav                    = zero
egrav0                    = zero
ekin                      = zero
nu_strc                   = zero
nu_stre                   = zero
g_force_c                 = zero
g_pot_c                   = zero
g_force_e                 = zero
g_pot_e                   = zero
e_nu_c                    = zero
f_nu_e                    = zero
lapse_c                   = one
lapse_e                   = one
u_edge                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_evh1_sweep_arrays
                                                                                
