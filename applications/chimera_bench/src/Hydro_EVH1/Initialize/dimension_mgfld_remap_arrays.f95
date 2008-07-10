SUBROUTINE dimension_mgfld_remap_arrays( max_12, nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_mgfld_remap_arrays
!    Module:       dimension_mgfld_remap_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the mgfld remap arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x (radial) array dimension
!  nez       : neutrino energy array dimension
!  nnu       : neutrino flavor array dimension
!  nnc       : composition array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, mgfld_remap_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE mgfld_remap_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: max_12        ! max(nx,ny,nz)+12
INTEGER, INTENT(in)              :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)              :: nnc           ! abundance array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' MGLFD remap arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in mgfld_remap_module_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ ALLOCATE MGFLD_REMAP_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initial values
!-----------------------------------------------------------------------

ALLOCATE (r0i(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r0i       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvoli(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvoli     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmi(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmi       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t0i(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't0i       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ei0i(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei0i      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dei(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dei       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Values after lagrangian step
!-----------------------------------------------------------------------

ALLOCATE (r0l(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r0l       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvoll(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvoll     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dxl(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dxl       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xal(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xal       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dml(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dml       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u0l(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u0l       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t0l(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't0l       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ei0l(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei0l      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye0l(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye0l      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Final values
!-----------------------------------------------------------------------

ALLOCATE (r0e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r0e       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvole(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvole     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dxe(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dxe       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xae(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xae       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dme(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dme       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u0e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u0e       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t0e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't0e       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ei0e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei0e      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye0e(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye0e      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Working psi0 array
!-----------------------------------------------------------------------

ALLOCATE (psi0_re(max_12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_re   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Working composition arrays
!-----------------------------------------------------------------------

ALLOCATE (comp(max_12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Working coordinate and therodynamic state arrays
!-----------------------------------------------------------------------

ALLOCATE (r(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (temp(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'temp      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xa(max_12+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xa        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dx(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dx        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xa0(max_12+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xa0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dx0(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dx0       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Binding energy arrays
!-----------------------------------------------------------------------

ALLOCATE (eb(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eb        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxbe(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxbe    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_bind_zn0(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_bind_zn0 '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Composition ye flux arrays
!-----------------------------------------------------------------------

ALLOCATE (fluxye_comp(max_12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxye_comp'; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ INITIALIZE MGFLD_REMAP_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initial values
!-----------------------------------------------------------------------

r0i                       = zero
dvoli                     = zero
dmi                       = zero
t0i                       = zero
ei0i                      = zero
dei                       = zero

!-----------------------------------------------------------------------
!  Values after lagrangian step
!-----------------------------------------------------------------------

r0l                       = zero
dvoll                     = zero
dxl                       = zero
xal                       = zero
dml                       = zero
u0l                       = zero
t0l                       = zero
ei0l                      = zero
ye0l                      = zero

!-----------------------------------------------------------------------
!  Final values
!-----------------------------------------------------------------------

r0e                       = zero
dvole                     = zero
dxe                       = zero
xae                       = zero
dme                       = zero
u0e                       = zero
t0e                       = zero
ei0e                      = zero
ye0e                      = zero

!-----------------------------------------------------------------------
!  Working psi0 array
!-----------------------------------------------------------------------

psi0_re                   = zero

!-----------------------------------------------------------------------
!  Working composition arrays
!-----------------------------------------------------------------------

comp                      = zero

!-----------------------------------------------------------------------
!  Working coordinate and therodynamic state arrays
!-----------------------------------------------------------------------

r                         = zero
xa                        = zero
dx                        = zero
xa0                       = zero
dx0                       = zero

!-----------------------------------------------------------------------
!  Binding energy arrays
!-----------------------------------------------------------------------

eb                        = zero
fluxbe                    = zero
e_bind_zn0                = zero

!-----------------------------------------------------------------------
!  Composition ye flux arrays
!-----------------------------------------------------------------------

fluxye_comp               = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_mgfld_remap_arrays
