SUBROUTINE dimension_nucbrn_arrays( nx, nnc, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_nucbrn_arrays
!    Module:       dimension_nucbrn_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the nuclear arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array dimension
!  nnc        : composition array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, nucbrn_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE nucbrn_module
USE parallel_module, ONLY : myid, ierr

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: nnc           ! abundance array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

  101 FORMAT (' Nucburn arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_nucbrn_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE NUCBRN_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Nuclear statistical equilibrium flag
!-----------------------------------------------------------------------

ALLOCATE (nse(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Temperature after nuclear burn
!-----------------------------------------------------------------------

ALLOCATE (T_burn(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'T_burn    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Temperature change due to nuclear burning
!-----------------------------------------------------------------------

ALLOCATE (dt_burn(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dT_burn   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Abundance parameters and burning energies
!-----------------------------------------------------------------------

ALLOCATE (a_name(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_name    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (xn(nx,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dudt_nuc(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_nuc  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uburn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc_rep(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_rep(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_rep(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (m_ex_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm_ex_nuc  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Screening corrections
!-----------------------------------------------------------------------

ALLOCATE (fescrn(nx,nnc,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fescrn    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fascrn(nx,nnc,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fascrn    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrnc(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrnc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrnc(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrnc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrnc(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrnc     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Eos cube index
!-----------------------------------------------------------------------

ALLOCATE (idty(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idty      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Reaction rate table
!-----------------------------------------------------------------------

ALLOCATE (rrdata(17,nx,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rrdata    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE NUCBRN_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

inuc                 = 0
nprint               = 0
ncycle               = 0
nse                  = 1
jdTmax               = 0
jdynmax              = 0

dTmax                = zero
dynmax               = zero
t_cntl_burn          = zero
dtime_burn           = zero
T_burn               = zero
dT_burn              = zero

xn                   = zero
dudt_nuc             = zero
uburn                = zero
be_nuc_rep           = zero
a_nuc_rep            = 56.d00
z_nuc_rep            = 28.d00
a_nuc                = zero
z_nuc                = zero
m_ex_nuc             = zero
be_nuc               = zero

itnuc                = 0
ttolnuc              = zero
ytolnuc              = zero
ynmin                = zero

rdpg                 = zero
rhegp                = zero
r3a                  = zero
rg3a                 = zero
rcag                 = zero
roga                 = zero
roag                 = zero
rnega                = zero
rneag                = zero
rmgga                = zero
rmgag                = zero
rsiga                = zero
rcaag                = zero
rtiga                = zero
r1212                = zero
r1216                = zero
r1616                = zero

fescrn               = zero
fascrn               = zero

idrnc                = 0
itrnc                = 0
iyrnc                = 0

idty                 = 0
dgrid                = zero
tgrid                = zero
ygrid                = zero
rhoes                = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_nucbrn_arrays
