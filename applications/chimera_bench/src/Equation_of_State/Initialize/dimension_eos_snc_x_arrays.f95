SUBROUTINE dimension_eos_snc_x_arrays( nx, ij_ray_dim, ik_ray_dim, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_eos_snc_x_arrays
!    Module:       dimension_eos_snc_x_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to eos_snc_x_module arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  nnc        : composition aray extent
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, eos_snc_x_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE eos_snc_x_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nnc           ! abundance array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Equation of state arrays for radial rays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_eos_snc_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ ALLOCATE EOS_SNC_X_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Nuclear statistical equilibrium
!-----------------------------------------------------------------------

ALLOCATE (nse(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated equation of state variables
!-----------------------------------------------------------------------

ALLOCATE (aesv(nx,12,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesv      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvd(nx,12,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvt(nx,12,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvt     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvy(nx,12,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvy     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (gam1(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gam1      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gam2(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gam2      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gam3(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gam3      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (duesrc(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'duesrc    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (idty(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idty      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idr(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itr(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyr(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyr       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  EOS table
!-----------------------------------------------------------------------

ALLOCATE (estble(12,nx,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'estble    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (escnst(12,nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'escnst    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------

ALLOCATE (a_name(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_name    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (xn(nx,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,1001) var_name; END IF

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
ALLOCATE (be_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!           \\\\\ INITIALIZE EOS_SNC_X_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Nuclear statistical equilibrium
!-----------------------------------------------------------------------

nse                       = 1

!-----------------------------------------------------------------------
!  Interpolated equation of state variables
!-----------------------------------------------------------------------

aesv                      = zero
aesvd                     = zero
aesvt                     = zero
aesvy                     = zero
gam1                      = zero
gam2                      = zero
gam3                      = zero
duesrc                    = zero
idty                      = 0

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idr                       = 0
itr                       = 0
iyr                       = 0

!-----------------------------------------------------------------------
!  EOS table
!-----------------------------------------------------------------------

estble                    = zero
escnst                    = zero

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------

xn                        = zero
be_nuc_rep                = 492.3d0
a_nuc_rep                 = 56.d0
z_nuc_rep                 = 28.d0
a_nuc                     = zero
z_nuc                     = zero
be_nuc                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_eos_snc_x_arrays
