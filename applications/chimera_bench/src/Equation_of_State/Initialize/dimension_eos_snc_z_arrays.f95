SUBROUTINE dimension_eos_snc_z_arrays( nz, ij_ray_dim, k_ray_dim, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_eos_snc_z_arrays
!    Module:       dimension_eos_snc_z_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To allocate dimensions to certain of the eos_snc_z arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nz         : z-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!  nnc        : composition aray extent
!
!    Output arguments:
!        none
!
!    Include files:
!      numerical_module
!      edit_module, eos_snc_z_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE eos_snc_z_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nz            ! z (azimuthal) array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: nnc           ! abundance array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Equation of state arrays for z (azimuthal) rays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_eos_snc_z_arrays')

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

ALLOCATE (nse(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated equation of state variables
!-----------------------------------------------------------------------

ALLOCATE (aesv(nz,12,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesv      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvd(nz,12,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvt(nz,12,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvt     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvy(nz,12,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvy     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (gam1(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gam1      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gam2(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gam2      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gam3(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gam3      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (duesrc(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'duesrc    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (idty(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idty      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idr(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itr(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyr(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyr       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  EOS table
!-----------------------------------------------------------------------

ALLOCATE (estble(12,nz,2,2,2,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'estble    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (escnst(12,nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'escnst    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------

ALLOCATE (xn(nz,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (be_nuc_rep(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_rep(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_rep(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!           \\\\\ INITIALIZE EOS_SNC_Y_MODULE ARRAYS /////
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
END SUBROUTINE dimension_eos_snc_z_arrays
