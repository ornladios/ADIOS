SUBROUTINE dimension_nu_dist_arrays( nx, ny, nz, nez, nnu, nnc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_nu_dist_arrays
!    Module:       dimension_nu_dist_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the neutrino distribution arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x-arrays dimension
!  ny         : y (angular) array extent
!  nz         : z (azimuthal) array extent
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  nnc        : composition array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  j_ray_dim  : number of radial zones on a processor after swapping with y
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, nu_dist_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero, one

USE edit_module, ONLY : nlog
USE nu_dist_module
USE parallel_module, ONLY : myid

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx         ! radial array extent
INTEGER, INTENT(in)              :: ny         ! y (angular) array extent
INTEGER, INTENT(in)              :: nz         ! z (azimuthal) array extent
INTEGER, INTENT(in)              :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu        ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc        ! composition array dimension
INTEGER, INTENT(in)              :: ij_ray_dim ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: j_ray_dim  ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: k_ray_dim  ! number of radial zones on a processor after swapping with z

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat      ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Neutrino distribution arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_nu_dist_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE NU_DIST_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------

ALLOCATE (unu(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unu       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunu(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unub(nx,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unub      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (unue(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unue      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunue(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunue     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unube(nx,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unube     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Functions of the neutrino group energies
!-----------------------------------------------------------------------

ALLOCATE (ncoefa(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ncoefa    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ecoefa(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefa    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ecoefae(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefae   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Momeonts of the neutrino occupation distribution
!-----------------------------------------------------------------------

ALLOCATE (psi0(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (stwt(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'stwt      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino number density
!-----------------------------------------------------------------------

ALLOCATE (n_nu(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n_nu      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Arrays associated with neutrino diffusion
!-----------------------------------------------------------------------

ALLOCATE (rhs1(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhs1      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tmfp(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tmfp      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tmfp_j(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tmfp_j    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_mfp_inv(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_mfp_inv '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dc(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dcr(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxlm(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxlm    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  State variables used in the transport subroutines
!-----------------------------------------------------------------------

ALLOCATE (rjmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rjmh      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gamgr_nu(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gamgr_nu  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agr_nu(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_nu    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrjmh_nu(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrjmh_nu '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (area(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'area      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (areajmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'areajmh   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (vol(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vol       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (drjmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'drjmh     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (drjmh_inv(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'drjmh_inv '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Thermodynamic quantities derived from the neutrino occupation
!   distribution
!-----------------------------------------------------------------------

ALLOCATE (aunu(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aunu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (apnu(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'apnu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (apnun(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'apnun     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (apnur(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'apnur     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (apmnun(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'apmnun    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (apmnur(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'apmnur    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rhonun(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhonun    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rhonur(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhonur    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino  fluxws
!-----------------------------------------------------------------------

ALLOCATE (fnu(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fnu       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxnuk(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxnuk   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxnu(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxnu    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino  stresses
!-----------------------------------------------------------------------

ALLOCATE (strsnu_x(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'strsnu_x  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (stress_x(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'stress_x  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_str_ex(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_ex '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_str_cx(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_cx '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (strsnu_y(ny,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'strsnu_y  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (stress_y(ny,nnu,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'stress_y  '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (strsnu_z(nz,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'strsnu_z  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (stress_z(nz,nnu,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'stress_z  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Nwutrinospheres
!-----------------------------------------------------------------------

ALLOCATE (j_sphere(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'j_sphere  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r_sphere(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r_sphere  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_sphere(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_sphere  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_sphere(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_sphere  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (m_sphere(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm_sphere  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Transport parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

ALLOCATE (dnurad(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dnurad    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unulum(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unulum    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_rad(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_rad     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unurad(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unurad    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unucr(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucr     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unuinfty(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unuinfty  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unucrt(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrt    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unukrad(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unukrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujrad(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujcr(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujcr    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujinfty(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujinfty '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnulum(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnulum    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (elec_rad(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'elec_rad  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnurad(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnurad    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnucr(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnucr     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnukrad(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnukrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnujrad(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnujrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnujcr(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnujcr    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujtdt(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujtdt  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujt     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dudt_nu(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_nu   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Absorption and emission parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

ALLOCATE (unucrea(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrea   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujeadt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujeadt '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujea(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujea    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-electron scattering parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

ALLOCATE (unucrnis(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrnis  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujnisdt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujnisdt'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujnis(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujnis   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation parameters derived from psi0
!   and psi1
!-----------------------------------------------------------------------

ALLOCATE (unucrpa(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrpa   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujpadt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujpadt '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujpa(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujpa    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Bremsstrahlung parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

ALLOCATE (unucrba(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrba   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujbadt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujbadt '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujba(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujba    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-nucleus elastic scattering parameters derived from psi0
!   and psi1
!-----------------------------------------------------------------------

ALLOCATE (unucrns(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrns   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujnsdt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujnsdt '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujns(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujns    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering parameters derived from psi0
!   and psi1
!-----------------------------------------------------------------------

ALLOCATE (unucrnns(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrnns  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujnnsdt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujnnsdt'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujnns(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujnns   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Parameterization of the neutrino distribution
!-----------------------------------------------------------------------

ALLOCATE (enuta(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'enuta     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (enucpa(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'enucpa    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Composition mass fractions
!-----------------------------------------------------------------------

ALLOCATE (xn(nx,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ INITIALIZE NU_DIST_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------

unu                       = zero
dunu                      = zero
unub                      = zero

unue                      = zero
dunue                     = zero
unube                     = zero

!-----------------------------------------------------------------------
!  Functions of the neutrino group energies
!-----------------------------------------------------------------------

ncoefa                    = zero
ecoefa                    = zero
ecoefae                   = zero

!-----------------------------------------------------------------------
!  Momeonts of the neutrino occupation distribution
!-----------------------------------------------------------------------

psi0                      = zero
psi1                      = zero
stwt                      = zero

!-----------------------------------------------------------------------
!  Neutrino number density
!-----------------------------------------------------------------------

n_nu                      = zero

!-----------------------------------------------------------------------
!  Arrays associated with neutrino diffusion
!-----------------------------------------------------------------------

rhs1                      = zero
tmfp                      = zero
tmfp_j                    = zero
e_mfp_inv                 = zero
dc                        = zero
dcr                       = zero
fluxlm                    = zero

!-----------------------------------------------------------------------
!  State variables used in the transport subroutines
!-----------------------------------------------------------------------

rjmh                      = zero
gamgr_nu                  = one
agr_nu                    = one
agrjmh_nu                 = one
area                      = zero
areajmh                   = zero
vol                       = zero
drjmh                     = zero
drjmh_inv                 = zero

!-----------------------------------------------------------------------
!  Thermodynamic quantities derived from the neutrino occupation
!   distribution
!-----------------------------------------------------------------------

aunu                      = zero
apnu                      = zero
apnun                     = zero
apnur                     = zero
apmnun                    = zero
apmnur                    = zero
rhonun                    = zero
rhonur                    = zero

!-----------------------------------------------------------------------
!  Neutrino  fluxws
!-----------------------------------------------------------------------

fnu                       = zero
fluxnuk                   = zero
fluxnu                    = zero

!-----------------------------------------------------------------------
!  Neutrino  stresses
!-----------------------------------------------------------------------

strsnu_x                  = zero
stress_x                  = zero
strsnu_y                  = zero
stress_y                  = zero

!-----------------------------------------------------------------------
!  Nwutrinospheres
!-----------------------------------------------------------------------

j_sphere                  = 1
r_sphere                  = zero
d_sphere                  = zero
t_sphere                  = zero
m_sphere                  = zero

!-----------------------------------------------------------------------
!  Transport parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

dnurad                    = zero
unulum                    = zero
e_rad                     = zero
unurad                    = zero
unucr                     = zero
unuinfty                  = zero
unucrt                    = zero
unukrad                   = zero
unujrad                   = zero
unujcr                    = zero
unujinfty                 = zero
nnulum                    = zero
elec_rad                  = zero
nnurad                    = zero
nnucr                     = zero
nnukrad                   = zero
nnujrad                   = zero
nnujcr                    = zero
dunujtdt                  = zero
unujt                     = zero
dudt_nu                   = zero

!-----------------------------------------------------------------------
!  Absorption and emission parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

unucrea                   = zero
dunujeadt                 = zero
unujea                    = zero

!-----------------------------------------------------------------------
!  Neutrino-electron scattering parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

unucrnis                  = zero
dunujnisdt                = zero
unujnis                   = zero

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation parameters derived from psi0
!   and psi1
!-----------------------------------------------------------------------

unucrpa                   = zero
dunujpadt                 = zero
unujpa                    = zero

!-----------------------------------------------------------------------
!  Bremsstrahlung parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

unucrba                   = zero
dunujbadt                 = zero
unujba                    = zero

!-----------------------------------------------------------------------
!  Neutrino-nucleus elastic scattering parameters derived from psi0
!   and psi1
!-----------------------------------------------------------------------

unucrns                   = zero
dunujnsdt                 = zero
unujns                    = zero

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering parameters derived from psi0
!   and psi1
!-----------------------------------------------------------------------

unucrnns                  = zero
dunujnnsdt                = zero
unujnns                   = zero

!-----------------------------------------------------------------------
!  Parameterization of the neutrino distribution
!-----------------------------------------------------------------------

enuta                     = zero
enucpa                    = zero

!-----------------------------------------------------------------------
!  Composition mass fractions
!-----------------------------------------------------------------------

xn                        = zero

WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_nu_dist_arrays
