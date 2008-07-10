SUBROUTINE dimension_e_advct_arrays( nx, ny, nz, nez, nnu, ij_ray_dim, &
& ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_e_advct_arrays
!    Module:       dimension_e_advct_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate the dimensions to certain of the energy advection arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array dimension
!  ny         : y (angular) array dimension
!  nz         : z (azimuthal) array dimension
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  e_advct_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE e_advct_module
USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=15)               :: var_name

INTEGER                          :: istat         ! allocation status
INTEGER                          :: ndim          ! array extent

!-----------------------------------------------------------------------
!        Formate
!-----------------------------------------------------------------------

  101 FORMAT (' e_advection arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a15,' in dimension_e_advct_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set ndim to the maximum of nx and ny
!-----------------------------------------------------------------------

ndim                      = MAX( nx, ny, nz )

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE E_ADVCT_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------

ALLOCATE (jdt_nu_e_advct(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'jdt_nu_e_advct '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_cntrl_e_advct(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_cntrl_e_advct'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi0nmax(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi0nmax      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psivmin(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psivmin        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsivmx(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsivmx        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dt_nu_e_advct(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dt_nu_e_advct  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Quantities needed to compute neutrino energy advection rate
!-----------------------------------------------------------------------

ALLOCATE (rjmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rjmh           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rajmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rajmh          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (yjmh(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yjmh           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (yajmh(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yajmh          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zjmh(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zjmh           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zajmh(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zajmh          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (adot(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'adot           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (adot_a(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'adot_a         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (bdot_b(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'bdot_b         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ddot(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddot           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ddot_d(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddot_d         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rdot(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rdot           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rdot_r(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rdot_r         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v_e1(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_e1           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v_e2(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_e2           '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (v_e(ndim,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_e            '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi(nez+1+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xk(nez+1+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xk             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dk(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dk             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xk_0(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xk_0           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dk_0(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dk_0           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xk_1(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xk_1           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dk_1(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dk_1           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (vk(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vk             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvolk(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvolk          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvolk_0(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvolk_0        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvolk_1(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvolk_1        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Quantities derived from the neutrino energy zoning
!-----------------------------------------------------------------------

ALLOCATE (ncoefa(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ncoefa         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ncoefaa(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ncoefaa        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ecoefa(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefa         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ecoefaa(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefaa        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  State variables
!-----------------------------------------------------------------------

ALLOCATE (u(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u              '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rho(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rhoa(ndim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhoa           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r              '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ra(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ra             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (y(ny+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y              '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ya(ny+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ya             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z(nz+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z              '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (za(nz+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'za             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrjmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrjmh         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrajmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrajmh        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rstmss(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rstmss         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmrst(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmrst          '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino occupation distribution
!-----------------------------------------------------------------------

ALLOCATE (psi0(ndim,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1(ndim,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1           '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi0_a(ndim,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_a         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flxf(nx,nez+1,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flxf           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (Ef(nx,nez+1,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'Ef             '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Energy advection parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

ALLOCATE (unucrv(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unucrv         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunujvdt(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunujvdt       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujv(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujv          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dndt_v(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dndt_v         '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE E_ADVCT_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------

jdt_nu_e_advct            = 0
t_cntrl_e_advct           = 1.d-01
dpsi0nmax                 = zero
psivmin                   = 1.0d-01
dpsivmx                   = zero
dt_nu_e_advct             = 1.0d+20

!-----------------------------------------------------------------------
!  Quantities needed to compute neutrino energy advection rate
!-----------------------------------------------------------------------

rjmh                      = zero
rajmh                     = zero
yjmh                      = zero
yajmh                     = zero
adot                      = zero
adot_a                    = zero
bdot_b                    = zero
ddot                      = zero
ddot_d                    = zero
rdot                      = zero
rdot_r                    = zero
v_e1                      = zero
v_e2                      = zero
v_e                       = zero

psi                       = zero
xk                        = zero
dk                        = zero
xk_0                      = zero
dk_0                      = zero
xk_1                      = zero
dk_1                      = zero
vk                        = zero
dvolk                     = zero
dvolk_0                   = zero
dvolk_1                   = zero
  
!-----------------------------------------------------------------------
!  Quantities derived from the neutrino energy zoning
!-----------------------------------------------------------------------

ncoefa                    = zero
ncoefaa                   = zero
ecoefa                    = zero
ecoefaa                   = zero

!-----------------------------------------------------------------------
!  State variables
!-----------------------------------------------------------------------

u                         = zero
rho                       = zero
rhoa                      = zero
r                         = zero
ra                        = zero
y                         = zero
ya                        = zero
z                         = zero
za                        = zero
agrjmh                    = zero
agrajmh                   = zero
rstmss                    = zero
dmrst                     = zero

!-----------------------------------------------------------------------
!  Neutrino occupation distribution
!-----------------------------------------------------------------------

psi0                      = zero
psi1                      = zero
psi0_a                    = zero
flxf                      = zero
Ef                        = zero

!-----------------------------------------------------------------------
!  Energy advection parameters derived from psi0 and psi1
!-----------------------------------------------------------------------

unucrv                    = zero
dunujvdt                  = zero
unujv                     = zero
dndt_v                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_e_advct_arrays
