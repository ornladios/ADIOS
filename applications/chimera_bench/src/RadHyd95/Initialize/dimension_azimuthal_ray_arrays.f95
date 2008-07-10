SUBROUTINE dimension_azimuthal_ray_arrays( nz, ij_ray_dim, k_ray_dim, &
& nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_azimuthal_ray_arrays
!    Module:       dimension_azimuthal_ray_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/13/05
!
!    Purpose:
!      To allocate the dimensions of the radhyd angular arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nz         : y (angular) array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  nnc        : composition array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  angular_ray_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE azimuthal_ray_module
USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nz            ! z-array extent
INTEGER, INTENT(in)               :: k_ray_dim     ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)               :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)               :: nnc           ! composition array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=12)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Azimuthal arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_azimuthal_ray_arrays')

!-----------------------------------------------------------------------
!
!                \\\\\ ALLOCATE ANGULAR_RAY_MODULE /////
!                /////            ARRAYS           \\\\\
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!           \\\\\ TIME STEPS AND TIME STEP CONSTRAINTS /////
!
!-----------------------------------------------------------------------

ALLOCATE (jdt_z(3,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'jdt_z       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dt_z(3,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dt_z        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!         \\\\\ PHYSICAL AZIMUTHAL COORDINATES DIFFERENCES /////
!
!-----------------------------------------------------------------------

ALLOCATE (dzphys_c(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dzphys_c    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------

ALLOCATE (rho_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_z       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_z         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_z        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ei_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_z        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (p_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_z         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gc_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gc_z        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ge_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ge_z        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (u_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_z         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_z         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_z         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (w_e(nz+1,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_e         '; WRITE  (*,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - INITIAL VALUES /////
!
!-----------------------------------------------------------------------

ALLOCATE (rho_zi(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_zi      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_zi(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_zi        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_zi(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_zi       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ei_zi(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_zi       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (u_zi(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_zi        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v_zi(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_zi        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w_zi(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_zi        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                     \\\\\ GR VARIABLES /////
!
!-----------------------------------------------------------------------
!  agr_z : lapse function
!-----------------------------------------------------------------------

ALLOCATE (agr_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_z       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!               \\\\\ SHOCK STABILIZING VARIABLES /////
!
!-----------------------------------------------------------------------
!  flat_z   : variables aligned along the y-axis indicating the presence
!   of an angular shock
!
!  flat_x_z : variables aligned along the y-axis indicating the presence
!   of a radial shock
!-----------------------------------------------------------------------

ALLOCATE (flat_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_x      '; WRITE  (*,1001) var_name; END IF
ALLOCATE (flat_x_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_z_x    '; WRITE  (*,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ GRAVITATION ACCELERATIONS /////
!
!-----------------------------------------------------------------------
!  grav_y_cy  : zone-centered y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_y_ey  : zone-edged y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!-----------------------------------------------------------------------

ALLOCATE (grav_z_cz(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_z_cz   '; WRITE  (*,1001) var_name; END IF
ALLOCATE (grav_z_ez(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_z_ez   '; WRITE  (*,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                  ||||| COMPOSITION VARIABLES |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------
!  xn(j,nc,i)      : mass fraction of the ith nucleus.
!
!  uburn(j,i)      : cumulative energy generated in zone j by nuclear
!   reactions (ergs/gm).
!
!  be_nuc_rep(j,i) : binding energy of the representative heavy nucleus
!   (MeV).
!
!  a_nuc_rep(j,j)  : mass number of the representative heavy nucleus.
!
!  z_nuc_rep(j,i)  : charge number of the representative heavy nucleus.
!-----------------------------------------------------------------------

ALLOCATE (xn_z(nz,nnc,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn_z        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uburn_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn_z     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc_rep_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep_z'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_rep_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep_z '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_rep_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep_z '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Nuclear statistical equilibrium
!-----------------------------------------------------------------------
!  nse_z(j,i) : a nuclear statistical equilibrium flag for angular zone j
!
!     nse_z(j,i) = 0 : material not in nuclear statistical equilibrium;
!      nuclear reaction network must be turned on to evolve the matter
!      composition.
!     nse_z(j,i) = 1 : material in nuclear statistical equilibrium;
!      nuclear reaction network turned off.
!-----------------------------------------------------------------------

ALLOCATE (nse_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse_z       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Cube grid
!-----------------------------------------------------------------------
!  idty_z(j,i) : the density regime (i.e., 1, 2, or 3) of radial zone j as
!   given by the above 
!  inequalities.
!     regime 1:             rho < rhoes(1)
!     regime 2:        rhoes(1) < rho < rhoes(2)
!     regime 3:             rhoes(2) < rho
!-----------------------------------------------------------------------

ALLOCATE (idty_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idty_z     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Moments of the neutrino distribution function
!-----------------------------------------------------------------------
!  psi0(j,k,n,i) : the zero moment of the occupation distribution for
!   neutrinos at the midpoint of angular zone j, of energy zone k, and
!   of type n.
!
!  psi1(j,k,n,i) : the first moment of the occupation distribution for
!   neutrinos at the outer boundary angular zone j, of energy zone k,
!   and of type n.
!
!  nu_str_cy           : y-component of zone-centered neutrino stress
!   (dynes g^{-1})
!
!  nu_str_ey           : y-component of zone-edge neutrino stress
!   (dynes g^{-1})
!-----------------------------------------------------------------------

ALLOCATE (psi0_z(nz,nez,nnu,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_z      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1_z(nz,nez,nnu,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1_z      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_str_cz(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_cy   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_str_ez(nz+1,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_ey   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_nu_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu_z      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                  \\\\\ INITIALIZE VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Time step determinations
!-----------------------------------------------------------------------

jdt_z                  = 0
dt_z                   = 1.d+20
dt_z_state             = 1.d+20

!-----------------------------------------------------------------------
!  Physical coordinate difference of z_ei(i+1) - z_ei(i)
!-----------------------------------------------------------------------

dzphys_c               = zero

!-----------------------------------------------------------------------
!  State variables
!-----------------------------------------------------------------------

rho_z                  = zero
t_z                    = zero
ye_z                   = zero
ei_z                   = zero
p_z                    = zero
gc_z                   = zero
ge_z                   = zero
u_z                    = zero
v_z                    = zero
w_z                    = zero
w_e                    = zero

!-----------------------------------------------------------------------
!  State variables - initial values
!-----------------------------------------------------------------------

rho_zi                 = zero
t_zi                   = zero
ye_zi                  = zero
ei_zi                  = zero
u_zi                   = zero
v_zi                   = zero
w_zi                   = zero

!-----------------------------------------------------------------------
!  Shock stabilizing variables
!-----------------------------------------------------------------------

flat_z                 = zero
flat_x_z               = zero

!-----------------------------------------------------------------------
!  Gravitational accelerations
!-----------------------------------------------------------------------

grav_z_cz              = zero
grav_z_ez              = zero

!-----------------------------------------------------------------------
!  Composition variables
!-----------------------------------------------------------------------

xn_z                   = zero

nse_z                  = 1
uburn_z                = zero
be_nuc_rep_z           = 492.3d0
a_nuc_rep_z            = 56.d0
z_nuc_rep_z            = 28.d0

!-----------------------------------------------------------------------
!  EOS cube regime indices
!-----------------------------------------------------------------------

idty_z                 = 0

!-----------------------------------------------------------------------
!  Moments of the neutrino distribution function
!-----------------------------------------------------------------------

psi0_z                 = zero
psi1_z                 = zero
nu_str_cz              = zero
nu_str_ez              = zero
e_nu_z                 = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_azimuthal_ray_arrays
                                                                                
