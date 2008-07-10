SUBROUTINE dimension_angular_ray_arrays( ny, j_ray_dim, ik_ray_dim, nez, &
& nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_angular_ray_arrays
!    Module:       dimension_angular_ray_arrays
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
!  ny         : y (angular) array extent
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  j_ray_dim  : number of radial zones on a processor after swapping with y
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

USE numerical_module, ONLY : zero, one

USE angular_ray_module
USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: ny            ! y-array extent
INTEGER, INTENT(in)               :: j_ray_dim     ! number of y-zones on a processor after swapping with y
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
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

  101 FORMAT (' Radhyd angular arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_angular_ray_arrays')

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

ALLOCATE (jdt_y(3,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'jdt_y       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dt_y(3,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dt_y        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ PHYSICAL ANGULAR COORDINATES /////
!
!-----------------------------------------------------------------------

ALLOCATE (dyphys_c(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dyphys_c    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------

ALLOCATE (rho_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_y       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_y         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_y        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ei_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_y        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (p_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_y         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gc_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gc_y        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ge_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ge_y        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (u_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_y         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_y         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_y         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (v_e(ny+1,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_e         '; WRITE  (*,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - INITIAL VALUES /////
!
!-----------------------------------------------------------------------

ALLOCATE (rho_yi(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_yi      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_yi(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_yi        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_yi(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_yi       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ei_yi(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_yi       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (u_yi(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_yi        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v_yi(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_yi        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w_yi(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_yi        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                     \\\\\ GR VARIABLES /////
!
!-----------------------------------------------------------------------
!  agr_y : lapse function
!-----------------------------------------------------------------------

ALLOCATE (agr_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_y       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!               \\\\\ SHOCK STABILIZING VARIABLES /////
!
!-----------------------------------------------------------------------
!  flat_y   : variables aligned along the y-axis indicating the presence
!   of an angular shock
!
!  flat_x_y : variables aligned along the y-axis indicating the presence
!   of a radial shock
!-----------------------------------------------------------------------

ALLOCATE (flat_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_x      '; WRITE  (*,1001) var_name; END IF
ALLOCATE (flat_x_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_x_y    '; WRITE  (*,1001) var_name; END IF

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

ALLOCATE (grav_y_cy(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_y_cy   '; WRITE  (*,1001) var_name; END IF
ALLOCATE (grav_y_ey(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_y_ey   '; WRITE  (*,1001) var_name; END IF

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

ALLOCATE (xn_y(ny,nnc,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn_y        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uburn_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn_y     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc_rep_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep_y'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_rep_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep_y '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_rep_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep_y '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Nuclear statistical equilibrium
!-----------------------------------------------------------------------
!  nse_y(j,i) : a nuclear statistical equilibrium flag for angular zone j
!
!     nse_y(j,i) = 0 : material not in nuclear statistical equilibrium;
!      nuclear reaction network must be turned on to evolve the matter
!      composition.
!     nse_y(j,i) = 1 : material in nuclear statistical equilibrium;
!      nuclear reaction network turned off.
!-----------------------------------------------------------------------

ALLOCATE (nse_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse_y       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Cube grid
!-----------------------------------------------------------------------
!  idty_y(j,i) : the density regime (i.e., 1, 2, or 3) of radial zone j as
!   given by the above 
!  inequalities.
!     regime 1:             rho < rhoes(1)
!     regime 2:        rhoes(1) < rho < rhoes(2)
!     regime 3:             rhoes(2) < rho
!-----------------------------------------------------------------------

ALLOCATE (idty_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idty_y     '; WRITE (nlog,1001) var_name; END IF

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

ALLOCATE (psi0_y(ny,nez,nnu,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_y      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1_y(ny,nez,nnu,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1_y      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_str_cy(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_cy   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_str_ey(ny+1,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_ey   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_nu_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu_y      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                  \\\\\ INITIALIZE VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Time step determinations
!-----------------------------------------------------------------------

jdt_y                  = 0
dt_y                   = 1.d+20
dt_y_state             = 1.d+20

!-----------------------------------------------------------------------
!  Physical angular coordinate differences
!-----------------------------------------------------------------------

dyphys_c               = zero

!-----------------------------------------------------------------------
!  State variables
!-----------------------------------------------------------------------

rho_y                  = zero
t_y                    = zero
ye_y                   = zero
ei_y                   = zero
p_y                    = zero
gc_y                   = zero
ge_y                   = zero
u_y                    = zero
v_y                    = zero
w_y                    = zero
v_e                    = zero

!-----------------------------------------------------------------------
!  State variables - initial values
!-----------------------------------------------------------------------

rho_yi                 = zero
t_yi                   = zero
ye_yi                  = zero
ei_yi                  = zero
u_yi                   = zero
v_yi                   = zero
w_yi                   = zero

!-----------------------------------------------------------------------
!  Lapse function
!-----------------------------------------------------------------------

agr_y                  = one

!-----------------------------------------------------------------------
!  Shock stabilizing variables
!-----------------------------------------------------------------------

flat_y                 = zero
flat_x_y               = zero

!-----------------------------------------------------------------------
!  Gravitational accelerations
!-----------------------------------------------------------------------

grav_y_cy              = zero
grav_y_ey              = zero

!-----------------------------------------------------------------------
!  Composition variables
!-----------------------------------------------------------------------

xn_y                   = zero

nse_y                  = 1
uburn_y                = zero
be_nuc_rep_y           = 492.3d0
a_nuc_rep_y            = 56.d0
z_nuc_rep_y            = 28.d0

!-----------------------------------------------------------------------
!  EOS cube regime indices
!-----------------------------------------------------------------------

idty_y                 = 0

!-----------------------------------------------------------------------
!  Moments of the neutrino distribution function
!-----------------------------------------------------------------------

psi0_y                 = zero
psi1_y                 = zero
nu_str_cy              = zero
nu_str_ey              = zero
e_nu_y                 = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_angular_ray_arrays
                                                                                
