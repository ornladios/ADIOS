SUBROUTINE dimension_radial_ray_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_radial_ray_arrays
!    Module:       dimension_radial_ray_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/04
!
!    Purpose:
!      To allocate the dimensions of the radhyd arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array extent
!  ny         : y (angular) array extent
!  nz         : z (azimuthal) array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  j_ray_dim  : number of radial zones on a processor after swapping with y
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  nnc        : composition array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero, one

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE radial_ray_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: ny            ! y-array extent
INTEGER, INTENT(in)               :: nz            ! z-array extent
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)               :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)               :: k_ray_dim     ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)               :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=14)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Radhyd_ray arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_radial_ray_arrays')

!-----------------------------------------------------------------------
!
!                \\\\\ ALLOCATE RADIALRAY_MODULE /////
!                /////          ARRAYS           \\\\\
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!           \\\\\ TIME STEPS AND TIME STEP CONSTRAINTS /////
!
!-----------------------------------------------------------------------
!  jdt          : radial zone for which dt(j,i_ray) us a minimum
!
!  j_radial_dt  : radial ray responsible for each dt_process
!
!  j_angular_dt : angular ray responsible for each dt_process
!
!  dt_process   : minimum time step for process i
!
!  dt           : minimum time step for process i along radial ray i_ray
!
!  dtime_trans  : minimum time step source and transport from radial ray
!   i_ray
!-----------------------------------------------------------------------

ALLOCATE (jdt(50,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'jdt           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (j_radial_dt(50), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'j_radial_dt   '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (j_angular_dt(50), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'j_angular_dt  '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (j_azimuthal_dt(50), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'j_azimuthal_dt'; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dt_process(50), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dt_process    '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dt(50,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dt            '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dtime_trans(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtime_trnas   '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                \\\\\ COMPOSITION BOUNDARY VALUES /////
!
!-----------------------------------------------------------------------
!  comp_bcl : right boundary composition
!
!  comp_bcr : right boundary composition
!-----------------------------------------------------------------------

ALLOCATE (comp_bcl(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp_bcl      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (comp_bcr(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp_bcr      '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ SOLID ANGLES SUBTENDED BY RAYS /////
!
!-----------------------------------------------------------------------
!  d_omega : solid angles subtended by radial rays
!
!  cos_theta : cos( theta )
!
!  sin_theta : sin( theta )
!
!  cos_phi   : cos ( phi )
!
!  sin_phi   : sin( phi )
!-----------------------------------------------------------------------

ALLOCATE (d_omega(ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_omega       '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (cos_theta(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cos_theta     '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (sin_theta(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sin_theta     '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (cos_phi(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cos_phi       '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (sin_phi(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sin_phi       '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!         \\\\\ COORDINATES - VALUES AT CYCLE BEGINNING /////
!
!-----------------------------------------------------------------------
!  x_ei  : x grid zone face locations at beginning of cycle
!
!  y_ei  : y grid zone face locations at beginning of cycle
!
!  z_ei  : z grid zone face locations at beginning of cycle
!
!  dx_ci : x_ei(i+1) - x_ei(i)
!
!  dy_ci : y_ei(i+1) - y_ei(i)
!
!  dz_ci : z_ei(i+1) - z_ei(i)
!
!  x_ci  : x grid zone midpoint locations at beginning of cycle
!
!  y_ci  : y grid zone midpoint locations at beginning of cycle
!
!  z_ci  : z grid zone midpoint locations at beginning of cycle
!-----------------------------------------------------------------------

ALLOCATE (x_ei(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'x_ei          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (y_ei(ny+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y_ei          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_ei(nz+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_ei          '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (dx_ci(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dx_ci         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dy_ci(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dy_ci         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dz_ci(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dz_ci         '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (x_ci(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'x_ci          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (y_ci(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y_ci          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_ci(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_ci          '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!       \\\\\ COORDINATES - VALUES AFTER LAGRANGIAN UPDATE /////
!
!-----------------------------------------------------------------------
!  x_el  : x grid zone face locations after Lagrangian update
!
!  y_el  : y grid zone face locations after Lagrangian update
!
!  z_el  : z grid zone face locations after Lagrangian update
!
!  dx_cl : x_el(i+1) - x_el(i)
!
!  dy_cl : y_el(i+1) - y_el(i)
!
!  dz_cl : z_el(i+1) - z_el(i)
!
!  x_cl  : x grid zone midpoint locations after Lagrangian update
!
!  y_cl  : y grid zone midpoint locations after Lagrangian update
!
!  z_cl  : z grid zone midpoint locations after Lagrangian update
!-----------------------------------------------------------------------

ALLOCATE (x_el(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'x_el          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (y_el(ny+1,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y_el          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_el(nz+1,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_el          '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (dx_cl(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dx_cl         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dy_cl(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dy_cl         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dz_cl(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dz_cl         '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (x_cl(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'x_cl          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (y_cl(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y_cl          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_cl(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_cl          '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ COORDINATES - VALUES AT CYCLE END /////
!
!-----------------------------------------------------------------------
!  x_ef  : x grid zone face locations at cycle end
!
!  y_ef  : y grid zone face locations at cycle end
!
!  z_ef  : z grid zone face locations at cycle end
!
!  dx_cf : x_ef(i+1) - x_ef(i)
!
!  dy_cf : y_ef(i+1) - y_ef(i)
!
!  dz_cf : z_ef(i+1) - z_ef(i)
!
!  x_cf  : x grid zone midpoint locations at cycle end
!
!  y_cf  : y grid zone midpoint locations at cycle end
!
!  z_cf  : z grid zone midpoint locations at cycle end
!-----------------------------------------------------------------------

ALLOCATE (x_ef(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'x_ef          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (y_ef(ny+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y_ef          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_ef(nz+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_ef          '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (dx_cf(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dx_cf         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dy_cf(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dy_cf         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dz_cf(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dz_cf         '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (x_cf(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'x_cf          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (y_cf(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y_cf          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_cf(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_cf          '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  rho_c    : density: zone average (g cm^{-3})
!
!  t_c      : temperature: zone average (K)
!
!  ye_c     : electron fraction: zone average
!
!  ei_c     : internal energy: zone average (ergs g^{-1})
!
!  e_v_c    : total energy minus kinetic energy (ergs g^{-1})
!
!  p_c      : pressure: zone average (ergs cm^{-3})
!
!  gc_c     : 1st adiabatic index
!
!  ge_c     : 1 + p/(ei*rho)
!
!  u_c      : zone centered average velocity x direction (cm s^{-1})
!
!  u_e      : x grid edge velocity for moving grid option
!
!  v_c      : zone centered average velocity y direction (cm s^{-1})
!
!  w_c      : zone centered average velocity z direction (cm s^{-1})
!-----------------------------------------------------------------------

ALLOCATE (rho_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_c         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (t_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_c           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (ye_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_c          '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (ei_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_c          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (e_v_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_v_c         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (p_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_c           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (gc_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gc_c          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (ge_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ge_c          '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (u_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_c           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (u_e(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_e           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (v_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_c           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (w_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_c           '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!         \\\\\ STATE VARIABLES - AFTER LAGRANGIAN UPDATE/////
!
!-----------------------------------------------------------------------
!  u_l      : x grid edge velocity for moving grid option
!
!  v_l      : zone centered average velocity y direction (cm s^{-1})
!
!  w_l      : zone centered average velocity z direction (cm s^{-1})
!
!  rho_l    : density: zone average (g cm^{-3})
!-----------------------------------------------------------------------

ALLOCATE (u_l(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_l           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (v_l(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_l           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (w_l(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_l           '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (rho_l(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_l         '; WRITE  (nlog,1001) var_name; END IF


!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - ANGULAR AVERAGES /////
!
!-----------------------------------------------------------------------
!  rhobar(j)           : angular averaged density (g cm^{-3})
!
!  p_bar (j)           : angular averaged matter pressure (ergs cm^{-3})
!
!  e_bar (j)           : angular averaged matter internal energy (ergs
!   cm^{-3})
!
!  vx_bar(j)           : angular averaged radial velocity (cm s^{-1})
!
!  v2_bar(j)           : angular averaged velocity square (cm^{2} s^{-2})
!
!  e_nu_c_bar(j)       : angular averaged neutrino energy density
!   (ergs cm^{-3}).
!
!  f_nu_e_bar(j)       : angular averaged neutrino energy flux
!   (ergs cm^{-2} s^{-1}).
!-----------------------------------------------------------------------

ALLOCATE (rhobar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhobar        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (p_bar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_bar         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (e_bar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_bar         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (vx_bar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vx_bar        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (v2_bar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v2_bar        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (e_nu_c_bar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu_c_bar    '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (f_nu_e_bar(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f_nu_e_bar    '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - INITIAL VALUES /////
!
!-----------------------------------------------------------------------
!  rho_ci : density: zone average (g cm^{-3}) at beginning of cycle
!
!  t_ci   : temperature: zone average (K) at beginning of cycle
!
!  ye_ci  : electron fraction: zone average at beginning of cycle
!
!  ei_ci  : internal energy: zone average (ergs g^{-1}) at beginning of
!   cycle
!
!  u_ci   : zone centered average velocity x direction (cm s^{-1}) at
!   beginning of cycle
!
!  v_ci   : zone centered average velocity y direction (cm s^{-1}) at
!   beginning of cycle
!
!  w_ci   : zone centered average velocity z direction (cm s^{-1}) at
!   beginning of cycle
!-----------------------------------------------------------------------

ALLOCATE (rho_ci(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_ci        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (t_ci(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_ci          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (ye_ci(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_ci         '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (ei_ci(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_ci         '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (u_ci(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_ci          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (v_ci(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_ci          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (w_ci(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_ci          '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                     \\\\\ GR VARIABLES /////
!
!-----------------------------------------------------------------------
!  agr_e   : zone-edged value of the lapse function after gravitational
!   update
!
!  agr_c   : zone-centered value of the lapse function after
!   gravitational update
!
!  agr_e_r : zone-edged value of the lapse function priot to
!   gravitational update
!
!  agr_c_r : zone-centered value of the lapse function prior to
!   gravitational update
!-----------------------------------------------------------------------

ALLOCATE (agr_e(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_e         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (agr_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_c         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (agr_e_r(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_e_r       '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (agr_c_r(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr_c_r       '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ SHOCK STABILIZING VARIABLES /////
!
!-----------------------------------------------------------------------
!  flat_x   : variables indicating the presence of a radial shock
!
!  flat_y_x : variables indicating the presence of an angular shock
!   aligned along the shock
!-----------------------------------------------------------------------

ALLOCATE (flat_x(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_x        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (flat_y_x(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_y_x      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (flat_z_x(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flat_z_x      '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!         \\\\\ GRAVITATION ACCELERATIONS AND POTENTIALS /////
!
!-----------------------------------------------------------------------
!  grav_x_c   : zone-centered x-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_y_c   : zone-centered y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_z_c   : zone-centered z-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_pot_c : zone-centered gravitational potential (energy required
!   to move a unit mass from x to infinity (erg g^{-1})
!
!  grav_pot_c_i : zone-centered gravitational potential (energy required
!   to move a unit mass from x to infinity (erg g^{-1}) saved for
!   subtraction from grav_pot_c computed at the end of the remap to
!   calculate dPhi/dt
!
!  gtot_pot_c : zone-centered gravitational potential (energy required
!   to move a unit mass from x to infinity with no mass overlying it
!   (erg g^{-1})
!
!  grav_x_e   : zone-edged x-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_y_e   : zone-edged y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_z_e   : zone-edged y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_pot_e : zone-edged gravitational potential (energy required
!   to move a unit mass from x to infinity (erg g^{-1})
!
!  gtot_pot_e : zone-edged gravitational potential (energy required
!   to move a unit mass from x to infinity with no mass overlying it
!   (erg g^{-1})
!-----------------------------------------------------------------------

ALLOCATE (grav_x_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_x_c      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_y_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_y_c      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_z_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_z_c      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_pot_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_pot_c    '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_pot_c_i(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_pot_c_i  '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (gtot_pot_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gtot_pot_c    '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_x_e(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_x_e      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_y_e(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_y_e      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_z_e(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_z_e      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (grav_pot_e(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'grav_pot_e    '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (gtot_pot_e(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gtot_pot_e    '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!          \\\\\ COMPOSITION VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  xn_c(i,ic,j,k)        : mass fraction of the ith nucleus.
!
!  uburn_c(i,j,k)        : cumulative energy generated in zone j by
!   nuclear reactions (ergs/gm).
!
!  be_nuc_rep_c(i,j,k)   : binding energy of the representative heavy
!   nucleus (MeV).
!
!  a_nuc_rep_c(i,j,k)    : mass number of the representative heavy nucleus.
!
!  z_nuc_rep_c(i,j,k)    : charge number of the representative heavy nucleus.
!
!  a_nuc_c(n)            : mass number of the nth nuclear species.
!
!  z_nuc_c(n)            : charge number of the nth nuclear species.
!
!  be_nuc_c(n)           : binding energy of the nth nuclear species (MeV).
!
!  e_bind_zn_c(j,k)      : total nuclear binding energy initially in 
!                           each mass shell (ergs)
!
!  fluxbe_c(j,k)         : total nuclear binding energy transferred
!                           during remap (ergs)
!
!  eb_c(j,k)             : mean nuclear binding energy (ergs g^{-1})
!-----------------------------------------------------------------------

ALLOCATE (xn_c(nx,nnc,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn_c          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (uburn_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn_c       '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_rep_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep_c   '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_rep_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep_c   '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (be_nuc_rep_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep_c  '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_c(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_c       '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_c(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_c       '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (be_nuc_c(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_c      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (e_bind_zn_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_bind_zn_c   '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (fluxbe_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxbe_c      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (eb_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eb_c          '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                      \\\\\ NSE FLAG /////
!
!-----------------------------------------------------------------------
!  nse(i,j,k) : a nuclear statistical equilibrium flag for zone i,j,k.
!
!     nse(i,j,k) = 0 : material not in nuclear statistical equilibrium;
!      nuclear reaction network must be turned on to evolve the matter
!      composition.
!     nse(i,j,k) = 1 : material in nuclear statistical equilibrium;
!      nuclear reaction network turned off.
!-----------------------------------------------------------------------

ALLOCATE (nse_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse_c         '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                  \\\\\ EOS AND OPACITY GRID /////
!
!-----------------------------------------------------------------------
!  idty(j) : the density regime (i.e., 1, 2, or 3) of radial zone j as
!   given by the following inequalities.
!          regime 1:             rho < rhoes(1)
!          regime 2:        rhoes(1) < rho < rhoes(2)
!          regime 3:             rhoes(2) < rho
!-----------------------------------------------------------------------

ALLOCATE (idty(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idty          '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!        \\\\\ INTERPOLATED EQUATION OF STATE VARIABLES /////
!
!-----------------------------------------------------------------------
!  aesv(j,i,i_ray)  : equation of state dependent variable i at radial
!   zone j.
!
!     i = 1   : pressure
!     i = 2   : energy
!     i = 3   : entropy
!     i = 4   : neutron chemical potential
!     i = 5   : proton chemical potential
!     i = 6   : electron chemical potential
!     i = 7   : free neutron mass fraction
!     i = 8   : free proton mass fraction
!     i = 9   : heavy nucleus mass fraction
!     i = 10  : heavy nucleus mass number
!     i = 11  : heavy nucleus charge number
!     i = 12  : gamma1
!-----------------------------------------------------------------------

ALLOCATE (aesv_c(nx,12,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesv_c        '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                   \\\\\ NEUTRINO ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (psi0_c(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_c        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (psi1_e(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1_e        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nu_str_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_c      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nu_str_e(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_str_e      '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (rhs1_c(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhs1_c        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dc_e(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_e          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (e_nu_c(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu_c        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (f_nu_e(nx+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f_nu_e        '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                   \\\\\ NEUTRINO ENERGIES /////
!
!-----------------------------------------------------------------------

ALLOCATE (unu_c(nx,nez,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unu_c         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (unub_c(nx,nez+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unub_c        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dunu_c(nx,nez,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunu_c        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (unue_e(nx,nez,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unue_e        '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (unube_e(nx,nez+1,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unube_e       '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (dunue_e(nx,nez,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunue_e       '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                     \\\\\ EDIT COUNTERS /////
!
!-----------------------------------------------------------------------

ALLOCATE (nedc(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedc          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nede(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nede          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedmi(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedmi         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedma(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedma         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedh(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedh          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedps(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedps         '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedu(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedu          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedy(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedy          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedsc(20,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedsc         '; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (nedn(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedn          '; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (nedng(100,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedng         '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                  \\\\\ INITIALIZE VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Current cycle number and maximum cycle number
!-----------------------------------------------------------------------

ncycle                 = 0
ncymax                 = 9000000

!-----------------------------------------------------------------------
!  Time amnd time steps
!-----------------------------------------------------------------------

time                   = zero
t_start                = zero
t_bounce               = zero
t_stop                 = zero
tb_stop                = zero

dtnph                  = 1.d+20
dtnmh                  = 1.d+20

jdt                    = 0
j_radial_dt            = 0
j_angular_dt           = 0
j_azimuthal_dt         = 0
dt                     = 1.d+20

!-----------------------------------------------------------------------
!  Geometry and boundary conditions
!-----------------------------------------------------------------------

ndim                   = -1
ngeomx                 = -1
ngeomy                 = -1
ngeomz                 = -1
nleftx                 = -1
nrightx                = -1
nlefty                 = -1
nrighty                = -1
nleftz                 = -1
nrightz                = -1

!-----------------------------------------------------------------------
!  Minimum and maximum logical dimensions
!-----------------------------------------------------------------------

imin                   = 1
imax                   = 1
jmin                   = 1
jmax                   = 1
kmin                   = 1
kmax                   = 1

!-----------------------------------------------------------------------
!  Neutron star parameters
!-----------------------------------------------------------------------

vel_ns                 = zero
mass_ns                = zero

!-----------------------------------------------------------------------
!  Problem dimensions
!-----------------------------------------------------------------------

xmin                   = zero
xmax                   = zero
ymin                   = zero
ymax                   = 1.d0
zmin                   = zero
zmax                   = 2.d0

xmin_i                 = zero
xmax_i                 = zero
ymin_i                 = zero
ymax_i                 = 1.d0
zmin_i                 = zero
zmax_i                 = 2.d0

!-----------------------------------------------------------------------
!  Lagrangina or Eulerian
!-----------------------------------------------------------------------

lagr                   = 'no'
t_bounce_lagr_chg      = 1.d-03

!-----------------------------------------------------------------------
!  Regridding option
!-----------------------------------------------------------------------

rezn                   = 'no'
t_bounce_mgrd_chg      = 1.d-03

!-----------------------------------------------------------------------
!  Angles - initial values
!-----------------------------------------------------------------------

omega                  = zero
d_omega                = zero
cos_theta              = zero
sin_theta              = zero
cos_phi                = zero
sin_phi                = zero

!-----------------------------------------------------------------------
!  Coordinate - initial values
!-----------------------------------------------------------------------

x_ei                   = zero
dx_ci                  = zero
x_ci                   = zero

y_ei                   = zero
dy_ci                  = zero
y_ci                   = zero

z_ei                   = zero
dz_ci                  = zero
z_ci                   = zero

!-----------------------------------------------------------------------
!  Coordinates after Lagrangian update
!-----------------------------------------------------------------------

x_el                   = zero
dx_cl                  = zero
x_cl                   = zero

y_el                   = zero
dy_cl                  = zero
y_cl                   = zero

z_el                   = zero
dz_cl                  = zero
z_cl                   = zero

!-----------------------------------------------------------------------
!  Coordinate - final values
!-----------------------------------------------------------------------

x_ef                   = zero
dx_cf                  = zero
x_cf                   = zero

y_ef                   = zero
dy_cf                  = zero
y_cf                   = zero

z_ef                   = zero
dz_cf                  = zero
z_cf                   = zero

u_e                    = zero

!-----------------------------------------------------------------------
!  State variables - current values
!-----------------------------------------------------------------------

rho_c                  = zero
t_c                    = zero
ye_c                   = zero
ei_c                   = zero
e_v_c                  = zero
p_c                    = zero
gc_c                   = zero
ge_c                   = zero
u_c                    = zero
v_c                    = zero
w_c                    = zero
nu_str_c               = zero
nu_str_e               = zero
rhs1_c                 = zero
dc_e                   = zero
e_nu_c                 = zero
f_nu_e                 = zero

!-----------------------------------------------------------------------
!  State variables - after Lagrangian update
!-----------------------------------------------------------------------

u_l                    = zero
v_l                    = zero
w_l                    = zero
rho_l                  = zero

!-----------------------------------------------------------------------
!  State variables - angular averages
!-----------------------------------------------------------------------

rhobar                 = zero
p_bar                  = zero
e_bar                  = zero
vx_bar                 = zero
v2_bar                 = zero
e_nu_c_bar             = zero
f_nu_e_bar             = zero

!-----------------------------------------------------------------------
!  State variables - initial values
!-----------------------------------------------------------------------

rho_ci                 = zero
t_ci                   = zero
ye_ci                  = zero
ei_ci                  = zero
u_ci                   = zero
v_ci                   = zero
w_ci                   = zero

!-----------------------------------------------------------------------
!  GR variables
!-----------------------------------------------------------------------

agr_e                  = one
agr_c                  = one
agr_e_r                = one
agr_c_r                = one

!-----------------------------------------------------------------------
!  Shock stabilizing variables
!-----------------------------------------------------------------------

flat_x                 = zero
flat_y_x               = zero
flat_z_x               = zero

!-----------------------------------------------------------------------
!  Gravitation acceleration and potentials
!-----------------------------------------------------------------------

grav_x_c               = zero
grav_y_c               = zero
grav_z_c               = zero
grav_pot_c             = zero
grav_pot_c_i           = zero
gtot_pot_c             = zero
grav_x_e               = zero
grav_y_e               = zero
grav_z_e               = zero
grav_pot_e             = zero
gtot_pot_e             = zero

!-----------------------------------------------------------------------
!  Composition variables
!-----------------------------------------------------------------------

xn_c                   = zero

nse_c                  = 1
uburn_c                = zero
a_nuc_rep_c            = 56.d0
z_nuc_rep_c            = 28.d0
be_nuc_rep_c           = 492.3d0
a_nuc_c                = zero
z_nuc_c                = zero
be_nuc_c               = zero
e_bind_zn_c            = zero
fluxbe_c               = zero
eb_c                   = zero

!-----------------------------------------------------------------------
!  EOS cube regime indices
!-----------------------------------------------------------------------

idty                   = 0

!-----------------------------------------------------------------------
!  Interpolated EOS variables
!-----------------------------------------------------------------------

aesv_c                 = zero

!-----------------------------------------------------------------------
!  Moments of the neutrino distribution function
!-----------------------------------------------------------------------

psi0_c                 = zero
psi1_e                 = zero

!-----------------------------------------------------------------------
!  Neutrino energies
!-----------------------------------------------------------------------

unu_c                  = zero
unub_c                 = zero
dunu_c                 = zero
unue_e                 = zero
unube_e                = zero
dunue_e                = zero

!-----------------------------------------------------------------------
!  Edit counters
!-----------------------------------------------------------------------

nedc                   = 0
nede                   = 0
nedmi                  = 0
nedma                  = 0
nedh                   = 0
nedps                  = 0
nedu                   = 0
nedy                   = 0
nedsc                  = 0

nedn                   = 0
nedng                  = 0

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_radial_ray_arrays
