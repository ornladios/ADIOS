SUBROUTINE problem_setup( nx, ny, nz, nez, nnu, ij_ray_dim, ik_ray_dim, &
& nnc )
!-----------------------------------------------------------------------
!
!    File:         problem_setup
!    Module:       problem_setup
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  K. R. DeNisco, Dept. of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/10/04
!
!    Purpose:
!      To initialize the problem.
!
!
!    Subprograms called:
!  setup_initial_hydro    : initializes hydro varialbles from initial model
!  setup_restart_hydro    : initializes hydro varialbles from restart files
!  angular_ave            : computes the angular average of certain quantities
!  radhyd_to_poisson      : computes Legendre polynomilas in the Poisson solver
!  setup_initial_GR       : initializes GR varialbles from initial model
!  setup_restart_GR       : initializes GR varialbles from restart files
!  setup_initial_trans    : initializes transport varialbles from initial model
!  setup_restart_trans    : initializes transport varialbles from restart files
!  load_radial_ray_arrays : transfers computed variables to radial_ray_module
!
!    Input arguments:
!  nx                     : x_array extent
!  ny                     : y_array extent
!  nz                     : z_array extent
!  nez                    : neutrino energy array extent
!  nnu                    : neutrino flavor array extent
!  ij_ray_dim             : number of y-zones on a processor before swapping
!  ik_ray_dim             : number of z-zones on a processor before swapping
!  nnc                    : abundance array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  cycle_module, edit_module, parallel_module, radial_ray_module
! 
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE cycle_module, ONLY: nrst
USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid, ierr
USE radial_ray_module, ONLY : imin, imax, rho_c, t_c, ye_c, ye_ci, x_ei, &
& dx_ci, x_ci, u_c, u_e, v_c, w_c, agr_e, agr_c, agr_e_r, agr_c_r, xn_c, &
& be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, psi0_c, psi1_e, aesv_c, regrid, &
& rhobar, time, t_bounce, rho_regrid, nse_c
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-zone extent
INTEGER, INTENT(in)              :: ny            ! y-zone extent
INTEGER, INTENT(in)              :: nz            ! z-zone extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nnc           ! abundance array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: l_regrid      ! regrid flag

INTEGER                          :: ij_ray        ! j-index of a radial ray
INTEGER                          :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Hydro variables are being initialized on processor',i6,' for ray',2i4,' with nrst = 0')
  103 FORMAT (' Hydro variables have been being initialized on processor',i6,' for ray',2i4,' with nrst = 0')
  105 FORMAT (' Hydro variables have been initialized for all rays on processor',i6,' with nrst = 0')
  107 FORMAT (' Hydro variables are being initialized on processor',i6,' for ray',2i4,' with nrst /= 0')
  109 FORMAT (' Hydro variables have been initialized on processor',i6,' for ray',2i4,' with nrst /= 0')
  111 FORMAT (' Hydro variables have been initialized for all rays on processor',i6,' with nrst /= 0')
  113 FORMAT (' Grid is being initialized on processor',i6)
  115 FORMAT (' Grid has been initialized on processor',i6)
  121 FORMAT (' Angular averages of selected quantities are being computed on processor',i6)
  123 FORMAT (' Angular averages of selected quantities have been computed on processor',i6)
  125 FORMAT (' Legendre polynomials and gravity are being computed on processor',i6)
  127 FORMAT (' Legendre polynomials and gravity have been computed on processor',i6,' and on all processors')
  129 FORMAT (' Neutrino group energies are being computed in problem_setup on processor',i6)
  131 FORMAT (' Neutrino group energies have been computed in problem_setup on processor',i6)
  133 FORMAT (' GR variables are being initialized on processor',i6,' for ray',2i4,' with nrst = 0')
  135 FORMAT (' GR variables have been initialized on processor',i6,' for ray',2i4,' with nrst = 0')
  137 FORMAT (' GR variables have been initialized for all rays on processor',i6,' with nrst = 0')
  139 FORMAT (' GR variables are being initialized on processor',i6,' for ray',2i4,' with nrst /= 0')
  141 FORMAT (' GR variables have been initialized on processor',i6,' for ray',2i4,' with nrst /= 0')
  143 FORMAT (' GR variables have been initialized for all rays on processor',i6,' with nrst /= 0')
  145 FORMAT (' Transport variables are being initialized on processor',i6,' for ray',2i4,' with nrst = 0')
  147 FORMAT (' Transport variables have been initialized on processor',i6,' for ray',2i4,' with nrst = 0')
  149 FORMAT (' Transport variables have been initialized for all rays on processor',i6,' with nrst = 0')
  151 FORMAT (' Transport variables are being initialized on processor',i6,' for ray',2i4,' with nrst /= 0')
  153 FORMAT (' Transport variables have been initialized on processor',i6,' for ray',2i4,' with nrst /= 0')
  155 FORMAT (' Transport variables have been initialized for all rays on processor',i6,' with nrst /= 0')
  157 FORMAT (' Radial ray arrays are being loaded on processor',i6,' for ray',2i4)
  159 FORMAT (' Radial ray arrays have been loaded on processor',i6,' for ray',2i4)
  161 FORMAT (' Radial ray arrays have been loaded for all rays on processor',i6)
  167 FORMAT (' Regridding is being performed on processor',i6)
  169 FORMAT (' Regridding has been performed on processor',i6)

!-----------------------------------------------------------------------
!  Proceed with problem setup
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        (1) Compute rest masses and radii
!        (2) Fill variable arrays at time n-1 and time n+1 are filled
!        (3) Modify problem (if appropriate)
!        (4) Evaluate equation of state variables
!        (5) Compute the GR variables from the initial data
!        (6) The central and edge values of the neutrino energy bins are computed;
!        (7) The neutrino absorption, emission, scattering, pair production, and bremsstrahlung
!         pair production function tables are initialized
!
!        Calls
!          pblmst1
!          esrgnz_x
!          eqstz_x
!          gammaz_x
!          setup_rel
!          gamgr_nu_cal
!          agr_nu_cal
!          e_zone
!          enu_cal
!          pre_trans
!          pblmst2
!          gennur
!          abemset
!          scataset
!          scateset
!          scatiset
!          pairset
!          bremset
!          scatnset
!          scatnnset
!          scatnAset
!          abemrate
!          sctirate
!          sctarate
!          scterate
!          pairrate
!          bremrate
!          sctnrate
!          sctnnrate
!          nu_number
!          mfp_cal
!          nu_sphere
!          diffc
!          nu_stress_x
!          eddington
!          nu_U
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!       \\\\\ INITIALIZE HYDRO AND THERMODYNAMIC VARIABLES /////
!
!-----------------------------------------------------------------------

IF ( nrst == 0 ) THEN

!-----------------------------------------------------------------------
!                ||||| BEGIN SUM OVER RADIAL RAYS |||||
!-----------------------------------------------------------------------

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Initialize hydro and thermodynamic variables
!-----------------------------------------------------------------------

      WRITE (nlog,101) myid, ij_ray, ik_ray
      CALL setup_initial_hydro( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
&      ik_ray_dim, nx, ny, nz, nnc, rho_c, t_c, ye_c, ye_ci, x_ei, dx_ci, &
&      u_c, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, aesv_c )
      WRITE (nlog,103) myid, ij_ray, ik_ray
      IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,105) myid

!-----------------------------------------------------------------------
!                 ||||| END SUM OVER RADIAL RAYS |||||
!-----------------------------------------------------------------------

    END DO ! ij_ray
  END DO ! ik_ray

!-----------------------------------------------------------------------
!  Compute the angular averages of certian quantities
!-----------------------------------------------------------------------

  WRITE (nlog,121) myid
  CALL angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )
  WRITE (nlog,123) myid

!-----------------------------------------------------------------------
!  Regrid model to constant drho/rho in the NSE material
!-----------------------------------------------------------------------

  WRITE (nlog,167) myid
  CALL regridder( imin, imax, nx, ij_ray_dim, ny, ik_ray_dim, nz, nez, &
&  nnu, nnc, nse_c, x_ei, dx_ci, x_ci, rhobar, rho_c, t_c, ye_c, u_c, &
&  u_e, v_c, w_c, psi0_c, regrid, rho_regrid, time, t_bounce, l_regrid )
  WRITE (nlog,169) myid

ELSE ! nrst /= 0

!-----------------------------------------------------------------------
!                ||||| BEGIN SUM OVER RADIAL RAYS |||||
!-----------------------------------------------------------------------

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Reinitialize hydro and thermodynamic variables
!-----------------------------------------------------------------------

      WRITE (nlog,107) myid, ij_ray, ik_ray
      CALL setup_restart_hydro( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
&      ik_ray_dim, nx, ny, nz, nnc, rho_c, t_c, ye_c, ye_ci, x_ei, dx_ci, &
&      u_c, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, aesv_c )
      WRITE (nlog,109) myid, ij_ray, ik_ray
      IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,111) myid

!-----------------------------------------------------------------------
!                 ||||| END SUM OVER RADIAL RAYS |||||
!-----------------------------------------------------------------------

    END DO ! ij_ray
  END DO ! ik_ray

!-----------------------------------------------------------------------
!  Compute the angular averages of certian quantities
!-----------------------------------------------------------------------

  WRITE (nlog,121) myid
  CALL angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )
  WRITE (nlog,123) myid

!-----------------------------------------------------------------------
!  Initialize hydro and thermodynamic variables complete
!-----------------------------------------------------------------------

END IF ! nrst == 0

!-----------------------------------------------------------------------
!  Set up the grid
!-----------------------------------------------------------------------

WRITE (nlog,113) myid
CALL initialize_grid
WRITE (nlog,115) myid

!-----------------------------------------------------------------------
!  Evaluate Legendre polynomials and compute gravity.
!-----------------------------------------------------------------------

WRITE (nlog,125) myid
CALL radhyd_to_poisson( 1, nx, ij_ray_dim, ik_ray_dim, ny, nz )
WRITE (nlog,127) myid

!-----------------------------------------------------------------------
!  Neutrino group energies at infinity
!-----------------------------------------------------------------------

WRITE (nlog,129) myid
CALL e_zone
WRITE (nlog,131) myid

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE GRAVITY AND TRANSPORT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                ||||| BEGIN SUM OVER RADIAL RAYS |||||
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Initialize relativistic variables
!-----------------------------------------------------------------------

    IF ( nrst == 0 ) THEN

      WRITE (nlog,133) myid, ij_ray, ij_ray
      CALL setup_initial_GR( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
&      ik_ray_dim, nx, ny, nz, nez, nnu, rho_c, x_ei, dx_ci, agr_e, agr_c )
      agr_e_r         = agr_e
      agr_c_r         = agr_c
      WRITE (nlog,135) myid, ij_ray, ij_ray
      IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) &
&      WRITE (nlog,137) myid

    ELSE ! nrst /= 0

      WRITE (nlog,139) myid, ij_ray, ij_ray
      CALL setup_restart_GR( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
&      ik_ray_dim, nx, ny, nz, nez, nnu, rho_c, x_ei, dx_ci, agr_e, agr_c )
      agr_e_r         = agr_e
      agr_c_r         = agr_c
      WRITE (nlog,141) myid, ij_ray, ij_ray
      IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) &
&      WRITE (nlog,143) myid

    END IF ! nrst == 0

!-----------------------------------------------------------------------
!  Initialize neutrino transport variables
!-----------------------------------------------------------------------

    IF ( nrst == 0 ) THEN

      WRITE (nlog,145) myid, ij_ray, ij_ray
      CALL setup_initial_trans( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
&      ik_ray_dim, nx, nez, nnu, rho_c, t_c, ye_c, x_ei, dx_ci, u_c, &
&      agr_e, agr_c, psi0_c, psi1_e, rhobar )
      WRITE (nlog,147) myid, ij_ray, ij_ray
      IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) &
&      WRITE (nlog,149) myid

    ELSE ! nrst /= 0

      WRITE (nlog,151) myid, ij_ray, ij_ray
      CALL setup_restart_trans( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
&      ik_ray_dim, nx, nez, nnu, rho_c, t_c, ye_c, x_ei, dx_ci, u_c, &
&      agr_e, agr_c, psi0_c, psi1_e, rhobar )
      WRITE (nlog,153) myid, ij_ray, ij_ray
      IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) &
&      WRITE (nlog,155) myid

    END IF ! nrst == 0

!-----------------------------------------------------------------------
!  Load thermodynamic and neutrino variables into radhyd_ray arrays
!-----------------------------------------------------------------------

    WRITE (nlog,157) myid, ij_ray, ij_ray
    CALL load_radial_ray_arrays( nx, nez, nnu, nnc, ij_ray, ik_ray )
    WRITE (nlog,159) myid, ij_ray, ij_ray
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) &
&    WRITE (nlog,161) myid

!-----------------------------------------------------------------------
!                 ||||| END SUM OVER RADIAL RAYS |||||
!-----------------------------------------------------------------------

  END DO ! ij_ray
END DO ! ik_ray

RETURN
END SUBROUTINE problem_setup

