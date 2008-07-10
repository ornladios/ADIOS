SUBROUTINE radhyd_x_W( nx, ny, nz, nez, nnu, nnc, ndim, ij_ray_min, &
& ij_ray_max, ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, i_edit, &
& time, dtnph, nu_equil )
!-----------------------------------------------------------------------
!
!    File:         radhyd_x_W
!    Module:       radhyd_x_W
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/23/07
!
!    Purpose:
!      To direct the computation of the x-evolution.
!
!    Subprograms called:
!  parameter_reset            : change rezn or m_grid parameters according to criteria
!  store_int_grid             : store the initial grid
!  store_int_radial_var       : store the initial radial ray variables
!  radhyd_to_nu_energy_flux   : transfer variables and computes the neutrino energy density and flux
!  angular_ave                : calculates the angular average of quantities
!  radhyd_to_grav_i           : transfer variables and compute the gravitational potential and accelerations at beginning of cycle
!  radhyd_to_nu_agr_e_advct_x : update neutrino distributions from lapse update
!  radhyd_to_monitor          : transfer variables and monitor energy and lepton number conservation
!  radhyd_to_eos_x_reset      : transfer variables and reset the EOS tables along the radial rays
!  radhyd_to_nu_stress_x      : transfer variables and compute the x-component of the neutrino stress
!  radhyd_to_evh1_x_lagr      : transfer variables and perform the x-Lagrangian hydro
!  set_final_radial_grid      : set the final radial grid
!  radhyd_to_grav_l           : transfer variables and compute the gravitational potential and accelerations after Lagrangian update
!  eramp                      : deposit energy to simulate an explosion
!  radhyd_to_nuclear          : transfer variables and perform the nuclear abundance update
!  radhyd_to_nu_e_advct_x     : transfer variables and perform the neutrino advection in energy for the x-hydro
!  radhyd_to_equilibrate_x    : equilibrates neutrinos to matter at high densities
!  radhyd_to_transport        : transfer variables and perform the neutrino source and transport step
!  radhyd_to_remap_x          : transfer variables and perform the x-remap
!  radhyd_to_remap_x_e        : transfer variables and perform the x-remap of the total energy
!  radhyd_to_edit             : transfer variables and perform the 1D edits
!  radhyd_to_edit_2D          : transfer variables and perform the 2D edits
!  radhyd_to_edit_Global      : transfer variables and perform the Global edits
!  radhyd_to_edit_HDF         : transfer variables and write HDF files
!  radhyd_to_restart          : transfer variables and write restart files
!  radhyd_to_terminate        : transfer variables and terminate simulation according to criteria
!
!    Input arguments:
!  nx                         : x-array extent
!  ny                         : y-array extent
!  nz                         : z-array extent
!  nez                        : neutrino energy array extent
!  nnu                        : neutrino flavor array extent
!  nnc                        : composition array extent
!  ndim                       : number of spatial dimensions of the siulation
!  ij_ray_min                 : minimum j-index of a radial ray
!  ij_ray_max                 : maximum j-index of a radial ray
!  ij_ray_dim                 : number of y-zones on a processor before swapping with y
!  ik_ray_min                 : minimum k-index of a radial ray
!  ik_ray_max                 : maximum k-index of a radial ray
!  ik_ray_dim                 : number of z-zones on a processor before swapping with z
!  i_edit                     : edit parameter
!  time                       : elapsed time
!  dtnph                      : time step
!  nu_equil                   : neutrino equilibration flag
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  bomb_module, edit_module, evh1_sweep, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE bomb_module, ONLY : bomb_time, t_start_bomb, e_bomb, jexpl_min, jexpl_max
USE edit_module, ONLY : nlog
USE evh1_sweep, ONLY : nmin, nmax, sweep
USE prb_cntl_module, ONLY : i_bomb

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(in)    :: nu_equil   ! neutrino equilibration flag

INTEGER, INTENT(in)              :: nx         ! x-array extent
INTEGER, INTENT(in)              :: ny         ! y-array extent
INTEGER, INTENT(in)              :: nz         ! z-array extent
INTEGER, INTENT(in)              :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu        ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc        ! composition array extent
INTEGER, INTENT(in)              :: ndim       ! number of spatial dimensions
INTEGER, INTENT(in)              :: ij_ray_min ! minimum j-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_max ! maximum j-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_min ! minimum k-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray_max ! maximum k-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray_dim ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_edit     ! edit parameter

REAL(KIND=double), INTENT(in)    :: time       ! elapsed time
REAL(KIND=double), INTENT(in)    :: dtnph      ! time step
REAL(KIND=double)                :: tmax_bomb  ! time to start adding energy

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER(len=2)                 :: reset_comp_eos  ! composition EOS reset flag

INTEGER                          :: ij_ray          ! j-index of a radial ray
INTEGER                          :: ik_ray          ! k-index of a radial ray

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                     \\\\\ RESET PARAMETERS /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling parameter_reset, ij_ray_dim, ik_ray_dim=',2i4)") &
& ij_ray_dim, ik_ray_dim

CALL parameter_reset( ij_ray_dim, ik_ray_dim )

!-----------------------------------------------------------------------
!
!                    \\\\\ STORE INITIAL GRID /////
!
!-----------------------------------------------------------------------

WRITE (nlog,*) ' Calling store_int_grid'

CALL store_int_grid

!-----------------------------------------------------------------------
!
!       \\\\\ STORE INITIAL RAIDAL VARIABLES FOR LATER USE /////
!
!-----------------------------------------------------------------------

WRITE (nlog,*) ' Calling store_int_radial_var'

CALL store_int_radial_var

!-----------------------------------------------------------------------
!
!           \\\\\ NEUTRINO ENERGY DENSITIES AND FLUXES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_nu_energy_flux, ij_ray_dim, &
&  ik_ray_dim, nx, nez, nnu=',5i4)") &
&  ij_ray_dim, ik_ray_dim, nx, nez, nnu

CALL radhyd_to_nu_energy_flux( ij_ray_dim, ik_ray_dim, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE ANGULAR AVERAGES OF QUANTITIES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling angular_ave, nx, ny, nz, ij_ray_dim, &
&  ik_ray_dim=',5i4)") &
&  nx, ny, nz, ij_ray_dim, ik_ray_dim

CALL angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )

!-----------------------------------------------------------------------
!
!               \\\\\ MONITOR CONSERVED QUANTITIES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_monitor, ij_ray_dim, &
&  ik_ray_dim, nx, ny, nz, nez, nnu, nnc=',8i4)") &
&  ij_ray_dim, ik_ray_dim, nx, ny, nz, nez, nnu, nnc

CALL radhyd_to_monitor( ij_ray_dim, ik_ray_dim, nx, ny, nz, nez, &
& nnu, nnc )

!-----------------------------------------------------------------------
!
!       ||||| THE LOOP 1 OVER THE ij_ray_dim * ik_ray_dim |||||
!       |||||  RADIAL RAYS ON EACH PROCESSOR BEGINS HERE  |||||
!
!-----------------------------------------------------------------------

sweep                   = 'x'
DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim

    WRITE (nlog,"(' Starting first ij_ray, ik_ray loop, ik_ray, ik_ray_dim, &
&    ij_ray, ij_ray_dim',4i4)") ik_ray, ik_ray_dim, ij_ray, ij_ray_dim

!........Artificially increase internal energy to promote explosion.....

    IF ( i_bomb == 1 ) THEN
      IF ( ( time >= t_start_bomb)  .AND. ( time <= tmax_bomb ) ) THEN
        CALL eramp( e_bomb, bomb_time, dtnph, jexpl_min, jexpl_max, &
&        ij_ray, ik_ray )
      END IF ! time >= t_start_bomb
    END IF ! i_bomb == 1

!-----------------------------------------------------------------------
!
!                    \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_eos_x_reset 1, &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

    reset_comp_eos      = 'ye'
    CALL radhyd_to_eos_x_reset( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!                    \\\\\ NUCLEAR BURN STEP /////
!
!  Load abundance variables and perform the nuclear burn. Return updated
!   variables to radial_ray_module.
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_nuclear, nx, ij_ray_dim, ik_ray_dim, &
&    ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

    CALL radhyd_to_nuclear( nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, &
&    nnc )

!-----------------------------------------------------------------------
!
!                 \\\\\ NEUTRINO TRANSPORT STEP /////
!
!  Load radiation variables and the state variables and perform the
!   neutrino transport step. Return updated variables to radial_ray_module.
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_transport, nx, ij_ray_dim, ik_ray_dim, &
&    ij_ray, ik_ray, nez, nnu, nnc=',8i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nez, nnu, nnc

    CALL radhyd_to_transport( nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, &
&    nez, nnu, nnc )

!-----------------------------------------------------------------------
!
!                    \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_eos_x_reset 3, nx, ij_ray_dim, &
&    ik_ray_dim, ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP 1 OVER THE ij_ray_dim * ik_ray_dim |||||
!        |||||   RADIAL RAYS ON EACH PROCESSOR ENDS HERE   |||||
!
!-----------------------------------------------------------------------

  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

!-----------------------------------------------------------------------
!
!           \\\\\ NEUTRINO ENERGY DENSITIES AND FLUXES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_nu_energy_flux, ij_ray_dim, &
&  ik_ray_dim, nx, nez, nnu=',5i4)") &
&  ij_ray_dim, ik_ray_dim, nx, nez, nnu

CALL radhyd_to_nu_energy_flux( ij_ray_dim, ik_ray_dim, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE ANGULAR AVERAGES OF QUANTITIES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling angular_ave, nx, ny, nz, ij_ray_dim, ik_ray_dim=',5i4)") &
& nx, ny, nz, ij_ray_dim, ik_ray_dim

CALL angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )

!-----------------------------------------------------------------------
!
!    \\\\\ COMPUTE GRAVITY AND GENERAL RELATIVISTIC VARIABLES /////
!
!  Gravity may have changed from call to radhyd_to_grav_f in preceding
!   cycle due to y and z sweeps.
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_grav_l, nx, ny, nz, ij_ray_dim, &
& ik_ray_dim=',5i4)") nx, ny, nz, ij_ray_dim, ik_ray_dim

CALL radhyd_to_grav_i( nx, ny, nz, ij_ray_dim, ik_ray_dim )

!-----------------------------------------------------------------------
!
!    \\\\\ UPDATE NEUTRINO DISTRIBUTION FROM CHAINGING LAPSE /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_nu_agr_e_advct_x, &
& nx, ij_ray_dim, ik_ray_dim, nez, nnu =',5i4)") &
& nx, ij_ray_dim, ik_ray_dim, nez, nnu

CALL radhyd_to_nu_agr_e_advct_x( nx, ij_ray_dim, ik_ray_dim, nez, nnu )

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP 2 OVER THE ij_ray_dim * ik_ray_dim |||||
!        |||||  RADIAL RAYS ON EACH PROCESSOR BEGINS HERE  |||||
!
!-----------------------------------------------------------------------

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!
!                    \\\\\ NEUTRINO X-STRESS /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_nu_stress_x, ij_ray, ik_ray, &
&    ij_ray_dim, ik_ray_dim, nx, nez, nnu=',7i4)") &
&    ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, nx, nez, nnu

    CALL radhyd_to_nu_stress_x( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&    nx, nez, nnu )

!-----------------------------------------------------------------------
!
!                 \\\\\ X-LAGRANGIAN HYDRO STEP /////
!
!  Load EVH1 Variables and perform EVH1 Lagrangian hydro along
!   x-direction.
!  Return updated variables to radial_ray_module.
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_evh1_x_lag, nx, ny, nz, ij_ray, &
&    ik_ray, ij_ray_dim, ik_ray_dim, nez, nnu=',9i4)") &
&    nx, ny, nz, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, nez, nnu

    CALL radhyd_to_evh1_x_lagr( nx, ny, nz, ij_ray, ik_ray, ij_ray_dim, &
&    ik_ray_dim, nez, nnu )

!-----------------------------------------------------------------------
!
!              \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!  Load radiation variables and the state variables before and after the 
!   Lagrangian hydro steo, and perform the neutrino energy advection step.
!  Return updated variables to radial_ray_module.
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_nu_e_advct_x, nx, ij_ray_dim, &
&    ik_ray_dim, ij_ray, ik_ray, nez, nnu=',7i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nez, nnu

    CALL radhyd_to_nu_e_advct_x( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nez, nnu )

!-----------------------------------------------------------------------
!
!                    \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_eos_x_reset 2, nx, ij_ray_dim, &
&    ik_ray_dim, ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!               \\\\\ NEUTRINO EQUILIBRATION STEP /////
!
!  Equilibrate neutrinos to matter above a specified density conserving
!   entropy
!-----------------------------------------------------------------------

    IF ( nu_equil == 'ye' ) THEN
    WRITE (nlog,"(' Calling radhyd_to_equilibrate_x, nx, nez, nnu, &
&    ij_ray_dim, ik_ray_dim, ij_ray, ik_ray=',7i4)") &
&    nx, nez, nnu, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray

      CALL radhyd_to_equilibrate_x( nx, nez, nnu, ij_ray_dim, ik_ray_dim, &
&      ij_ray, ik_ray )
    END IF ! nu_equil == 'ye'

!-----------------------------------------------------------------------
!
!       ||||| THE LOOP 3 OVER THE ij_ray_dim * ik_ray_dim |||||
!       |||||  RADIAL RAYS ON EACH PROCESSOR BEGINS HERE  |||||
!
!-----------------------------------------------------------------------

  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

!-----------------------------------------------------------------------
!
!                      \\\\\ SET FINAL GRID /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling set_final_radial_grid, ij_ray_dim, ik_ray_dim, nx=',3i4)") &
& ij_ray_dim, ik_ray_dim, nx

CALL set_final_radial_grid( ij_ray_dim, ik_ray_dim, nx )

WRITE (nlog,"(' Calling radhyd_to_regrid_final_grid, nx, ij_ray_dim, ny, & 
& ik_ray_dim, nz=',5i4)") &
& nx, ij_ray_dim, ny, ik_ray_dim, nz

CALL radhyd_to_regrid_final_grid( nx, ij_ray_dim, ny, ik_ray_dim, nz )

!-----------------------------------------------------------------------
!
!       ||||| THE LOOP 3 OVER THE ij_ray_dim * ik_ray_dim |||||
!       |||||  RADIAL RAYS ON EACH PROCESSOR BEGINS HERE  |||||
!
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!
!                    \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_eos_x_reset 2, nx, ij_ray_dim, &
&    ik_ray_dim, ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!                       \\\\\ X-REMAP STEP /////!                       
!
!  Load updated variables and the x-Lagrangian altered grid the, x-remap
!   all variables except the energy to the final grid. Energy is remapped
!   after the gravitational potential has been recalculated
!  Return updated variables to radial_ray_module.
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_remap_x, nx, ij_ray_dim, &
&    ik_ray_dim, ij_ray, ik_ray, nez, nnu, nnc=',8i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nez, nnu, nnc

    CALL radhyd_to_remap_x( nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, &
&    nez, nnu, nnc )

    WRITE (nlog,"(' Calling radhyd_to_eos_x_reset 4, nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!               \\\\\ NEUTRINO EQUILIBRATION STEP /////
!
!  Equilibrate neutrinos to matter above a specified density conserving
!   entropy
!-----------------------------------------------------------------------

    IF ( nu_equil == 'ye' ) THEN
    WRITE (nlog,"(' Calling radhyd_to_equilibrate_x, nx, nez, nnu, &
&    ij_ray_dim, ik_ray_dim, ij_ray, ik_ray=',7i4)") &
&    nx, nez, nnu, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray

      CALL radhyd_to_equilibrate_x( nx, nez, nnu, ij_ray_dim, ik_ray_dim, &
&      ij_ray, ik_ray )
    END IF

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP 3 OVER THE ij_ray_dim * ik_ray_dim |||||
!        |||||   RADIAL RAYS ON EACH PROCESSOR ENDS HERE   |||||
!
!-----------------------------------------------------------------------

  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

!-----------------------------------------------------------------------
!
!           \\\\\ NEUTRINO ENERGY DENSITIES AND FLUXES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_nu_energy_flux, ij_ray_dim, &
&  ik_ray_dim, nx, nez, nnu=',5i4)") &
&  ij_ray_dim, ik_ray_dim, nx, nez, nnu

CALL radhyd_to_nu_energy_flux( ij_ray_dim, ik_ray_dim, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE ANGULAR AVERAGES OF QUANTITIES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling angular_ave, nx, ny, nz, ij_ray_dim, &
& ik_ray_dim=',5i4)") &
& nx, ny, nz, ij_ray_dim, ik_ray_dim

CALL angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE GENERAL RELATIVISTIC VARIABLES /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_grav_f, nx, ny, nz, ij_ray_dim, &
& ik_ray_dim =',5i4)") nx, ny, nz, ij_ray_dim, ik_ray_dim

CALL radhyd_to_grav_f( nx, ny, nz, ij_ray_dim, ik_ray_dim )

!-----------------------------------------------------------------------
!
!    \\\\\ UPDATE NEUTRINO DISTRIBUTION FROM CHAINGING LAPSE /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_nu_agr_e_advct_x, &
& nx, ij_ray_dim, ik_ray_dim, nez, nnu =',5i4)") &
& nx, ij_ray_dim, ik_ray_dim, nez, nnu

CALL radhyd_to_nu_agr_e_advct_x( nx, ij_ray_dim, ik_ray_dim, nez, nnu )

!-----------------------------------------------------------------------
!
!       ||||| THE LOOP 3 OVER THE ij_ray_dim * ik_ray_dim |||||
!       |||||  RADIAL RAYS ON EACH PROCESSOR BEGINS HERE  |||||
!
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim


!-----------------------------------------------------------------------
!
!                 \\\\\ X-REMAP THE TOTAL ENERGY /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_remap_x_e, nx, ij_ray_dim, &
&    ik_ray_dim, ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

    CALL radhyd_to_remap_x_e( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nnc )

!-----------------------------------------------------------------------
!
!       ||||| THE LOOP 3 OVER THE ij_ray_dim * ik_ray_dim |||||
!       |||||   RADIAL RAYS ON EACH PROCESSOR ENDS HERE   |||||
!
!-----------------------------------------------------------------------

  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT 1D /////
!
!  Load updated variables and test criteria for a 1D edit. Perform an
!   edit in accordance with the criteria satisfied.
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_edit, ij_ray_min, ij_ray_max, ij_ray_dim, &
& ik_ray_min, ik_ray_max, ik_ray_dim, i_edit, nx, ny, nz, nez, nnu, nnc=',13i4)") &
& ij_ray_min, ij_ray_max, ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, &
& i_edit, nx, ny, nz, nez, nnu, nnc

CALL radhyd_to_edit( ij_ray_min, ij_ray_max, ij_ray_dim, ik_ray_min, &
& ik_ray_max, ik_ray_dim, i_edit, nx, ny, nz, nez, nnu, nnc )

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT MD /////
!
!  Load updated variables and test criteria for a MD edit. Perform an
!   edit in accordance with the criteria satisfied.
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_edit_2D, nx, nez, nnu, ij_ray_dim, ny, &
& ik_ray_dim, nz, nnc=',8i4)") &
& nx, nez, nnu, ij_ray_dim, ny, ik_ray_dim, nz, nnc

CALL radhyd_to_edit_MD( nx, nez, nnu, ij_ray_dim, ny, ik_ray_dim, &
& nz, nnc )

!-----------------------------------------------------------------------
!
!                       \\\\\ EDIT GLOBAL /////
!
!  Load updated variables and test criteria for a Global edit. Perform an
!   edit in accordance with the criteria satisfied.
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_edit_Global, nx, nez, nnu, nnc, &
& ij_ray_dim, ny, ik_ray_dim, nz=',8i4)") &
& nx, nez, nnu, nnc, ij_ray_dim, ny, ik_ray_dim, nz

CALL radhyd_to_edit_Global( nx, nez, nnu, nnc, ij_ray_dim, ny, &
& ik_ray_dim, nz )

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT HDF /////
!
!  Load updated variables and test criteria for a MD edit. Perform an
!   edit in accordance with the criteria satisfied.
!-----------------------------------------------------------------------

!  CALL radhyd_to_edit_HDF( i_ray_min, i_ray_max, i_ray_dim, nx, ny, nez, nnu, nnc )

!-----------------------------------------------------------------------
!
!                         \\\\\ RESTART /////
!
!  Examine criteria for writing a restart file, and do so if the
!   criteria are met.
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_restart, ndim, nx, nez, nnu=',4i4)") &
& ndim, nx, nez, nnu

CALL radhyd_to_restart( ndim, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!                        \\\\\ TERMINATE/////
!
!  Examine criteria for terminating the simulation, and write edit and
!   restart files if the criteria are met.
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_terminate, ij_ray_min, ij_ray_max, &
& ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, nx, nez, nnu, nnc, ndimm=',11i4)") &
& ij_ray_min, ij_ray_max, ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, &
& nx, nez, nnu, nnc, ndim

CALL radhyd_to_terminate( ij_ray_min, ij_ray_max, ij_ray_dim, ik_ray_min, &
& ik_ray_max, ik_ray_dim, nx, nez, nnu, nnc, ndim )

RETURN
END SUBROUTINE radhyd_x_W
