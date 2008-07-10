SUBROUTINE radhyd_z( kmin, kmax, nx, ny, nz, nez, nnu, nnc, ndim, &
& ij_ray_dim, ik_ray_dim, k_ray_dim, n_proc, n_proc_y, n_proc_z, dtnph, &
& nu_equil, sub_cy_yz )
!-----------------------------------------------------------------------
!
!    File:         radhyd_z
!    Module:       radhyd_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/28/07
!
!    Purpose:
!      To direct the computation of the z-evolution.
!
!    Subprograms called:
!  radhyd_to_transpose_x_z     : transpose variables for the z-evolution
!  radhyd_to_eos_z_reset       : transfer variables and reset the EOS tables along the z (angular) rays
!  store_int_azimuthal_var     : store the initial azimuthal variables
!  radhyd_to_hydro_z_time_step : compute the z-hydro time step
!  radhyd_to_nu_stress_z       : transfer variables and compute the y-component of the neutrino stress
!  radhyd_to_evh1_z_lagr       : transfer variables and perform the z-Lagrangian hydro
!  radhyd_to_nu_e_advct_z      : transfer variables and perform the neutrino advection in energy for the z-hydro
!  radhyd_to_remap_z           : transfer variables and perform the z-remap
!  radhyd_to_equilibrate_z     : equilibrates neutrinos to matter at high densities
!  radhyd_to_equilibrate_x     : equilibrates neutrinos to matter at high densities
!  radhyd_to_transpose_z_x     : transpose variables for the x-evolution
!  radhyd_to_eos_x_reset       : transfer variables and reset the EOS tables along the x (radial) rays
!  set_final_azimuthal_grid    : set the angular grid at the end of the time step
!
!    Input arguments:
!  kmin                        : minimum z-index
!  kmax                        : maximum z-index
!  nx                          : x-array extent
!  ny                          : y-array extent
!  nz                          : z-array extent
!  nez                         : neutrino energy array extent
!  nnu                         : neutrino flavor array extent
!  nnc                         : composition array extent
!  ndim                        : number of spatial dimensions of the siulation
!  ij_ray_dim                  : number of y-zones on a processor before swapping with y
!  ik_ray_dim                  : number of z-zones on a processor before swapping with z
!  k_ray_dim                   : number of radial zones on a processor after swapping with z
!  n_proc                      : number of processors assigned to the run
!  n_proc_y                    : the number of processors assigned to the y-zones
!  n_proc_z                    : the number of processors assigned to the z-zones
!  dtnph                       : time step
!  nu_equil                    : neutrino equilibration flag
!  sub_cy_yz                   : subcycle flag
!
!    Output arguments:
!        none
!      
!    Include files:
!  kind_module, numerical_module
!  edit_module, evh1_sweep, parallel_module
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE edit_module, ONLY : nlog, nprint
USE evh1_sweep, ONLY : sweep
USE parallel_module, ONLY : myid, myid_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(in)    :: nu_equil   ! neutrino equilibration flag
CHARACTER (len=2), INTENT(in)    :: sub_cy_yz  ! subcycle flag

INTEGER, INTENT(in)              :: kmin            ! minimum z-index
INTEGER, INTENT(in)              :: kmax            ! maximum z-index
INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: ny              ! y-array extent
INTEGER, INTENT(in)              :: nz              ! z-array extent
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc             ! composition array extent
INTEGER, INTENT(in)              :: ndim            ! number of spatial dimensions
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: k_ray_dim       ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: n_proc          ! number of processors assigned to the run
INTEGER, INTENT(in)              :: n_proc_y        ! the number of processors assigned to the y-zones
INTEGER, INTENT(in)              :: n_proc_z        ! the number of processors assigned to the z-zones

REAL(KIND=double)                :: dtnph           ! time step

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER(len=2)                 :: reset_comp_eos  ! composition EOS reset flag

INTEGER                          :: ij_ray          ! j-index of a radial ray
INTEGER                          :: ik_ray          ! k-index of a radial ray
INTEGER                          :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER                          :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray
INTEGER                          :: i_radial        ! the unshifted radial zone (angular ray) corresponding to ji_ray
INTEGER                          :: j_radial        ! the shifted radial zone (angular ray) corresponding to ji_ray

INTEGER                          :: i_time          ! subcycling index
INTEGER                          :: i_time_max      ! maximum subcycling index for subcycling
INTEGER                          :: i_time_lim      ! maximum subcycling index

REAL(KIND=double)                :: dt_z_hydro      ! time step for the z-sweep hydro
REAL(KIND=double)                :: dtnph_s         ! temporary storage for dtnph
REAL(KIND=double)                :: d_time          ! accumulated time for angular array subcycling
REAL(KIND=double), PARAMETER     :: dt_tol = 1.d-10 ! time step for the z-sweep hydro

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Need MPI version of radhyd_z since n_proc=',i4,' > 1')
  201 FORMAT (' i_time=',i4,' = i_time_max=',i4,',  dtnph=',es11.3,' dtnph_s=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Stop if n_proc > 1.
!-----------------------------------------------------------------------

IF ( n_proc > 1 ) THEN
  WRITE (nprint,101) n_proc
  WRITE (nlog,101) n_proc
  STOP
END IF ! n_proc > 1

!-----------------------------------------------------------------------
!
!              \\\\\ TRANSPOSE VARIABLES FOR Z-SWEEP /////
!
!-----------------------------------------------------------------------

CALL radhyd_to_transpose_x_z( nx, ij_ray_dim, ik_ray_dim, nz, &
& k_ray_dim, ny, nez, nnu, nnc, n_proc, n_proc_y, n_proc_z )

!-----------------------------------------------------------------------
!
!               \\\\\ SET THE FINAL AZIMUTHAL GRID /////
!
!-----------------------------------------------------------------------

CALL set_final_azimuthal_grid

!-----------------------------------------------------------------------
!
!      |||||    THE LOOP OVER THE ij_ray_dim * k_ray_dim    |||||
!      |||||  AZIMUTHAL RAYS ON EACH PROCESSOR BEGINS HERE  |||||
!
!-----------------------------------------------------------------------

sweep                 = 'z'
DO kj_ray = 1,ij_ray_dim
  DO ki_ray = 1,k_ray_dim
    i_radial          = k_ray_dim * myid_z + ki_ray
    j_radial          = i_radial + 1

!-----------------------------------------------------------------------
!
!               \\\\\ UPDATE Z-ARRAY EOS TABLES /////
!
!-----------------------------------------------------------------------

    reset_comp_eos    = 'ye'
    CALL radhyd_to_eos_z_reset( nz, ki_ray, kj_ray, ij_ray_dim, &
&    k_ray_dim, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!               \\\\\ COMPUTE TIME STEP FOR Z-SWEEP /////
!
!-----------------------------------------------------------------------

    CALL radhyd_to_hydro_z_time_step( kmin, kmax, ki_ray, kj_ray, &
&    ij_ray_dim, k_ray_dim, dt_z_hydro )

!-----------------------------------------------------------------------
!
!          \\\\\ SUBCYCLE Z-SWEEP IF DTIME_HYDRO < DTNPH /////
!
!-----------------------------------------------------------------------

    dtnph_s           = dtnph
    d_time            = zero
    IF ( dt_z_hydro >= dtnph ) THEN
      i_time_max      = 1
    ELSE
      i_time_max      = 1000
      dtnph           = dt_z_hydro
    END IF

    IF ( sub_cy_yz == 'ye' ) THEN
      i_time_lim      = i_time_max
    ELSE
      i_time_lim      = 1
    END IF

    DO i_time = 1,i_time_lim
      dtnph           = DMIN1( dtnph, dtnph_s - d_time )
      d_time          = d_time + dtnph

!-----------------------------------------------------------------------
!
!     \\\\\ STORE INITIAL AZIMUTHAL VARIABLES FOR LATER USE /////
!
!-----------------------------------------------------------------------

      CALL store_int_azimuthal_var( ki_ray, kj_ray, nz )

!-----------------------------------------------------------------------
!
!                  \\\\\ NEUTRINO Z-STRESS /////
!
!-----------------------------------------------------------------------

      CALL radhyd_to_nu_stress_z( ki_ray, kj_ray, ij_ray_dim, k_ray_dim, &
&      i_radial, j_radial, nx, ny, nz, nez, nnu  )

!-----------------------------------------------------------------------
!
!                \\\\\ Z-LAGRANGIAN HYDRO STEP /////
!
!  Load EVH1 Variables and perform EVH1 Lagrangian hydro along        
!   z-direction. Return updated variables to azimuthal_ray_module.
!-----------------------------------------------------------------------

      CALL radhyd_to_evh1_z_lagr( nx, ny, nz, ki_ray, kj_ray, &
&      ij_ray_dim, k_ray_dim, i_radial, nez, nnu )

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

      reset_comp_eos    = 'no'
      CALL radhyd_to_eos_z_reset( nz, ki_ray, kj_ray, ij_ray_dim, &
&      k_ray_dim, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!  Load radiation variables and the state variables before and after
!   the Lagrangian hydro steo, and perform the neutrino energy advection
!   step.
!  Return updated variables to azimuthal_ray_module.
!-----------------------------------------------------------------------

      CALL radhyd_to_nu_e_advct_z( nx, nz, ki_ray, kj_ray, ij_ray_dim, &
&      k_ray_dim, i_radial, j_radial, nez, nnu )

!-----------------------------------------------------------------------
!
!                    \\\\\ Z-REMAP STEP /////
!
!  Load updated variables and the grid before and after the z-Lagrangian
!   step, and perform the z-remap step. Return updated variables to 
!   azimuthal_ray_module.
!-----------------------------------------------------------------------

      CALL radhyd_to_remap_z( nx, nz, ki_ray, kj_ray, ij_ray_dim, &
&      k_ray_dim, i_radial, nez, nnu, nnc )

      reset_comp_eos    = 'no'
      CALL radhyd_to_eos_z_reset( nz, ki_ray, kj_ray, ij_ray_dim, &
&      k_ray_dim, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!              \\\\\ NEUTRINO EQUILIBRATION STEP /////
!
!-----------------------------------------------------------------------

      IF ( nu_equil == 'ye'  .and.  i_time > 1 ) THEN
        CALL radhyd_to_equilibrate_z( nx, nz, nez, nnu, ki_ray, kj_ray, &
&        ij_ray_dim, k_ray_dim, i_radial )
      END IF

!-----------------------------------------------------------------------
!
!              \\\\\ END SUBCYCLE Z-SWEEP LOOP /////
!
!-----------------------------------------------------------------------

      IF ( d_time > dtnph_s - dt_tol ) EXIT
      IF ( i_time == i_time_max ) THEN
        print  *,' i_time, i_time_max, dtnph, dtnph_s',             &
&        i_time, i_time_max, dtnph, dtnph_s
        WRITE (nlog,201) i_time,i_time_max,dtnph,dtnph_s
        STOP
      END IF

    END DO ! i_time

    dtnph               = dtnph_s

!-----------------------------------------------------------------------
!
!       |||||    THE LOOP OVER THE ij_ray_dim * k_ray_dim   |||||
!       |||||    ANGULAR RAYS ON EACH PROCESSOR ENDS HERE   |||||
!
!-----------------------------------------------------------------------

  END DO ! kj_ray
END DO ! ki_ray

!-----------------------------------------------------------------------
!
!           \\\\\ TRANSPOSE VARIABLES BACK FOR X-SWEEP /////
!
!-----------------------------------------------------------------------

CALL radhyd_to_transpose_z_x( nx, ij_ray_dim, ik_ray_dim, nz, k_ray_dim, &
& nez, nnu, nnc, n_proc, n_proc_z )

!-----------------------------------------------------------------------
!
!         ||||| THE LOOP OVER THE ij_ray_dim * ik_ray_dim |||||
!         ||||| RADIAL RAYS ON EACH PROCESSOR BEGINS HERE |||||
!
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    reset_comp_eos          = 'no'

    CALL radhyd_to_eos_x_reset( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
&    ik_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!              \\\\\ NEUTRINO EQUILIBRATION STEP /////
!
!  Equilibrate neutrinos to matter above a specified density conserving
!   entropy
!-----------------------------------------------------------------------

    IF ( nu_equil == 'ye' ) THEN
      CALL radhyd_to_equilibrate_x( nx, nez, nnu, ij_ray_dim, ik_ray_dim, &
&      ij_ray, ik_ray )
    END IF ! nu_equil == 'ye'

!-----------------------------------------------------------------------
!
!         ||||| THE LOOP OVER THE ij_ray_dim * ik_ray_dim |||||
!         |||||  RADIAL RAYS ON EACH PROCESSOR ENDS HERE  |||||
!
!-----------------------------------------------------------------------

  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

RETURN
END SUBROUTINE radhyd_z

