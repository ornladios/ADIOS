SUBROUTINE radhyd_y_W( jmin, jmax, nx, ny, nz, nez, nnu, nnc, ndim, &
& ij_ray_dim, ik_ray_dim, j_ray_dim, n_proc, n_proc_y, n_proc_z, dtnph, &
& nu_equil, sub_cy_yz )
!-----------------------------------------------------------------------
!
!    File:         radhyd_y_W
!    Module:       radhyd_y_W
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/28/07
!
!    Purpose:
!      To direct the computation of the y-evolution.
!
!    Subprograms called:
!  radhyd_to_transpose_x_y     : transpose variables for the y-evolution
!  radhyd_to_eos_y_reset       : transfer variables and reset the EOS tables along the y (angular) rays
!  store_int_angular_var       : store the initial angular variables
!  radhyd_to_hydro_y_time_step : compute the y-hydro time step
!  radhyd_to_nu_stress_y       : transfer variables and compute the y-component of the neutrino stress
!  radhyd_to_evh1_y_lagr       : transfer variables and perform the y-Lagrangian hydro
!  radhyd_to_nu_e_advct_y      : transfer variables and perform the neutrino advection in energy for the y-hydro
!  radhyd_to_remap_y           : transfer variables and perform the y-remap
!  radhyd_to_equilibrate_x     : equilibrates neutrinos to matter at high densities
!  radhyd_to_equilibrate_y     : equilibrates neutrinos to matter at high densities
!  radhyd_to_transpose_y_x     : transpose variables for the x-evolution
!  radhyd_to_eos_x_reset       : transfer variables and reset the EOS tables along the x (radial) rays
!  set_final_angular_grid      : set the angular grid at the end of the time step
!
!    Input arguments:
!  jmin                        : minimum y-index
!  jmax                        : maximum y-index
!  nx                          : x-array extent
!  ny                          : y-array extent
!  nz                          : z-array extent
!  nez                         : neutrino energy array extent
!  nnu                         : neutrino flavor array extent
!  nnc                         : composition array extent
!  ndim                        : number of spatial dimensions of the siulation
!  ij_ray_dim                  : number of y-zones on a processor before swapping with y
!  ik_ray_dim                  : number of z-zones on a processor before swapping with z
!  j_ray_dim                   : number of radial zones on a processor after swapping with y
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

USE edit_module, ONLY : nlog
USE evh1_sweep, ONLY : sweep
USE parallel_module, ONLY : myid, myid_y, ierr

USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(in)    :: nu_equil        ! neutrino equilibration flag
CHARACTER (len=2), INTENT(in)    :: sub_cy_yz       ! subcycle flag

INTEGER, INTENT(in)              :: jmin            ! minimum y-index
INTEGER, INTENT(in)              :: jmax            ! maximum y-index
INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: ny              ! y-array extent
INTEGER, INTENT(in)              :: nz              ! z-array extent
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc             ! composition array extent
INTEGER, INTENT(in)              :: ndim            ! number of spatial dimensions
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: j_ray_dim       ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: n_proc          ! number of processors assigned to the run
INTEGER, INTENT(in)              :: n_proc_y        ! the number of processors assigned to the y-zones
INTEGER, INTENT(in)              :: n_proc_z        ! the number of processors assigned to the z-zones

REAL(KIND=double)                :: dtnph           ! time step

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER(len=2)                 :: reset_comp_eos  ! composition EOS reset flag

LOGICAL                          :: print_y = .true.

INTEGER                          :: ij_ray          ! j-index of a radial ray
INTEGER                          :: ik_ray          ! k-index of a radial ray
INTEGER                          :: ji_ray          ! x (radial) index of a specific y (angular) ray
INTEGER                          :: jk_ray          ! z (azimuthal) index of a specific y (angular) ray
INTEGER                          :: i_radial        ! the unshifted radial zone (angular ray) corresponding to ji_ray
INTEGER                          :: j_radial        ! the shifted radial zone (angular ray) corresponding to ji_ray

INTEGER                          :: i_time          ! subcycling index
INTEGER                          :: i_time_max      ! maximum subcycling index for subcycling
INTEGER                          :: i_time_lim      ! maximum subcycling index

REAL(KIND=double)                :: dt_y_hydro      ! time step for the y-sweep hydro
REAL(KIND=double)                :: dtnph_s         ! temporary storage for dtnph
REAL(KIND=double)                :: d_time          ! accumulated time for angular array subcycling
REAL(KIND=double), PARAMETER     :: dt_tol = 1.d-10 ! time step for the y-sweep hydro

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  201 FORMAT (' i_time=',i4,' = i_time_max=',i4,',  dtnph=',es11.3,' dtnph_s=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!              \\\\\ TRANSPOSE VARIABLES FOR Y-SWEEP /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling radhyd_to_transpose_x_y,  nx, ij_ray_dim, &
& ik_ray_dim, ny, j_ray_dim, nz, nez, nnu, nnc, n_proc, n_proc_y, &
& n_proc_z=',12i4)") &
& nx, ij_ray_dim, ik_ray_dim, ny, j_ray_dim, nz, nez, nnu, nnc, n_proc, &
& n_proc_y, n_proc_z

CALL radhyd_to_transpose_x_y( nx, ij_ray_dim, ik_ray_dim, ny, &
& j_ray_dim, nz, nez, nnu, nnc, n_proc, n_proc_y, n_proc_z )

!-----------------------------------------------------------------------
!
!                \\\\\ SET THE FINAL ANGULAR GRID /////
!
!-----------------------------------------------------------------------

WRITE (nlog,"(' Calling set_final_angular_grid')") 

CALL set_final_angular_grid

!-----------------------------------------------------------------------
!
!       |||||    THE LOOP OVER THE j_ray_dim * ik_ray_dim   |||||
!       |||||   ANGULAR RAYS ON EACH PROCESSOR BEGINS HERE  |||||
!
!-----------------------------------------------------------------------

sweep                 = 'y'
DO jk_ray = 1,ik_ray_dim
  DO ji_ray = 1,j_ray_dim
    i_radial          = j_ray_dim * myid_y + ji_ray
    j_radial          = i_radial + 1

!-----------------------------------------------------------------------
!
!               \\\\\ UPDATE Y-ARRAY EOS TABLES /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_eos_y_reset 1, ny, ji_ray, &
&    jk_ray, j_ray_dim, ik_ray_dim, nnc=',6i4)") &
&    ny, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, nnc

    reset_comp_eos    = 'ye'
    CALL radhyd_to_eos_y_reset(  ny, ji_ray, jk_ray, j_ray_dim, &
&    ik_ray_dim, nnc, reset_comp_eos  )

!-----------------------------------------------------------------------
!
!               \\\\\ COMPUTE TIME STEP FOR Y-SWEEP /////
!
!-----------------------------------------------------------------------

    WRITE (nlog,"(' Calling radhyd_to_hydro_y_time_step, jmin, jmax, &
&    ji_ray, jk_ray, j_ray_dim, ik_ray_dim, dt_y_hydro=',6i4,es11.3)") &
&    jmin, jmax, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, dt_y_hydro

    CALL radhyd_to_hydro_y_time_step( jmin, jmax, ji_ray, jk_ray, &
&    j_ray_dim, ik_ray_dim, dt_y_hydro )

!-----------------------------------------------------------------------
!
!          \\\\\ SUBCYCLE Y-SWEEP IF DTIME_HYDRO < DTNPH /////
!
!-----------------------------------------------------------------------

    dtnph_s           = dtnph
    d_time            = zero
    IF ( dt_y_hydro >= dtnph ) THEN
      i_time_max      = 1
    ELSE
      i_time_max      = 1000
      dtnph           = dt_y_hydro
    END IF

    IF ( sub_cy_yz == 'ye' ) THEN
      i_time_lim      = i_time_max
    ELSE
      i_time_lim      = 1
    END IF

    DO i_time = 1,i_time_lim
      dtnph           = DMIN1( dtnph, dtnph_s - d_time )
      d_time          = d_time + dtnph

      WRITE (nlog,"(' i_time=',i4)") i_time

!-----------------------------------------------------------------------
!
!      \\\\\ STORE INITIAL ANGULAR VARIABLES FOR LATER USE /////
!
!-----------------------------------------------------------------------

      WRITE (nlog,"(' Calling store_int_angular_var, ( ji_ray, jk_ray, ny )=',3i4)") &
&      ji_ray, jk_ray, ny

      CALL store_int_angular_var( ji_ray, jk_ray, ny )

!-----------------------------------------------------------------------
!
!                  \\\\\ NEUTRINO Y-STRESS /////
!
!-----------------------------------------------------------------------

      WRITE (nlog,"(' Calling radhyd_to_nu_stress_y, ji_ray, jk_ray, &
&      j_ray_dim, ik_ray_dim, i_radial, j_radial, nx, ny, nez, nnu=',10i4)") &
&      ji_ray, jk_ray, j_ray_dim, ik_ray_dim, i_radial, j_radial, nx, ny, nez, nnu

      CALL radhyd_to_nu_stress_y( ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
&      i_radial, j_radial, nx, ny, nez, nnu )

!-----------------------------------------------------------------------
!
!                \\\\\ Y-LAGRANGIAN HYDRO STEP /////
!
!  Load EVH1 Variables and perform EVH1 Lagrangian hydro along
!   y-direction.
!  Return updated variables to radial_ray_module.
!-----------------------------------------------------------------------

      WRITE (nlog,"(' Calling radhyd_to_evh1_y_lagr, nx, ny, nz, ji_ray, &
&      jk_ray, j_ray_dim, ik_ray_dim, i_radial, nez, nnu=',10i4)") &
&      nx, ny, nz, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, i_radial, nez, nnu

      CALL radhyd_to_evh1_y_lagr( nx, ny, nz, ji_ray, jk_ray, &
&      j_ray_dim, ik_ray_dim, i_radial, nez, nnu )

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------
 
      WRITE (nlog,"(' Calling radhyd_to_eos_y_reset 2, ny, ji_ray, &
&      jk_ray, j_ray_dim, ik_ray_dim, nnc=',6i4)") &
&      ny, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, nnc

      reset_comp_eos    = 'no'
      CALL radhyd_to_eos_y_reset( ny, ji_ray, jk_ray, j_ray_dim, &
&      ik_ray_dim, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!  Load radiation variables and the state variables before and after
!   the Lagrangian hydro steo, and perform the neutrino energy advection
!   step.
!  Return updated variables to radial_ray_module.
!-----------------------------------------------------------------------

      WRITE (nlog,"(' Calling radhyd_to_nu_e_advct_y, nx, ny, ji_ray, &
&      jk_ray, j_ray_dim, ik_ray_dim, i_radial, j_radial, nez, nnu=',10i4)") &
&      nx, ny, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, i_radial, j_radial, nez, nnu

      CALL radhyd_to_nu_e_advct_y( nx, ny, ji_ray, jk_ray, j_ray_dim, &
&      ik_ray_dim, i_radial, j_radial, nez, nnu )

!-----------------------------------------------------------------------
!
!                    \\\\\ Y-REMAP STEP /////
!
!  Load updated variables and the grid before and after the y-Lagrangian
!   step, and perform the y-remap step. Return updated variables to
!   angular_ray_module.
!-----------------------------------------------------------------------
 
      WRITE (nlog,"(' Calling radhyd_to_remap_y, nx, ny, ji_ray, jk_ray, &
&      j_ray_dim, ik_ray_dim, i_radial, nez, nnu, nnc=',10i4)") &
&      nx, ny, ji_ray, jk_ray, j_ray_dim, i_radial, ik_ray_dim, nez, nnu, nnc

      CALL radhyd_to_remap_y( nx, ny, ji_ray, jk_ray, j_ray_dim, &
&      ik_ray_dim, i_radial, nez, nnu, nnc )
 
      WRITE (nlog,"(' Calling radhyd_to_eos_y_reset 3, ny, ji_ray, jk_ray, &
&      j_ray_dim, ik_ray_dim, nnc=',6i4)") &
&      ny, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, nnc

      reset_comp_eos    = 'no'
      CALL radhyd_to_eos_y_reset( ny, ji_ray, jk_ray, j_ray_dim, &
&      ik_ray_dim, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!              \\\\\ NEUTRINO EQUILIBRATION STEP /////
!
!-----------------------------------------------------------------------

      IF ( nu_equil == 'ye'  .and.  i_time > 1 ) THEN
        WRITE (nlog,"(' Calling radhyd_to_equilibrate_y nx, ny, nez, &
&        nnu, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, i_radial=',9i4)") &
&        nx, ny, nez, nnu, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, i_radial

        CALL radhyd_to_equilibrate_y( nx, ny, nez, nnu, ji_ray, &
&        jk_ray, j_ray_dim, ik_ray_dim, i_radial )
      END IF ! nu_equil == 'ye'

!-----------------------------------------------------------------------
!
!              \\\\\ END SUBCYCLE Y-SWEEP LOOP /////
!
!-----------------------------------------------------------------------

      IF ( d_time > dtnph_s - dt_tol ) EXIT
      IF ( i_time == i_time_max ) THEN
        print  *,' i_time,i_time_max,dtnph,dtnph_s',i_time,i_time_max,dtnph,dtnph_s
        WRITE (nlog,201) i_time,i_time_max,dtnph,dtnph_s
        STOP
      END IF

    END DO ! i_time

    dtnph               = dtnph_s

!-----------------------------------------------------------------------
!
!       |||||    THE LOOP OVER THE j_ray_dim * ik_ray_dim   |||||
!       |||||    ANGULAR RAYS ON EACH PROCESSOR ENDS HERE   |||||
!
!-----------------------------------------------------------------------

  END DO ! ji_ray
END DO ! jk_ray

!-----------------------------------------------------------------------
!
!           \\\\\ TRANSPOSE VARIABLES BACK FOR X-SWEEP /////
!
!-----------------------------------------------------------------------
 
 WRITE (nlog,"(' Calling radhyd_to_transpose_y_x, nx, ij_ray_dim, &
& ik_ray_dim, ny, j_ray_dim, nez, nnu, nnc, n_proc, n_proc_y=',10i4)") &
& nx, ij_ray_dim, ik_ray_dim, ny, j_ray_dim, nez, nnu, nnc, n_proc, &
&  n_proc_y

CALL radhyd_to_transpose_y_x( nx, ij_ray_dim, ik_ray_dim, ny, &
& j_ray_dim, nez, nnu, nnc, n_proc, n_proc_y )

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

    WRITE (nlog,"(' Calling radhyd_to_eos_x_reset 5, nx, ij_ray_dim, &
&    ik_ray_dim, ij_ray, ik_ray, nnc=',6i4)") &
&    nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, nnc

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
      WRITE (nlog,"(' Calling radhyd_to_equilibrate_x, nx, nez, nnu, &
&      ij_ray_dim, ik_ray_dim, ij_ray, ik_ray=',7i4)") &
&      nx, nez, nnu, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray

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
END SUBROUTINE radhyd_y_W

