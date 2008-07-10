!=======================================================================
!
!    \\\\\\\\\\        B E G I N   P R O G R A M          //////////
!    //////////               R A D H Y D                 \\\\\\\\\\
!
!                            Developed by
!                         Stephen W. Bruenn
!             Florida Atlantic University, Boca Raton, FL
!
!=======================================================================

PROGRAM radhyd

!-----------------------------------------------------------------------
!
!    File:         radhyd
!    Module:       radhyd
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/04
!
!    Purpose:
!      To initiate and run the simulation.
!
!    Subprograms called:
!  read_pack_array_dimensions   : reads and packs array dimensions
!  unpack_array_dimenisons      : unpacks array dimensions and distributes to all processors
!  load_array_module            : loads array dimensions into array_module
!  initialize                   : initialize problem
!  initialize_grid              : set the initial values of the grid
!  rhobar_cal                   : calculate the mean density in an angular sector
!  radhyd_to_edit               : transfer variables and call edit routines
!  radhyd_to_monitor            : transfer variables and call the energy and lepton number conservation check
!  cycle                        : update the cycle number
!  parameter_reset              : change rezn or m_grid parameters according to criteria
!  store_int_grid               : store the initial grid
!  store_int_radial_var         : store the initial radial ray variables
!  radhyd_to_eos_x_reset        : transfer variables and reset the EOS tables along the radial rays
!  radhyd_to_evh1_x_lagr        : transfer variables and perform the x-Lagrangian hydro
!  set_final_radial_grid        : set the final radial grid
!  eramp                        : deposit energy to simulate an explosion
!  radhyd_to_nuclear            : transfer variables and perform the nuclear abundance update
!  radhyd_to_nu_e_advct_x       : transfer variables and perform the neutrino advection in energy for the x-hydro
!  radhyd_to_transport          : transfer variables and perform the neutrino source and transport step
!  radhyd_to_remap_x            : transfer variables and perform the x-remap
!  radhyd_to_edit_2D            : transfer variables and perform the 2D edits
!  radhyd_to_edit_HDF           : transfer variables and write HDF files
!  radhyd_to_restart            : transfer variables and write restart files
!  radhyd_to_terminate          : transfer variables and terminate simulation according to criteria
!  radhyd_to_transpose_x_y      : transpose variables for the y-evolution
!  radhyd_to_eos_y_reset        : transfer variables and reset the EOS tables along the angular rays
!  store_int_angular_var        : store the initial angular variables
!  radhyd_to_hydro_y_time_step  : compute the y-hydro time step
!  radhyd_to_evh1_y_lagr        : transfer variables and perform the y-Lagrangian hydro
!  radhyd_to_remap_y            : transfer variables and perform the y-remap
!  radhyd_to_transpose_y_x      : transpose variables for the x-evolution
!  radyhd_to_time_step_select   : transpose variables and compute the next time step
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  angular_ray_module, bomb_module, edit_module, eos_snc_x_module,
!  evh1_global, evh1_sweep, incrmnt_module, mdl_cnfg_module,
!  parallel_module, prb_cntl_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE angular_ray_module, ONLY : v_y, rho_y, t_y, ye_y
USE bomb_module, ONLY : bomb_time, t_start_bomb, e_bomb, jexpl_min, jexpl_max
USE edit_module, ONLY : data_path, log_path, nlog
USE eos_snc_x_module, ONLY: nse
USE evh1_global, ONLY : svel
USE evh1_sweep, ONLY : nmin, nmax, sweep
USE incrmnt_module, ONLY : dtmpmn, dtmpnn, dye, dyecnvt
USE mdl_cnfg_module, ONLY : jr_max
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY : i_bomb
USE radial_ray_module, ONLY : jmin, jmax, ndim, dtnph, time, sub_cy_yz, &
& nu_equil

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER(len=128), DIMENSION(3) :: c_path_data     ! character array of data path
CHARACTER(len=128)               :: log_file        ! character array of data path
CHARACTER(len=2)                 :: reset_comp_eos  ! composition EOS reset flag

INTEGER                          :: nx              ! x-array extent
INTEGER                          :: ny              ! y-array extent
INTEGER                          :: nz              ! z-array extent
INTEGER                          :: nez             ! neutrino energy array extent
INTEGER                          :: nnu             ! neutrino flavor array extent
INTEGER                          :: nnc             ! composition array extent
INTEGER                          :: max_12          ! max(nx,ny,nz)+12
INTEGER                          :: n_proc          ! number of processors assigned to the run
INTEGER                          :: i_ray_dim       ! number of x-rays to put on a processor
INTEGER                          :: j_ray_dim       ! number of y-rays to put on a processor

INTEGER, DIMENSION(20)           :: n_dim_data      ! integer array of array dimensions

INTEGER                          :: i               ! d_temperature index, x-coordinate index
INTEGER                          :: j               ! zone index, y-coordinate index
INTEGER                          :: i_ray           ! radial ray index
INTEGER                          :: i_ray_min       ! minimum index denoting a specific radial ray
INTEGER                          :: i_ray_max       ! maximum index denoting a specific radial ray
INTEGER                          :: j_ray           ! angular ray index
INTEGER                          :: n               ! neutrino flavor index
INTEGER                          :: i_radial        ! the unshifted radial zone (angular ray) corresponding to j_ray
INTEGER                          :: j_radial        ! the shifted radial zone (angular ray) corresponding to j_ray

INTEGER                          :: i_time          ! subcycling index
INTEGER                          :: i_time_max      ! maximum subcycling index for subcycling
INTEGER                          :: i_time_lim      ! maximum subcycling index

INTEGER                          :: i_edit          ! edit parameter
INTEGER                          :: istat           ! open-close file flag

REAL(KIND=double)                :: tmax_bomb       ! time to start adding energy
REAL(KIND=double)                :: dt_y_hydro      ! time step for the y-sweep hydro
REAL(KIND=double)                :: dtnph_s         ! temporary storage for dtnph
REAL(KIND=double)                :: d_time          ! accumulated time for angular array subcycling
REAL(KIND=double), PARAMETER     :: dt_tol = 1.d-10 ! time step for the y-sweep hydro

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Array dimensions have been read')
  105 FORMAT (' Array dimensions have been broadcast and unpacked')
  107 FORMAT (' Output data path is ',a128)
  109 FORMAT (' Simulation log path is ',a128)
  111 FORMAT (' Array modules have been loaded')
  201 FORMAT (' i_time=',i4,' = i_time_max=',i4,',  dtnph=',es11.3,' dtnph_s=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

myid              = 0

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZATION AND SET UP /////
!
!  Read in array dimensions and number of processors to be used,
!   broadcast to all processors.
!
!  Open the file to which the simulation log is to be directed.
!
!  Read in problem control keys and initial model, broadcast to
!   all processors. Initialize the problem.
!        
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN
  CALL read_pack_array_dimensions( c_path_data, n_dim_data )
END IF
WRITE (nlog,101)

CALL unpack_array_dimenisons( c_path_data, n_dim_data, nx, ny, nz, nez, &
& nnu, nnc, n_proc, n_proc_y, n_proc_z, ij_ray_dim, ik_ray_dim, j_ray_dim, &
& k_ray_dim, data_path, log_path )

IF ( log_path(1:1) == ' ' ) THEN
  nlog               = 6
ELSE
  nlog               = 49
  WRITE (log_file,'(a17,i4.4,a2)') '/Log_File/sim_log',myid,'.d'
  log_file           = TRIM(log_path)//TRIM(log_file)
  OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='old')
END IF

WRITE (nlog,105)
WRITE (nlog,107) data_path
WRITE (nlog,109) log_path

CALL load_array_module( nx, ny, nz, nez, nnu, nnc, n_proc, n_proc_y, &
& n_proc_z, ij_ray_dim, ik_ray_dim, j_ray_dim, k_ray_dim, max_12 )
WRITE (nlog,111)

CALL initialize( nx, ny, nz, nez, nnu, nnc, max_12, n_proc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )

i_edit               = 0
i_ray_min            = 1
i_ray_max            = i_ray_dim

CALL initialize_grid

CALL rhobar_cal

CALL radhyd_to_edit( i_ray_min, i_ray_max, i_ray_dim, i_edit, nx, nez, nnu, nnc )

CALL radhyd_to_monitor( i_ray_dim, nx, nez, nnu, nnc )

tmax_bomb   = t_start_bomb + bomb_time

!-----------------------------------------------------------------------
!
!             \\\\\ MAIN EVOLUTION LOOP /////
!             /////     BEGINS HERE     \\\\\
!
!-----------------------------------------------------------------------

WRITE (nlog,*) 'Beginning problem loop'

DO

  IF ( log_path(1:1) /= ' ' ) THEN
    CLOSE (unit=nlog,status='keep')
    OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='old',POSITION='append')
  END IF ! log_path(1:1) /= ' '

  WRITE (nlog,*) 'time =',time, 'dtnph =',dtnph
 
!........Update cycle number............................................

  CALL cycle( time, dtnph )

 
!........Check for parameter change.....................................

  CALL parameter_reset( i_ray_dim )

!........Initialize increment arrays....................................

  DO i_ray = 1,i_ray_dim
    DO i = 1,10
      DO j = 1,jr_max
        dtmpmn(j,i,i_ray) = zero
        dtmpnn(j,i,i_ray) = zero
      END DO
    END DO
  END DO

  DO i_ray = 1,i_ray_dim
    DO n = 1,nnu
      DO j = 1,jr_max
        dye(j,n,i_ray)    = zero
      END DO
    END DO
  END DO

  DO i_ray = 1,i_ray_dim
    DO j = 1,jr_max
      dyecnvt(j,i_ray)    = zero
    END DO
  END DO

  svel                    = zero

!-----------------------------------------------------------------------
!
!                  \\\\\ STORE INITIAL GRID /////
!
!-----------------------------------------------------------------------

  CALL store_int_grid

!-----------------------------------------------------------------------
!
!      \\\\\ STORE INITIAL RAIDAL VARIABLES FOR LATER USE /////
!
!-----------------------------------------------------------------------

  CALL store_int_radial_var( i_ray_dim, nx )

!-----------------------------------------------------------------------
!
!                  \\\\\ COMPUTE MEAN DENSITIES /////
!
!-----------------------------------------------------------------------

  CALL rhobar_cal

!-----------------------------------------------------------------------
!
!               \\\\\ MONITOR CONSERVED QUANTITIES /////
!
!-----------------------------------------------------------------------

  CALL radhyd_to_monitor( i_ray_dim, nx, nez, nnu, nnc )

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP 1 OVER THE i_ray_dim RADIAL RAYS |||||
!        |||||      ON EACH PROCESSOR BEGINS HERE        |||||
!
!-----------------------------------------------------------------------

  sweep                   = 'x'
  DO i_ray = 1,i_ray_dim

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    reset_comp_eos        = 'ye'
    CALL radhyd_to_eos_x_reset( nx, i_ray_dim, i_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!                  \\\\\ NEUTRINO X-STRESS /////
!
!-----------------------------------------------------------------------

    CALL radhyd_to_nu_stress_x( i_ray, i_ray_dim, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!                \\\\\ X-LAGRANGIAN HYDRO STEP /////
!
!        Load EVH1 Variables and perform EVH1 Lagrangian hydro
!         along x-direction.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

    CALL radhyd_to_evh1_x_lagr( nx, i_ray, i_ray_dim, ny, nz, nez, nnu )

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP 1 OVER THE i_ray_dim RADIAL RAYS |||||
!        |||||       ON EACH PROCESSOR ENDS HERE         |||||
!
!-----------------------------------------------------------------------

  END DO

!-----------------------------------------------------------------------
!
!                     \\\\\ SET FINAL GRID /////
!
!-----------------------------------------------------------------------

  CALL set_final_radial_grid( i_ray_dim, nx, ny )

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP 2 OVER THE i_ray_dim RADIAL RAYS |||||
!        |||||      ON EACH PROCESSOR BEGINS HERE        |||||
!
!-----------------------------------------------------------------------

  DO i_ray = 1,i_ray_dim

!........Artificially increase internal energy to promote explosion.....

    IF (i_bomb .eq. 1) THEN
      IF ((time >= t_start_bomb) .AND. (time <= tmax_bomb)) THEN
        CALL eramp( e_bomb, bomb_time, dtnph, jexpl_min, jexpl_max, i_ray )
      END IF
    END IF  

!-----------------------------------------------------------------------
!
!                   \\\\\ NUCLEAR BURN STEP /////
!
!        Load abundance variables and perform the nuclear burn.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

    CALL radhyd_to_nuclear( nx, i_ray_dim, i_ray, nnc )

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!        Load radiation variables and the state variables before and
!         after the Lagrangian hydro steo, and perform the neutrino
!         energy advection step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

    CALL radhyd_to_nu_e_advct_x( nx, i_ray_dim, i_ray, nez, nnu )

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, i_ray_dim, i_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!              \\\\\ NEUTRINO EQUILIBRATION STEP /////
!
!-----------------------------------------------------------------------

    IF ( nu_equil == 'ye' ) THEN
      CALL radhyd_to_equilibrate_x( nx, nez, nnu, i_ray_dim, i_ray )
    END IF

!-----------------------------------------------------------------------
!
!               \\\\\ NEUTRINO TRANSPORT STEP /////
!
!        Load radiation variables and the state variables and perform
!         the neutrino transport step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

    CALL radhyd_to_transport( nx, i_ray_dim, i_ray, nez, nnu )

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, i_ray_dim, i_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!                    \\\\\ X-REMAP STEP /////
!
!        Load updated variables and the grid before and after the
!         x-Lagrangian step, and perform the x-remap step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

    CALL radhyd_to_remap_x( nx, i_ray, i_ray_dim, nez, nnu, nnc )

    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, i_ray_dim, i_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP 2 OVER THE i_ray_dim RADIAL RAYS |||||
!        |||||       ON EACH PROCESSOR ENDS HERE         |||||
!
!-----------------------------------------------------------------------

  END DO

!-----------------------------------------------------------------------
!
!                       \\\\\ EDIT 1D /////
!
!        Load updated variables and test criteria for a 1D edit.
!        Perform an edit in accordance with the criteria satisfied.
!
!-----------------------------------------------------------------------

  CALL radhyd_to_edit( i_ray_min, i_ray_max, i_ray_dim, i_edit, nx, nez, nnu, nnc )

!-----------------------------------------------------------------------
!
!                       \\\\\ EDIT 2D /////
!
!        Load updated variables and test criteria for a 2D edit.
!        Perform an edit in accordance with the criteria satisfied.
!
!-----------------------------------------------------------------------

  CALL radhyd_to_edit_2D( nx, nez, nnu, i_ray_dim, ny, nnc )

!-----------------------------------------------------------------------
!
!                       \\\\\ EDIT GLOBAL /////
!
!        Load updated variables and test criteria for a Global edit.
!        Perform an edit in accordance with the criteria satisfied.
!
!-----------------------------------------------------------------------

  CALL radhyd_to_edit_Global( nx, nez, nnu, i_ray_dim, ny )

!-----------------------------------------------------------------------
!
!                       \\\\\ EDIT HDF /////
!
!        Load updated variables and test criteria for a MD edit.
!        Perform an edit in accordance with the criteria satisfied.
!
!-----------------------------------------------------------------------

!  CALL radhyd_to_edit_HDF( i_ray_min, i_ray_max, i_ray_dim, nx, ny, nez, nnu, nnc )

!-----------------------------------------------------------------------
!
!                      \\\\\ RESTART /////
!
!        Examine criteria for writing a restart file, and do so if the
!         criteria are met.
!
!-----------------------------------------------------------------------

  CALL radhyd_to_restart( ndim, nez, nnu, i_ray_dim )

!-----------------------------------------------------------------------
!
!                     \\\\\ TERMINATE/////
!
!        Examine criteria for terminating the simulation, and write
!         edit and restart files if the criteria are met.
!
!-----------------------------------------------------------------------

  CALL radhyd_to_terminate( i_ray_min, i_ray_max, i_ray_dim, nx, nez, nnu, &
& nnc, ndim )

!-----------------------------------------------------------------------
!
!            \\\\\ CONTINUE WITH Y-SWEEP IF NDIM > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1 ) THEN

!-----------------------------------------------------------------------
!
!              \\\\\ TRANSPOSE VARIABLES FOR Y-SWEEP /////
!
!-----------------------------------------------------------------------

    CALL radhyd_to_transpose_x_y( nx, i_ray_dim, ny, j_ray_dim, nez, nnu, &
&    nnc, n_proc )

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP OVER THE j_ray_dim ANGULAR RAYS |||||
!        |||||     ON EACH PROCESSOR BEGINS HERE        |||||
!
!-----------------------------------------------------------------------

    sweep                 = 'y'
    DO j_ray = 1,j_ray_dim
      i_radial            = j_ray_dim * myid + j_ray
      j_radial            = i_radial + 1

!-----------------------------------------------------------------------
!
!               \\\\\ UPDATE Y-ARRAY EOS TABLES /////
!
!-----------------------------------------------------------------------

      reset_comp_eos      = 'ye'
      CALL radhyd_to_eos_y_reset( ny, j_ray_dim, j_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!      \\\\\ STORE INITIAL ANGULAR VARIABLES FOR LATER USE /////
!
!-----------------------------------------------------------------------

      CALL store_int_angular_var( j_ray, ny )

!-----------------------------------------------------------------------
!
!               \\\\\ COMPUTE TIME STEP FOR Y-SWEEP /////
!
!-----------------------------------------------------------------------

      CALL radhyd_to_hydro_y_time_step( jmin, jmax, j_ray, j_ray_dim, &
&      dt_y_hydro )

!-----------------------------------------------------------------------
!
!          \\\\\ SUBCYCLE Y-SWEEP IF DTIME_HYDRO < DTNPH /////
!
!-----------------------------------------------------------------------

      dtnph_s             = dtnph
      d_time              = zero
      IF ( dt_y_hydro >= dtnph ) THEN
        i_time_max        = 1
      ELSE
        i_time_max        = 1000
        dtnph             = dt_y_hydro
      END IF

      IF ( sub_cy_yz == 'ye' ) THEN
        i_time_lim        = i_time_max
      ELSE
        i_time_lim        = 1
      END IF

      DO i_time = 1,i_time_lim
        dtnph             = DMIN1( dtnph, dtnph_s - d_time )
        d_time            = d_time + dtnph

!-----------------------------------------------------------------------
!
!                  \\\\\ NEUTRINO Y-STRESS /////
!
!-----------------------------------------------------------------------

        CALL radhyd_to_nu_stress_y( j_ray, j_ray_dim, i_radial, j_radial, &
&        nx, ny, nez, nnu )

!-----------------------------------------------------------------------
!
!                \\\\\ Y-LAGRANGIAN HYDRO STEP /////
!
!        Load EVH1 Variables and perform EVH1 Lagrangian hydro
!         along y-direction.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

        CALL radhyd_to_evh1_y_lagr( ny, j_ray, j_ray_dim, i_radial, &
&        j_radial, nx, i_ray_dim, nz, nez, nnu )

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

        reset_comp_eos    = 'no'
        CALL radhyd_to_eos_y_reset( ny, j_ray_dim, j_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!        Load radiation variables and the state variables before and
!         after the Lagrangian hydro steo, and perform the neutrino
!         energy advection step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

        CALL radhyd_to_nu_e_advct_y( nx, ny, j_ray_dim, j_ray, i_radial, &
&        j_radial, nez, nnu )

!-----------------------------------------------------------------------
!
!                    \\\\\ Y-REMAP STEP /////
!
!        Load updated variables and the grid before and after the
!         x-Lagrangian step, and perform the x-remap step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

        CALL radhyd_to_remap_y( nx, ny, j_ray, j_ray_dim, nez, nnu, nnc )

        reset_comp_eos    = 'no'
        CALL radhyd_to_eos_y_reset( ny, j_ray_dim, j_ray, nnc, reset_comp_eos )

!-----------------------------------------------------------------------
!
!              \\\\\ NEUTRINO EQUILIBRATION STEP /////
!
!-----------------------------------------------------------------------

        IF ( nu_equil == 'ye' ) THEN
          CALL radhyd_to_equilibrate_y( nx, ny, nez, nnu, j_ray_dim, j_ray, &
&          i_radial )
        END IF

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

      END DO

      dtnph               = dtnph_s

!-----------------------------------------------------------------------
!
!        ||||| THE LOOP OVER THE j_ray_dim ANGULAR RAYS |||||
!        |||||       ON EACH PROCESSOR ENDS HERE        |||||
!
!-----------------------------------------------------------------------

    END DO

!-----------------------------------------------------------------------
!
!           \\\\\ TRANSPOSE VARIABLES BACK FOR X-SWEEP /////
!
!-----------------------------------------------------------------------

    CALL radhyd_to_transpose_y_x( nx, i_ray_dim, ny, j_ray_dim, nez, nnu, &
&    nnc, n_proc )

  END IF ! ndim > 1

!-----------------------------------------------------------------------
!
!                  \\\\\ UPDATE EOS TABLES /////
!
!-----------------------------------------------------------------------

  DO i_ray = 1,i_ray_dim
    reset_comp_eos        = 'no'
    CALL radhyd_to_eos_x_reset( nx, i_ray_dim, i_ray, nnc, reset_comp_eos )
  END DO

!-----------------------------------------------------------------------
!
!                   \\\\\ SELECT TIME STEP /////
!
!        Load updated variables and the grid before and after the
!         x-Lagrangian step, and perform the x-remap step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

  CALL radyhd_to_time_step_select( nx, i_ray_dim, ny, ndim )

END DO

END PROGRAM radhyd

