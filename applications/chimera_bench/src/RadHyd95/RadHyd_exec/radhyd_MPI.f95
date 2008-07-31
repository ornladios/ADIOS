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
!    File:         radhyd_MPI
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
!  radhyd_to_edit               : transfer variables and call edit routines
!  radhyd_to_monitor            : transfer variables and call the energy and lepton number conservation check
!  cycle                        : update the cycle number
!  radhyd_x                     : directs the computation of the x-sweep
!  radhyd_y                     : directs the computation of the y-sweep
!  radhyd_z                     : directs the computation of the z-sweep
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
!  cycle_module, edit_module, parallel_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : data_path, log_path, reset_path, nlog
USE parallel_module, ONLY : myid, myid_y, myid_z, ierr, MPI_COMM_ROW, &
& MPI_COMM_COL
USE radial_ray_module, ONLY : jmin, jmax, kmin, kmax, ndim, dtnph, time, &
& sub_cy_yz, nu_equil

USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER(len=128), DIMENSION(3) :: c_path_data     ! character array of data path
CHARACTER(len=128)               :: log_file        ! character array of data path
CHARACTER(len=2)                 :: reset_comp_eos  ! composition EOS reset flag

LOGICAL                          :: l_regrid        ! regrid flag

INTEGER                          :: i               ! x-array extent
INTEGER                          :: nx              ! x-array extent
INTEGER                          :: ny              ! y-array extent
INTEGER                          :: nz              ! z-array extent
INTEGER                          :: nez             ! neutrino energy array extent
INTEGER                          :: nnu             ! neutrino flavor array extent
INTEGER                          :: nnc             ! composition array extent
INTEGER                          :: max_12          ! max(nx,ny,nz)+12
INTEGER                          :: n_proc          ! number of processors assigned to the run
INTEGER                          :: n_proc_y        ! the number of processors assigned to the y-zones
INTEGER                          :: n_proc_z        ! the number of processors assigned to the z-zones
INTEGER                          :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
INTEGER                          :: ik_ray_dim      ! number of z-zones on a processor before swapping with z
INTEGER                          :: j_ray_dim       ! number of radial zones on a processor after swapping with y
INTEGER                          :: k_ray_dim       ! number of radial zones on a processor after swapping with z
INTEGER                          :: j_block_col     ! j block on a given processor
INTEGER                          :: k_block_row     ! k block on a given processor

INTEGER, DIMENSION(20)           :: n_dim_data      ! integer array of array dimensions

INTEGER                          :: ij_ray          ! j-index of a radial ray
INTEGER                          :: ij_ray_min      ! minimum j-index of radial ray
INTEGER                          :: ij_ray_max      ! maximum j-index of radial ray
INTEGER                          :: ik_ray          ! k-index of a radial ray
INTEGER                          :: ik_ray_min      ! minimum k-index of radial ray
INTEGER                          :: ik_ray_max      ! maximum k-index of radial ray
INTEGER                          :: ji_ray          ! x (radial) index of a specific y (angular) ray
INTEGER                          :: jk_ray          ! z (azimuthal) index of a specific y (angular) ray
INTEGER                          :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER                          :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray
INTEGER                          :: i_radial        ! the unshifted radial zone (angular ray) corresponding to ji_ray
INTEGER                          :: j_radial        ! the shifted radial zone (angular ray) corresponding to ji_ray

INTEGER                          :: i_time          ! subcycling index
INTEGER                          :: i_time_max      ! maximum subcycling index for subcycling
INTEGER                          :: i_time_lim      ! maximum subcycling index

INTEGER                          :: i_edit          ! edit parameter
INTEGER                          :: istat           ! open-close file flag
INTEGER                          :: j_ray_bndl      ! polar index of the radial ray bundle on processor myid
INTEGER                          :: k_ray_bndl      ! azimuthal index of the radial ray bundle on processor myid

INTEGER                          :: num_procs       ! number of processors assigned to the run
INTEGER                          :: num_procs_y     ! number of processors assigned to the y index of the rays
INTEGER                          :: num_procs_z     ! number of processors assigned to the z index of the rays

#ifdef ADIOS_MODEL || ADIOS_KEYS
INTEGER                          :: adios_err       ! error flag for ADIOS
#endif

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Array dimensions have been read')
  103 FORMAT (' The number of processors granted=',i6,' is not equal to the number assumed=',i6)
  105 FORMAT (' Array dimensions have been broadcast and unpacked')
  107 FORMAT (' Output data path is ',a128)
  109 FORMAT (' Simulation log path is ',a128)
  111 FORMAT (' Reset keys path is ',a128)
  113 FORMAT (' Array modules have been loaded')
  115 FORMAT (' The number of processors granted for the y block=',i4,' is not equal to the number assumed=',i4)
  117 FORMAT (' The number of processors granted for the z block=',i4,' is not equal to the number assumed=',i4)
  119 FORMAT (' n_proc_y * n_proc_z=',i4,' is not equal to the number of processors assigned=',i4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                \\\\\ INITIALIZATION AND SET UP /////
!
!  Initialize MPI
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

!-----------------------------------------------------------------------
!  Initialize MPI
!-----------------------------------------------------------------------

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid , ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, ierr )

!-----------------------------------------------------------------------
!
!                        \\\\\ ADIOS INIT /////
!
!-----------------------------------------------------------------------
#ifdef ADIOS_MODEL || ADIOS_KEYS
CALL adios_init ('config.xml'//char(0), adios_err)

CALL init_prof('./log'//char(0),20,401,MPI_COMM_WORLD,num_procs,myid) 

! craypat trace
CALL PAT_tracing_state(0, adios_err)

#endif

!-----------------------------------------------------------------------
!                ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Read in and pack array dimensions for simulation if myid = 0
!-----------------------------------------------------------------------

  CALL read_pack_array_dimensions( c_path_data, n_dim_data )
  WRITE (nlog,101)

!-----------------------------------------------------------------------
!  Check the number of processors assigned against those granted
!-----------------------------------------------------------------------

  IF ( n_dim_data(7) /= num_procs ) THEN
    WRITE (nlog,103) num_procs, n_dim_data(7)
    WRITE (*,103) num_procs, n_dim_data(7)
    STOP
  END IF

!-----------------------------------------------------------------------
!                 ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast array dimensions
!-----------------------------------------------------------------------

CALL MPI_BCAST( c_path_data, 3*128, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST( n_dim_data , 20   , MPI_INTEGER  , 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Umpack array dimensions
!-----------------------------------------------------------------------

CALL unpack_array_dimenisons( c_path_data, n_dim_data, nx, ny, nz, nez, &
& nnu, nnc, n_proc, n_proc_y, n_proc_z, ij_ray_dim, ik_ray_dim, j_ray_dim, &
& k_ray_dim, data_path, log_path, reset_path )

!-----------------------------------------------------------------------
!  Ray bundle coordinates
!-----------------------------------------------------------------------

j_ray_bndl           = MOD( myid, n_proc_y ) * ij_ray_dim + 1
k_ray_bndl           = ( myid/n_proc_y ) * ik_ray_dim + 1

!-----------------------------------------------------------------------
!  Set up paths for log and edit files
!-----------------------------------------------------------------------

IF ( log_path(1:1) == ' ' ) THEN
  nlog               = 6
ELSE
  nlog               = 49
  WRITE (log_file,'(a17,i4.4,a1,i4.4,a2)') '/Log_File/sim_log',j_ray_bndl, &
&  '_',k_ray_bndl,'.d'
  log_file           = TRIM(log_path)//TRIM(log_file)
  OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='old')
END IF

!-----------------------------------------------------------------------
!  Set up default path for restart keys file 'reset.d
!-----------------------------------------------------------------------

IF ( reset_path(1:1) == ' ' ) THEN
  WRITE (reset_path,'(a13)') 'Data3/Restart'
END IF ! reset_path(1:1) == ' '

!-----------------------------------------------------------------------
!  Write data paths
!-----------------------------------------------------------------------

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
WRITE (nlog,105)
WRITE (nlog,107) data_path
WRITE (nlog,109) log_path
WRITE (nlog,111) reset_path

!-----------------------------------------------------------------------
!  Load array dimensions into array_module
!-----------------------------------------------------------------------

CALL load_array_module( nx, ny, nz, nez, nnu, nnc, n_proc, n_proc_y, &
& n_proc_z, ij_ray_dim, ik_ray_dim, j_ray_dim, k_ray_dim, max_12 )
IF ( myid == 0 ) WRITE (nlog,113)

!-----------------------------------------------------------------------
!  Define j_block_col and k_block_row, the j and k blocks for each processor
!-----------------------------------------------------------------------

j_block_col          = MOD(myid,n_proc_y)
k_block_row          = myid / n_proc_y

!-----------------------------------------------------------------------
!  Create communicators for each value of k_block_row (for XY transpose)
!-----------------------------------------------------------------------

CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, k_block_row, myid, MPI_COMM_ROW, ierr)
CALL MPI_COMM_RANK(MPI_COMM_ROW, myid_y, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_ROW, num_procs_y, ierr)

!-----------------------------------------------------------------------
!  Create communicators for each value of j_block_col (for XZ transpose)
!-----------------------------------------------------------------------

CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, j_block_col, myid, MPI_COMM_COL, ierr)
CALL MPI_COMM_RANK(MPI_COMM_COL, myid_z, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_COL, num_procs_z, ierr)

!-----------------------------------------------------------------------
!  Check the number of processors assigned against those assumed
!-----------------------------------------------------------------------

  IF ( MOD( num_procs_y, n_proc_y ) /= 0 ) THEN
!  IF ( num_procs_y /= n_proc_y ) THEN
    WRITE (*,115) num_procs_y, n_proc_y
    WRITE (nlog,115) num_procs_y, n_proc_y
    STOP
  END IF

  IF ( MOD( num_procs_z, n_proc_z ) /= 0 ) THEN
!  IF ( num_procs_z /= n_proc_z ) THEN
    WRITE (*,117) num_procs_z, n_proc_z
    WRITE (nlog,117) num_procs_z, n_proc_z
    STOP
  END IF

  IF ( MOD( n_proc_y * n_proc_z, n_proc ) /= 0 ) THEN
!  IF ( n_proc_y * n_proc_z /= n_proc ) THEN
    WRITE (*,119) n_proc_y * n_proc_z, n_proc
    WRITE (nlog,119) n_proc_y * n_proc_z, n_proc
    STOP
  END IF

!-----------------------------------------------------------------------
!  Read in and initialize the simulation
!-----------------------------------------------------------------------

CALL initialize( nx, ny, nz, nez, nnu, nnc, max_12, n_proc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )

i_edit               = 0
ij_ray_min           = 1
ij_ray_max           = ij_ray_dim
ik_ray_min           = 1
ik_ray_max           = ik_ray_dim

!-----------------------------------------------------------------------
!  Initial edit
!-----------------------------------------------------------------------

CALL radhyd_to_edit( ij_ray_min, ij_ray_max, ij_ray_dim, ik_ray_min, &
& ik_ray_max, ik_ray_dim, i_edit, nx, ny, nz, nez, nnu, nnc )

!-----------------------------------------------------------------------
!  Initialize for monitoring energy and lepton conservation
!-----------------------------------------------------------------------

CALL radhyd_to_monitor_MPI( ij_ray_dim, ik_ray_dim, nx, ny, nz, nez, nnu, &
& nnc )

!-----------------------------------------------------------------------
!
!                  \\\\\ MAIN EVOLUTION LOOP /////
!                  /////     BEGINS HERE     \\\\\
!
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------

WRITE (nlog,*) 'Beginning problem loop'

DO

  IF ( log_path(1:1) /= ' ' ) THEN
    CLOSE (unit=nlog,status='keep')
    OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nlog,FILE=TRIM(log_file),STATUS='old',POSITION='append')
  END IF ! log_path(1:1) /= ' '

!-----------------------------------------------------------------------
!  Update cycle number
!-----------------------------------------------------------------------
 
  CALL cycle( time, dtnph )

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_start(ncycle)
#endif

  CALL MPI_Barrier( MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!
!                         \\\\\ X-SWEEP /////
!
!-----------------------------------------------------------------------

  CALL radhyd_x_MPI( nx, ny, nz, nez, nnu, nnc, ndim, ij_ray_min, &
& ij_ray_max, ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, i_edit, &
& time, dtnph, nu_equil )

!-----------------------------------------------------------------------
!
!            \\\\\ CONTINUE WITH Y-SWEEP IF NDIM > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  ny > 1 ) THEN
    CALL radhyd_y_MPI( jmin, jmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, j_ray_dim, n_proc, n_proc_y, n_proc_z, &
&     dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  ny > 1

!-----------------------------------------------------------------------
!
!            \\\\\ CONTINUE WITH Z-SWEEP IF NDIM > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  nz > 1 ) THEN
    CALL radhyd_z_MPI( kmin, kmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, k_ray_dim, n_proc, n_proc_y, n_proc_z, &
&    dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  nz > 1

!-----------------------------------------------------------------------
!
!                     \\\\\ SELECT TIME STEP /////
!
!        Load updated variables and the grid before and after the
!         x-Lagrangian step, and perform the x-remap step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

  CALL radyhd_to_time_step_select_MPI( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, ndim )

!-----------------------------------------------------------------------
!
!                         \\\\\ RESTART /////
!
!        Examine criteria for writing a restart file, and do so if the
!         criteria are met.
!
!-----------------------------------------------------------------------

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL radhyd_to_restart_adios( ndim, nx, nez, nnu )
#else
  CALL radhyd_to_restart( ndim, nx, nez, nnu )
#endif

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_end(ncycle)
#endif

!-----------------------------------------------------------------------
!  Update cycle number
!-----------------------------------------------------------------------

  CALL cycle( time, dtnph )

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_start(ncycle)
#endif

!-----------------------------------------------------------------------
!
!              \\\\\ Z-SWEEP IF NDIM > 1 AND NZ > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  nz > 1 ) THEN
    CALL radhyd_z_MPI( kmin, kmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, k_ray_dim, n_proc, n_proc_y, n_proc_z, &
&    dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  nz > 1

!-----------------------------------------------------------------------
!
!            \\\\\ CONTINUE WITH Y-SWEEP IF NDIM > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  ny > 1 ) THEN
    CALL radhyd_y_MPI( jmin, jmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, j_ray_dim, n_proc, n_proc_y, n_proc_z, &
&    dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  ny > 1

!-----------------------------------------------------------------------
!
!                  \\\\\ CONTINUE WITH X-SWEEP /////
!
!-----------------------------------------------------------------------

  CALL radhyd_x_MPI( nx, ny, nz, nez, nnu, nnc, ndim, ij_ray_min, &
& ij_ray_max, ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, i_edit, &
& time, dtnph, nu_equil )

!-----------------------------------------------------------------------
!
!                     \\\\\ SELECT TIME STEP /////
!
!        Load updated variables and the grid before and after the
!         x-Lagrangian step, and perform the x-remap step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

  CALL radyhd_to_time_step_select_MPI( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, ndim )

!-----------------------------------------------------------------------
!
!                         \\\\\ RESTART /////
!
!        Examine criteria for writing a restart file, and do so if the
!         criteria are met.
!
!-----------------------------------------------------------------------

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL radhyd_to_restart_adios( ndim, nx, nez, nnu )
#else
  CALL radhyd_to_restart( ndim, nx, nez, nnu )
#endif

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_end(ncycle)
#endif

!-----------------------------------------------------------------------
!  Update cycle number
!-----------------------------------------------------------------------

  CALL cycle( time, dtnph )

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_start(ncycle)
#endif

  CALL MPI_Barrier( MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!
!                        \\\\\ X-SWEEP /////
!
!-----------------------------------------------------------------------

  CALL radhyd_x_MPI( nx, ny, nz, nez, nnu, nnc, ndim, ij_ray_min, &
& ij_ray_max, ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, i_edit, &
& time, dtnph, nu_equil )

!-----------------------------------------------------------------------
!
!            \\\\\ CONTINUE WITH Z-SWEEP IF NDIM > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  nz > 1 ) THEN
    CALL radhyd_z_MPI( kmin, kmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, k_ray_dim, n_proc, n_proc_y, n_proc_z, &
&    dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  nz > 1

!-----------------------------------------------------------------------
!
!            \\\\\ CONTINUE WITH Y-SWEEP IF NDIM > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  ny > 1 ) THEN
    CALL radhyd_y_MPI( jmin, jmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, j_ray_dim, n_proc, n_proc_y, n_proc_z, &
&     dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  ny > 1

!-----------------------------------------------------------------------
!
!                     \\\\\ SELECT TIME STEP /////
!
!        Load updated variables and the grid before and after the
!         x-Lagrangian step, and perform the x-remap step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

  CALL radyhd_to_time_step_select_MPI( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, ndim )

!-----------------------------------------------------------------------
!
!                         \\\\\ RESTART /////
!
!        Examine criteria for writing a restart file, and do so if the
!         criteria are met.
!
!-----------------------------------------------------------------------

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL radhyd_to_restart_adios( ndim, nx, nez, nnu )
#else
  CALL radhyd_to_restart( ndim, nx, nez, nnu )
#endif

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_end(ncycle)
#endif

!-----------------------------------------------------------------------
!  Update cycle number
!-----------------------------------------------------------------------

  CALL cycle( time, dtnph )

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_start(ncycle)
#endif

!-----------------------------------------------------------------------
!
!            \\\\\ CONTINUE WITH Y-SWEEP IF NDIM > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  ny > 1 ) THEN
    CALL radhyd_y_MPI( jmin, jmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, j_ray_dim, n_proc, n_proc_y, n_proc_z, &
&    dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  ny > 1

!-----------------------------------------------------------------------
!
!              \\\\\ Z-SWEEP IF NDIM > 1 AND NZ > 1 /////
!
!-----------------------------------------------------------------------

  IF ( ndim > 1  .and.  nz > 1 ) THEN
    CALL radhyd_z_MPI( kmin, kmax, nx, ny, nz, nez, nnu, nnc, ndim, &
&    ij_ray_dim, ik_ray_dim, k_ray_dim, n_proc, n_proc_y, n_proc_z, &
&    dtnph, nu_equil, sub_cy_yz )
  END IF ! ndim > 1  .and.  nz > 1

!-----------------------------------------------------------------------
!
!                  \\\\\ CONTINUE WITH X-SWEEP /////
!
!-----------------------------------------------------------------------

  CALL radhyd_x_MPI( nx, ny, nz, nez, nnu, nnc, ndim, ij_ray_min, &
& ij_ray_max, ij_ray_dim, ik_ray_min, ik_ray_max, ik_ray_dim, i_edit, &
& time, dtnph, nu_equil )

!-----------------------------------------------------------------------
!
!                     \\\\\ SELECT TIME STEP /////
!
!        Load updated variables and the grid before and after the
!         x-Lagrangian step, and perform the x-remap step.
!        Return updated variables to radial_ray_module.
!
!-----------------------------------------------------------------------

  CALL radyhd_to_time_step_select_MPI( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, ndim )

!-----------------------------------------------------------------------
!
!                         \\\\\ RESTART /////
!
!        Examine criteria for writing a restart file, and do so if the
!         criteria are met.
!
!-----------------------------------------------------------------------

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL radhyd_to_restart_adios( ndim, nx, nez, nnu )
#else
  CALL radhyd_to_restart( ndim, nx, nez, nnu )
#endif

#ifdef ADIOS_MODEL || ADIOS_KEYS
  CALL cycle_end(ncycle)
#endif

END DO

!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------

#ifdef ADIOS_MODEL || ADIOS_KEYS
CALL adios_finalize (myid)
#endif

CALL MPI_FINALIZE(ierr)

END PROGRAM radhyd

