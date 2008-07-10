SUBROUTINE initialize( nx, ny, nz, nez, nnu, nnc, max_12, n_proc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         initialize
!    Module:       initialize
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/04
!
!    Purpose:
!      To set up the the problem for execution.
!
!    Subprograms called:
!      dimension_arrays, initialize_variables, problem_read, genst_hy,
!      data_check, rezone, load_radhyd_arrays, load_radhyd_arrays,
!      problem_setup, svel_init
!
!    Input-output arguments:
!  nx         : x-array extent
!  ny         : y-array extent
!  nz         : z-array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  nnc        : composition array extent
!  n_proc     : number of processors assigned to the run
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  j_ray_dim  : number of radial zones on a processor after swapping with y
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, evh1_global, nucbrn_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nprint, data_path, nlog
USE evh1_global, ONLY : small, smallp, smallr
USE nucbrn_module, ONLY : inuc
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input - Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(inout)              :: nx             ! x-array extent
INTEGER, INTENT(inout)              :: ny             ! y-array extent
INTEGER, INTENT(inout)              :: nz             ! z-array extent
INTEGER, INTENT(inout)              :: nez            ! neutrino energy array extent
INTEGER, INTENT(inout)              :: nnu            ! neutrino flavor array extent
INTEGER, INTENT(inout)              :: nnc            ! composition array extent
INTEGER, INTENT(inout)              :: max_12         ! max(nx,ny,nz)+12
INTEGER, INTENT(inout)              :: n_proc         ! number of processors assigned to the run
INTEGER, INTENT(inout)              :: ij_ray_dim     ! number of y-zones on a processor before swapping
INTEGER, INTENT(inout)              :: ik_ray_dim     ! number of z-zones on a processor before swapping
INTEGER, INTENT(inout)              :: j_ray_dim      ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(inout)              :: k_ray_dim      ! number of radial zones on a processor after swapping with z

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)                 :: rst_mod_file   ! character string containing name of restart_model file
CHARACTER (len=128)                 :: rst_tmp_file1  ! character string containing name of temporary restart file
CHARACTER (len=128)                 :: rst_tmp_file2  ! character string containing name of temporary restart file
CHARACTER (len=128)                 :: rst_fnl_file   ! character string containing name of final restart file

CHARACTER(len=128), DIMENSION(1)    :: c_init_data    ! character array of initial data
CHARACTER (len=2), DIMENSION(20)    :: c_radhyd_data  ! character array of radhyd keys
CHARACTER (len=1), DIMENSION(1)     :: c_eos_data     ! character array of edit keys
CHARACTER (len=5), DIMENSION(nnc)   :: c_nuc_data     ! character array of nuclei

LOGICAL, DIMENSION(1)               :: l_rezone_data  ! lagrangian flag
LOGICAL                             :: l_multi_d      ! multiD flag

INTEGER, DIMENSION(2)               :: i_init_data    ! integer array of initial data
INTEGER, DIMENSION(50)              :: i_radhyd_data  ! integer array of radhyd keys
INTEGER, DIMENSION(40+2*nnu)        :: i_trans_data   ! integer array of transport keys
INTEGER, DIMENSION(5)               :: i_e_advct_data ! integer array of e_advect keys
INTEGER, DIMENSION(1200+3*40*nnu)   :: i_edit_data    ! integer array of edit keys
INTEGER, DIMENSION(30)              :: i_hydro_data   ! integer array of transport keys
INTEGER, DIMENSION(10+nx,ij_ray_dim,ik_ray_dim) :: i_nuc_data ! integer array of edit keys
INTEGER, DIMENSION(2)               :: i_model_data   ! integer array of initial model data
INTEGER                             :: nrst           ! cycle number to start simulation
INTEGER                             :: nouttmp        ! reatart read flag
INTEGER                             :: istat          ! open file flag
INTEGER                             :: iskip          ! print flag for echoing the reading of data
INTEGER                             :: nuc_number     ! number of nuclear species (not counting representative heavy nucleus)

INTEGER, PARAMETER                                          :: n_restart = 23

REAL(KIND=double), DIMENSION(50)                            :: d_radhyd_data  ! 64 bit real array of radhyd keys
REAL(KIND=double), DIMENSION(14)                            :: d_eos_data     ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION((110+3*nnu+3*nez+1))           :: d_trans_data   ! 64 bit real array of transport keys
REAL(KIND=double), DIMENSION(5+2*nnu)                       :: d_e_advct_data ! 64 bit real array of e_advect keys
REAL(KIND=double), DIMENSION(50)                            :: d_edit_data    ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION((30+nx))                       :: d_hydro_data   ! 64 bit real array of transport keys
REAL(KIND=double), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data     ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)         :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), DIMENSION(2,nx)                               :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), DIMENSION(2,nx+1)                             :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), DIMENSION(8,nez,nnu,ij_ray_dim,ik_ray_dim)    :: d_psi_data1   ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: d_psi_data2   ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION(2,nx,nnu,ij_ray_dim,ik_ray_dim)     :: d_psi_data3   ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION(2,nnu,ij_ray_dim,ik_ray_dim)        :: d_psi_data4   ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION(2,ij_ray_dim,ik_ray_dim)            :: d_psi_data5   ! 64 bit real array of edit keys
REAL(KIND=double), DIMENSION(20+3*ny+3*nz+2)                :: d_rezone_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (' Need MPI version of initialize since n_proc=',i4,' > 1')
  101 FORMAT (' Arrays are now being dimensioned on processor',i6)
  103 FORMAT (' All arrays have now been dimensioned on processor',i6)
  105 FORMAT (' Variables are being initialized on processor',i6)
  107 FORMAT (' All variables have now been initialized on processor',i6)
  109 FORMAT (' Problem data is now being read on processor',i6)
  111 FORMAT (' Problem data has now been read on processor',i6)
  113 FORMAT (' Data are now checked for consistency on processor',i6)
  115 FORMAT (' Data have been checked for consistency on processor',i6)
  117 FORMAT (' y- and x-zoning is now being performed, and the model rezoned, if needed on processor',i6)
  119 FORMAT (' y- and x-zoning has been performed, and the model rezoned, if needed on processor',i6)
  121 FORMAT (' Arrays are being distributed, unpacked, and loaded on processor',i6)
  123 FORMAT (' Arrays have been distributed, unpacked, and loaded on processor',i6)
  125 FORMAT (' Model extents are now checked for consistency on processor',i6)
  127 FORMAT (' Model extents have been checked for consistency on processor',i6)
  129 FORMAT (' NSE array is being set on processor',i6)
  131 FORMAT (' NSE array has been set on processor',i6)
  133 FORMAT (' Electron fraction is being adjusted to be consistent with the composition on processor',i6)
  135 FORMAT (' Electron fraction has bee adjusted to be consistent with the composition on processor',i6)
  137 FORMAT (' Solid angles are being computed on all processors on processor',i6)
  139 FORMAT (' Solid angles have been computed on all processors on processor',i6)
  141 FORMAT (' Problem is being initialized on all processors on processor',i6)
  143 FORMAT (' Problem has been initialized on all processors on processor',i6)
  145 FORMAT (' Nuclear reaction rate and partition function data will be read in and broadcast on processor',i6)
  147 FORMAT (' Nuclear reaction rate and partition function data have been read in and broadcast on processor',i6)
  149 FORMAT (' Nuclear reactions are turned off on processor',i6)
  151 FORMAT (' EVH1 arrays are loaded on all processors on processor',i6)
  153 FORMAT (' EVH1 arrays have been loaded on all processors on processor',i6)
  155 FORMAT (' Time step check is being performed on all processors on processor',i6)
  157 FORMAT (' Time step check has been performed on all processors on processor',i6)
 1001 FORMAT (' nouttmp=',i3,' in initialize, which is not in the range of 1 - 4; nrst=',i8, &
 &              ' ny=',i4,' nz=',i4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Stop if n_proc > 1.
!-----------------------------------------------------------------------

IF ( n_proc > 1 ) THEN
  WRITE (nprint,1) n_proc
  WRITE (nlog,1) n_proc
  STOP
END IF ! n_proc > 1

!-----------------------------------------------------------------------
!  Set some floor values for pressure and density
!-----------------------------------------------------------------------

small                = 1.d-15 
smallp               = 1.d-15 
smallr               = 1.d-15 

!-----------------------------------------------------------------------
!
!                 \\\\\ DIMENSION ARRAYS /////
!
!  dimension_arrays
!  Calls
!   dimension_radhyd_arrays
!   dimension_hydro_arrays
!   dimension_tov_pot_arrays
!   dimension_mgfld_arrays
!   dimension_edit_arrays
!   dimension_eos_arrays
!   dimension_eos_drv_arrays
!   dimension_nucbrn_arrays
!   dimension_e_advct_arrays
!-----------------------------------------------------------------------

WRITE (nlog,101) myid
CALL dimension_arrays( nx, ny, nz, nez, nnu, nnc, max_12, n_proc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )
WRITE (nlog,103) myid

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE VARIABLES /////
!
!  initialize_variables
!  Calls
!   initialize_global_var
!   initialize_cycle_arrays
!   initialize_it_tol_arrays
!   initialize_bomb_arrays
!   initialize_rezone_arrays
!-----------------------------------------------------------------------

WRITE (nlog,105) myid
CALL initialize_variables
WRITE (nlog,107) myid

!-----------------------------------------------------------------------
!
!                        \\\\\ DATA READ /////
!
!  Open files reset.d, superdump.d, rstdmp1.d, rstdmp2.d, and read in
!   initial data or restart data. Pack data for broadcasting to the
!   assigned processors.
!
!  Calls problem_read
!        ------------
!
!    Calls read_pack                      : orchestrates the read-in of the problem
!          ---------
!
!      Calls read_pack_init               : reads header, cycle number (nrst), and flag nouttmp
!
!      IF nrst = 0 (initializing a simulation)
!      Calls radhyd_read                  : opens and closes radhyd_keys.d
!        Calls read_pack_radhyd_keys      : reads and packs the data in radhyd_keys.d
!      Calls  eos_read                    : opens and closes eos_keys.d
!        Calls read_pack_eos_keys         : reads and packs the data in eos_keys.d
!      Calls  transport_read              : opens and closes transport_keys.d
!        Calls read_pack_transport_keys   : reads and packs the data in transport_keys.d
!      Calls  e_advct_read                : opens and closes e_advct_keys.d
!        Calls read_pack_e_advct_keys     : reads and packs the data in e_advct_keys.d
!      Calls  edit_read                   : opens and closes edit_keys.d
!        Calls read_pack_edit_keys        : reads and packs the data in edit_keys.d
!      Calls  hydro_read                  : opens and closes hydro_keys.d
!        Calls read_pack_hydro_keys       : reads and packs the data in hydro_keys.d
!      Calls  nuclear_read                : opens and closes nuclear_keys.d
!        Calls read_pack_nuclear_keys     : reads and packs the data in nuclear_keys.d
!      Calls  model_read                  : opens and closes initial_model.d
!        Calls read_pack_initial_model    : reads and packs the data in initial_model.d
!
!      IF nrst > 0 (restarting a simulation)
!        Opens restart files depending on the value of the nouttmp key
!        Calls read_pack_radhyd_keys      : reads and packs the radhyd_keys data from common keys file
!        Calls read_pack_eos_keys         : reads and packs the eos_keys data from common keys file
!        Calls read_pack_transport_keys   : reads and packs the transport_keys data from common keys file
!        Calls read_pack_e_advct_keys     : reads and packs the e_advct_keys data from common keys file
!        Calls read_pack_edit_keys        : reads and packs the edit_keys data from common keys file
!        Calls read_pack_hydro_keys       : reads and packs the hydro_keys data from common keys file
!        Calls read_pack_nuclear_keys     : reads and packs the muclear_keys data from common keys file
!
!      IF ny = 1 and nz = 1 (1-dimensional simultion)
!        Opens restart model file depending on the value of the nouttmp key
!        Calls read_pack_restart_model    : reads in the model configuration
!
!      IF ny > 1 or nz > 1 (multi-dimensional simultion)
!        Opens restart model file depending on the value of the nouttmp key
!        Calls read_pack_restart_model    : reads in the model configuration for selected data
!
!       Saves initial file reset.d to reset_initial.d
!
!  Calls data_check                       : checks keys, array values, processor
!                                            requests, etc. for consistency
!  Calls rezone                           : sets up y and z grids
!
!  Calls unpack_arrays                    : Orchestrates the unpacking the receive buffers
!        -------------                       and the loading of modules on each node
!
!    Calls unpack_init                    : unpacks the header, cycle number, and nouttmp key
!    Calls unpack_radial_ray_keys         : unpacks the radhyd ray keys
!    Calls unpack_rezone_arrays           : unpacks cpprdinae and grid data
!    Calls unpack_eos_keys                : unpacks the equation of state keys
!    Calls unpack_transport_keys          : unpacks the transport keys
!    Calls unpack_e_advct_keys            : unpacks the neutrino energy advection keys
!    Calls upack_edit_keys                : unpacks the edit keys
!    Calls unpack_hydro_keys              : unpacks the hydro keys
!    Calls unpack_nuclear_data            : unpacks the nuclear keys
!    Calls unpack_initial_model           : unpacks the initial model, or parts of the reatart model
!    Calls unpack_init                    : unpacks the header, cycle number, and nouttmp key
!
!  IF ( nrst /= 0 and  ny > 1 or nz > 1 ) and nouttmp /= 5 ) THEN
!    Opens restart model file depending on the value of the nouttmp key
!    Calls read_restart_model_to_node     : reads in model data to each node
!
!  Calls ye_fix                           : Ensures that the value of ye and the ye given
!                                            by the nuclei are consistent
!
!  Calls problem_setup                    : Orchestrates the initiation or restart of the simulation
!        -------------
!    Calls setup_initial_hydro            : sets up all the hydro and EOS variables (nrst = 0)
!    Calls setup_restart_hydro            : sets up all the hydro and EOS variables (nrst > 0)
!    Calls solid_angle                    : computes the solid angle subtended by each ray
!    Calls angular_ave                    : computes the angular averages of selected quantities
!    Calls radhyd_to_poisson              : computes the gravitational forces
!    Calls setup_initial_GR               : initializes the GR variables (nrst = 0)
!    Calls setup_restart_GR               : initializes the GR variables (nrst > 0)
!    Calls setup_initial_trans            : initializes neutrino opacities and transport quantities (nrst = 0)
!    Calls setup_restart_trans            : initializes neutrino opacities and transport quantities (nrst > 0)
!    Calls load_radial_ray_arrays         : loads initialized quantities into RADIAL_RAY_MODULE
!
!  Calls network_setup                    : initializes the nuclear reaction rates and quantities
!
!  Calls nse_set                          : sets nse flag for all x-array zones outside of model to the outermost
!                                            value of the flag in the model domain
!
!  Calls load_evh1_arrays                 : loads keys and initialized quantities into the hydro modules
!
!  Calls time_step_check                  : Warns if the input time step exceeds the Courant time step
!
!  Calls solid_angle                      : Computes once more the solid angle subtended by each ray
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize buffers
!-----------------------------------------------------------------------

i_radhyd_data      = 0
i_trans_data       = 0
i_e_advct_data     = 0
i_edit_data        = 0
i_hydro_data       = 0
i_nuc_data         = 0
i_model_data       = 0

d_radhyd_data      = zero
d_eos_data         = zero
d_trans_data       = zero
d_e_advct_data     = zero
d_edit_data        = zero
d_hydro_data       = zero
d_nuc_data         = zero
d_model_data1      = zero
d_model_data2      = zero
d_model_data3      = zero
d_psi_data1        = zero
d_psi_data2        = zero
d_psi_data3        = zero
d_psi_data4        = zero
d_psi_data5        = zero
d_rezone_data      = zero

!-----------------------------------------------------------------------
!  Input data
!-----------------------------------------------------------------------

WRITE (nlog,109) myid
  CALL problem_read( c_init_data, c_radhyd_data, c_eos_data, c_nuc_data, &
& i_init_data, i_radhyd_data, i_trans_data, i_e_advct_data, i_edit_data, &
& i_hydro_data, i_nuc_data, i_model_data, d_radhyd_data, d_eos_data,     &
& d_trans_data, d_e_advct_data, d_edit_data, d_hydro_data, d_nuc_data,   &
& d_model_data1, d_model_data2, d_model_data3, d_psi_data1, d_psi_data2, &
& d_psi_data3, d_psi_data4, d_psi_data5, nrst, nouttmp )
WRITE (nlog,111) myid

!-----------------------------------------------------------------------
!
!                \\\\\ DATA CONSISTENCY CHECK /////
!
!-----------------------------------------------------------------------

WRITE (nlog,113) myid
CALL data_check( c_radhyd_data, i_radhyd_data, d_radhyd_data, &
& d_hydro_data, nx, ny, nz, nnu, n_proc )
WRITE (nlog,115) myid

!-----------------------------------------------------------------------
!
!             \\\\\ Y-, Z- GRIDDING, X-REGRIDDING /////
!
!  rezone
!  Calls
!   lagregrid
!   eulregrid
!-----------------------------------------------------------------------

WRITE (nlog,117) myid
CALL rezone( c_radhyd_data, i_radhyd_data, d_radhyd_data,             &
& l_rezone_data, d_rezone_data )
WRITE (nlog,119) myid

!-----------------------------------------------------------------------
!
!                     \\\\\ UNPACK ARRAYS /////
!
!  unpack_arrays
!  Calls
!   unpack_init
!   unpack_radial_ray_keys
!   unpack_eos_keys
!   upack_transport_keys
!   unpack_e_advct_keys
!   upack_edit_keys
!   unpack_hydro_keys
!   unpack_nuclear_data
!   unpack_initial_model
!-----------------------------------------------------------------------

l_multi_d          = .false.
WRITE (nlog,121) myid
CALL unpack_arrays( c_init_data, c_radhyd_data, c_eos_data, c_nuc_data,  &
& i_init_data, i_radhyd_data, i_trans_data, i_e_advct_data, i_edit_data, &
& i_hydro_data, i_nuc_data, i_model_data, l_rezone_data, d_radhyd_data,  &
& d_eos_data, d_trans_data, d_e_advct_data, d_edit_data, d_hydro_data,   &
& d_nuc_data, d_model_data1, d_model_data2, d_model_data3, d_psi_data1,  &
& d_psi_data2, d_psi_data3, d_psi_data4, d_psi_data5, d_rezone_data,     &
& l_multi_d )
WRITE (nlog,123) myid

!-----------------------------------------------------------------------
!
!            \\\\\ READ MODEL DATA TO EACH NODE /////
!            \\\\\  IF NDIM > 1 AND NRST > 0    /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Unpack cycle number and nouttmp key
!-----------------------------------------------------------------------

nrst               = i_init_data(1)
nouttmp            = i_init_data(2)

!-----------------------------------------------------------------------
!  Open files
!-----------------------------------------------------------------------

IF ( nrst /= 0  .and.  ( ny > 1  .or.  nz > 1 ) ) THEN

  IF ( nouttmp == 1 ) THEN
    WRITE (rst_tmp_file1,'(a18,i4.4,a2)') '/Restart/rst_tmp1_', &
&    myid,'.d'
    rst_tmp_file1 = TRIM(data_path)//TRIM(rst_tmp_file1)
    OPEN (UNIT=n_restart,FILE=TRIM(rst_tmp_file1), STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(rst_tmp_file1), &
&    FORM='unformatted', STATUS='old')

  ELSE IF ( nouttmp == 2 ) THEN
    WRITE (rst_tmp_file2,'(a18,i4.4,a2)') '/Restart/rst_tmp2_', &
&    myid,'.d'
    rst_tmp_file2 = TRIM(data_path)//TRIM(rst_tmp_file2)
    OPEN (UNIT=n_restart,FILE=TRIM(rst_tmp_file2),STATUS='new', &
&    FORM='unformatted',  IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(rst_tmp_file2), &
&    FORM='unformatted', STATUS='old')
    
  ELSE IF ( nouttmp == 3 ) THEN
    WRITE (rst_mod_file,'(a22,i4.4,a1,i7.7,a2)') '/Restart/restart_model', &
&    myid,'_',nrst,'.d'
    rst_mod_file  = TRIM(data_path)//TRIM(rst_mod_file)
    OPEN (UNIT=n_restart,FILE=TRIM(rst_mod_file),STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(rst_mod_file), &
&    FORM='unformatted', STATUS='old')
    
  ELSE IF ( nouttmp == 4 ) THEN
    WRITE (rst_fnl_file,'(a26,i4.4,a2)') '/Restart/restart_final_mod', &
&    myid,'.d'
    rst_fnl_file  = TRIM(data_path)//TRIM(rst_fnl_file)
    OPEN (UNIT=n_restart,FILE=TRIM(rst_fnl_file),STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(rst_fnl_file), &
&    FORM='unformatted', STATUS='old')

  ELSE
    WRITE (nprint,1001) nouttmp, nrst, ny, nz
    WRITE (nlog,1001) nouttmp, nrst, ny, nz
    STOP

  END IF ! nouttmp == 1

!-----------------------------------------------------------------------
!  Commemce reading
!-----------------------------------------------------------------------

  iskip             = 0
  nuc_number        = i_nuc_data(3,1,1)

  REWIND n_restart

  CALL read_restart_model_to_node( n_restart, nprint, iskip, nx, nez, nnu, &
&  nnc, nrst, ij_ray_dim, ik_ray_dim, i_nuc_data, nuc_number )

!-----------------------------------------------------------------------
!  Close restart file
!-----------------------------------------------------------------------

  CLOSE (UNIT=n_restart,STATUS='keep')

END IF ! nrst /= 0  .and.  ny > 1

!-----------------------------------------------------------------------
!  Check consistency of declared extents and model extents.
!-----------------------------------------------------------------------

WRITE (nlog,125) myid
CALL model_check( i_radhyd_data, i_model_data )
WRITE (nlog,127) myid

!-----------------------------------------------------------------------
!  Set nse flag for all x-array zones outside of model.
!-----------------------------------------------------------------------

WRITE (nlog,129) myid
CALL nse_set( nx, ij_ray_dim, ik_ray_dim )
WRITE (nlog,131) myid

!-----------------------------------------------------------------------
!  Set the value of ye for nonNSE material consistent with the
!   composition
!-----------------------------------------------------------------------

WRITE (nlog,133) myid
CALL ye_fix( nx, ij_ray_dim, ik_ray_dim )
WRITE (nlog,135) myid

!-----------------------------------------------------------------------
!  Compute the solid angle subtended by each radial ray.
!-----------------------------------------------------------------------

WRITE (nlog,137) myid
CALL solid_angle
WRITE (nlog,139) myid

!-----------------------------------------------------------------------
!              \\\\\ FINISH SETTING UP THE PROBLEM /////
!
!  problem_setup
!
!  Calls
!   setup_initial_hydro
!   setup_restart_hydro
!   angular_ave
!   radhyd_to_poisson
!   setup_initial_GR
!   setup_restart_GR
!   setup_initial_trans
!   setup_restart_trans
!   load_radial_ray_arrays
!-----------------------------------------------------------------------

WRITE (nlog,141) myid
CALL problem_setup( nx, ny, nz, nez, nnu, ij_ray_dim, ik_ray_dim, nnc )
WRITE (nlog,143) myid

!-----------------------------------------------------------------------
!  Read in and set nuclear reaction rate data if inuc = 1.
!-----------------------------------------------------------------------

WRITE (nlog,145) myid
IF ( inuc == 1 ) THEN
  CALL network_setup
  WRITE (nlog,147) myid
ELSE
  WRITE (nlog,149) myid
END IF

!-----------------------------------------------------------------------
!  Load evh1 arrays.
!-----------------------------------------------------------------------

WRITE (nlog,151) myid
CALL load_evh1_arrays
WRITE (nlog,153) myid

!-----------------------------------------------------------------------
!  Compute sound crossing time.
!-----------------------------------------------------------------------

WRITE (nlog,155) myid
CALL time_step_check( ij_ray_dim, ik_ray_dim )
WRITE (nlog,157) myid

RETURN
END SUBROUTINE initialize
