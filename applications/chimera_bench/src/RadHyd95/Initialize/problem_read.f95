SUBROUTINE problem_read( c_init_data, c_radhyd_data, c_eos_data,        &
& c_nuc_data, i_init_data, i_radhyd_data, i_trans_data, i_e_advct_data, &
& i_edit_data, i_hydro_data, i_nuc_data, i_model_data, d_radhyd_data,   &
& d_eos_data, d_trans_data, d_e_advct_data, d_edit_data, d_hydro_data,  &
& d_nuc_data, d_model_data1, d_model_data2, d_model_data3, d_psi_data1, &
& d_psi_data2, d_psi_data3, d_psi_data4, d_psi_data5, nrst, nouttmp )
!-----------------------------------------------------------------------
!
!    File:         problem_read
!    Module:       problem_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  K. R. DeNisco, Dept. of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/10/02
!
!    Purpose:
!      To open files reset.d, superdump.d, rstdmp1.d, rstdmp2.d, and
!       to read in initial value or restart data.
!
!
!    Subprograms called:
!  read_pack      : directs the reading and packing of the initial or restart data
!
!    Input arguments:
!        none
!
!    Output arguments:
!  c_init_data    : character array of initial data
!  c_radhyd_data  : character array of radhyd keys
!  c_eos_data     : character array of edit keys
!  c_nuc_data     : character array of nuclei
!  i_init_data    : integer array of initial data
!  i_radhyd_data  : integer array of radhyd keys
!  i_trans_data   : integer array of transport keys
!  i_e_advct_data : integer array of e_advect keys
!  i_edit_data    : integer array of transport keys
!  i_hydro_data   : integer array of transport keys
!  i_nuc_data     : integer array of edit keys
!  i_model_data   : integer array of initial model data
!  d_radhyd_data  : 64 bit real array of radhyd keys
!  d_eos_data     : 64 bit real array of edit keys
!  d_trans_data   : 64 bit real array of transport keys
!  d_e_advct_data : 64 bit real array of e_advect keys
!  d_edit_data    : 64 bit real array of edit keys
!  d_hydro_data   : 64 bit real array of hydro keys
!  d_nuc_data     : 64 bit real array of nuclear keys
!  d_model_data1  : 64 bit real array of initial model data
!  d_model_data2  : 64 bit real array of initial model data
!  d_model_data3  : 64 bit real array of initial model data
!  d_psi_data1    : 64 bit real array of neutrino data
!  d_psi_data2    : 64 bit real array of neutrino data
!  d_psi_data3    : 64 bit real array of neutrino data
!  d_psi_data4    : 64 bit real array of neutrino data
!  d_psi_data5    : 64 bit real array of neutrino data
!  nrst           : cycle number to start simulation
!  nouttmp        : restart read flag
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nez, nnu, nnc, ij_ray_dim, ik_ray_dim
USE numerical_module, ONLY : zero

USE edit_module, ONLY: nread, nprint, nlog, nrstd1, nrstd2, nrrst, data_path

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER(len=128), INTENT(out), DIMENSION(1)    :: c_init_data    ! character array of initial data
CHARACTER (len=2), INTENT(out), DIMENSION(20)    :: c_radhyd_data  ! character array of radhyd keys
CHARACTER (len=1), INTENT(out), DIMENSION(1)     :: c_eos_data     ! character array of edit keys
CHARACTER (len=5), INTENT(out), DIMENSION(nnc)   :: c_nuc_data     ! character array of nuclei

INTEGER, INTENT(out), DIMENSION(2)               :: i_init_data    ! integer array of initial data
INTEGER, INTENT(out), DIMENSION(50)              :: i_radhyd_data  ! integer array of radhyd keys
INTEGER, INTENT(out), DIMENSION(40+2*nnu)        :: i_trans_data   ! integer array of transport keys
INTEGER, INTENT(out), DIMENSION(5)               :: i_e_advct_data ! integer array of e_advect keys
INTEGER, INTENT(out), DIMENSION(1200+3*40*nnu)   :: i_edit_data    ! integer array of edit keys
INTEGER, INTENT(out), DIMENSION(30)              :: i_hydro_data   ! integer array of transport keys
INTEGER, INTENT(out), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim)   :: i_nuc_data     ! integer array of edit keys
INTEGER, INTENT(out), DIMENSION(2)               :: i_model_data   ! integer array of initial model data
INTEGER, INTENT(out)                             :: nrst           ! cycle number to start simulation
INTEGER, INTENT(out)                             :: nouttmp        ! reatart read flag

REAL(KIND=double), INTENT(out), DIMENSION(50)                            :: d_radhyd_data  ! 64 bit real array of radhyd keys
REAL(KIND=double), INTENT(out), DIMENSION(14)                            :: d_eos_data     ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION((110+3*nnu+3*nez+1))           :: d_trans_data   ! 64 bit real array of transport keys
REAL(KIND=double), INTENT(out), DIMENSION(5+2*nnu)                       :: d_e_advct_data ! 64 bit real array of e_advect keys
REAL(KIND=double), INTENT(out), DIMENSION(50)                            :: d_edit_data    ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION((30+nx))                       :: d_hydro_data   ! 64 bit real array of transport keys
REAL(KIND=double), INTENT(out), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data     ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)                :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx)                                      :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx+1)                                    :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(8,nez,nnu,ij_ray_dim,ik_ray_dim)           :: d_psi_data1   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim)        :: d_psi_data2   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,nx,nnu,ij_ray_dim,ik_ray_dim)            :: d_psi_data3   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,nnu,ij_ray_dim,ik_ray_dim)               :: d_psi_data4   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,ij_ray_dim,ik_ray_dim)                   :: d_psi_data5   ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                                        :: istat          ! open-close file flag

   11 FORMAT (a128)
   21 FORMAT (10x,i10)
   23 FORMAT ()
   25 FORMAT ('!-----------------------------------------------------------------------')
   27 FORMAT ('!')
   29 FORMAT ('! Set nrst to 0 to read initial model')
   31 FORMAT ('! If nrst /= 0, set noutput to 1 to read in rst_tmp1_ files')
   33 FORMAT ('! If nrst /= 0, set noutput to 2 to read in rst_tmp2_ files')
   35 FORMAT ('! If nrst /= 0, set noutput to 3 to read in restart_model files at cycle number nrst')
   37 FORMAT ('! If nrst /= 0, set noutput to 4 to read in restart_final_mod files')
   39 FORMAT ('! If nrst /= 0, set noutput to 5 to read 1-D restart_model files at cycle number nrst and begin MD run')
 1001 FORMAT (' Error in closing reset.d in subroutine problem_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                  ||||| MYID == 0 THROUGHOUT |||||
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

c_init_data        = ' '
c_eos_data         = ' '
c_nuc_data         = ' '

i_init_data        = 0
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

!-----------------------------------------------------------------------
!  Set unit numbers for nread and nprint
!-----------------------------------------------------------------------

nread              = 11
nprint             = 41

!-----------------------------------------------------------------------
!
!             \\\\\ OPEN FILES RESET.D AMD SUPERDUMP.D /////
!
!        File reset.d locates the initialization or restart files.
!
!        File superdump.d is a scratch file for writing dagnostics.
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/reset.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/reset.d',STATUS='old')
OPEN (UNIT=nprint,FILE=TRIM(data_path)//'/Run_Log/superdump.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nprint,FILE=TRIM(data_path)//'/Run_Log/superdump.d',STATUS='old')

!-----------------------------------------------------------------------
!  Set unit numbers for nrrst, nrstd1, and nrstd2
!-----------------------------------------------------------------------

nrrst              = nread
nrstd1             = 19
nrstd2             = 21

!-----------------------------------------------------------------------
!
!            \\\\\ READ INITIAL MODEOL AND PROBLEM KEYS /////
!
!         read_pack
!
!          If nrst = 0 Calls
!           read_init
!           radhyd_read
!            Calls read_pack_radhyd_keys
!           eos_read
!            Calls read_pack_eos_keys
!           transport_read
!            Calls read_pack_transport_keys
!           e_advct_read
!            Calls read_pack_e_advct_keys
!           edit_read
!            Calls read_pack_edit_keys
!           hydro_read
!            Calls read_pack_hydro_keys
!           nuclear_read
!            Calls read_pack_initial_nuclear_data
!           model_read
!            Calls read_pack_initial_model
!
!          If nrst > 0 Calls
!            read_pack_radhyd_keys
!            read_pack_eos_keys
!            read_pack_transport_keys
!            read_pack_edit_keys
!            read_pack_hydro_keys
!            read_pack_nuclear_keys
!            read_pack_initial_model
!
!           If ny = 1 and nz = 0
!             read_pack_restart_model
!           IF ny > 1 or nz > 0
!             read_pack_restart_model
!             read_restart_model_to_node
!-----------------------------------------------------------------------

CALL read_pack( c_init_data, c_radhyd_data, c_eos_data, c_nuc_data,      &
& i_init_data, i_radhyd_data, i_trans_data, i_e_advct_data, i_edit_data, &
& i_hydro_data, i_nuc_data, i_model_data, d_radhyd_data, d_eos_data,     &
& d_trans_data, d_e_advct_data, d_edit_data, d_hydro_data, d_nuc_data,   &
& d_model_data1, d_model_data2, d_model_data3, d_psi_data1, d_psi_data2, &
& d_psi_data3, d_psi_data4, d_psi_data5, nrst, nouttmp )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
IF ( istat /= 0 ) WRITE (nlog,1001)

!-----------------------------------------------------------------------
!  Save initial file reset.d
!-----------------------------------------------------------------------

IF ( nrst == 0 ) THEN
  OPEN (UNIT=nread,FILE='Data3/Initial_Data/reset_initial.d',STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/reset_initial.d',STATUS='old')
  WRITE (nread,11) c_init_data
  WRITE (nread,21) nrst
  WRITE (nread,23) 
  WRITE (nread,23) 
  WRITE (nread,23) 
  WRITE (nread,25) 
  WRITE (nread,27) 
  WRITE (nread,29) 
  WRITE (nread,31) 
  WRITE (nread,33) 
  WRITE (nread,35) 
  WRITE (nread,37) 
  WRITE (nread,39) 
  WRITE (nread,27) 
  WRITE (nread,25) 
  CLOSE (UNIT=nread,STATUS='keep')
END IF
      
RETURN
END SUBROUTINE problem_read

