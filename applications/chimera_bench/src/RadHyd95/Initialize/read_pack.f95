SUBROUTINE read_pack( c_init_data, c_radhyd_data, c_eos_data, c_nuc_data, &
& i_init_data, i_radhyd_data, i_trans_data,i_e_advct_data, i_edit_data,   &
& i_hydro_data, i_nuc_data, i_model_data, d_radhyd_data, d_eos_data,      &
& d_trans_data, d_e_advct_data, d_edit_data, d_hydro_data, d_nuc_data,    &
& d_model_data1, d_model_data2, d_model_data3, d_psi_data1, d_psi_data2,  &
& d_psi_data3, d_psi_data4, d_psi_data5, nrst, nouttmp )
!-----------------------------------------------------------------------
!
!    File:         read_pack
!    Module:       read_pack
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To direct the read-in of the problem configuration.
!
!    Subprograms called:
!  read_pack_init           : reads in nrst and noutput from file reset.d
!  radhyd_read              : reads in basic simulation keys from file radhyd_keys.d
!  eos_read                 : reads in equation of state keys from file eos_keys.d
!  transport_read           : reads in transport keys from file transport_keys.d
!  e_advct_read             : reads in energy advection keys from file e_advct_keys.d
!  edit_read                : reads in edit keys from file edit_keys.d
!  hydro_read               : reads in hydro keys from file hydro_keys.d
!  nuclear_read             : reads in nuclear keys from file nuclear_keys.d
!  model_read               : reads in model configuration from file initial_model.d
!  read_pack_radhyd_keys    : reads and packs basic simulation keys from a restart file
!  read_pack_eos_keys       : reads and packs equation of state keys from a restart file
!  read_pack_transport_keys : reads and packs transport keys from a restart file
!  read_pack_e_advct_keys   : reads and packs energy advection keys from a restart file
!  read_pack_edit_keys      : reads and packs edit keys from a restart file
!  read_pack_hydro_keys     : reads and packs hydro keys from a restart file
!  read_pack_nuclear_keys   : reads and packs nuclear keys from a restart file
!  read_pack_restart_model  : reads and packs configuration from a restart file
!
!
!    Input arguments:
!        none
!
!    Output arguments:
!  c_init_data              : character array of initial data
!  c_radhyd_data            : character array of radhyd keys
!  c_eos_data               : character array of edit keys
!  c_nuc_data               : character array of nuclei
!  i_init_data              : integer array of initial data
!  i_radhyd_data            : integer array of radhyd keys
!  i_trans_data             : integer array of transport keys
!  i_e_advct_data           : integer array of e_advect keys
!  i_edit_data              : integer array of transport keys
!  i_hydro_data             : integer array of transport keys
!  i_nuc_data               : integer array of edit keys
!  i_model_data             : integer array of initial model data
!  d_radhyd_data            : 64 bit real array of radhyd keys
!  d_eos_data               : 64 bit real array of edit keys
!  d_trans_data             : 64 bit real array of transport keys
!  d_e_advct_data           : 64 bit real array of e_advect keys
!  d_edit_data              : 64 bit real array of edit keys
!  d_hydro_data             : 64 bit real array of hydro keys
!  d_nuc_data               : 64 bit real array of nuclear keys
!  d_model_data1            : 64 bit real array of initial model data
!  d_model_data2            : 64 bit real array of initial model data
!  d_model_data3            : 64 bit real array of initial model data
!  d_psi_data1              : 64 bit real array of neutrino data
!  d_psi_data2              : 64 bit real array of neutrino data
!  d_psi_data3              : 64 bit real array of neutrino data
!  d_psi_data4              : 64 bit real array of neutrino data
!  d_psi_data5              : 64 bit real array of neutrino data
!  nrst                     : cycle number to start simulation
!
!    Include files:
!  array_module, kind_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, ny, nz, nez, nnu, nnc, nezp1, ij_ray_dim, &
& ik_ray_dim

USE edit_module, ONLY : nrrst, nprint, nlog, data_path
USE parallel_module, ONLY : myid, ierr

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)                              :: rst_file       ! character string containing name of restart_keys file
CHARACTER (len=128)                              :: rst_keys_file  ! character string containing name of restart_keys file
CHARACTER (len=128)                              :: rst_mod_file   ! character string containing name of restart_model file
CHARACTER (len=128)                              :: rst_tmp_file1  ! character string containing name of temporary restart file
CHARACTER (len=128)                              :: rst_tmp_file2  ! character string containing name of temporary restart file
CHARACTER (len=128)                              :: rst_fnl_file   ! character string containing name of final restart file

CHARACTER (len=128), INTENT(out), DIMENSION(1)   :: c_init_data    ! character array of initial data
CHARACTER (len=2  ), INTENT(out), DIMENSION(20)  :: c_radhyd_data  ! character array of radhyd keys
CHARACTER (len=1  ), INTENT(out), DIMENSION(1)   :: c_eos_data     ! character array of edit keys
CHARACTER (len=5  ), INTENT(out), DIMENSION(nnc) :: c_nuc_data     ! character array of nuclei

INTEGER, INTENT(out), DIMENSION(2)               :: i_init_data    ! integer array of initial data
INTEGER, INTENT(out), DIMENSION(50)              :: i_radhyd_data  ! integer array of radhyd keys
INTEGER, INTENT(out), DIMENSION(40+2*nnu)        :: i_trans_data   ! integer array of transport keys
INTEGER, INTENT(out), DIMENSION(5)               :: i_e_advct_data ! integer array of e_advect keys
INTEGER, INTENT(out), DIMENSION(1200+3*40*nnu)   :: i_edit_data    ! integer array of edit keys
INTEGER, INTENT(out), DIMENSION(30)              :: i_hydro_data   ! integer array of transport keys
INTEGER, INTENT(out), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim)   :: i_nuc_data ! integer array of edit keys
INTEGER, INTENT(out), DIMENSION(2)               :: i_model_data   ! integer array of initial model data
INTEGER, INTENT(out)                             :: nrst           ! cycle number to start simulation

REAL(KIND=double), INTENT(out), DIMENSION(50)                            :: d_radhyd_data  ! 64 bit real array of radhyd keys
REAL(KIND=double), INTENT(out), DIMENSION(14)                            :: d_eos_data     ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION((110+3*nnu+3*nez+1))           :: d_trans_data   ! 64 bit real array of transport keys
REAL(KIND=double), INTENT(out), DIMENSION(5+2*nnu)                       :: d_e_advct_data ! 64 bit real array of e_advect keys
REAL(KIND=double), INTENT(out), DIMENSION(50)                            :: d_edit_data    ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION((30+nx))                       :: d_hydro_data   ! 64 bit real array of hydro keys
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

INTEGER, PARAMETER                :: n0 = 0
INTEGER, PARAMETER                :: n1 = 1
INTEGER, PARAMETER                :: n_restart = 23
INTEGER                           :: iskip         ! print flag for echoing the reading of data
INTEGER                           :: istat         ! open-close file flag
INTEGER                           :: nuc_number    ! number of nuclear species (not counting representative heavy nucleus)
INTEGER                           :: nouttmp       ! reatart read flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

   61 FORMAT (' return error in subroutine readst')
  101 FORMAT (' Radhyd keys have been read and packed')
  103 FORMAT (' Eos keys have been read and packed')
  105 FORMAT (' Transport keys have been read and packed')
  107 FORMAT (' E_Advection keys have been read and packed')
  109 FORMAT (' Edit keys have been read and packed')
  111 FORMAT (' Hydro keys have been read and packed')
  113 FORMAT (' Nuclear keys have been read and packed')
  115 FORMAT (' Initial model has been read and packed')
 1001 FORMAT (' nouttmp=',i4,' in read_pack, which is not in the range of 1 - 6; nrst=',i4, &
& ' ny=',i4,' nz=',i4)
 1003 FORMAT (' nouttmp=',i4,' is not allowed with ny=',i4,' and nz=',i4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                  ||||| MYID == 0 THROUGHOUT |||||
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ READ NRST, NOUTTMP, AND HEADER /////
!
!            ||||| ONLY MYID = 0 SHOULD BE IN HERE! |||||
!
!-----------------------------------------------------------------------

CALL read_pack_init( nrrst, c_init_data, i_init_data, nrst, nouttmp )

!-----------------------------------------------------------------------
!
!                     \\\\\ NRST AND NOUTTMP /////
!
!                              nrst = 0
!                              --------
!
!   If nrst = 0, (i.e., a new simulation is being initiated) parameters
!    are initialized and the problem configuration is read from data
!    files in data directory "Initial_Data." If ndim > 1 (ny > 1  or 
!    nz > 1 ), the problem keys and configuration are broadcast to all
!    processors.
!
!                     nrst > 0, ny = 1 and nz = 1
!                     ---------------------------
!
!   If nrst > 0 and ndim = 1 (ny = 1  and nz = 1) the problem keys are
!    read from one of the following files in the "Restart" directory
!    as directed by the value of "nouttmp":
!
!    rst_tmp1_keys.d                    : nouttmp = 1
!    rst_tmp2_keys.d                    : nouttmp = 2
!    restart_keys_XXXXXXX.d             : nouttmp = 3  (XXXXXXX=nrst)
!    restart_final.d                    : nouttmp = 4
!    error                              : nouttmp = 5
!
!   After this the problem configuration is read in from
!
!    rst_tmp1_model.d                   : nouttmp = 1
!    rst_tmp2_model.d                   : nouttmp = 2
!    restart_modelXXXXXXX.d             : nouttmp = 3  (XXXXXXX=nrst)
!    restart_final_mod.d                : nouttmp = 4
!    error                              : nouttmp = 5
!
!                     nrst > 0, ny > 0 or nz > 0
!                     ---------------------------
!
!   If nrst > 0 and ndim > 1 (ny > 1  or  nz > 1) problem keys are read
!    (by processor 0) from one of the following files in the "Restart"
!    directory as directed by the value of "nouttmp," and broadcast to
!    all processors.
!
!    rst_tmp1_keys.d                    : nouttmp = 1
!    rst_tmp2_keys.d                    : nouttmp = 2
!    restart_keys_XXXXXXX.d             : nouttmp = 3  (XXXXXXX=nrst)
!    restart_final_keys.d               : nouttmp = 4
!    restart_keys_XXXXXXX.d             : nouttmp = 5  (XXXXXXX=nrst)
!    restart_keys_XXXXXXX.d             : nouttmp = 6  (XXXXXXX=nrst)
!
!   After this, a few of the quantities common to all nodes are read in
!    (by processor 0) from one of the following files in the "Restart"
!    directory as directed by the value of "nouttmp," and broadcast to
!    all processors. If nouttmp = 5. the entire model configuration is
!    passed to all nodes.
!
!    rst_tmp1_0000.d                    : nouttmp = 1
!    rst_tmp2_0000.d                    : nouttmp = 2
!    restart_model0000_XXXXXXX.d        : nouttmp = 3  (XXXXXXX=nrst)
!    restart_final_mod0000.d            : nouttmp = 4
!    restart_modelXXXXXXX.d             : nouttmp = 5  (XXXXXXX=nrst)
!
!   Finally, execution is passed back to "initlalize_MPI" and the
!    the nuclear abundance data and model configuration are read by all
!    processors from files
!
!    rst_tmp1_YYYY.d                    : nouttmp = 1  (YYYY=myid)
!    rst_tmp2_YYYY.d                    : nouttmp = 2  (YYYY=myid)
!    restart_modelYYYY_XXXXXXX.d        : nouttmp = 3  (YYYY=myid, XXXXXXX=nrst)
!    restart_final_modYYYY.d            : nouttmp = 4  (YYYY=myid)
!    Nothing done                       : nouttmp = 5
!
!-----------------------------------------------------------------------

IF ( nrst == 0 ) THEN

!-----------------------------------------------------------------------
!
!                   \\\\\ READ RADHYD_KEYS /////
!
!-----------------------------------------------------------------------

  CALL radhyd_read( c_radhyd_data, i_radhyd_data, d_radhyd_data, nrst )

!-----------------------------------------------------------------------
!
!                    \\\\\ READ EOS_KEYS /////
!
!-----------------------------------------------------------------------

  CALL eos_read( c_eos_data, d_eos_data, nrst )

!-----------------------------------------------------------------------
!
!                  \\\\\ READ TRANSPORT_KEYS /////
!
!-----------------------------------------------------------------------

  CALL transport_read( i_trans_data, d_trans_data, nrst )

!-----------------------------------------------------------------------
!
!                   \\\\\ READ E_ADVCT_KEYS /////
!
!-----------------------------------------------------------------------

  CALL e_advct_read( i_e_advct_data, d_e_advct_data, nrst )

!-----------------------------------------------------------------------
!
!                    \\\\\ READ EDIT_KEYS /////
!
!-----------------------------------------------------------------------

  CALL edit_read( i_edit_data, d_edit_data, nrst )

!-----------------------------------------------------------------------
!
!                   \\\\\ READ HYDRO_KEYS /////
!
!-----------------------------------------------------------------------

  CALL hydro_read( i_hydro_data, d_hydro_data, nrst )

!-----------------------------------------------------------------------
!
!                  \\\\\ READ NUCLEAR_KEYS /////
!
!-----------------------------------------------------------------------

  CALL nuclear_read( c_nuc_data, i_nuc_data, d_nuc_data, nrst )

!-----------------------------------------------------------------------
!
!                 \\\\\ READ INITIAL_MODEL /////
!
!-----------------------------------------------------------------------

  CALL model_read( i_model_data, d_model_data1, d_model_data2,          &
&  d_model_data3, d_psi_data2, nrst )

!-----------------------------------------------------------------------
!
!                 \\\\\ NRST /= 0, READ MODEL KEYS /////
!
!  Restart problem keys are read from file
!     rst_tmp1_keys.d      : if nouttmp = 1
!     rst_tmp2_keys.d      : if nouttmp = 2
!     restart_keysxxxxx.d  : if nouttmp = 3
!     restart_final_keys.d : if nouttmp = 4
!     restart_keysxxxxx.d  : if nouttmp = 5
!     restart_keysxxxxx.d  : if nouttmp = 6
!  in directory 'Restart'.
!
!-----------------------------------------------------------------------

ELSE ! nrst > 0

!-----------------------------------------------------------------------
!  Open restart_keys file
!-----------------------------------------------------------------------

  IF ( nouttmp == 1 ) THEN
    OPEN (UNIT=n_restart, FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', &
&    STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart, FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', &
&    STATUS='old')

  ELSE IF ( nouttmp == 2 ) THEN
    OPEN (UNIT=n_restart, FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', &
&    STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart, FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', &
&    STATUS='old')
    
  ELSE IF ( nouttmp == 3  .or.  nouttmp == 5  .or.  nouttmp == 6 ) THEN
    WRITE (rst_keys_file,'(a21,i7.7,a2)') '/Restart/restart_keys', nrst,'.d'
    rst_keys_file  = TRIM(data_path)//TRIM(rst_keys_file)
    OPEN (UNIT=n_restart,FILE=TRIM(rst_keys_file), STATUS='new', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart, FILE=TRIM(rst_keys_file), STATUS='old')
    
  ELSE IF ( nouttmp == 4 ) THEN
    OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', &
&    STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', &
&    STATUS='old')

  ELSE
    WRITE (nprint,1001) nouttmp, nrst, ny, nz
    WRITE (nlog,1001) nouttmp, nrst, ny, nz
    STOP

  END IF ! nouttmp == 1

  REWIND n_restart


!-----------------------------------------------------------------------
!  Commemce reading
!-----------------------------------------------------------------------

  iskip            = 0

  CALL read_pack_radhyd_keys( n_restart, nprint, iskip, c_radhyd_data, &
& i_radhyd_data, d_radhyd_data, nrst )
  IF ( myid == 0 ) WRITE (nlog,101)

  CALL read_pack_eos_keys( n_restart, nprint, iskip, c_eos_data, d_eos_data, &
&  nrst )
  IF ( myid == 0 ) WRITE (nlog,103)
  
  CALL read_pack_transport_keys( n_restart, nprint, iskip, nez, nezp1, nnu, &
& i_trans_data, d_trans_data, nrst )
  IF ( myid == 0 ) WRITE (nlog,105)

  CALL read_pack_e_advct_keys( n_restart, nprint, iskip, nnu, i_e_advct_data, &
& d_e_advct_data, nrst )
  IF ( myid == 0 ) WRITE (nlog,107)

  CALL read_pack_edit_keys( n_restart, nprint, iskip, nez, nnu, i_edit_data, &
& d_edit_data, nrst )
  IF ( myid == 0 ) WRITE (nlog,109)

  CALL read_pack_hydro_keys( n_restart, nprint, iskip, nx, i_hydro_data, &
& d_hydro_data, nrst )
  IF ( myid == 0 ) WRITE (nlog,111)

  CALL read_pack_nuclear_keys( n_restart, nprint, iskip, nx, ij_ray_dim, &
&  ik_ray_dim, nnc, c_nuc_data, i_nuc_data, d_nuc_data, nrst, nuc_number )
  IF ( myid == 0 ) WRITE (nlog,113)

!-----------------------------------------------------------------------
!  Cloes restart_keys file
!-----------------------------------------------------------------------

  CLOSE (UNIT=n_restart,STATUS='keep')

!-----------------------------------------------------------------------
!
!                \\\\\ READ MODEL IF PROBLEM IS 1D /////
!
!  If problem is 1D, restart model is read from file
!
!     rst_tmp1_model.d            : if nouttmp = 1
!     rst_tmp2_model.d            : if nouttmp = 2
!     restart_modelxxxxxxx.d      : if nouttmp = 3  (xxxxxxx=nrst)
!     restart_final_mod.d         : if nouttmp = 4
!     error                       : if nouttmp = 5 or 6
!
!   If problem is MD, restart model is read in by each processor
!    from files
!
!     rst_tmp1_xxxx.d             : if nouttmp = 1
!     rst_tmp2_xxxx.d             : if nouttmp = 2
!     restart_modelxxxx_xxxxxxx.d : if nouttmp = 3
!     restart_final_modxxxx.d     : if nouttmp = 4
!     restart_modelxxxxxxx.d      : if nouttmp = 5  (xxxxxxx=nrst)
!
!-----------------------------------------------------------------------

  IF ( ny == 1  .and.  nz == 1 ) THEN

!-----------------------------------------------------------------------
!  Open restart model file
!-----------------------------------------------------------------------

    IF ( nouttmp == 1 ) THEN
      OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/rst_tmp1_model.d', &
&      FORM='unformatted', STATUS='new', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/rst_tmp1_model.d', &
&      FORM='unformatted', STATUS='old')

    ELSE IF ( nouttmp == 2 ) THEN
      OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/rst_tmp2_model.d', &
&      FORM='unformatted', STATUS='new',IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/rst_tmp2_model.d', &
&      FORM='unformatted', STATUS='old')

    ELSE IF ( nouttmp == 3 ) THEN
      WRITE (rst_file,'(a22,i7.7,a2)') '/Restart/restart_model', nrst,'.d'
      rst_file       = TRIM(data_path)//TRIM(rst_file)
      OPEN (UNIT=n_restart, FILE=TRIM(rst_file), FORM='unformatted', &
&      STATUS='new', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=n_restart, FILE=TRIM(rst_file), &
&      FORM='unformatted', STATUS='old')
    
    ELSE IF ( nouttmp == 4 ) THEN
      OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/restart_final_mod.d', &
&      FORM='unformatted', STATUS='new',IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=n_restart,FILE=TRIM(data_path)//'/Restart/restart_final_mod.d', &
&      FORM='unformatted', STATUS='old')
    
    ELSE IF ( nouttmp == 5  .or.  nouttmp == 6 ) THEN
      WRITE (nprint,1003) nouttmp, ny, nz
      WRITE (nlog,1003) nouttmp, ny, nz
      STOP

    ELSE
      WRITE (nprint,1001) nouttmp, nrst, ny, nz
      WRITE (nlog,1001) nouttmp, nrst, ny, nz
      STOP

    END IF ! nouttmp == 1

    REWIND n_restart

!-----------------------------------------------------------------------
!  Commemce reading
!-----------------------------------------------------------------------

    iskip          = 0

    CALL read_pack_restart_model( n_restart, nprint, iskip, nx, nez, nnu,     &
&    nnc, ij_ray_dim, ik_ray_dim, i_model_data, d_model_data1, d_model_data2, &
&    d_model_data3, i_nuc_data, d_nuc_data, d_psi_data1, d_psi_data2,         &
&    d_psi_data3, d_psi_data4, d_psi_data5, nuc_number, nrst )
    IF ( myid == 0 ) WRITE (nlog,115)

!-----------------------------------------------------------------------
!  Cloes restart_mod file
!-----------------------------------------------------------------------

    CLOSE (UNIT=n_restart,STATUS='keep')

!-----------------------------------------------------------------------
!
!            \\\\\ READ MODEL TO NODE 0 IF PROBLEM IS MD /////
!
!  If problem is MD, restart model is read from the file specified
!   in directory 'Restart' in order to get a few of the quantities
!   common to all nodes.
!
!-----------------------------------------------------------------------

  ELSE ! ny /= 0 and/or nz /= 0

    IF ( nouttmp == 5 ) THEN
    
      WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&      nrst,'.d'
      rst_mod_file   = TRIM(data_path)//TRIM(rst_mod_file)
      OPEN (UNIT=n_restart, FILE=TRIM(rst_mod_file), STATUS='new', &
&      FORM='unformatted', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=n_restart, FILE=TRIM(rst_mod_file), &
&      FORM='unformatted', STATUS='old')

      REWIND n_restart

!-----------------------------------------------------------------------
!  Commemce reading
!-----------------------------------------------------------------------

      iskip            = 0

      CALL read_pack_restart_model( n_restart, nprint, iskip, nx, nez, nnu,     &
&      nnc, ij_ray_dim, ik_ray_dim, i_model_data, d_model_data1, d_model_data2, &
&      d_model_data3, i_nuc_data, d_nuc_data, d_psi_data1, d_psi_data2,         &
&      d_psi_data3, d_psi_data4, d_psi_data5, nuc_number, nrst )

!-----------------------------------------------------------------------
!  Cloes restart_mod file
!-----------------------------------------------------------------------

      CLOSE (UNIT=n_restart,STATUS='keep')

    END IF ! nouttmp == 5

  END IF ! ny == 1  and  nz == 1

END IF ! nrst == 0

RETURN
END SUBROUTINE read_pack
