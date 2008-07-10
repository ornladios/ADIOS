SUBROUTINE radhyd_to_restart_adios( ndim, nx, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_restart
!    Module:       radhyd_to_restart
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/27/05
!
!    Purpose:
!        To load variables arrays for writing a restart file.
!
!    Subprograms called:
!  restart1D_write : Writes the restart files
!
!    Input arguments:
!  ndim            : number of spatial dimensions of the simulation
!  nx              : x-array extent
!  nez             : neutrino energy array extent
!  nnu             : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, edit_module, radial_ray_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE cycle_module, ONLY: ncycle, ncymax
USE edit_module, ONLY: nnrst, ncyrst, nprm, intprm, noutpmt, intrst, nrstd2, &
& nrstd1, nrstfl, dtcomplot, nprint, data_path, reset_path, nlog
USE parallel_module, ONLY : myid
USE t_cntrl_module, ONLY: time, t_bounce, dtnph

     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ndim             ! number of spatial dimensions of the simulation
INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=128)              :: rstfile          ! character string containing name of restart file
CHARACTER (len=128)              :: rst_key_file     ! character string containing name of restart_keys file
CHARACTER (len=128)              :: rst_mod_file     ! character string containing name of restart_model file
CHARACTER (len=128)              :: rst_tmp_file1    ! character string containing name of temporary restart file
CHARACTER (len=128)              :: rst_tmp_file2    ! character string containing name of temporary restart file

INTEGER                          :: i                ! do index
INTEGER                          :: nnrst_save       ! temporarily stored value of nnrst
INTEGER                          :: ndump            ! unit number to write temporary restart dumps
INTEGER                          :: istat            ! open file flag
INTEGER                          :: itime            ! used to determine permanent restart file
INTEGER                          :: itime_prev       ! used to determine permanent restart file

REAL(KIND=double)                :: tmult            ! factor to determine when to write permanent restart dump
REAL(KIND=double)                :: t_tb             ! time from core bounce

#ifdef ADIOS_KEYS || ADIOS_MODEL
character (len=128) :: restart_filename
#else
#endif

 1001 FORMAT (' Problem closing file superdump.d in subroutine mgfld_edit')
 1003 FORMAT (' Problem reopening file superdump.d in subroutine mgfld_edit')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Increment edit counters
!-----------------------------------------------------------------------

nprm                 = nprm  + 1
nnrst                = nnrst + 1

!-----------------------------------------------------------------------
!
!             \\\\\ WRITE A PERMANENT RESTART FILE /////
!
!  Write to permanent restart file if ncycle = ncyrst(i) for some i.
!         
!  nnrst (the temporary restart dump counter) is temporarily set to zero
!   to be recorded as such in the restart file.
!
!-----------------------------------------------------------------------

DO i = 1,100
  IF ( ncycle == ncyrst(i) ) THEN
    nnrst_save       = nnrst
    nnrst            = 0 

    IF ( nprm >= intprm ) nprm = 0      ! Don't write restart file twice

    IF ( ndim == 1 ) THEN

#ifdef ADIOS_KEYS

      WRITE (rstfile,'(a21,i7.7,a2)') '/Restart/restart_keys',ncycle,'.d'
      rstfile        = TRIM(data_path)//TRIM(rstfile)

!-----------------------------------------------------------------------
!  Print the restart keys file
!-----------------------------------------------------------------------

      CALL restart_keys_write_adios ( rstfile, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart keys file
!-----------------------------------------------------------------------

      WRITE (rstfile,'(a21,i7.7,a2)') '/Restart/restart_keys',ncycle,'.d'
      rstfile        = TRIM(data_path)//TRIM(rstfile)
      OPEN (UNIT=noutpmt, FILE=TRIM(rstfile), STATUS='new', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rstfile), STATUS='old')

!-----------------------------------------------------------------------
!  Print the restart keys file
!-----------------------------------------------------------------------

      CALL restart_keys_write( noutpmt, nez, nnu )

!-----------------------------------------------------------------------
!  Close the restart keys file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

#ifdef ADIOS_MODEL
      WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&      ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write_adios( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart model file
!-----------------------------------------------------------------------

      WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&      ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
      OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), STATUS='new', &
&      FORM='unformatted', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), &
&      FORM='unformatted', STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write( noutpmt, nx, nnu )

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

    ELSE ! ndim > 1

!-----------------------------------------------------------------------
!                  ||||| MYID = 0 STARTS HERE |||||
!-----------------------------------------------------------------------

      IF ( myid == 0 ) THEN

#ifdef ADIOS_KEYS

        WRITE (rst_key_file,'(a21,i7.7,a2)') '/Restart/restart_keys', ncycle,'.d'
        rst_key_file = TRIM(data_path)//TRIM(rst_key_file)
 
!-----------------------------------------------------------------------
!  Print the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CALL restart_keys_write_adios ( rst_key_file, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart keys file if myid = 0
!-----------------------------------------------------------------------

        WRITE (rst_key_file,'(a21,i7.7,a2)') '/Restart/restart_keys', ncycle,'.d'
        rst_key_file = TRIM(data_path)//TRIM(rst_key_file)
        OPEN (UNIT=noutpmt, FILE=TRIM(rst_key_file), STATUS='new', IOSTAT=istat)
        IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_key_file), STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CALL restart_keys_write( noutpmt, nez, nnu )

!-----------------------------------------------------------------------
!  Close the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

!-----------------------------------------------------------------------
!                   ||||| MYID = 0 ENDS HERE |||||
!-----------------------------------------------------------------------

      END IF ! myid = 0

#ifdef ADIOS_MODEL

      WRITE (rst_mod_file,'(a22,i4.4,a1,i7.7,a2)') '/Restart/restart_model', &
&      myid,'_',ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write_adios( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart model file
!-----------------------------------------------------------------------

      WRITE (rst_mod_file,'(a22,i5.5,a1,i7.7,a2)') '/Restart/restart_model', &
&      myid,'_',ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
      OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), STATUS='new', &
&       FORM='unformatted', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), &
&       FORM='unformatted', STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write( noutpmt, nx, nnu )

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

    END IF ! ndim == 1
    nnrst            = nnrst_save
  END IF !  ncycle = ncyrst(i)
END DO ! i = 1,100

!-----------------------------------------------------------------------
!
!             \\\\\ WRITE A PERMANENT RESTART FILE /////
!
!  Write to permanent restart file if nprm = intprm
!
!  nnrst (the temporary restart dump counter) is temporarily set to zero
!   to be recorded as such in the restart file.
!
!-----------------------------------------------------------------------

IF ( nprm >= intprm ) THEN
  nnrst_save         = nnrst
  nnrst              = 0
  nprm               = 0

  IF ( ndim == 1 ) THEN

#ifdef ADIOS_KEYS

    WRITE (rstfile,'(a21,i7.7,a2)') '/Restart/restart_keys',ncycle,'.d'
    rstfile        = TRIM(data_path)//TRIM(rstfile)

!-----------------------------------------------------------------------
!  Print the restart keys file
!-----------------------------------------------------------------------

    CALL restart_keys_write_adios ( rstfile, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart keys file
!-----------------------------------------------------------------------

    WRITE (rstfile,'(a21,i7.7,a2)') '/Restart/restart_keys',ncycle,'.d'
    rstfile        = TRIM(data_path)//TRIM(rstfile)
    OPEN (UNIT=noutpmt, FILE=TRIM(rstfile), STATUS='new', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rstfile), STATUS='old')

!-----------------------------------------------------------------------
!  Print the restart keys file
!-----------------------------------------------------------------------

    CALL restart_keys_write( noutpmt, nez, nnu )

!-----------------------------------------------------------------------
!  Close the restart keys file
!-----------------------------------------------------------------------

    CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

#ifdef ADIOS_MODEL

    WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&    ncycle,'.d'
    rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write_adios( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart model file
!-----------------------------------------------------------------------

    WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&    ncycle,'.d'
    rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
    OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), &
&    FORM='unformatted', STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write( noutpmt, nx, nnu )

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

    ELSE ! ndim > 1

!-----------------------------------------------------------------------
!                  ||||| MYID = 0 STARTS HERE |||||
!-----------------------------------------------------------------------

      IF ( myid == 0 ) THEN

#ifdef ADIOS_KEYS

        WRITE (rst_key_file,'(a21,i7.7,a2)') '/Restart/restart_keys', ncycle,'.d'
        rst_key_file = TRIM(data_path)//TRIM(rst_key_file)
 
!-----------------------------------------------------------------------
!  Print the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CALL restart_keys_write_adios( rst_key_file, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart keys file if myid = 0
!-----------------------------------------------------------------------

        WRITE (rst_key_file,'(a21,i7.7,a2)') '/Restart/restart_keys', ncycle,'.d'
        rst_key_file = TRIM(data_path)//TRIM(rst_key_file)
        OPEN (UNIT=noutpmt, FILE=TRIM(rst_key_file), STATUS='new', IOSTAT=istat)
        IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_key_file), STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CALL restart_keys_write( noutpmt, nez, nnu )

!-----------------------------------------------------------------------
!  Close the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

!-----------------------------------------------------------------------
!                   ||||| MYID = 0 ENDS HERE |||||
!-----------------------------------------------------------------------

      END IF ! myid = 0

#ifdef ADIOS_MODEL

      WRITE (rst_mod_file,'(a22,i4.4,a1,i7.7,a2)') '/Restart/restart_model', &
&      myid,'_',ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)

!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write_adios( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart model file
!-----------------------------------------------------------------------

      WRITE (rst_mod_file,'(a22,i5.5,a1,i7.7,a2)') '/Restart/restart_model', &
&      myid,'_',ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
      OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), STATUS='new', &
&      FORM='unformatted', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), &
&      FORM='unformatted', STATUS='old')

!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write( noutpmt, nx, nnu )

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

    END IF ! ndim == 1
  nnrst              = nnrst_save
END IF ! nprm >= intprm

!-----------------------------------------------------------------------
!
!             \\\\\ WRITE A PERMANENT RESTART FILE /////
!
!  Write to permanent restart file if t_bounce is a multiple dtcomplot.
!-----------------------------------------------------------------------

IF ( t_bounce /= zero ) THEN
  t_tb               = time - t_bounce
  tmult              = 1.d+3/dtcomplot
  itime              = int( t_tb * tmult )
  itime_prev         = int( ( t_tb - dtnph ) * tmult )
  IF ( itime /= itime_prev ) THEN
    nnrst_save       = nnrst
    nnrst            = 0 

    IF ( ndim == 1 ) THEN

#ifdef ADIOS_KEYS

      WRITE (rstfile,'(a21,i7.7,a2)') '/Restart/restart_keys',ncycle,'.d'
      rstfile        = TRIM(data_path)//TRIM(rstfile)

!-----------------------------------------------------------------------
!  Print the restart keys file
!-----------------------------------------------------------------------

      CALL restart_keys_write_adios( rstfile, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart keys file
!-----------------------------------------------------------------------

      WRITE (rstfile,'(a21,i7.7,a2)') '/Restart/restart_keys',ncycle,'.d'
      rstfile        = TRIM(data_path)//TRIM(rstfile)
      OPEN (UNIT=noutpmt, FILE=TRIM(rstfile), STATUS='new', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rstfile), STATUS='old')

!-----------------------------------------------------------------------
!  Print the restart keys file
!-----------------------------------------------------------------------

      CALL restart_keys_write( noutpmt, nez, nnu )

!-----------------------------------------------------------------------
!  Close the restart keys file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

#ifdef ADIOS_MODEL

      WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&      ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write_adios( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart model file
!-----------------------------------------------------------------------

      WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&      ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
      OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), STATUS='new', &
&      FORM='unformatted', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), &
&      FORM='unformatted', STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write( noutpmt, nx, nnu )

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

    ELSE ! ndim > 1

!-----------------------------------------------------------------------
!                  ||||| MYID = 0 STARTS HERE |||||
!-----------------------------------------------------------------------

      IF ( myid == 0 ) THEN

#ifdef ADIOS_KEYS

        WRITE (rst_key_file,'(a21,i7.7,a2)') '/Restart/restart_keys', ncycle,'.d'
        rst_key_file = TRIM(data_path)//TRIM(rst_key_file)
 
!-----------------------------------------------------------------------
!  Print the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CALL restart_keys_write_adios( rst_key_file, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart keys file if myid = 0
!-----------------------------------------------------------------------

        WRITE (rst_key_file,'(a21,i7.7,a2)') '/Restart/restart_keys', ncycle,'.d'
        rst_key_file = TRIM(data_path)//TRIM(rst_key_file)
        OPEN (UNIT=noutpmt, FILE=TRIM(rst_key_file), STATUS='new', IOSTAT=istat)
        IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_key_file), STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CALL restart_keys_write( noutpmt, nez, nnu )

!-----------------------------------------------------------------------
!  Close the restart keys file if myid = 0
!-----------------------------------------------------------------------

        CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

!-----------------------------------------------------------------------
!                   ||||| MYID = 0 ENDS HERE |||||
!-----------------------------------------------------------------------

      END IF ! myid = 0

#ifdef ADIOS_MODEL

      WRITE (rst_mod_file,'(a22,i7.7,a2)') '/Restart/restart_model', &
&      ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write_adios( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the restart model file
!-----------------------------------------------------------------------

      WRITE (rst_mod_file,'(a22,i5.5,a1,i7.7,a2)') '/Restart/restart_model', &
&      myid,'_',ncycle,'.d'
      rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
      OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), STATUS='new', &
&      FORM='unformatted', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=noutpmt, FILE=TRIM(rst_mod_file), &
&      FORM='unformatted', STATUS='old')

!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

      CALL restart_model_write( noutpmt, nx, nnu )

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

      CLOSE (UNIT=noutpmt,STATUS='keep')

#endif

    END IF ! ndim == 1
    nnrst            = nnrst_save
  END IF ! itime /= itime_prev
END IF ! t_bounce /= zero

!-----------------------------------------------------------------------
!
!             \\\\\ WRITE A TEMPORARY RESTART FILE /////
!
!  Write to restart file if nnrst = intrst.
!
!-----------------------------------------------------------------------

IF ( nnrst >= intrst ) THEN
  nnrst              = 0

  IF ( ndim == 1 ) THEN

#ifdef ADIOS_KEYS

    restart_filename = TRIM(data_path)//'/Restart/rst_tmp1_keys.d'
    !restart_filename = TRIM(data_path)//'/Restart/rst_tmp2_keys.d'

!-----------------------------------------------------------------------
!  Print the temporary restart keys file
!-----------------------------------------------------------------------

    CALL restart_keys_write_adios ( restart_filename, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the temporary restart keys file
!-----------------------------------------------------------------------

    OPEN (UNIT=nrstd1,FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', &
&    STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nrstd1,FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', &
&    STATUS='old')
    OPEN (UNIT=nrstd2,FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', &
&    STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nrstd2,FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', &
&    STATUS='old')

    REWIND nrstd1
    REWIND nrstd2

    IF ( nrstfl == 2 ) THEN
      IF ( ndump == nrstd1 ) THEN
        ndump        = nrstd2
      ELSE
        ndump        = nrstd1
      END IF ! ndump = nrstd1
    END IF ! nrstfl = 2
 
!-----------------------------------------------------------------------
!  Print the temporary restart keys file
!-----------------------------------------------------------------------

    CALL restart_keys_write( ndump, nez, nnu )

!-----------------------------------------------------------------------
!  Close the temporary restart keys file
!-----------------------------------------------------------------------

    CLOSE (UNIT=nrstd1,STATUS='keep')
    CLOSE (UNIT=nrstd2,STATUS='keep')

#endif

#ifdef ADIOS_MODEL
    WRITE (rst_tmp_file1,'(a25)') '/Restart/rst_tmp1_model.d'
    rst_tmp_file1  = TRIM(data_path)//TRIM(rst_tmp_file1)
    WRITE (rst_tmp_file2,'(a25)') '/Restart/rst_tmp2_model.d'
    rst_tmp_file2  = TRIM(data_path)//TRIM(rst_tmp_file2)

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

    CALL restart_model_write_adios( rst_tmp_file1, nx, nnu )
    CALL restart_model_write_adios( rst_tmp_file2, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the temporary restart model file
!-----------------------------------------------------------------------

    WRITE (rst_tmp_file1,'(a25)') '/Restart/rst_tmp1_model.d'
    rst_tmp_file1  = TRIM(data_path)//TRIM(rst_tmp_file1)
    WRITE (rst_tmp_file2,'(a25)') '/Restart/rst_tmp2_model.d'
    rst_tmp_file2  = TRIM(data_path)//TRIM(rst_tmp_file2)

    OPEN (UNIT=nrstd1, FILE=rst_tmp_file1, STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nrstd1, FILE=TRIM(rst_tmp_file1), &
&    FORM='unformatted', STATUS='old')
    OPEN (UNIT=nrstd2, FILE=rst_tmp_file2, STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nrstd2, FILE=TRIM(rst_tmp_file2), &
&    FORM='unformatted', STATUS='old')

    REWIND nrstd1
    REWIND nrstd2
 
!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

    CALL restart_model_write( ndump, nx, nnu )

!-----------------------------------------------------------------------
!  Close the temporary restart model file
!-----------------------------------------------------------------------

    CLOSE (UNIT=nrstd1,STATUS='keep')
    CLOSE (UNIT=nrstd2,STATUS='keep')

#endif

  ELSE ! ndim > 1

!-----------------------------------------------------------------------
!  Alternate restart files if nrstfl = 2
!-----------------------------------------------------------------------

    IF ( nrstfl == 2 ) THEN
      IF ( ndump == nrstd1 ) THEN
        ndump        = nrstd2
      ELSE
        ndump        = nrstd1
      END IF ! ndump = nrstd1
    END IF ! nrstfl = 2

!-----------------------------------------------------------------------
!                  ||||| MYID = 0 STARTS HERE |||||
!-----------------------------------------------------------------------

    IF ( myid == 0 ) THEN

#ifdef ADIOS_KEYS

      restart_filename = TRIM(data_path)//'/Restart/rst_tmp1_keys.d'
      !restart_filename = TRIM(data_path)//'/Restart/rst_tmp2_keys.d'
 
!-----------------------------------------------------------------------
!  Print the temporary restart keys file if myid = 0
!-----------------------------------------------------------------------

      CALL restart_keys_write_adios( restart_filename, nez, nnu )

!-----------------------------------------------------------------------
!  Close the temporary restart keys file if myid = 0
!-----------------------------------------------------------------------

#else

!-----------------------------------------------------------------------
!  Name and open the temporary restart keys file if myid = 0
!-----------------------------------------------------------------------

      OPEN (UNIT=nrstd1,FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', &
&      STATUS='new',IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=nrstd1,FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', &
&      STATUS='old')
      OPEN (UNIT=nrstd2,FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', &
&      STATUS='new',IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=nrstd2,FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', &
&      STATUS='old')

      REWIND nrstd1
      REWIND nrstd2
 
!-----------------------------------------------------------------------
!  Print the temporary restart keys file if myid = 0
!-----------------------------------------------------------------------

      CALL restart_keys_write( ndump, nez, nnu )

!-----------------------------------------------------------------------
!  Close the temporary restart keys file if myid = 0
!-----------------------------------------------------------------------

      CLOSE (UNIT=nrstd1,STATUS='keep')
      CLOSE (UNIT=nrstd2,STATUS='keep')

#endif

!-----------------------------------------------------------------------
!                  ||||| MYID = 0 STARTS HERE |||||
!-----------------------------------------------------------------------

    END IF ! myid = 0

#ifdef ADIOS_MODEL

    WRITE (rst_tmp_file1,'(a18,i4.4,a2)') '/Restart/rst_tmp1_', &
&    myid,'.d'
    rst_tmp_file1  = TRIM(data_path)//TRIM(rst_tmp_file1)
    WRITE (rst_tmp_file2,'(a18,i4.4,a2)') '/Restart/rst_tmp2_', &
&    myid,'.d'
    rst_tmp_file2  = TRIM(data_path)//TRIM(rst_tmp_file2)

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

    CALL restart_model_write_adios( rst_tmp_file1, nx, nnu )
    CALL restart_model_write_adios( rst_tmp_file2, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the temporary restart model file
!-----------------------------------------------------------------------

    WRITE (rst_tmp_file1,'(a18,i5.5,a2)') '/Restart/rst_tmp1_', &
&    myid,'.d'
    rst_tmp_file1  = TRIM(data_path)//TRIM(rst_tmp_file1)
    WRITE (rst_tmp_file2,'(a18,i5.5,a2)') '/Restart/rst_tmp2_', &
&    myid,'.d'
    rst_tmp_file2  = TRIM(data_path)//TRIM(rst_tmp_file2)

    OPEN (UNIT=nrstd1, FILE=rst_tmp_file1, STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nrstd1, FILE=TRIM(rst_tmp_file1), &
&    FORM='unformatted', STATUS='old')
    OPEN (UNIT=nrstd2, FILE=rst_tmp_file2, STATUS='new', &
&    FORM='unformatted', IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nrstd2, FILE=TRIM(rst_tmp_file2), &
&    FORM='unformatted', STATUS='old')

    REWIND nrstd1
    REWIND nrstd2
 
!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

    CALL restart_model_write( ndump, nx, nnu )

!-----------------------------------------------------------------------
!  Close the temporary restart model file
!-----------------------------------------------------------------------

    CLOSE (UNIT=nrstd1,STATUS='keep')
    CLOSE (UNIT=nrstd2,STATUS='keep')

#endif

  END IF ! ndim == 1

!-----------------------------------------------------------------------
!                  ||||| MYID = 0 STARTS HERE |||||
!-----------------------------------------------------------------------

  IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Record the location of the restart file in reset.d if myid = 0.
!-----------------------------------------------------------------------

    OPEN (UNIT=ndump,FILE=TRIM(reset_path)//'/reset.d', STATUS='new', &
&    IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=ndump, FILE=TRIM(reset_path)//'/reset.d', &
&    STATUS='old')
    REWIND ndump
    CALL init_write( ndump )
    CLOSE (UNIT=ndump, STATUS='keep')

!-----------------------------------------------------------------------
!  Close and reopen Data3/superdump.d
!-----------------------------------------------------------------------

    CLOSE (UNIT=nprint,STATUS='keep',IOSTAT=istat)
    IF ( istat /= 0 ) WRITE (nlog,1001)

    OPEN (UNIT=nprint,FILE=TRIM(data_path)//'/Run_Log/superdump.d', &
&    STATUS='old', POSITION='append', IOSTAT=istat)
    IF ( istat /= 0 ) WRITE (nlog,1003)

!-----------------------------------------------------------------------
!                   ||||| MYID = 0 ENDS HERE |||||
!-----------------------------------------------------------------------

  END IF ! myid = 0

END IF ! nnrst >= intrst

RETURN
END SUBROUTINE radhyd_to_restart_adios
