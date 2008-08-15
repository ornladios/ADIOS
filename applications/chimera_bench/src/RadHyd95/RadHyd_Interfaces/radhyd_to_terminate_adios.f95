SUBROUTINE radhyd_to_terminate_adios ( ij_ray_min, ij_ray_max, ij_ray_dim, ik_ray_min, &
& ik_ray_max, ik_ray_dim, nx, nez, nnu, nnc, ndim )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_terminate
!    Module:       radhyd_to_terminate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To load variables into 3-d arrays for porting into MGFLD edit
!         via subroutine edit_in.
!
!    Subprograms called
!  edit_term_in     : sets up variables and prints the final edit
!  restart1D_write  : prints the final restart file
!
!    Input arguments
!  ij_ray_min : minimum j-index of radial ray
!  ij_ray_max : maximum j-index of radial ray
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_min : minimum k-index of radial ray
!  ik_ray_max : maximum k-index of radial ray
!  ik_ray_dim : number of z-zones on a processor before swapping
!  nx               : x_array extent
!  nez              : neutrino energy array extent
!  nnu              : neutrino flavor array extent
!  nnc              : neutrino abundance array extent
!  ndim             : number of spatial dimensions of the simulation
!
!    Output arguments:
!        none
!
!    Include files:
!      kind_module, numerical_module
!      edit_module, parallel_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE edit_module, ONLY : ncychg, data_path, nlog
USE parallel_module, ONLY : myid
USE mpi
USE radial_ray_module, ONLY : imin, imax, rho_ci, rho_c, t_ci, t_c, ye_ci, &
& ye_c, u_c, x_ef, psi0_c, psi1_e, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, &
& nse_c, dtnph, time, t_bounce, ncycle, ncymax, t_stop, tb_stop, nprint
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray_min       ! minimum j-index of radial ray
INTEGER, INTENT(in)              :: ij_ray_max       ! maximum j-index of radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_min       ! minimum k-index of radial ray
INTEGER, INTENT(in)              :: ik_ray_max       ! maximum k-index of radial ray
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc              ! composition array extent
INTEGER, INTENT(in)              :: ndim             ! number of spatial dimensions of the simulation

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=128)              :: rst_mod_file     ! character string containing name of restart_model file

INTEGER                          :: is               ! minimum x-array index
INTEGER                          :: ie               ! maximum x-array index
INTEGER                          :: idim             ! x-array extent
INTEGER                          :: ij_ray           ! j-index of a radial ray
INTEGER                          :: ik_ray           ! k-index of a radial ray
INTEGER, PARAMETER               :: n_dump = 23      ! unit number for final restart dump
INTEGER                          :: istat            ! open file flag
INTEGER                          :: ierr             ! mpi flag
INTEGER                          :: adios_err        ! adios error flag

#ifdef ADIOS_KEYS || ADIOS_MODEL
character (len=128) :: restart_filename
#endif

 1001 FORMAT (' Job Done')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                  \\\\\ TEST FOR TERMINATION /////
!
!-----------------------------------------------------------------------

IF ( t_bounce == zero ) THEN
  IF ( ncycle < ncymax  .and.  time < t_stop ) RETURN
ELSE
  IF ( ncycle < ncymax  .and.  time < t_stop  .and.  time - t_bounce < tb_stop ) RETURN
END IF

!-----------------------------------------------------------------------
!
!                       \\\\\ FIMAL EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Load integer scalars
!-----------------------------------------------------------------------

is                    = imin
ie                    = imax
idim                  = nx

!-----------------------------------------------------------------------
!              ||||| Begin loop over radial rays |||||
!-----------------------------------------------------------------------

DO ik_ray = ik_ray_min,ik_ray_max
  DO ij_ray = ij_ray_min,ij_ray_max

!-----------------------------------------------------------------------
!  Transfer variables to for editing
!-----------------------------------------------------------------------

    CALL edit_term_in( is, ie, idim, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&    nx, nez, nnu, nnc, rho_ci, rho_c, t_ci, t_c, ye_ci, ye_c,  x_ef, u_c, &
&    psi0_c, psi1_e, dtnph, time, ncycle, xn_c, be_nuc_rep_c, a_nuc_rep_c, &
&    z_nuc_rep_c, nse_c )

!-----------------------------------------------------------------------
!               ||||| End loop over radial rays |||||
!-----------------------------------------------------------------------

  END DO ! ij_ray = ij_ray_min,ij_ray_max
END DO ! ik_ray = ik_ray_min,ik_ray_max

!-----------------------------------------------------------------------
!
!                   \\\\\ FIMAL RESTART DUMP /////
!
!-----------------------------------------------------------------------

IF ( ndim == 1 ) THEN

#ifdef ADIOS_KEYS

  restart_filename = TRIM(data_path)//'/Restart/restart_final_keys.d'

  IF ( ncycle == ncymax) ncymax = ncymax + ncychg

!-----------------------------------------------------------------------
!  Print the final restart keys file
!-----------------------------------------------------------------------

  CALL restart_keys_write_adios( restart_filename, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the final restart keys dump file
!-----------------------------------------------------------------------

  OPEN (UNIT=n_dump,FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', &
&  STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=n_dump,FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', &
&  STATUS='old')

  IF ( ncycle == ncymax) ncymax = ncymax + ncychg

!-----------------------------------------------------------------------
!  Print the final restart keys file
!-----------------------------------------------------------------------

  CALL restart_keys_write( n_dump, nez, nnu )

!-----------------------------------------------------------------------
!  Close the final restart keys dump file
!-----------------------------------------------------------------------

  CLOSE (UNIT=n_dump,STATUS='keep')

#endif

#ifdef ADIOS_MODEL

  WRITE (rst_mod_file,'(a28)') '/Restart/restart_final_mod.d'
   rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

  CALL restart_model_write_adios ( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the final restart model file
!-----------------------------------------------------------------------

  WRITE (rst_mod_file,'(a28)') '/Restart/restart_final_mod.d'
   rst_mod_file = TRIM(data_path)//TRIM(rst_mod_file)
  OPEN (UNIT=n_dump, FILE=TRIM(rst_mod_file), STATUS='new', &
&   FORM='unformatted', IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=n_dump, FILE=TRIM(rst_mod_file), &
&   FORM='unformatted', STATUS='old')
 
!-----------------------------------------------------------------------
!  Print the restart model file
!-----------------------------------------------------------------------

  CALL restart_model_write( n_dump, nx, nnu )

!-----------------------------------------------------------------------
!  Close the restart model file
!-----------------------------------------------------------------------

  CLOSE (UNIT=n_dump,STATUS='keep')

#endif

ELSE ! ndim > 1

  IF (myid == 0) THEN

#ifdef ADIOS_KEYS

    restart_filename = TRIM(data_path)//'/Restart/restart_final_keys.d'

    IF ( ncycle == ncymax) ncymax = ncymax + ncychg

!-----------------------------------------------------------------------
!  Print the final restart filele
!-----------------------------------------------------------------------

    CALL restart_keys_write_adios ( restart_filename, nez, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the final restart keys file
!-----------------------------------------------------------------------

    OPEN (UNIT=n_dump,FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', &
&  STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=n_dump,FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', &
&  STATUS='old')

    IF ( ncycle == ncymax) ncymax = ncymax + ncychg

!-----------------------------------------------------------------------
!  Print the final restart filele
!-----------------------------------------------------------------------

    CALL restart_keys_write( n_dump, nez, nnu )

!-----------------------------------------------------------------------
!  Close the final restart keys file
!-----------------------------------------------------------------------

    CLOSE (UNIT=n_dump,STATUS='keep')

#endif

  END IF ! myid == 0

#ifdef ADIOS_MODEL

  WRITE (rst_mod_file,'(a26,i4.4,a2)') '/Restart/restart_final_mod', &
&  myid,'.d'
  rst_mod_file      = TRIM(data_path)//TRIM(rst_mod_file)

  IF ( ncycle == ncymax) ncymax = ncymax + ncychg

!-----------------------------------------------------------------------
!  Print the final restart model file
!-----------------------------------------------------------------------

  CALL restart_model_write_adios( rst_mod_file, nx, nnu )

#else

!-----------------------------------------------------------------------
!  Name and open the final restart model file
!-----------------------------------------------------------------------

  WRITE (rst_mod_file,'(a26,i4.4,a2)') '/Restart/restart_final_mod', &
&  myid,'.d'
  rst_mod_file      = TRIM(data_path)//TRIM(rst_mod_file)
  OPEN (UNIT=n_dump, FILE=TRIM(rst_mod_file), STATUS='new', &
&  FORM='unformatted', IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=n_dump, FILE=TRIM(rst_mod_file), &
&  FORM='unformatted', STATUS='old')

  IF ( ncycle == ncymax) ncymax = ncymax + ncychg

!-----------------------------------------------------------------------
!  Print the final restart model file
!-----------------------------------------------------------------------

  CALL restart_model_write( n_dump, nx, nnu )

!-----------------------------------------------------------------------
!  Close the final restart model file
!-----------------------------------------------------------------------

  CLOSE (UNIT=n_dump,STATUS='keep')

#endif

END IF ! ndim == 1

!-----------------------------------------------------------------------
!  Job done
!-----------------------------------------------------------------------

WRITE (nprint,1001)
WRITE (nlog,1001)

#ifdef ADIOS_MODEL || ADIOS_KEYS
CALL adios_finalize (myid, adios_err)
#endif

#ifdef ADIOS_PROFILING
  CALL cycle_end(ncycle)
  CALL finalize_prof()
#endif


CALL MPI_FINALIZE(ierr)
STOP

END SUBROUTINE radhyd_to_terminate_adios
