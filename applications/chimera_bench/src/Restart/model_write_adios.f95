SUBROUTINE model_write_adios( ndump, nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         model_write
!    Module:       model_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/25/05
!
!    Purpose:
!      To dump the model configuration.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ndump       : unit number to dump the data
!  nx          : x-array extent
!  nnu         : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  edit_modulee, eos_snc_x_module, nu_dist_module, nu_energy_grid_module,
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE array_module, ONLY : nez, nnc, ij_ray_dim, ik_ray_dim

USE edit_module, ONLY : nlog, nu_r, nu_rt, nu_rho, nu_rhot, psi0dat, psi1dat

USE eos_snc_x_module, ONLY : duesrc

USE nu_dist_module, ONLY : dnurad, unukrad, unujrad, nnukrad, nnujrad, &
& e_rad, unurad, elec_rad, nnurad

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

USE radial_ray_module, ONLY : imin, imax, ncycle, nprint, nse_c, rho_c, &
& t_c, ye_c, u_c, v_c, w_c, x_ef, dx_cf, psi0_c, xn_c, a_nuc_rep_c, &
& z_nuc_rep_c, be_nuc_rep_c, uburn_c, e_nu_c_bar, f_nu_e_bar

!old
USE edit_module, ONLY : nlog

USE nu_dist_module, ONLY : dnurad, unukrad, unujrad, nnukrad, nnujrad, &
& e_rad, unurad, elec_rad, nnurad

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

USE radial_ray_module, ONLY : imin, imax, ncycle, nprint, nse_c, rho_c, &
& t_c, ye_c, u_c, v_c, w_c, x_ef, dx_cf, psi0_c, xn_c, a_nuc_rep_c, &
& z_nuc_rep_c, be_nuc_rep_c, uburn_c, e_nu_c_bar, f_nu_e_bar

! added by zf2
USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

character(len=128), INTENT(in)   :: ndump           ! unit number to write restart file
INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

INTEGER                          :: io_count
LOGICAL                          :: io_initialized = .FALSE.
INTEGER                          :: error            
INTEGER                          :: adios_err       ! ADIOS error flag 

integer*8 :: io_type, handle
#define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b,adios_err)

!-----------------------------------------------------------------------
!        Formats: Document the dump
!-----------------------------------------------------------------------

 1001 FORMAT (' ***Model dump written at cycle      ',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!
!               \\\\\ WRITE THE MODEL RESTART FILES /////
!
!-----------------------------------------------------------------------

IF (io_initialized == .FALSE.) THEN
  io_count = 0
  io_initialized = .TRUE.
END IF

! added by zf: sync all procs before io
    CALL MPI_Barrier( MPI_COMM_WORLD, error )

! open start
CALL open_start(ncycle, io_count)

!CALL adios_get_group (io_type, 'restart.model'//char(0))
!CALL adios_open (handle, io_type, trim(ndump)//char(0))

CALL adios_open (handle, 'restart.model'//char(0), trim(ndump)//char(0), 'w'//char(0),adios_err)

! open end
CALL open_end(ncycle, io_count)

! write start
CALL write_start(ncycle, io_count)

ADIOS_WRITE(handle,nx)
ADIOS_WRITE(handle,nx+1)
ADIOS_WRITE(handle,nez)
ADIOS_WRITE(handle,nnu)
ADIOS_WRITE(handle,nnc)
ADIOS_WRITE(handle,ij_ray_dim)
ADIOS_WRITE(handle,ik_ray_dim)

!-----------------------------------------------------------------------
!  Radial index bounds for MGFLD shifted radial arrays
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,imin)
ADIOS_WRITE(handle,imax)

!-----------------------------------------------------------------------
!  Independent thermodynamic variables
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,rho_c)
ADIOS_WRITE(handle,t_c)
ADIOS_WRITE(handle,ye_c)

!-----------------------------------------------------------------------
!  Independent mechanical variables
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,u_c)
ADIOS_WRITE(handle,v_c)
ADIOS_WRITE(handle,w_c)
ADIOS_WRITE(handle,dx_cf)
ADIOS_WRITE(handle,x_ef)

!-----------------------------------------------------------------------
!  Independent radiation variables and bookkeeping arrays
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,psi0_c)
ADIOS_WRITE(handle,dnurad)
ADIOS_WRITE(handle,unukrad)
ADIOS_WRITE(handle,unujrad)
ADIOS_WRITE(handle,e_rad)
ADIOS_WRITE(handle,unurad)
ADIOS_WRITE(handle,nnukrad)
ADIOS_WRITE(handle,nnujrad)
ADIOS_WRITE(handle,elec_rad)
ADIOS_WRITE(handle,nnurad)
ADIOS_WRITE(handle,e_nu_c_bar)
ADIOS_WRITE(handle,f_nu_e_bar)

!-----------------------------------------------------------------------
!  Net number of neutrinos radiated from density rho_nurad and radius
!   r_nurad
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,nu_r)
ADIOS_WRITE(handle,nu_rt)
ADIOS_WRITE(handle,nu_rho)
ADIOS_WRITE(handle,nu_rhot)

!-----------------------------------------------------------------------
!  Time integrated psi0 and psi1
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,psi0dat)
ADIOS_WRITE(handle,psi1dat)

!-----------------------------------------------------------------------
!  nse - non-bse flag
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,nse_c)

!-----------------------------------------------------------------------
!  Nuclear abundances
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,xn_c)

!-----------------------------------------------------------------------
!  Auxiliary heavy nucleus
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,a_nuc_rep_c)
ADIOS_WRITE(handle,z_nuc_rep_c)
ADIOS_WRITE(handle,be_nuc_rep_c)

!-----------------------------------------------------------------------
!  Nuclear energy released
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,uburn_c)

!-----------------------------------------------------------------------
!  Energy offsets
!-----------------------------------------------------------------------

ADIOS_WRITE(handle,duesrc)

! write end
CALL write_end(ncycle, io_count)

! close start
CALL close_start(ncycle, io_count)

CALL adios_close (handle, adios_err)

! cloe end
CALL close_end(ncycle, io_count)

io_count = io_count+1

!-----------------------------------------------------------------------
!
!                    \\\\\ RECORD THE DUMP /////
!
!-----------------------------------------------------------------------

WRITE (nlog,1001) ncycle,ndump

RETURN
END SUBROUTINE model_write_adios
