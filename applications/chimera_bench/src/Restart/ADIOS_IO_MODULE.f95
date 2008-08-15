!-----------------------------------------------------------------------
!    File:         IO_MODULE
!    Module:       io_module
!    Type:         Module
!    Author:       Reuben D Budiardja, Dept of Physics, UTK,
!                  Knoxville, TN 37996
!
!    Date:         03/21/08
!
!    Purpose:
!      Contains the subroutines necessary for IO
!
!-----------------------------------------------------------------------

MODULE adios_io_module

  USE kind_module, ONLY : double
  USE array_module, ONLY: n_proc
  USE numerical_module, ONLY : zero, third, half, one, epsilon, frpi, frpith
  USE physcnst_module, ONLY : ergmev, rmu

  USE MPI
  
  USE parallel_module, ONLY : myid, myid_y, myid_z, &
      MPI_COMM_ROW, MPI_COMM_COL
  USE edit_module, ONLY : nlog, nu_r, nu_rt, nu_rho, nu_rhot, psi0dat, & 
      psi1dat, data_path
  USE eos_snc_x_module, ONLY : duesrc, aesv, nse_e=>nse
  USE mdl_cnfg_module, ONLY : jr_min, jr_max
  USE nu_dist_module, ONLY : dnurad, unukrad, unujrad, nnukrad, nnujrad, &
      e_rad, unurad, elec_rad, nnurad, dudt_nu, dunujeadt
  USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
  USE nucbrn_module, ONLY: dudt_nuc, nse_n=>nse
  USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, ncycle, &
      time, nprint, nse_c, rho_c, t_c, ye_c, rhobar, u_c, v_c, w_c, &
      x_ei, x_ef, y_ef, z_ef, dx_ci, dx_cf, dy_cf, dz_cf, x_ci, x_cf, y_cf, z_cf, &
      psi0_c, psi1_e, xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, uburn_c, &
      e_nu_c_bar, f_nu_e_bar, e_nu_c, unu_c, dunu_c, unue_e, dunue_e, &
      grav_x_c, grav_y_c, grav_z_c, gtot_pot_c, agr_e, agr_c
  USE shock_module, ONLY : pq_x, pqy_x, j_shk_radial_p
  
  IMPLICIT none
  
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: filename = 'RadHyd3D_output_'
  
  
  LOGICAL, PRIVATE                    :: io_initialized = .FALSE.
  INTEGER, PRIVATE                    :: nx              ! x-array extent
  INTEGER, PRIVATE                    :: ny              ! y-array extent globaly
  INTEGER, PRIVATE                    :: nz              ! z-array extent globaly
  INTEGER, PRIVATE                    :: nez             ! Number of Neutrino Energy-Zone 
  INTEGER, PRIVATE                    :: nnu             ! neutrino flavor array extent
  INTEGER, PRIVATE                    :: nnc             ! Composition array extent

  
integer*8, private :: io_type, adios_handle,adios_groupsize, adios_totalsize, adios_err, adios_buf_size


  !-----------------------------------------------------------------------
  !        Variables related to the domain decomposition 
  !
  !  my_j_ray_dim : Number of ray owned by local processor in j direction
  !  my_k_ray_dim : Number of ray owned by local processor in k direction
  !  my_j_ray     : Index of ray in j dimension local to a processor
  !  my_k_ray     : Index of ray in k dimension local to a processor
  !  j_ray        : Global index of j ray
  !  k_ray        : Global index of k ray
  !  j_ray_min    : Min of global index of j ray for this processor
  !  k_ray_min    : Min of global index of k ray for this processor
  !  j_ray_max    : Max of global index of j ray for this processor
  !  k_ray_max    : Max of global index of k ray for this processor
  !  nproc_y      : Number of processors assigned to the y dimension
  !  nproc_z      : Number of processors assigned to the z dimension
  !  io_count     : Counter for the number of calls to model_write_hdf5
  !-----------------------------------------------------------------------
  
  INTEGER, PRIVATE                    :: my_j_ray_dim
  INTEGER, PRIVATE                    :: my_k_ray_dim
  INTEGER, PRIVATE                    :: my_j_ray
  INTEGER, PRIVATE                    :: my_k_ray
  INTEGER, PRIVATE                    :: j_ray
  INTEGER, PRIVATE                    :: k_ray
  INTEGER, PRIVATE                    :: j_ray_min
  INTEGER, PRIVATE                    :: k_ray_min
  INTEGER, PRIVATE                    :: j_ray_max
  INTEGER, PRIVATE                    :: k_ray_max
  INTEGER, PRIVATE                    :: nproc_y
  INTEGER, PRIVATE                    :: nproc_z
  INTEGER, PRIVATE                    :: io_count
  
  !--------------------------------------------------------------------------
  !  nz_hyperslabs    : Number of separate hyperslabs output (ie. files output) 
  !                     the whole domain is divided into. Set this to 1 to get
  !                     one file output per dump, set to "nz" to get nz files 
  !                     per dump.
  !                     nz / nz_hyperslabs should be an even number
  ! nz_hyperslab_width: nz / nz_hyperslabs (computed, do not set)
  ! k_hyperslab_min   : Min of index k in the hyperslab group
  !--------------------------------------------------------------------------
  INTEGER, PRIVATE, PARAMETER         :: nz_hyperslabs = 1
  
  INTEGER, PRIVATE                    :: nz_hyperslab_width
  INTEGER, PRIVATE                    :: my_hyperslab_group
  INTEGER, PRIVATE                    :: k_hyperslab_min
  INTEGER, PRIVATE                    :: mpi_comm_per_hyperslab_group
  LOGICAL, PRIVATE                    :: hyperslab_group_master
  
  
    
  REAL(kind=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: pMD
  REAL(kind=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: sMD
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:),   PRIVATE :: r_shock      ! radius of shock maximum
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:),   PRIVATE :: r_shock_mn   ! minimum estimateed shock radius
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:),   PRIVATE :: r_shock_mx   ! maximum estimateed shock radius
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: rsphere_mean ! mean neutrinosphere radius
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: dsphere_mean ! mean neutrinosphere density
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: tsphere_mean ! mean neutrinosphere temperature
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: msphere_mean ! mean neutrinosphere enclosed mass
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: esphere_mean ! mean neutrinosphere energy
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE :: r_gain       ! gain radius
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:),   PRIVATE :: tau_adv      ! advection time scale (s)
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:),   PRIVATE :: tau_heat_nu  ! neutrino heating time scale (s)
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:),   PRIVATE :: tau_heat_nuc ! nuclear heating time scale (s)
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:),   PRIVATE :: r_nse        ! radius of NSE-nonNSE boundary
  
  REAL(KIND=double), PRIVATE                                :: io_walltime  

  
  PUBLIC      :: model_write_hdf5_adios
  PUBLIC      :: model_read_hdf5_adios
  
  PRIVATE     :: compute_plot_variables
  PRIVATE     :: mean
  

  CONTAINS

  
  SUBROUTINE initialized_io()
    
    INTEGER :: error
    INTEGER :: nproc_per_hyperslab_group
    INTEGER :: mpi_world_group
    INTEGER :: mpi_group_per_hyperslab_group
    INTEGER :: iproc
    INTEGER, DIMENSION(:), ALLOCATABLE :: hyperslab_group_member
    
    
    !------------------------------------------------------------------------
    !       Initialize Variables related to the decomposition
    !------------------------------------------------------------------------
    
    !-- FIXME: Check if this is the best way to do it
  
    if(.NOT.io_initialized)THEN
      nx = size(x_cf)
      ny = size(y_cf)
      nz = size(z_cf)
      nez = size(psi0_c, dim=2)
      nnu = size(psi0_c, dim=3)
      nnc = size(xn_c, dim=2)

      CALL MPI_COMM_SIZE(MPI_COMM_ROW, nproc_y, error)
      CALL MPI_COMM_SIZE(MPI_COMM_COL, nproc_z, error)
      my_j_ray_dim = ny/nproc_y
      my_k_ray_dim = nz/nproc_z
      
      j_ray_min = MOD(myid, nproc_y) * my_j_ray_dim + 1
      k_ray_min = (myid/nproc_y) * my_k_ray_dim + 1
      j_ray_max = MOD(myid, nproc_y) * my_j_ray_dim + my_j_ray_dim
      k_ray_max = (myid/nproc_y) * my_k_ray_dim + my_k_ray_dim
      
      ALLOCATE(pMD(nx,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(sMD(nx,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(r_shock(my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(r_shock_mn(my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(r_shock_mx(my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(rsphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(dsphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(tsphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(msphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(esphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(r_gain(nnu+1,my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(tau_adv(my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(tau_heat_nu(my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(tau_heat_nuc(my_j_ray_dim,my_k_ray_dim))
      ALLOCATE(r_nse(my_j_ray_dim,my_k_ray_dim))
      
      !-- Sanity check for nz_hyperslabs
      nz_hyperslab_width = nz / nz_hyperslabs
      IF(mod(nz, nz_hyperslabs) /= 0)THEN
        PRINT*, 'nz should be evenly divisible by nz_hyperslabs'
        PRINT*, 'nz = ', nz, ' nz_hyperslabs = ', nz_hyperslabs
        CALL MPI_ABORT(MPI_COMM_WORLD, error)
      END IF
      IF(nz_hyperslabs > nproc_z)THEN
        PRINT*, 'nz_hyperslabs cannot be more than nproc_z'
        PRINT*, 'nproc_z = ', nproc_z, ' nz_hyperslabs = ', nz_hyperslabs
        CALL MPI_ABORT(MPI_COMM_WORLD, error)
      END IF
      
      !-- Set the "master" (ie. first processor) of each hyperslab group
      my_hyperslab_group = myid_z / nz_hyperslab_width + 1
      hyperslab_group_master = .false.
      nproc_per_hyperslab_group = nproc_y * nproc_z / nz_hyperslabs
      IF((my_hyperslab_group-1) * nproc_per_hyperslab_group == myid) &
        hyperslab_group_master = .true.
      
      !-- Set the k index in the hyperslab group
      k_hyperslab_min = mod(k_ray_min-1, nz_hyperslab_width) + 1
      
      !-- Create MPI communicator for each group of writers
      ALLOCATE(hyperslab_group_member(nproc_per_hyperslab_group))
      hyperslab_group_member &
        = (/(iproc, iproc = (my_hyperslab_group-1)*nproc_per_hyperslab_group, &
                            my_hyperslab_group*nproc_per_hyperslab_group)/)
    
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD, mpi_world_group, error)
      CALL MPI_GROUP_INCL(mpi_world_group, nproc_per_hyperslab_group, &
             hyperslab_group_member, mpi_group_per_hyperslab_group, error)
      call MPI_COMM_CREATE(MPI_COMM_WORLD, mpi_group_per_hyperslab_group, &
             mpi_comm_per_hyperslab_group, error)
      call MPI_GROUP_FREE(mpi_group_per_hyperslab_group, error)
      deALLOCATE(hyperslab_group_member)

      io_walltime = zero
      io_count    = 0

      io_initialized = .TRUE.
    END IF
  
  END SUBROUTINE initialized_io


  SUBROUTINE model_write_hdf5_adios()
    !-----------------------------------------------------------------------
    !
    !    File:         model_write_hdf5
    !    Type:         Subprogram
    !    Author:       Reuben Budiardja, Dept of Physics, UTK
    !
    !    Date:         3/20/08
    !
    !    Purpose:
    !      To dump the model configuration in HDF5 file.
    !
    !    Subprograms called:
    !        write_ray_hyperslab
    !        write_1d_slab
    !
    !    Input arguments:
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


    !-----------------------------------------------------------------------
    !       File, group, dataset, and dataspace Identifier 
    !-----------------------------------------------------------------------
    !CHARACTER(LEN=29)                :: suffix
    CHARACTER(LEN=35)                :: suffix
    
    REAL(KIND=double)                :: io_startime
    REAL(KIND=double)                :: io_endtime
    
! added by zf
    INTEGER, dimension(3)            :: array_dimensions
    INTEGER, dimension(2)            :: radial_index_bound
    INTEGER, dimension(2)            :: theta_index_bound
    INTEGER, dimension(2)            :: phi_index_bound

    CALL initialized_io()    
    
    !------------------------------------------------------------------------
    !       Compute Derived Physical Values for Plotting
    !------------------------------------------------------------------------
    
    pMD(imin:imax,:,:) = aesv(imin+1:imax+1,1,:,:)
    sMD(imin:imax,:,:) = aesv(imin+1:imax+1,3,:,:)
    CALL compute_plot_variables(imin, imax, nx, nez, nnu, jmin, jmax, my_j_ray_dim, &
               ny, kmin, kmax, my_k_ray_dim, nz, x_ef, x_cf, y_cf, z_cf, u_c, &
               v_c, w_c, rho_c, t_c, ye_c, rhobar, psi0_c, psi1_e, &
               e_nu_c, unu_c, dunu_c, unue_e, dunue_e, nse_c, gtot_pot_c, &
               r_shock, r_shock_mn, r_shock_mx, &
               rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, &
               esphere_mean, r_gain, r_nse, tau_adv, tau_heat_nu, &
               tau_heat_nuc)
    
! added by zf
    array_dimensions = (/nx,ny,nz/)
    radial_index_bound = (/imin,imax/)
    theta_index_bound = (/jmin,jmax/)
    phi_index_bound = (/kmin,kmax/)

    !------------------------------------------------------------------------
    !       Create and Initialize File using Default Properties
    !------------------------------------------------------------------------

    WRITE(suffix, fmt='(i9.9,a7,i4.4,a3)') ncycle,'_group_',my_hyperslab_group,'.bp'
    !WRITE(suffix, fmt='(i9.9,a7,i4.4,a1,i6.6,a3)') ncycle,'_group_',my_hyperslab_group,'_',myid,'.bp'

#ifdef ADIOS_PROFILING 
! added by zf: sync all procs before io
    CALL MPI_Barrier( MPI_COMM_WORLD, error )
#endif

    io_startime = MPI_WTIME()

#if ADIOS_PROFILING 
! open start
    CALL open_start(ncycle, io_count)
#endif

    CALL adios_open (adios_handle, 'restart.model'//char(0), TRIM(data_path)//'/Restart/'//filename//TRIM(suffix)//char(0), 'w'//char(0),adios_err)

#if ADIOS_PROFILING
! open end
    CALL open_end(ncycle, io_count)
#endif

#if ADIOS_PROFILING
! write start
    CALL write_start(ncycle, io_count)
#endif

    #include "gwrite_restart.model.fh"
    WRITE(*,*)"Open ok:",myid

    CALL adios_close (adios_handle, adios_err)
WRITE(*,*)"adios close ok ",myid

#if ADIOS_PROFILING
! cloe end
    CALL close_end(ncycle, io_count)
#endif

    io_endtime = MPI_WTIME()
    io_walltime = io_walltime + (io_endtime-io_startime)
    io_count = io_count+1
    
    IF(myid == 0)THEN
      WRITE(nlog, 'a30,i9,a10,i5,a20,f10.5')'*** HDF5 Model dump at cycle ', ncycle, &
        'IO Count:', io_count, 'IO elapsed time: ', io_walltime 
    ENDIF
    
    RETURN

  END SUBROUTINE model_write_hdf5_adios
  
  
  SUBROUTINE model_read_hdf5_adios(read_cycle, nuc_number, i_nuc_data, nprint, &
&  nlog, i_model_data )
    !-----------------------------------------------------------------------
    !
    !    File:         model_read_hdf5
    !    Type:         Subprogram
    !    Author:       Reuben Budiardja, Dept of Physics, UTK
    !
    !    Date:         3/20/08
    !
    !    Purpose:
    !      To read the model configuration from HDF5 file.
    !
    !    Subprograms called:
    !        read_ray_hyperslab
    !        read_1d_slab
    !
    !    Input arguments:
    !  read_cycle  : the restart file cycle number to read 
    !
    !    Output arguments:
    !        none
    !
    !    Include files:
    !  edit_modulee, eos_snc_x_module, nu_dist_module, nu_energy_grid_module,
    !  radial_ray_module
    !
    !-----------------------------------------------------------------------
    
    INTEGER, INTENT(in)              :: read_cycle
    INTEGER, INTENT(in)              :: nprint
    INTEGER, INTENT(in)              :: nlog
    INTEGER, INTENT(in)              :: nuc_number    ! number of nuclear species (not counting representative heavy nucleus)
    INTEGER, INTENT(in), DIMENSION(:,:,:) :: i_nuc_data ! integer array of edit keys
    
    !-----------------------------------------------------------------------
    !        Output variables.
    !-----------------------------------------------------------------------

    INTEGER, INTENT(out), DIMENSION(2) :: i_model_data  ! integer array of initial model data

    !-----------------------------------------------------------------------
    !       File, group, dataset, and dataspace Identifier 
    !-----------------------------------------------------------------------
    
    CHARACTER(LEN=35)                :: suffix
    INTEGER                          :: imin_read       
    INTEGER                          :: imax_read
    INTEGER                          :: n
    INTEGER                          :: iadjst
    INTEGER, dimension(2)            :: radial_index_bound
    INTEGER                          :: ij_ray        ! radial ray j index
    INTEGER                          :: ik_ray        ! radial ray k index
    INTEGER                          :: j,k
    INTEGER                          :: fileUnit
    INTEGER                          :: error 

    REAL(KIND=double)                :: xn_tot        ! sum of themass fractions in a zone
    
    !-----------------------------------------------------------------------
    !        Formats
    !-----------------------------------------------------------------------

     3001 FORMAT (' jr_max + 1 =',i10,' > nx =',i4)
     3003 FORMAT (' jr_max,n_proc =',2i6,' MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 ) /= 0')

    
    CALL initialized_io()
    
    x_ei  = zero
    dx_ci = zero
    
    rho_c = zero
    t_c   = zero
    ye_c  = zero
    u_c   = zero
    v_c   = zero
    w_c   = zero
    
    psi0_c       = zero
    nse_c        = 0
    xn_c         = zero
    be_nuc_rep_c = zero
    a_nuc_rep_c  = zero
    z_nuc_rep_c  = zero
    uburn_c      = zero
    
    
    WRITE(suffix, fmt='(i9.9,a7,i4.4,a3)') read_cycle,'_group_',my_hyperslab_group,'.bp'
    !WRITE(suffix, fmt='(i9.9,a7,i4.4,a1,i6.6,a3)') ncycle,'_group_',my_hyperslab_group,'_',myid,'.bp'

    if(hyperslab_group_master)THEN
      PRINT*, '***Read Model Dump', TRIM(data_path)//'/Restart/'//filename//suffix
      WRITE(nlog, 'a50'), '***Read Model Dump', TRIM(data_path)//'/Restart/'//filename//suffix
    END IF

    ! open bp file for read
    CALL adios_open (adios_handle, 'restart.model'//char(0), TRIM(data_path)//'/Restart/'//filename//TRIM(suffix)//char(0), 'r'//char(0),error)
    
    if(error /= 0)then
      PRINT*, '***ERROR in trying to open ', TRIM(data_path)//'/Restart/'//filename//suffix
      WRITE(nlog, 'a50'), '***ERROR in trying to open ', TRIM(data_path)//'/Restart/'//filename//suffix
      CALL MPI_ABORT(MPI_COMM_WORLD, error)
    end if

    #include "gread_restart.model.fh"

    CALL adios_close(adios_handle,adios_err)

    imin_read = radial_index_bound(1)
    imax_read = radial_index_bound(2)

    
    !-----------------------------------------------------------------------
    !     Set values not read directly
    !-----------------------------------------------------------------------
    x_ci(1:nx) = x_ei(1:nx) + half * dx_ci(1:nx)
    nse_e(2:nx,:,:) = nse_c(1:nx-1,:,:)
    nse_n(2:nx)     = nse_c(1:nx-1,1,1)
    
    !-----------------------------------------------------------------------
    !
    !             \\\\\ CHECK DATA EXTENTS FOR CONSISTENCY /////
    !
    !        jr_max must be compatible with n_proc, so the y_arrays and
    !         z-arrays can be distributed uniformly over the processors.
    !
    !-----------------------------------------------------------------------
    jr_min                 = imin_read + 1
    jr_max                 = imax_read + 1
    
    i_model_data(1)        = imin_read
    i_model_data(2)        = imax_read

    IF ( jr_max + 1 > nx ) THEN
      WRITE (nlog,3001) jr_max+1,nx
      WRITE (nprint,3001) jr_max+1,nx
      STOP
    END IF ! jr_max + 1 > nx

    IF ( MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 ) /= 0 ) THEN
      WRITE (nlog,3003) jr_max,n_proc
      WRITE (nprint,3003) jr_max,n_proc
      STOP
    END IF ! MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 )

    !-----------------------------------------------------------------------
    !
    !                \\\\\ ADJUST NUCLEAR ABUNDANCES /////
    !
    !-----------------------------------------------------------------------

    iadjst                 = i_nuc_data(4,1,1)

    IF ( iadjst == 1 ) THEN
      DO ik_ray = 1,my_k_ray_dim
        DO ij_ray = 1,my_j_ray_dim
          DO j = 1,nx
            IF ( nse_c(j,ij_ray,ik_ray) == 0 ) THEN
              xn_tot         = zero
              DO n = 1,nuc_number+1
                xn_tot       = xn_tot + xn_c(j,n,ij_ray,ik_ray)
              END DO
              DO n = 1,nuc_number+1
                xn_c(j,n,ij_ray,ik_ray) = xn_c(j,n,ij_ray,ik_ray)/( xn_tot + epsilon )
              END DO
            END IF ! nse(j,ij_ray,ik_ray) == 0
          END DO ! j
        END DO ! ij_ray
      END DO ! ik_ray
    END IF ! iadjst == 1
    
    RETURN

  END SUBROUTINE model_read_hdf5_adios
  

  SUBROUTINE compute_plot_variables(imin, imax, nx, nez, nnu, jmin, jmax, ij_ray_dim, &
               ny, kmin, kmax, ik_ray_dim, nz, x_e_in, x_c_in, y_in, z_in, uMD_in, &
               vMD_in, wMD_in, rhoMD_in, tMD_in, yeMD_in, rhobar, psi0_in, psi1_in, &
               e_nu_MD_in, unu_in, dunu_in, unue_in, dunue_in, nse_in, gtot_pot_c, &
               r_shock_in, r_shock_in_mn, r_shock_in_mx, &
               rsphere_mean_in, dsphere_mean_in, tsphere_mean_in, msphere_mean_in, &
               esphere_mean_in, r_gain_in, r_nse_in, tau_adv_in, tau_heat_nu_in, &
               tau_heat_nuc_in)
    
    !    Subroutine:   compute_plot_variables
    !    Type:         Subprogram
    !
    !    Date:         03/26/08
    !
    !    Purpose:
    !      To compute variables needed for output files and plot
    !
    !    Subprograms called:
    !
    !    Input arguments:
    !
    !  imin         : minimum x-array index for the edit
    !  imax         : maximum x-array index for the edit
    !  jmin         : minimum y-array index for the edit
    !  jmax         : maximum y-array index for the edit
    !  kmin         : minimum z-array index for the edit
    !  kmax         : maximum z-array index for the edit
    !  nx           : x-array extent
    !  ny           : y-array extent
    !  nz           : z-array extent
    !  ij_ray_dim   : number of y-zones on a processor before swapping with y
    !  ik_ray_dim   : number of z-zones on a processor before swapping with z
    !  nez          : neutrino energy array extent
    !  nnu          : neutrino flavor array extent
    !  x_e_in       : radial edge of zone (cm)
    !  x_c_in       : radial midpoint of zone (cm)
    !  y_in         : y (angular) midpoint of zone
    !  z_in         : z (azimuthal) midpoint of zone
    !  uMD_in       : x (radial) velocity of zone (cm s^{-1})
    !  vMD_in       : y (angular) velocity of zone (cm s^{-1})
    !  wMD_in       : z (azimuthal) velocity of zone (cm s^{-1})
    !  rhoMD_in     : density of zone (g cm^{-3})
    !  tMD_in       : temperature of zone (K)
    !  yeMD_in      : electron fraction of zone
    !  rhobar       : mean density at a given radius (g cm^{-3})
    !  psi0_in      : zero moment of the neutrino distribution
    !  psi1_in      : first moment of the neutrino distribution
    !  e_nu_MD_in   : neutrino energy density (ergs cm^{-3})
    !  unu_in       : radial zone-centered neutrino energy
    !  dunu_in      : radial zone-centered neutrino energy zone width
    !  unue_in      : radial zone-edged neutrino energy
    !  dunue_in     : radial zone-edged neutrino energy zone width
    !  nse_in       : NSE flag
    !  gtot_pot_c   : unshifted zone-centered gravitational potential energy [ergs]
    
    !-----------------------------------------------------------------------
    !        Input variables
    !-----------------------------------------------------------------------
    
    INTEGER, INTENT(in)                        :: imin            ! minimum x-array index for the edit
    INTEGER, INTENT(in)                        :: imax            ! maximum x-array index for the edit
    INTEGER, INTENT(in)                        :: jmin            ! minimum y-array index for the edit
    INTEGER, INTENT(in)                        :: jmax            ! number of radial rays assigned to a processor
    INTEGER, INTENT(in)                        :: kmin            ! minimum y-array index for the edit
    INTEGER, INTENT(in)                        :: kmax            ! maximum z-array index for the edit
    INTEGER, INTENT(in)                        :: nx              ! x-array extent
    INTEGER, INTENT(in)                        :: ny              ! y_array extent
    INTEGER, INTENT(in)                        :: nz              ! z_array extent
    INTEGER, INTENT(in)                        :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
    INTEGER, INTENT(in)                        :: ik_ray_dim      ! number of z-zones on a processor before swapping with z
    INTEGER, INTENT(in)                        :: nez             ! neutrino energy array extent
    INTEGER, INTENT(in)                        :: nnu             ! neutrino flavor array extent
    INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: nse_in ! NSE flag
    REAL(KIND=double), INTENT(in), DIMENSION(nx+1)  :: x_e_in          ! radial midpoint of zone (cm)
    REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: x_c_in          ! radial edge of zone (cm)
    REAL(KIND=double), INTENT(in), DIMENSION(ny)    :: y_in            ! y (angular) midpoint of zone
    REAL(KIND=double), INTENT(in), DIMENSION(nz)    :: z_in            ! z (azimuthal) midpoint of zone
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: uMD_in   ! x (radial) velocity of zone (cm s^{-1})
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: vMD_in   ! y (angular) velocity of zone (cm s^{-1})
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: wMD_in   ! z (azimuthal) velocity of zone (cm s^{-1})
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: rhoMD_in ! density of zone (g cm^{-3})
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: tMD_in   ! temperature of zone (K)
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: yeMD_in  ! entropy of zone
    REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rhobar          ! mean density at a given radius (g cm^{-3})
    REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_in  ! zero moment of the neutrino distribution
    REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1_in  ! first moment of the neutrino distribution
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: e_nu_MD_in ! neutrino energy density (ergs cm^{-3})
    REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: unu_in       ! radial zone-centered neutrino energy
    REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: dunu_in      ! radial zone-centered neutrino energy zone width
    REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: unue_in      ! radial zone-edged neutrino energy
    REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: dunue_in     ! radial zone-edged neutrino energy zone width
    REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: gtot_pot_c   ! unshifted zone-centered gravitational potential energy (ergs g^{-1}

    !-----------------------------------------------------------------------
    !        Output variables
    !-----------------------------------------------------------------------
    REAL(KIND=double), INTENT(OUT), DIMENSION(ij_ray_dim,ik_ray_dim) :: r_shock_in      ! radius of shock maximum
    REAL(KIND=double), INTENT(OUT), DIMENSION(ij_ray_dim,ik_ray_dim) :: r_shock_in_mn   ! minimum estimateed shock radius
    REAL(KIND=double), INTENT(OUT), DIMENSION(ij_ray_dim,ik_ray_dim) :: r_shock_in_mx   ! maximum estimateed shock radius
    REAL(KIND=double), INTENT(OUT), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)   :: rsphere_mean_in ! mean neutrinosphere radius
    REAL(KIND=double), INTENT(OUT), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)   :: dsphere_mean_in ! mean neutrinosphere density
    REAL(KIND=double), INTENT(OUT), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)   :: tsphere_mean_in ! mean neutrinosphere temperature
    REAL(KIND=double), INTENT(OUT), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)   :: msphere_mean_in ! mean neutrinosphere enclosed mass
    REAL(KIND=double), INTENT(OUT), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)   :: esphere_mean_in ! mean neutrinosphere energy
    REAL(KIND=double), INTENT(OUT), DIMENSION(nnu+1,ij_ray_dim,ik_ray_dim) :: r_gain_in       ! gain radius
    REAL(KIND=double), INTENT(OUT), DIMENSION(ij_ray_dim,ik_ray_dim)       :: tau_adv_in      ! advection time scale (s)
    REAL(KIND=double), INTENT(OUT), DIMENSION(ij_ray_dim,ik_ray_dim)       :: tau_heat_nu_in  ! neutrino heating time scale (s)
    REAL(KIND=double), INTENT(OUT), DIMENSION(ij_ray_dim,ik_ray_dim)       :: tau_heat_nuc_in ! nuclear heating time scale (s)
    REAL(KIND=double), INTENT(OUT), DIMENSION(ij_ray_dim,ik_ray_dim)       :: r_nse_in  ! radius of NSE-nonNSE boundary

    !-----------------------------------------------------------------------
    !        Local variables
    !-----------------------------------------------------------------------

    INTEGER                                         :: i               ! radial zone index
    INTEGER                                         :: j               ! y (angular) zone index
    INTEGER                                         :: k               ! z (azimuthal) zone index
    INTEGER                                         :: n               ! neutrino energy index
    INTEGER                                         :: jr_min          ! shifted miminim radial zone index
    INTEGER                                         :: jr_max          ! shifted maximum radial zone index
    INTEGER                                         :: jr_maxp         ! jr_max + 1
    INTEGER                                         :: j_shock         ! shifted radial zone index of shock maximum
    INTEGER                                         :: j_shock_mn      ! shifted radial zone index of minimum estimated shock radius
    INTEGER                                         :: j_shock_mx      ! shifted radial zone index of maximum estimated shock radius
    INTEGER                                         :: i_shock         ! radial zone index of shock maximum
    INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)       :: i_shock_in      ! radial zone index of shock maximum
    INTEGER                                         :: i_shockm        ! radial zone index of shock next to maximum
    INTEGER                                         :: i_shock_mn      ! radial zone index of minimum estimated shock radius
    INTEGER                                         :: i_shock_mx      ! radial zone index of maximum estimated shock radius
    INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)       :: j_shock_in      ! shifted radial zone index of shock maximum
    INTEGER, DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: j_sphere        ! zone index of the neutrinospheres
    INTEGER, DIMENSION(nnu,ij_ray_dim,ik_ray_dim)   :: jsphere_mean_in ! shifted radial zone <= mean n-neutrinosphere
    INTEGER                                         :: j_gain          ! shifted radial index of gain radius
    INTEGER                                         :: i_gain          ! unshifted radial index of gain radius
    INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)       :: i_gain_in       ! index of the gain radius
    INTEGER                                         :: nwt             ! standard deviation flag for linear fitting

    
    REAL(KIND=double), DIMENSION(nx)                :: dr_c            ! radial zone thickness (cm)
    REAL(KIND=double), DIMENSION(nx)                :: dvol            ! radial zone volume (cm^{3})
    REAL(KIND=double), DIMENSION(nx)                :: m_neut_e        ! unshifted zone_edged Newtonian mass (g)
    REAL(KIND=double), DIMENSION(nx)                :: m_neut_c        ! unshifted zone_centered Newtonian mass (g)
    
    REAL(KIND=double), PARAMETER                    :: UTOT0 = 8.9d0   ! change in the zero of energy (MeV)
    REAL(KIND=double), PARAMETER                    :: ku = ergmev/rmu ! ( # nucleons/gram )( erg/mev )
    REAL(KIND=double), PARAMETER                    :: e_bar0   = ku * UTOT0

    REAL(KIND=double), PARAMETER                    :: pqmin = 1.d0    ! shock criterion
    REAL(KIND=double)                               :: q_shock_in1     ! strength of shock maximum
    REAL(KIND=double)                               :: q_shock_in2     ! strength of shock second to maximum
    REAL(KIND=double)                               :: r_shock_in1     ! radius of shock maximum
    REAL(KIND=double)                               :: r_shock_in2     ! radius of shock second to maximum
    REAL(KIND=double)                               :: r_shk           ! radius of shock maximum
    REAL(KIND=double)                               :: m_shk           ! mass enclosed by shock maximum
    
    REAL(KIND=double), DIMENSION(nx)                :: rho             ! shifted density (cm^{-3})
    REAL(KIND=double), DIMENSION(nx)                :: t               ! shifted temperature (K)
    REAL(KIND=double), DIMENSION(nx)                :: rstmss          ! enclosed rest mass (g)
    
    REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: lum_in ! neutrino luminosity (foes)
    REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: e_rms_stat_in   ! SQRT( SUM psi0 * w5dw/SUM w3ww ) (MeV)
    REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: e_rms_trns_in   ! SQRT( SUM psi1 * w5dw/SUM w3ww ) (MeV)
    REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: r_sphere        ! neutrinosphere radius
    REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: d_sphere        ! neutrinosphere density
    REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: t_sphere        ! neutrinosphere temperature
    REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: m_sphere        ! neutrinosphere enclosed mass
    REAL(KIND=double)                               :: e               ! energy per unit mass (ergs g^{-1})
    REAL(KIND=double)                               :: E_env           ! total binding energy of heating region
    REAL(KIND=double)                               :: Q_nu            ! total neutrino heating rate in heating region
    REAL(KIND=double)                               :: Q_nuc           ! total nuclear heating rate in heating region
    REAL(KIND=double), DIMENSION(nx)                :: sig             ! set of standard deviations for linear fits
    REAL(KIND=double)                               :: a_x             ! parameter a_x in the straight line fit y = a_x x + a_c
    REAL(KIND=double)                               :: a_c             ! parameter a_c in the straight line fit y = a_x x + a_c

    !-----------------------------------------------------------------------
    !  MGFLD shifted indices
    !-----------------------------------------------------------------------

    jr_min                    = imin + 1
    jr_max                    = imax + 1
    jr_maxp                   = jr_max + 1

    !-----------------------------------------------------------------------
    !  Mean enclosed mass
    !-----------------------------------------------------------------------

    dr_c(1:imax)              = x_e_in(2:imax+1) - x_e_in(1:imax)
    dvol(1:imax)              = frpi * dr_c(1:imax) * ( x_e_in(1:imax) * ( x_e_in(1:imax) + dr_c(1:imax) ) &
    &                         + dr_c(1:imax) * dr_c(1:imax) * third ) 

    i                         = imin
    IF ( x_e_in(i) == zero ) THEN
      m_neut_e(i)             = zero
    ELSE
      m_neut_e(i)             = rhobar(i) * frpith * x_e_in(i)**3
    END IF

    DO  i = imin,imax
      m_neut_e(i+1)           = m_neut_e(i) + rhobar(i) * dvol(i)
    END  DO

    m_neut_c(imin:imax)       = half * ( m_neut_e(imin:imax ) + m_neut_e(imin+1:imax+1) )
    
    !-----------------------------------------------------------------------
    !  Shock location
    !-----------------------------------------------------------------------

    DO k = 1,ik_ray_dim
      DO j = 1,ij_ray_dim
        CALL findshock( imin+1, imax, j, k, pqmin, j_shock, j_shock_mx, j_shock_mn, &
    &    nx, x_c_in, m_neut_c, r_shk, m_shk )
        i_shock                 = MAX( j_shk_radial_p(j,k) - 1, 1 )
        i_shock_in(j,k)         = i_shock
        i_shock_mx              = MAX( j_shock_mx - 1, 1 )
        i_shock_mn              = MAX( j_shock_mn - 1, 1 )
        j_shock_in(j,k)         = i_shock + 1
        r_shock_in1             = x_c_in(i_shock)
        r_shock_in_mn(j,k)      = x_c_in(i_shock_mn)
        r_shock_in_mx(j,k)      = x_c_in(i_shock_mx)
        IF ( i_shock_mn == 1 ) r_shock_in_mn(j,k) = 0.d0
        IF ( i_shock_mx == 1 ) r_shock_in_mx(j,k) = 0.d0
        IF ( i_shock    == 1 ) THEN
          r_shock_in (j,k)      = 0.d0
        ELSE
          IF ( ( pq_x(j_shock+1,j,k) + pqy_x(j_shock+1,j,k) )/( aesv(j_shock+1,1,j,k) + epsilon )  &
    &        > ( pq_x(j_shock-1,j,k) + pqy_x(j_shock-1,j,k) )/( aesv(j_shock-1,1,j,k) + epsilon ) ) THEN
            i_shockm            = i_shock + 1
          ELSE
            i_shockm            = i_shock - 1
          END IF ! pq_x/aesv
          r_shock_in2           = x_c_in(i_shockm)
          q_shock_in1           = ( pq_x(i_shock+1 ,j,k) + pqy_x(i_shock+1 ,j,k) )/( aesv(i_shock+1 ,1,j,k) + epsilon )
          q_shock_in2           = ( pq_x(i_shockm+1,j,k) + pqy_x(i_shockm+1,j,k) )/( aesv(i_shockm+1,1,j,k) + epsilon )
          r_shock_in(j,k)       = mean( q_shock_in1, q_shock_in2, r_shock_in1, r_shock_in2 )
        END IF ! i_shock   -1 /= 1
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim

    !-----------------------------------------------------------------------
    !  Neutrino luminosities
    !-----------------------------------------------------------------------

    DO k = 1,ik_ray_dim
      DO j = 1,ij_ray_dim
        CALL luminosity_MD( imin, imax, j, k, ij_ray_dim, ik_ray_dim, nx, nez, &
    &    nnu, x_e_in, psi1_in, unue_in, dunue_in, lum_in )
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim

    !-----------------------------------------------------------------------
    !  Neutrino rms energies
    !-----------------------------------------------------------------------

    DO k = 1,ik_ray_dim
      DO j = 1,ij_ray_dim
        CALL e_rms_MD( imin, imax, j, k, ij_ray_dim, ik_ray_dim, nx, nez, nnu, &
    &    psi0_in, psi1_in, unu_in, dunu_in, e_rms_stat_in, e_rms_trns_in )
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim

    !-----------------------------------------------------------------------
    !  Mean neutrinospheres
    !-----------------------------------------------------------------------

    rstmss(1)                 = zero
    DO i = imin+1,imax+1
      rstmss(i)               = rstmss(i-1) + frpith * ( x_e_in(i)**3 - x_e_in(i-1)**3 ) * rhobar(i-1)
    END DO

    DO k = 1,ik_ray_dim
      DO j = 1,ij_ray_dim

        rho(jr_min:jr_max)    = rhoMD_in(imin:imax,j,k)
        t  (jr_min:jr_max)    = tMD_in  (imin:imax,j,k)

        CALL nu_sphere( jr_min, jr_maxp, j, k, ij_ray_dim, ik_ray_dim, x_e_in, &
    &    rho, t, rstmss, nx, nez, nnu, j_sphere, r_sphere, d_sphere, t_sphere, &
    &    m_sphere )

        CALL nu_sphere_mean( jr_min, jr_maxp, j, k, ij_ray_dim, ik_ray_dim,    &
    &    e_rms_trns_in, j_sphere, r_sphere, d_sphere, t_sphere, m_sphere, nx,  &
    &    nez, nnu, rsphere_mean_in, dsphere_mean_in, tsphere_mean_in,          &
    &    msphere_mean_in, esphere_mean_in, jsphere_mean_in )

      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim

    !-----------------------------------------------------------------------
    !  Gain radius
    !
    !   Determine r_gain(1:2,i_ray) by working inward from the shock and
    !    locating the radius where dunujeadt changes sign
    !
    !   Determine r_gain(5,i_ray)  by working inward from the shock and
    !    locating the radius where two adjecent dudt_nu are positive
    !    followed inward by two adjecent negative dudt_nu
    !
    !   Set gain radius to zero if either
    !    (1) neutrinosphere == 0
    !    (2) no shock
    !-----------------------------------------------------------------------

    DO n = 1,2
      DO j = 1,ij_ray_dim
        DO k = 1,ik_ray_dim
          j_gain              = jr_min
          i_gain              = jr_min - 1
          IF ( jsphere_mean_in(n,j,k) == jr_min  .or.  r_shock_in(j,k) == zero ) THEN
            r_gain_in(n,j,k)  = zero
          ELSE
            DO i = jr_max,jsphere_mean_in(n,j,k)+1,-1
              IF ( nse_in(I+1,j,k) == 0 ) CYCLE
              IF ( dunujeadt(i,n,j,k) > zero  .and.  dunujeadt(i+1,n,j,k) < zero ) THEN
                j_gain        = i
                i_gain        = i - 1
                EXIT
              END IF ! dunujeadt(i,n,j,k) > zero  .and.  dunujeadt(i+1,n,j,k) < zero
            END DO ! i = jr_max,jsphere_mean_in(n,j,k)+1,-1
            IF ( j_gain == jr_min ) THEN
              r_gain_in(n,j,k)  = zero
              CYCLE
            END IF ! j_gain == jr_min
            r_gain_in(n,j,k)  = mean( dunujeadt(j_gain,n,j,k), -dunujeadt(j_gain+1,n,j,k), &
    &                                 x_c_in(i_gain+1), x_c_in(i_gain) )
          END IF ! jsphere_mean_in(n,j) == jr_min  .or.  r_shock_in(j) == zero
        END DO ! k = 1,ik_ray_dim
      END DO ! j = 1,ij_ray_dim
    END DO ! n = 1,2

    DO j = 1,ij_ray_dim
      DO k = 1,ik_ray_dim
        j_gain                = jr_min
        i_gain                = jr_min - 1
        i_gain_in(j,k)        = i_gain
        IF ( jsphere_mean_in(1,j,k) == jr_min  .or.                        &
    &        jsphere_mean_in(2,j,k) == jr_min  .or.                        &
    &        r_shock_in(j,k) == zero )                                      THEN
          r_gain_in(5,j,k)    = zero
        ELSE
          DO i = jr_max,jsphere_mean_in(2,j,k)+1,-1
            IF ( nse_in(I+1,j,k) == 0 ) CYCLE
            IF ( dudt_nu(i-1,j,k) < zero  .and.                            &
    &            dudt_nu(i,j,k) < zero    .and.                            &
    &            dudt_nu(i+1,j,k) > zero  .and.                            &
    &            dudt_nu(i+2,j,k) > zero )                                  THEN
              j_gain          = i
              i_gain          = i - 1
              i_gain_in(j,k)  = i_gain
              EXIT
            END IF ! dudt_nu(i,j,k) > zero  .and.  dudt_nu(i+1,j,k) < zero 
          END DO ! i = jr_max,jsphere_mean_in(2,j,k)+1,-1
          IF ( j_gain == jr_min ) THEN
            r_gain_in(5,j,k)  = zero
            CYCLE
          ELSE
            r_gain_in(5,j,k)  = mean( -dudt_nu(j_gain,j,k), dudt_nu(j_gain+1,j,k), &
    &                                  x_c_in(i_gain+1), x_c_in(i_gain) )
          END IF ! j_gain == jr_min
        END IF ! jsphere_mean_in(1,j,k) == jr_min  .or.  jsphere_mean_in(2,j,k) == jr_min  .or.  r_shock_in(j,k) == zero
      END DO ! k = 1,ik_ray_dim
    END DO ! j = 1,ij_ray_dim
    
    !-----------------------------------------------------------------------
    !  Advection time scale
    !
    !   Determine the advection time scale by integrating dr/(-u(r)) from
    !    the neutrinosphere to the shock
    !-----------------------------------------------------------------------

    DO j = 1,ij_ray_dim
      DO k = 1,ik_ray_dim
        IF ( jsphere_mean_in(1,j,k) == jr_min              .or.             &
    &        jsphere_mean_in(2,j,k) == jr_min              .or.             &
    &        r_shock_in(j,k)        == zero                .or.             &
    &        i_shock_in(j,k) - 2    <= i_gain_in(j,k) + 2  .or.             &
    &        i_gain_in(j,k)         == 1 )                   THEN
          tau_adv_in(j,k)     = 1.d+100
        ELSE
          nwt                 = 0
          sig                 = 1.d0
          CALL linear_fit( x_c_in, uMD_in(:,j,k), nx, i_gain_in(j,k), i_shock_in(j,k)-2, sig, nwt, a_c, a_x )
          tau_adv_in(j,k)     = zero
          DO i = i_shock_in(j,k) - 2, i_gain_in(j,k) + 1, -1
            tau_adv_in(j,k)   = tau_adv_in(j,k) + dr_c(i)/( DABS( a_x * x_c_in(i) + a_c ) + epsilon )
          END DO ! i = i_shock_in(j,k),i_gain_in(j,k),-1
        END IF ! jsphere_mean_in(1,j,k) == jr_min. etc.
      END DO ! k = 1,ik_ray_dim
    END DO ! j = 1,ij_ray_dim

    !-----------------------------------------------------------------------
    !  Neutrino heating time scale
    !
    !   Determine the heating time scale by integrating E/dEdt from
    !    the neutrinosphere to the shock
    !-----------------------------------------------------------------------

    DO j = 1,ij_ray_dim
      DO k = 1,ik_ray_dim
        IF ( tau_adv_in(j,k) == 1.d+100 ) THEN
          tau_heat_nu_in(j,k) = 1.d+100
          tau_heat_nuc_in(j,k) = 1.d+100
        ELSE
          E_env               = zero
          Q_nu                = zero
          Q_nuc               = zero
          DO i = i_shock_in(j,k) - 2, i_gain_in(j,k) + 1, -1
            e                 = half * ( uMD_in(i,j,k) * uMD_in(i,j,k)      &
    &                         + vMD_in(i,j,k) * vMD_in(i,j,k) )             &
    &                         + aesv(i+1,2,j,k) - e_bar0 + gtot_pot_c(i,j,k)
            E_env             = E_env + e * ( m_neut_e(i+1) - m_neut_e(i) )
            Q_nu              = Q_nu  + dudt_nu(i+1,j,k)  * ( m_neut_e(i+1) - m_neut_e(i) )
            Q_nuc             = Q_nuc + dudt_nuc(i+1,j,k) * ( m_neut_e(i+1) - m_neut_e(i) )
          END DO ! i = i_shock_in(j,k),i_gain_in(j,k)+1,-1
          tau_heat_nu_in(j,k) = E_env/( Q_nu + epsilon )
          tau_heat_nuc_in(j,k) = E_env/( Q_nuc + epsilon )
        END IF ! tau_adv_in(j,k) == 1.d+100
      END DO ! k = 1,ik_ray_dim
    END DO ! j = 1,ij_ray_dim
    
    
    !-----------------------------------------------------------------------
    !  NSE-nonNSE boundary
    !-----------------------------------------------------------------------

      DO j = 1,ij_ray_dim
        DO k = 1,ik_ray_dim
          DO i = imin,imax
            IF ( nse_in(i,j,k) == 0 ) THEN
              r_nse_in(j,k) = x_e_in(i)
              EXIT
            END IF ! nse_in(i,j,k) == 0
          END DO ! i = imin,imax
        END DO ! k = kmin,kmax
      END DO ! j = jmin,jmax


    
  END SUBROUTINE compute_plot_variables
  
  
  FUNCTION mean(q1,q2,x1,x2)
  
    REAL (KIND=double) :: mean
    REAL (KIND=double) :: q1
    REAL (KIND=double) :: q2
    REAL (KIND=double) :: x1
    REAL (KIND=double) :: x2

    mean                 = ( q1 * x1 + q2 * x2 )/( q1 + q2 )

  END FUNCTION mean

END MODULE adios_io_module
