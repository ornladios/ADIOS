!! Simulation parameter and control module
module sml_module
  integer :: sml_mype !! processor index
  integer :: sml_totalpe !! total number of processors
  integer :: sml_mstep  !! total number of time steps for simulation
  integer :: sml_istep  !! current time step number
  integer :: sml_ipc    !! 2nd order Runge-Kutta index, or predictor-corrector index
  integer :: sml_stdout !! NOT USING NOW. standard output number 
  integer :: sml_iptl   !! NOT USING NOW. implicit arguments for debugging
  real (kind=8) :: sml_time !! simulation time 
  real (kind=8) :: sml_dt   !! time step size

  integer :: sml_deltaf  !! delta-f switch. 0 for off, 1 for on. delta-f simulation is not verified yet, especially for output file.
  integer :: sml_dwdt_exb_only !! delta-f weight evolution switch. 0 for whole (grad_b, curl_b,and ExB), 1 for exb only
  logical :: sml_dwdt_fix_bg  !! delta-f weight evolution switch. false for (1-w) factor true for (1) factor. default is .false.
  integer :: sml_no_para_nonlinear !! parallel nonlinearity switch 0 for usual full simulation 1 for no nonlinearity
  integer :: sml_canonical_maxwell  !! NOT IMPLEMENTED YET. Set initial or f0 distribution to canonical Maxwellian. 0 for off, 1 for on.
  real (kind=8) :: sml_marker_temp_factor !! Apply virtual temperature for initial loading. The virtual temperature will be (temp_factor)*(real temperature).
  
  integer :: sml_concentric !! Simulation for concentric circular geometry
  integer :: sml_read       !! Read B-field profile from a file
  !   real (kind=8), parameter :: sml_2pi=6.28318530717959D0  !! 2 Pi
  !  real (kind=8), parameter :: sml_pi =3.14159265358979D0 !! Pi
  real (kind=8) :: sml_2pi !! 2 Pi
  real (kind=8) :: sml_pi  !!  Pi
  real (kind=8) :: sml_sqrtpi  !! Sqrt(Pi)
  real (kind=8) :: sml_sqrt2pi !! Sqrt(2Pi)
  
  real (kind=8), parameter :: sml_e_charge=1.6022D-19  !! electron charge (MKS)
  real (kind=8), parameter :: sml_epsilon0=8.8542D-12  !! permittivity of free space (MKS)
  real (kind=8), parameter :: sml_prot_mass=1.6720D-27 !! proton mass (MKS)
!  real (kind=8), parameter :: sml_elec_mass=9.1094D-31 !! electron mass (MKS)
  real (kind=8), parameter :: sml_elec_mass=0.01D0*sml_prot_mass !! electron mass (MKS)

  real (kind=8), parameter :: sml_ev2j=sml_e_charge, sml_j2ev=1D0/sml_e_charge ! e, 1/e
 
  integer :: sml_bounce         !! Bounce routine switch  0 for off, 1 for inner boundary, 2 for both boundaries
  integer :: sml_restart        !! Restarting simulation switch.  0 for original run, 1 for restart
  integer :: sml_minusB         !! Change the direction of whole B-field
  integer :: sml_bt_sign        !! Change the direction of toroidal B-field only
  integer :: sml_monte_num      !! Number of sample for Monte-Carlo volume calculation. It should be greater than particle number for correct simulation.
  real (kind=8) :: sml_bd_min_r, sml_bd_max_r,sml_bd_min_z,sml_bd_max_z !! simulation boundary in R-Z space
  real (kind=8) :: sml_marker_den !! Marker density. Used for loading and initial weight calculation.
  real (kind=8) :: sml_inpsi, sml_outpsi  !! Inner and outer boundary for initial loading. 
!  integer :: sml_pot_diag_period
  integer :: sml_restart_write_period !! Number of steps between particle data dumps for restart.

  real (kind=8) :: sml_en_order_kev  !! Characteristic ion energy (keV)
  real (kind=8) :: sml_en_order      !! Characteristic ion energy in normalized unit
  real (kind=8) :: sml_tran !! Torodial transit time of an ion with the characteristic energy. 
  
  integer :: sml_machine  !! Machine type 
  integer :: sml_limiter  !! Switch for limiter
  integer :: sml_special  !! Switch for special simulation. eg. single particle simulation.

  real (kind=8) :: sml_bp_mult !! Artificial multiplication of poloidal B-field
  real (kind=8) :: sml_bt_mult !! Artificial multiplication of toroidal B-field

  integer :: sml_angle_stored  !! Flag for real angle storage for each particle. 1 for valid data, 0 for invalid data.

  !for gyrokinetic
  integer, parameter  :: sml_nlarmor=4
  integer :: sml_nphi_total
  integer :: sml_turb_efield=0
  integer :: sml_no_00_efield=0

  !for multi speces
!  real (kind=8) ::   sml_c2_2m_rel(2), sml_c_m_rel(2),sml_mass_rel(2), sml_charge_rel(2)

  !for ExB supress
  integer :: sml_exb_suppress
  real (kind=8) :: sml_exb_on(2)=(/1D0, 0D0/), sml_exb_suppress_time, sml_exb_suppress_width
  
  ! for initial distribution relaxation
  integer :: sml_relax_dist
  integer :: sml_relax_dist_num

  ! for delta-f simulation f0 function
  real (kind=8) :: sml_f0den_edge, sml_f0den_out, sml_f0den_ped_c,sml_f0den_ped_width
  integer :: sml_deltaf_f0_mode
  real (kind=8), allocatable :: sml_f0_n(:), sml_f0_gradn(:,:)
  real (kind=8) :: sml_f0_1_Ln, sml_f0_1_Lt  ! for sml_deltaf_f0_mode==-1  1/Ln , 1/Lt - artificial f0

  ! grid file input
  character (len=65) :: sml_node_file, sml_ele_file, sml_bfollow_file
  character (len=65) :: sml_diag_node_file, sml_diag_ele_file

  ! time averaging potential
  real (kind=8) :: sml_tavg_factor
  integer :: sml_tavg_on
  
  !multi PE for one plane
  integer :: sml_plane_per_pe, sml_pe_per_plane
#ifdef ADIOS
  integer :: sml_pe_group,sml_pe_mype,sml_pe_comm
#endif
  integer :: sml_plane_group, sml_plane_comm, sml_plane_totalpe, sml_plane_mype,sml_plane_index
  ! electron on/off switch
  integer :: sml_electron_on
  
  integer :: sml_fem_matrix, sml_smtrc_ad_elec

  !for p-c.
  real (kind=8), allocatable :: sml_pc_coeff(:,:)
  integer :: sml_push_mode !! 1 for RK2, 2 for P-C
  integer :: sml_pc_order
  
  ! grid guess table size
  integer :: sml_guess_table_size

  ! delta-f loading parameter
  real (kind=8) :: sml_initial_deltaf_noise

  ! zero inner boundary of 00 solver
  integer :: sml_zero_inner_bd

  ! toroidal mode select for linear GK simulation
  integer :: sml_mode_select_on, sml_mode_select_n
  integer :: sml_mode_initial_n

  !
  real (kind=8) :: sml_bd_ext_delta1, sml_bd_ext_delta2, sml_bd_ext_delta3, sml_bd_ext_delta4

  integer :: sml_bfollow, sml_bfollow_read
  ! input file directory
  character (len=65) :: sml_input_file_dir

  integer :: sml_bounce_zero_weight

  logical :: sml_use_pade=.true.
  logical :: sml_use_simple00=.false.
  integer :: sml_simple00_nsmth=0

  integer :: sml_add_pot0  ! additional potential 0 for no pot, 1 for read-in, 2 for neoclassical pot0 (simple)
  character (len=256) :: sml_add_pot0_file
  integer :: sml_flat_electron_density ! smoothing of background electron density for ion only simulation.   0 for no smoothing (default) 1 for smoothing.
  logical :: sml_supress_weight_growth  ! if this value is true then weight of deltaf will be between -sml_weight_max and sml_weight_max. default value is .false.
  real (kind=8) :: sml_weight_max   ! maximum weight(= phase(6)) for delta-f simulation. default value is 10.

  ! for restore temperature gradient
  logical :: sml_restore_temp  ! restore temperature gradient. Work for delta-f=1, sml_deltaf_f0_mode=-1, and rk2 method only
  integer :: sml_restore_temp_period ! how frequently restore temperature
  real (kind=8) :: sml_restore_temp_avg_scale ! time average scale - dt base, default=0.01 (means 100 dt scale)
  integer :: sml_max_mat_width

  ! for cold electron at the boundary
  integer :: sml_bd_Te_mode    ! 0 for no effect, 1 for innner only, 2 for both boundary, 3 for outer only
  real (kind=8) :: sml_bd_Te_width !radial decay length of electron temperature at the boundary 

  logical :: sml_zero_out_total_charge ! total charge is zero always
  
  real (kind=8),parameter :: sml_boundary_diagonal=1D-6  ! diagonal value of boundary node point

  logical :: sml_exclude_private

  integer :: sml_sheath_mode
  real (kind=8) :: sml_sheath_init_pot_factor

  logical :: sml_rgn1_pot0_only

  logical :: sml_heat_on
end module


!! Equilibrium module
module eq_module
  character (len=80):: eq_header, eq_filename
  real (kind=8) :: eq_x_psi, eq_x_z, eq_x_r, eq_min_r, eq_max_r, eq_min_z, eq_max_z,eq_axis_r,eq_axis_z,eq_axis_b
  integer :: eq_mr, eq_mz, eq_mpsi,eq_sep
  real (kind=8), allocatable :: eq_I(:), eq_psi_grid(:),eq_rgrid(:), eq_zgrid(:),eq_psi_rz(:,:)
  
  !not in file
!! Initial ion temperature profile variable.  Tanh profile.
!!    - edge : inside separatrix
!!    - out : outside separatrix
!!    - ped_c : pedestal center 
!!    - ped_width : pedestal width.
  real (kind=8) :: eq_tempi_ev_edge, eq_tempi_ev_out,eq_tempi_ped_c,eq_tempi_ped_width  
!! Initial electron temperature profile variable.  Tanh profile.
!!    - edge : inside separatrix
!!    - out : outside separatrix
!!    - ped_c : pedestal center 
!!    - ped_width : pedestal width.
  real (kind=8) :: eq_tempe_ev_edge, eq_tempe_ev_out, eq_tempe_ped_c, eq_tempe_ped_width
!! Initial ion density profile variable.  Tanh profile.
!!    - edge : inside separatrix
!!    - out : outside separatrix
!!    - ped_c : pedestal center 
!!    - ped_width : pedestal width.
  real (kind=8) :: eq_den_edge, eq_den_out, eq_den_ped_c, eq_den_ped_width


  integer, parameter :: eq_mid_r_npsi=50  !!  For poloidal flux -> midplane R function. Number of evaluation points 
  real (kind=8) :: eq_mid_r_psi(eq_mid_R_npsi) !! For poloidal flux -> midplane R function. Function data
  real (kind=8) :: eq_mid_r_dp  !! For poloidal flux -> midplane R function. Delta psi 
end module


!! Particle module
module ptl_module
  type species_type
     integer :: num  !! number of particle for each processor
     integer :: maxnum
     integer :: type !! species type number 1 for ion, 2 for electron
     real (kind=8) :: mass   !! Ion mass (Atomic Unit)
     real (kind=8) :: charge !! Ion charge (electron charge unit)
     real (kind=8) :: c_m     
     integer :: nphase=9 !! Number of phase variable + 3(weight)
     integer :: nphase2=18 !! 2*ptl_nphase
#ifdef SMALL_GID
     integer,  pointer :: gid(:) !! Global ID number of each particle. >0 for inside simulation region, <0 for outside simulation region.
     integer  ::  maxgid !! Largest global ID number used.
#else
     integer (kind=8), pointer :: gid(:) !! Global ID number of each particle. >0 for inside simulation region, <0 for outside simulation region.
     integer (kind=8) ::  maxgid !! Largest global ID number used.
#endif

!     integer, pointer :: reset_derivs(:)
     real (kind=8), pointer :: phase(:,:) !! Ion phase variable - phase(var,ptl)
     real (kind=8), pointer :: phase0(:,:) !! Ion phase variable - storage for old phse variable  
     real (kind=8), pointer :: derivs(:,:,:) !! Storage for time derivative of phase variable (time history) derivs(var,time,ptl)
     ! variables for ion only
     real (kind=8), pointer :: angle(:)
     real (kind=8), pointer :: weight0(:)
     ! lost particle save
     integer :: lost_num, lost_nummax
     real (kind=8), pointer :: lost_index(:)
     ! particle-grid position save
     integer, pointer :: tr_save(:,:,:)  !size of electron and ion are different
     real (kind=8), pointer :: p_save(:,:,:,:) ! size of electron and ion are different
     real (kind=8), pointer :: pos_save(:,:)
     integer , pointer :: pc_index(:)
     integer :: stored_derivs
#ifdef CIRCULAR_SPECIAL
     integer, pointer :: rindex(:,:,:), tindex(:,:,:,:)
     real (kind=8), pointer :: rweight(:,:,:), tweight(:,:,:,:)
#endif

  end type species_type
  type ptl_type
     type(species_type) :: elec, ion


  end type ptl_type
  real (kind=8) :: ptl_mass(2), ptl_charge(2),  ptl_c_m(2), ptl_c2_2m(2)
  integer :: ptl_rz_outside, ptl_nphase_max=0
  !special simulation purpose : such as single particle simulation
  real (kind=8) :: ptl_special_r, ptl_special_z,  ptl_special_phi, ptl_special_en_ev, ptl_special_pitch
  integer :: ptl_special_n
  integer :: ptl_maxnum, ptl_e_maxnum,ptl_num, ptl_e_num, ptl_lost_nummax
  real (kind=8) :: ptl_mass_au, ptl_charge_eu, ptl_e_mass_au, ptl_e_charge_eu

#ifdef COREVV2
  real (kind=8), allocatable :: ptl_exb(:,:)
#endif

  contains
    subroutine ptl_mem_allocation( sp, sp_type, maxnum,mass, charge, pc_order, nlarmor, lost_nummax)
      integer, intent(in) :: sp_type, pc_order, maxnum, nlarmor, lost_nummax
      real (kind=8), intent(in) :: mass, charge
      type(species_type) :: sp
      integer :: np

      sp%maxnum=maxnum
      np=sp%nphase
      sp%mass=mass
      sp%charge=charge
      sp%c_m=charge/mass

      ptl_nphase_max=max(ptl_nphase_max,sp%nphase)

      sp%lost_nummax=lost_nummax
      sp%lost_num=0 

      allocate( sp%phase(np, maxnum), sp%phase0(np,maxnum), &
           sp%gid(maxnum),sp%derivs(np,pc_order,maxnum), &
           sp%lost_index(lost_nummax) )
           
      
      if(sp_type==1) then !for ion 
         allocate( sp%angle(maxnum), sp%weight0(maxnum) )
         allocate( sp%tr_save(nlarmor,2,maxnum),sp%p_save(3,nlarmor,2,maxnum) )
         sp%type=1
#ifdef COREVV2
         allocate(ptl_exb(2,maxnum))
#endif
#ifdef CIRCULAR_SPECIAL
         allocate( sp%rindex(nlarmor,2,maxnum), sp%tindex(2,nlarmor,2,maxnum), &
              sp%rweight(nlarmor,2,maxnum ), sp%tweight(2,nlarmor,2,maxnum) )
#endif
      else
         allocate( sp%tr_save(1,1,maxnum),sp%p_save(3,1,1,maxnum) )
#ifdef CIRCULAR_SPECIAL
         allocate( sp%rindex(1,2,maxnum), sp%tindex(2,1,2,maxnum), &
              sp%rweight(1,2,maxnum), sp%tweight(2,1,2,maxnum) )
#endif
         sp%type=2
      endif

      allocate(sp%pos_save(5,maxnum))

      allocate(sp%pc_index(pc_order))
      do i=1,pc_order
         sp%pc_index(i)=i
      enddo
      sp%stored_derivs=0
!      sp%reset_derivs=0
    end subroutine ptl_mem_allocation
end module ptl_module


!! Module for Field value
Module fld_module
  type fld_type
    real (kind=8) :: r, z, phi  !! NOT USED. R,Z coordinate variable
    real (kind=8) :: I         !! R*Bt
    real (kind=8) :: epot      !! Electric potential
    real (kind=8) :: br, bz, bphi  ! B_r, B_z, B_phi
    real (kind=8) :: dbrdr, dbrdz, dbrdp, dbzdr, dbzdz, & 
         dbzdp, dbpdr, dbpdz, dbpdp !! dB_x / dx , x=R,Z,phi
    real (kind=8)   :: Er, Ez, Ephi  !! E_r, E_z, E_phi
    !for weight calculation
    real (kind=8) :: dIdpsi, dpsidr,dpsidz !! dI/dpsi, dpsi/dr, dpsi/dz for weight calculation
    real (kind=8) :: psi  !! Poloidal Flux
    !for weight calculation : f0 information
    real (kind=8) :: f0_den,f0_gradn(2),f0_temp,f0_gradt(2)
 end type fld_type
end module


!! Interpolation module
module itp_module
#if defined(PSPLINE)
  use EZspline_obj
  type(EZspline2_r8) :: spl(0:2,0:2)
  type(EZspline1_r8) :: spl_psi
#endif


  integer :: itp_mpsi, itp_mr, itp_mz
  integer, parameter :: itp_mr2=100, itp_mz2=100
  integer,parameter :: itp_korder=3, itp_korder_rz=5, itp_korder_rz2=5
  real (kind=8) :: itp_min_r, itp_max_r, itp_min_z, itp_max_z
  real (kind=8) :: itp_min_psi, itp_max_psi
  
  real (kind=8), allocatable :: itp_I_cscoef(:,:), itp_I_break(:)
  real (kind=8), allocatable :: itp_psi_knot(:), itp_r_knot(:), itp_z_knot(:)
  real (kind=8), allocatable :: itp_psi_bscoef(:)

  real (kind=8) :: itp_rgrid2(itp_mr2), itp_zgrid2(itp_mz2)
  real (kind=8) :: itp_r_knot2(itp_mr2+itp_korder_rz2), itp_z_knot2(itp_mz2+itp_korder_rz2)

  real (kind=8) :: itp_psi_bscoef00(itp_mr2*itp_mz2), itp_psi_bscoef01(itp_mr2*itp_mz2), &
       itp_psi_bscoef10(itp_mr2*itp_mz2), itp_psi_bscoef02(itp_mr2*itp_mz2), &
       itp_psi_bscoef11(itp_mr2*itp_mz2), itp_psi_bscoef20(itp_mr2*itp_mz2)       
end module

#ifdef NETCDF
module netcdf_out
  use netcdf
  type netcdf1d
     ! netCDF file ID
     integer :: ncfileid
     ! netCDF dimension IDs
     integer :: psi_ncdim, time_ncdim
     ! program variables for dimesions
     integer :: time_ncstart(1), ncstart(2)
     integer :: time_nccount(1), nccount(2)
     ! netCDF variable IDs
     integer :: time_ncvarid, psi_ncvarid
     integer, allocatable :: ncvarid(:,:)
     ! netCDF variable shapes
     integer :: time_ncdims(1), psi_ncdims(1), ncdims(2)
  end type netcdf1d
  type netcdf2d
     ! netCDF file ID
     integer :: ncfileid
     ! netCDF dimension IDs
     integer :: psi_ncdim, time_ncdim
     ! program variables for dimesions
     integer :: time_ncstart(1), ncstart(2)
     integer :: time_nccount(1), nccount(2)
     ! netCDF variable IDs
     integer :: time_ncvarid, psi_ncvarid
     integer, allocatable :: ncvarid(:,:)
     ! netCDF variable shapes
     integer :: time_ncdims(1), psi_ncdims(1), ncdims(2)     
  end type netcdf2d
  type netcdf0d
     ! netCDF file ID
     integer :: ncfileid
     ! netCDF dimension IDs
     integer :: time_ncdim
     ! program variaables for dimensions
     integer :: time_ncstart(1), ncstart(2)
     integer :: time_nccount(1), nccount(2)
     ! netCDF variable IDs
     integer :: time_ncvarid
     integer, allocatable :: ncvarid(:)
     ! netCDF variable shpes
     integer :: time_ncdims(1), ncdims(1)
  end type netcdf0d
  contains
    subroutine netcdf1d_finalize(data1d)
      implicit none
      type(netcdf1d) :: data1d
      ! netCDF error status return
      integer :: iret

      ! sete netCDF attribute for ElVis to stop monitoring netCDF file
      iret = nf90_redef(data1d%ncfileid)
      call check_err(iret)
      iret = nf90_put_att(data1d%ncfileid, NF90_GLOBAL, 'running', 'false')
      call check_err(iret)

      ! leave netCDF define mode and close file
      iret = nf90_enddef(data1d%ncfileid)
      call check_err(iret)
      iret = nf90_close(data1d%ncfileid)
      call check_err(iret)
      
      deallocate(data1d%ncvarid)
    end subroutine netcdf1d_finalize     
end module netcdf_out
#endif

!! Diagnosis module
module diag_module
#ifdef NETCDF
  use netcdf_out
#endif
  integer, parameter :: diag_max_sp_num=2  ! two species

  integer :: diag_efld_period  !! Period for efiled diagnosis
  integer :: diag_flow_period  !! Period for flow diagnosis

  integer :: diag_tracer_period !! Period for tracer
  integer :: diag_tracer_n     !! Particle index for tracer routine
  integer :: diag_tracer_sp     ! Species index for tracer 
  
  integer :: diag_pw_period  !! Period for particle weight diagnosis
  integer :: diag_pw_on      !! on-off switch for particle weight diagnosis


  ! for flow
  !! Number of data point (shell) for flow diagnosis
  integer, parameter :: diag_flow_npsi=50 , diag_flow_nvars1=7, diag_flow_nvars2=diag_flow_nvars1*2+11 
  real (kind=8) :: diag_flow_vol(diag_flow_npsi) !! Volume for each flow diagnosis shell
  real (kind=8) :: diag_flow_pin,diag_flow_pout,diag_flow_dp !! Inner boundary psi, Outer boundary psi, delta psi for flow diagnosis

  integer :: diag_pot_period      !! frequency of ASCII dumps of pot data
  integer :: diag_binout_period   !! frequency of binpack dumps of 2d data

  !avg diagnosis - flux 
  integer :: diag_avg_on
  integer :: diag_avg_outperiod !! Write-out period for time averging diagnosis
  integer, parameter :: diag_avg_npsi=50 !! Number of data point (shell) for time averaing diagnosis. 
  real (kind=8) :: diag_avg_pin, diag_avg_pout, diag_avg_dp !! Inner boundary psi, outer boundary psi, delta psi for time averaging diagnosis
  integer, parameter :: diag_avg_nv1=10, diag_avg_nv2=5 
  real (kind=8) :: diag_avg_vol(diag_avg_npsi)
  real (kind=8) :: diag_avg_flux(diag_avg_npsi,diag_avg_nv1,diag_max_sp_num)

  !diagnosis - f
  integer,parameter :: diag_f_nv=30, diag_f_npsi=20
  integer :: diag_f_on,diag_f_count, diag_f_mod, diag_f_skip
  real (kind=8) :: diag_f_vmax, diag_f_pin, diag_f_pout, diag_f_dpsi, diag_f_angle1, diag_f_angle2, diag_f_slope1, diag_f_slope2
  real (kind=8) ,dimension(diag_f_nv,-diag_f_nv:diag_f_nv,diag_f_npsi) :: diag_f, diag_f0,diag_f_c,diag_f0_c
  real (kind=8) , allocatable :: diag_f_n0(:), diag_f_vt(:) , diag_f_vcs(:),diag_f_dvol(:)

  !diagnosis GAM
  integer :: diag_gam_node_begin, diag_gam_on, diag_gam_node_num

  !diagnosis particle dump
  integer :: diag_ptl_on, diag_ptl_begin, diag_ptl_end, diag_ptl_num
  integer (kind=8) :: diag_ptl_ifile, diag_ptl_efile
  real (kind=8), allocatable :: diag_ptl_times(:), diag_ptl_data1(:), &
       diag_ptl_data2(:)

#ifdef NETCDF


  ! netCDF file ID for flow diagnostics --- diag_flow
  integer :: diag_flow_ncfileid
  ! netCDF dimension IDs
  integer :: diag_psi_ncdim, diag_time_ncdim
  ! program variables for dimensions
  integer :: diag_time_ncstart(1), diag_flow_ncstart(2)
  integer :: diag_time_nccount(1), diag_flow_nccount(2)
  ! netCDF variable IDs
  integer :: diag_time_ncvarid, diag_psi_ncvarid
  integer :: diag_flow_ncvarid(diag_flow_nvars2,diag_max_sp_num) ! this number should be NN2 in flow_diagnosis
  ! netCDF variable shapes
  integer :: diag_time_ncdims(1), diag_psi_ncdims(1), diag_flow_ncdims(2)
  !
  ! netCDF file ID for flow diagnostics --- diag_avg

  integer :: diag_avg_ncfileid
  ! netCDF dimension IDs
  integer :: diag_avg_psi_ncdim, diag_avg_time_ncdim
  ! program variables for dimensions
  integer :: diag_avg_time_ncstart(1), diag_avg_ncstart(2)
  integer :: diag_avg_time_nccount(1), diag_avg_nccount(2)
  ! netCDF variable IDs
  integer :: diag_avg_time_ncvarid, diag_avg_psi_ncvarid 
  integer :: diag_avg_ncvarid(diag_avg_nv1+diag_avg_nv2,diag_max_sp_num) ! this number should be NN2 in flow_diagnosis
  ! netCDF variable shapes
  integer :: diag_avg_time_ncdims(1), diag_avg_psi_ncdims(1), diag_avg_ncdims(2)

  type(netcdf1d) :: diag_efld_nc
  
#endif

end module diag_module


!****************************************************************************
! simulation parameter for collision module definition
!
! first created : 2000/10/26
! last modified : 2002/05/15
!****************************************************************************
!! Collision module
Module col_module
  integer :: col_mode
  integer :: col_period
  integer, parameter :: col_imp=1 ! imp=1 impurity species, 0=none    

  real (kind=8)  :: col_accel  ! artificial collision amplifying factor
  
  !2002/05/15 -- col mode 2 variables
  integer,parameter :: col_2_m=30 !! Number of radial slice for conserving collision
  integer,parameter :: col_2_mtheta=8 !! Number of poloidal slice for conserving collision
  real (kind=8) :: col_2_pin, col_2_pout !! Inner and outer boundary  
  real (kind=8) :: col_2_dp      !! Delta psi of a cell for conserving collision
  real (kind=8) :: col_2_dtheta  !! Delta theta of a cell for conserving collision 
  real (kind=8), dimension(col_2_m,col_2_mtheta)  :: col_2_vol !! Volume of each cell for conserving collision
  real (kind=8), dimension(col_2_m,col_2_mtheta) :: col_2_dvp_sum, &
       col_2_dv2_sum, col_2_dw_sum, col_2_w0, col_2_org_weight_sum, col_2_org_v2_avg, col_2_org_vp_avg  !! Momentum and energy values or their change for conserving collision.
  

  !2004/01/16 -- impurity collision
  integer :: col_imp_on  !!Switch for impurity collision
  real (kind=8) :: col_imp_charge, col_imp_mass  !! Charge(electorn) and Mass(AU) of impurity
  real (kind=8) ::   col_den_imp_edge,col_den_imp_out,col_den_imp_ped_c, col_den_imp_ped_width  !! Profile information of impurity ions. Tanh model.

  
  !2002/08/30 -- varying background
  integer :: col_varying_bg !!Switch for background plasma update
  integer :: col_vb_m       !! Number of background update points
  integer :: col_vb_period  !! Period of background updating
  real (kind=8) :: col_vb_pin, col_vb_pout, col_vb_dp !! Inner and outer boundary and delta psi for background update
  real (kind=8) , allocatable :: col_vb_den(:), col_vb_temp(:),col_vb_vol(:) !! Density, Temperature and volume of each shell for backgournd update
  
  !2002/10/09
  integer :: col_en_col_on   !! Switch for energy collision
end module col_module


!****************************************************************************
! ripple parameter 
!
! first created : 2000/11/01
! last modified : 2000/11/01
!! ripple module
!! ripple strenght = ratio*B*exp(tau/tau0)
!! tau= (R-R0)**2 + elong*z**2
!****************************************************************************
Module rpl_module

    integer :: rpl_mode	!! Ripple switch
    real (kind=8) :: rpl_N_coil  !! Number of toroidal coil
    real (kind=8) :: rpl_R0 !! Eliptic ripple model. tau=(R-R0)**2 + elong*z**2
    real (kind=8) :: rpl_tau0 !!  Eliptic ripple model. increase rate.  ripple strength = ratio*B*exp(tau/tau0)
    real (kind=8) :: rpl_ratio !! Eliptic ripple model. Realative strength of ripple field to B at the center of elipse. ripple strength = ratio*B*exp(tau/tau0)
    real (kind=8) :: rpl_elong !! Eliptic ripple model. Elongation of elipise. tau=(R-R0)**2 + elong*z**2
    !r0 : center of elipse
    !elong : elongation of elipse
    !tau0 : increase rate
    !ratio : relative strength to B at the center of elipse
    !ripple model
    ! ripple strenght = ratio*B*exp(tau/tau0)
    ! tau= (R-R0)**2 + elong*z**2
	
end module

!! Radial E-field module
module efld_module
  integer :: efld_mode !! Mode selection for E-field model : 0 for no E-field, 1 for static modeled E-field, 2 for self-consistent E-field
  
  !static efield
  real (kind=8) :: efld_1_psi1, efld_1_psie, efld_1_psim, efld_1_pot0, efld_1_potm , efld_1_dpsi !! Varialble for static E-field 

  ! varying efield variables
  integer,parameter :: efld_npsi=50 !! Number of data point for varying E-field evaluation.

  real (kind=8) :: efld_pout, efld_pin,efld_dp !! Inner and outer boundary and delta psi for varying E-field
  real (kind=8) :: efld_pol_factor !! A constant parameter for radial E-field evaluation (epsilon_0/m_i)
  real (kind=8) :: efld_dpdp(efld_npsi) !! d(Epot) / dpsi
  real (kind=8) :: efld_vol(efld_npsi) !! volume of each flux shell
  real (kind=8) :: efld_d2pdpdt(efld_npsi) !! d^2(Epot) / dpsi dt
  real (kind=8) :: efld_pot(efld_npsi) !! potential
  real (kind=8) :: efld_t0 !! time at E-field evaluation 
  real (kind=8) :: efld_dpdp0_outside 

  ! efield reset
  integer :: efld_reset_step !! the E-field will be set to zero when reset_step == istep 
  integer :: efld_reset !! Switch for E_field reset
  integer :: efld_start !! Time step when self-consistent E-field calculation starts
  
  !efield cut-off
  integer ::efld_cutoff !! Switch for E-field cut-off : Set zero E-field when psi < cutoff_psi
  real (kind=8) :: efld_cutoff_psi !! Set zero E-field when psi < cutoff_psi
 
  !efield maximum 2003/10/17
  real (kind=8) :: efld_max_dpdp   !! Maximum E-field allowed. Preventing too rapid growth of E-field. 

  !initiail efield
  !real (kind=8) :: efld_init_dpdp
  integer :: efld_read_dpdp  !! Switch for reading in initial E-field

  real (kind=8) :: efld_dpdr_mid(efld_npsi) !! for diagnosis only - dpsi/dr value at the midplane as a function of psi
end module


!last modified 2002/05/18
module rf_module
  integer :: rf_on
  real (kind=8) :: rf_r_res, rf_za,rf_lap,rf_e0
  real (kind=8),allocatable :: rf_r_save(:)
end module

#if defined(IMSL)
module arg_module ! argument passing for IMSL routine
  real (kind=8) ::  arg_bnc_psi, arg_bnc_r,arg_bnc_z,arg_bnc_E,arg_bnc_L,arg_bnc_mu
end module arg_module
#endif

module neu_module
  character (LEN=64) :: neu_sepfile
  integer ::neu_col_mode, neu_adjust_n0, neu_cx_period, &
       neu_ion_period, neu_varying_mfp
  integer ::neu_col_period
  real (kind=8) ::neu_t0
  
  integer , parameter :: neu_sep_mtheta=100
  integer :: neu_sep_mtheta_file
  real (kind=8) :: neu_sep_r(neu_sep_mtheta), neu_sep_z(neu_sep_mtheta)
  real (kind=8), allocatable :: neu_sep_r_file(:), neu_sep_z_file(:)
  real (kind=8) :: neu_n0, neu_delta_n, neu_delta_theta, neu_mfp, &
       neu_theta_x, neu_recycle_rate, neu_mfp0
  real (kind=8) :: neu_old_inside_weight, neu_weight_sum_lost, &
       neu_actual_ionized, neu_ionize2_psi
  integer :: neu_ionize_mode, neu_start_time
  real (kind=8) :: neu_temp_factor

  !---------- neutral2 ----------------
  integer, parameter :: neu_grid_mpsi=50, neu_grid_mtheta=10
  integer,parameter :: neu_mtheta=120
  real (kind=8) :: neu_r0(neu_mtheta),neu_z0(neu_mtheta), &
       neu_jacob(neu_mtheta),neu_dtheta,neu_psi_edge,neu_temp0
  real (kind=8) :: neu_dt
  integer :: neu_istep_max, neu_num,neu_mpol,neu_monte_num
  real (kind=8) :: neu_grid_den(neu_grid_mpsi,neu_grid_mtheta),&
       neu_grid_temp(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_wsum(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_esum(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_vol(neu_grid_mpsi,neu_grid_mtheta)
  
  real (kind=8) :: neu_grid_dtheta, neu_grid_dpsi, &
       neu_grid_min_psi, neu_grid_max_psi
  integer :: neu_mode2_period
  
end module neu_module

module tbl_module
  integer :: tbl_diffusion_on, tbl_period, tbl_stop_time
  real (kind=8) :: tbl_D_coeff,tbl_inpsi
  integer :: tbl_mult_max,tbl_mult_period
  real (kind=8), allocatable :: tbl_mult(:)

end module

module heat_module
  integer :: heat_period, heat_mult_max, heat_mult_period
  real (kind=8) :: heat_inpsi, heat_outpsi, heat_power
  real (kind=8), allocatable :: heat_mult(:)
end module

module lim_module
  character (LEN=64) :: lim_filename  !! Limiter file name
  real (kind=8), allocatable :: lim_r(:,:), lim_z(:,:), lim_weight(:,:), lim_en(:,:)&
       ,lim_org_r(:),lim_org_z(:),lim_ds(:,:)
  real (kind=8):: lim_dz, lim_zmin, lim_zmax,lim_r0_down,lim_r0_up,lim_psi_min
  integer :: lim_zindex(2), lim_store_mz,lim_mdata
end module


module bnc_module
  integer, parameter :: bnc_nr=100
  real (kind=8) :: bnc_min_r, bnc_max_r, bnc_dr
  real (kind=8) :: bnc_z_psi_min(bnc_nr)
  real (kind=8) :: bnc_arg_pass_r
end module bnc_module

!----------------------PETSc STUFF-------------------------------------------------

module perf_monitor
  integer, parameter :: mon_NX=59
  integer, parameter :: mon_N=23    ! limit on number of PETSc log events?
  integer, parameter :: mon_N2=59
  integer, parameter :: TOTAL_           = 1
  integer, parameter :: INIT1_           = 2
  integer, parameter :: LOAD_            = 3
  integer, parameter :: INIT2_           = 4
  integer, parameter :: MAIN_LOOP_       = 5
  integer, parameter :: FINALIZE_SPEC3_  = 6
  integer, parameter :: FINALIZE_        = 7
  integer, parameter :: IPC_LOOP_        = 8
  integer, parameter :: COLL_SNAP_       = 9
  integer, parameter :: COLLISION_       =10
  integer, parameter :: NEUTRAL_         =11
  integer, parameter :: DIFFUSION_       =12
  integer, parameter :: HEATING_         =13
  integer, parameter :: RESTART_WRITE_   =14
  integer, parameter :: CHARGEI_         =15
  integer, parameter :: CHARGEE_         =16
  integer, parameter :: POISSON_         =17
  integer, parameter :: DERIVSI_         =18
  integer, parameter :: DERIVSE_         =19
  integer, parameter :: DIAGNOSIS_       =20
  integer, parameter :: PUSHI_           =21
  integer, parameter :: PUSHE_           =22
  integer, parameter :: SHIFT_           =23
  integer, parameter :: NEUTRAL_STEP_    =24
  integer, parameter :: NEUTRAL_SUM_     =25
  integer, parameter :: NEUTRAL_GATHER_  =26
  integer, parameter :: CHARGEI_ALL_PART_=27
  integer, parameter :: CHARGEI_SR_      =28
  integer, parameter :: Z_Z_MODE_EXT_    =29
  integer, parameter :: CHARGEI_ACCUM_   =30
  integer, parameter :: CHARGEE_ALL_PART_=31
  integer, parameter :: CHARGEE_DENSITY_ =32
  integer, parameter :: PHASE0_          =33
  integer, parameter :: PHASE1_          =34
  integer, parameter :: SOLVER1_         =35
  integer, parameter :: PHASE2_          =36
  integer, parameter :: SOLVER2_         =37
  integer, parameter :: PHASE3_          =38
  integer, parameter :: SOLVER3_         =39
  integer, parameter :: PHASE4_          =40
  integer, parameter :: PHASE5_          =41
  integer, parameter :: PUSH_ALL_PART_   =42
  integer, parameter :: SHIFT_IE_RED_    =43
  integer, parameter :: SHIFT_IE_SR_R_   =44
  integer, parameter :: SHIFT_IE_SR_L_   =45
  integer, parameter :: SMOOTH_          =46
  integer, parameter :: ANGLE_           =47
  integer, parameter :: COL_SNAP_RED_    =48
  integer, parameter :: CONS_COL_RED_    =49
  integer, parameter :: CONS_COL_RED2_   =50
  integer, parameter :: DIAG_F_RED_      =51
  integer, parameter :: DIAG_TOR_RED_    =52
  integer, parameter :: FLOW_DIAG_RED_   =53
  integer, parameter :: DIAG_TIME_AVG_   =54
  integer, parameter :: TIME_AVG_RED_    =55
  integer, parameter :: PWEIGHT_RED_     =56
  integer, parameter :: DIAG2D_RED_      =57
  integer, parameter :: HEAT_RED_        =58
  integer, parameter :: GET_VOL_RED_     =59

#if !defined(NO_PETSC)
  integer :: event(mon_N)
#else
  real (kind=8) :: mon_time(mon_N)
  real (kind=8) :: mon_sum(mon_N)
#endif
  character (len=21) :: mon_str (2*mon_NX)
  logical            :: mon_sync(2*mon_NX)

contains
  subroutine init_perf_monitor( ierr ) 
#if defined(CAM_TIMERS)
    use sml_module
    use perf_mod, only: t_barrier_onf, t_initf
#endif
    implicit none
    integer ierr,i
    logical masterproc
#if ( defined TIMING_BARRIERS )
    logical :: timing_barrier = .true.
#else
    logical :: timing_barrier = .false.
#endif
    character (len=21) :: ctemp1, ctemp2
    include 'mpif.h'

#if defined(CAM_TIMERS)
    if (sml_mype == 0) then
       masterproc = .true.
    else
       masterproc = .false.
    endif
    call t_initf("perf_in", LogPrint=masterproc, Mpicom=MPI_COMM_WORLD, &
                 MasterTask=masterproc)
    timing_barrier = t_barrier_onf()
#endif

    mon_str(TOTAL_)            = 'TOTAL                '
    mon_str(INIT1_)            = 'INIT1                '
    mon_str(LOAD_)             = 'LOAD                 '
    mon_str(INIT2_)            = 'INIT2                '
    mon_str(MAIN_LOOP_)        = 'MAIN_LOOP            '
    mon_str(FINALIZE_SPEC3_)   = 'FINALIZE_SPEC3       '
    mon_str(FINALIZE_)         = 'FINALIZE             '
    mon_str(IPC_LOOP_)         = 'IPC_LOOP             '
    mon_str(COLL_SNAP_)        = 'COLL_SNAP            '
    mon_str(COLLISION_)        = 'COLLISION            '
    mon_str(NEUTRAL_)          = 'NEUTRAL              '
    mon_str(DIFFUSION_)        = 'DIFFUSION            '
    mon_str(HEATING_)          = 'HEATING              '
    mon_str(RESTART_WRITE_)    = 'RESTART_WRITE        '
    mon_str(CHARGEI_)          = 'CHARGEI              '
    mon_str(CHARGEE_)          = 'CHARGEE              '
    mon_str(POISSON_)          = 'POISSON              '
    mon_str(DERIVSI_)          = 'DERIVSI              '
    mon_str(DERIVSE_)          = 'DERIVSE              '
    mon_str(DIAGNOSIS_)        = 'DIAGNOSIS            '
    mon_str(PUSHI_)            = 'PUSHI                '
    mon_str(PUSHE_)            = 'PUSHE                '
    mon_str(SHIFT_)            = 'SHIFT                '
    mon_str(NEUTRAL_STEP_)     = 'NEUTRAL_STEP         '
    mon_str(NEUTRAL_SUM_)      = 'NEUTRAL_SUM          '
    mon_str(NEUTRAL_GATHER_)   = 'NEUTRAL_GATHER       '
    mon_str(CHARGEI_ALL_PART_) = 'CHARGEI_ALL_PART     '
    mon_str(CHARGEI_SR_)       = 'CHARGEI_SENDRECV     '
    mon_str(Z_Z_MODE_EXT_)     = '0_0_MODE_EXT         '
    mon_str(CHARGEI_ACCUM_)    = 'CHARGEI_ACCUM        '
    mon_str(CHARGEE_ALL_PART_) = 'CHARGEE_ALL_PART     '
    mon_str(CHARGEE_DENSITY_)  = 'CHARGEE_DENSITY      '
    mon_str(PHASE0_)           = 'PHASE0               '
    mon_str(PHASE1_)           = 'PHASE1               '
    mon_str(SOLVER1_)          = 'SOLVER1              '
    mon_str(PHASE2_)           = 'PHASE2               '
    mon_str(SOLVER2_)          = 'SOLVER2              '
    mon_str(PHASE3_)           = 'PHASE3               '
    mon_str(SOLVER3_)          = 'SOLVER3              '
    mon_str(PHASE4_)           = 'PHASE4               '
    mon_str(PHASE5_)           = 'PHASE5               '
    mon_str(PUSH_ALL_PART_)    = 'PUSH_ALL_PART        '
    mon_str(SHIFT_IE_RED_)     = 'SHIFT_IE_RED         '
    mon_str(SHIFT_IE_SR_R_)    = 'SHIFT_IE_SR_R        '
    mon_str(SHIFT_IE_SR_L_)    = 'SHIFT_IE_SR_L        '
    mon_str(SMOOTH_)           = 'SMOOTH               '
    mon_str(ANGLE_)            = 'ANGLE                '
    mon_str(COL_SNAP_RED_)     = 'COL_SNAP_RED         '
    mon_str(CONS_COL_RED_)     = 'CONS_COL_RED         '
    mon_str(CONS_COL_RED2_)    = 'CONS_COL_RED2        '
    mon_str(DIAG_F_RED_)       = 'DIAG_F_RED           '
    mon_str(DIAG_TOR_RED_)     = 'DIAG_TOR_RED         '
    mon_str(FLOW_DIAG_RED_)    = 'FLOW_DIAG_RED        '
    mon_str(DIAG_TIME_AVG_)    = 'DIAG_TIME_AVG        '
    mon_str(TIME_AVG_RED_)     = 'TIME_AVG_RED         '
    mon_str(PWEIGHT_RED_)      = 'PWEIGHT_RED          '
    mon_str(DIAG2D_RED_)       = 'DIAG2D_RED           '
    mon_str(HEAT_RED_)         = 'HEAT_RED             '
    mon_str(GET_VOL_RED_)      = 'GET_VOL_RED          '

    mon_sync(:) = .false.
    mon_sync(INIT1_)           = timing_barrier
    mon_sync(INIT2_)           = timing_barrier
    mon_sync(CHARGEI_)         = timing_barrier
    mon_sync(CHARGEE_)         = timing_barrier
    mon_sync(POISSON_)         = timing_barrier
    mon_sync(DERIVSI_)         = timing_barrier
    mon_sync(DERIVSE_)         = timing_barrier
    mon_sync(DIAGNOSIS_)       = timing_barrier
    mon_sync(PUSHI_)           = timing_barrier
    mon_sync(PUSHE_)           = timing_barrier
    mon_sync(SHIFT_)           = timing_barrier
    mon_sync(COLL_SNAP_)       = timing_barrier
    mon_sync(COLLISION_)       = timing_barrier
    mon_sync(NEUTRAL_)         = timing_barrier
    mon_sync(DIFFUSION_)       = timing_barrier
    mon_sync(HEATING_)         = timing_barrier
    mon_sync(RESTART_WRITE_)   = timing_barrier
    mon_sync(FINALIZE_)        = timing_barrier
    mon_sync(SOLVER1_)         = timing_barrier
    mon_sync(SOLVER2_)         = timing_barrier
    mon_sync(SOLVER3_)         = timing_barrier
    mon_sync(PHASE4_)          = timing_barrier
    mon_sync(PHASE5_)          = timing_barrier
    mon_sync(CHARGEI_SR_)      = timing_barrier
    mon_sync(Z_Z_MODE_EXT_)    = timing_barrier
    mon_sync(CHARGEE_DENSITY_) = timing_barrier

    do i=1,mon_NX
       ctemp1 = mon_str(i)
       ctemp2(1: 5) = 'sync_'
       ctemp2(6:21) =  ctemp1(1:16)
       mon_str(i+mon_NX) = ctemp2
    enddo

#if !defined(NO_PETSC)
    do i=1, mon_N
       call PetscLogEventRegister( event(i), mon_str(i), 0, ierr )
       !call xgc_event_register( event(i) ,mon_str(i), ierr)
    enddo
#else
    mon_time(:)=0D0
    mon_sum(:)=0D0
#endif

  end subroutine init_perf_monitor

  subroutine finish_perf_monitor( ierr ) 
    use sml_module
#if defined(CAM_TIMERS)
    use perf_mod, only: t_finalizef, t_prf
#endif
    implicit none
    integer ierr,i
    include 'mpif.h'

#if defined(NO_PETSC)  
    ! deallocat psn%values ...
    if(sml_mype==0) then
       do i=1, mon_N
          write(*,1002) mon_str(i), mon_sum(i), mon_sum(i)/mon_sum(0)*100.
       enddo
    end if
1002 format(a10,e8.3,' ',f6.2,'%')
#else
    call petsc_end( ierr )
#endif

#if defined (CAM_TIMERS)
    call t_prf('timing_all', MPI_COMM_WORLD )
    call t_finalizef ()
#endif

  end subroutine finish_perf_monitor

end module perf_monitor

!---------------------- PETSC SOLVER ---------------------------------------------
#if !defined(NO_PETSC)
! 
!
!
module xgc_solver_module
#define PETSC_AVOID_DECLARATIONS
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h" 
#include "include/finclude/petscksp.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscpc.h"
#undef PETSC_AVOID_DECLARATIONS
  type xgc_solver
     Vec  xVec,bVec
     Mat  mat, rhs_mat, rhs2_mat
     KSP  ksp
     VecScatter scat
     IS petsc_xgc_is
     integer :: n_rhs_mat
     integer :: comm, mype, totalpe,prefix
  end type xgc_solver
  !
contains
  subroutine new_solver( this, nnodes, maxn1, maxn2, maxn3, ierr )
    !-----[--.----+----.----+----.-----------------------------------------]
    !     
    !     
    !     input:
    !     output:
    !     ierr: error code
    !     
    !-----[--.----+----.----+----.-----------------------------------------]
    use sml_module
    implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
    !#include "include/finclude/petscviewer.h"
    integer comm
    type(xgc_solver) :: this
    integer :: nnodes, maxn1,maxn2,maxn3, ierr
    real (kind=8) ::  tol
    Mat AA
    Vec x,y
    KSP ksp

    this%scat = 0 !PETSC_NULL_OBJECT
    this%petsc_xgc_is = 0 !PETSC_NULL_OBJEC

    this%mat=0
    this%rhs_mat=0
    this%rhs2_mat=0
    !
!    if( prefix .eq. 2 ) then
!       comm = sml_plane_comm
!    else
!       comm = petsc_comm_world
!    endif

    comm = this%comm


    ! LHS --------------
    call MatCreateMPIAIJ( comm, PETSC_DECIDE, PETSC_DECIDE, nnodes, nnodes, &
         maxn1, PETSC_NULL_INTEGER, maxn1/2, PETSC_NULL_INTEGER, &
         AA, ierr )
    CHKERRQ(ierr)  
    this%mat = AA
    call MatSetFromOptions(AA,ierr)
    CHKERRQ(ierr)
    call MatSetOption(AA,MAT_IGNORE_OFF_PROC_ENTRIES,ierr)
    CHKERRQ(ierr)
    call MatSetOption(AA,MAT_KEEP_ZEROED_ROWS,ierr)
    CHKERRQ(ierr)

    ! RHS -------------
    if(this%n_rhs_mat>=1) then
       !  call  MatCreate(comm,AA,ierr)
       call MatCreateMPIAIJ( comm, PETSC_DECIDE, PETSC_DECIDE, nnodes, nnodes, &
            maxn2, PETSC_NULL_INTEGER, maxn2/2, PETSC_NULL_INTEGER, &
            AA, ierr )
       CHKERRQ(ierr)
       this%rhs_mat = AA
       call MatSetFromOptions(AA,ierr)
       CHKERRQ(ierr)
       call MatSetOption(AA,MAT_IGNORE_OFF_PROC_ENTRIES,ierr)
       CHKERRQ(ierr)
    endif

    ! RHS2 -------------
    if(this%n_rhs_mat>=2) then
       !  call  MatCreate(comm,AA,ierr)
       call MatCreateMPIAIJ( comm, PETSC_DECIDE, PETSC_DECIDE, nnodes, nnodes, &
            maxn3, PETSC_NULL_INTEGER, maxn3/2, PETSC_NULL_INTEGER, &
            AA, ierr )
       CHKERRQ(ierr)
       this%rhs2_mat = AA
       call MatSetFromOptions(AA,ierr)
       CHKERRQ(ierr)
       call MatSetOption(AA,MAT_IGNORE_OFF_PROC_ENTRIES,ierr)
       CHKERRQ(ierr)
    endif

    ! KSP -------------
    call KSPCreate(comm, ksp, ierr )
    CHKERRQ(ierr)
    if( this%prefix .eq. 2 ) then
       call KSPSetOptionsPrefix( ksp, 's2_', ierr )
       CHKERRQ(ierr)
    end if
    this%xVec = 0
    this%bVec = 0

    this%ksp = ksp
    tol = 1.d-6
    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
    CHKERRQ(ierr)
    call KSPSetFromOptions( ksp, ierr )
    CHKERRQ(ierr)

    call MatZeroEntries(AA,ierr)
    CHKERRQ(ierr)

  end subroutine new_solver

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine delete_solver( this, ierr )
    implicit none
    type(xgc_solver) :: this
    integer ierr
    if( this%bVec .ne. 0 ) then
       call VecDestroy(this%xVec,ierr)
       ! if(ierr .ne. 0) call MPI_ABORT(PETSC_COMM_WORLD,1,ierr)
       !call VecDestroy(this%bVec,ierr)
       !if(ierr .ne. 0) call MPI_ABORT(PETSC_COMM_WORLD,1,ierr)

       !       call VecScatterDestroy( this%scat, ierr )
       ! if(ierr .ne. 0) call MPI_ABORT(PETSC_COMM_WORLD,1,ierr)

       ! call ISDestroy( this%petsc_xgc_is, ierr )
    endif

    if( this%mat .ne. 0 ) then
       call MatDestroy( this%mat,ierr)
       ! if(ierr .ne. 0) call MPI_ABORT(PETSC_COMM_WORLD,1,ierr)
    endif

    if( this%rhs_mat .ne. 0) then
       call MatDestroy( this%rhs_mat,ierr)
    endif

    if( this%rhs2_mat .ne. 0) then
       call MatDestroy( this%rhs2_mat,ierr)
    endif

    if( this%ksp .ne. 0 ) then
       call KSPDestroy(this%ksp,ierr)
       ! if(ierr .ne. 0) call MPI_ABORT(PETSC_COMM_WORLD,1,ierr)            
    endif
  end subroutine delete_solver
end module xgc_solver_module


#endif

module mat_class
  type mat_type
     integer :: N,width
     real (kind=8), pointer :: value(:,:)
     integer, pointer :: eindex(:,:),nelement(:)
  end type mat_type
contains
  ! create a empty matrix
  subroutine new_mat(mat,n,w)
    implicit none
    type(mat_type) :: mat
    integer,intent(in) :: n,w
    mat%n=n
    mat%width=w
    allocate(mat%value(w,n), mat%eindex(w,n), mat%nelement(n))
    mat%nelement=0
    mat%value=0D0
  end subroutine new_mat
  
  !delete a matrix
  subroutine del_mat(mat)
    implicit none
    type(mat_type) :: mat

    if(associated(mat%value)) then
       deallocate(mat%value,mat%eindex,mat%nelement)
       nullify(mat%value)
       nullify(mat%eindex)
       nullify(mat%nelement)
    endif
    mat%n=-1
    mat%width=-1
  end subroutine del_mat

    ! set value
  subroutine set_value(mat,i,j,value,flag)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: i,j,flag
    real (kind=8) :: value
    integer :: l

    !check same node
    do l=1,mat%nelement(i)
       ! redundant point
       if(j==mat%eindex(l,i))then
          if(flag==0) then
             mat%value(l,i)=value
          else
             mat%value(l,i)=mat%value(l,i)+value
          endif
          return
       endif
    enddo

    ! new point
    if(mat%width <= mat%nelement(i) ) then
       print *, 'Error : not enough memory space for matrix - too many off diagnonal comp.'
       print *, '        Increase width ', mat%width, mat%nelement(i),i
       return
    endif
    mat%nelement(i)=mat%nelement(i)+1
    l=mat%nelement(i)
    mat%eindex(l,i)=j
    mat%value(l,i)=value
  end subroutine set_value
  subroutine get_max_width(mat,mwidth)
    implicit none
    type(mat_type) :: mat
    integer, intent(out) :: mwidth
    integer :: i

    mwidth=-1
    do i=1, mat%n
       if(mat%nelement(i)> mwidth) then
          mwidth=mat%nelement(i)
       endif
    enddo
  end subroutine get_max_width
  subroutine output_matrix(mat,filename)
    implicit none
    integer :: i,j
    type(mat_type) :: mat
    character (len=20) :: filename 
    open(unit=201,file=filename,status='replace')
    do i=1, mat%n
     do j=1, mat%nelement(i)
        write(201,*) i,mat%eindex(j,i),mat%value(j,i)
     enddo
    enddo
    close(201)

  end subroutine output_matrix

end module mat_class

module boundary_class
  type range_type
     integer :: start,end     
  end type range_type
  type boundary_type
     type(range_type) :: in,out1,out2
  end type boundary_type
contains
  logical function is_inside(i,bd) 
    ! check if inside of boundary
    ! return .true. if ith node point is inside of simulation boundary (extended)  
    ! type == 1 check only outer boundary
    ! type == 2 check both inner and outer boundary 
    implicit none
    type(boundary_type) :: bd
    integer , intent(in) ::i
    
    if( (bd%in%start<= i .and. i<= bd%in%end) .or. & 
          (bd%out1%start <= i .and. i<= bd%out1%end) .or. &
          (bd%out2%start <= i .and. i<= bd%out2%end) ) then 
       is_inside=.false.
    else
       is_inside=.true.
    endif
    
  end function is_inside

  subroutine set_boundary_values(arr,value,bd)
    implicit none
    type(boundary_type) :: bd
    real (kind=8) :: arr(bd%out2%end), value


    if(bd%in%end>0) then
       arr(bd%in%start:bd%in%end)=value
    endif
    arr(bd%out1%start:bd%out1%end)=value
    arr(bd%out2%start:bd%out2%end)=value
    
  end subroutine set_boundary_values

end module boundary_class
module psn_class
  use mat_class
#if !defined(NO_PETSC)
  use xgc_solver_module
  use boundary_class
#endif
  type psn_type
#if !defined(NO_PETSC)
     type(xgc_solver) :: solver00,solverH
#endif
     !field variables
     ! move to proper class later
     type(mat_type) :: mat
     real (kind=8), pointer :: idensity(:,:),idensity0(:),tempi_ev(:),factor1(:)
     real (kind=8), pointer :: tite(:),dpot(:,:),rhoi2(:),tempe(:)  !! ion density
     real (kind=8), pointer :: edensity0(:),ddensity0(:),density0_full(:),f0_density0(:)
     real (kind=8), pointer :: pot0(:), add_pot0(:)     
     
     real (kind=8), pointer :: E_para(:,:)
     real (kind=8), pointer :: E_perp_node(:,:,:), E_perp0_node(:,:)
     real (kind=8), pointer :: E_perp_tr(:,:,:), E_perp0_tr(:,:)
     
     real (kind=8), pointer :: bfollow_p(:,:,:)
     integer, pointer :: bfollow_tr(:,:)
     real (kind=8), pointer :: bfollow_1_dx(:)
#ifdef CIRCULAR_SPECIAL2
     real (kind=8), pointer :: pot_cc2(:,:,:), epara_cc2(:,:,:)
     real (kind=8) :: dt_cc2
     integer :: nt_cc2
#endif
     integer :: nwall, wall_start
     real (kind=8), pointer :: sheath_pot(:)

     !boundary
     type(boundary_type) :: cbd0, pbd0, cbdh, pbdh ! charge / potential boundary of High-k solver / 0 solver
  end type psn_type
  
contains

  subroutine psn_mem_alloc(psn,n,m,ntr)
!    use sml_module
    implicit none
    type(psn_type) :: psn
    integer :: n,m, ntr !node, nphi, ntriangle
    allocate( psn%idensity(n,0:m), psn%idensity0(n), psn%tempi_ev(n), psn%factor1(n), &
         psn%tite(n),psn%rhoi2(n),psn%tempe(n), &
         psn%dpot(n,-1:m+1), psn%edensity0(n), psn%ddensity0(n),psn%density0_full(n),psn%pot0(n), &
         psn%E_para(n,0:m),psn%E_perp_tr(2,ntr,0:m),psn%E_perp0_tr(2,ntr) , &
         psn%E_perp_node(2,n,0:m),psn%E_perp0_node(2,n) )
    allocate( psn%bfollow_p(3,2,n), psn%bfollow_tr(2,n), psn%bfollow_1_dx(n) )

    
    psn%idensity=0D0
    psn%idensity0=0D0
    psn%tempi_ev=0D0
    psn%factor1=0D0
    psn%tite=0D0
    psn%tempe=0D0
    psn%dpot=0D0
    psn%edensity0=0D0
    psn%ddensity0=0D0
    psn%density0_full=0D0
    psn%pot0=0D0
    psn%E_para=0D0
    psn%E_perp_tr=0D0
    psn%E_perp0_tr=0D0
    psn%E_perp_node=0D0
    psn%E_perp0_node=0D0
    psn%bfollow_p=0D0
    psn%bfollow_tr=0
    psn%bfollow_1_dx=0D0
    psn%rhoi2=0D0
    
  end subroutine psn_mem_alloc

  subroutine psn_sheath_alloc(psn)
    implicit none
    type(psn_type) :: psn
    allocate(psn%sheath_pot(psn%nwall))
  end subroutine psn_sheath_alloc
end module psn_class

module smooth_module
  use mat_class
  type smooth_type     
     real (kind=8), pointer :: weight(:)
     integer :: n, mode, type
  end type smooth_type
  type smooth_r_type
     type(mat_type) :: mat
     real(kind=8) :: d0 ! basic distance
     integer :: n ! number of point
     integer :: type 
  end type smooth_r_type
  type(smooth_type) :: smooth00,smooth0L,smoothH,smoothdiag
  type(smooth_r_type) :: smooth_r1
  integer :: smooth_n_in, smooth_mode_in,smooth_type_in
  integer :: smooth_H_n_in, smooth_H_mode_in, smooth_H_type_in
  integer :: smooth_diag_n_in, smooth_diag_mode_in, smooth_diag_type_in  
  ! input parameters for smooth_r
  real (kind=8) :: smooth_r1_d0_in
  integer :: smooth_r1_n_in, smooth_r1_type_in
  contains 
    ! initialize smooth_init    
    subroutine smooth_init(smooth)
      implicit none
      type(smooth_type) :: smooth
      integer :: i
      
      if(smooth%mode==-1) return  ! No smooth
      if(smooth%mode== 0) return  ! flux average
      
      allocate(smooth%weight(smooth%n))
      do i=1, smooth%n
         if(smooth%type==1) then
            smooth%weight(i)=exp(- real(i-1)**2/real(smooth%n)**2*4D0 )
         else
            smooth%weight(i)=1D0
         endif
      enddo
      smooth%weight(1)=0.5D0*smooth%weight(1)
    end subroutine smooth_init
    !
    ! delete smooth_pol object
    subroutine smooth_pol_delete(smooth)
      implicit none
      type(smooth_type) :: smooth
      if(associated(smooth%weight)) then
         deallocate(smooth%weight)
         nullify(smooth%weight) 
      endif
    end subroutine smooth_pol_delete
    !
    ! initialize smooth_r
    subroutine smooth_r_init1(smooth_r,d0,n,type,nnode)
      implicit none
      type(smooth_r_type) :: smooth_r
      real (kind=8) :: d0
      integer :: n, type,nnode
      
      smooth_r%n=n
      smooth_r%d0=d0
      smooth_r%type=type

      call new_mat(smooth_r%mat,nnode,(2*n+1)*3)     

      ! call smooth_r_init2 with grid information

    end subroutine smooth_r_init1
    !
    ! delete smooth_r object
    subroutine smooth_r_delete(smooth_r)
      implicit none
      type(smooth_r_type) :: smooth_r
!      if(allocated(smooth_r%mat)) then
         call del_mat(smooth_r%mat)
!      endif
    end subroutine smooth_r_delete
end module


module rtemp_module
  type rtemp_type
     integer :: np
     real (kind=8) :: pin, pout, dp, inv_dp
     real (kind=8),pointer ::dtemp(:)
  end type rtemp_type
  type(rtemp_type) :: rtempi
  
  contains
    subroutine retore_temp_setup(rtemp,n,pin,pout)
      use sml_module
      implicit none
      type(rtemp_type) :: rtemp
      integer :: n
      real (kind=8) :: pin, pout

      !check validity of simulation parameters
      if(sml_deltaf/=1 .and. sml_mype==0) then
         print *, '***************** *     Warning       *****************************'
         print *, 'restore temperature is designed only for delta-f simulation'
         print *, '*******************************************************************'
      endif
      if(sml_deltaf_f0_mode/=-1 .and. sml_mype==0) then
         print *, '***************** *     Warning       *****************************'
         print *, 'restore temperature is designed only for deltaf_f0_mode==-1 simulation'
         print *, '*******************************************************************'
      endif
      if(sml_push_mode/=1 .and. sml_mype==0) then
         print *, '***************** *     Warning       *****************************'
         print *, 'restore temperature is designed only for RK2 method'
         print *, 'Error of weight evolution maybe happens with P-C particle pushing'
         print *, '*******************************************************************'
      endif

     ! set variables from input, allocate memory, initialize and set dpsi variable 
      rtemp%np=n
      rtemp%pin=pin
      rtemp%pout=pout
      rtemp%dp= (pout-pin)/real(n)
      rtemp%inv_dp=1D0/rtemp%dp
      allocate(rtemp%dtemp(n))
      rtemp%dtemp=0D0 ! delta temp = 0.

    end subroutine retore_temp_setup

    subroutine restore_temp(sp, rtemp)
      use sml_module
      use ptl_module
      use eq_module
      implicit none
      type(species_type) :: sp
      type(rtemp_type) :: rtemp
      integer :: i,j
      real (kind=8) :: dtem(rtemp%np), wsum(rtemp%np)
      real (kind=8), allocatable :: e(:)
      real (kind=8) :: r,z,psi,phi, b
      real (kind=8) :: c2_2m, inv_temp, alpha
      real (kind=8), external :: b_interpol
      ! allocate temp array (2,ptl num)
      allocate(e(sp%num))
      
      ! initialization
      dtem=0D0
      wsum=0D0
      c2_2m=ptl_c2_2m(sp%type)
      inv_temp=1D0/(eq_tempi_ev_out*sml_ev2j) !unit: 1 over joule

      ! get energy of each particle and flux sum (delta-f energy and full-f weight)
      do i=1,sp%num
         ! find index
         psi=sp%phase(9,i)
         j=(psi-rtemp%pin)*rtemp%inv_dp + 1
         if(1<=j .and. j<=rtemp%np) then
            ! get B-field        
            r=sp%phase(1,i)
            z=sp%phase(2,i)
            phi=sp%phase(3,i)
            b=b_interpol(r,z,phi)
            ! get E
            e(i)=c2_2m*(sp%phase(4,i)*b)**2 + sp%phase(5,i)*b   
            ! store it
            dtem(j)=dtem(j) + e(i)* sp%phase(6,i) * sp%phase(8,i)
            wsum(j)=wsum(j) + sp%phase(8,i)
         endif
      enddo
      
      ! allreduce 
      call my_mpi_allreduce(dtem,dtem,rtemp%np)
      call my_mpi_allreduce(wsum,wsum,rtemp%np)
      ! get delta temperature, decaying rate
      wsum=max(wsum,1D0) ! to prevent NaN.
      dtem=dtem/wsum*inv_temp  ! get normalized delta T_i relative to T0
      alpha=sml_restore_temp_avg_scale*real(sml_restore_temp_period)
      rtemp%dtemp = (1D0-alpha)*rtemp%dtemp + alpha*dtem
      
      ! Change weight of each particle
      do i=1, sp%num
         psi=sp%phase(9,i)
         j=(psi-rtemp%pin)*rtemp%inv_dp + 1
         if(1<=j .and. j<=rtemp%np) then
            sp%phase(6,i)=sp%phase(6,i) - (e(i)*inv_temp - 1.5d0)*rtemp%dtemp(j)
         endif
      enddo
      ! deallocate
      deallocate(e)
#ifdef XGC_DEBUG_RESTORE_TEMP
      if(mod(sml_istep,sml_restore_temp_period*10) .and. sml_mype==0) then
         open(unit=2345,file='fort.rtemp',position='append')
         do j=1, rtemp%np
            write(2345,*) rtemp%pin+rtemp%dp*(real(j)-0.5), rtemp%dtemp(j),dtem(j)
         enddo
         write(2345,*) ' '
         close(2345)
      endif
#endif
    end subroutine restore_temp

  end module rtemp_module
  
