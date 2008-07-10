subroutine setup(ptl)
  use sml_module
  use ptl_module
  use itp_module
  use eq_module
  use heat_module
  use tbl_module
  use neu_module
  use efld_module
  use lim_module
  use diag_module
  use col_module
  use bnc_module
  use rf_module
  use rpl_module
  use grid_class
  use smooth_module
  use rtemp_module
  implicit none
  type(ptl_type) :: ptl
  integer :: i,j,ij,ierror,sp
  integer, parameter :: nspace=2
  real (kind=8) :: vol, meanr, r1,r2,z1,z2
  real (kind=8) :: v_the,psi
  real (kind=8), external :: tempe_ev,init_tempi_ev,init_den
  character (len=20) :: ptlfile, varname
  
  namelist /sml_param/ sml_restart, sml_bp_mult, sml_bt_mult,&
       sml_special, sml_machine,sml_mstep, sml_en_order_kev,sml_dt,sml_deltaf,&
       sml_concentric,&
       sml_canonical_maxwell,sml_deltaf_f0_mode,sml_marker_temp_factor,&
       sml_restart_write_period,sml_monte_num,sml_inpsi,sml_outpsi,&
       sml_bounce,sml_bounce_zero_weight,sml_bt_sign, sml_minusB,&
       sml_read,sml_limiter,&
       sml_bd_min_r,sml_bd_max_r,sml_bd_min_z,sml_bd_max_z,&
       sml_relax_dist, sml_relax_dist_num, &
       sml_f0den_edge, sml_f0den_out, sml_f0den_ped_c, sml_f0den_ped_width,&
       sml_f0_1_Ln, sml_f0_1_Lt, &
       sml_node_file, sml_ele_file, sml_bfollow_file,&       
       sml_turb_efield, sml_no_00_efield, & 
       sml_bfollow, sml_bfollow_read,&
       sml_tavg_factor, sml_tavg_on,&
       sml_electron_on, sml_nphi_total, sml_fem_matrix, sml_smtrc_ad_elec, &
       sml_push_mode, sml_pc_order, &
       sml_dwdt_exb_only,&
       sml_guess_table_size,&
       sml_initial_deltaf_noise,&
       sml_zero_inner_bd,&
       sml_mode_select_on, sml_mode_select_n, sml_mode_initial_n,&
       sml_bd_ext_delta1, sml_bd_ext_delta2, sml_bd_ext_delta3, sml_bd_ext_delta4, &
       sml_input_file_dir,&
       sml_dwdt_fix_bg, sml_no_para_nonlinear, sml_add_pot0, sml_add_pot0_file, &
       sml_flat_electron_density,&
       sml_supress_weight_growth, sml_weight_max,&
       sml_restore_temp, sml_restore_temp_period,sml_restore_temp_avg_scale, sml_use_simple00,&
       sml_simple00_nsmth, &
       sml_max_mat_width, &
       sml_bd_Te_mode, sml_bd_Te_width, sml_zero_out_total_charge, &
       sml_sheath_mode, sml_sheath_init_pot_factor,sml_exclude_private,sml_rgn1_pot0_only, &
       sml_use_pade, sml_heat_on

  namelist /ptl_param/  ptl_mass_au,ptl_charge_eu,ptl_num,& 
       ptl_e_mass_au, ptl_e_charge_eu, ptl_e_num, ptl_lost_nummax, ptl_maxnum,ptl_e_maxnum, &
       ptl_special_r,ptl_special_z,ptl_special_phi,&
       ptl_special_en_ev,ptl_special_pitch


  namelist /eq_param/ eq_filename,eq_den_edge, eq_den_out, eq_den_ped_c, &
       eq_den_ped_width, &
       eq_tempi_ev_edge,eq_tempi_ev_out,eq_tempi_ped_c,eq_tempi_ped_width,&
       eq_tempe_ev_edge,eq_tempe_ev_out,eq_tempe_ped_c,eq_tempe_ped_width
       
  namelist /col_param/ col_mode, col_period, col_en_col_on, col_accel,&
       col_2_pin, col_2_pout, col_2_dtheta,&
       col_varying_bg, col_vb_m, col_vb_period, col_vb_pin, col_vb_pout,&
       col_imp_on, col_imp_charge,col_imp_mass,& 
       col_den_imp_edge,col_den_imp_out,col_den_imp_ped_c,col_den_imp_ped_width
  
  namelist /diag_param/ diag_tracer_period, diag_tracer_n, diag_tracer_sp, &
       diag_efld_period, diag_pw_on, diag_pw_period, &
       diag_flow_period, diag_flow_pin,diag_flow_pout, &
       diag_pot_period, diag_binout_period, &
       diag_avg_on,diag_avg_outperiod,diag_avg_pin,diag_avg_pout,&
       diag_f_on,diag_f_skip,diag_f_mod,&
       diag_f_vmax,diag_f_pin,diag_f_pout,diag_f_angle1,diag_f_angle2,&
       diag_gam_on, diag_gam_node_begin, diag_gam_node_num, &
       diag_ptl_on, diag_ptl_begin, diag_ptl_end, diag_ptl_num
  
  namelist /efld_param/ efld_mode, efld_1_psie,efld_1_dpsi,efld_1_psim,&
       efld_1_pot0,efld_1_potm,&
       efld_pin,efld_pout,efld_dpdp0_outside,&
       efld_cutoff,efld_cutoff_psi,efld_max_dpdp,&
       efld_reset,efld_reset_step,efld_start,efld_read_dpdp

!  namelist /rf_param/ rf_on, rf_r_res, rf_za, rf_lap, rf_e0

  namelist /neu_param/ neu_sepfile,neu_col_mode,neu_adjust_n0,&
       neu_ionize_mode,neu_ionize2_psi,&
       neu_cx_period,neu_ion_period,neu_col_period,&
       neu_start_time, neu_recycle_rate,&
       neu_varying_mfp,&
       neu_t0,neu_n0,neu_delta_n,neu_delta_theta,neu_mfp0,&
       neu_theta_x, neu_temp_factor,&
       neu_monte_num,neu_mode2_period,neu_mpol,&
       neu_temp0,neu_dt,neu_istep_max,neu_psi_edge,neu_grid_max_psi,&
       neu_grid_min_psi,&
       neu_sepfile
       
  namelist /lim_param/ lim_store_mz,lim_filename
  namelist /smooth_param/ smooth_mode_in, smooth_n_in, smooth_type_in, &
       smooth_H_mode_in, smooth_H_n_in, smooth_H_type_in, &
       smooth_diag_mode_in, smooth_diag_n_in, smooth_diag_type_in,&
       smooth_r1_d0_in, smooth_r1_n_in, smooth_r1_type_in
  namelist /tbl_param/ tbl_diffusion_on, tbl_stop_time, tbl_period, tbl_d_coeff
  namelist /heat_param/ heat_period, heat_inpsi, heat_outpsi, heat_power

  ! set processor ID 
  call check_point2('sml_param')
  sml_nphi_total=16 ! default 16

  ! read sml_parameters
  open(unit=14,file='input',action='read')
  read(14,nml=sml_param)
  close(14)
  call reset_nphi_total

  ! set domain/processor relations
  if( sml_nphi_total >= sml_totalpe ) then
       sml_pe_per_plane=1
       sml_plane_per_pe=sml_nphi_total/sml_totalpe
  else
       sml_plane_per_pe=1 
       sml_pe_per_plane=sml_totalpe/sml_nphi_total
  endif

  if(sml_totalpe*sml_plane_per_pe/=sml_nphi_total*sml_pe_per_plane) then
     print *, 'Invalid PE number and sml_nphi_total'
     stop
  endif

  !for now
!  if(sml_nphi_total < sml_totalpe) then
!     print *, 'sml_nphi_total cannot be smaller than sml_totalpe'
!     print *, 'This feature is under developing now'
!     print *, 'Set sml_nphi_total = number of processor in input'
!     stop
!  endif

  ! Make a list of te processes in the new communicator 
  call new_communicator
  
  ! 0. constant - change to parameter later
  sml_pi=4D0*datan(1D0)
  sml_2pi=8D0*datan(1D0)
  sml_sqrtpi=sqrt(sml_pi)
  sml_sqrt2pi=sqrt(sml_2pi)

  ! 1.  set deafault value to phase1 variables
  sml_input_file_dir='./'
  sml_restart=0
  sml_bp_mult=1D0
  sml_bt_mult=1D0
  sml_special=0  ! 0 : nomral simulation
                 ! 1 : single particel simulation
                 ! 2 : undefined
                 ! 3 : velocity hole calculation
  sml_machine=0  !                0 : circular 
                 ! macine number  1: D3D
                 !                2: NSTX
                 !                3: alcator-CMOD
                 !                4: ASDEX-U    
  sml_mstep=3000          ! total time step. (simulation time) = dt*mstep
  sml_en_order_kev = 0.2D0  ! energy order in KeV

  sml_dt=0.001D0           !unit of tran --> code unit
  sml_deltaf=0             !0: full f, 1: delta f - delta f simulation is not verified yet
  sml_canonical_maxwell=0   !  cannonical distribution - not implimented yet
  sml_deltaf_f0_mode=0
  sml_marker_temp_factor=1D0   ! Marker temperature ratio -- (Marker when weight =1 )
  sml_restart_write_period=10000000
  sml_inpsi=0.9             ! eq_x_psi unit --> code unit
  sml_outpsi=1.05           ! eq_x_psi unit --> code unit
  ! if current I_T is in -phi direction, sml_minusB=0
  ! if current I_T is in the phi direction, sml_minusB=1
  ! phi is the cylindrical angle
  sml_bt_sign=-1            !sml_minusB=1 , vec B -> - vec B
  sml_minusB=0              ! 0 if cocurrent, 1 if counter current

  sml_electron_on=0

  sml_relax_dist=0
  sml_relax_dist_num=1000
  sml_turb_efield=0


  sml_tavg_on=0
  sml_tavg_factor= 1D-3

  sml_fem_matrix=1
  sml_smtrc_ad_elec=0
  sml_push_mode=1
  sml_pc_order=3
  sml_dwdt_exb_only=0
  sml_dwdt_fix_bg=.false.
  sml_no_para_nonlinear=0
  sml_guess_table_size=1000
  sml_initial_deltaf_noise=1D-3
  sml_zero_inner_bd=0
  sml_mode_select_on=0
  sml_mode_select_n=1
  sml_mode_initial_n=0
  sml_f0_1_Lt=6.92D0
  sml_f0_1_Ln=sml_f0_1_Lt/3.114D0
  sml_bounce_zero_weight=0

  sml_bd_ext_delta1=0.01 ! inner boundary extension for poisson eq
  sml_bd_ext_delta2=0.02 ! inner boundary extension for charge zeroing out
  sml_add_pot0=0
  sml_add_pot0_file='pot0.dat'
  sml_flat_electron_density=0
  sml_heat_on=.false.

  sml_supress_weight_growth=.false.  ! if this value is true then weight of deltaf will be between -sml_weight_max and sml_weight_max. default value is .false.
  sml_weight_max=10.   ! maximum weight(= phase(6)) for delta-f simulation. default value is 10.

  ! restore temperature gradient
  sml_restore_temp = .false.  ! restore temperature gradient. Work for delta-f=1, sml_deltaf_f0_mode=-1, and rk2 method only
  sml_restore_temp_period=5 ! how frequently restore temperature
  sml_restore_temp_avg_scale=0.01 ! time average scale - dt base, default=0.01 (means 100 dt scale)

  sml_bd_Te_mode=0
  sml_bd_Te_width=0.01D0

  sml_use_pade=.true.

  ! for concentric circular, 00-mode poission equation is solved analytically when we use poisson_two_solvers
  ! it has smoothing routine in it and the number of smoothing operation of charge density is give by sml_simple00_nsmth
  sml_use_simple00=.false.
  sml_simple00_nsmth=0

  sml_zero_out_total_charge=.true.

  ! read sml_parameters
  open(unit=14,file='input',action='read')
  read(14,nml=sml_param)
  close(14)
  call reset_nphi_total
  
  sml_bd_ext_delta3=sml_bd_ext_delta1 ! outer boundary extension for poisson eq
  sml_bd_ext_delta4=sml_bd_ext_delta2 ! outer boundary extension for charge zeroing out

  
  ! 2 for core simulation, 1 for edge simulation 
  if(sml_outpsi>1D0) then
     sml_bounce=1
  else
     sml_bounce=2
  endif

  if (sml_machine==0) then
!     sml_concentric=1
     sml_concentric=0
     sml_read=1
     sml_limiter=0
     eq_filename='circular.eqd'
  else
     sml_concentric=0
     sml_read=1
     sml_limiter=0
     if(sml_machine==1) then
        eq_filename="g096333.eqd"
     elseif(sml_machine==2) then
        eq_filename="g104316_nstx.eqd"
     elseif(sml_machine==3) then
        print *, 'not prepared yet'
        stop
     elseif(sml_machine==4) then
        eq_filename="g17151.ASDEXU.eqd"
     endif
  endif
 
  if(sml_machine==1) then
     ! simulation boundary for D3D
     sml_bd_min_r=0.7  ! inner boundary
     sml_bd_max_r=2.5  ! outer boundary
     sml_bd_min_z=-1.8 ! lower boundary
     sml_bd_max_z=1.8  ! upper boundary
  elseif(sml_machine==2) then
                                                                                
     ! simulation boundary for NSTX
     sml_bd_min_r=0.18 ! inner boundary
     sml_bd_max_r=1.6  ! outer boundary
     sml_bd_min_z=-1.8 ! lower boundary
     sml_bd_max_z=1.5  ! upper boundary
  elseif(sml_machine==3) then
                                                                                
     ! simulation boundary for alcator-CMOD
     sml_bd_min_r=0.4  ! inner boundary
     sml_bd_max_r=1.   ! outer boundary
     sml_bd_min_z=-.7  ! lower boundary
     sml_bd_max_z=0.65 ! upper boundary
  elseif(sml_machine==4) then
     ! simulation boundary for ASDEX-U
     sml_bd_min_r=1D0
     sml_bd_max_r=2.25D0 
     sml_bd_min_z= -1.5D0 
     sml_bd_max_z=1.3 
  else
     sml_bd_min_r=0.0001
     sml_bd_max_r=1000.
     sml_bd_min_z=-1000.
     sml_bd_max_z=1000.
  endif


!  sml_node_file="coarsed3d.1.node"
!  sml_ele_file="coarsed3d.1.ele"
!  sml_bfollow_file="coarsed3d.bf.dat"
  sml_node_file="neo10.1.node"
  sml_ele_file="neo10.1.ele"
  sml_bfollow_file="neo10.bf.dat"
  sml_bfollow=0
  sml_bfollow_read=0
  sml_max_mat_width=100
  
  sml_exclude_private=.true.
  sml_rgn1_pot0_only= .false.

  sml_sheath_mode=0
  sml_sheath_init_pot_factor=2D0
  ! read sml_parameters
  open(unit=14,file='input',action='read')
  read(14,nml=sml_param)
  close(14)
  call reset_nphi_total
  


  !---------------------------------------------------------------
  call check_point2('ptl_param')
  ptl_mass_au=2D0  !mass ratio to proton
  ptl_e_mass_au=2D-2
  ptl_charge_eu=1D0  ! charge number
  ptl_e_charge_eu=-1D0

  ! special simulation set definition
  if(sml_special==1) then !single particle simulation
     ptl_maxnum=10
     ptl_num=1
     ptl_special_r=1.95D0
     ptl_special_z=0.
     print *, ptl_special_r, ptl_special_z
     ptl_special_phi=0D0
     ptl_special_en_ev=7205
     ptl_special_pitch=0.276D0
  else if(sml_special==3) then
     ptl_num=5000
     ptl_special_n=100
     ptl_special_r=1.95
     ptl_special_z=0.
     ptl_special_phi=0D0
     ptl_special_en_ev=1000D0
  else
     ptl_maxnum=4000
     ptl_num=3000   !Particle number of each PE
  endif

  open(unit=14,file='input',action='read')
  read(14,nml=ptl_param)
  close(14)
  ptl_e_maxnum=ptl_maxnum
  ptl_e_num=ptl_num


  ptl_lost_nummax=max(ptl_num/10,100)

  open(unit=14,file='input',action='read')
  read(14,nml=ptl_param)
  close(14)

  if(sml_special==1) then
     if(ptl_num/=1) print *, 'Warning. ptl_num is not 1 for single ptl simulation'
  endif

  sml_monte_num=max(5*ptl_num,10000)
  if(sml_turb_efield==0) then
     ptl_maxnum=max(ptl_num,ptl_maxnum)
     ptl_e_maxnum=max(ptl_e_num,ptl_e_maxnum)
  else
     ptl_maxnum=max(ptl_num+ptl_num/5,ptl_maxnum)
     ptl_e_maxnum=max(ptl_e_num+ptl_e_num/5,ptl_e_maxnum)
  endif
  ptl%ion%num=ptl_num
  ptl%elec%num=ptl_e_num

  ptl_mass(1)=ptl_mass_au*sml_prot_mass
  ptl_charge(1)=ptl_charge_eu*sml_e_charge
  ptl_mass(2)=ptl_e_mass_au*sml_prot_mass
  ptl_charge(2)=ptl_e_charge_eu*sml_e_charge

  ptl_c_m(:)=ptl_charge(:)/ptl_mass(:)
  ptl_c2_2m(:)=0.5D0*ptl_charge(:)**2/ptl_mass(:)


  !-------------READ Equilbirum file------------------------------
  call check_point2('eq_param')
  open(unit=14,file='input',action='read') ! for eq_filename
  read(14,nml=eq_param)
  close(14)

  if(sml_mype==0) print *, 'reading......'
  call add_dir_name1
  call read
  
  !---------------------------------------------------------------
  eq_den_edge=  5D19   ! m^-3 
  eq_den_out =0.51D19  ! m^-3
  eq_den_ped_c= 0.5*(sml_inpsi+sml_outpsi) !x_psi unit -> code unit   
  eq_den_ped_width=0.001  ! x_psi unit -> code unit
  
  !read eq param 1st
  open(unit=14,file='input',action='read')
  read(14,nml=eq_param)
  close(14)

  
  ! tempi - ion temperature
  eq_tempi_ev_edge=1D3   !ev
  eq_tempi_ev_out=0.1D3  !ev
  eq_tempi_ped_c=eq_den_ped_c
  eq_tempi_ped_width=eq_den_ped_width

  !read eq param 2nd
  open(unit=14,file='input',action='read')
  read(14,nml=eq_param)
  close(14)
  
  !electron temperature
  eq_tempe_ev_edge=eq_tempi_ev_edge
  eq_tempe_ev_out=eq_tempi_ev_out
  eq_tempe_ped_c=eq_tempi_ped_c
  eq_tempe_ped_width=eq_tempi_ped_width

  ! read eq_parameters
  open(unit=14,file='input',action='read')
  read(14,nml=eq_param)
  close(14)

  ! set dependant sml_parameters
  sml_f0den_edge=eq_den_edge
  sml_f0den_out=eq_den_out
  sml_f0den_ped_c=eq_den_ped_c
  sml_f0den_ped_width=eq_den_ped_width

  ! re-read sml_parameters
  open(unit=14,file='input',action='read')
  read(14,nml=sml_param)
  close(14)
  call reset_nphi_total



  !-----------------------------------------------------
  ! Collision parameters
  call check_point2('col_param')

  col_mode=1  ! 0: off 1 : monte-carlo (non-conserving)
              ! 2: conserving collision (monte-carlo) 
  col_period=3  ! collsion each 3 time step
  col_en_col_on=0  ! energy collision on / off
  col_accel=1D0    ! artifitial increase of collision

  col_2_pin=sml_inpsi  ! x_psi unit -> code unit
  col_2_pout=min(sml_outpsi,0.999D0)  ! x_psi unit -> code unit
  
  col_varying_bg=0  ! change background profile for collision as time goes
  col_vb_m=30  
  col_vb_period=20000  ! period of change
  col_vb_pin=sml_inpsi        ! x_psi unit -> code unit
  col_vb_pout=sml_outpsi      ! x_psi unit -> code unit
  
  col_imp_on=0        ! impurity collision on/off 
  col_imp_charge=4.5  ! impurity charge number
  col_imp_mass=12.    ! impurity mass (AU)
  col_den_imp_edge=eq_den_edge*0.027D0  ! m-3 -> code unit
  col_den_imp_out =eq_den_out *0.027D0  ! m^-3 -> code unit
  col_den_imp_ped_c =  eq_den_ped_c     ! x_psi unit -> code unit
  col_den_imp_ped_width=eq_den_ped_width ! x_psi unit -> code unit

  ! read col_parameters
  open(unit=14,file='input',action='read')
  read(14,nml=col_param)
  close(14)

  col_vb_period=max(col_vb_period,1)

  !----------------------------------------------------
  ! diagnosis
  call check_point2('diag_param')
  ! tracer : record time history of a given particle
  diag_tracer_period=100 
  diag_tracer_n=1  ! index number of traced particle
  diag_tracer_sp=1 ! species index number of traced particle

  ! efield change output
  diag_efld_period=sml_mstep/50 

  ! particle weight ("noise") diagnostic for delta-f simulation 
  diag_pw_on=0
  diag_pw_period=max(1,sml_mstep/50000)


  ! density and flow 
  diag_flow_period=sml_mstep/50
  diag_flow_pin=sml_inpsi         !x_psi unit -> code unit
  diag_flow_pout=sml_outpsi       !x_psi unit -> code unit

  ! 2d potential and binpack data dumps
  diag_pot_period=sml_mstep/5
  diag_binout_period=sml_mstep/5

  ! time averaged calculation - transport coefficient
  diag_avg_on=0
  diag_avg_outperiod=sml_mstep/50
  diag_avg_pin=sml_inpsi          !x_psi unit -> code unit
  diag_avg_pout=sml_outpsi        !x_psi unit -> code unit

  ! Distribution function on midplane (between small angle gap)
!  diag_f_mod=100       ! # of step averaged
!  diag_f_skip=sml_mstep/5/diag_f_mod     ! skip*mod is period.
  diag_f_on=0
  diag_f_skip=5
  diag_f_mod=sml_mstep/5

  diag_f_vmax = 3 ! realtive to thermal velocity of ion
  diag_f_pin=sml_inpsi             ! x_psi unit -> code unit
  diag_f_pout=sml_outpsi           ! x_psi unit -> code unit
  diag_f_angle1=10                 ! degree -> radian
  diag_f_angle2=-10                 ! degree -> radian

  diag_gam_on=0
  diag_gam_node_begin=0
  diag_gam_node_num=15

  diag_ptl_on=0
  diag_ptl_begin=1
  diag_ptl_end=sml_mstep
  diag_ptl_num=ptl_num
  
  ! read diag_parameters
  open(unit=14,file='input',action='read')
  read(14,nml=diag_param)
  close(14)

  diag_avg_outperiod=max(diag_avg_outperiod,1)
  diag_flow_period=max(diag_flow_period,1)
  diag_efld_period=max(diag_efld_period,1)
  diag_pot_period=max(diag_pot_period,1)
  diag_binout_period=max(diag_binout_period,1)

  if(sml_mype/=0) then
     diag_tracer_n=0
  endif

  if (sml_mype==0 .and. diag_ptl_on==1) then
     allocate(diag_ptl_times((diag_ptl_end-diag_ptl_begin)+1),&
              diag_ptl_data1(diag_ptl_num),diag_ptl_data2(diag_ptl_num))
  endif
#ifdef BINPACK
     ptlfile="xgc.ion.bp"
     call bpopenfile(ptlfile,diag_ptl_ifile)
! dummy mesh data
     varname = "nspace"
     call bpwritescalarint(diag_ptl_ifile,nspace,varname)
     varname = "coordinates"
     call bpbeginbasicgroup(diag_ptl_ifile,varname)
     call bpendbasicgroup(diag_ptl_ifile)
! begin ion time history data
     call bpbegintimestepgroup(diag_ptl_ifile)
     if (sml_electron_on==1) then
        ptlfile="xgc.electron.bp"
        call bpopenfile(ptlfile,diag_ptl_efile)
! dummy mesh data
        varname = "nspace"
        call bpwritescalarint(diag_ptl_efile,nspace,varname)
        varname = "coordinates"
        call bpbeginbasicgroup(diag_ptl_efile,varname)
        call bpendbasicgroup(diag_ptl_efile)
! begin ion time history data
        call bpbegintimestepgroup(diag_ptl_efile)
     endif
     allocate(diag_ptl_times((diag_ptl_end-diag_ptl_begin)+1),&
              diag_ptl_data1(diag_ptl_num),diag_ptl_data2(diag_ptl_num))
  endif
#endif

  !-----------------------------------------------------
  call check_point2('efld_param')
  ! -1 : zero efield, no poisson eq
  !  0 : zero efield, solve poisson eq for diagnosis purpose
  !  1 : static efield, load efield profile from file 
  !  2 : self-consistent varying efield 
  ! -2 : static efield of modeled efeild. debug purpose. efld_mode will be set to 1 after initialization
  efld_mode=2   
  
  !static e-field paramter (mode 1) -- DIII-D parameter
  efld_1_psie=(1D0+(0.56105-0.559)/0.41)  ! x_psi unit -> code unit
  efld_1_dpsi=0.0082D0/0.41D0             ! x_psi unit -> code unit
  efld_1_psim=0.54/0.559                  ! x_psi unit -> code unit 
  efld_1_pot0=0.2D3                       ! V core potential -> code unit


  efld_pin=max(sml_inpsi,0.01d0)            ! x_psi unit -> code unit
  efld_pout=min(1.005D0,sml_outpsi)       ! x_psi unit -> code unit
  efld_dpdp0_outside=0D0                  ! Volt/m -> code unit
  efld_cutoff=1  ! enforce zero efield inside of some point
  efld_cutoff_psi=0.9*sml_inpsi + 0.1*sml_outpsi  ! x_psi unit -> code unit
  
  efld_reset=0
  efld_reset_step=2000
  efld_start=0  ! 0 - always on, N - start at N*dt time.

  efld_read_dpdp=0

  open(unit=14,file='input',action='read')
  read(14,nml=efld_param)
  close(14)

  !neutral ----------------------------------------------------------------
  
  !neutral collision-----------------------------------------------------------
  call check_point2('neu_param')
  ! neutral collision switch
  ! 1 : simple neutral collision, 0: no neutral collision 2 : consistent neutral calculation
  neu_col_mode=0  
  neu_adjust_n0=1         ! 0 fixed n0, 1 varying n0 - neu_recycle_rate should be specified.
  neu_ionize_mode=1       ! 1 - limiter hit recycling, 2 - sep flow recycleing
  neu_ionize2_psi=1D0      ! x_psi unit -> code unit
  neu_cx_period=300       !charge exchange period
  neu_ion_period=300      ! ionization freq. (period)
  neu_col_period=10       ! neutral elastic collision period
  neu_start_time=0*col_vb_period !Start the neutral effect after some background adjustment
  neu_recycle_rate=0.9    ! recycle rate of neutral
  neu_varying_mfp=1       ! 1: varying mean free path  0 : constant mean free path -- for mode 1 only
  neu_t0= eq_tempi_ev_out ! ev -> normalized unit, temperature -- for mode 1 only
  neu_n0=1.0D17           ! m^-3 -> nomalized unit , base neutral density
  neu_delta_n=  1.0D0     ! relative neut den over neu_n0 at poloidal peak
  neu_delta_theta=  sml_pi/6D0  ! poloidal peak witdh
  ! mean free path = v_0 * tau_ionize = sqrt(2 T_0) / (n_e0 * <sigma V> ionize )
  !neu_mfp0=2D0*sqrt(2D0 * neu_t0) / ( &
  !     eq_den_edge * 0.8D-8 *sqrt(0.75D0*eq_tempe_ev_edge*1D3) &
  !     * (exp(-13.56/eq_tempe_ev_edge/0.75D3)/(1D0+0.01D0*eq_tempe_ev_edge*0.75D3) ) * sml_norm_t &
  !     )
  neu_mfp0=0.05           ! cm -> code unit 5cm - initial mean free path 
  neu_theta_x= sml_2pi - acos(  (eq_x_r-eq_axis_r)/dsqrt((eq_x_r-eq_axis_r)**2+(eq_x_z-eq_axis_z)**2)  ) ! X-point theta
  neu_temp_factor=0.1
  ! neutral 2 only-------------------
  ! neutral 2 gets neutral temp and density consistently with boundary conditions
  neu_monte_num=10*sml_monte_num
  neu_mode2_period=col_vb_period
  neu_num=min(10,50/sml_totalpe + 1)
  neu_mpol=200/sml_totalpe + 1
  neu_temp0=3D0 !neutral temp. boundary condition  eV
  neu_istep_max=4000
  neu_sepfile=''

  neu_psi_edge=1.04
  neu_grid_max_psi=neu_psi_edge
  neu_grid_min_psi=0.9D0

  open(unit=14,file='input',action='read')
  read(14,nml=neu_param)
  close(14)

    
  ! limiter ---------------------------------------------------------------
  call check_point2('lim_param')
  lim_store_mz=40
  lim_filename=''
  open(unit=14,file='input',action='read')
  read(14,nml=lim_param)
  close(14)

  ! smooth -----------------------------------------------------------------
  call check_point2('smooth_param')
  smooth_mode_in=2
  smooth_n_in=30
  smooth_type_in=1  !gaussian type

  smooth_H_mode_in=2
  smooth_H_n_in=3
  smooth_H_type_in=1

  smooth_r1_n_in=0
  smooth_r1_d0_in=0.01
  smooth_r1_type_in=1

  open(unit=14,file='input',action='read')
  read(14,nml=smooth_param)
  close(14)

  smooth_diag_mode_in=smooth_mode_in
  smooth_diag_n_in=smooth_n_in
  smooth_diag_type_in=smooth_type_in

  open(unit=14,file='input',action='read')
  read(14,nml=smooth_param)
  close(14)


  ! turbulence -----------------------------------------------------------
  call check_point2('tbl_param')
  tbl_diffusion_on=0
  tbl_stop_time=sml_mstep
  tbl_period=100
  tbl_d_coeff=1D0
  open(unit=14,file='input',action='read')
  read(14,nml=tbl_param)
  close(14)



  ! heat ------------------------
  if(sml_heat_on) then
     call check_point2('heat_param')
     heat_period=100
     heat_inpsi=sml_inpsi
     heat_outpsi=sml_inpsi + 0.2*(sml_outpsi-sml_inpsi)
     heat_power=1D6 ! 1 MW
     open(unit=14,file='input',action='read')
     read(14,nml=heat_param)
     close(14)
  endif
  !***************************************************************************
  
  ! Enforce some input parameter 
  if(sml_push_mode==1) sml_pc_order=2


  

  ! Write out whole parameters
  if(sml_mype==0) then
     open(unit=15,file='fort.input.used',status='replace')
     write(15,nml=sml_param)
     write(15,nml=ptl_param)
     write(15,nml=eq_param)
     write(15,nml=col_param)
     write(15,nml=diag_param)
     write(15,nml=efld_param)
     !  write(15,nml=rf_param)
     write(15,nml=neu_param)
     write(15,nml=lim_param)
     write(15,nml=smooth_param)
     write(15,nml=tbl_param)
     write(15,nml=heat_param)
     close(15)
  endif
  ! 6.  Process  initializing with changing unit  

  ! Second Definition part (before memory allocation and interpolation)---
  
  if(sml_mype==0) print *,'ptl_num for each CPU =', ptl_num
  
  ! maximum number of particle per processor
  if(ptl_maxnum < ptl_num ) then
     print *, 'too small memory space. Increase ptl_maxnum in module.f90 or decrease ptl_num. '
     stop
  endif
  
  !append dir to filename.
  call add_dir_name2

  ! memory allocation
  call ptl_mem_allocation( ptl%ion,  1, ptl_maxnum  ,ptl_mass(1)  , ptl_charge(1),   sml_pc_order, sml_nlarmor, ptl_lost_nummax)
  call ptl_mem_allocation( ptl%elec, 2, ptl_e_maxnum,ptl_mass(2)  , ptl_charge(2),   sml_pc_order, 1          , ptl_lost_nummax)    
  call mem_allocation
  if(sml_mype==0) print *, 'init_interpolation'
  call init_interpolation
  if(sml_mype==0) print *, 'third difinition part'
  
  !-third definition part ----------------------------------------------
  if(sml_push_mode==2) then
     allocate(sml_pc_coeff(sml_pc_order,2))
     
     if(sml_pc_order==3) then
        sml_pc_coeff(:,1)=(/ 23D0/12D0 , -16D0/12D0, 5D0/12D0 /)
        sml_pc_coeff(:,2)=(/ 5D0/12D0 , 8D0/12D0, -1D0/12D0 /)
     elseif(sml_pc_order==4) then
        sml_pc_coeff(:,1)=(/ 55D0/24D0 , -59D0/24D0, 37D0/24D0, -9D0/24D0 /)
        sml_pc_coeff(:,2)=(/ 9D0/24D0 , 19D0/24D0, -5D0/24D0, 1D0/24D0 /)
     elseif(sml_pc_order==6) then
        sml_pc_coeff(:,1)=(/ 4277D0, -7923D0, 9982D0, -7298D0, 2877D0, -475D0 /) /1440D0
        sml_pc_coeff(:,2)=(/ 475D0, 1427D0, -798D0, 482D0, -173D0, 27D0 /) /1440D0
     else
        print *, 'Invalid pc order'
        stop
     endif
  endif
  ! energy order in MKS unit
  sml_en_order = sml_en_order_kev*1D3 * sml_e_charge 
  ! transit time for a particle (sp=1) at the mag axis with pitch1(??) , R=1 
  sml_tran = sml_2pi *eq_axis_r/ sqrt(2.D0 * sml_en_order/ptl_mass(1))
  ! delta t of one time step
  sml_dt= sml_dt*sml_tran ! Delta t - step size

  !simulation boundary
  sml_inpsi=sml_inpsi* eq_x_psi  !inner boundary of simulation region 
  sml_outpsi=sml_outpsi * eq_x_psi ! outter boundary of simulation region 
 
  sml_exb_suppress=0
  sml_exb_suppress_time=0.5 * sml_tran
  sml_exb_suppress_width=0.25 * sml_tran

  sml_f0den_ped_c= sml_f0den_ped_c * eq_x_psi   
  sml_f0den_ped_width=sml_f0den_ped_width * eq_x_psi

  sml_f0_1_Lt=sml_f0_1_Lt/eq_axis_r
  sml_f0_1_Ln=sml_f0_1_Ln/eq_axis_r

  sml_bd_Te_width=sml_bd_Te_width*eq_x_psi ! radial decay length of electron temperature at the boundary 
  !density and temperature : tanh function --------------------------
  ! edge : inside separatrix 
  ! out  : outside separatrix
  ! ped_c : center of tanh function
  ! ped_width : wdith of tanh function
  ! den - ion(electron) density
  eq_den_ped_c= eq_den_ped_c * eq_x_psi   
  eq_den_ped_width=eq_den_ped_width * eq_x_psi
  
  ! tempi - ion temperature
  eq_tempi_ped_c=eq_tempi_ped_c * eq_x_psi
  eq_tempi_ped_width=eq_tempi_ped_width * eq_x_psi
  
  !electron temperature
  eq_tempe_ped_c=eq_tempe_ped_c * eq_x_psi
  eq_tempe_ped_width=eq_tempe_ped_width * eq_x_psi


  !-----------------------------------------------------------------------
  
    
  !collision input------------------------------------------------------------
  col_2_pin=col_2_pin*eq_x_psi
  col_2_pout=col_2_pout*eq_x_psi
  col_2_dp=(col_2_pout-col_2_pin)/real(col_2_m)
  col_2_dtheta= sml_2pi/real(col_2_mtheta)
  
  col_vb_pin=col_vb_pin*eq_x_psi
  col_vb_pout=col_vb_pout*eq_x_psi
  col_vb_dp=(col_vb_pout-col_vb_pin)/real(col_vb_m-1)
  
  col_den_imp_ped_c = col_den_imp_ped_c * eq_x_psi
  col_den_imp_ped_width=col_den_imp_ped_width * eq_x_psi

  ! diagnosis ---------------------------------------------------------------------

  diag_flow_pin=diag_flow_pin * eq_x_psi
  diag_flow_pout=diag_flow_pout * eq_x_psi
  diag_flow_dp=(diag_flow_pout-diag_flow_pin)/real(diag_flow_npsi)

  diag_avg_pin=diag_avg_pin * eq_x_psi
  diag_avg_pout=diag_avg_pout * eq_x_psi
  diag_avg_dp= (diag_avg_pout-diag_avg_pin)/real(diag_avg_npsi-1)

  diag_f_pin=diag_f_pin * eq_x_psi
  diag_f_pout=diag_f_pout * eq_x_psi
  diag_f_dpsi= (diag_f_pout-diag_f_pin)/real(diag_f_npsi)
  diag_f_angle1=diag_f_angle1*sml_pi/180.0
  diag_f_angle2=diag_f_angle2*sml_pi/180.0
  diag_f_slope1=tan(diag_f_angle1)
  diag_f_slope2=tan(diag_f_angle2)

  ! efield------------------------------------------------------------------------

  efld_1_psie=efld_1_psie * eq_x_psi
  efld_1_dpsi=efld_1_dpsi * eq_x_psi !
  efld_1_psim=efld_1_psim * eq_x_psi ! 
  efld_1_potm=-2.0D0*efld_1_pot0 
  efld_1_psi1=efld_1_psim/sqrt(1D0-efld_1_potm/efld_1_pot0/(1D0-dexp(-((efld_1_psim-efld_1_psie)/efld_1_dpsi)**2)))

  !varying e-field parameter - global variable
  efld_dpdp=0D0 ! initialize
  efld_d2pdpdt=0D0 !initialize
  efld_t0=0D0

  efld_pin=efld_pin*eq_x_psi
  efld_pout=efld_pout*eq_x_psi
  efld_dp=(efld_pout-efld_pin)/real(efld_npsi-1)
  efld_dpdp0_outside=0.
  efld_pol_factor=8.8542D-12/(ptl_mass(1)) ! epsilon zero = 8.854D-12
  efld_cutoff_psi=efld_cutoff_psi * eq_x_psi
  efld_max_dpdp=5.0D5 ! sec^-1  - maximum dpot/dpsi
  !outside

  

  !rf : not implimented yet ----------------------------------------------
  rf_on=0
  rf_r_res=1.5
  rf_za=0D0
  rf_lap=1
  rf_e0= 1E5 

  !neutral collision-----------------------------------------------------------
  neu_ionize2_psi=neu_ionize2_psi*eq_x_psi
  ! mean free path = v_0 * tau_ionize = sqrt(2 T_0) / (n_e0 * <sigma V> ionize )
  !neu_mfp0=2D0*sqrt(2D0 * neu_t0) / ( &
  !     eq_den_edge * 0.8D-8 *sqrt(0.75D0*eq_tempe_ev_edge*1D3) &
  !     * (exp(-13.56/eq_tempe_ev_edge/0.75D3)/(1D0+0.01D0*eq_tempe_ev_edge*0.75D3) ) &
  !     )
  neu_mfp0=neu_mfp0 ! 5cm
  neu_mfp=neu_mfp0
  if(sml_mype==0 .and. neu_col_mode/=0) then
     print *, 'vth=', sqrt(2D0*neu_t0/ptl_mass(1))
     ! n_e0 should be 'col_denb_edge/2 , otherwise change_neu_mfp routine should be modified.
     print *,'mfp =' , neu_mfp
     print *, 'theta_x =' , neu_theta_x/sml_pi
  endif
  ! neutral 2 only-------------------
  ! neutral 2 gets neutral temp and density consistently with boundary conditions
  neu_dt=0.5D-4*sml_tran/sqrt(neu_temp0*1D-3) !heuristic
 

  neu_psi_edge=neu_psi_edge*eq_x_psi
  neu_grid_max_psi=neu_grid_max_psi*eq_x_psi
  neu_grid_min_psi=neu_grid_min_psi*eq_x_psi
  neu_grid_dtheta=sml_2pi/real(neu_grid_mtheta)
  neu_grid_dpsi=(neu_grid_max_psi-neu_grid_min_psi)/real(neu_grid_mpsi-1)
  neu_dtheta=sml_2pi/real(neu_mtheta)  



  ! turbulence -----------------------------------------------------------------
  

  tbl_d_coeff=tbl_d_coeff
  tbl_inpsi=sml_inpsi  !Should be greater than "efld_cutoff_psi"
  !Multiplication number for tbl_d_coeff, 10/23/2003
  tbl_mult_max=10
  tbl_mult_period=sml_mstep/tbl_mult_max
  allocate(tbl_mult(tbl_mult_max))
  tbl_mult(1:5)=1.d0
  tbl_mult(6:tbl_mult_max)=1D0  

  ! heating --------------------------------------------------------------------
  ! heating on/off
  !Multiplication number for heat_power, 09/23/2003
  heat_mult_max=10
  heat_mult_period=sml_mstep/heat_mult_max
  allocate(heat_mult(heat_mult_max))
  heat_mult(1)=1.d0
  heat_mult(2)=1.d0
  heat_mult(3)=1.d0
  heat_mult(4)=1.d0
  heat_mult(5)=1.d0
  heat_mult(6)=1.d0
  heat_mult(7)=1.d0
  heat_mult(8)=1.d0
  heat_mult(9)=1.d0
  heat_mult(10:heat_mult_max)=1.d0
  
  heat_inpsi = heat_inpsi*eq_x_psi
  heat_outpsi = heat_outpsi*eq_x_psi

  ! smooth --------------------------------------------------------------------



  ! ---------------------- others --------------------------------------------

  sml_bd_min_r=max(sml_bd_min_r,itp_min_r)
  sml_bd_max_r=min(sml_bd_max_r,itp_max_r)
  sml_bd_min_z=max(sml_bd_min_z,itp_min_z)
  sml_bd_max_z=min(sml_bd_max_z,itp_max_z)

  if(sml_mype==0) print *, 'R range :', sml_bd_min_r,sml_bd_max_r,&
       'Z range :',sml_bd_min_z,sml_bd_max_z
  bnc_min_r=sml_bd_min_r 
  bnc_max_r=sml_bd_max_r
  bnc_dr=(bnc_max_r-bnc_min_r)/real(bnc_nr-1)
  
  rpl_mode=0
  rpl_N_coil=24D0 ! D3D 24, real number for ripple strenght expresion
  rpl_R0= 1.14 ! D3D 1.14, TFTR 2.23 , ITER 6.75 - 0.00034*z**2
  rpl_ratio= 5.D-9  ! D3D 1D-9, TFTR 1.4E-5, ITER 3.75E-6
  rpl_elong=0.25 ! D3D 0.25, TFTR 1.1 , ITER 0.268 
  rpl_tau0 = 0.072 ! D3D 0.072, TFTR 0.18, ITER 0.535 meter
  

!-------------------------pre processing ------------------------

  sml_angle_stored=0
  
  col_2_dw_sum=0D0
  col_2_dvp_sum=0D0
  col_2_dv2_sum=0D0

  ptl_rz_outside=0
  
  ! initialize
  diag_avg_flux=0D0


! file description
#ifndef ADIOS
  if(sml_mype==0) then
     open(unit=22, file='fort.tracer', status='replace')
     write(22,'(a400)') '# [1]time [2]r [3]z [4]phi [5]rho_|| [6] E_kin(ev) [7] E_tot(ev) [8]pitch [9]psi [10]canonical P [11]mu [12]weight(df) [13-15] derivs [16] species'
     close(22)
     open(unit=71, file='fort.gam',status='replace')
     write(71,*) '# [1] Time(tran) [2] pot0 [3] dpot [4] idensity0 [5] edensity0'
     !close(71)
     open(unit=72, file='fort.fenergy',status='replace')
     write(72,*) '# [1] Time(tran) [2] Field Energy(J)'
     !close(72)
     write(30,*) '#[1] psi [2] density [3] flow_toro [4] flow_polo [5] flow_para [6] flow_ExB_pol [7] flow_raidal'

     if(col_mode==2)write(40,*) '#[1] psi [2] col_2_dw_sum  [3] col_2_dvp_sum  [4] col_2_dv2_sum '
     write(50,*) 
     
     write(39,*) '# [1]psi [2] dpdp  [3] E-field at midplane' ! E-field diagnosis
     write(300,*) '# special simulation 3 : [1] energy [2] pitch [3] mass [4] region [5] ptl_index'
     write(301,*) '# specical 3 - reamined particle : [1] K_|| [2] K_perp [3] energy [4] pitch [5] ptl_index'
     write(299,*) '# specical 3 - lost particel : [1] K_|| [2] K_perp [3] energy [4] pitch [5] ptl_index'
  
     write(103,*) '#bounce setup : '  
     write(10,*) 'sml_totalpe', sml_totalpe
     write(10,*) 'sml_dt(s)=', sml_dt
     write(10,*) 'sml_dt(tau)=',sml_dt/sml_tran
     write(10,*) 'sml_dt(normalized)=',sml_dt
     write(10,*) 'sml_tran(s)=',sml_tran
     write(10,*) 'ptl_num(total)=',ptl_num*sml_totalpe
     write(10,*) 'eq_axis_r(m)=', eq_axis_r
     write(10,*) 'eq_axis_z(m)=', eq_axis_z
     write(10,*) 'eq_axis_B(T)=', eq_axis_b
     write(10,*) 'psi_x(code unit)=', eq_x_psi
     write(10,*) 'psi_x(MKS)=', eq_x_psi
!     write(10,*) 'neu_mean_free_path(m)', neu_mfp
!     write(10,*) 'neu_den (mks)', neu_n0 
!     write(10,run_parameters)     
     close(10)
     open(unit=987,file='fort.pot_zero',status='replace')
     write(987,*) '# r, phi00, den00, grid%psi(grid%itheta0(i))'
     close(987)
  endif
#endif
  ! bounce setup 
!  ptl_nbounce=0D0
  if(sml_concentric/=1 .and. sml_bounce/=0) then
     if(sml_mype==0) print *, 'bounce setup'
     call bounce_setup    
  endif

  !limiter setup 
  if(sml_limiter==1) then  
     if(sml_mype==0) print *, 'limiter setup'
     call limiter_setup 
  endif 

  !read initial e-field
  if(efld_read_dpdp==1) then
     efld_dpdp=1D-3
  endif

  !col_vb initialize
  allocate(col_vb_vol(col_vb_m), col_vb_den(col_vb_m), col_vb_temp(col_vb_m))
  allocate( diag_f_vcs(diag_f_npsi) ,diag_f_n0(diag_f_npsi), diag_f_vt(diag_f_npsi))
  allocate( diag_f_dvol(diag_f_npsi))
  diag_f_count=0
  diag_f=0D0
  diag_f0=0D0
  
  do i = 1, diag_f_npsi
     psi=diag_f_pin+(i-0.5)*diag_f_dpsi
!     diag_f_dvol(i) = ! infinitesimal volume 
     ! set critical slowing down velocity
     ! vcs = v_the*((3*PI^0.5*m_e*Z_i)/(4*m_i))^(1/3)

     v_the= sqrt(2D0* tempe_ev(psi)*sml_ev2j*ptl_mass(1)/ptl_mass(2) )
     diag_f_vcs(i) =  v_the * ( (3.*sml_sqrtpi*ptl_charge(1))/(4.*ptl_mass(1)) )**(1D0/3D0)
       
     diag_f_n0(i) = init_den(psi)! set density
!     diag_f_vt(i) = sqrt( 2D0* init_tempi_ev(psi) *sml_ev2j/ptl_mass(1))! set thermal velocity
     diag_f_vt(i) = sqrt( 2D0* eq_tempi_ev_edge *sml_ev2j/ptl_mass(1))! set thermal velocity
!     diag_f_vt(i) = sqrt(2D0 * sml_en_order /ptl_mass(1)) ! set 200 eV
  enddo
  
!  neu_weight_sum_lost=0D0
!  if(neu_col_mode/=0 .and. (sml_concentric/=1 .or.sml_read/=0) ) then
!     open(unit=112, file=neu_sepfile, action='read')
!     do i=1, neu_sep_mtheta
!        read(112,*) neu_sep_r(i),neu_sep_z(i)
!     enddo
!     close(112)
!  endif
!1000 format(e19.13, 1x, e19.13)
  neu_weight_sum_lost=0D0
  if(neu_col_mode==1 ) then
     do i=1, neu_sep_mtheta
        !! incompete  
        ! Don't use neu_mode=1 until it is fixed
        ! neu_sep_r should be obtained from linear interpolation of neu_sep_*_file
        neu_sep_r(i)=neu_sep_r_file(1)
        neu_sep_z(i)=neu_sep_z_file(1)
     enddo
  endif

  if(sml_mype==0) print *, 'get_mid_r setup'
  call get_mid_r_setup
  if(sml_mype==0) print *, 'diag_efld setup'
  call diag_efld_setup  ! store efld_dpdr_mid value

  if(neu_col_mode == 2 ) then
     if(sml_mype==0) print *, 'neutral2_setup'
     call neutral2_setup
  endif


! Smooth initialization
  if(sml_tavg_on==1) then
     smooth00%mode=0
     smooth00%n=1 ! default
     smooth00%type=1 ! default     
  else
     smooth00%mode=smooth_mode_in
     smooth00%n=smooth_n_in
     smooth00%type=smooth_type_in
  endif
  call smooth_init(smooth00)

  !need for sml_tavg_on/=0 and diagnosis
  smooth0L%mode=smooth_mode_in
  smooth0L%n=smooth_n_in
  smooth0L%type=smooth_type_in
  call smooth_init(smooth0L)

  ! for high mode smoothing
  smoothH%mode=smooth_H_mode_in
  smoothH%n=smooth_H_n_in
  smoothH%type=smooth_H_type_in
  call smooth_init(smoothH)

  ! for diagnosis smoothing
  smoothdiag%mode=smooth_diag_mode_in
  smoothdiag%n=smooth_diag_n_in
  smoothdiag%type=smooth_diag_type_in
  call smooth_init(smoothdiag)

  ! setup restore temperature module
  if(sml_restore_temp) call  retore_temp_setup(rtempi,25,sml_inpsi,sml_outpsi)

  !-------------CHECKING INPUT VAILIDITY----------------------------------------------
  if(sml_totalpe*ptl_maxnum > huge(ptl%ion%maxgid)) then
     if(sml_mype==0) print *, 'You need 64-bit GID variable for sml_total_pe(',sml_totalpe,') X ptl_maxnum(',ptl_maxnum,').',huge(ptl%ion%maxgid)
     stop
  endif


end subroutine setup

subroutine mem_allocation
  use rf_module
  use ptl_module
  use itp_module
  use sml_module
  use eq_module, only : eq_mpsi, eq_mr, eq_mz
  implicit none


  ! ptl module allocation
  allocate(rf_r_save(ptl_maxnum))

  ! itp module allocation
  itp_mpsi=eq_mpsi
  itp_mr=eq_mr
  itp_mz=eq_mz
  
  allocate(itp_psi_knot(itp_mpsi+itp_korder), itp_I_cscoef(4,itp_mpsi),itp_I_break(itp_mpsi))
  allocate(itp_r_knot(itp_mr+itp_korder_rz),itp_z_knot(itp_mz+itp_korder_rz), &
       itp_psi_bscoef(itp_mr*itp_mz))


end subroutine mem_allocation
  


subroutine get_mid_r_setup
  use sml_module
  use eq_module
  use itp_module, only : itp_max_r
  implicit none
  real (kind=8) :: r1,r2,r0,rend,dr
  integer :: i,j,find
  integer, parameter :: NN=500
  real(kind=8), external :: psi_interpol
  real (kind=8) :: sqrt_out, sqrt_in, sqrt_p1, sqrt_p2, sqrt_dest

  sqrt_out=sqrt(sml_outpsi)
  sqrt_in=sqrt(sml_inpsi)
  eq_mid_r_dp= (sqrt_out-sqrt_in)/real(eq_mid_r_npsi-1)
  r0=eq_axis_r
  rend=itp_max_r
  dr= (itp_max_r-r0)/real( NN )
  j=1

  do i=1, eq_mid_r_npsi
     sqrt_dest= real(i-1) * eq_mid_r_dp + sqrt_in
     find=0
     ! find proper j ( R-index )
     j=1
     do while( find/=1 )
        r1= real(j-1) * dr + r0
        r2= real(j)* dr + r0
        sqrt_p1=sqrt(psi_interpol(r1,eq_axis_z,0,0))
        sqrt_p2=sqrt(psi_interpol(r2,eq_axis_z,0,0))
!        if( sqrt_p1 < sqrt_dest .and. sqrt_dest <= sqrt_p2 ) then
!           find=1
!        else if ( sqrt_dest <= sqrt_p1 ) then
!           print *, 'Error in get_mid_R_setup'
!           stop
!           find=1
        if( sqrt_dest < sqrt_p2) then
           find=1
        else
           j=j+1
        endif
     enddo
     eq_mid_r_psi(i)= ((sqrt_dest-sqrt_p1)*r2 + (sqrt_p2-sqrt_dest)*r1) / (sqrt_p2-sqrt_p1)
     eq_mid_r_psi(i)=max(eq_mid_r_psi(i),r0)
!     print *, i, eq_mid_r_psi(i), r2, r1, psi1, psi_dest, psi2
  end do
  
end subroutine get_mid_r_setup

real (kind=8) function get_mid_r(psi)
  use eq_module
  use sml_module
  implicit none
  real  (kind=8),intent(in) :: psi
  real (kind=8) :: aa,bb, sqrt_psi, sqrt_in
  integer :: i
  
  if( psi < sml_inpsi .or. psi > sml_outpsi) print * , 'Warning : psi in get_mid_r is out of interpolation range', psi/eq_x_psi,psi 

  sqrt_psi=sqrt(psi)
  sqrt_in=sqrt(sml_inpsi)
  

  i=int((sqrt_psi-sqrt_in)/eq_mid_r_dp)+1 
  i=min(eq_mid_r_npsi-1,max(1,i))
  bb= (sqrt_psi-sqrt_in)/eq_mid_r_dp + 1D0 - real(i)
  aa=1D0-bb

  get_mid_r = eq_mid_r_psi(i)*aa + eq_mid_r_psi(i+1) * bb

end function get_mid_r

subroutine diag_efld_setup
  use efld_module
  implicit none
  integer :: i
  real (kind=8) :: r,dpdr,dpdz,psi
  real (kind=8), external :: get_mid_R,psi_interpol

  do i=1, efld_npsi
     psi=efld_dp*real(i-1)+efld_pin
     r=get_mid_R(psi)
     dpdr=psi_interpol(r,0.d0,1,0)
     dpdz=psi_interpol(r,0.d0,0,1)
     efld_dpdr_mid(i)=sqrt(dpdr**2+dpdz**2)
!     print *, psi, r, dpdr, dpdz, efld_dpdr_mid(i)
  enddo
end subroutine diag_efld_setup

subroutine reset_nphi_total
  use sml_module
  
  if(sml_turb_efield/=1) then
     sml_nphi_total=1
  endif  

end subroutine reset_nphi_total

subroutine add_dir_name1
  use sml_module, only : sml_input_file_dir
  use eq_module, only : eq_filename

  implicit none
  integer :: endpos
  
  endpos=len_trim(sml_input_file_dir)
  eq_filename       = sml_input_file_dir(1:endpos) // eq_filename

!!$  print *, eq_filename

end subroutine add_dir_name1
subroutine add_dir_name2
  use sml_module, only : sml_node_file, sml_ele_file, sml_bfollow_file,&
       sml_input_file_dir, sml_add_pot0_file
  use eq_module, only : eq_filename
  use neu_module , only : neu_sepfile
  use lim_module, only : lim_filename
  implicit none
  integer :: endpos
  
  endpos=len_trim(sml_input_file_dir)
  sml_node_file     = sml_input_file_dir(1:endpos) // sml_node_file
  sml_ele_file      = sml_input_file_dir(1:endpos) // sml_ele_file
  sml_bfollow_file  = sml_input_file_dir(1:endpos) // sml_bfollow_file
  neu_sepfile       = sml_input_file_dir(1:endpos) // neu_sepfile
  lim_filename      = sml_input_file_dir(1:endpos) // lim_filename
  sml_add_pot0_file = sml_input_file_dir(1:endpos) // sml_add_pot0_file

!!$  print *, 'node',sml_node_file
!!$  print *, 'ele',sml_ele_file
!!$  print *, 'bf',sml_bfollow_file
!!$  print *, 'neu',neu_sepfile
!!$  print *, 'lim',lim_filename

end subroutine add_dir_name2

subroutine check_point2(str)
  use sml_module
  implicit none
  character (len=*) :: str
  integer :: n,ierr
  include 'mpif.h'

  call mpi_barrier(MPI_COMM_WORLD,ierr)  
  if(sml_mype==0) print *, str
end subroutine check_point2
