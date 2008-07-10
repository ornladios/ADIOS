!********************************************************
!! XGC 2 
!!@author Seunghoe Ku
!!@version 2.02
!********************************************************
program xgc2
#ifdef NETCDF
  use netcdf
#endif
  use sml_module
  use heat_module
  use tbl_module
  use neu_module
  use ptl_module
  use rf_module
  use col_module
  use efld_module
  use grid_class
  use psn_class
  use diag_module
  use smooth_module
  use rtemp_module
#if !defined(IMSL)
  use Ecuyer_random
#endif 
  use perf_monitor

  implicit none
  include 'mpif.h'
  
  integer :: istep,ierr,ipc, save_efld_mode,i
  character (len=10) :: ic(0:15)
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(ptl_type) ,target :: ptl
  interface
     subroutine restart_write(ptl)
       use ptl_module
       implicit none
       type(ptl_type),target :: ptl  
     end subroutine restart_write
#ifdef ELMFIRE_COOL
     subroutine elmfire_cool(ptl)
       use ptl_module
       implicit none
       type(ptl_type),target :: ptl
     end subroutine elmfire_cool
#endif
  end interface
 

  ! MPI initialize 
  call my_mpi_init 
#if !defined(NO_PETSC)
  call petsc_init( ierr )
#endif
  !ADIOS INITIALIZE
#ifdef ADIOS
  call adios_init ("config.xml"//char(0))!, MPI_COMM_WORLD, MPI_COMM_SELF, MPI_INFO_NULL)
#endif
  call init_perf_monitor( ierr )
  call monitor_start (TOTAL_)

  call monitor_start (INIT1_)
#if defined(IMSL)
  call rnset(sml_mype) !random number - seed initialize
#else
!  call init_seeds(sml_mype,sml_mype+1,sml_mype+2) !random number - seed initialize
!  call init_seeds(sml_mype,sml_mype+1,sml_mype+4) !random number - seed initialize
   call init_seeds(1234*sml_mype,2345*sml_mype+6789,4321*sml_mype+10) !random number - seed initialize
#endif

  if(sml_mype==0) print *, 'Total simulation processors : ',sml_totalpe
  ! Setup parameters and poisson array
  if(sml_mype==0) print *, 'setting up.........',sml_mype
  call setup(ptl)
  if(sml_mype==0) print *, 'setting up grid system.....',sml_mype
  call setup_grid(grid,psn)
#ifdef CIRCULAR_SPECIAL2
  call init_cc2_parameters(grid,psn)
#endif
  call smooth_r_init1(smooth_r1,smooth_r1_d0_in,smooth_r1_n_in,smooth_r1_type_in,grid%nnode)
  call smooth_r_init2(smooth_r1,grid)
  call init_f0(grid)
  call efield_init(grid,psn)
  call monitor_stop (INIT1_)

  call monitor_start (LOAD_)
  ! particle loading
  call load(grid,psn,ptl)
  if(sml_sheath_mode/=0) then
     if(sml_mype==0)  print *, 'sheath_pot_init',sml_mype
     call sheath_pot_init(grid,psn)
  endif
  call monitor_stop (LOAD_)
  call monitor_start (INIT2_)
  if(sml_mype==0) then
     print *, 'initial diagnosis'
     call output_initial
     call output_bfield
!     call tracer(diag_tracer_n,diag_tracer_sp,grid,psn,ptl)
!ADIOS ENABLE
!#ifdef BINPACK
! Dump bfield component data in binary file for post-processing and viz
     call dump_bfield(grid)
!#endif
#ifdef NETCDF
     print *, 'NETCDF (flow) init'
     call diag_flow_ncinit
     print *, 'NETCDF (avg) init'
     call diag_avg_ncinit
!        print *, 'NETCDF (efld) init'
!        call diag_efld_ncinit
#endif
!     call output_psi_der
  endif
  if(col_varying_bg==1 ) call col_snapshot(ptl)

  if(sml_electron_on/=1) call chargee_background(grid,psn,ptl)
  ! initial distribution ---------------------
  if(sml_relax_dist==1 .and. sml_electron_on==1) then
     save_efld_mode=efld_mode
     efld_mode=0
     do istep=1, sml_relax_dist_num
        do ipc=1,2
           call chargee(grid,psn,ptl%elec)
           call push(istep,ipc,grid,psn,ptl%elec)
        enddo
     enddo
     efld_mode=save_efld_mode !restore
     ptl%ion%phase(1:3,1:ptl%ion%num)=ptl%elec%phase(1:3,1:ptl%ion%num)
     ptl%ion%phase(6:9,1:ptl%ion%num)=ptl%elec%phase(6:9,1:ptl%ion%num)
     ptl%ion%gid(1:ptl%ion%num)=ptl%elec%gid(1:ptl%elec%num)
     !temp?
  endif
  !initialize local variable for poisson subroutine

  call poisson(grid,psn,0)

  call shift_ie(grid,ptl,0)
  call conserving_collision(ptl,0) ! if col_mode==2 ??
  ! main loop --------------------------
  if(sml_mype==0) print *, 'main loop started', sml_mype
  call monitor_stop (INIT2_)

  do istep=1 ,sml_mstep
     call monitor_start (MAIN_LOOP_)
     sml_istep=sml_istep+1
     if(sml_mype==0) print *, istep
     if(sml_exb_suppress==1) call set_exb_suppress(sml_time)
     do ipc=1,2 ! 1 for first RK2 step or Predicter, 2 for 2nd RK2 step or Corrector
        call monitor_start (IPC_LOOP_)
        sml_ipc=ipc

        call monitor_start (CHARGEI_)
!        print *, 'chargei', sml_mype
        call chargei(grid,psn,ptl%ion)
        call monitor_stop (CHARGEI_)
!        print *, 'chargee', sml_mype
        if(sml_electron_on==1) then
           call monitor_start (CHARGEE_)
           call chargee(grid,psn,ptl%elec)
           call monitor_stop (CHARGEE_)
        endif

!        print *, 'poisson', sml_mype
        call monitor_start (POISSON_)
        call poisson(grid,psn,sml_istep)
        call monitor_stop (POISSON_)

        call monitor_start (DERIVSI_)
        call phase_derivs(istep,ipc,grid,psn,ptl%ion)
        call monitor_stop (DERIVSI_)

        if(sml_electron_on==1) then
           call monitor_start (DERIVSE_)
           call phase_derivs(istep,ipc,grid,psn,ptl%elec)
           call monitor_stop (DERIVSE_)
        endif

!        print *, 'diagnosis', sml_mype        
        call monitor_start (DIAGNOSIS_)
        call diagnosis(sml_istep,ipc,grid,psn,ptl)
        if(ipc==1 .and. mod(istep-1,diag_pot_period)==0) then
           if(sml_mype==0) &
                print *, int(real(istep)/real(sml_mstep)*100.) ,'%'
        endif
        call monitor_stop (DIAGNOSIS_)

!        print *, 'pushi', sml_mype
        call monitor_start (PUSHI_)
        call push(istep,ipc,grid,psn,ptl%ion)
        call monitor_stop (PUSHI_)

!        print *, 'pushe', sml_mype
        if(sml_electron_on==1) then
           call monitor_start (PUSHE_)
           call push(istep,ipc,grid,psn,ptl%elec)
           call monitor_stop (PUSHE_)
        endif

        if(sml_special/=3 ) then
           call monitor_start (SHIFT_)
           call memory_cleaning(ptl)
#if defined(XGC_DEBUG) 
           call mem_clean_check(ptl)
#endif
           !        print *, 'shift_ie', sml_mype,ipc
           call shift_ie(grid,ptl,1)
#if defined(XGC_DEBUG)
           if( mod(sml_istep,1000)==0 ) call diagnosis_toroidal(grid,ptl,316)
#endif
           call monitor_stop (SHIFT_)
        endif
        
        call monitor_stop (IPC_LOOP_)
     enddo  !end ipc-loop


     !Calculate f0 for ion-ion collision 
!    if(sml_mype==0) print *, 'collision'
     call monitor_start (COLL_SNAP_)
     if(col_varying_bg==1 .and. col_mode/=0 ) then
        if(mod(sml_istep,col_vb_period)==0) then
           call col_snapshot(ptl) 
        endif
     endif
     call monitor_stop (COLL_SNAP_)

     if(sml_restore_temp .and. mod(sml_istep,sml_restore_temp_period)==0) then
        call restore_temp(ptl%ion,rtempi)
     endif

     ! particle-plasma collision
     call monitor_start (COLLISION_)
     call collision(sml_istep,ptl)
     call monitor_stop (COLLISION_)

     ! particle-neutral collision
     call monitor_start (NEUTRAL_)
     if(neu_col_mode/=0) then
        call neutral_col(sml_istep,ptl)
        !        else if(sml_injection_on==1) then 
        !           call injection(sml_istep)
     endif
     call monitor_stop (NEUTRAL_)

     ! Simply modeled turbulence diffusion - random walk
     call monitor_start (DIFFUSION_)
     if(tbl_diffusion_on==1) then
        call tbl_diffusion(sml_istep,ptl)
     endif
     call monitor_stop (DIFFUSION_)

     ! Inner boundary heating
     call monitor_start (HEATING_)
     if(sml_heat_on .AND. mod(sml_istep,heat_period)==0 ) then
        call heating(sml_istep,ptl)
     endif
     call monitor_stop (HEATING_)

#ifdef ELMFIRE_COOL
     call elmfire_cool(ptl)
#endif
    
     ! enforce zero efield for efld_mode=2
     if(efld_reset==1.AND.efld_reset_step>sml_istep) then
        efld_dpdp=0D0
        efld_d2pdpdt=0D0
        efld_reset=0
     endif

     sml_time=sml_time+sml_dt
! periodic restart dump
    if (mod(sml_istep,sml_restart_write_period)==0) then
       call monitor_start (RESTART_WRITE_)
       call restart_write(ptl)
       call monitor_stop (RESTART_WRITE_)
    endif 
    call monitor_stop (MAIN_LOOP_)
  enddo  ! end of main loop

  call monitor_start (FINALIZE_)    
#ifdef NETCDF
  if(sml_mype==0) then
     call diag_flow_ncfinal
     call diag_avg_ncfinal
  endif
#endif
  ! output for loss hole calcuation -- special simulation mode 3
  if(sml_special==3 .and. sml_mype==0) then
     call monitor_start (FINALIZE_SPEC3_)
     call special3_finalize
     call monitor_stop (FINALIZE_SPEC3_)
  endif
  
  if(sml_mype==0) then
     close(30)
     close(31)
     close(71)
     close(22)
#ifdef BINPACK
! finalize particle data binary dump files
      write(0,*)"dump finalized"

 !    call particle_dump_finalize(ptl)
#endif
  endif
  call monitor_stop (FINALIZE_)

  call monitor_stop (TOTAL_)
  call finish_perf_monitor(ierr)
#if !defined(NO_PETSC)  
  ! call delete_solver( solver, ierr )
#endif
  ! MPI_finalize
#ifdef ADIOS
  call adios_finalize (sml_mype)
#endif
  call my_mpi_finalize
  
  ! free memories
  call smooth_pol_delete(smooth0L)
  call smooth_pol_delete(smoothH)
  call smooth_pol_delete(smoothdiag)
  call smooth_r_delete(smooth_r1)
  ! free particle memories
  ! ----??

  ! free solver memories
  ! ----??

1000 format (8(e10.3))
1001 format (8(f9.3,'%'))
end program xgc2


!routines for debug-----------
!! NaN check
subroutine nan_check(tag,ptl)
  use ptl_module
  use sml_module
  implicit none
  type(ptl_type) :: ptl
  integer :: i,j
  real (kind=8) :: var
  character (len=30) :: tag

  do i=1, ptl%ion%num
     if(ptl%ion%gid(i)>0) then
        do j=1, ptl%ion%nphase
           var=ptl%ion%phase(j,i)
           if(var > 1D0 .or. var < 2D0 )  then
           else
              print *,'NaN found in',tag, sml_mype,i,j,' : ',ptl%ion%phase(:,i) 
              stop
           endif
        enddo
     endif
  enddo
end subroutine nan_check

!! subroutine for mu<0 check, only for debug
subroutine minus_mu_check(tag,ptl)
  use ptl_module
  use sml_module
  implicit none
  type(ptl_type) :: ptl
  integer :: i,j
  real (kind=8) :: var
  character (len=30) :: tag
  do i=1, ptl%ion%num
     if(ptl%ion%gid(i)>0) then
        var=ptl%ion%phase(5,i)
        if(var <0D0 )  then
           print *,'minus mu found ',tag, sml_mype,i,' : ',ptl%ion%phase(:,i) 
           print *,'-----------------'
           stop
        endif
     endif
  enddo
end subroutine 

subroutine setup_grid(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use diag_module, only : diag_gam_node_begin, diag_gam_node_num
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: ierr,i
  integer, parameter :: wall_num=100

  grid%guess_n=(/sml_guess_table_size,sml_guess_table_size/)
  grid%guess_min=(/sml_bd_min_r,sml_bd_min_z/)
  grid%guess_max=(/sml_bd_max_r,sml_bd_max_z/)
  if(sml_mype==0)  then
     print *,'Guess min', grid%guess_min
     print *,'Guess max', grid%guess_max
  endif
  
  if(sml_mype==0) print *, 'init_grid'
  call init_grid(grid,sml_node_file,sml_ele_file)
  if(sml_mype==0) print *, 'init_guess_table'
  call init_guess_table(grid)

     

  call psn_mem_alloc(psn,grid%nnode,grid%nphi,grid%ntriangle)
  if(sml_sheath_mode/=0) then
     !find wall point
     do i=1, grid%nnode
        if(grid%p(i)==wall_num) exit
     enddo
     psn%wall_start=i
     psn%nwall=grid%nnode-i+1
     call psn_sheath_alloc(psn)
  endif
  if(sml_turb_efield==1) then
     call init_bfollow(grid,psn)
  endif
  if(sml_mype==0) then 
     print *, 'init_poisson' 
     print *, 'original inner bd : ',grid%bd%in%start,grid%bd%in%end
     print *, 'original outer bd : ',grid%bd%out1%start,grid%bd%out1%end,grid%nnode
  endif
  
  call extend_boundary2(grid,psn)

!ADIOS ENABLE
!#ifdef BINPACK
! Dump grid description data in binary file for post-processing and viz
  if (sml_mype==0) call dump_grid(grid,psn)
!#endif
     
  ! set gam parameter
  if(diag_gam_node_begin<1) then
     diag_gam_node_begin=(psn%cbd0%in%end+min(psn%cbd0%out1%start,psn%cbd0%out2%start))/2-diag_gam_node_num/2
  endif

  call init_poisson(grid,psn,3)
  
  !  call get_matrix00(grid,psn)
#ifdef XGC_DEBUG1
  if(sml_mype==0) call output_matrix(grid,psn)
#endif
  
  if(sml_deltaf==1) call psn_deltaf_f0(grid,psn)
  
end subroutine setup_grid


!-------------------------------------------------------------------------------------
!!$subroutine extend_boundary(grid,psn)
!!$  use grid_class
!!$  use psn_class
!!$  use sml_module
!!$  use eq_module
!!$  implicit none
!!$  type(grid_type) :: grid
!!$  type(psn_type) :: psn
!!$  integer :: i,si
!!$  real (kind=8) , parameter :: epsilon=1D-8   
!!$  !  real (kind=8) ,parameter :: delta_psi=0.01D0, delta_psi2=0.02D0
!!$  ! ext_delta1 : psi margin for inner boundary for poisson equation -- pbd%in%end   ! inner_bd(2) 
!!$  ! ext_delta2 : psi margin for inner boundary for charge zero region -- cbd%in%end  !inner_bd(4)
!!$  ! ext_delta3 : psi margin for outer boundary for poisson equation -- pbd%out1%start, pbd%out2%start !outer_bd(1) , bd(3)
!!$  ! ext_delta4 : psi margin for outer boundary for charge zero region -- cbd%out1%start, cbd%out2%start ! outer_bd(5), bd(7)
!!$  
!!$
!!$
!!$  !*********** First, set H-solver boundary *******************************
!!$
!!$  ! set default boundary as intrinsic grid boundary
!!$  psn%cbdH=grid%bd
!!$  psn%pbdH=grid%bd
!!$
!!$  !  grid%outer_bd(4)=-1
!!$  psn%pbdH%out2%end=-1
!!$  !potential boundary setup
!!$  do si=1, grid%nnode     
!!$     i=si !grid%inv_sort(si)
!!$     if( grid%psi(i) - epsilon < sml_inpsi + sml_bd_ext_delta1*eq_x_psi) then  ! -epsilon is added for including magnetic axis node when small error exist at  (0. < 0). evaluation
!!$        !grid%inner_bd(2)=max(grid%inner_bd(2),si)
!!$        psn%pbdH%in%end=max(psn%pbdH%in%end,si) 
!!$     endif
!!$     if(sml_sheath_mode==0) then ! this valid only with no sheath simulation
!!$
!!$        ! Excluding private region
!!$        if(sml_exclude_private) then
!!$           if( grid%psi(i) > sml_outpsi - sml_bd_ext_delta3*eq_x_psi) then 
!!$!              grid%outer_bd(1)=min(grid%outer_bd(1),si)
!!$              psn%pbdH%out1%start=min(psn%pbdH%out1%start,si)
!!$           endif
!!$        else
!!$           ! Including private region
!!$           if( grid%psi(i) > sml_outpsi - sml_bd_ext_delta3*eq_x_psi ) then
!!$!              grid%outer_bd(3)=min(grid%outer_bd(3),si)
!!$              psn%pbdH%out2%start=min(psn%pbdH%out2%start,si)
!!$           endif
!!$           if( grid%rgn(i) == 2 .or. grid%rgn(i) == 1) then
!!$!              grid%outer_bd(4)= max(grid%outer_bd(4),si)   !maxium of region 1 and 2 node
!!$              psn%pbdH%out2%end = max(psn%pbdH%out2%end,si)   !maxium of region 1 and 2 node
!!$           endif
!!$        endif
!!$     endif
!!$  enddo
!!$
!!$  ! when there is no region3 area
!!$  if(psn%pbdH%out2%end<0) psn%pbdH%out2%end = psn%pbdH%out1%end 
!!$  ! end of potential boundary
!!$
!!$  ! charge boundary setup
!!$  !grid%outer_bd(8)=-1
!!$  psn%cbdH%out2%end=-1
!!$  do si=1, grid%nnode 
!!$     i=si !grid%inv_sort(si)
!!$     if( grid%psi(i) - epsilon < sml_inpsi + sml_bd_ext_delta2*eq_x_psi) then 
!!$        psn%cbdH%in%end=max(psn%cbdH%in%end,si)
!!$     endif
!!$     ! Excluding private region
!!$     if(sml_exclude_private) then
!!$        if( grid%psi(i)> sml_outpsi -sml_bd_ext_delta4*eq_x_psi) then 
!!$           psn%cbdH%out1%start=min(psn%cbdH%out1%start,si)
!!$        endif
!!$     else
!!$        ! Including private region
!!$        if( grid%psi(i) > sml_outpsi - sml_bd_ext_delta4*eq_x_psi ) then
!!$           psn%cbdH%out2%start=min(psn%cbdH%out2%start,si)
!!$        endif
!!$        if( grid%rgn(i) == 2 .or. grid%rgn(i) == 1) then
!!$           psn%cbdH%out2%end = max(psn%cbdH%out2%end,si)   !maxium of region 1 and 2 node
!!$        endif
!!$     endif
!!$  enddo
!!$  if(psn%cbdH%out2%end<0) psn%cbdH%out2%end=psn%cbdH%out1%end 
!!$
!!$
!!$
!!$  if(sml_mype==0) then
!!$!     print *, 'extended inner bd    : ',grid%inner_bd(1),grid%inner_bd(2)
!!$!     print *, 'extended inner bd-2  : ',grid%inner_bd(3),grid%inner_bd(4)
!!$!     print *, 'extended outer bd    : ',grid%outer_bd(1),grid%outer_bd(2),grid%nnode
!!$!     print *, 'extended outer bd2   : ',grid%outer_bd(3),grid%outer_bd(4),grid%nnode
!!$!     print *, 'extended outer bd-2  : ',grid%outer_bd(5),grid%outer_bd(6),grid%nnode
!!$!     print *, 'extended outer bd2-2 : ',grid%outer_bd(7),grid%outer_bd(8),grid%nnode
!!$
!!$     print *, 'extended inner bd (pot)     : ',psn%pbdH%in%start,psn%pbdH%in%end
!!$     print *, 'extended inner bd (charge)  : ',psn%cbdH%in%start,psn%cbdH%in%end
!!$     print *, 'extended outer bd  (pot)    : ',psn%pbdH%out1%start,psn%pbdH%out1%end
!!$     print *, 'extended outer bd2 (pot)    : ',psn%pbdH%out2%start,psn%pbdH%out2%end ,grid%nnode
!!$     print *, 'extended outer bd (charge)  : ',psn%cbdH%out1%start,psn%cbdH%out1%end ,grid%nnode
!!$     print *, 'extended outer bd2(charge)  : ',psn%cbdH%out2%start,psn%cbdH%out2%end ,grid%nnode
!!$  endif
!!$
!!$
!!$end subroutine extend_boundary

subroutine extend_boundary2(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8) :: in_bd_psi, out_bd_psi
  logical :: is_in, is_out, is_private
  real (kind=8), parameter :: epsil=0.005

  ! nonzero n-mode,  potential boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta1*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta3*eq_x_psi
  is_in=.true.
  is_out=.true.
  is_private=.false.
  if(sml_mype==0) print *, 'nonzero n-mode, potential boundary'
  call extend_bd_single(grid,psn%pbdH,in_bd_psi,out_bd_psi,is_in, is_out, is_private)

  ! nonzero n-mode, charge boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta2*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta4*eq_x_psi
  is_in=.true.
  is_out=.true.
  is_private=.false.
  if(sml_mype==0) print *, 'nonzero n-mode, charge boundary'
  call extend_bd_single(grid,psn%cbdH,in_bd_psi,out_bd_psi,is_in, is_out, is_private)

  ! n=0-mode, potential boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta1*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta3*eq_x_psi
  if(sml_rgn1_pot0_only) out_bd_psi=eq_x_psi !*(1D0-epsil)
  is_in=(sml_zero_inner_bd==1)
  is_out=(sml_sheath_mode<=0)
  is_private=.NOT. sml_exclude_private
  if(sml_mype==0) print *, 'n=0-mode, potential boundary'
  call extend_bd_single(grid,psn%pbd0,in_bd_psi,out_bd_psi,is_in, is_out, is_private)

  ! n=0-mode, charge boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta2*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta4*eq_x_psi
  if(sml_rgn1_pot0_only) out_bd_psi=eq_x_psi*(1D0-sml_bd_ext_delta4)
  is_in=.true.
  is_out=.true.
  is_private=.NOT. sml_exclude_private
  if(sml_mype==0) print *, 'n=0-mode, charge boundary'
  call extend_bd_single(grid,psn%cbd0,in_bd_psi,out_bd_psi,is_in, is_out, is_private)

end subroutine


subroutine extend_bd_single(grid,bd,in_bd_psi,out_bd_psi,is_in,is_out,is_private)
  use grid_class
  use sml_module
  use boundary_class
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(boundary_type) :: bd
  real (kind=8) :: in_bd_psi, out_bd_psi
  logical :: is_in, is_out, is_private
  integer :: i
  real (kind=8) , parameter :: epsilon=1D-8

  bd=grid%bd

  bd%out2%end=-1
  if(.NOT. is_in) bd%in%end=-1
  !potential boundary setup
  do i=1, grid%nnode
     if(is_in) then
        if( grid%psi(i) - epsilon < in_bd_psi) then  ! -epsilon is added for including magnetic axis node when small error exist at  (0. < 0). evaluation
           bd%in%end=max(bd%in%end,i)
        endif
     endif
     if(is_out) then 
        ! Excluding private region
        if(.NOT. is_private) then
           if( grid%psi(i) > out_bd_psi) then
              bd%out1%start=min(bd%out1%start,i)
           endif
        else
           ! Including private region
           if( grid%psi(i) > out_bd_psi ) then
              bd%out2%start=min(bd%out2%start,i)
           endif
           if( grid%rgn(i) == 2 .or. grid%rgn(i) == 1) then
              bd%out2%end = max(bd%out2%end,i)   !maxium of region 1 and 2 node
           endif
        endif
     endif
  enddo
  ! when there is no region3 area
  if(bd%out2%end<0) bd%out2%end = bd%out1%end

  if(sml_mype==0) then
     print *, 'extended inner bd  : ',bd%in%start,bd%in%end
     print *, 'extended outer bd (rgn2-rgn3)  : ',bd%out1%start,bd%out1%end
     print *, 'extended outer bd2(rgn3-wall)  : ',bd%out2%start,bd%out2%end ,grid%nnode
  endif

end subroutine
!-------------------------------------------------------------------------------------
subroutine monitor_start(i)
  use perf_monitor
#if defined(CAM_TIMERS)
  use perf_mod, only: t_barrierf, t_startf
#endif
  implicit none
  include 'mpif.h'

  integer, intent(in) :: i

  integer ierr

#if defined(CAM_TIMERS)
  if (i <= mon_N2) then
    if (mon_sync(i)) then
      call t_barrierf(mon_str(i+mon_NX), MPI_COMM_WORLD)
    endif

    call t_startf(mon_str(i))
  endif
#else
  if ((i <= mon_N) .and. (mon_sync(i))) call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

#if !defined(NO_PETSC)
  if (i <= mon_N) call PetscLogEventBegin( event(i), ierr )
#else
  if (i <= mon_N) call cpu_time( mon_time(i) )
#endif

end subroutine monitor_start

!-------------------------------------------------------------------------------------
subroutine monitor_stop(i)
  use perf_monitor
#if defined(CAM_TIMERS)
  use perf_mod, only: t_stopf
#endif
  implicit none
  integer, intent(in) :: i
  integer ierr
#if defined(NO_PETSC)
  real (kind=8) :: t2
#endif
  
#if !defined(NO_PETSC)
  if (i <= mon_N) call PetscLogEventEnd( event(i), ierr )
#else
  if (i <= mon_N) then
    call cpu_time(t2)
    mon_sum(i)=mon_sum(i)+t2-mon_time(i)
  endif
#endif

#if defined(CAM_TIMERS)
  if (i <= mon_N2) call t_stopf(mon_str(i))
#endif

end subroutine monitor_stop

!-------------------------------------------------------------------------------------
subroutine set_exb_suppress(time)
  use sml_module
  implicit none
  real (kind=8) :: time

!  sml_exb_sup = 1D0 + tanh( (time - sml_exb_suppress_time )/sml_exb_suppress_width )
  sml_exb_on(1)=1D0
  sml_exb_on(2)=0D0
end subroutine set_exb_suppress

!-------------------------------------------------------------------------------------
subroutine psn_deltaf_f0(grid,psn)
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i
  real (kind=8) :: z, psi
  real (kind=8), external :: f0_den
  
  allocate( psn%f0_density0(grid%nnode))
  
  do i=1, grid%nnode
     z=grid%x(2,i)
     psi=grid%psi(i)
     psn%f0_density0(i)=f0_den(psi,z) 
  enddo
  
end subroutine psn_deltaf_f0

subroutine new_communicator
   use sml_module
   implicit none
   include 'mpif.h'
   integer :: i,ierr, plane_0_pe
   integer :: MPI_GROUP_WORLD
   integer :: pe_ranks(0:sml_pe_per_plane-1)
#ifdef ADIOS
   integer :: plane_ranks(0:sml_totalpe/sml_pe_per_plane-1)
#endif
!   return
   plane_0_pe=int(sml_mype/sml_pe_per_plane)*sml_pe_per_plane
   do i=0, sml_pe_per_plane-1
      pe_ranks(i)=plane_0_pe+i
   enddo

  !get the group underlying MPI_COMM_WORLD
   call mpi_comm_group(MPI_COMM_WORLD, MPI_GROUP_WORLD,ierr)
  ! Create the new group
   call mpi_group_incl(MPI_GROUP_WORLD, sml_pe_per_plane, pe_ranks, sml_plane_group, ierr)
  ! Create the new communicator
   call mpi_comm_create(MPI_COMM_WORLD, sml_plane_group, sml_plane_comm,ierr)
   call mpi_comm_size(sml_plane_comm,sml_plane_totalpe,ierr)
   call mpi_comm_rank(sml_plane_comm,sml_plane_mype,ierr)
   sml_plane_index=sml_mype/sml_plane_totalpe

#ifdef ADIOS   
   plane_0_pe=int(sml_mype/(sml_totalpe/sml_pe_per_plane))*(sml_totalpe/sml_pe_per_plane)
   do i=0, (sml_totalpe/sml_pe_per_plane)-1
      plane_ranks(i)=plane_0_pe+i
   enddo
   call mpi_group_incl(MPI_GROUP_WORLD, sml_totalpe/sml_pe_per_plane, plane_ranks, sml_pe_group, ierr)
   call mpi_comm_create(MPI_COMM_WORLD, sml_pe_group, sml_pe_comm,ierr)
   call mpi_comm_rank(sml_pe_comm,sml_pe_mype,ierr)
#endif
   !print *,sml_pe_per_plane , plane_0_pe, sml_plane_mype, sml_totalpe, sml_mype ,sml_plane_index
   !stop

end subroutine

subroutine mem_clean_check(ptl)
  use ptl_module
  implicit none
  type(ptl_type) :: ptl
  integer :: i
  do i=1, ptl%ion%num
     if( ptl%ion%gid(i)<0) then
        print *, 'error-ion', ptl%ion%gid(i)
        call err_count
     endif
  enddo
  do i=1, ptl%elec%num
     if( ptl%elec%gid(i)<0) then
        print *, 'error-ion', ptl%elec%gid(i)
        call err_count
     endif
  enddo
end subroutine mem_clean_check

subroutine range_check(ptl)
  use sml_module
  use ptl_module
  use eq_module
  implicit none
  type(ptl_type) ,target :: ptl
  integer :: i,sp_type,rz_outside
  type(species_type) ,pointer :: sp
  real (kind=8) :: r,z
  do sp_type=1,1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif

     do i=1, sp%num
        r=sp%phase(1,i)
        z=sp%phase(2,i)
        rz_outside=0
        if(r<eq_min_r) then
           r=eq_min_r
           rz_outside=1
        else if (r>eq_max_r)then
           r=eq_max_r
           rz_outside=1
        endif
        if(z<eq_min_z) then
           z=eq_min_z
           rz_outside=1
        else if (z>eq_max_z)then
           z=eq_max_z
           rz_outside=1
        endif
     enddo
     if(rz_outside==1)  then
        print *, 'Outside',r,z
        call err_count
     endif
  enddo
end subroutine range_check

  
subroutine err_count
  implicit none
  integer :: count=0
  integer :: sml_max_error_msg_allow
  save count

  sml_max_error_msg_allow=3000
  count=count+1
  if(count>sml_max_error_msg_allow) then
     print *, '###############################################################'
     print *, '# error count is exceeded maximum number allowed',sml_max_error_msg_allow,count
     print *, '###############################################################'
     stop
  endif
end subroutine

#ifdef XGC_DEBUG_LOAD_GRID
subroutine load_grid(grid,psn,ptl)
  use ptl_module
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(ptl_type),target :: ptl
  type(psn_type) :: psn
  type(species_type),pointer :: sp 
  integer :: sp_type, i
  real (kind=8), parameter :: epsil=1D-7
  real (kind=8), external :: ranx

  do sp_type=1,1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif


     ! set particle number to grid number
     sp%num=grid%nnode
     if(sp%num >= sp%maxnum) then
        print *, 'increase maxnum', sp%num, sp%maxnum
        stop
     endif
     
     do i=1, grid%nnode
        ! set particle position 
        sp%phase(1,i)= grid%x(1,i) + epsil*ranx()
        sp%phase(2,i)= grid%x(2,i) + epsil*ranx()

        ! set particle velocity
        sp%phase(3,i)=(grid%phimin + grid%phimax)*0.5
        sp%phase(4,i)=0D0
        sp%phase(5,i)=0D0
        sp%phase(6:7,i)=0D0
        ! set particle weight
        sp%phase(8,i)=grid%node_vol(i)
     enddo

  enddo


end subroutine load_grid
#endif

#ifdef ELMFIRE_COOL
subroutine elmfire_cool(ptl)
  use ptl_module
  use sml_module
  use eq_module
  implicit none
  type(ptl_type),target :: ptl
  type(species_type),pointer :: sp 
  integer :: i
  real (kind=8) :: r_minor, dneu, r0
  real (kind=8) :: alpha ,prob
  real (kind=8) , external :: ranx
  dneu=0.01
  r0=0.5D0 * 0.358D0 * eq_axis_r  ! half minor radius
  alpha=1D0-1D-3

  sp=> ptl%ion
  do i=1, sp%num
     r_minor=sqrt( (sp%phase(1,i) - eq_axis_r)**2 + (sp%phase(2,i)-eq_axis_z)**2 )
     prob=1/cosh( (r_minor-r0)/dneu )**2
     if(ranx() < prob) then
        sp%phase(4,i)=alpha*sp%phase(4,i)
        sp%phase(5,i)=alpha**2 * sp%phase(5,i)
     endif
  enddo
end subroutine elmfire_cool
#endif
