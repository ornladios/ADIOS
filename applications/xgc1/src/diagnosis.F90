!*************************************************
! diagnosis main 
!************************************************
subroutine diagnosis(istep,irk,grid,psn,ptl)
  use sml_module
  use diag_module
  use efld_module
  use col_module
  use grid_class
  use psn_class
  use ptl_module
  implicit none
  integer, intent(in) :: istep,irk
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(ptl_type) :: ptl
  integer :: i
  interface
     subroutine tracer(n,sp_type,grid,psn,ptl)
       use ptl_module
       use grid_class
       use psn_class
       implicit none
       type(ptl_type) :: ptl
       type(grid_type) :: grid
       type(psn_type) :: psn
       integer, intent(in) :: n,sp_type
     end subroutine tracer
     subroutine flow_diagnosis(grid,psn,ptl,istep)
       use ptl_module
       use grid_class
       use psn_class
       implicit none
       type(grid_type) :: grid
       type(psn_type) :: psn
       type(ptl_type),target :: ptl       
       integer :: istep      
     end subroutine flow_diagnosis
     subroutine diagnosis_time_avg(istep,ptl)
       use ptl_module
       implicit none
       type(ptl_type),target :: ptl  
       integer, intent(in) :: istep
     end subroutine diagnosis_time_avg
     subroutine pweight_diagnosis(istep,ptl)
       use ptl_module
       implicit none
       type(ptl_type),target :: ptl  
       integer, intent(in) :: istep
     end subroutine pweight_diagnosis
!!ADIOS JC
!#ifdef BINPACK
     subroutine particle_dump(istep,ptl)
       use ptl_module
       implicit none
       type(ptl_type),target :: ptl  
       integer, intent(in) :: istep
     end subroutine particle_dump
!#endif
!!ADIOS JC
     subroutine diag_2d(istep,grid,psn,ptl)
       use grid_class
       use psn_class
       use ptl_module
       implicit none
       integer :: istep
       type(grid_type) :: grid
       type(psn_type) :: psn
       type(ptl_type),target :: ptl
     end subroutine diag_2d
     subroutine diag_gam(istep,grid,psn,ptl)
       use grid_class
       use psn_class
       use ptl_module
       implicit none
       integer :: istep
       type(grid_type) :: grid
       type(psn_type) :: psn
       type(ptl_type),target :: ptl
     end subroutine diag_gam
  end interface
  select case( irk )
  case(1) ! IRK = 1

     ! diagnosis sumation 
     if(diag_avg_on==1) then
        call diagnosis_time_avg(istep,ptl)
        if(mod(istep,diag_avg_outperiod)==0) then
           call time_avg_output(istep)
        endif
     endif

     ! f diagnosis 
     if(mod((istep-1)/diag_f_mod,diag_f_skip)==0 .and. diag_f_on==1) then
        call diagnosis_f(ptl,istep)
     endif
     
     ! tracer --------------------
     if(diag_tracer_n/=0) then
        if(mod(istep,diag_tracer_period)==0) then
           call tracer(diag_tracer_n,diag_tracer_sp,grid,psn,ptl)
        endif
     endif
     
#ifdef XGC_DEBUG4 
     !charging routine test - with single particle
     call charging_test(istep,grid,psn)     
#endif
  
     if(mod(istep,diag_pot_period)==0 .or. &
        mod(istep,diag_binout_period)==0) then
        if(sml_mype==0)write(0,*)"diag_pot"
        call diag_pot(istep,grid,psn)
        if(sml_mype==0)write(0,*)"diag_2d"
        call diag_2d(istep,grid,psn,ptl)
     endif
     if(diag_gam_on==1) then
	call diag_gam(istep,grid,psn,ptl)
     endif

     ! flow diagnosis
     if(mod(istep,diag_flow_period)==0) then
        call flow_diagnosis(grid,psn,ptl,istep)
     endif

     ! dump rms deviation of ion particle weight 
     if (mod(istep,diag_pw_period)==0 .and. diag_pw_on==1) then
        call pweight_diagnosis(istep,ptl)
     endif

!!ADIOS JC
!#ifdef BINPACK
     ! dump 2d fluid field and particle data in binpack format - redundant
!     if(mod(istep,diag_pot_period)==0 .or. &
!        mod(istep,diag_binout_period)==0) then
!        call diag_2d(istep,grid,psn,ptl)
!     endif
     if (istep>=diag_ptl_begin .and. istep<=diag_ptl_end) then
        if (diag_ptl_on==1) then
           call particle_dump(istep,ptl)
        endif
     endif
!#endif
!!ADIOS JC

  case(2) ! IRK = 2
     

  end select

end subroutine

subroutine diag_gam(istep,grid,psn,ptl)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(ptl_type), target :: ptl
  integer , parameter :: N=6, NT=17
  integer :: i,tmp
  integer, intent(in) :: istep
  integer :: iphi
  real (kind=8) :: sum_total(NT), suma(NT)
  integer :: counter=0, init=0
  integer :: sp_type
  real (kind=8) , allocatable :: b(:,:)
  real (kind=8) :: en, bptl,ephi
  real (kind=8), external :: b_interpol
  type(species_type), pointer :: sp
  save counter, init
#ifdef ADIOS
  real (kind=8) :: gam(100)
  real (kind=8) :: fenergy(NT+1)
  integer :: ncol_gam,ncol_fenergy
  integer*8:: grp_id,buf_id
  character(len=50) :: filename
#endif

#ifdef ADIOS

  if(sml_mype==0) then
!!NEW OPEN
!     call adios_get_group(grp_id,"fort.gam"//char(0))
     filename="fort.bp"//char(0)
!     call adios_open_append(buf_id,grp_id,filename)
     call adios_open(buf_id,"fort.gam"//char(0),filename,"a"//char(0))
!     call adios_set_path(grp_id,"gam"//char(0))
     call adios_set_path(buf_id,"gam"//char(0))
  endif 
#endif
  counter=mod(counter+1,10)
  if(sml_mype==0) then    

#ifdef ADIOS
     gam(1)=istep*sml_dt/sml_tran
     do i=0,diag_gam_node_num-1
        gam(2+4*i)=psn%pot0(diag_gam_node_begin+i)
        gam(3+4*i)=psn%dpot(diag_gam_node_begin+i,0)
        gam(4+4*i)=psn%idensity0(diag_gam_node_begin+i)
        gam(5+4*i)=psn%edensity0(diag_gam_node_begin+i)
     enddo
     ncol_gam=diag_gam_node_num*4+1
#else
     tmp=1+3*diag_gam_node_num

     write(71,1000) istep*sml_dt/sml_tran,&
          (psn%pot0(diag_gam_node_begin+i),&
          psn%dpot(diag_gam_node_begin+i,0),& 
          psn%idensity0(diag_gam_node_begin+i),&
          psn%edensity0(diag_gam_node_begin+i),&
          i=0,diag_gam_node_num-1)
     
     if(counter==0) then
        close(71)
        open(unit=71, file='fort.gam', position='append')
     endif
#endif
  endif
  !****end of gam diagnosis***********8

  !************************* field energy diagnosis***************************************
  
!!$  if(init==0) then ! initialization for B-field at node point
!!$     allocate(b(4,grid%nnode))
!!$     do i=1, grid%nnode        
!!$        call bvec_interpol(grid%x(1,i),grid%x(2,i),0D0,b(1,i),b(2,i),b(3,i))
!!$        b(4,i)=sqrt(b(1,i)**2+b(2,i)**2+b(3,i)**2)
!!$     enddo
!!$     init=1
!!$  endif
  suma(:)=0D0

  do iphi=1, grid%nphi
     do i=1, grid%nnode
        suma(1)=suma(1)+psn%E_para(i,iphi)**2*grid%node_vol(i)
        suma(2)=suma(2)+(psn%E_perp_node(1,i,iphi)**2+psn%E_perp_node(2,i,iphi)**2)*grid%node_vol(i)
        Ephi=(psn%E_para(i,iphi)*grid%bfield(4,i) &
             - psn%E_perp_node(1,i,iphi)*grid%bfield(1,i) &
             - psn%E_perp_node(2,i,iphi)*grid%bfield(2,i) &
             )/grid%bfield(3,i)
        suma(2)=suma(2)+Ephi**2*grid%node_vol(i)
        suma(3)=suma(3)+(psn%idensity(i,iphi)**2)*grid%node_vol(i)
        suma(4)=suma(4)+psn%dpot(i,iphi)**2*grid%node_vol(i)
        suma(5)=suma(5)+(psn%dpot(i,iphi)+psn%pot0(i))**2*grid%node_vol(i)
     enddo
     do i=1, grid%ntriangle
        suma(1)=suma(1)+(psn%E_perp_tr(1,i,iphi)**2 + psn%E_perp_tr(2,i,iphi)**2)*grid%tr_vol(i)
     enddo     
  enddo

  ! particle weight diagnosis 
  do sp_type=1, 1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif
        
     do i=1, sp%num
        suma(6+N*sml_electron_on)=suma(6+N*sml_electron_on) + sp%phase(6,i)*sp%phase(8,i)  ! deltaf weight sum
        suma(7+N*sml_electron_on)=suma(7+N*sml_electron_on) + (sp%phase(6,i)*sp%phase(8,i))**2  ! deltaf weight**2
        suma(8+N*sml_electron_on)=suma(8+N*sml_electron_on) + sp%phase(8,i)              ! full-f weight sum        
        suma(9+N*sml_electron_on)=suma(9+N*sml_electron_on) + sp%phase(6,i)**2
        bptl=b_interpol(sp%phase(1,i),sp%phase(2,i),sp%phase(3,i))
        en=(sp%phase(5,i)+ptl_c2_2m(sp_type)*sp%phase(4,i)**2*bptl)*bptl
        suma(11+N*sml_electron_on)=suma(11+N*sml_electron_on) + en*sp%phase(6,i)*sp%phase(8,i)
     enddo
     suma(10+N*sml_electron_on)=real(sp%num)
  enddo
  
  call my_mpi_reduce(suma,sum_total,NT)
  
  if(sml_mype==0) then
#ifdef ADIOS
     ncol_fenergy=6+(1+sml_electron_on)*N
     fenergy(1)=istep*sml_dt/sml_tran
     fenergy(2)=0.5D0*8.8541878D-12*sum_total(1)/real(sml_pe_per_plane)
     fenergy(3)=0.5D0*8.8541878D-12*sum_total(2)/real(sml_pe_per_plane)
     fenergy(4)=sum_total(3)/real(sml_pe_per_plane)
     fenergy(5)=sum_total(4)/real(sml_pe_per_plane)
     fenergy(6)=sum_total(5)/real(sml_pe_per_plane)
     fenergy(7:ncol_fenergy)=sum_total(6:5+(1+sml_electron_on)*N)
     call adios_gwrite(buf_id,"fort.gam")
     call adios_close(buf_id)
#else
     write(72,1000) istep*sml_dt/sml_tran,&
          0.5D0*8.8541878D-12*sum_total(1)/real(sml_pe_per_plane),&
          0.5D0*8.8541878D-12*sum_total(2)/real(sml_pe_per_plane),&
          sum_total(3)/real(sml_pe_per_plane),&
          sum_total(4)/real(sml_pe_per_plane),&
          sum_total(5)/real(sml_pe_per_plane),&
          (sum_total(i),i=6,5+(1+sml_electron_on)*N)
     if(counter==0) then
        close(72)
        open(unit=72, file='fort.fenergy', position='append')
     endif
#endif
  endif
  
1000 format(100(e19.13,1x))
end subroutine diag_gam
  
subroutine diag_pot(istep,grid,psn)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  implicit none
  include 'mpif.h'

  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,istep,iphi,nspace=2
  character (len=30) :: filename, filename2
#ifdef ADIOS 
  integer*8 :: buf_id, grp_id
#endif
  real (kind=8) , allocatable :: potential(:), pot3d(:,:)
  integer   :: nnode,nphi

  if (sml_mype==0) then
     iphi=1
     if (mod(istep,diag_pot_period)==0) then
#ifdef ASCII
        write(filename,'("fort.p",".",i2.2)') (istep)/(diag_pot_period)
        open(500, file=filename, status='replace', form='formatted') 
        do i=1, grid%nnode
           write(500,1000) grid%x(1,i),grid%x(2,i) ,&
             psn%idensity(i,iphi),psn%dpot(i,iphi),&
             psn%idensity0(i),psn%edensity0(i),&
             psn%density0_full(i),psn%pot0(i),&
#ifdef XGC_DEBUG10
             grid%node_vol(i)*real(sml_nphi_total),&
             grid%rtmp3(i), grid%rtmp4(i),&
             grid%rtmp5(i), grid%rtmp6(i),&
             grid%rtmp7(i), grid%rtmp8(i)
#else
	     grid%node_vol(i)*real(sml_nphi_total)
#endif
        enddo
        close(500)
      write(filename,'("timestep_",".",i2.2)') (istep)/(diag_pot_period)
      filename=trim(filename)//char(0)
!! NEW OPEN
!      call adios_get_group(grp_id,"fort.p"//char(0))
!      call adios_set_path(grp_id,filename) 
!      call adios_open_append(buf_id,grp_id,"diag_pot.p.bp"//char(0))
      call adios_open(buf_id,"fort.p"//char(0),"diag_pot.p.bp"//char(0),"a"//char(0))
      call adios_set_path(buf_id,filename)
      call adios_gwrite(buf_id,"fort.p") 
      call adios_close(buf_id)
#endif

     endif !istep%diag_pot_period
  endif !if(sml_mype==0)
  if (sml_mype==0) then
     iphi=1
#ifdef ADIOS
! NEW OPEN
!     call adios_get_group(grp_id,"fort.ef"//char(0))
     write(filename,'("dir_",i2.2,"_",i2.2)') sml_mype/sml_pe_per_plane,(istep)/(diag_pot_period)
     filename=trim(filename)//char(0)
!     call adios_set_path(grp_id,filename)
!     call adios_open_append(buf_id,grp_id,"diag_pot.ef.bp"//char(0))
     call adios_open(buf_id,"fort.ef"//char(0),"diag_pot.ef.bp"//char(0),"a"//char(0))
     call adios_set_path(buf_id,filename)
     call adios_gwrite(buf_id,"fort.ef") 
     call adios_close(buf_id)
#else
     open(500, file=filename, status='replace', form='formatted') 
     do i=1, grid%nnode
        write(500,1000) grid%x(1,i),grid%x(2,i),&
             psn%dpot(i,iphi),&
             psn%E_para(i,iphi)
     enddo
     close(500)
#endif
  endif

1000 format(19(e19.13,' '))

#ifdef ADIOS
     nnode=grid%nnode
     nphi=grid%nphi
     if(sml_pe_per_plane==1)then
        allocate(pot3d(grid%nnode,grid%nphi))
        pot3d(:,1:grid%nphi) = psn%dpot(:,1:grid%nphi)
     else
        allocate(pot3d(grid%nnode,1))
        nphi=1
        pot3d(:,1) = psn%dpot(:,1)
     endif
     write(filename2,'("diag_pot_fieldp",".",i4.4,".bp")') (istep)/(diag_binout_period)
     filename2=trim(filename2)//char(0)
!! NEW OPEN
!     call adios_get_group(grp_id,"diag_pot.fieldp"//char(0))
     call adios_open(buf_id,"diag_pot.fieldp"//char(0),filename2,"w"//char(0))

     write(filename2,'("node_data[",i0,"]")')sml_pe_mype+1 
     filename2=trim(filename2)//char(0)
!     call adios_set_path(grp_id,filename2)
     call adios_set_path(buf_id,filename2)
     call adios_gwrite(buf_id,"diag_pot.fieldp")
     call adios_close(buf_id)
     deallocate(pot3d)
     if (sml_mype==0) then
         allocate(potential(grid%nnode))
         potential(:) = psn%pot0(:)  ! normalized potential (Volts)
         write(filename2,'("fieldp",".",i4.4,".bp")') (istep)/(diag_binout_period)
         filename2=trim(filename2)//char(0)
!! NEW OPEN
!         call adios_get_group(grp_id,"diag_pot.fieldp.header"//char(0))
!         call adios_open_append(buf_id,"diag_pot.fieldp.header"//char(0),filename2,"a"//char(0))
         call adios_open(buf_id,"diag_pot.fieldp.header"//char(0),filename2,"a"//char(0))
         call adios_gwrite(buf_id,"diag_pot.fieldp.header")
         call adios_close(buf_id)
         deallocate(potential)
     endif
  !endif
#endif
end subroutine diag_pot

subroutine diagnosis_toroidal(grid,ptl, fnum)
  use sml_module
  use grid_class
  use ptl_module
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(ptl_type),target :: ptl
  type(species_type),pointer :: sp
  integer :: fnum
  integer :: iphi,i,sp_type,ierror
  integer :: num(sml_nphi_total),dum(sml_nphi_total),num2(sml_nphi_total,2)
 
!!$  if(sml_mype==0)then       
!!$     do i=1, sml_nphi_total
!!$        num(i)=i
!!$     enddo
!!$     print *, 'numbers/ions/electrons :'
!!$     print 1000, num(1:sml_nphi_total)
!!$  endif
!!$  
!!$  
  num2=0
  do sp_type=1, 1+sml_electron_on
     num=0
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif
     
     do i=1, sp%num 
        iphi= FLOOR(sp%phase(3,i)/grid%delta_phi) +1
!        if(iphi>grid%nphi+grid%iphi_offset .or. iphi < grid%iphi_offset) then
!           print *, 'error in diagnosis_toroidal', iphi, phase(i,3),sp,sml_mype,fnum
!        else
        if(iphi<=0 .or. iphi>sml_nphi_total) print *, iphi,sp%phase(3,i),grid%delta_phi
           num(iphi)=num(iphi)+1
!        endif
     enddo
!     call my_mpi_reduce(num,dum,sml_nphi_total)
     call monitor_start (DIAG_TOR_RED_)
     call mpi_reduce(num,dum,sml_nphi_total,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
     call monitor_stop (DIAG_TOR_RED_)

     if(sml_mype==0) then
        !        print 1000, dum(1:sml_nphi_total) 
        num2(1:sml_nphi_total,sp_type)=dum(1:sml_nphi_total)
     endif
  enddo
  
  if (sml_mype==0) then
     do i=1, sml_nphi_total
        write(fnum,*) i, num2(i,1),num2(i,2)
     enddo
     write(fnum,*) ' '
  end if
1000 format (64(I5))

end subroutine diagnosis_toroidal

!weight, momentum, energy change due to collision
subroutine col_2_diagnosis
  use col_module
  implicit none
  integer :: i,l
#ifdef ADIOS
  integer*8 :: grp_id, buf_id
  integer, allocatable :: idx(:,:)
  real(kind=8), allocatable :: psi(:,:)
#else
  real (kind=8) :: psi
#endif

#ifdef ADIOS
  allocate(idx(col_2_m,col_2_dp),psi(col_2_m,col_2_dp)) 
#endif

  do i=1, col_2_m
#ifdef ADIOS
     psi(i,:)= col_2_dp*(real(i) - 0.5) + col_2_pin
     do l=1, col_2_mtheta
        idx(i,:)=l
     enddo
#else
     psi= col_2_dp*(real(i) - 0.5) + col_2_pin
     do l=1, col_2_mtheta
        write(40,*) psi, l, col_2_dw_sum(i,l),col_2_dvp_sum(i,l),col_2_dv2_sum(i,l)
     enddo
     write(40,*) ' '
#endif
  enddo

#ifdef ADIOS
!! NEW OPEN
!  call adios_get_group(grp_id,"fort.40"//char(0))
!  call adios_open_append(buf_id,grp_id,"fort.bp"//char(0))
  call adios_open(buf_id,"fort.40"//char(0),"fort.bp"//char(0),"a"//char(0))
  call adios_set_path(grp_id,"fort.40"//char(0))
  call adios_gwrite(buf_id,"fort.40")
  call adios_close(buf_id)
  deallocate (idx,psi)
#else
   write(40,*) ' '
#endif
end subroutine col_2_diagnosis

!*********************************************************
! single particle tracer
!*********************************************************
subroutine tracer(n,sp_type,grid,psn,ptl)
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(ptl_type) :: ptl
  type(fld_type) :: fld    
  integer, intent(in) :: n,sp_type
  type(species_type) , pointer :: sp
  real (kind=8) :: r,z,phi,rho,mu,b,psi,psi_c,b_interpol,ev,pitch,br,bz,bphi
  real (kind=8) :: exb1, exb2
#ifdef SMALL_GID
  integer :: gid
#else
  integer (kind=8) :: gid
#endif
  target ptl
#ifdef ADIOS
  real (kind=8) :: tracerdata(18)
  integer :: ncol_tracer
  integer*8:: grp_id,buf_id
  character(len=50) :: filename
#endif

  if(sp_type==1) then
     sp=>ptl%ion
  else
     sp=>ptl%elec
  endif
  
  gid=sp%gid(n)
  if (gid>0) then
     fld%r=sp%phase(1,n)
     fld%z=sp%phase(2,n)
     fld%phi=sp%phase(3,n)
     fld%psi=sp%phase(9,n)

     call field(fld,sml_time)
     call efield(grid,psn,sp,n,fld,sml_time)
     if (sml_deltaf==1) call f0_info(grid,n,sp,fld)
     r=fld%r
     z=fld%z
     phi=fld%phi
     rho=sp%phase(4,n)
     mu=sp%phase(5,n)
     psi=fld%psi
     
     b=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
     bphi=fld%bphi
     br=fld%br
     bz=fld%bz
     
     call rho_mu_to_ev_pitch2(rho,mu,b,ev,pitch,sp_type)
     if(sml_minusB/=1) then
        psi_c=psi + rho*bphi*r
     else
        psi_c=psi - rho*bphi*r
     endif
     
     exb1= (fld%bz*fld%Ephi - fld%Bphi * fld%Ez)/B/B
     exb2= (fld%bphi * fld%Er - fld%br * fld%Ephi )/B/B  
#ifdef ADIOS
!! NEW OPEN
!     call adios_get_group(grp_id,"fort.tracer"//char(0))
     filename="fort.bp"//char(0)
!     call adios_open_append(buf_id,grp_id,filename)
     call adios_open(buf_id,"fort.tracer"//char(0),filename,"a"//char(0))
!    call adios_set_path(grp_id,"tracer"//char(0))
     call adios_set_path(buf_id,"tracer"//char(0))
      tracerdata(1)=sml_time/sml_tran
      tracerdata(2)=r
      tracerdata(3)=z
      tracerdata(4)=phi
      tracerdata(5)=rho
      tracerdata(6)=ev
      tracerdata(7)=ev+fld%Epot
      tracerdata(8)=pitch
      tracerdata(9)=psi/eq_x_psi
      tracerdata(10)=psi_c
      tracerdata(11)=mu
      tracerdata(12)=sp%phase(6,n)
      tracerdata(13)=sp%derivs(1,sp%pc_index(1),n)
      tracerdata(14)=sp%derivs(2,sp%pc_index(1),n)
      tracerdata(15)=sp%derivs(4,sp%pc_index(1),n)
      tracerdata(16)=real(gid)
      tracerdata(17)=(exb1*fld%dpsidr + exb2*fld%dpsidz)/sqrt(fld%dpsidr**2+fld%dpsidz**2)
      tracerdata(18)=(sp%derivs(1,sp%pc_index(1),n)*fld%Er+sp%derivs(2,sp%pc_index(1),n)*fld%Ez+sp%derivs(3,sp%pc_index(1),n)*fld%r*fld%Ephi)/fld%f0_temp
     !  write(22,*) sml_time,ev,psi,psi_c
     ncol_tracer=18

     call adios_gwrite(buf_id,"fort.tracer")
     call adios_close(buf_id)
#else
     open(unit=22, file='fort.tracer', position='append')
     write(22,1000) sml_time/sml_tran,r, z, phi, rho,&
          ev,ev+fld%Epot,pitch,psi/eq_x_psi,psi_c,mu,sp%phase(6,n), & 
          sp%derivs(1,sp%pc_index(1),n),sp%derivs(2,sp%pc_index(1),n),sp%derivs(4,sp%pc_index(1),n),& 
          real(gid), (exb1*fld%dpsidr + exb2*fld%dpsidz)/sqrt(fld%dpsidr**2+fld%dpsidz**2), &
          (sp%derivs(1,sp%pc_index(1),n)*fld%Er+&
          sp%derivs(2,sp%pc_index(1),n)*fld%Ez+&
          sp%derivs(3,sp%pc_index(1),n)*fld%r*fld%Ephi)/fld%f0_temp
     !  write(22,*) sml_time,ev,psi,psi_c
     close(22)
#endif
  endif
1000 format(18(e19.13,1x))
end subroutine tracer

!**********************************************************
! special3 finalize 
!**********************************************************
subroutine special3_finalize(ptl)
  use sml_module
  use ptl_module
  implicit none
  type(ptl_type):: ptl
  real (kind=8) :: en,pitch,v_perp,v_para
  integer :: i,mass,region
 
#ifndef ADIOS
  do i=1, ptl%ion%num
    ! sp=ptl_species(i)
     mass=ptl%ion%mass
     !     region=isign(1,ptl%ion%gid(i))
     if(ptl%ion%gid(i)> 0) then
        region=1
     else
        region=-1
     endif
     en=ptl_special_en_ev*real(i/ptl_special_n)/real(ptl_num/ptl_special_n)
     pitch=-1D0 + 2* real(mod(i,ptl_special_n))/real(ptl_special_n-1)
     write(300,*) en,pitch, mass,region,i
     write(300+region,*) en*(pitch**2)*dsign(1D0,pitch), en*(1.-pitch**2),en,pitch,i
  enddo
#endif
end subroutine special3_finalize


subroutine flow_diagnosis(grid,psn,ptl,istep)
#ifdef ADIOS
  use eq_module, only : eq_x_psi
#endif
  use psn_class ! XGC-1
  use ptl_module
  use eq_module
  use sml_module
  use efld_module
  use diag_module
  use fld_module
  use grid_class
  use perf_monitor
  implicit none
  type(psn_type) :: psn
  type(grid_type) :: grid
  type(fld_type) :: fld
  type(ptl_type), target :: ptl
  type(species_type), pointer :: sp

  integer,parameter :: NN=diag_flow_nvars1 ! When you change this number, change array size of diag_flow_ncvarid in module.F90
  integer,parameter :: NN2=diag_flow_nvars2 ! When you change this number, change array size of diag_flow_ncvarid in module.F90
  integer :: i,j,k,sp_type,istep !,sp
  real (kind=8) :: b, bp, rho, ww(2), y(ptl_nphase_max), psi, time_now
  real (kind=8) :: flow(diag_flow_nvars1,diag_flow_npsi,2),weight(diag_flow_npsi,2), &
       mean(diag_flow_nvars1,2) , v(diag_flow_nvars1),dum(diag_flow_nvars1,diag_flow_npsi,2),wdum(diag_flow_npsi,2)
  real (kind=8) :: en_sum(2,diag_flow_npsi,2),en(2),dum2(2,diag_flow_npsi,2),temp(2,2)
!  real (kind=8) :: count_in(diag_flow_npsi)
  real (kind=8) ::  yprime(ptl_nphase_max) ,yprime2(ptl_nphase_max)
  real (kind=8), external :: f0_den, f0_temp_ev
  real (kind=8) :: c_m, c2_2m
  real (kind=8) :: out_var(diag_flow_nvars2) ! When you change array size of this, change array size of diag_flow_ncvarid in module.F90
  !m3d coupling
  integer, parameter :: N_m3d=9
  integer :: itemp, numqr
  integer, parameter :: numqrmax=300
  real (kind=8) :: m3d_couple(N_m3d,diag_flow_npsi),norm_pot,get_mid_r,I_interpol
  real (kind=8) :: inipsi, endpsi, temppsi,a1e,a1egrad,a1efield
  real (kind=8) :: electemp(diag_flow_npsi),k13(diag_flow_npsi),k23(diag_flow_npsi),&
               Eradial(diag_flow_npsi),epspar(diag_flow_npsi),poloB(diag_flow_npsi),&
               elecden(diag_flow_npsi),electempgrad(diag_flow_npsi),kc(diag_flow_npsi), &
               elecdengrad(diag_flow_npsi),paraflow(diag_flow_npsi),totB(diag_flow_npsi),&
               iontemp(diag_flow_npsi),nuste(diag_flow_npsi),nusti(diag_flow_npsi),&
               elecdenraw(diag_flow_npsi),electempraw(diag_flow_npsi)
  real (kind=8) :: elecbootj(diag_flow_npsi),ionbootj(diag_flow_npsi) 
  real (Kind=8) :: neoparaflow(diag_flow_npsi),elecbootej(diag_flow_npsi), &
                   elecbootgradj(diag_flow_npsi),iontempgrad(diag_flow_npsi), &
                   neoparagradflow(diag_flow_npsi),neoparaeflow(diag_flow_npsi), &
                   neoionbootj(diag_flow_npsi),neobootj(diag_flow_npsi)
  real (kind=8) :: qprofile(diag_flow_npsi),qprofileout(numqrmax),qpsidata(numqrmax), &  
                   qrdata(numqrmax),psigrid(diag_flow_npsi),figrid(diag_flow_npsi), &
                   potgrid(diag_flow_npsi),paraBflow(diag_flow_npsi),&
                   potgradgrid(diag_flow_npsi)
  real (kind=8) :: dpdr
  real (kind=8) :: B2_sum(diag_flow_npsi),dumB2(diag_flow_npsi)
  real (kind=8) :: ffderiv(diag_flow_npsi),pderiv(diag_flow_npsi)
  real (kind=8) :: aveJdotB_EFIT(diag_flow_npsi)
  real (kind=8) :: r, z
  real (kind=8) , external :: psi_interpol,dinit_den_dpsi , init_den
  character (len=26) :: filename
  character (len=32) :: var_names(diag_flow_nvars2)
  save  B2_sum, aveJdotB_EFIT
  
#ifdef ADIOS
  real (kind=8) :: simtime(1)
  real (kind=8) :: netcdf_flow(diag_flow_npsi,NN2)
  character (len=5) :: sp_name(2)=(/'ion__', 'elec_'/)
  character (len=30) ::filename2
  integer*8 :: buf_id,grp_id
  integer*8 :: buf_id_fort,grp_id_fort
  real (kind=8), dimension(diag_flow_npsi) :: psi1
  var_names(1) ='psi'
  var_names(2) ='density(df)'
  var_names(3) ='toroidal_flow(df)'
  var_names(4) ='poloidal_flow(df)'
  var_names(5) ='parallel_flow(df)'
  var_names(6) ='poloidal_comp._of_ExB(df)'
  var_names(7) ='Radial_flow_times_grad_psi(df)'
  var_names(8) ='Radial_current_density(df)'
  var_names(9) ='v_para_times_B(df)'
  var_names(10)='parallel_mean_energy(df)'
  var_names(11)='perp._temperature(df)'
  var_names(12)='parallel_temperature(df)'
  var_names(13)='density'
  var_names(14)='toroidal_flow'
  var_names(15)='poloidal_flow'
  var_names(16)='parallel_flow'
  var_names(17)='poloidal_comp._of_ExB'
  var_names(18)='Radial_flow_times_grad_psi'
  var_names(19)='Radial_current_density'
  var_names(20)='v_para_times_B'
  var_names(21)='parallel_mean_energy'
  var_names(22)='perp._temperature'
  var_names(23)='parallel_temperature'
  var_names(24)='delta-f_ratio'
  var_names(25)='flow_shell_volume'
#endif 

  do sp_type=1, 1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif
     c_m=ptl_c_m(sp_type) 
     c2_2m=ptl_c2_2m(sp_type) 

     flow=0D0
     en_sum=0D0
     weight=1D-20
     if (istep==0) then 
        B2_sum=0D0
     endif

     time_now=sml_time
#ifdef ADIOS
     call adios_open(buf_id_fort,"fort.30"//char(0),"fort.bp"//char(0),"a"//char(0))
     write(filename2,'("/istep_",i5.5)')istep
     filename2=trim(filename2)//char(0)
     call adios_set_path(buf_id_fort,"fort.30"//filename2)
#else
     if(sml_mype==0) then
        write(30+sp_type-1,*) '#',diag_flow_npsi, sml_time/sml_tran ! insert empty line for gnuplot  - time separate
        write(30+sp_type-1,*)
     endif
#endif     
     do i=1, sp%num
        y(:)=sp%phase(:,i)
        if(sp%gid(i)>0) then
           fld%r=y(1)
           fld%z=y(2)
           fld%phi=y(3)
           fld%psi=y(9)
           call field(fld,time_now)
           ! obtaining E-field for each particle 
           call efield(grid,psn,sp,i,fld,time_now) !! diff. from XGC-0
           ! get time derivatives
           call diagnosis_derivs(fld,time_now,y,yprime,yprime2,sp_type)
           yprime=sp%derivs(:,sp%pc_index(1),i) !! diff. from XGC-0

           if(diag_flow_pin<fld%psi .AND.fld%psi<diag_flow_pout .AND. fld%z>eq_x_z) then
              rho=y(4)
              j=int((fld%psi-diag_flow_pin)/diag_flow_dp)+1
              j=min(diag_flow_npsi,max(1,j))

              bp=sqrt(fld%bz**2+fld%br**2) + 1D-30 ! B_p strength , 1D-30 added avoid divide by zero
              B=sqrt(fld%bz**2+fld%br**2+fld%bphi**2) ! B-field strength

              v(1)=yprime(3)*fld%r ! toroidal velocity
              v(2)=(yprime(1)*fld%br+yprime(2)*fld%bz)/bp ! poloidal velocity
              v(3)=c_m*rho*B ! v_|| parallel velocity
              v(4)=(yprime2(1)*fld%br+yprime2(2)*fld%bz)/bp  ! poloidal component of ExB
              v(5)=(yprime(1)*fld%bz-yprime(2)*fld%br)*fld%r   !radial component of velocity * grad psi
              v(6)=sqrt(fld%bz**2+fld%br**2)*fld%r           ! for radial current density
              v(7)=c_m*rho*(B**2) ! v_|| B 
              
              !for temperature calculation
              en(1)= c2_2m* (rho*B)**2
              en(2)= y(5)*B
              
              if(sml_minusb==1) v(5)=-v(5)

              if(sml_deltaf==1) then
                 ww(1)=sp%phase(8,i)*sp%phase(6,i)
                 ww(2)=sp%phase(8,i)
                 !              count_in(j,sp)=count_in(j,sp)+1D0
              else
                 ww(1:2)=sp%phase(8,i)             
              endif
              do k=1, 2
                 flow(:,j,k)=v(:)*ww(k)+flow(:,j,k)
                 en_sum(:,j,k)=en(:)*ww(k)+en_sum(:,j,k)
                 weight(j,k)=ww(k) + weight(j,k)
              enddo
              if (istep==0) B2_sum(j)=(B**2)*ww(1) + B2_sum(j)
           endif
        endif !if(sp%gid(i)>0) then
     enddo !do i=1, sp%num

     call monitor_start (FLOW_DIAG_RED_)
     !MPI reduce using my_mpi_reduce
     !  call my_mpi_reduce(count_in,wdum,diag_flow_npsi*ptl_sp_num)
     !  count_in=wdum
     call my_mpi_reduce(flow,dum,NN*diag_flow_npsi*2)
     flow=dum
     call my_mpi_reduce(en_sum,dum2,2*diag_flow_npsi*2)
     en_sum=dum2
     call my_mpi_reduce(weight,wdum,diag_flow_npsi*2)
     if (istep==0) then
        call my_mpi_reduce(B2_sum,dumB2,diag_flow_npsi)
        B2_sum=dumB2
     endif
     call monitor_stop(FLOW_DIAG_RED_)

     if (sml_mype==0 .and. (mod(istep,diag_flow_period)==0)) then  
        if(sml_deltaf==1) then
           do i=1,diag_flow_npsi
              psi=diag_flow_dp*(real(i)-0.5) + diag_flow_pin
              weight(i,2)=diag_flow_vol(i)*init_den(psi)
              en_sum(1,i,1)=0.5D0*weight(i,1)*f0_temp_ev(psi,0D0)*sml_ev2j
              en_sum(2,i,1)=2D0*en_sum(1,i,1)
           enddo
           weight(:,1)=weight(:,1) + wdum(:,1)
           en_sum(:,:,1)=en_sum(:,:,1) + dum2(:,:,1)
           weight(:,2)=wdum(:,2)
           en_sum(:,:,2)=dum2(:,:,2)
        else
           weight(:,:)=wdum(:,:)
           en_sum(:,:,:)=dum2(:,:,:)
        endif
        do i=1, diag_flow_npsi
           psi=diag_flow_pin+diag_flow_dp*(real(i)-0.5D0)
           if(weight(i,1)==0) then 
              mean=0D0
           else
              do k=1,2
                 mean(:,k)=flow(:,i,k)/weight(i,k)
                 temp(:,k)=en_sum(:,i,k)/weight(i,k)
              enddo
           endif
           do k=1,2
              mean(5,k)=mean(5,k)*weight(i,k)/(diag_flow_vol(i))* &
                   sp%charge*sml_e_charge/(mean(6,k)+1D-40)
           enddo
           
           !output variables
           out_var(1)         =psi/eq_x_psi ! normalized psi
           out_var(2)         =weight(i,1)/(diag_flow_vol(i))  ! density - delta-f
           out_var(3:NN+2)    = mean(1:NN,1)    ! flows - delta-f
           out_var(NN+2)      =out_var(NN+2)                ! v_para*B 
           out_var(NN+3:NN+4) =temp(1:2,1)*sml_j2ev*1D-3        ! mean energy - delta-f
           out_var(NN+5)      = 2D0*(temp(1,1) - 0.5D0*ptl_mass(sp_type)*mean(3,1)**2)*sml_j2ev*1D-3 ! temp - delta-f
           out_var(NN+6)      = weight(i,2)/(diag_flow_vol(i)) ! density - full-f
           out_var(NN+7:2*NN+6) = mean(1:NN,2)  ! flows - full-f
           out_var(2*NN+6)    =out_var(2*NN+6)             ! v_para*B
           out_var(2*NN+7:2*NN+8)=temp(1:2,2)*sml_j2ev*1D-3    ! mean energy - full-f
           out_var(2*NN+9)    =2D0*(temp(1,2) - 0.5D0*ptl_mass(sp_type)*mean(3,2)**2)*sml_j2ev*1D-3 ! temp - full-f
           out_var(2*NN+10)   =wdum(i,1)/weight(i,1)
           out_var(NN2)   =diag_flow_vol(i)           
#ifdef ADIOS 
           netcdf_flow(i,:)=out_var
#else
           write(30+sp_type-1,2000) (out_var(j),j=1,NN2)
#endif
        enddo !do i=1, diag_flow_npsi
#ifdef ADIOS
       call adios_gwrite(buf_id_fort,"fort.30")
       call adios_close(buf_id_fort)
#endif
#ifdef ADIOS

     if(istep==diag_flow_period) then
        do j=1, diag_flow_npsi
           psi1(j) = (diag_flow_pin + diag_flow_dp * (real(j)-0.5D0)) / eq_x_psi
        enddo
     endif
     
     !diag_time_ncstart(1) = istep/diag_flow_period + 1

     filename2="flowdiag.bp"//char(0)
     call adios_open(buf_id,"diagnosis.flow"//char(0),filename2,"a"//char(0))
     call adios_gwrite(buf_id,"diagnosis.flow")
     call adios_close(buf_id)

#endif
     endif !if (sml_mype==0 .and. (mod(istep,diag_flow_period)==0)) then  
  enddo !965: do sp_type=1, 1+sml_electron_on

2000 format(100(e19.13,' '))
2100 format(30(e19.13,' '))
3000 format(5(e16.9,' '))  
3100 format(3(e16.9,' '))
end subroutine flow_diagnosis


!!****************************************************************************
!!> derivatives for diagnosis
!!
!!<***************************************************************************
subroutine diagnosis_derivs(fld,T, Y,YPRIME,yprime2,sp_type)
    use fld_module
    use sml_module
    use ptl_module
    use rpl_module
    implicit none
    type(fld_type), intent(in) :: fld
    real (kind=8), intent(in) :: T, Y(ptl_nphase_max) 
    real (kind=8), intent(out) :: YPRIME(ptl_nphase_max),YPRIME2(ptl_nphase_max)
    integer, intent(in) :: sp_type !! sp = 2 for elec, sp =1 for ion
    real (kind=8):: D, B, nb_curl_nb, dbdr, dbdz, fr,fp,fz,r,z,phi,rho,mu,dbdphi
    real (kind=8):: ripp,dripp_dr,dripp_dz
    integer :: i
    real (kind=8) :: mass, charge,c_m  !relative charge and mass to normalizing charge and mass for multi-species
    !! for weigth calculation
    real (kind=8) :: exb1,exb2,ddendpsi,dtempdpsi,energy,temp,den,gradn1,gradn2,&
         gradt1,gradt2,exbgradn, exbgradt, exbgradb
    !! variables for optimization. 1/B, 1/B^2, rho^2, mu*rho^2*B, B^2
    real (kind=8):: over_B, over_B2,cmrho2,cmrho,murho2b,murho2b_c,b2
    !! initial ion temp, initial ion density, dT/dpsi, dn/dpsi
    real (kind=8), external :: init_tempi_ev, init_den
    
    

    mass=ptl_mass(sp_type)
    charge=ptl_charge(sp_type) !-1.D0/ptl_charge
    c_m=charge/mass
 
    r=y(1)
    z=y(2)
    phi=y(3)
    rho=y(4)
    mu=y(5)    

    B = sqrt( fld%br**2 + fld%bz**2 + fld%bphi **2 )
    b2=b**2
    over_B=1/B
    over_B2=1/B2
       
    ! normalized b dot curl of normalized b
    nb_curl_nb= 1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - (fld%bphi/r  + fld%dbpdr)*fld%bz )
    D=1.D0/ ( 1.D0 + rho * nb_curl_nb )
       
    dbdphi=0D0
!    if(RPL_mode==1) then  ! don't use with delta-f , not implimented yet
!       call get_ripple(r,z,ripp,dripp_dr,dripp_dz)
!       dbdphi=B*rpl_N_coil*dsin(rpl_N_coil*phi)*ripp
!    endif
       
    dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
    dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
    ! bug fix -- rho**2 B grad B term added 2002/1/15
    !! optimization variables
    cmrho2=c_m*rho**2
    cmrho =c_m*rho
    murho2b=(mu+charge*cmrho2 *B)
    murho2b_c=murho2b/charge
    fr =  fld%Er - (murho2B_c) * dbdr ! fr is force over q . -gradient of  Hamiltonian
    fp =  fld%Ephi - (murho2B_c) * dbdphi/r  ! modified by shlee 5/30/2001 -- /r added
    fz =  fld%Ez - (murho2B_c) * dbdz

    ! F is grad H/q


    yprime(1)= D*( (fld%bz*Fp - fld%Bphi * Fz) * over_B2         &
         +  cmrho * fld%br                       &
         +  cmrho2 * (fld%dbzdp/r - fld%dbpdz ) )
    yprime(2)= D*( (fld%bphi * fr - fld%br * fp ) * over_B2      &
         +  cmrho * fld%bz                       &
         +  cmrho2 * (fld%bphi/r + fld%dbpdr-fld%dbrdp/r) )
    yprime(3)= D*( (fld%br * fz - fld%bz * fr) * over_B2         &
         +  cmrho * fld%bphi                     &
         +  cmrho2 * ( fld%dbrdz - fld%dbzdr) )
    yprime(3)=yprime(3)/r  ! modified by shlee 5/30/2001
       ! modified equation from lagrangian 2001/09  , 2002/01/23 -- 1/D corrected to D
    yprime(4)=D*over_B2 *( &
         fld%br*fr + fld%bz*fz + fld%bphi*fp &
         + rho*( fr*(fld%dbzdp/r-fld%dbpdz) + fz*(fld%bphi/r+fld%dbpdr-fld%dbrdp/r) + fp*(fld%dbrdz-fld%dbzdr)) &
         )

    yprime(5)=0D0  ! constant magnetic moment
    yprime(9)=0D0
    yprime(8)=0D0


    ! yprime2 : ExB force only-----------------------------------------------

    fr =  fld%Er ! fr is force over q . -gradient of  Hamiltonian
    fp =  fld%Ephi   ! modified by shlee 5/30/2001 -- /r added
    fz =  fld%Ez 
    yprime2(1)= D* (fld%bz*Fp - fld%Bphi * Fz) * over_B2         
    yprime2(2)= D* (fld%bphi * fr - fld%br * fp ) * over_B2      
    yprime2(3)= D* (fld%br * fz - fld%bz * fr) * over_B2         
    yprime2(3)=yprime(3)/r  ! modified by shlee 5/30/2001
       ! modified equation from lagrangian 2001/09  , 2002/01/23 -- 1/D corrected to D
    yprime2(4)=D*over_B2 *( &
         fld%br*fr + fld%bz*fz + fld%bphi*fp &
         + rho*( fr*(fld%dbzdp/r-fld%dbpdz) + fz*(fld%bphi/r+fld%dbpdr-fld%dbrdp/r) + fp*(fld%dbrdz-fld%dbzdr)) &
         )


end subroutine

subroutine diagnosis_time_avg(istep,ptl)
  use eq_module , only : eq_x_psi, eq_x_z
  use sml_module
  use sml_module, only : sml_mype,sml_totalpe, sml_minusb, sml_time, sml_e_charge
  use ptl_module
  use diag_module
  use perf_monitor
  implicit none
  integer, intent(in) :: istep
  type(ptl_type),target :: ptl
  integer :: i,j,st,jp
  type(species_type) , pointer :: sp
!  real (kind=8) :: diag_avg_flux(diag_avg_npsi,ptl_sp_num), eflux(diag_avg_npsi,ptl_sp_num),&
!       gradpsi(diag_avg_npsi,ptl_sp_num), weight(diag_avg_npsi,ptl_sp_num),dum(diag_avg_npsi,ptl_sp_num)
!  real (kind=8) :: t_para(diag_avg_npsi,ptl_sp_num), t_perp(diag_avg_npsi,ptl_sp_num)
!  real (kind=8) :: dum(diag_avg_npsi,ptl_sp_num)
  real (kind=8) :: psi,br,bz,bphi,b2,dpsi_dr,dpsi_dz, gradpsi2,psidot,pweight,pweight2
  real (kind=8) :: k_para, k_perp,v_para
  real (kind=8), external :: b_interpol
  real (kind=8) :: B,r,z,phi,charge,rho,mu
  real (kind=8) :: norm_tmp
  real (kind=8) :: aa(0:1)
  real (kind=8) :: c2_2m
  real (kind=8) :: c_m 
!  save dum
 
  call monitor_start (DIAG_TIME_AVG_)

  do st=1, 1+sml_electron_on
     if(st==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif

     c2_2m=ptl_c2_2m(st)
     c_m=ptl_c_m(st)
     

     do i=1,sp%num
        if(sp%gid(i) > 0) then
           r=sp%phase(1,i)
           z=sp%phase(2,i)
           phi=sp%phase(3,i)
           rho=sp%phase(4,i)
           mu=sp%phase(5,i)
           psi=sp%phase(9,i)
           !           if(efld_pin < psi .AND. psi < efld_pout .AND.( y(2) > NC_xD_z .or. psi > NC_x_psi )) then ! exclude region 3 only 2003/11/04
           if(diag_avg_pin < psi .AND. psi < diag_avg_pout .AND. z > eq_x_z ) then
              j=int((psi-diag_avg_pin)/diag_avg_dp)+1
              j=min(diag_avg_npsi-1,max(1,j))
              aa(1)= (psi-diag_avg_pin)/diag_avg_dp +1 -j
              aa(0)=1D0-aa(1)
              !              print *,i,efld_m_psi, j !debug
              call bvec_interpol(r,z,phi,br,bz,bphi)
              b2=br**2+bz**2+bphi**2
              b=sqrt(b2)
              dpsi_dr= bz*r ! get grad psi from bvec
              dpsi_dz= -br*r
              if(sml_minusb==1) then 
                 dpsi_dr=-dpsi_dr
                 dpsi_dz=-dpsi_dz
              endif
              gradpsi2=dpsi_dr**2+dpsi_dz**2
!              psidot=sp%derivs(1,i)*dpsi_dr+sp%derivs(2,i)*dpsi_dz ! radial velocity  * grad psi  - FOR XGC0
              psidot=sp%derivs(1,sp%pc_index(1),i)*dpsi_dr+sp%derivs(2,sp%pc_index(1),i)*dpsi_dz ! radial velocity  * grad psi - FOR XGC1
#ifdef COREVV2
              psidot=ptl_exb(1,i)*dpsi_dr+ptl_exb(2,i)*dpsi_dz
#endif
              k_para=c2_2m*rho**2*b2     ! parallel energy
              v_para=c_m*rho*b           ! parallel velocity
              k_perp=mu*b                           ! perp energy
              
              pweight2=sp%phase(8,i)
              if(sml_deltaf==1) then
                 pweight=pweight2*sp%phase(6,i)
              else
                 pweight=pweight2  ! full-f weight
              endif
              
              do jp=0,1
                 diag_avg_flux(j+jp,1,st) = aa(jp)               *pweight + diag_avg_flux(j+jp,1,st)   ! weight itself
                 diag_avg_flux(j+jp,2,st) = aa(jp)*psidot        *pweight + diag_avg_flux(j+jp,2,st)   ! radial particle flux
                 diag_avg_flux(j+jp,3,st) = aa(jp)*psidot*k_para *pweight + diag_avg_flux(j+jp,3,st)   ! parallel Energy flux
                 diag_avg_flux(j+jp,4,st) = aa(jp)*psidot*k_perp *pweight + diag_avg_flux(j+jp,4,st)   ! perp.    Energy flux
                 diag_avg_flux(j+jp,5,st) = aa(jp)*v_para        *pweight + diag_avg_flux(j+jp,5,st)   ! parallel flow              
                 diag_avg_flux(j+jp,6,st) = aa(jp)*k_para        *pweight + diag_avg_flux(j+jp,6,st)   ! Parallel energy
                 diag_avg_flux(j+jp,7,st) = aa(jp)*k_perp        *pweight + diag_avg_flux(j+jp,7,st)   ! Perp energy
                 diag_avg_flux(j+jp,8,st) = aa(jp)*psidot*v_para *pweight + diag_avg_flux(j+jp,8,st)   ! v_r v_|| correlation
                 diag_avg_flux(j+jp,9,st) = aa(jp)*gradpsi2      *pweight2+ diag_avg_flux(j+jp,9,st)   ! grad_psi **2
                 diag_avg_flux(j+jp,10,st)= aa(jp)               *pweight2+ diag_avg_flux(j+jp,10,st)  ! full f weight  

              enddo
           end if
        endif
     enddo
  enddo
  
  call monitor_stop (DIAG_TIME_AVG_)
          
end subroutine diagnosis_time_avg

subroutine time_avg_output(istep) !1359:1548
  use diag_module
  use sml_module
  use eq_module
  use perf_monitor
  use ptl_module
  implicit none
  integer, intent(in) :: istep
  integer, parameter :: NN=diag_avg_nv1, NA=diag_avg_nv2
  real (kind=8) :: flux_total(diag_avg_npsi,NN,1+sml_electron_on)
  integer :: array_size
  real (kind=8) :: composite_value(diag_avg_npsi,NA,1+sml_electron_on)
  real (kind=8) :: psi, weight0(diag_avg_npsi,1+sml_electron_on), en_sum0(2,diag_avg_npsi,1+sml_electron_on)
  real (kind=8) :: tmp(diag_avg_npsi)
  real (kind=8),external :: f0_den,f0_temp_ev
  integer :: i,j,st

#ifdef ADIOS
  real (kind=8) :: simtime(1)
  integer*8 :: grp_id,buf_id
  real (kind=8), allocatable:: psi_eq_x_psi(:)
  character (len=5) :: sp_name(2)=(/'ion__', 'elec_'/)
  character (len=30) :: filename2
  real (kind=8), dimension(diag_flow_npsi) :: psi1
  character (len=32) :: var_names(diag_flow_nvars2)
  ! flux name definition
  var_names(1 )='psi'
  var_names(2 )='radial_flux(avg)'
  var_names(3 )='radial_E_para_flux(avg)'
  var_names(4 )='radial_E_perp_flux(avg)'
  var_names(5 )='parallel_flow(avg)'
  var_names(6 )='parallel_mean_E(avg)'
  var_names(7 )='perp_temperature(avg)'
  var_names(8 )='v_r_v_para_corel(avg)'
  var_names(9 )='average_gradpsi_sqr(avg)'
  var_names(10)='full_weight(avg)'
  var_names(11)='density(avg)'
  var_names(12)='total_E_flux(avg)'
  var_names(13)='grad_T(avg)'
  var_names(14)='Thermal_conductivity(psi:avg)'
  var_names(15)='Thermal_conductivity(m:avg)'

  allocate(psi_eq_x_psi(diag_avg_npsi))
#endif
!  MPI reduce using my_mpi_reduce

  call monitor_start (TIME_AVG_RED_)

  array_size = NN*diag_avg_npsi*(1+sml_electron_on)
  call my_mpi_reduce(diag_avg_flux,flux_total,array_size)
  
  call monitor_stop (TIME_AVG_RED_)
  
  if(sml_mype==0) then
     ! delta-f correction
     do st=1, 1+sml_electron_on

        if(sml_deltaf==1) then
           do i=1,diag_avg_npsi
              psi=diag_avg_dp*(real(i)-0.5) + diag_avg_pin
              weight0(i,st)=diag_avg_vol(i)*f0_den(psi,0D0)
              en_sum0(1,i,st)=0.5D0*weight0(i,st)*f0_temp_ev(psi,0D0)*sml_ev2j
              en_sum0(2,i,st)=2D0*en_sum0(1,i,st)
           enddo
           flux_total(:,1,st)=flux_total(:,1,st) + weight0(:,st)*real(diag_avg_outperiod)
           flux_total(:,6,st)=flux_total(:,6,st) + en_sum0(1,:,st)*real(diag_avg_outperiod)
           flux_total(:,7,st)=flux_total(:,7,st) + en_sum0(2,:,st)*real(diag_avg_outperiod)
        endif

     enddo

     ! get average value
     do st=1, 1+sml_electron_on
        
        flux_total(:,2,st) = flux_total(:,2,st) / flux_total(:,1,st)  ! Flux/n - radial flow
        flux_total(:,3,st) = flux_total(:,3,st) / flux_total(:,1,st)  ! E_para_Flux / n - radial e_para flow
        flux_total(:,4,st) = flux_total(:,4,st) / flux_total(:,1,st)  ! E_perp_flux / n - radial e_perp flow 
        flux_total(:,5,st) = flux_total(:,5,st) / flux_total(:,1,st)  ! Parallel mean flow
        flux_total(:,6,st) = flux_total(:,6,st) / flux_total(:,1,st)  ! Paralllel mean energy
        flux_total(:,7,st) = flux_total(:,7,st) / flux_total(:,1,st)  ! Perp Temp.
        flux_total(:,8,st) = flux_total(:,8,st) / flux_total(:,1,st)  ! coorelation
        flux_total(:,9,st) = flux_total(:,9,st) / flux_total(:,10,st)  ! average gradpsi**2 
        
        ! get composite value - density, thermal conductivity
        ! ** density
        composite_value(:,1,st) = flux_total(:,1,st) / diag_avg_vol(:)  /real(diag_avg_outperiod)  ! density
        ! ** thermal flux
        composite_value(:,2,st) = (  flux_total(:,3,st) + flux_total(:,4,st)  & ! total energy  flow
             - flux_total(:,2,st)*( flux_total(:,6,st) + flux_total(:,7,st) ) & ! - <v_r> *( v^2)
             + ptl_mass(st) * flux_total(:,2,st) * &
             ( flux_total(:,2,st)*flux_total(:,5,st) - flux_total(:,8,st) ) & ! m * u_|| * ( <v_r> u_|| - <v_r v||> )
             )
        ! ** gradient
        tmp(:)=flux_total(:,6,st) + flux_total(:,7,st) - 0.5D0*ptl_mass(st)*flux_total(:,5,st)**2
        composite_value(2:diag_avg_npsi-1,3,st) = (tmp(3:diag_avg_npsi) - tmp(1:diag_avg_npsi-2))/diag_avg_dp*0.5D0   !temp grad.(psi)
        composite_value(1,3,st)=composite_value(2,3,st)
        composite_value(diag_avg_npsi,3,st)=composite_value(diag_avg_npsi-1,3,st)
        ! ** thermal conductivity (psi-space)
        composite_value(:,4,st)=composite_value(:,2,st)/composite_value(:,3,st)
        ! ** thermal conductivity (real-space)
        composite_value(:,5,st)=composite_value(:,4,st)/flux_total(:,9,st)

     enddo

     !write
     
     do st=1, 1+sml_electron_on
        ! normalize
!!$        flux_total(:,2,st) = flux_total(:,2,st)    ! radial flow psi/s * 
!!$        flux_total(:,3:4,st) = flux_total(:,3:4,st)  ! J psi/s
!!$        flux_total(:,5,st) = flux_total(:,5,st)     ! Parallel mean flow m/s       
        flux_total(:,6:7,st) = flux_total(:,6:7,st)*sml_j2ev       ! eV
!!$        flux_total(:,8,st) = flux_total(:,8,st)  ! m^2/s^2
!!$        flux_total(:,9,st) = flux_total(:,9,st)  ! dimension psi^2/m^2
!!$        
!!$
!!$        composite_value(:,1,st)=composite_value(:,1,st) 
!!$        composite_value(:,2,st)=composite_value(:,2,st)  ! J psi/s
!!$        composite_value(:,3,st)=composite_value(:,3,st)  ! J / psi
!!$        composite_value(:,4,st)=composite_value(:,4,st)  ! psi^2/s
!!$        composite_value(:,5,st)=composite_value(:,5,st)  ! m^2/s
        
#ifdef ADIOS       
     if(sml_mype==0)write(0,*)"write fort.50" 
!! NEW OPEN
     !call adios_get_group(grp_id,"fort.50"//char(0))
     !call adios_open_append(buf_id,grp_id,"fort.bp"//char(0))
     call adios_open(buf_id,"fort.50"//char(0),"fort.bp"//char(0),"a"//char(0))

     !if(st==1)call adios_set_path(grp_id,"fort.50/ion"//char(0))
     !if(st==2)call adios_set_path(grp_id,"fort.50/electron"//char(0))
     if(st==1)call adios_set_path(buf_id,"fort.50/ion"//char(0))
     if(st==2)call adios_set_path(buf_id,"fort.50/electron"//char(0))
#else
        write(50+st-1,*) ' '
        write(50+st-1,*) '# t(tau)', (istep-1)*sml_dt/sml_tran 

#endif 
        do i=1, diag_avg_npsi
           psi= diag_avg_pin + diag_avg_dp*real(i-0.5)
#ifdef ADIOS
           psi_eq_x_psi(i)= diag_avg_pin + diag_avg_dp*real(i-0.5)
#else
           write(50+st-1,1000) psi/eq_x_psi, (flux_total(i,j,st),j=2,NN), (composite_value(i,j,st),j=1,NA)           
#endif
        enddo
     enddo
#ifdef ADIOS       
     call adios_gwrite(buf_id,"fort.50") 
     call adios_close(buf_id)
     deallocate(psi_eq_x_psi)
#endif

#ifdef ADIOS
     if(istep==diag_avg_outperiod) then
        do j=1, diag_avg_npsi
           psi1(j) = (diag_avg_pin + diag_avg_dp * (real(j)-0.5D0)) / eq_x_psi
        enddo
     endif
     filename2 ="fluxdiag.bp" //char(0)
!! NEW OPEN
!     call adios_get_group(grp_id,"diagnosis.flux"//char(0))
!     call adios_open_append(buf_id,grp_id,filename2)
     call adios_open(buf_id,"diagnosis.flux"//char(0),filename2,"a"//char(0))
     simtime(1) = sml_time/sml_tran
     if(sml_mype==0)write(0,*)"istep: ",istep/diag_avg_outperiod
     call adios_gwrite(buf_id,"diagnosis.flux")
     !do st=1, 1+sml_electron_on
     !   do j=2,NN
     !      filename2=sp_name(st)//trim(var_names(j))//char(0)
     !      call adios_write(buf_id,filename2,flux_total(:,j,st))
     !   enddo
     !   do j=NN+1, NN+NA
     !      filename2=sp_name(st)//trim(var_names(j))//char(0)
     !      call adios_write(buf_id,filename2,composite_value(:,j,st))
     !   enddo
     !enddo
     call adios_close(buf_id)
#endif
  endif

  ! reset variables
  diag_avg_flux=0D0

1000 format(100(e19.13,1x))
end subroutine time_avg_output


subroutine M3D_coupling_tail
  use lim_module
  use eq_module
  use sml_module
  implicit none
  integer :: i,rl,j
  real (kind=8) :: dr,dz,psi0
  character (LEN=10) :: date


  ! Limiter
  write(444,*) lim_zindex(1)+lim_zindex(2)
  write(444,2002) (((lim_r(i,rl)),(lim_z(i,rl)),i=1,lim_zindex(rl)),rl=1,2)

  ! Mesh point
  dr=(eq_max_r-eq_min_r)/real(eq_mr-1)
  dz=(eq_max_z-eq_min_z)/real(eq_mz-1)
  write(444,*) eq_mr
  write(444,2004) ((eq_min_r+dr*(i-1)), i=1, eq_mr)
  write(444,*) eq_mz
  ! Psi value
  write(444,2004) ((eq_min_z+dz*(j-1)), j=1, eq_mz)
  write(444,2004) ((eq_psi_rz(i,j), i=1, eq_mr),j=1,eq_mz)

2002 format(2(e19.13,' '))
2004 format(2(e19.13,' '))

  !Tail information
  write(444,*) eq_header,eq_filename
  call date_and_time(date)
  write(444,*) "Created by XGC2 ",date,'.'


end subroutine M3D_coupling_tail


subroutine pweight_diagnosis(istep,ptl)
  use ptl_module
  use diag_module, only : diag_flow_npsi, diag_flow_dp, &
     diag_flow_pin, diag_flow_pout, diag_pw_period
  use eq_module, only : eq_x_z
  use sml_module, only : sml_mype
  use perf_monitor
  implicit none
  type(ptl_type),target :: ptl  
  type(species_type) , pointer:: sp
  integer, intent(in) :: istep
  real (kind=8), dimension(diag_flow_npsi) :: pweight, pwdev, wdum, count_in
  real (kind=8) :: psi, z, ww, wdev
  integer :: i, j
  character(len=12) :: filename
#ifdef ADIOS
  integer*8:: grp_id,buf_id
  real (kind=8), dimension(diag_flow_npsi) :: psi_arr,pw_arr
#endif 
  sp=> ptl%ion

  pweight=0D0
  pwdev=0D0
  count_in=0D0
  do i=1, sp%num
     if(sp%gid(i)>0) then
        psi=sp%phase(9,i)
        z=sp%phase(2,i)
        if(diag_flow_pin<psi .AND. psi<diag_flow_pout .AND. z>eq_x_z) then
           ww = (1D0 + sp%phase(6,i)) * sp%phase(8,i)
           wdev = ww - sp%weight0(i)
           j=int((psi-diag_flow_pin)/diag_flow_dp)+1
           j=min(diag_flow_npsi,max(1,j))
           count_in(j) = count_in(j) + 1D0
           pweight(j) = pweight(j) + ww
           pwdev(j) = pwdev(j) + wdev*wdev
        endif
     endif
  enddo

  call monitor_start (PWEIGHT_RED_)
! MPI reduce using my_mpi_reduce
  call my_mpi_reduce(count_in,wdum,diag_flow_npsi)
  count_in=wdum
  call my_mpi_reduce(pweight,wdum,diag_flow_npsi)
  pweight=wdum
  call my_mpi_reduce(pwdev,wdum,diag_flow_npsi)
  pwdev=wdum
  call monitor_stop (PWEIGHT_RED_)

  if (sml_mype==0) then
     write(filename,'("pw",i5.5,".dat")') istep/diag_pw_period
     open(555, file=filename, status='replace', form='formatted')
     do j=1, diag_flow_npsi
        psi = diag_flow_pin + diag_flow_dp * (real(j)-0.5D0)
        write(555,2005) psi,sqrt(pwdev(j)/count_in(j))/(pweight(j)/count_in(j))
#ifdef ADIOS
        psi_arr(i)=psi
        pw_arr(i)=sqrt(pwdev(j)/count_in(j))/(pweight(j)/count_in(j))     
#endif 
     enddo
     close(555)
#ifdef ADIOS
!!NEW OPEN
     !call adios_get_group(grp_id,"diagnosis.particleweight"//char(0))
     !call adios_open_append(buf_id,grp_id,"diagnosis_pw.bp"//char(0))
     call adios_open(buf_id,"diagnosis.particleweight"//char(0),"diagnosis_pw.bp"//char(0),"a"//char(0))
     !call adios_set_path(grp_id,"particleweight"//char(0))
     call adios_set_path(buf_id,"particleweight"//char(0))
     call adios_gwrite(buf_id,"diagnosis.particleweight")
     call adios_close(buf_id)
#endif
  endif

2005 format(2(e19.13,' '))  
end subroutine pweight_diagnosis
subroutine particle_dump(istep,ptl)
  use ptl_module 
  use diag_module, only : diag_ptl_begin, diag_ptl_num, &
    diag_ptl_ifile, diag_ptl_efile, diag_ptl_times, &
    diag_ptl_data1, diag_ptl_data2
  use sml_module, only : sml_mype, sml_time,  &
    sml_deltaf, sml_electron_on  
  implicit none
  type(ptl_type) , target:: ptl
  type(species_type), pointer :: sp
  integer, intent(in) :: istep
  integer :: i, sp_type, pnum, index
  real (kind=8) :: c_m, inv_m, r, z, phi, B
  real (kind=8), external :: b_interpol
  character (len=50) :: grpname 
#ifdef ADIOS
  integer*8 :: buf
  character (len=3) :: sp_name(2)=(/'ion', 'ele'/)
#endif

#ifdef ADIOS
  do sp_type=1, 1+sml_electron_on

     call adios_open(buf,"diagnosis.particle"//char(0),"particle.bp"//char(0),"a"//char(0))
     call adios_set_path(buf,sp_name)
     c_m=ptl_c_m(sp_type)
     inv_m=1.D0/ptl_mass(sp_type)
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif    
     if (sml_mype==0) then
        pnum=diag_ptl_num
! save current time for this particle dump in array
        index = istep - diag_ptl_begin
        diag_ptl_times(index+1) = sml_time 
! Write current timestamp
        call adios_write(buf,"sml_time"//char(0))
        write(grpname,'("time_node_data[",i5.5,"]")') index
        call adios_set_path(buf,grpname//char(0)) 
! write global ids
        call adios_write(buf,"gid"//char(0),"sp%gid(1:pnum)"//char(0))
! write position data
        diag_ptl_data1(1:pnum)=sp%phase(1,1:pnum)
        diag_ptl_data2(1:pnum)=sp%phase(2,1:pnum)
        call adios_write(buf,"position_r"//char(0),"diag_ptl_data1"//char(0))
        call adios_write(buf,"position_z"//char(0),"diag_ptl_data2"//char(0))
        
! compute particle velocity components
        do i=1,pnum
           if (sp%gid(i)>0) then
              r = sp%phase(1,i)
              z = sp%phase(2,i)
              phi=sp%phase(3,i)
              B=b_interpol(r,z,phi)
              diag_ptl_data1(i) = c_m*sp%phase(4,i)*B 
              diag_ptl_data2(i) = sqrt(2.D0*inv_m*sp%phase(5,i)*B)
           else
              diag_ptl_data1(i) = 0.D0
              diag_ptl_data2(i) = 0.D0
           endif
        enddo
! write particle velocity components
        call adios_write(buf,"vpar"//char(0),"diag_ptl_data1"//char(0))
        call adios_write(buf,"vperp"//char(0),"diag_ptl_data2"//char(0))
! write particle weight data
        if(sml_deltaf==1) then
           diag_ptl_data1(1:pnum)=sp%phase(8,1:pnum)*sp%phase(6,1:pnum)
        else
           diag_ptl_data1(1:pnum)=sp%phase(8,1:pnum) ! for full f simulation only
        endif
        call adios_write(buf,"weight"//char(0),"diag_ptl_data1"//char(0))
     endif
     call adios_close(buf)
   enddo
#endif
end subroutine particle_dump

subroutine particle_dump_finalize
  use diag_module, only : diag_ptl_on, diag_ptl_begin, diag_ptl_end, &
    diag_ptl_ifile, diag_ptl_efile, diag_ptl_times, &
    diag_ptl_data1, diag_ptl_data2
  use sml_module, only : sml_electron_on
  implicit none
  integer :: sp, nsteps
  integer (kind=8) :: fileptr
  integer, parameter :: vardim=1
  character (len=30) :: varname
#if BINPACK
  if (diag_ptl_on==1) then
    nsteps = diag_ptl_end - diag_ptl_begin + 1
    varname = "time"
    do sp=1, 1+sml_electron_on
      if(sp==1) then
        fileptr=diag_ptl_ifile
      else
        fileptr=diag_ptl_efile
      endif
! end of timestep series
      call bpendtimestepgroup(fileptr)
! write out timestamp array
      call bpwritevardouble(fileptr,diag_ptl_times,vardim,nsteps,varname)
! close binpack file
      call bpclosefile(fileptr)
    enddo

! deallocate particle diagnostic arrays
    deallocate(diag_ptl_data1,diag_ptl_data2,diag_ptl_times)
  endif
#endif
end subroutine particle_dump_finalize
!! ADIOS JC
!!#endif

subroutine diag_2d(istep,grid,psn,ptl)
  use diag_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  use sml_module
  use smooth_module
  use perf_monitor
  implicit none
  integer :: istep
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  type(ptl_type),target :: ptl
  type(species_type), pointer :: sp
  integer :: sp_type,i,nd1,nd2,j
  real (kind=8) , allocatable :: weight(:),flow(:),flowe(:),enpa(:),enpe(:),&
     flower(:),flowez(:),dum(:)
  real (kind=8) :: c_m,c2_m,B,Bpol,vp,en1,en2,particle_weight,p(3),ppw,psi,&
     D,nb_curl_nb,over_B2,ver,vez,ve
  real (kind=8),external :: f0_den,f0_temp_ev
  integer :: nodes(3),node,itr,n_n,l,nlarmor
  character (len=30) :: filename, filename2

#ifdef ADIOS 
  integer :: nspace
  integer*8        :: buf_id 
#endif


  character (len=5) :: sp_name(2)=(/'ion__', 'elec_'/)
  n_n = grid%nnode
  allocate(weight(n_n),flow(n_n),flowe(n_n),&
       enpa(n_n),enpe(n_n),flower(n_n),flowez(n_n),dum(n_n))

  do sp_type=1, 1+sml_electron_on

     if(sp_type==1) then
        sp=>ptl%ion
        nlarmor=sml_nlarmor
     else
        sp=>ptl%elec
        nlarmor=1
     endif

     weight=1D-100
     flow=0D0
     enpa=0D0
     enpe=0D0
     flowe=0D0
     flower=0D0
     flowez=0D0

     c_m=ptl_c_m(sp_type) 
     c2_m=ptl_charge(sp_type)*ptl_c_m(sp_type) 

     do i=1, sp%num
        if(sp%gid(i)>0) then ! Assume B-field variation within gyro radius is neglegible.
           fld%r = sp%phase(1,i)
           fld%z = sp%phase(2,i)
           fld%phi=sp%phase(3,i)
           fld%psi=sp%phase(9,i)
           call field(fld,sml_time)
           B = sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
           vp = c_m*sp%phase(4,i)*B
           en1 = c2_m* (sp%phase(4,i)*B)**2
           en2 = sp%phase(5,i)*B

           ! poloidal ExB drift
           call efield(grid,psn,sp,i,fld,sml_time)  ! Get gyro averaged E-field
           ! normalized b dot curl of normalized b
           over_B2 = 1.D0 / B**2
           nb_curl_nb = over_B2 * ( fld%dbpdz * fld%br + &
              (fld%dbzdr - fld%dbrdz) * fld%bphi - &
              (fld%bphi / fld%r + fld%dbpdr) * fld%bz )
           D = 1.D0 / ( 1.D0 + sp%phase(4,i) * nb_curl_nb )
           ver = D * (fld%bz * fld%Ephi - fld%Bphi * fld%Ez) * over_B2
           vez = D * (fld%bphi * fld%Er - fld%br * fld%Ephi) * over_B2
           Bpol = sqrt(fld%br**2 + fld%bz**2) + 1.D-30
           ve = (ver*fld%br + vez*fld%bz) / Bpol

           if(sml_deltaf==1) then
              particle_weight=sp%phase(8,i)*sp%phase(6,i)
           else
              particle_weight=sp%phase(8,i) ! for full f simulation only
           endif

           ! find gyrocenter position in triangle mesh
           !call search_tr(grid,sp%phase(1:2,i),itr,p)
           
           ! Assign values to gyro-ring of ion / center of electron           
           do l=1, nlarmor
              if( sp%tr_save(l,1,i)>0 ) then           ! incomplete  - iphi_frac index summation is needed for exact calculation
                 nodes=grid%nd(:,sp%tr_save(l,1,i))
                 p=sp%p_save(:,l,1,i)
                 do j=1,3
                    node=nodes(j)
                    ppw=p(j)*particle_weight/real(nlarmor)
                    weight(node)= weight(node)+ppw
                    flow(node)=flow(node)+ppw*vp
                    flowe(node)=flowe(node)+ppw*ve
                    enpa(node)=enpa(node)+ppw*en1
                    enpe(node)=enpe(node)+ppw*en2
                    flower(node)=flower(node)+ppw*ver
                    flowez(node)=flowez(node)+ppw*vez
                 enddo
              endif
           enddo
        endif
     enddo
     
     call monitor_start (DIAG2D_RED_)
     !reduce    
     call my_mpi_reduce(weight,dum,n_n)
     weight=dum
     call my_mpi_reduce(flow,dum,n_n)
     flow=dum
     call my_mpi_reduce(flowe,dum,n_n)
     flowe=dum
     call my_mpi_reduce(enpa,dum,n_n)
     enpa=dum
     call my_mpi_reduce(enpe,dum,n_n)
     enpe=dum
     call my_mpi_reduce(flower,dum,n_n)
     flower=dum
     call my_mpi_reduce(flowez,dum,n_n)
     flowez=dum
     call monitor_stop (DIAG2D_RED_)
     
     if(sml_mype==0) then

        if(sml_deltaf==1) then           
           do i=1, grid%npsi
              nd1=grid%itheta0(i)
              nd2=grid%itheta0(i)+grid%ntheta(i)-1              
              psi=grid%psi(nd1)
!              dum(grid%inv_sort(nd1:nd2))=grid%node_vol(grid%inv_sort(nd1:nd2))*sml_nphi_total*f0_den(psi,0D0)
!              enpa(grid%inv_sort(nd1:nd2))=enpa(grid%inv_sort(nd1:nd2))+dum(grid%inv_sort(nd1:nd2))*f0_temp_ev(psi,0D0)*sml_ev2j
!              enpe(grid%inv_sort(nd1:nd2))=enpe(grid%inv_sort(nd1:nd2))+dum(grid%inv_sort(nd1:nd2))*f0_temp_ev(psi,0D0)*sml_ev2j
               dum(nd1:nd2)=grid%node_vol(nd1:nd2)*sml_nphi_total*f0_den(psi,0D0)
               enpa(nd1:nd2)=enpa(nd1:nd2)+dum(nd1:nd2)*f0_temp_ev(psi,0D0)*sml_ev2j
               enpe(nd1:nd2)=enpe(nd1:nd2)+dum(nd1:nd2)*f0_temp_ev(psi,0D0)*sml_ev2j
           enddo
           weight(:)=dum+weight(:)           
        endif

        flow=flow/weight
        flowe=flowe/weight
        enpa=enpa/weight
        enpe=enpe/weight
        flower=flower/weight
        flowez=flowez/weight
!!$        
!!$        do i=grid%inner_bd(1), grid%inner_bd(2)
!!$           flow(i)=0D0
!!$           flowe(i)=0D0
!!$           enpa(i)=0D0
!!$           enpe(i)=0D0
!!$           flower(i)=0D0
!!$           flowez(i)=0D0
!!$        enddo
!!$        do i=grid%outer_bd(1), grid%outer_bd(2)
!!$           flow(i)=0D0
!!$           flowe(i)=0D0
!!$           enpa(i)=0D0
!!$           enpe(i)=0D0
!!$           flower(i)=0D0
!!$           flowez(i)=0D0
!!$        enddo
        call set_boundary_values(flow,0D0,psn%cbdH)
        call set_boundary_values(flowe,0D0,psn%cbdH)
       
        call smooth_pol0(grid,flow,smoothdiag)
        call smooth_pol0(grid,flowe,smoothdiag)
        call smooth_pol0(grid,enpa,smoothdiag)
        call smooth_pol0(grid,enpe,smoothdiag)
        call smooth_pol0(grid,flower,smoothdiag)
        call smooth_pol0(grid,flowez,smoothdiag)

! normalize diagnostic quantities and convert into desired physical units
        weight(:) = weight(:) * grid%inv_node_vol(:) / &
                    ( real(sml_nphi_total)) ! m^-3
        enpa(:) = enpa(:) * sml_j2ev * 1.D-3           ! keV
        enpe(:) = enpe(:) * sml_j2ev * 1.D-3

#ifdef ASCII
        if (mod(istep,diag_pot_period)==0) then
        ! dump 2d fields in ASCII tabular format for Matlab plotting
           if(sp_type==1) then
              write(filename,'("fort.i",".",i2.2)') &
                (istep)/(diag_pot_period)
           else
              write(filename,'("fort.e",".",i2.2)') &
                (istep)/(diag_pot_period)
           endif

           open(501, file=filename, status='replace', form='formatted') 
           do i=1, n_n
              write(501,3001) grid%x(1,i), &
                              grid%x(2,i), &
                              weight(i), flow(i), flowe(i), &
                              enpa(i), enpe(i)
           enddo
           close(501)
        endif
#endif

#ifdef ADIOS 
        if (mod(istep,diag_binout_period)==0) then
           ! dump 2d field diagnostics in binpack format
           if(sml_mype==0)write(0,*)"my diagnosis.F90"
           if(sp_type==1) then
              write(filename2,'("fieldi",".",i4.4,".bp")') (istep)/(diag_binout_period)
              filename2=trim(filename2)//char(0)
              if(sml_mype==0)write(0,*)"my diagnosis.F90:",filename2
              call adios_open(buf_id,"diagnosis.fieldi"//char(0),filename2,"w"//char(0))
              call adios_gwrite(buf_id,"diagnosis.fieldi")
              call adios_close(buf_id)
              
              call adios_open(buf_id,"diagnosis.fieldi.1"//char(0),filename2,"a"//char(0))
              call adios_gwrite(buf_id,"diagnosis.fieldi.1")
              call adios_close(buf_id)
              
              call adios_open(buf_id,"diagnosis.fieldi.2"//char(0),filename2,"a"//char(0))
              call adios_gwrite(buf_id,"diagnosis.fieldi.2")
              call adios_close(buf_id)

              call adios_open(buf_id,"diagnosis.fieldi.3"//char(0),filename2,"a"//char(0))
              call adios_gwrite(buf_id,"diagnosis.fieldi.3")
              call adios_close(buf_id)

              call adios_open(buf_id,"diagnosis.fieldi.4"//char(0),filename2,"a"//char(0))
              call adios_gwrite(buf_id,"diagnosis.fieldi.4")
              call adios_close(buf_id)

              call adios_open(buf_id,"diagnosis.fieldi.5"//char(0),filename2,"a"//char(0))
              call adios_gwrite(buf_id,"diagnosis.fieldi.5")
              call adios_close(buf_id)
             if(sml_mype==0)write(*,*)"end my diagnosis.F90"
           endif ! sp_type==1
        endif  !istep%binout_period
#endif
     endif  !sml_mype==0
  enddo  !sp_type=1,1+sml_electron_on

  deallocate(weight,flow,flowe,enpa,enpe,flower,flowez,dum)
3001 format(7(e19.13,' '))
end subroutine diag_2d

!!ADIOS JC
!#ifdef BINPACK
#ifndef ADIOS 
subroutine restart_write(ptl)
  use sml_module
  use ptl_module
  implicit none
  type(ptl_type),target :: ptl
  type(species_type),pointer :: sp
  integer :: i, sp_type
  character (len=50) :: filename, varname, varname2, grpname
  integer (kind=8) :: fileptr
  integer, parameter :: vardim=1

  write(filename,'("xgc.restart.",i4.4,".bp")') sml_mype
  call bpopenfile(filename,fileptr)

! timestep index
  varname = "timestep"
  call bpwritescalarint(fileptr,sml_istep,varname)

! current time
  varname = "timeval"
  call bpwritescalardouble(fileptr,sml_time,varname)

  do sp_type=1, 1+sml_electron_on
     if (sp_type==1) then
        grpname = "restart_data_i"
        varname = "ptl_num_i"
        varname2 = "ptl_maxgid_i"
        sp=>ptl%ion
     else
        grpname = "restart_data_e"
        varname = "ptl_num_e"
        varname2 = "ptl_maxgid_e"
        sp=>ptl%elec
     endif

     call bpwritescalarint(fileptr,sp%num,varname)
#ifdef SMALL_GID
     call bpwritescalarint(fileptr,sp%maxgid,varname2)
#else
     call bpwritescalarlong(fileptr,sp%maxgid,varname2)
#endif
     call bpbeginbasicgroup(fileptr,grpname)
! write global ids
     varname = "gid"
#ifdef SMALL_GID
     call bpwritevarint(fileptr,sp%gid(1:sp%num),vardim,sp%num,varname)
#else
     call bpwritevarlong(fileptr,sp%gid(1:sp%num),vardim,sp%num,varname)
#endif
! write phase space data and particle weights
     do i=1,sp%nphase
        write(varname,'("phase[",i1.1,"]")') i
        call bpwritevardouble(fileptr,sp%phase(i,1:sp%num),vardim,sp%num,varname)
     enddo
     call bpendbasicgroup(fileptr)
  enddo

  call bpclosefile(fileptr)
  return
end subroutine restart_write
#endif

#ifndef ADIOS
subroutine restart_read(ptl)
  use sml_module
  use ptl_module
  implicit none
  type(ptl_type),target :: ptl
  type(species_type),pointer :: sp
  integer :: i, sp_type, varsize, varnum
  character (len=30) :: filename, tmpname, varname, varname2, grpname
  integer (kind=8) :: fileptr
  integer, parameter :: vardim=1

  write(filename,'("xgc.restart.",i4.4,".bp")') sml_mype
  call bprdopenfile(filename,fileptr)

! timestep index
  varname = "timestep"
  call bpreadscalarattr(fileptr,sml_istep,tmpname)
! sanity check: tmpname should be the same as varname
     if (tmpname(1:8) /= varname(1:8)) then
        print*,'Restart:  Expected variable name ',varname
        print*,'          Found variable name ',tmpname
        stop
     endif

! current time
  varname = "timeval"
  call bpreadscalarattr(fileptr,sml_time,tmpname)
! sanity check: tmpname should be the same as varname
     if (tmpname(1:7) /= varname(1:7)) then
        print*,'Restart:  Expected variable name ',varname
        print*,'          Found variable name ',tmpname
        stop
     endif

  do sp_type=1, 1+sml_electron_on
     if (sp_type==1) then
        grpname = "restart_data_i"
        varname = "ptl_num_i"
        varname2 = "ptl_maxgid_i"
        sp=>ptl%ion
     else
        grpname = "restart_data_e"
        varname = "ptl_num_e"
        varname2 = "ptl_maxgid_e"
        sp=>ptl%elec
     endif

     call bpreadscalarattr(fileptr,sp%num,tmpname)
! sanity check: tmpname should be the same as varname
     if (tmpname(1:9) /= varname(1:9)) then
        print*,'Restart:  Expected variable name ',varname
        print*,'          Found variable name ',tmpname
        stop
     endif
! maximum number of particle per processor
     if (sp%maxnum < sp%num) then
        print *,'Not enough memory allocated. Increase ptl_maxnum in module.f90.'
        stop
     endif
     call bpreadscalarattr(fileptr,sp%maxgid,tmpname)
! sanity check: tmpname should be the same as varname2
     if (tmpname(1:12) /= varname2(1:12)) then
        print*,'Restart:  Expected variable name ',varname2
        print*,'          Found variable name ',tmpname
        stop
     endif

     call bprdbeginbasicgroup(fileptr,tmpname)
! sanity check: tmpname should be the same as grpname
     if (tmpname(1:14) /= grpname(1:14)) then
        print*,'Restart:  Expected group name ',grpname
        print*,'          Found group name ',tmpname
        stop
     endif
! read global ids
     varname = "gid"
     call bpreadvarinfo(fileptr,varsize,varnum,tmpname)
#ifdef SMALL_GID
! sanity check: varsize should be size of integer
     if (varsize /= 4) then
        print*,'Restart:  Expected variable size of 4'
        print*,'          Found variable size of ',varsize
        stop
     endif
#else
! sanity check: varsize should be size of integer*8
     if (varsize /= 8) then
        print*,'Restart:  Expected variable size of 8'
        print*,'          Found variable size of ',varsize
        stop
     endif
#endif
! sanity check: varnum should be equal to sp%num
     if (varnum /= sp%num) then
        print*,'Restart:  Expected variable number of ',sp%num
        print*,'          Found variable number of ',varnum
        stop
     endif
! sanity check: tmpname should be the same as varname
     if (tmpname(1:3) /= varname(1:3)) then
        print*,'Restart:  Expected variable name ',varname
        print*,'          Found variable name ',tmpname
        stop
     endif
     call bpreadvardata(fileptr,varsize,varnum,sp%gid(1:sp%num))
! read phase space data and particle weights
     do i=1,sp%nphase
        write(varname,'("phase[",i1.1,"]")') i
        call bpreadvarinfo(fileptr,varsize,varnum,tmpname)
! sanity check: varsize should be size of double
        if (varsize /= 8) then
           print*,'Restart:  Expected variable size of 8'
           print*,'          Found variable size of ',varsize
        stop
        endif
! sanity check: varnum should be equal to sp%num
        if (varnum /= sp%num) then
           print*,'Restart:  Expected variable number of ',sp%num
           print*,'          Found variable number of ',varnum
        stop
        endif
! sanity check: tmpname should be the same as varname
        if (tmpname(1:8) /= varname(1:8)) then
           print*,'Restart:  Expected variable name ',varname
           print*,'          Found variable name ',tmpname
        stop
        endif
        call bpreadvardata(fileptr,varsize,varnum,sp%phase(i,1:sp%num))
     enddo
     call bprdendbasicgroup(fileptr)
  enddo

  call bprdclosefile(fileptr)
  return
end subroutine restart_read
#endif

subroutine dump_grid(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: m, n, dir, itr, n_n, n_t, n_psi, dims(1)
  integer, parameter :: ndim=1
  real (kind=8), pointer :: coord(:,:)
  integer, pointer :: nodeid(:,:), n_itheta0(:), nextn(:)
  character (len=30) :: meshbp_file, varname
  integer (kind=8) :: meshbp_fileptr
#ifdef ADIOS
  integer (kind=8) :: grp, buf 
#endif
  n_n = grid%nnode
  n_t = grid%ntriangle
  allocate(coord(2,n_n),nodeid(3,n_t))
! node coords scaled to physical length unit of meters
  coord(:,:) = grid%x(:,:) 
! node id's with zero-base indexing for C code
  nodeid(:,:) = grid%nd(:,:) - 1

  n_psi = grid%npsi
  allocate(n_itheta0(n_psi))
! node id's with zero-base indexing for C code
  n_itheta0(:) = grid%itheta0(:) - 1

  allocate(nextn(n_n))
  if(sml_turb_efield==1) then
! use bfollow data arrays to determine node connectivity from one plane to next
    dir = 1
    if (sml_bt_sign < 0) dir=2
    do n=1,n_n
      itr = psn%bfollow_tr(dir,n)
      m = 1
      if (psn%bfollow_p(2,dir,n) > psn%bfollow_p(m,dir,n)) m=2
      if (psn%bfollow_p(3,dir,n) > psn%bfollow_p(m,dir,n)) m=3
! node id's with zero-base indexing for C code
      nextn(n) = grid%nd(m,itr) - 1
    enddo
  else
! just assume toroidal symmetry in mesh
    do n=1,n_n
! node id's with zero-base indexing for C code
      nextn(n) = n - 1
    enddo
  endif
#ifdef ADIOS
  write(*,*)"write mesh"
  meshbp_file = "mesh.bp"//char(0)
  call adios_open(buf,"diagnosis.mesh"//char(0),meshbp_file,"w"//char(0))
  call adios_gwrite(buf,"diagnosis.mesh")
  call adios_close(buf)
  write(*,*)"end write mesh"
#else
! Dump mesh information with binpack
  meshbp_file = "xgc.mesh.bp"
  call bpopenfile(meshbp_file,meshbp_fileptr)
  call bpwritetrimesh(meshbp_fileptr,n_n,n_t,coord,nodeid)
! Dump additional info on mesh point ordering and connectivity between planes
  dims(1) = n_psi
  varname = "itheta0"
  call bpwritevarint(meshbp_fileptr,n_itheta0,ndim,dims,varname)
  dims(1) = n_n
  varname = "nextnode"
  call bpwritevarint(meshbp_fileptr,nextn,ndim,dims,varname)
  call bpclosefile(meshbp_fileptr)
  deallocate(coord,nodeid,n_itheta0,nextn)
#endif
end subroutine dump_grid

subroutine dump_bfield(grid)
  use grid_class
  use eq_module
  use sml_module
  implicit none
  type(grid_type) :: grid
  integer :: i,n_n,dims(2)
  real (kind=8) :: r,z,phi,br,bz,bphi
  real (kind=8), allocatable :: bfld(:,:)
  integer, parameter :: nvar=2, veclen=3, flddim=2, nspace=2, veclen2=1,flddim2=1
  character (len=30) :: filename, fldname, units, varname
  integer (kind=8) :: fileptr
#ifdef ADIOS
  integer(kind=8) :: grp,buf
#endif
  n_n = grid%nnode
  allocate(bfld(3,n_n))
  phi = 3.141592/2.
  do i=1,n_n
     r = grid%x(1,i)
     z = grid%x(2,i)
     call bvec_interpol(r,z,phi,br,bz,bphi)
     bfld(1,i) = br 
     bfld(2,i) = bz 
     bfld(3,i) = bphi 
  enddo

#ifndef ADIOS
! Dump bfield components with binpack
  filename = "xgc.bfield.bp"
  call bpopenfile(filename,fileptr)
! dummy mesh data
  varname = "nspace"
  call bpwritescalarint(fileptr,nspace,varname)
  varname = "coordinates"
  call bpbeginbasicgroup(fileptr,varname)
  call bpendbasicgroup(fileptr)
! begin field data
  call bpbeginnodegroup(fileptr,n_n,nvar)
  varname(:) = " "
  fldname = "bfield"
  units = "tesla"
  dims(1) = n_n
  dims(2) = 3
  call bpwritefielddouble(fileptr,veclen,bfld,flddim,dims,&
                          fldname,units,varname)
  varname(:) = " "
  fldname = "poloidalflux"
  units = "Tm^2"
  dims(1) = n_n
  dims(2) = 1
  call bpwritefielddouble(fileptr,veclen2,grid%psi,flddim2,dims,&
       fldname,units,varname)

  call bpendnodegroup(fileptr)
  call bpclosefile(fileptr)
#endif
#ifdef ADIOS
  filename="bfield.bp"//char(0)
!! NEW OPEN
!  call adios_get_group(grp,"diagnosis.bfield"//char(0))
!  call adios_open(buf,grp,filename)
  call adios_open(buf,"diagnosis.bfield"//char(0),filename,"w"//char(0))
  call adios_gwrite(buf,"diagnosis.bfield")
  call adios_close(buf)

!  call adios_get_group(grp,"diagnosis.bfield.1"//char(0))
!  call adios_open_append(buf,grp,filename)
  call adios_open(buf,"diagnosis.bfield.1"//char(0),filename,"a"//char(0))
  call adios_gwrite(buf,"diagnosis.bfield.1")
  call adios_close(buf)

#endif

  deallocate(bfld)
end subroutine dump_bfield
!#endif
!!ADIOS JC

#ifdef ADIOS
!! restart_write_adios
subroutine restart_write(ptl)
  use sml_module
  use ptl_module
  implicit none
  type(ptl_type),target :: ptl
  type(species_type),pointer :: sp
  integer :: i,j
  character (len=30) :: filename, dirname
  integer*8 :: buf_id
  integer*8 :: icomm
  integer   :: imype
  include 'mpif.h'
  
  if(sml_mype==0)then
    open(1,FILE="timestep.dat",FORM="unformatted")
    write(1)sml_istep
    close(1) 
  endif 
  icomm=MPI_COMM_WORLD
  imype=sml_mype 

  write(filename,'("restart.",i4.4,".bp")') sml_istep
  filename=trim(filename)//char(0)
  call adios_open(buf_id,"restart"//char(0),filename,"w"//char(0))
  write(dirname,'("node_",i4.4)') sml_mype
  dirname=trim(dirname)//char(0)
  call adios_set_path(buf_id,dirname)
  sp=>ptl%ion
  call adios_gwrite(buf_id,"restart")
  call adios_close(buf_id)
!endif
end subroutine restart_write
#endif

#ifdef ADIOS
!restart_read_adios
subroutine restart_read(ptl)
  use sml_module
  use ptl_module
  implicit none
  interface
  subroutine adios_read_int(a,b,c)
      integer*8 :: a 
      character::b*(*)
      integer,dimension(:),pointer::c 
  end subroutine adios_read_int 
!  end interface  

!  interface
  subroutine adios_read_long(a,b,c)
      integer*8 :: a 
      character::b*(*)
      integer*8,dimension(:),pointer::c 
  end subroutine adios_read_long
!  end interface  

!  interface
  subroutine adios_read_double(a,b,c)
      integer*8 :: a 
      character::b*(*)
      real*8,dimension(:,:),pointer::c 
  end subroutine adios_read_double 
  end interface  

  type(ptl_type),target :: ptl
  type(species_type),pointer :: sp
  integer :: ierr
  character (len=30) :: filename,varname
  integer, pointer :: pnum
#ifdef SMALL_GID
  integer, pointer :: gid(:)
#else
  integer*8, pointer :: gid(:)
#endif
  real*8 ,pointer:: temp(:)
  real*8 :: arr(9,2)
  integer*8 :: grp_id, buf_id,icomm=0
  integer   :: imype=0,imaxnum=0
  include "mpif.h" 

  if(sml_mype==0)then
    open(1,FILE="timestep.dat",FORM="unformatted")
    read(1)sml_istep
    close(1)
  endif
  call MPI_BCAST(sml_istep, 1,MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  write(filename,'("restart.",i4.4,".bp")') sml_istep
  filename=trim(filename)//char(0)
  !if(sml_mype==10)write(*,*)"restart read file: ",filename 
  call adios_open(buf_id,"restart"//char(0),filename,"r"//char(0))
  write(filename,'("node_",i4.4)') sml_mype
  filename=trim(filename)//char(0)
  call adios_set_path(buf_id,filename)
  icomm=MPI_COMM_WORLD
  sp=>ptl%ion
  call adios_gread(buf_id,"restart")
  !#ifdef SMALL_GID
  !call adios_read_int(buf_id,"igid"//char(0),sp%gid)
  !#else
  !call adios_read_long(buf_id,'igid'//char(0),sp%gid)
  !#endif
  !call adios_read_double(buf_id,'iphase'//char(0),sp%phase)
  call adios_close(buf_id)
  !if(sml_mype==0)write(*,*)"0, restart read end: sp%phase",sp%phase(:,2)
  !if(sml_mype==11)write(*,*)"1, restart read end: sp%phase",sp%phase(:,2)
  !if(sml_mype==22)write(*,*)"2, restart read end: sp%phase",sp%phase(:,2)
  !if(sml_mype==31)write(*,*)"4, restart read end: sp%phase",sp%phase(:,2)
end subroutine restart_read
#endif
 
#ifdef ADIOS
subroutine adios_read_int(a,b,c)
   integer*8 :: a 
   character::b*(*)
   integer,dimension(:),pointer::c 
   call adios_read(a,b,c)
end subroutine adios_read_int

subroutine adios_read_long(a,b,c)
   integer*8 :: a 
   character::b*(*)
   integer*8,dimension(:),pointer::c 
   call adios_read(a,b,c)
end subroutine adios_read_long

subroutine adios_read_double(a,b,c)
   integer*8 :: a 
   character::b*(*)
   !real*8,dimension(:),pointer::c
   real*8,dimension(:,:) ,pointer::c
   call adios_read(a,b,c)
end subroutine adios_read_double
#endif

!!$subroutine diagnosis_tot_cam(fnum)
!!$  use sml_module
!!$  use ptl_module
!!$  implicit none
!!$  include 'mpif.h'
!!$  integer :: fnum
!!$  integer :: i,sp,ierror
!!$  integer, pointer :: pnum,gid(:)
!!$  real (kind=8), pointer :: phase(:,:)
!!$  real (kind=8) :: tot_cam(2),dum(2),psi,psi_c,r,z,phi,rho,mu,br,bz,bphi
!!$
!!$  tot_cam=0D0
!!$
!!$  do sp=1, 1+sml_electron_on
!!$     if(sp==1) then
!!$        phase=>ptl_phase
!!$        pnum=>ptl_num
!!$        gid=>ptl_gid
!!$     else
!!$        phase=>ptl_electron
!!$        pnum=>ptl_e_num
!!$        gid=>ptl_e_gid
!!$     endif
!!$     
!!$     
!!$     do i=1, pnum
!!$        if(gid(i)>0) then
!!$           r=phase(i,1)
!!$           z=phase(i,2)
!!$           phi=phase(i,3)
!!$           rho=phase(i,4)
!!$           mu=phase(i,5)
!!$           psi=phase(i,9)
!!$
!!$           call bvec_interpol(r,z,phi,br,bz,bphi)
!!$                      
!!$           if(sml_minusB/=1) then
!!$              psi_c=psi + rho*bphi*r
!!$           else
!!$              psi_c=psi - rho*bphi*r
!!$           endif
!!$           
!!$           tot_cam(sp)=tot_cam(sp)+psi_c
!!$
!!$
!!$
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$  call my_mpi_reduce(tot_cam,dum,2)
!!$
!!$  if(sml_mype==0) then
!!$     write(fnum,*) dum
!!$  end if
!!$
!!$end subroutine diagnosis_tot_cam
!!$
!!$
!!$subroutine diagnosis_tot_mu(fnum)
!!$  use sml_module
!!$  use ptl_module
!!$  implicit none
!!$  include 'mpif.h'
!!$  integer :: fnum
!!$  integer :: i,sp,ierror
!!$  integer, pointer :: pnum,gid(:)
!!$  real (kind=8), pointer :: phase(:,:)
!!$  real (kind=8) :: tot_mu(2),dum(2),psi,psi_c,r,z,phi,rho,mu,br,bz,bphi
!!$
!!$  tot_mu=0D0
!!$
!!$  do sp=1, 1+sml_electron_on
!!$     if(sp==1) then
!!$        phase=>ptl_phase
!!$        pnum=>ptl_num
!!$        gid=>ptl_gid
!!$     else
!!$        phase=>ptl_electron
!!$        pnum=>ptl_e_num
!!$        gid=>ptl_e_gid
!!$     endif
!!$     
!!$     
!!$     do i=1, pnum
!!$        if(gid(i)>0) then
!!$           mu=phase(i,5)
!!$           tot_mu(sp)=tot_mu(sp)+mu
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$  call my_mpi_reduce(tot_mu,dum,2)
!!$
!!$  if(sml_mype==0) then
!!$     write(fnum,*) dum
!!$  end if
!!$
!!$end subroutine diagnosis_tot_mu
!!$
!!$subroutine diagnosis_tot_pnum(fnum)
!!$  use sml_module
!!$  use ptl_module
!!$  implicit none
!!$  include 'mpif.h'
!!$  integer :: fnum
!!$  integer :: i,sp,ierror
!!$  integer, pointer :: pnum,gid(:)
!!$  real (kind=8), pointer :: phase(:,:)
!!$  real (kind=8) :: tot_pnum(2),dum(2)
!!$
!!$  tot_pnum=0D0
!!$
!!$  do sp=1, 1+sml_electron_on
!!$     if(sp==1) then
!!$        phase=>ptl_phase
!!$        pnum=>ptl_num
!!$        gid=>ptl_gid
!!$     else
!!$        phase=>ptl_electron
!!$        pnum=>ptl_e_num
!!$        gid=>ptl_e_gid
!!$     endif
!!$     
!!$     
!!$     do i=1, pnum
!!$        if(gid(i)>0) then
!!$           tot_pnum(sp)=tot_pnum(sp)+real(1D0)
!!$        else
!!$           print *, ' minus region',fnum,sml_mype
!!$           stop
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$  call my_mpi_reduce(tot_pnum,dum,2)
!!$
!!$  if(sml_mype==0) then
!!$     write(fnum,*) dum
!!$  end if
!!$
!!$end subroutine diagnosis_tot_pnum
!!$
!!$


! For debug only
subroutine diag_efield(istep,grid,psn)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,istep,iphi,nd(3)
  character (len=30) :: filename
  real (kind=8) :: x(2), E(2)
#ifdef ADIOS
  integer*8 :: buf_id,grp_id
#endif
  if (mod(istep,diag_pot_period)==0) then
     ! perp efield
     write(filename,'("fort.ef",".",i2.2)') (istep)/(diag_pot_period)
     open(500, file=filename, status='replace', form='formatted') 
     do i=1, grid%ntriangle
        nd=grid%nd(:,i)
        x= 1D0/3D0 * ( grid%x(:,nd(1)) + grid%x(:,nd(2)) + grid%x(:,nd(3)) ) !center of a triangle
        E(1:2)= psn%E_perp_tr(:,i,1)
        write(500,1000) x(1), x(2),&
           E(1), &
           E(2), &
           sqrt( E(1)**2 + E(2)**2 ) 
     enddo
     close(500)
  endif
1000 format(5(e19.13,' '))

end subroutine diag_efield

subroutine charging_test(istep,grid,psn)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,istep,iphi
  real (kind=8) :: sum1(sml_nphi_total),tmp(sml_nphi_total)

  if(sml_mype==0) then
     do i=1, grid%nnode
        if(psn%idensity0(i)/=0D0) then
           write(1222,1000) sml_time/sml_tran,grid%x(1,i),grid%x(2,i), psn%idensity0(i)
        endif
     enddo
  endif
  write(1222,*) ' '

  do i=1, grid%nphi
     sum1(i+grid%iphi_offset)= sum(psn%idensity(:,i))
  enddo
  
  call my_mpi_reduce(sum1,tmp,grid%nphi)

  do i=1, sml_nphi_total
     write(1333,2000) sml_time/sml_tran,real(i),tmp(i)
  enddo

  write(1333,*) ' '

1000 format(4(e19.13,' '))
2000 format(3(e19.13,' '))
  
end subroutine charging_test
