!! Obtain charge density from gyro-center particle
subroutine chargei(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use smooth_module
  use perf_monitor
  implicit none
  include 'mpif.h'

  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
!  integer, intent(in) :: iflag
  integer :: i,j,larmor,itr,k
  real (kind=8) :: rho,x_ring(2),p(3)
  real (kind=8) :: dx_unit(2,sml_nlarmor)
  real (kind=8) :: particle_weight !! weight of single particle
  real (kind=8) :: phi,phi_weight(2),sumtmp
  integer :: nodes(3)
  integer :: iphi,iphi_frac
  real (kind=8), external :: gyro_radius, init_den
  integer :: icount, idest,isource,isendtag,irecvtag,ierror
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (kind=8) :: psi,den, x(2)
  real (kind=8) :: inv_delta_phi,inv_nphi_total,inv_nlarmor, angle, ranx, dx_unit2(2),cosa,sina
  integer :: init,count
  real (kind=8) :: dpdr, dpdz, dp, psi_interpol
  logical, parameter :: USE_SEARCH_TR2 = .true.

#ifdef CIRCULAR_SPECIAL
  call chargei_cc(grid,psn,sp)
  return
#endif

  inv_delta_phi=1D0/grid%delta_phi
  inv_nphi_total=1D0/real(sml_nphi_total)
  inv_nlarmor=1D0/real(sml_nlarmor)

  !dx_unit(:,1)=(/1.,0./)
  !dx_unit(:,2)=(/0.,1./)
  !dx_unit(:,3)=(/0.,-1./)
  !dx_unit(:,4)=(/-1.,0./)
  do larmor=1, sml_nlarmor
     angle=sml_2pi/real(sml_nlarmor)*real(larmor-1)
     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
  enddo

  psn%idensity=0D0
  !for all particle
  call monitor_start(CHARGEI_ALL_PART_)
  do i=1, sp%num
     if(sp%gid(i)<0) cycle
  ! for 4-point gyro ring
     !get proper toroidal angle index and weight
     if(sml_deltaf==1) then
        particle_weight=sp%phase(6,i)*sp%phase(8,i)
     else
        particle_weight=sp%phase(8,i) ! for full f simulation only
     endif
     phi=sp%phase(3,i)
     iphi= FLOOR(phi*inv_delta_phi)    - grid%iphi_offset
     if(iphi<0) then 
        print *, 'chargei Warning : iphi exceeded angle range-', iphi,phi, i,sp%gid(i),sml_mype
        call err_count
        iphi=0
     else if(iphi>=grid%nphi) then
        print *, 'chargei Warning : iphi exceeded angle range+', iphi,phi, i,sp%gid(i),sml_mype
        call err_count
        iphi=grid%nphi-1
     endif

!     iphi= min(grid%nphi-1,max(iphi,0))
     phi_weight(2)= (phi*inv_delta_phi)   - iphi - grid%iphi_offset ! larger index weight
     phi_weight(1)=1D0 - phi_weight(2)  ! smaller index weight
#ifdef XGC_DEBUG5
     if(phi_weight(1)>1D0 .or. phi_weight(1)<0D0) then
        print *, phi_weight(1), grid%iphi_offset, iphi, phi*inv_delta_phi, sml_mype
        stop
     endif
#endif     
     do iphi_frac=1, sml_bfollow+1
        ! get field following posision
        
        if(sml_bfollow==1) then
!           call field_following_pos(sp%phase(1:2,i),phi, iphi_frac,phi_weight(1),grid%delta_phi,x)
           call field_following_pos(sp%phase(1:2,i),phi,iphi_frac,phi_weight(iphi_frac),grid%delta_phi,x)
        else
           x=sp%phase(1:2,i)
        endif

        rho=gyro_radius(x,sp%phase(5,i))  !gyro radius
#ifdef GYRO_RADIAL_AVG
        dpdr=psi_interpol(x(1),x(2),1,0)
        dpdz=psi_interpol(x(1),x(2),0,1)
        dp=sqrt(dpdr**2 + dpdz**2)
        cosa=dpdr/dp
        sina=dpdz/dp
#else
        angle=sml_2pi*ranx()
        cosa=cos(angle)
        sina=sin(angle)
#endif
        do larmor=1, sml_nlarmor
           
           dx_unit2(1)= dx_unit(1,larmor)*cosa+dx_unit(2,larmor)*sina
           dx_unit2(2)=-dx_unit(1,larmor)*sina+dx_unit(2,larmor)*cosa

           x_ring = x + rho* dx_unit2(:)
           
           ! find position for 
           if (USE_SEARCH_TR2) then
              call search_tr2(grid,x_ring,itr,p)
           else
              if(larmor==1) then
                 call search_tr(grid,x_ring,itr,p)
                 init=itr
              else
                 call search_tr_with_guess(grid,x_ring,init,itr,p,count)
              endif

           endif
           sp%tr_save(larmor,iphi_frac,i)=itr
           sp%p_save(:,larmor,iphi_frac,i)=p(:)
        enddo
     enddo

     if(sml_bfollow==0) then
        sp%tr_save(:,2,i)=sp%tr_save(:,1,i)
        sp%p_save(:,:,2,i)=sp%p_save(:,:,1,i)
     endif

     if(minval(sp%tr_save(:,:,i)) > 0 ) then
        do iphi_frac=1, 2
           do larmor=1, sml_nlarmor
              itr=sp%tr_save(larmor,iphi_frac,i)
              p=sp%p_save(:,larmor,iphi_frac,i)
              nodes=grid%nd(:,itr)
              
              do j=1, 3
                 psn%idensity(nodes(j),iphi+iphi_frac-1)= psn%idensity(nodes(j),iphi+iphi_frac-1)+&
                      p(j)*particle_weight*phi_weight(iphi_frac)
                 !              min_node=min(nodes(j),min_node)
                 !              max_node=max(nodes(j),max_node)
              enddo
           enddo
        enddo
     else
        !eliminate particle
        call remove_particle(sp,i,-1) 
     endif
  enddo
  call monitor_stop(CHARGEI_ALL_PART_)
  !1. toroidal
  
  ! send and receive data from other PE
  call monitor_start (CHARGEI_SR_)
  grid%rtmp1=psn%idensity(:,0)
  grid%rtmp2=0D0
  icount=grid%nnode
  idest=mod(sml_mype-sml_pe_per_plane+sml_totalpe,sml_totalpe)
  isource=mod(sml_mype+sml_pe_per_plane,sml_totalpe)
  isendtag=sml_mype
  irecvtag=isource
  
  
  
  call mpi_sendrecv(grid%rtmp1,icount,MPI_REAL8,idest,isendtag,&
       grid%rtmp2,icount,MPI_REAL8,isource,irecvtag,MPI_COMM_WORLD,istatus,ierror)
  
  psn%idensity(:,grid%nphi) = psn%idensity(:,grid%nphi) + grid%rtmp2(:)
  call monitor_stop (CHARGEI_SR_)

!***** moved to 'after-divede-by-vol'
  !Poloidal smoothing of ion density
!  do i=1, grid%nphi
!     call smooth_pol0(grid,psn%idensity(:,i),smoothH)
!  enddo


  !zero-zero mode extraction -- How?
  ! sum up all density value in a PE
  call monitor_start (Z_Z_MODE_EXT_)
  do i=1, grid%nnode
     psn%idensity0(i)=sum(psn%idensity(i,1:grid%nphi))
  enddo
  ! sum-up
  call my_mpi_allreduce(psn%idensity0,grid%rtmp1,grid%nnode) !rtmp1 is used for a temporory purpose. variable
  ! get density from charge summation
  psn%idensity0=grid%rtmp1*inv_nphi_total*grid%inv_node_vol(:)*inv_nlarmor
  call monitor_stop (Z_Z_MODE_EXT_)
  
  !average out density0 following the field line
  !do i=1, grid%npsi
  !   sumtmp=0D0
  !   do j=1, grid%ntheta(i)
  !      sumtmp=sumtmp+psn%idensity0(grid%itheta0(i)+j-1)
  !   enddo
  !   do j=1, grid%ntheta(i)
  !      psn%idensity0(grid%itheta0(i)+j-1)=sumtmp / grid%surf_vol(i)
  !   enddo
  !enddo


  ! get density from charge accumulation  
  call monitor_start (CHARGEI_ACCUM_)
  if(sml_pe_per_plane==1) then
     do i=1, grid%nphi
        psn%idensity(:,i)=psn%idensity(:,i)*grid%inv_node_vol(:)*inv_nlarmor
     enddo
  else
     call mpi_allreduce(psn%idensity(:,1),grid%rtmp1,grid%nnode,MPI_DOUBLE_PRECISION,mpi_sum,sml_plane_comm,ierror)
#ifdef XGC_DEBUG5
     grid%rtmp2=psn%idensity(:,1)*inv_nlarmor
#endif
     psn%idensity(:,1)=grid%rtmp1(:)*grid%inv_node_vol(:)*inv_nlarmor
  end if

 
  !Poloidal smoothing of ion density
  do i=1, grid%nphi
     call smooth_pol0(grid,psn%idensity(:,i),smoothH)
  enddo
 
  ! radial smoothing
  do i=1, grid%nphi
     call smooth_r(psn%idensity(:,i),grid%rtmp1,smooth_r1,grid)
     psn%idensity(:,i)=grid%rtmp1
  enddo
  call smooth_r(psn%idensity0(:),grid%rtmp1,smooth_r1,grid)
  psn%idensity0(:)=grid%rtmp1
  
  call monitor_stop (CHARGEI_ACCUM_)
end subroutine chargei


!! Obtain charge density from gyro-center particle
subroutine chargee(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use smooth_module
  use ptl_module
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
!  integer, intent(in) :: iflag
  integer :: i,j,itr
  real (kind=8) :: x(2),p(3)
  real (kind=8) :: particle_weight !! weight of single particle
  integer :: nodes(3)
  integer :: iphi
  real (kind=8), external ::  init_den
!  integer :: icount, idest,isource,isendtag,irecvtag,istatus,ierror
  real (kind=8) :: psi,den, inv_nphi_total

  logical, parameter :: USE_SEARCH_TR2 = .true.
  integer, parameter :: wall_num = 100

  inv_nphi_total=1D0/real(sml_nphi_total)
  psn%edensity0=0D0

  !for all particle
  call monitor_start (CHARGEE_ALL_PART_)
  do i=1, sp%num
     if(sp%gid(i)<0) cycle
     x = sp%phase(1:2,i)
        
     ! find position and save it.
     if (USE_SEARCH_TR2) then
       call search_tr2(grid,x,itr,p)
     else
       call search_tr(grid,x,itr,p)
     endif

     sp%tr_save(1,1,i)=itr
     sp%p_save(:,1,1,i)=p(:)

     ! particle wall hit check

     if(itr<0) then
        if(sml_sheath_mode==0) then
           call remove_particle(sp,i,-1)
           cycle ! return to loop start for next particle
        else
           ! sheath calculation 
           call sheath_calculation(grid,psn,sp,i,1,itr,p)
        endif
        if(sp%gid(i)<0) then
           cycle  ! return to loop start for next particle
        else
           sp%tr_save(1,1,i)=itr
           sp%p_save(:,1,1,i)=p(:)
        endif                  
     endif

     if(sml_deltaf==1) then
        particle_weight=sp%phase(8,i)*sp%phase(6,i)
     else
        particle_weight=sp%phase(8,i) ! for full f simulation only
     endif
     
     nodes=grid%nd(:,itr)
     
     do j=1, 3
        psn%edensity0(nodes(j))= psn%edensity0(nodes(j))+p(j)*particle_weight
     enddo
  enddo
  call monitor_stop (CHARGEE_ALL_PART_)

  call monitor_start (CHARGEE_DENSITY_)
  ! sum-up
  call my_mpi_allreduce(psn%edensity0,grid%rtmp1,grid%nnode) !sendl is used for a temporory purpose. variable
  ! get density from charge accumulation

  psn%edensity0=grid%rtmp1*grid%inv_node_vol(:)*inv_nphi_total

  ! radial smoothing
  call smooth_r(psn%edensity0(:),grid%rtmp1,smooth_r1,grid)
  psn%edensity0(:)=grid%rtmp1
  call monitor_stop (CHARGEE_DENSITY_)
  
end subroutine chargee


subroutine smooth_pol0(grid,data,smooth)
  use grid_class  
  use sml_module
  use smooth_module
  use eq_module, only : eq_x_psi
  use perf_monitor
  implicit none
  type(smooth_type) :: smooth
  type(grid_type) :: grid
  real (kind=8) ,intent(inout) :: data(grid%nnode)
  integer :: nd1,nd2, nd,rl,i,cnd
  real (kind=8) :: dsum(2),wsum(2),s,a
  
  ! mode | smoothing
  !  -1     no smoothing
  !  0      flux average for all
  !  1      flux average for inside separatrix
  !  2      smoothing with smoothing function which is decided by smoothing type

  ! no smooth for -1 mode
  if(smooth%mode==-1) then
     return
  endif

  !flux averaging -- for mode 0 and 1
  ! for mode 1, psi <= eq_x_psi
  call monitor_start(SMOOTH_)
  if(smooth%mode==0 .or. smooth%mode==1) then !flux average for mode==0 .or. mode==1
     do i=1, grid%npsi - 1 ! -1 : not for wall boundary 
        nd1=grid%itheta0(i) 
        nd2=grid%itheta0(i)+grid%ntheta(i)-1


        if(smooth%mode==0 .or. grid%psi(nd1) < 0.999D0*eq_x_psi .and. smooth%mode==1) then
!           Volume (flux) average        
!           s=sum( data(nd1:nd2)*grid%node_vol(nd1:nd2) )
!           data(nd1:nd2)=s/grid%surf_vol(i)*real(sml_nphi_total)
!           Area average           
           s=sum( data(nd1:nd2)*grid%node_area(nd1:nd2) )
           a=sum( grid%node_area(nd1:nd2) )
           data(nd1:nd2)=s/a
        else
           exit
        endif
     enddo
  else  ! set start node for smooth
     nd1=1
  endif

  if(smooth%mode==1 .or. smooth%mode==2) then
     do nd=nd1, grid%nnode
        dsum=0D0
        wsum=0D0
        do rl=1, 2
           cnd=nd
           NODES : do i=1,smooth%n
              dsum(rl)=dsum(rl)+smooth%weight(i)*data(cnd)*grid%node_area(cnd)
              wsum(rl)=wsum(rl)+smooth%weight(i)*grid%node_area(cnd)
              !find next node
              cnd=grid%nn(rl,cnd)
              if(cnd <=0 .or. cnd > grid%nnode) exit NODES
           enddo NODES
        enddo
        grid%rtmp1(nd)=0.5D0*(dsum(1)/wsum(1) + dsum(2)/wsum(2))
     enddo
     
     data(nd1:grid%nnode)=grid%rtmp1(nd1:grid%nnode)
  endif

  call monitor_stop (SMOOTH_)
end subroutine

subroutine init_bfollow(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  include 'mpif.h'
  integer :: i,j
  type(psn_type) :: psn
  type(grid_type) :: grid
  real (kind=8), allocatable :: xp(:,:),xm(:,:)
  integer , allocatable :: nold(:)
  integer :: itr,dir,min_dist_node,nd,node_new,node_old,dum
  real (kind=8) :: min_dist,p(3),x(2),xtmp(4),dist,dphi,phi
  integer :: iseg_offset,segnum, segend, segpe, ierror

  logical, parameter :: USE_SEARCH_TR2 = .true.

  !allocate memeory
  allocate(nold(grid%nnode), xp(2,grid%nnode), xm(2,grid%nnode) )
  nold=0 !safty
  xp=0D0 !safty
  xm=0D0 !safty
  
!  print *, 'reading b-following file'
  ! file open
  
  if(sml_bfollow==1 .and. sml_bfollow_read==1) then
     open(unit=3000,file=sml_bfollow_file,action='read')
     
     ! set default value
     nold=-1
     ! read line until hit -1
     loop1: do while(.true.)
        read(3000,1000) dum, node_old, xtmp(1),xtmp(2),xtmp(3),xtmp(4)
        if(dum==-1 ) then
           exit loop1
        endif
        
        node_new=grid%old2new(node_old)
        if (node_new .lt. 1 .or. node_new .gt. grid%nnode) then
           print*,'bogus node from old2new map: ',node_old,' -> ',node_new
        else
           xp(:,node_new)=xtmp(1:2)
           xm(:,node_new)=xtmp(3:4)
           nold(node_new)=node_old
        endif
     enddo loop1
     close(3000)
  elseif(sml_bfollow==1 .and. sml_bfollow_read /= 1) then
     ! Generates bfollow
     nold=0
     phi=0D0 ! Cannot specifiy toroidal angle - symmetry is assumed.
     do i=1, grid%nnode
        dir=1
        if(sml_bt_sign==-1) dir=-dir
        if(sml_minusB==1) dir=-dir

        dphi=real(dir)*grid%delta_phi
        ! field following direction
        call field_following_pos2(grid%x(:,i),phi,dphi,xp(:,i))

        ! opposite direction to B
        dphi=-dphi
        call field_following_pos2(grid%x(:,i),phi,dphi,xm(:,i))
        
     enddo
  else
     xp(:,:)=grid%x(:,:)
     xm(:,:)=grid%x(:,:)
     nold(:)=0
  endif
  
  !  print *, 'end b-following file'

  do i=1, grid%nnode
     if(nold(i)==-1) then
        xp(:,i)=grid%x(:,i)
        xm(:,i)=grid%x(:,i)
     endif
  enddo
#if defined(XGC_DEBUG2)
  call bfollow_test(xp,xm,grid)
#endif 
  !Parallizing below routine
  ! algorithm - each pe find each segment
  
  segnum = grid%nnode/sml_totalpe
  if(grid%nnode > sml_totalpe*segnum ) then
     segnum=segnum+1
  endif
  iseg_offset = sml_mype*segnum
  segend = min(grid%nnode, iseg_offset + segnum)
  

  
  do i=1+iseg_offset, segend
     !     print *, i
     psn%bfollow_1_dx(i)=1D0/sqrt( (xp(1,i)-xm(1,i))**2 + (xp(2,i)-xm(2,i))**2 + (grid%x(1,i)*grid%delta_phi*2D0)**2)
     ! What is exact formula for curved derivative of curved system ??
     ! Here we assumed that it is the same with XYZ derivative of XYZ coordinate
     
     ! for two directions
     do dir=1,2
        if(dir==1) then
           x=xp(:,i)
        else
           x=xm(:,i)
        endif
        !        print *, dir , 'search tr',x
        
        if (USE_SEARCH_TR2) then
           call search_tr2(grid,x,itr,p)
        else
           call search_tr(grid,x,itr,p)
        endif
        
        !find  node and p
        !        print *, dir, 'search tr end'
        if(itr<=0) then
           !find nearest node point
           !           print *, 'find_nearest',dir
           min_dist=1D50
           do nd=1, grid%nnode
              dist=(grid%x(1,nd)-x(1))**2 + (grid%x(2,nd)-x(2))**2
              if(min_dist > dist) then
                 min_dist=dist
                 min_dist_node=nd
              endif
           enddo
           !           print *, 'min',min_dist, min_dist_node
           !choose one triangle that has this node point
           itr=grid%tr_node(1,min_dist_node)
           do j=1,3
              if(min_dist_node==grid%nd(j,itr)) then
                 p(j)=1D0
              else
                 p(j)=0D0
              endif
           enddo
        endif
        if(itr<=0 .or. itr> grid%ntriangle) then
           print *, 'Wrong itr', itr,grid%ntriangle
           call err_count
        endif
        psn%bfollow_tr(dir,i) = itr
        psn%bfollow_p(:,dir,i) = p
     enddo
  enddo
  
  !propagate result
  do segpe=0, sml_totalpe-1
     iseg_offset = segpe*segnum
     segend = min(grid%nnode, iseg_offset + segnum)
     if(segend-iseg_offset >0 ) then
         call MPI_BCAST(psn%bfollow_tr(:,iseg_offset+1:segend) , (segend-iseg_offset)*2  ,&
              MPI_INTEGER, segpe, MPI_COMM_WORLD, ierror)
         call MPI_BCAST(psn%bfollow_p(:,:,iseg_offset+1:segend) , (segend-iseg_offset)*6  ,&
              MPI_REAL8, segpe, MPI_COMM_WORLD, ierror)
         call MPI_BCAST(psn%bfollow_1_dx(iseg_offset+1:segend), (segend-iseg_offset) ,&
              MPI_REAL8, segpe, MPI_COMM_WORLD, ierror) 
     endif
  enddo
#if defined(XGC_DEBUG2)
!  call bfollow_test(xp,xm,grid)
#endif 
  !deallocate memeory
  deallocate(nold,xp,xm) 

1000 format (2I8,1x,4(e19.13,1x))

end subroutine init_bfollow

subroutine bfollow_test(xp,xm,grid)
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  real (kind=8) :: xp(2,grid%nnode), xm(2,grid%nnode)
  integer :: i
  real (kind=8) :: dx(3), r,z,phi,b(3),cost1, cost2, cost3
  
  if(sml_mype==0) then

     do i=1, grid%nnode
        r=grid%x(1,i)
        z=grid%x(2,i)
        dx(1:2)= xp(:,i)-xm(:,i)
        dx(3)=- 2D0*r *   grid%delta_phi ! r * delta phi
        call bvec_interpol(r,z,phi,b(1),b(2),b(3))
        
        cost1 = (dx(1)*b(1)+dx(2)*b(2)+dx(3)*b(3))/sqrt( (dx(1)**2+dx(2)**2+dx(3)**2) * (b(1)**2+b(2)**2+b(3)**2) )
        cost2 = (-dx(1)*b(1)-dx(2)*b(2)+dx(3)*b(3))/sqrt( (dx(1)**2+dx(2)**2+dx(3)**2) * (b(1)**2+b(2)**2+b(3)**2) )
        cost3 = (dx(3)*b(3))/sqrt( (dx(3)**2) * (b(1)**2+b(2)**2+b(3)**2) )
        write(1357,1000) i, r,z,1D0-cost1, 1D0-cost2, 1D0-cost3
        write(1358,1000) i, r,z, xp(1,i), xp(2,i), xm(1,i), xm(2,i)
     enddo
     close(1357)
  endif

1000 format (I6,1x,6(e19.13,1x))

end subroutine

subroutine chargee_background(grid,psn,ptl)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use smooth_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(ptl_type) :: ptl
  integer :: i
  real (kind=8) :: psi,init_den

  if(sml_deltaf==1) then
     psn%edensity0=0D0
  else
     call chargei(grid,psn,ptl%ion)
     psn%edensity0(:)=psn%idensity0(:)
     if(sml_flat_electron_density==1) then
        call smooth_pol0(grid, psn%edensity0, smooth00)  !???? how? call it or not
     endif
!     do i=1, grid%nnode
!        if( (grid%inner_bd(1)> i .or. i> grid%inner_bd(2)) .and. & 
!             (grid%outer_bd(1) > i .or. i> grid%outer_bd(2)) ) then
!           psi=grid%psi(i)
!           psn%edensity0(i)=init_den(psi,grid%x(2,i))
!        else
!           psn%edensity0(i)=0D0
!        endif
!     enddo
  endif
     
end subroutine chargee_background


subroutine field_following_pos_org(x_org,phi,dir,w,delta_phi,x_dest)
  implicit none
  real (kind=8), intent(in) :: x_org(2), w, delta_phi,phi
  real (kind=8), intent(out) :: x_dest(2)
  integer :: dir
  real (kind=8) :: B(2),Bphi,dphi, x_mid(2)
  real (kind=8) :: hh,h6,dx1(2),dx2(2),dx3(2),x_tmp(2)

  if(dir==2) then   ! dir==2 -> larger index. To larger phi plane. dphi >0
     dphi=delta_phi*(1D0-w)
  else              ! dir==1 -> smaller index. To smaller phi plane. dphi <0
     dphi=delta_phi*(w-1D0)
  endif
  
  call bvec_interpol(x_org(1),x_org(2),0D0,b(1),b(2),bphi)

  if( .false. ) then ! first order calculation
     x_dest(:)=x_org(:) + b(:)/bphi*(x_org(1)*dphi)
  else if( .true. ) then  
     ! second order calculation - rk2
     ! obtain mid point
     x_mid(:)=x_org(:) + b(:)/bphi*(x_org(1)*dphi*0.5D0) 
     
     ! get new derivative
     call bvec_interpol(x_mid(1),x_mid(2),0D0,b(1),b(2),bphi)

     ! advance one step using mid-point derivative
     x_dest(:)=x_org(:) + b(:)/bphi*(x_mid(1)*dphi) 
  else
     ! 4th Order Calculation - rk4
     !
     hh=dphi*0.5D0
     h6=dphi/6D0     
     ! derivative 1 (from x_org)- obtained already
     dx1=b/bphi*x_org(1)     ! dydx from y       : y-> x_org , dx1 -> dydx
     x_tmp=x_org + hh*dx1    ! yt=y+hh*dydx      : yt -> x_tmp
     !
     ! derivative 2 (from x_tmp) 
     call bvec_interpol(x_tmp(1),x_tmp(2),0D0,b(1),b(2),bphi)
     dx2=b/bphi*x_tmp(1)     ! dyt from yt       : dyt -> dx2
     x_tmp=x_org + hh*dx2    ! yt=y+hh*dyt       : yt -> x_tmp
     !
     ! derivative 3 (from x_tmp)
     call bvec_interpol(x_tmp(1),x_tmp(2),0D0,b(1),b(2),bphi)
     dx3=b/bphi*x_tmp(1)     ! dym from yt       : dym -> dx3 
     x_tmp=x_org + dphi*dx3     ! yt=y + h*dym      : yt -> x_tmp
     dx3 = dx2 + dx3         ! dym = dyt + dym   : dym -> dx3 , dyt -> dx2
     !
     ! derivative 4 (from x_tmp)
     call bvec_interpol(x_tmp(1),x_tmp(2),0D0,b(1),b(2),bphi)
     dx2=b/bphi*x_tmp(1)    ! dyt from yt       : dyt -> dx2, yt -> x_tmp
     x_dest=x_org+h6* ( dx1 + dx2 + 2D0*dx3)     ! yout = y + h6* (dydx+dyt+2D0*dym)  : yout -> x_dest, dydx -> dx1, dyt-> dx2 , dym -> dx3
  endif
     
  !x_dest(:)=x_org(:)
end subroutine field_following_pos_org


subroutine field_following_pos(x_org,phi,dir,w,delta_phi,x_dest)
  implicit none
  real (kind=8), intent(in) :: x_org(2), w, delta_phi,phi
  real (kind=8), intent(out) :: x_dest(2)
  real (kind=8) :: dphi
  integer :: dir
  
  if(dir==2) then   ! dir==2 -> larger index. To larger phi plane. dphi >0
     dphi=delta_phi*(1D0-w)
  else              ! dir==1 -> smaller index. To smaller phi plane. dphi <0
     dphi=delta_phi*(w-1D0)
  endif
  
  call field_following_pos2(x_org,phi,dphi,x_dest)
  
end subroutine field_following_pos

subroutine field_following_pos2(x_org,phi0,dphi,x_dest)
  implicit none
  real (kind=8), intent(in) :: x_org(2), dphi,phi0
  real (kind=8), intent(out) :: x_dest(2)
  real (kind=8) :: B(2),Bphi, x_mid(2)
  real (kind=8) :: hh,h6,dx1(2),dx2(2),dx3(2),x_tmp(2)

  call bvec_interpol(x_org(1),x_org(2),phi0,b(1),b(2),bphi)

  if( .false. ) then ! first order calculation
     x_dest(:)=x_org(:) + b(:)/bphi*(x_org(1)*dphi)
  else if( .true. ) then  
     ! second order calculation - rk2
     ! obtain mid point
     hh=dphi*0.5D0
     x_mid(:)=x_org(:) + b(:)/bphi*(x_org(1)*hh) 
     
     ! get new derivative
     call bvec_interpol(x_mid(1),x_mid(2),phi0+hh,b(1),b(2),bphi)

     ! advance one step using mid-point derivative
     x_dest(:)=x_org(:) + b(:)/bphi*(x_mid(1)*dphi) 
  else
     ! 4th Order Calculation - rk4
     !
     hh=dphi*0.5D0
     h6=dphi/6D0     
     ! derivative 1 (from x_org)- obtained already
     dx1=b/bphi*x_org(1)     ! dydx from y       : y-> x_org , dx1 -> dydx
     x_tmp=x_org + hh*dx1    ! yt=y+hh*dydx      : yt -> x_tmp
     !
     ! derivative 2 (from x_tmp) 
     call bvec_interpol(x_tmp(1),x_tmp(2),phi0+hh,b(1),b(2),bphi)
     dx2=b/bphi*x_tmp(1)     ! dyt from yt       : dyt -> dx2
     x_tmp=x_org + hh*dx2    ! yt=y+hh*dyt       : yt -> x_tmp
     !
     ! derivative 3 (from x_tmp)
     call bvec_interpol(x_tmp(1),x_tmp(2),phi0+hh,b(1),b(2),bphi)
     dx3=b/bphi*x_tmp(1)     ! dym from yt       : dym -> dx3 
     x_tmp=x_org + dphi*dx3     ! yt=y + h*dym      : yt -> x_tmp
     dx3 = dx2 + dx3         ! dym = dyt + dym   : dym -> dx3 , dyt -> dx2
     !
     ! derivative 4 (from x_tmp)
     call bvec_interpol(x_tmp(1),x_tmp(2),phi0+dphi,b(1),b(2),bphi)
     dx2=b/bphi*x_tmp(1)    ! dyt from yt       : dyt -> dx2, yt -> x_tmp
     x_dest=x_org+h6* ( dx1 + dx2 + 2D0*dx3)     ! yout = y + h6* (dydx+dyt+2D0*dym)  : yout -> x_dest, dydx -> dx1, dyt-> dx2 , dym -> dx3
  endif
     
  !x_dest(:)=x_org(:)
end subroutine field_following_pos2



subroutine smooth_r_init2(smooth_r,grid)
  use smooth_module
  use grid_class
  implicit none
  type(smooth_r_type) :: smooth_r
  type(grid_type) :: grid
  real (kind=8) :: vring(smooth_r%n), fring(0:smooth_r%n)
  real (kind=8) :: weight,p(3),dpsi(2),r_norm(2),psi_interpol,x0(2),x(2),tmp
  integer :: i,j,kr,dir,itr,nodes(3)

  if(smooth_r%n<=0) then
     return ! no smoothing
!  else if(smooth_r%n==1) then
!     print *, 'Error : invalid smooth_r%n. Please use smooth_r%n=0 for no smoothing', smooth_r%n
!     stop
!  else if(smooth_r%n==2) then
!     vring(1)=-0.707107D0
!     vring(2)= 0.707107D0
!     fring(1)= 0.5D0
!     fring(2)= 0.5D0
  else if(smooth_r%n==1) then
     vring(1)=1.22474D0
     fring(0)=2D0/3D0
     fring(1)=1D0/6D0
!     vring(1)=-1.22474D0
!     vring(2)= 0D0
!     vring(3)= 1.22474D0
!     fring(2)= 2D0/3D0 
!     fring(1)= 1D0/6D0
!     fring(3)= 1D0/6D0
!  else if(smooth_r%n==4) then
!     vring(1)= -1.65068D0
!     vring(2)= -0.524648D0
!     vring(3)=  0.524648D0
!     vring(4)=  1.65068D0
!     fring(1)=  (3D0-sqrt(6D0))/12D0
!     fring(2)=  (3D0+sqrt(6D0))/12D0
!     fring(3)=  (3D0+sqrt(6D0))/12D0
!     fring(4)=  (3D0-sqrt(6D0))/12D0
  else if(smooth_r%n==2) then
     vring(1)= 0.958572D0
     vring(2)= 2.02018D0
     fring(0)= 8D0/15D0
     fring(1)= (7D0+2D0*sqrt(10D0))/60
     fring(2)= (7D0-2D0*sqrt(10D0))/60
!     vring(1)= -2.02018D0
!     vring(2)= -0.958572D0
!     vring(3)= 0D0
!    vring(4)= 0.958572D0
!     vring(5)= 2.02018D0
!     fring(1)= (7D0-2D0*sqrt(10D0))/60
!     fring(2)= (7D0+2D0*sqrt(10D0))/60
!     fring(3)= 8D0/15D0
!     fring(4)= (7D0+2D0*sqrt(10D0))/60
!     fring(5)= (7D0-2D0*sqrt(10D0))/60
  else
     !simple linear gauss smoothing
     do i=1, smooth_r%n
        vring(i)= real(2.2*i)/real(smooth_r%n)
        fring(i)= exp( - vring(i)**2 )
     enddo
     fring(0)=1D0
     tmp=1D0 + 2D0 * sum(fring(1:smooth_r%n))
     fring=fring/tmp
     !print *, 'error : invalid smooth_r%n, not implimented yet', smooth_r%n
     !stop
  endif

  do i=1, grid%nnode

     ! 1st point is the original grid point
     call set_value(smooth_r%mat,i,i,fring(0),1)
     
     ! position of grid points
     x0=grid%x(:,i)

     ! B-field vector
     dpsi(1)=psi_interpol(x0(1),x0(2),1,0)
     dpsi(2)=psi_interpol(x0(1),x0(2),0,1)
     r_norm=dpsi / sqrt(dpsi(1)**2+dpsi(2)**2)

     do kr=1,smooth_r%n
        do dir=-1,1,2
           x=x0+r_norm*real(dir)*smooth_r%d0
           call search_tr(grid,x,itr,p)
           if(itr>0) then
              nodes=grid%nd(:,itr)
              
              do j=1,3
                 weight=fring(kr)*p(j)
                 if(weight<0.00001) then
                    call set_value(smooth_r%mat,i,i,weight,1)  ! add to diagonal
                 else
                    call set_value(smooth_r%mat,i,nodes(j),weight,1)
                 endif
              enddo !end of 3-point interpolation loop
           else ! cannot find smoothing grid
              call set_value(smooth_r%mat,i,i,fring(kr),1)
           endif
        enddo !end of direction loop
     enddo ! end of n ring loop
  end do ! end of node loop
  
  
end subroutine smooth_r_init2


subroutine smooth_r(in,out,smoothr,grid)
  use grid_class
  use smooth_module
  implicit none
  type(smooth_r_type) :: smoothr
  type(grid_type) :: grid
  real (kind=8) , dimension(grid%nnode) :: in,out
  integer :: i,j,ierr

  if(smoothr%n<=0) then
     out=in
     return
  end if

  do i=1, grid%nnode
     out(i)=0D0
     do j=1, smoothr%mat%nelement(i)
        out(i)=out(i)+smoothr%mat%value(j,i)*in(smoothr%mat%eindex(j,i))
     enddo
  enddo

end subroutine smooth_r

!******************************************************************************************************
!********  Funtions for debugging ********************************************************************
!******************************************************************************************************
subroutine enforce_modeled_charge(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn

  

  psn%idensity=1D15
  psn%idensity0=1D15


end subroutine enforce_modeled_charge

#ifdef CIRCULAR_SPECIAL
subroutine chargei_cc(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use smooth_module
  use perf_monitor
  implicit none
  include 'mpif.h'

  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
!  integer, intent(in) :: iflag
  integer :: i,j,larmor,itr,k
  real (kind=8) :: rho,x_ring(2),p(3)
  real (kind=8) :: dx_unit(2,sml_nlarmor)
  
  real (kind=8) :: particle_weight !! weight of single particle
  real (kind=8) :: phi,phi_weight(2),sumtmp
  integer :: nodes(3)
  integer :: iphi,iphi_frac
  real (kind=8), external :: gyro_radius, init_den
  integer :: icount, idest,isource,isendtag,irecvtag,ierror
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (kind=8) :: psi,den, x(2)
  real (kind=8) :: inv_delta_phi,inv_nphi_total,inv_nlarmor, angle, ranx, dx_unit2(2),cosa,sina
  integer :: init,count
  real (kind=8) :: dpdr, dpdz, dp, psi_interpol
  logical, parameter :: USE_SEARCH_TR2 = .true.
  ! for r-theta interpolation
  real (kind=8) :: rw(2),tw(2),r,theta
  integer :: ri,ti,rindex,tindex
  integer :: node

  inv_delta_phi=1D0/grid%delta_phi
  inv_nphi_total=1D0/real(sml_nphi_total)
  inv_nlarmor=1D0/real(sml_nlarmor)

  !dx_unit(:,1)=(/1.,0./)
  !dx_unit(:,2)=(/0.,1./)
  !dx_unit(:,3)=(/0.,-1./)
  !dx_unit(:,4)=(/-1.,0./)
  do larmor=1, sml_nlarmor
     angle=sml_2pi/real(sml_nlarmor)*real(i-1)
     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
  enddo

  psn%idensity=0D0
  !for all particle
  call monitor_start(CHARGEI_ALL_PART_)
  do i=1, sp%num
  ! for 4-point gyro ring
     !get proper toroidal angle index and weight
     if(sml_deltaf==1) then
        particle_weight=sp%phase(6,i)*sp%phase(8,i)
     else
        particle_weight=sp%phase(8,i) ! for full f simulation only
     endif
     phi=sp%phase(3,i)
     iphi= FLOOR(phi*inv_delta_phi)    - grid%iphi_offset
     if(iphi<0) then 
        print *, 'chargei Warning : iphi exceeded angle range-', iphi,phi, i,sp%gid(i),sml_mype
        call err_count
        iphi=0
     else if(iphi>=grid%nphi) then
        print *, 'chargei Warning : iphi exceeded angle range+', iphi,phi, i,sp%gid(i),sml_mype
        call err_count
        iphi=grid%nphi-1
     endif

!     iphi= min(grid%nphi-1,max(iphi,0))
     phi_weight(2)= (phi*inv_delta_phi)   - iphi - grid%iphi_offset ! larger index weight
     phi_weight(1)=1D0 - phi_weight(2)  ! smaller index weight
#ifdef XGC_DEBUG5
     if(phi_weight(1)>1D0 .or. phi_weight(1)<0D0) then
        print *, phi_weight(1), grid%iphi_offset, iphi, phi*inv_delta_phi, sml_mype
        stop
     endif
#endif     
     do iphi_frac=1, sml_bfollow+1
        ! get field following posision
        
        if(sml_bfollow==1) then
!           call field_following_pos(sp%phase(1:2,i),phi, iphi_frac,phi_weight(1),grid%delta_phi,x)
           call field_following_pos(sp%phase(1:2,i),phi,iphi_frac,phi_weight(iphi_frac),grid%delta_phi,x)
        else
           x=sp%phase(1:2,i)
        endif

        rho=gyro_radius(x,sp%phase(5,i))  !gyro radius
#ifdef GYRO_RADIAL_AVG
        dpdr=psi_interpol(x(1),x(2),1,0)
        dpdz=psi_interpol(x(1),x(2),0,1)
        dp=sqrt(dpdr**2 + dpdz**2)
        cosa=dpdr/dp
        sina=dpdz/dp
#else
        angle=sml_2pi*ranx()
        cosa=cos(angle)
        sina=sin(angle)
#endif
        do larmor=1, sml_nlarmor
           
           dx_unit2(1)= dx_unit(1,larmor)*cosa+dx_unit(2,larmor)*sina
           dx_unit2(2)=-dx_unit(1,larmor)*sina+dx_unit(2,larmor)*cosa

           x_ring = x + rho* dx_unit2(:)
           
           ! get psi and theta
           call get_r_theta(x_ring,r,theta)
           if(r>=grid%radial(grid%npsi)) then
              call remove_particle(sp,i,-1)
              goto 1000
           endif
           ! find r-index
           call search_rindex(grid,r,ri,rw(1))
           rw(2)=1D0-rw(1)
           ! save index and weight
           sp%rindex(larmor,iphi_frac,i)=ri
           sp%rweight(larmor,iphi_frac,i)=rw(1)
           do rindex=1,2              
              ! find theta- theta+ using binary search              
              call search_tindex(grid,theta,ri+rindex-1,ti,tw(1))
              tw(2)=1D0-tw(1)
              ! save index and weight
              sp%tindex(rindex,larmor,iphi_frac,i)=ti
              sp%tweight(rindex,larmor,iphi_frac,i)=tw(1)
              do tindex=1,2
                 call find_node(grid,ri+rindex-1,ti+tindex-1,node)
                 ! assign charge
                 if(node>=0) then
                    psn%idensity(node,iphi+iphi_frac-1)=psn%idensity(node,iphi+iphi_frac-1) + &
                         particle_weight*phi_weight(iphi_frac)*rw(rindex)*tw(tindex)
                 else
                    !eliminate particle
                    call remove_particle(sp,i,-1)
                    goto 1000 ! break do right before  particle loop
                 endif
              enddo
           enddo
        enddo
     enddo
     
     if(sml_bfollow==0) then
        sp%rindex(:,2,i)=sp%rindex(:,1,i)
        sp%rweight(:,2,i)=sp%rweight(:,1,i)
        sp%tindex(:,:,2,i)=sp%tindex(:,:,1,i)
        sp%tweight(:,:,2,i)=sp%tweight(:,:,1,i)
     endif
1000 continue
  enddo
  call monitor_stop(CHARGEI_ALL_PART_)
  !1. toroidal
  
  ! send and receive data from other PE
  call monitor_start (CHARGEI_SR_)
  grid%rtmp1=psn%idensity(:,0)
  grid%rtmp2=0D0
  icount=grid%nnode
  idest=mod(sml_mype-sml_pe_per_plane+sml_totalpe,sml_totalpe)
  isource=mod(sml_mype+sml_pe_per_plane,sml_totalpe)
  isendtag=sml_mype
  irecvtag=isource
  
  
  
  call mpi_sendrecv(grid%rtmp1,icount,MPI_REAL8,idest,isendtag,&
       grid%rtmp2,icount,MPI_REAL8,isource,irecvtag,MPI_COMM_WORLD,istatus,ierror)
  
  psn%idensity(:,grid%nphi) = psn%idensity(:,grid%nphi) + grid%rtmp2(:)
  call monitor_stop (CHARGEI_SR_)

!***** moved to 'after-divede-by-vol'
  !Poloidal smoothing of ion density
!  do i=1, grid%nphi
!     call smooth_pol0(grid,psn%idensity(:,i),smoothH)
!  enddo


  !zero-zero mode extraction -- How?
  ! sum up all density value in a PE
  call monitor_start (Z_Z_MODE_EXT_)
  do i=1, grid%nnode
     psn%idensity0(i)=sum(psn%idensity(i,1:grid%nphi))
  enddo
  ! sum-up
  call my_mpi_allreduce(psn%idensity0,grid%rtmp1,grid%nnode) !rtmp1 is used for a temporory purpose. variable
  ! get density from charge summation
  psn%idensity0=grid%rtmp1*inv_nphi_total*grid%inv_node_vol(:)*inv_nlarmor
  call monitor_stop (Z_Z_MODE_EXT_)
  
  !average out density0 following the field line
  !do i=1, grid%npsi
  !   sumtmp=0D0
  !   do j=1, grid%ntheta(i)
  !      sumtmp=sumtmp+psn%idensity0(grid%itheta0(i)+j-1)
  !   enddo
  !   do j=1, grid%ntheta(i)
  !      psn%idensity0(grid%itheta0(i)+j-1)=sumtmp / grid%surf_vol(i)
  !   enddo
  !enddo


  ! get density from charge accumulation  
  call monitor_start (CHARGEI_ACCUM_)
  if(sml_pe_per_plane==1) then
     do i=1, grid%nphi
        psn%idensity(:,i)=psn%idensity(:,i)*grid%inv_node_vol(:)*inv_nlarmor
     enddo
  else
     call mpi_allreduce(psn%idensity(:,1),grid%rtmp1,grid%nnode,MPI_DOUBLE_PRECISION,mpi_sum,sml_plane_comm,ierror)
#ifdef XGC_DEBUG5
     grid%rtmp2=psn%idensity(:,1)*inv_nlarmor
#endif
     psn%idensity(:,1)=grid%rtmp1(:)*grid%inv_node_vol(:)*inv_nlarmor
  end if

 
  !Poloidal smoothing of ion density
  do i=1, grid%nphi
     call smooth_pol0(grid,psn%idensity(:,i),smoothH)
  enddo
 
  ! radial smoothing
  do i=1, grid%nphi
     call smooth_r(psn%idensity(:,i),grid%rtmp1,smooth_r1,grid)
     psn%idensity(:,i)=grid%rtmp1
  enddo
  call smooth_r(psn%idensity0(:),grid%rtmp1,smooth_r1,grid)
  psn%idensity0(:)=grid%rtmp1
  
  call monitor_stop (CHARGEI_ACCUM_)
end subroutine chargei_cc

subroutine find_node(grid,ri,ti,node)
  use grid_class
  implicit none
  type(grid_type) :: grid
  integer, intent(in) :: ri,ti
  integer, intent(out) :: node
  integer :: ti2

  ti2=ti
  if(ti2 >= grid%itheta0(ri)+grid%ntheta(ri)) ti2=ti2-grid%ntheta(ri)
!  node= grid%sort(ti2)
  node= ti2

end subroutine find_node

subroutine get_r_theta(x,radial,theta)
  use eq_module
  use sml_module
  implicit none
  real (kind=8) :: x(2), radial, theta

  radial=dsqrt( (x(1)-eq_axis_r)**2 + (x(2)-eq_axis_z)**2 )
  theta=dacos( (x(1)-eq_axis_r)/radial )
  if( x(2) < eq_axis_z ) then
     theta=sml_2pi-theta
  endif

end subroutine get_r_theta

subroutine get_cc_values(grid)
  use sml_module
  use grid_class
  use eq_module
  implicit none
  type(grid_type) :: grid
  real (kind=8) :: radial,theta
  real (kind=8) :: r,z
  integer :: i,j,node
  integer, external :: angle_compare

  do i=1, grid%npsi

     do j=1, grid%ntheta(i)
        node=grid%itheta0(i)+j-1
        r=grid%x(1,node)
        z=grid%x(2,node)
        radial=dsqrt( (r-eq_axis_r)**2 + (z-eq_axis_z)**2 )
        theta=dacos( (r-eq_axis_r)/radial )
        if( z < eq_axis_z ) then
           theta=sml_2pi-theta
        endif
        
        grid%radial(i)=radial
        grid%theta(node)=theta        
     enddo
  enddo
  

  !check if it is increasing function

  do i=1, grid%npsi-1
     if(grid%radial(i) >= grid%radial(i+1)) then
        print *, 'Error in get_cc_values (r)', grid%radial(i),grid%radial(i+1),i
        stop
     endif
  enddo
  do i=1, grid%npsi
     do j=1, grid%ntheta(i)-1
        node=grid%itheta0(i)+j-1
        if(angle_compare(grid%theta(node),grid%theta(node+1),grid%theta(grid%itheta0(i)))/=-1) then
           print *, 'Error in get_cc_values (theta)', grid%theta(node),grid%theta(node+1),grid%theta(grid%itheta0(i)),node,i
           stop
        endif
     enddo
  enddo

end subroutine get_cc_values

integer function angle_compare(a,b,base)
  use sml_module
  implicit none
  real (kind=8), intent(in):: a,b,base
  real (kind=8) :: ab,bb 
  
  ab=a-base
  if(ab<0D0) ab=ab+sml_2pi
  bb=b-base
  if(bb<0D0) bb=bb+sml_2pi

  if( ab>bb) then
     angle_compare = 1
  elseif( ab<bb) then
     angle_compare = -1
  else
     angle_compare = 0
  endif

end function angle_compare

subroutine search_rindex(grid,r,ri,rw)
  use grid_class
  implicit none
  type(grid_type) :: grid
  real (kind=8), intent(in):: r
  real (kind=8), intent(out):: rw
  integer, intent(out):: ri
  integer :: l,u ! lower and upper bound
  integer :: mid

  !binary search
  l=1
  u=grid%npsi
  do while ( u - l > 1 )
     mid=(l+u)/2
     if(grid%radial(mid) > r ) then
        u=mid
     else
        l=mid
     endif
  end do

  ri=l
  rw= (r-grid%radial(l))/(grid%radial(u)-grid%radial(l))
  rw=1D0-rw
#ifdef XGC_DEBUG11
  if(rw<0D0 .or. rw >1D0) print *, 'error in search_rindex', r,ri,rw
#endif

end subroutine search_rindex

subroutine search_tindex(grid,t,ri,ti,tw)
  use grid_class
  use sml_module, only: sml_2pi
  implicit none
  type(grid_type) :: grid
  real (kind=8), intent(in):: t
  integer, intent(in) :: ri
  real (kind=8), intent(out):: tw
  integer, intent(out):: ti
  integer :: l,u ! lower and upper bound
  integer :: mid
  real (kind=8) :: tl,ul ! detla theta from lower bound theta
  real (kind=8) :: base
  integer , external :: angle_compare

  !bineary search
  l=grid%itheta0(ri)
  base=grid%theta(l)

  u=grid%itheta0(ri)+grid%ntheta(ri)-1
  ! check the angle is in btween ntheta and 1
  if(angle_compare(grid%theta(u),t,base)==-1 ) then
     ti=u
     u=l
     l=ti
  else

     do while ( u - l > 1 )
        mid=(l+u)/2
        if(angle_compare(grid%theta(mid),t,base)==1) then
           u=mid
        else
           l=mid
        endif
     enddo
     ti=l
  endif

  tl=t-grid%theta(l)
  if(tl<0D0) tl=tl+sml_2pi
  ul=grid%theta(u)-grid%theta(l)
  if(ul<0D0) ul=ul+sml_2pi
  tw= 1D0 - tl/ ul

#ifdef XGC_DEBUG11
  if(tw<0D0 .or. tw >1D0) print *, 'error in search_tindex', t,ri,ti,tw
#endif

  
end subroutine search_tindex



#endif

subroutine sheath_calculation(grid,psn,sp,iptl,type,itrout,pout)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  integer, intent(in) :: type, iptl
  integer :: itrout
  integer :: i,l, wall,itr
  real (kind=8) :: pout(3)
  real (kind=8) :: p(3),psave(3)
  integer, parameter :: wall_num=100
  real (kind=8), parameter :: minus_val=-1D50
  real (kind=8) :: rho,b,en_para, x(2), time_now
  real (kind=8) , external :: b_interpol
  logical, parameter :: USE_SEARCH_TR2 = .true. 
  real (kind=8) :: new_phase(sp%nphase)

  ! find nearest wall point

  new_phase=sp%phase0(:,iptl)
  x = new_phase(1:2)
        
  ! find position of previous time step
  if (USE_SEARCH_TR2) then
     call search_tr2(grid,x,itr,p)
  else
     call search_tr(grid,x,itr,p)
  endif
  psave=p

  do i=1, 3
     l=maxloc(p,1)
     wall = grid%nd(l,itr)

     if(grid%p(wall)==wall_num) then
        exit
     else
        p(l)=minus_val
        if(i==3) then
           if(sml_mype==0) print *, 'No nearest wall', itr,new_phase(1),new_phase(2)
           call remove_particle(sp,iptl,-1)
           return
        endif
     endif
  enddo
  
  !check potential and energy
  
  rho=new_phase(4)
  b=b_interpol(x(1),x(2))
  
  en_para=ptl_c2_2m(sp%type) * (rho*b)**2
  
  if(en_para < psn%pot0(wall)) then
     ! reflection
     new_phase(4)=-new_phase(4)
     
     ! reset derivatives
     if(sml_push_mode/=1) then
        time_now=sml_time
        call reset_derivs_without_e(new_phase,sp,iptl,time_now)
     endif

     sp%phase(:,iptl)=new_phase
  else
     ! remove particle
     call remove_particle(sp,iptl,-1)
     return
  endif  

  if(sp%gid(iptl)>0) then
     itrout=itr
     pout=psave
  endif

end subroutine sheath_calculation
