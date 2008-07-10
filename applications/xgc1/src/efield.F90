subroutine efield(grid,psn,sp,i,fld,time)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use efld_module
  use eq_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  integer, intent(in) :: i
  type(species_type) :: sp
  real (kind=8) :: time


  if(efld_mode<=0 .or. efld_cutoff==1 .and.&
       ( fld%psi <sml_inpsi + 0.01*eq_x_psi .or. &
       fld%psi > sml_outpsi - 0.01*eq_x_psi .and. sml_bounce==2 ) ) then  ! change boundary value to global variable later
     fld%Er=0D0
     fld%Ez=0D0
     fld%Ephi=0D0
     fld%Epot=0D0
     return
  endif
  !Modeled E


  !Simple 1D E


  !Gyrokinetic E
#if !defined(CIRCULAR_SPECIAL)
  call efield_gk(grid,psn,sp,i,fld)
#else
  call efield_cc(grid,psn,sp,i,fld)
#endif


end subroutine efield


subroutine efield_gk(grid,psn,sp,i,fld) !require fld%(r,z,phi,br,bz,bphi)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  type(species_type) :: sp
  integer, intent(in) :: i
  integer :: ip,nodes(3),larmor,iphi,iphi_frac
  real (kind=8) :: E_para,pot,E(3),phi_weight(0:1),B!,b_unit(3)
!  real (kind=8), pointer :: E_perp0(:,:)!, pot0(:)
  real (kind=8) :: weight_tmp

  E_para=0D0
  pot=0D0
  E=0D0


  !get unit b-vector
  B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
!  b_unit(1)=fld%br/B
!  b_unit(2)=fld%bz/B
!  b_unit(3)=fld%bphi/B

  !get phi index
  iphi= FLOOR(fld%phi/grid%delta_phi)    - grid%iphi_offset
  if(iphi<0) then 
     print *, 'Warning : iphi exceeded angle range-', iphi,fld%phi, grid%iphi_offset,grid%nphi,sp%type,i,sp%gid(i)
     call err_count
     iphi=0
  else if(iphi>=grid%nphi) then
     print *, 'Warning : iphi exceeded angle range+', iphi,fld%phi, grid%iphi_offset,grid%nphi,sp%type,i,sp%gid(i)
     call err_count
     iphi=grid%nphi-1
  endif

  phi_weight(1)= (fld%phi/grid%delta_phi)   - iphi - grid%iphi_offset
  phi_weight(0)=1D0 - phi_weight(1)

  !get E
  if(sp%type==1) then !Ion ---------------------------------------------------------
     if(sml_turb_efield==1) then
        ! parallel field
        do iphi_frac=1,2
           do larmor=1, sml_nlarmor
              if(sp%tr_save(larmor,iphi_frac,i)>0) then
                 nodes=grid%nd( :, sp%tr_save(larmor,iphi_frac,i) )  
#if defined(LINEAR_E_PERP)
                 !parallel field
#if !defined(LINEAR_E_PARA)                 
                 do ip=1, 3
                    E_para=E_para + phi_weight(iphi_frac-1)*sp%p_save(ip,larmor,iphi_frac,i)* &
                         psn%E_para(nodes(ip),iphi+iphi_frac-1)
                 enddo
#else
                 if(iphi_frac==2) then
                    do ip=1,3
                       E_para=E_para + sp%p_save(ip,larmor,iphi_frac,i)*psn%E_para(nodes(ip),iphi+iphi_frac-1)
                    enddo
                 endif
#endif
                 
                 ! potential                 
                 do ip=1,3
                    pot=pot+ phi_weight(iphi_frac-1)*sp%p_save(ip,larmor,iphi_frac,i) * &
                         ( psn%pot0(nodes(ip))+psn%dpot(nodes(ip),iphi+iphi_frac-1) )
                 enddo

                 ! perp field
                 E(1:2)=E(1:2)+ phi_weight(iphi_frac-1)* &
                      psn%E_perp_tr(:,sp%tr_save(larmor,iphi_frac,i),iphi+iphi_frac-1)             
#else
                 do ip=1, 3
                    weight_tmp = phi_weight(iphi_frac-1)*sp%p_save(ip,larmor,iphi_frac,i)
                    
                    E_para = E_para + weight_tmp * &
                         psn%E_para(nodes(ip),iphi+iphi_frac-1)
                    
                    pot    = pot    + weight_tmp * &
                         ( psn%pot0(nodes(ip))+psn%dpot(nodes(ip),iphi+iphi_frac-1) )
                    
                    E(1:2) = E(1:2) + weight_tmp * &
                         psn%E_perp_node(:,nodes(ip),iphi+iphi_frac-1)
                   
!                    E(1:2) = E(1:2) + phi_weight(iphi_frac-1)* &
!                         psn%E_perp(:,sp%tr_save(larmor,iphi_frac,i),iphi+iphi_frac-1)          
                 enddo
#endif

              else
                 print *, 'E-field ion invaild tr', i,sp%tr_save(larmor,iphi_frac,i)
              endif                           
           enddo

        enddo
        
        if(sml_bfollow==1) then
           E(3)=(E_para*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi
        else
           E(3)=E_para*sign(1D0,fld%Bphi)
        endif
     else !*********** no turbulence field *********************
        ! perp field only
        do larmor=1, sml_nlarmor
           if(sp%tr_save(larmor,1,i)>0) then
              nodes=grid%nd( :, sp%tr_save(larmor,1,i) )  
#if defined(LINEAR_E_PERP)
              E(1:2)=E(1:2)+ psn%E_perp0_tr(:,sp%tr_save(larmor,1,i))
              
              do ip=1,3
                 pot=pot+ sp%p_save(ip,larmor,1,i) * psn%pot0(nodes(ip))
              enddo
#else
              do ip=1,3
                 pot=pot+ sp%p_save(ip,larmor,1,i) * psn%pot0(nodes(ip))
                 E(1:2)=E(1:2) + sp%p_save(ip,larmor,1,i) * psn%E_perp0_node(:,nodes(ip))
              enddo
#endif
           endif
        enddo
     endif
     
     E=E/real(sml_nlarmor)
     pot=pot/real(sml_nlarmor)
  else  !Electron --------------------------------------------------------------
     ! parallel field
!     pot0=>psn%pot0        !dummy variable. not needed now
!     E_perp0=>psn%E_perp0  !dummy variable. not needed now
     
     if(sp%tr_save(1,1,i)>0) then
        nodes=grid%nd(:,sp%tr_save(1,1,i))
        do ip=1,3
           pot=pot + sp%p_save(ip,1,1,i) *  psn%pot0(nodes(ip))
           !              E_para=E_para + phi_weight(j)*ptl_e_p_save(ip,i) * &
           !                   psn%E_para(nodes(ip),iphi+j)
#if !defined(LINEAR_E_PERP)
           E(1:2)=E(1:2) + psn%E_perp0_node(:,nodes(ip))
#endif           
        enddo
        
        !        E(:)=E_para*b_unit(:)
        
        !perp field
#if defined(LINEAR_E_PERP)
        E(1:2)=E(1:2) + psn%E_perp0_tr(:,sp%tr_save(1,1,i))
#endif      
     else
        print *, 'E-field electron invalid tr', i,sp%tr_save(1,1,i)
     endif
  endif
  
  fld%Er=E(1)
  fld%Ez=E(2)
  fld%Ephi=E(3)
  fld%Epot=pot
  
end subroutine efield_gk

subroutine efield_init(grid,psn)
  use grid_class
  use psn_class
  use efld_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  
  if(efld_mode/=1 .and. efld_mode/=-2) return 
  
  if(efld_mode==1) then
     call efield_load(grid,psn)
  endif
  
  if(efld_mode==-2) then
     call efield_model_load(grid,psn)
     efld_mode=1
  endif
  
  call get_potential_grad(grid,psn)

end subroutine efield_init

subroutine efield_load(grid,psn)
  use grid_class
  use psn_class
  use efld_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  
  print *, 'Not implimented yet - efield_load'
  stop
end subroutine efield_load

subroutine efield_model_load(grid,psn)
  use grid_class
  use psn_class
  use efld_module
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,j
  real (kind=8) :: psi, theta, r,z,dr,dz,l,cost,sint,mmode,nmode,phi

  mmode=14
  nmode=-7
  ! set potential value with a function of f(psi,real-theta) = f1(psi)*f2(theta)
  do i=1, grid%nnode
     r=grid%x(1,i)
     z=grid%x(2,i)
     dr=r-eq_axis_r
     dz=z
     l=sqrt(dr**2 + dz**2)
     theta=acos(dr/l)
     if(dz < 0D0) then
        theta= sml_2pi - theta
     endif
     cost=cos(mmode*theta)
     sint=sin(mmode*theta)
     psi=grid%psi(i)
     psn%pot0(i)=0D0 
     do j=-1, grid%nphi+1
        phi=real(j + grid%iphi_offset)*grid%delta_phi
        psn%dpot(i,j) = 2000D0*psi*sint*cos(nmode*phi)
     enddo
  enddo

end subroutine efield_model_load

#ifdef CIRCULAR_SPECIAL

subroutine efield_cc(grid,psn,sp,i,fld)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  type(species_type) :: sp
  integer, intent(in) :: i
  integer :: larmor, iphi, iphi_frac
  real (kind=8) :: E_para, pot, E(3), phi_weight(0:1),B
  real (kind=8) :: weight_tmp
  real (kind=8) :: rw(2),tw(2)
  real (kind=8) :: cost,sint, x(2),radial, theta
  integer :: ri,ti,rindex,tindex,node

#ifdef CIRCULAR_SPECIAL2
  call efield_cc2(grid,psn,sp,i,fld)
  return
#endif

  E_para=0D0
  pot=0D0
  E=0D0



    !get unit b-vector
  B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)

  !get phi index
  iphi= FLOOR(fld%phi/grid%delta_phi)    - grid%iphi_offset

  ! below block is only for debugging
  if(iphi<0) then 
     print *, 'Warning : iphi exceeded angle range-', iphi,fld%phi, grid%iphi_offset,grid%nphi,sp%type,i,sp%gid(i)
     call err_count
     iphi=0
  else if(iphi>=grid%nphi) then
     print *, 'Warning : iphi exceeded angle range+', iphi,fld%phi, grid%iphi_offset,grid%nphi,sp%type,i,sp%gid(i)
     call err_count
     iphi=grid%nphi-1
  endif

  phi_weight(1)= (fld%phi/grid%delta_phi)   - iphi - grid%iphi_offset
  phi_weight(0)=1D0 - phi_weight(1)

  !get E
  if(sp%type==1) then !Ion ---------------------------------------------------------
!     if(sml_turb_efield==1) then
        ! parallel field
        do iphi_frac=1,2
           do larmor=1, sml_nlarmor

              ri=sp%rindex(larmor,iphi_frac,i)
              rw(1)=sp%rweight(larmor,iphi_frac,i)
              rw(2)=1D0-rw(1)              
              do rindex=1,2
                 ti=sp%tindex(rindex,larmor,iphi_frac,i)
                 tw(1)=sp%tweight(rindex,larmor,iphi_frac,i)
                 tw(2)=1D0-tw(1)
                 do tindex=1,2
                    call find_node(grid,ri+rindex-1,ti+tindex-1,node)

                    weight_tmp = phi_weight(iphi_frac-1)*rw(rindex)*tw(tindex)                    
                    E_para = E_para + weight_tmp * &
                         psn%E_para(node,iphi+iphi_frac-1)
                    
                    pot    = pot    + weight_tmp * &
                         ( psn%pot0(node)+psn%dpot(node,iphi+iphi_frac-1) )
                    
                    E(1:2) = E(1:2) + weight_tmp * &
                         psn%E_perp_node(:,node,iphi+iphi_frac-1)
                 enddo
              enddo
           enddo
        enddo

        
        E=E/real(sml_nlarmor)
        pot=pot/real(sml_nlarmor)

        ! coordinate transform from C.C. to R.Z
        x(1)=fld%r
        x(2)=fld%z
        call get_r_theta(x,radial,theta)
        cost=cos(theta)
        sint=sin(theta)
        fld%Er= E(1)*cost - E(2)*sint
        fld%Ez= E(1)*sint + E(2)*cost
         
        ! get toroidal E-field from E_para
        if(sml_bfollow==1) then
           E(3)=(E_para*B- fld%Er*fld%Br - fld%Ez*fld%Bz)/fld%Bphi
        else
           E(3)=E_para*sign(1D0,fld%Bphi)
        endif
!     endif !*********** no turbulence field *********************

     
  endif
  
  ! coordinate transform from RZ to C.C.
  fld%Ephi=E(3)
  fld%EPot=Pot 
end subroutine efield_cc

subroutine get_potential_grad_cc(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  use smooth_module
  use perf_monitor
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,iphi,j,k
  integer :: iphi_p, iphi_m
  real (kind=8):: dpot_p, dpot_m
  integer :: nodes_p(3),nodes_m(3),nd,sgn
  real (kind=8) :: pot_r(2,0:grid%nphi), rdt, dtheta, tw(2), dr
  integer :: ri, ti, tindex, ntmp, nd_plus, nd_minus
  call monitor_start (PHASE5_)
  

  do i=1, grid%npsi, grid%npsi-1
     do j=1, grid%ntheta(i)
        nd=grid%itheta0(i)+j-1
        psn%E_perp_node(:,nd,:)=0D0
     enddo
  enddo

  do i=1+1, grid%npsi-1
     do j=1, grid%ntheta(i)
        nd=grid%itheta0(i)+j-1
       ! CC grid should be sorted grid (not reordered)
        
        ! radial
        dr=grid%radial(i+1)-grid%radial(i-1)
        pot_r=0D0
        do k=-1,1,2 ! r minus dr/2, r plus dr/2
           ri=i+k
           ! find ti and tw
           call search_tindex(grid,grid%theta(nd),ri,ti,tw(1))
           tw(2)=1D0-tw(1)
           do tindex=1,2
              call find_node(grid,ri,ti+tindex-1,ntmp)
              pot_r((k+3)/2,:)=pot_r((k+3)/2,:)+tw(tindex)*(psn%pot0(ntmp)+psn%dpot(ntmp,0:grid%nphi))
           enddo           
        enddo
        psn%E_perp_node(1,nd,:)=-(pot_r(2,:)-pot_r(1,:))/dr

        ! theta
        nd_plus=grid%itheta0(i)+mod(j,grid%ntheta(i))
        nd_minus=grid%itheta0(i)+mod(j-2+grid%ntheta(i),grid%ntheta(i))
        dtheta=modulo(grid%theta(nd_plus)-grid%theta(nd_minus),sml_2pi)
        rdt=grid%radial(i)*dtheta
        psn%E_perp_node(2,nd,:)=-( psn%pot0(nd_plus)-psn%pot0(nd_minus) +&
             psn%dpot(nd_plus,0:grid%nphi)-psn%dpot(nd_minus,0:grid%nphi) )/rdt
     enddo
  enddo

  if(sml_bt_sign>0) then
     sgn=1
  else
     sgn=-1
  endif

  ! store grad_para (phi0 + phi) in E_para
  if(sml_turb_efield==1) then
     do i=1, grid%nnode
        do iphi=0, grid%nphi
           iphi_p=iphi+sgn
           iphi_m=iphi-sgn
           nodes_p=grid%nd(:,psn%bfollow_tr(1,i))
           nodes_m=grid%nd(:,psn%bfollow_tr(2,i))
           dpot_p=0D0
           dpot_m=0D0
           do nd=1,3
              dpot_p=dpot_p + psn%dpot(nodes_p(nd),iphi_p)*psn%bfollow_p(nd,1,i)
              dpot_m=dpot_m + psn%dpot(nodes_m(nd),iphi_m)*psn%bfollow_p(nd,2,i)
           enddo
           psn%E_para(i,iphi) = -(dpot_p-dpot_m)*psn%bfollow_1_dx(i)    !(psn%pot0(i)
        enddo
     enddo
  else
     psn%E_para=0D0
  endif
  call monitor_stop (PHASE5_)


#ifdef CIRCULAR_SPECIAL2
  call prepare_table_cc2(grid,psn)
#endif  

end subroutine get_potential_grad_cc

#endif

#ifdef CIRCULAR_SPECIAL2
subroutine efield_cc2(grid,psn,sp,i,fld)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  type(species_type) :: sp
  integer, intent(in) :: i
  integer :: larmor, iphi, iphi_frac
  real (kind=8) :: E_para, pot, E(3), phi_weight(0:1),B
  real (kind=8) :: weight_tmp
  real (kind=8) :: rw(2),tw(2)
  real (kind=8) :: cost,sint, x(2),radial, theta
  integer :: ri,ti,rindex,tindex,node
  real (kind=8) :: r_orig, t_orig, pot1, e1(2), e_para1
  integer :: node1, node2
  
  E_para=0D0
  pot=0D0
  E=0D0


  !get unit b-vector
  B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)

  !get phi index
  iphi= FLOOR(fld%phi/grid%delta_phi)    - grid%iphi_offset

  ! below block is only for debugging
  if(iphi<0) then 
     print *, 'Warning : iphi exceeded angle range-', iphi,fld%phi, grid%iphi_offset,grid%nphi,sp%type,i,sp%gid(i)
     call err_count
     iphi=0
  else if(iphi>=grid%nphi) then
     print *, 'Warning : iphi exceeded angle range+', iphi,fld%phi, grid%iphi_offset,grid%nphi,sp%type,i,sp%gid(i)
     call err_count
     iphi=grid%nphi-1
  endif

  phi_weight(1)= (fld%phi/grid%delta_phi)   - iphi - grid%iphi_offset
  phi_weight(0)=1D0 - phi_weight(1)


  !get E
  if(sp%type==1) then !Ion ---------------------------------------------------------
!     if(sml_turb_efield==1) then
        ! parallel field
        do iphi_frac=1,2
           do larmor=1, sml_nlarmor

              ri=sp%rindex(larmor,iphi_frac,i)
              rw(1)=sp%rweight(larmor,iphi_frac,i)
              rw(2)=1D0-rw(1) 

              ! obtain original r***
              r_orig=grid%radial(ri)*rw(1)+grid%radial(ri+1)*rw(2)

              !              do rindex=1,2
              ti=sp%tindex(1,larmor,iphi_frac,i)
              tw(1)=sp%tweight(1,larmor,iphi_frac,i)
              tw(2)=1D0-tw(1)


              call find_node(grid,ri,ti,node1)
              call find_node(grid,ri,ti+1,node2)
              if(grid%theta(node2) > grid%theta(node1) ) then
                 t_orig=tw(1)*grid%theta(node1) + tw(2)*grid%theta(node2)
              else
                 t_orig=tw(1)*( grid%theta(node1) - sml_2pi ) + tw(2)*grid%theta(node2)
                 if(t_orig<0D0) t_orig=t_orig+sml_2pi
              endif

              !get potential and E-field from r_/t_original value

              call get_e_pot_cc2(grid,psn,r_orig,t_orig,iphi+iphi_frac-1,pot1,e1(1),e1(2),e_para1)
              
              weight_tmp= phi_weight(iphi_frac-1)/real(sml_nlarmor)
              
              E_para = E_para + weight_tmp *  E_para1
              
              pot    = pot    + weight_tmp *  pot1
              
              E(1:2) = E(1:2) + weight_tmp *  E1(1:2)
                    
           enddo
        enddo
        
        
        E=E/real(sml_nlarmor)
        pot=pot/real(sml_nlarmor)

        ! coordinate transform from C.C. to R.Z
        x(1)=fld%r
        x(2)=fld%z
        call get_r_theta(x,radial,theta)
        cost=cos(theta)
        sint=sin(theta)
        fld%Er= E(1)*cost - E(2)*sint
        fld%Ez= E(1)*sint + E(2)*cost
         
        ! get toroidal E-field from E_para
        if(sml_bfollow==1) then
           E(3)=(E_para*B- fld%Er*fld%Br - fld%Ez*fld%Bz)/fld%Bphi
        else
           E(3)=E_para*sign(1D0,fld%Bphi)
        endif
!     endif !*********** no turbulence field *********************

     
  endif
  
  ! coordinate transform from RZ to C.C.
  fld%Ephi=E(3)
  fld%EPot=Pot 


end subroutine efield_cc2

subroutine get_e_pot_cc2(grid,psn,r,t,iphi,pot,er,et,ep)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn 
  integer :: iphi
  real (kind=8) :: r, t, pot, er, et, ep
  integer :: ri,ti, tip1
  real (kind=8) :: rw, tw, o, a, b, c

  
  call search_rindex(grid,r,ri,rw)
  ti=t/psn%dt_cc2+1
  tw=1D0 - (t/psn%dt_cc2+1-ti)
  tip1=ti+1
  if(tip1>psn%nt_cc2) tip1=tip1-psn%nt_cc2
  

  ! Potential and er, et
  o=psn%pot_cc2(ri  ,ti  ,iphi)
  a=psn%pot_cc2(ri+1,ti  ,iphi)-o
  b=psn%pot_cc2(ri  ,tip1,iphi)-o
  c=psn%pot_cc2(ri+1,tip1,iphi)-o

  pot=a*rw + b*tw + (c-a-b)*rw*tw + o
  er = -(a + (c-a-b)*tw)/(grid%radial(ri+1)-grid%radial(ri))
  et = -(b + (c-a-b)*rw)/(r*psn%dt_cc2)

  ! e_para
  o=psn%epara_cc2(ri  ,ti,iphi)
  a=psn%epara_cc2(ri+1,ti,iphi)-o
  b=psn%epara_cc2(ri  ,tip1,iphi)-o
  c=psn%epara_cc2(ri+1,tip1,iphi)-o

  ep=a*rw + b*tw + (c-a-b)*rw*tw + o

end subroutine get_e_pot_cc2

subroutine prepare_table_cc2(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn 
  integer :: ri, ti, tindex, node, ti2
  real (kind=8) :: r,t,tw(2)
  psn%pot_cc2=0D0
  psn%epara_cc2=0D0

  call check_point(1000)
  do ri=1, grid%npsi
     do ti=1, psn%nt_cc2
        r=grid%radial(ri)
        t= real(ti-1)*psn%dt_cc2
        call search_tindex(grid,t,ri,ti2,tw(1))
        tw(2)=1D0-tw(1)
        
        do tindex=1,2
           call find_node(grid,ri,ti2+tindex-1,node)
           psn%pot_cc2(ri,ti,0:grid%nphi)=psn%pot_cc2(ri,ti,0:grid%nphi)+ &
                tw(tindex)*( psn%pot0(node)+psn%dpot(node,0:grid%nphi) )
           psn%epara_cc2(ri,ti,0:grid%nphi)=psn%epara_cc2(ri,ti,0:grid%nphi)+&
                tw(tindex)* psn%E_para(node,0:grid%nphi)
        enddo        
     enddo
  enddo
end subroutine prepare_table_cc2


subroutine init_cc2_parameters(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn 

  print *, 'init_cc2_parameters'
  psn%nt_cc2=grid%ntheta(grid%npsi-3)
  psn%dt_cc2=sml_2pi/real(psn%nt_cc2)
  allocate(psn%pot_cc2(grid%npsi,psn%nt_cc2,0:grid%nphi),&
       psn%epara_cc2(grid%npsi,psn%nt_cc2,0:grid%nphi))
!  print *, 'end of init_cc2_parameters'

end subroutine init_cc2_parameters
#endif
