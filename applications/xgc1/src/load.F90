!! Initialize Phase variabls of ions.
!!>For single particle simulation, call load_single
!!<For V-space hole simulation, call load_special3
subroutine load(grid,psn,ptl)
  use ptl_module
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(ptl_type) :: ptl
  type(psn_type) :: psn
  type(species_type),pointer :: sp 
  integer :: i, sp_type
  integer (kind=8) :: ioff, idum1,idum2
  real (kind=8) :: r,z,phi,rho,mu,b,psi
  real (kind=8) :: psi_interpol, b_interpol,init_den,ni,init_tempi_ev,ti_ev
  real (kind=8) :: ni_c, ti_ev_c, psi_c, en_ev,br,bz,bphi,ti_ev_vertual
  target ptl
  interface
     subroutine uniform_space_dist(grid,ptl)
       use grid_class
       use ptl_module
       type(grid_type), intent(in)::grid
       type(ptl_type),target  :: ptl
     end subroutine uniform_space_dist
     subroutine maxwell_v_dist(grid,ptl)
       use grid_class
       use ptl_module
       type(grid_type), intent(in)::grid
       type(ptl_type),target  :: ptl
     end subroutine maxwell_v_dist
     subroutine restart_write(ptl)
       use ptl_module
       implicit none
       type(ptl_type),target :: ptl  
     end subroutine restart_write
     subroutine canonical_maxwell_v_dist(grid,psn,ptl)
       use grid_class
       use ptl_module
       use psn_class
       type(grid_type) :: grid
       type(ptl_type) :: ptl
       type(psn_type) :: psn
     end subroutine canonical_maxwell_v_dist
  end interface

  ! loading for single particle simulation
  if(sml_special==1) then
     call load_single(ptl)
     return
  endif
  
  if(sml_special==3) then
     call load_special3(ptl)
     return
  endif

  ! restart from previous runs -----------------------------
  if(sml_restart /= 0)then
#if defined(ADIOS) 
     if(sml_mype==0)write(*,*)'restart read start'
     call restart_read(ptl)
     if(sml_mype==0)write(*,*)'restart read end'
     return
#else
     stop '** error ** Compile with adios or binpack support for restart!!'
#endif
  else
     sml_istep=0
     sml_time=0.D0
  endif
  !---------------------------------------------------------

  ! First, set particle position (r,z,phi) uniformly in simulation region
  call uniform_space_dist(grid,ptl)
  if(sml_deltaf/=1 .and. sml_electron_on==1 .and. ptl%ion%num==ptl%elec%num) then
     !debug
     ptl%elec%phase(1:3,1:ptl%ion%num) = ptl%ion%phase(1:3,1:ptl%ion%num)
     ptl%elec%phase(6:9,1:ptl%ion%num) = ptl%ion%phase(6:9,1:ptl%ion%num)
  endif

  if(sml_canonical_maxwell/=1) then
     call maxwell_v_dist(grid,ptl)
  else
     call canonical_maxwell_v_dist(grid,psn,ptl)
  endif

  do sp_type=1, 1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif

!     ioff = sml_mype*sp%num
     idum1=sml_mype
     idum2=sp%num 
     ioff = idum1 * idum2
     do i=1,sp%num
        sp%gid(i)=ioff+i
     enddo

!     sp%maxgid=sml_totalpe*sp%num
     idum1=sml_totalpe
     idum2=sp%num
     sp%maxgid=idum1*idum2
  enddo
end subroutine load

subroutine maxwell_v_dist(grid,ptl)
  use grid_class
  use sml_module
  use ptl_module
  use eq_module
  implicit none
  type(grid_type), intent(in)::grid
  type(ptl_type),target  :: ptl
  type(species_type),pointer :: sp
  real (kind=8) :: r,z,phi,b,psi,n,t_ev,t_ev_vertual,rho,mu,n_c,t_ev_c,bphi,br,bz,en_ev,psi_c
  integer :: i,sp_type
  real (kind=8),external :: init_den_wz,init_tempi_ev_wz,b_interpol,tempe_ev_wz,tempe_ev,init_den,init_tempi_ev,f0_den2
  real (kind=8) :: marker_den  
!  integer, pointer :: pnum
!  real (kind=8), pointer :: phase(:,:)
  real (kind=8), external :: ranx
#ifdef GAM_TEST
  real (kind=8), external :: gam_weight
#endif



  do sp_type=1, 1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
        marker_den=sml_marker_den
     else
        sp=>ptl%elec
        marker_den=sml_marker_den*real(ptl%elec%num)/real(ptl%ion%num)
     endif

     do i=1, sp%num
        !retrive position variable from the uniform_space_dist
        r=sp%phase(1,i)
        z=sp%phase(2,i)
        phi=sp%phase(3,i)
        b=b_interpol(r,z,phi)
        psi=sp%phase(9,i)
     
     
        !get density and temp
        if(sp_type==1) then
           n=init_den_wz(psi,z)  !get density
           t_ev=init_tempi_ev_wz(psi,z)
        else
           n=init_den_wz(psi,z)
           t_ev=tempe_ev_wz(psi,z)
        endif

        !setting particle weight          
        sp%phase(8,i)=n/marker_den !sml_marker_den=(ptl number per unit volume) -- background profile, 2002/05/30
        ! set very small or zero value for delta-f weight
        if(sml_deltaf==1) then
           if(sml_deltaf_f0_mode==1) then
              sp%phase(6:7,i)= (n-f0_den2(grid,r,z,psi))/n  
           else
              sp%phase(6:7,i)=0D0
           endif
           if(abs(sp%phase(6,i)) <= sml_initial_deltaf_noise ) then
              sp%phase(6:7,i) = sml_initial_deltaf_noise*2D0*(ranx()-0.5D0)*sp%phase(1,i)/eq_axis_r ! ballooning mode shape
              if(sml_mode_initial_n/=0) then
                 sp%phase(6:7,i)=abs(sp%phase(6,i))*cos(sp%phase(3,i)*sml_mode_initial_n)
              endif
           endif
#ifdef GAM_TEST
           sp%phase(6:7,i)=gam_weight(r,z,psi)
#endif
        else
           sp%phase(6:7,i)=0D0
        endif
     
        ! Velocity distribution
        t_ev_vertual=t_ev*sml_marker_temp_factor  !set marker particle temperature
        call maxwell_dist(t_ev_vertual,b,rho,mu,ptl_mass(sp_type),ptl_charge(sp_type))   !get rho and mu from Monte-Carlo method
        sp%phase(4,i)=rho                          ! set rho and mu to phase variable
        sp%phase(5,i)=mu 

        if(sml_marker_temp_factor/=1D0) then
           sp%phase(4,i)=rho                          ! set rho and mu to phase variable
           sp%phase(5,i)=mu
           en_ev=(ptl_c2_2m(sp_type)*(rho*b)**2 + mu*B)*sml_j2ev
           sp%phase(8,i)=sp%phase(8,i)*exp( - en_ev*(1D0/t_ev - 1D0/t_ev_vertual) ) * (t_ev_vertual/t_ev)**1.5D0
        endif
        
        ! canonical maxwellian distribution
        ! Adjust weight.
        ! Not good for even distribution - Algorithm need to be changed later
        if(sml_canonical_maxwell==1) then
           call bvec_interpol(r,z,phi,br,bz,bphi)
           b=sqrt(br**2+bz**2+bphi**2)

           if(sml_minusb/=1) then
              psi_c= psi + rho* r * bphi
           else
              psi_c= psi - rho* r * bphi
           endif

           !psi_c=min(sml_outpsi,max(sml_inpsi,psi_c))
           if(sp_type==1) then
              t_ev_c=init_tempi_ev_wz(psi_c,z)
              n_c=init_den_wz(psi_c,z)
           else
              t_ev_c=tempe_ev_wz(psi_c,z)
              n_c=init_den_wz(psi_c,z)
           endif

           en_ev=(ptl_c2_2m(sp_type)*(rho*b)**2 + mu*B)*sml_j2ev
           sp%phase(8,i) = sp%phase(8,i)* exp( - en_ev*(1D0/t_ev_c - 1D0/t_ev) ) * n_c/n * (t_ev/t_ev_c)**1.5D0
           !print *, exp( - en_ev*(1D0/ti_ev_c - 1D0/ti_ev) ) * ni_c/ni *(ti_ev/ti_ev_c),ni_c/ni
        endif

! Store initial (ion) particle weight for collisional dissipation diagnostic
        if (sp_type==1) then
           sp%weight0(i) = (1D0 + sp%phase(6,i)) * sp%phase(8,i)
        endif

     enddo

  enddo
     
end subroutine maxwell_v_dist

subroutine canonical_maxwell_v_dist(grid,psn,ptl)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use eq_module
  implicit none
  type(grid_type), intent(in)::grid
  type(psn_type), intent(in) :: psn
  type(ptl_type),target  :: ptl
  type(species_type),pointer :: sp
  integer :: sp_type, i,j,k
  real (kind=8) :: r, z, phi, psi, br, bz, bphi, b, psi_c,mass,charge, pot, pot_c
  real (kind=8) :: t_max, v_max, dv, v_phi, k_phi, t_c, t_ev_c,weight, marker_den
  real (kind=8) :: rnd, vr, vz, v_perp1,v_perp2,v_para,  zmax, zdum,rho,mu
  real (kind=8),external :: canonical_maxwell_psi, canonical_maxwell_f1, inverse_ftn,ranx, init_tempi_ev_wz, tempe_ev_wz
  real (kind=8), parameter :: v_max_factor=5D0
  integer, parameter :: nv=1000
  real (kind=8) :: a(nv),f1(0:2)
  
  

  do sp_type=1, 1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
        marker_den=sml_marker_den
        t_max=max(eq_tempi_ev_edge,eq_tempi_ev_out)*sml_ev2j
    else
        sp=>ptl%elec
        marker_den=sml_marker_den*real(ptl%elec%num)/real(ptl%ion%num)
        t_max=max(eq_tempe_ev_edge,eq_tempe_ev_out)*sml_ev2j
     endif

     mass=ptl_mass(sp_type)
     charge=ptl_charge(sp_type)
     v_max=v_max_factor*sqrt(2D0/mass*t_max)
     dv=2D0*v_max/real(nv-1)

     
     do i=1, sp%num
        !retrive position variable from the uniform_space_dist
        r=sp%phase(1,i)
        z=sp%phase(2,i)
        phi=sp%phase(3,i)
        psi=sp%phase(9,i)
        call bvec_interpol(r,z,phi,br,bz,bphi)
        b=sqrt(br**2+bz**2+bphi**2)
        
        pot=0D0  ! set initial potential 
        call get_add_pot0_psi(psi,pot,grid,psn)

        !2.  Evaluate f1(v_phi) and simpson integral
        
        a(1)=0D0

        v_phi=-v_max
        psi_c= canonical_maxwell_psi(r,psi,v_phi,mass,charge)
        pot_c=0D0  ! set canonical potential
        call get_add_pot0_psi(psi_c,pot_c,grid,psn)
        f1(0)= canonical_maxwell_f1(r,z,psi_c,v_phi,mass,pot,pot_c,sp_type)
        

        do j=2,nv
           do k=1,2
              v_phi=-v_max + real(j-2)*dv + real(k)*0.5D0*dv
              
              psi_c= canonical_maxwell_psi(r,psi,v_phi,mass,charge)
              pot_c=0D0  ! set canonical potential
              f1(k)= canonical_maxwell_f1(r,z,psi_c,v_phi,mass,pot,pot_c,sp_type)
           enddo
           
           a(j)=a(j-1) +  ( f1(0)+4D0*f1(1)+f1(2) ) * dv / 6D0  
           f1(0)=f1(2)
        enddo
           
        !4   get weight
        weight=a(nv)

        !5   normalize  -> a should be a increasing ftn between 0 and 1
        a=a/weight

        !6  random number and its inverse
        rnd=ranx()
        v_phi=inverse_ftn(rnd,a,nv,-v_max,dv)
        

        !7.  generate  v_perp1/
        psi_c= canonical_maxwell_psi(r,psi,v_phi,mass,charge) ! get real psi_c
        if(sp_type==1) then
           t_ev_c=init_tempi_ev_wz(psi_c,z)
        else
           t_ev_c=tempe_ev_wz(psi_c,z)
        endif
        t_c=t_ev_c*sml_ev2j

        
        zmax=1.d0 - dexp(-7.d0)
        zdum=zmax*ranx()
        v_perp1= sqrt(-2.d0/mass*dlog(1.d0-zdum)*t_c)   ! for unit temperature

        rnd=ranx()
        vr= v_perp1*cos(sml_2pi*rnd)
        vz= v_perp1*sin(sml_2pi*rnd)
        
        v_para=(vr*br + vz*bz + v_phi*bphi) /b
        v_perp2= vr**2 + vz**2 + v_phi**2 - v_para**2
        

        mu= 0.5D0*mass*v_perp2/b
        rho=mass*v_para/charge/b

        !---------- finallize

        sp%phase(4,i)=rho                          ! set rho and mu to phase variable
        sp%phase(5,i)=mu 
        sp%phase(8,i)=weight/marker_den
        ! Store initial (ion) particle weight for collisional dissipation diagnostic
        if (sp_type==1) then
           sp%weight0(i) = (1D0 + sp%phase(6,i)) * sp%phase(8,i)
        endif

     end do
  end do
end subroutine canonical_maxwell_v_dist

! get canonical angular momentum with 
real (kind=8) function canonical_maxwell_psi(r,psi,v_phi,mass,charge)
  use sml_module
  implicit none
  real (kind=8) :: r,psi,v_phi,mass,charge

  if(sml_minusb/=1) then
     canonical_maxwell_psi= psi + mass/charge* r * v_phi
  else
     canonical_maxwell_psi= psi - mass/charge* r * v_phi
  endif

end function canonical_maxwell_psi

real (kind=8) function canonical_maxwell_f1(r,z,psi_c,v_phi,mass,pot,pot_c,sp_type)
  use sml_module
  implicit none
  real (kind=8) :: r,z,psi_c,v_phi,mass,pot,pot_c
  integer :: sp_type
  real (kind=8) :: k_phi, t_c, t_ev_c
  real (kind=8) , external :: init_den_wz, tempe_ev_wz, init_tempi_ev_wz
  
  
  k_phi=0.5D0*mass*v_phi**2
  
  if(sp_type==1) then
     t_ev_c=init_tempi_ev_wz(psi_c,z)
  else
     t_ev_c=tempe_ev_wz(psi_c,z)
  endif
  t_c=t_ev_c*sml_ev2j
  
  
  canonical_maxwell_f1=init_den_wz(psi_c,z)*exp( (pot_c-pot)/t_ev_c )/sqrt(sml_2pi*t_c/mass)*exp(-k_phi/t_c)

end function canonical_maxwell_f1

real (kind=8) function inverse_ftn(r,a,nv,v0,dv)
  integer :: nv
  real (kind=8) :: a(nv), v0,dv
  real (kind=8) :: ind, r,ranx
  integer :: mid, left, right
  
  ! binary search
  left=1
  right=nv
  
  
  do while ( right > left+1)
     mid=(right+left)/2
     if(a(mid)<r ) then
        left=mid
     else
        right=mid
     endif
  enddo
  

  ! linear interpolation
  ind= left + (r-a(left))/(a(left+1)-a(left))
  inverse_ftn= v0 + (ind-1D0)*dv

end function inverse_ftn


!!  Distribute particle uniformly in space. (inpsi < psi < outpsi)
!!> The r,z,phi,psi of phase variable are set
!!< 
subroutine uniform_space_dist(grid,ptl)
  use grid_class
  use ptl_module
  use sml_module
  use eq_module
  implicit none
  type(ptl_type) :: ptl
  type(grid_type), intent(in) :: grid
  type(species_type),pointer :: sp
  integer :: sp_type
  integer :: valid  !! number of particle that generated inside the region inpsi<psi<outpsi
  integer :: total  !! number of total particle generated 
  real (kind=8) :: rdim, zdim   !! R,Z dimension of whole simulation region
  real (kind=8) :: roffset, zoffset !! offset (start position) of simulation region
  real (kind=8) :: r,z,psi,phi
  real (kind=8) , external :: ranx, psi_interpol
  real (kind=8) :: x(2),p(3)
  integer :: itr
  target ptl
  interface
       subroutine get_volume(ptl)
       use ptl_module
       type(ptl_type),target  :: ptl
     end subroutine get_volume
  end interface

! simulation boundary is imposed 2001/01/24
  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z
  
  do sp_type=1, 1+sml_electron_on
     valid=0
     total=0

     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif

     ! generate particle until # of valid particles become ptl_num 
     do while(valid<sp%num)
        !generate r,z in simulation region
        r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) ) !2002/05/27
        z=zdim*ranx() + zoffset
        psi=psi_interpol(r,z,0,0)
        total=total+1

        !check psi validity
        if(sml_inpsi < psi .AND. psi < sml_outpsi) then
           ! check inside wall validity
           x(1)=r
           x(2)=z
           call search_tr(grid,x,itr,p)
           if(itr>0)then
              valid=valid+1
              phi=(ranx()*grid%nphi+ grid%iphi_offset)*grid%delta_phi   ! set toroidal angle
              if(sml_2pi < phi) phi=sml_2pi   ! just for in case
              if(phi<0D0) phi=0D0             ! just for in case
              !set phase variables to global storage
              sp%phase(1,valid)=r
              sp%phase(2,valid)=z
              sp%phase(3,valid)=phi
              sp%phase(9,valid)=psi_interpol(r,z,0,0)
           end if
        endif
     enddo

  enddo

  call shift_check(grid,ptl%ion)
  if(sml_electron_on==1) call shift_check(grid,ptl%elec)
  
  ! Get space volume using Monte-Carlo method
  call get_volume(ptl)

end subroutine uniform_space_dist



subroutine maxwell_dist(ti_ev,b,rho,mu,mass,charge)
  use sml_module
  implicit none
  real (kind=8) , intent(in) :: b,ti_ev,mass,charge
  real (kind=8) , intent(out) :: rho,mu
  real (kind=8) :: theta,r,z,pitch,phi,zdum,zmax,v,ti
  real (kind=8), external :: ranx


  ti=ti_ev*sml_ev2j
  
  ! perp velocity
  zmax=1.d0 - dexp(-7.d0)
  zdum=zmax*ranx()
!  v= sqrt(-2.d0*dlog(1.d0-zdum)*ti/mass)
  v= sqrt(-2.d0/mass*dlog(1.d0-zdum)*ti)
  mu=0.5D0*mass*v**2/b   ! mu is indep to mass and charge for given temp
  ! parallel velocity
  zdum=zmax*ranx()
!  v= sqrt(-2.d0*dlog(1.d0-zdum)*ti/mass)
  v= sqrt(-2.d0/mass*dlog(1.d0-zdum)*ti)

  v= v*cos(sml_pi*ranx())
  rho=v/b*mass/charge !*mass/charge
  
endsubroutine


!!************************************************************
!!>Calculating the voulme of shells
!!
!! first created 2002/05/27
!! adopted from get_marker_den -> renaming and modifying
!! r dependant loading
!!<***********************************************************
subroutine get_volume(ptl)
  use diag_module
  use ptl_module
  use sml_module
  use eq_module
  use col_module
  use efld_module
  use perf_monitor
  implicit none
  type(ptl_type),target :: ptl
  integer ::  valid    !! number of valid particle 
  integer ::  total    !! number of total particle generated
  integer :: j,ierror, sp, i,l
  real (kind=8) :: rdim, zdim  !! R,Z dimension of whole simulation region
  real (kind=8) :: roffset, zoffset !! offset (start position) of simulation region
  real (kind=8) :: tvolume, vvolume !! Total volume and valid volume
  real (kind=8) :: r,z,psi,theta,phi,dum1,dum2
  real (kind=8) ::  dum_f(diag_flow_npsi),dum_a(diag_avg_npsi),&
       dum_c(col_2_m,col_2_mtheta),dum_e(efld_npsi),&
       dum_fv(diag_f_npsi),dum_vb(col_vb_m)  !! dummy variables for reduce operation
  real (kind=8) , external :: ranx, psi_interpol, init_den
  real (kind=8) :: aa,bb   !! Linear interpolation weight

! simulation boundary is imposed 2001/01/24
  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z

  valid=0
  total=0
  
  ! initialize volume
  diag_flow_vol=1D-50
  diag_avg_vol=1D-50
  col_2_vol=1D-50
  efld_vol=1D-50
  diag_f_dvol=1D-50
  col_vb_vol=1D-50
  
!  print *, 'get_vol'
  do while(valid<sml_monte_num)
     r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) ) 
     z=zdim*ranx() + zoffset
     psi=psi_interpol(r,z,0,0)
     total=total+1
     
     ! marker density
     if(sml_inpsi < psi .AND. psi < sml_outpsi .and. z>eq_x_z) then
        valid=valid+1
     endif
     
     !flow
     if(diag_flow_pin < psi .AND. psi < diag_flow_pout .and. z>eq_x_z) then
        j=int((psi-diag_flow_pin)/diag_flow_dp)+1
        j=min(diag_flow_npsi,max(1,j))
        diag_flow_vol(j)=diag_flow_vol(j)+1D0
     endif
     
     ! diag_avg
     if(diag_avg_pin < psi .AND. psi < diag_avg_pout .and. z>eq_x_z) then
        j=int((psi-diag_avg_pin)/diag_avg_dp)+1
        j=min(diag_avg_npsi,max(1,j))
        bb=(psi-diag_avg_pin)/diag_avg_dp + 1 -j
        aa= 1D0-bb
        diag_avg_vol(j)=diag_avg_vol(j)+1D0*aa
        diag_avg_vol(j+1)=diag_avg_vol(j+1)+1D0*bb
     endif
     
     ! mode 2 (conserving) collision
     if(col_2_pin < psi .AND. psi < col_2_pout .and. z>eq_x_z) then

        theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
        if(z < 0D0) then
           theta= sml_2pi - theta
        endif
                
        j=int((psi-col_2_pin)/col_2_dp)+1
        j=min(col_2_m,max(1,j))
        l=int(theta/col_2_dtheta)+1  
        l=min(col_2_mtheta,max(1,l))
        col_2_vol(j,l)=col_2_vol(j,l)+1D0
     endif     

     ! Background update
     if(col_vb_pin < psi .and. psi < col_vb_pout .and. z>eq_x_z) then
        j=int((psi-col_vb_pin)/col_vb_dp) +1
        j=min(col_vb_m-1,max(1,j))
        bb=(psi-col_vb_pin)/col_vb_dp + 1 - j
        aa= 1D0- bb
        col_vb_vol(j)=col_vb_vol(j) + 1D0 * aa
        col_vb_vol(j+1)=col_vb_vol(j+1) + 1D0 * bb
     endif

     !efield
     if(efld_pin < psi .AND. psi < efld_pout .and. z>eq_x_z) then
        j=int((psi-efld_pin)/efld_dp)+1
        j=min(efld_npsi-1,max(1,j))
        bb=(psi-efld_pin)/efld_dp + 1 - j
        aa= 1D0 - bb
        efld_vol(j)=efld_vol(j)+1D0*aa
        efld_vol(j+1)=efld_vol(j+1)+1D0*bb
     endif
     
     ! distribution function diagnosis
     if(diag_f_pin < psi .AND. psi < diag_f_pout .AND. z>eq_x_z &
          .and. r>1D0 .and. z< diag_f_slope1*(r-eq_axis_r) .and. z> diag_f_slope2*(r-eq_axis_r)) then
        j=int((psi-diag_f_pin)/diag_f_dpsi)+1
        j=min(diag_f_npsi,max(1,j))
        diag_f_dvol(j)=diag_f_dvol(j)+1D0
     endif

  enddo

  ! total sum
  dum1=real(total)
  
  call monitor_start (GET_VOL_RED_)
  call my_mpi_allreduce(dum1,dum2,1)
  call my_mpi_allreduce(diag_flow_vol,dum_f,diag_flow_npsi)
  diag_flow_vol=dum_f
  call my_mpi_allreduce(diag_avg_vol,dum_a,diag_avg_npsi)
  diag_avg_vol=dum_a
  call my_mpi_allreduce(col_2_vol, dum_c, col_2_m*col_2_mtheta)
  col_2_vol=dum_c
  call my_mpi_allreduce(efld_vol, dum_e, efld_npsi)
  efld_vol=dum_e
  call my_mpi_allreduce(diag_f_dvol,dum_fv,diag_f_npsi)
  diag_f_dvol=dum_fv
  call my_mpi_allreduce(col_vb_vol,dum_vb,col_vb_m)
  col_vb_vol=dum_vb
  call monitor_stop (GET_VOL_RED_)

  ! marker_den = density of simulation ptl = ptl_num*total_pe/(real volume * monte_num *totla_pe/sum of total )
  ! = ptl_num*sum of total /(real_volume*monte_num)
  sml_marker_den= real(ptl%ion%num)*dum2/ (  rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) * sml_monte_num)  ! normalized unit -  is normalization constant
  
  ! volume of each flux shell -- shell volume = (ptl number of each shell) / (sum of total) * volume
  !flow
  diag_flow_vol=diag_flow_vol/dum2*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) ) ! bug fix 2002/07/04
  !diag_avg
  diag_avg_vol=diag_avg_vol/dum2*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) ) ! bug fix 2002/07/04
  ! collision
  col_2_vol=col_2_vol/dum2*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) ) ! bug fix 2002/07/04
  !efield
  efld_vol=efld_vol/dum2*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) ) ! bug fix 2002/07/04
  !diag_f
  diag_f_dvol=diag_f_dvol/dum2*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) ) 
  !col_vb_vol
  col_vb_vol=col_vb_vol/dum2*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) )
  
  ! w0 calculation for local maxwell-----------
!  do sp=1, ptl_sp_num
!  print *, 'w0 calculation'
     do i=1,col_2_m
        psi=col_2_dp*(real(i)-0.5) + col_2_pin
        do l=1, col_2_mtheta
           col_2_w0(i,l)=col_2_vol(i,l)*init_den(psi)
        enddo
     enddo
!  enddo

end subroutine 

! for single particle simulation. Load a particle in a given position and velocity.
subroutine load_single(ptl)
  use ptl_module
  use sml_module
  implicit none
  type(ptl_type) :: ptl
  real (kind=8) :: rho,mu,en,pitch,b
  real (kind=8), external :: b_interpol,psi_interpol


  b=b_interpol(ptl_special_r,ptl_special_z,0D0)
  
  call ev_pitch_to_rho_mu2(ptl_special_en_ev,ptl_special_pitch,b,rho,mu,1)
  call rho_mu_to_ev_pitch2(rho,mu,b,en,pitch,1)
  print *, en,ptl_special_en_ev
  print *,rho,mu
  ptl%ion%phase(1,1)=ptl_special_r
  ptl%ion%phase(2,1)=ptl_special_z
  ptl%ion%phase(3,1)=ptl_special_phi
  ptl%ion%phase(4,1)=rho
  ptl%ion%phase(5,1)=mu
  ptl%ion%phase(6:7,1)=1D-10
  ptl%ion%phase(8,1)=1D0
  ptl%ion%phase(9,1)=psi_interpol(ptl_special_r,ptl_special_z,0,0)
  ptl%ion%gid(1)=1
  ptl%ion%maxgid=1
  ptl%ion%num=1

  if(sml_electron_on==1) then
     ptl%elec%phase(:,1)=ptl%ion%phase(:,1)
     ptl%elec%phase(4,1)=ptl%ion%phase(4,1)*sqrt(ptl_mass(2))/ptl_charge(2)
     ptl%elec%gid(1)=1
     ptl%elec%num=1
     ptl%elec%maxgid=1
  endif
end subroutine load_single

! for V-space hole calculation. sml_special==3
subroutine load_special3(ptl)
  use sml_module
  use ptl_module
  implicit none
  type(ptl_type) :: ptl
  real (kind=8) :: rho,mu,en,pitch,b
  real (kind=8), external :: b_interpol,psi_interpol
  integer :: i, ioff

  b=b_interpol(ptl_special_r,ptl_special_z,0D0)
  ptl%ion%phase(1,1:ptl%ion%num)=ptl_special_r
  ptl%ion%phase(2,1:ptl%ion%num)=ptl_special_z
  ptl%ion%phase(3,1:ptl%ion%num)=ptl_special_phi

  ioff = sml_mype*ptl%ion%num
  do i=1, ptl%ion%num
     ptl%ion%gid(i)=ioff+i
     en=ptl_special_en_ev*real(i/ptl_special_n)/real(ptl_num/ptl_special_n)
     pitch=-1D0 + 2* real(mod(i,ptl_special_n))/real(ptl_special_n-1)
     call ev_pitch_to_rho_mu2(en,pitch,b,rho,mu,1)
     ptl%ion%phase(4,i)=rho
     ptl%ion%phase(5,i)=mu
     ptl%ion%phase(6:7,i)=1D-10
     ptl%ion%phase(8,i)=1D0
     ptl%ion%phase(9,i)=psi_interpol(ptl_special_r,ptl_special_z,0,0)
  enddo
  ptl%ion%maxgid = sml_totalpe*ptl%ion%num

end subroutine load_special3

  
!! ******************************************************
!!>functions for energy/pitch <-> rho/mu conversion
!!<*****************************************************
subroutine rho_mu_to_ev_pitch2(rho,mu,b,ev,pitch,sp_type)
  use sml_module
  use ptl_module
  implicit none
  real (kind=8), intent(inout) :: rho,mu,b
  real (kind=8), intent(out) :: ev,pitch
  real (kind=8) :: enj,v_pal,v_pep
  integer :: sp_type

  if(mu<0) then
     print *, 'minus mu found :',rho,mu,b
     mu=0D0
  endif
  
  enj=(mu*b+ptl_c2_2m(sp_type)*(rho*b)**2)
  ev=enj*sml_j2ev
  v_pal=ptl_c_m(sp_type)*rho*b
  v_pep=SQRT(2.d0*mu*b/ptl_mass(sp_type))
  pitch=v_pal/SQRT(v_pal**2+v_pep**2)
  
end subroutine rho_mu_to_ev_pitch2

subroutine ev_pitch_to_rho_mu2(ev,pitch,b,rho,mu,sp_type)
  use sml_module
  use ptl_module
  implicit none
  real (kind=8), intent(in) :: ev,pitch,b
  real (kind=8), intent(out) :: rho,mu
  real (kind=8) :: en
  integer :: sp_type
  
  en=ev*sml_ev2j
  !       rho=pitch*DSQRT(2d0*en/mass)/b*mass/charge
  rho=pitch*SQRT(en/ptl_c2_2m(sp_type)) /B
  mu=max(0D0,(1d0-pitch**2)*en/b)
       
end subroutine ev_pitch_to_rho_mu2

! IMSL randum functions
function ranx()
#if !defined(IMSL)
  use Ecuyer_random
#endif
  implicit none
  real (kind=8) :: ranx
#if defined(IMSL)
  real (kind=8) , external :: drnunf ! IMSL function
  ranx=DRNUNF()
#else
  ranx=taus88()
#endif
end function ranx


real (kind=8) function init_den_wz(psi,z)
  use eq_module
  use sml_module
  implicit none
  real (kind=8), intent(in) :: psi,z
  real (kind=8), external :: init_den,ranx
  real (kind=8) :: tmp
  if(z>eq_x_z .or. psi > eq_x_psi) then
     tmp=psi
  else
     tmp=eq_x_psi
  endif
  init_den_wz=init_den(tmp)
end function init_den_wz

real (kind=8) function init_tempi_ev_wz(psi,z)
  use eq_module
  implicit none
  real (kind=8), intent(in) :: psi,z
  real (kind=8), external :: init_tempi_ev
  real (kind=8) :: tmp
  if(z>eq_x_z .or. psi > eq_x_psi) then
     tmp=psi
  else
     tmp=eq_x_psi
  endif
  init_tempi_ev_wz=init_tempi_ev(tmp)
  
end function init_tempi_ev_wz

real (kind=8) function tempe_ev_wz(psi,z)
  use eq_module
  implicit none
  real (kind=8), intent(in) :: psi,z
  real (kind=8), external :: tempe_ev
  real (kind=8) :: tmp
  if(z>eq_x_z .or. psi > eq_x_psi) then
     tmp=psi
  else
     tmp=eq_x_psi
  endif
  tempe_ev_wz=tempe_ev(tmp)

end function tempe_ev_wz

#ifdef GAM_TEST
real (kind=8) function gam_weight(r,z,psi)
  use sml_module
  use eq_module
  implicit none
  real (kind=8) :: r, z, psi
  real (kind=8) :: a,b,c,k, r_minor2

  a=300D0
  b=sml_inpsi
  k=0.137/3.6D-3
  r_minor2=(r-eq_axis_r)**2 + (z-eq_axis_z)**2
  gam_weight=sml_initial_deltaf_noise *sin( k * sqrt(r_minor2) )


!  gam_weight=sml_initial_deltaf_noise*sin(a*(psi-b))
!  gam_weight=sml_initial_deltaf_noise

end function gam_weight

#endif
