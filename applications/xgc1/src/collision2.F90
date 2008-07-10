! momentum conserving collision
! based on Plasma Phys. Control. Fusion, W X Wang, N nakjima, M Okamoto 41, 1091 (1999)
subroutine conserving_collision(ptl,flag)
  !1. store psi and B value for each ptls
  !2. store the sum of weight, energy, v_|| for each flux surface
  !2.5 get temparature and density as a function of psi
  !3. mode 2 collision 
  !4. store erf(y) value & coefficient a,b, and d
  !5. P 
  !6. P1
  !7. P2
  use sml_module
  use col_module
  use ptl_module
  use eq_module,only : eq_x_z, eq_x_psi
  use perf_monitor
  implicit none
  type(ptl_type):: ptl
  integer, intent(in) :: flag
  integer :: i,j,k,l,wi
  integer, allocatable :: org_j(:),org_l(:)
  real (kind=8) :: r,z,phi,b,rho,mu,psi,weight,theta,factor
  real (kind=8) :: two_therm_en, v2th,y, sqrty, dphidy, phiy!,deltapsi
  real (kind=8),allocatable :: org_b(:),coeff_a(:), coeff_b(:), coeff_d(:)
  real (kind=8) ,dimension(col_2_m,col_2_mtheta) :: new_weight_sum, new_v2_avg,new_vp_avg, &
       delta_weight, delta_v2, delta_vp,dum
  real (kind=8) ,external :: psi_interpol, b_interpol,derf,tempi_ev
!  real (kind=8) ,external :: initden

  save new_weight_sum, new_v2_avg, new_vp_avg, delta_weight, delta_v2, delta_vp, dum, wi
  save org_b,org_j,org_l, coeff_a, coeff_b, coeff_d
 
  if(flag==0) then
     allocate(org_j(ptl%ion%maxnum), org_l(ptl%ion%maxnum),org_b(ptl%ion%maxnum),&
          coeff_a(ptl%ion%maxnum), coeff_b(ptl%ion%maxnum), coeff_d(ptl%ion%maxnum))
     return
  endif

!  deltapsi=(col_2_pout-col_2_pin)/real(col_2_m)

  ! initialization
  col_2_org_weight_sum=1D-90 ! assign very small number to prevent divied by zero
  col_2_org_v2_avg=0D0
  col_2_org_vp_avg=0D0
  org_j(1:ptl%ion%num)=-1
  org_l(1:ptl%ion%num)=-1

  ! angle calculation
  !  if(sml_angle_stored==0) then
  call store_angle(ptl)
  !  endif

  ! 1, 2, 3. 4
  ! Get original weight, para_momentum, 2*kinetic_energy/m (=v^2)
  do i=1, ptl%ion%num
     if(ptl%ion%gid(i) > 0) then
        ! store psi and B ------------------------------------------------------
        weight=ptl%ion%phase(8,i)
        if(sml_deltaf==1) then
           weight=weight*ptl%ion%phase(6,i)
        endif
        r=ptl%ion%phase(1,i)
        z=ptl%ion%phase(2,i)
        phi=ptl%ion%phase(3,i)
        b=b_interpol(r,z,phi)
        rho=ptl%ion%phase(4,i)
        mu=ptl%ion%phase(5,i)
        psi=ptl%ion%phase(9,i)
        org_b(i)=b

        ! store the sum of weight, energy, and v_|| -----------------------------
        if(col_2_pin<psi .AND. psi<col_2_pout .AND. z>eq_x_z ) then
           j=int((psi-col_2_pin)/col_2_dp)+1
           j=min(col_2_m,max(1,j))
           org_j(i)=j
           theta=ptl%ion%angle(i)
           l=int(theta/col_2_dtheta) +1
           l=min(col_2_mtheta,max(1,l))
           org_l(i)=l
           col_2_org_weight_sum(j,l)=col_2_org_weight_sum(j,l)+weight
           col_2_org_v2_avg(j,l)=col_2_org_v2_avg(j,l)+ weight *  2D0*( mu*B+ptl_c2_2m(1)*(rho*B)**2 )/ptl_mass(1) ! v^2 sum
           col_2_org_vp_avg(j,l)=col_2_org_vp_avg(j,l)+ weight * ptl_c_m(1)*(rho*B) ! v_|| sum
        endif
     endif
  enddo

  ! parallel summation 

  call monitor_start (CONS_COL_RED_)
  call my_mpi_allreduce(col_2_org_weight_sum,dum,col_2_m*col_2_mtheta)
  if(sml_deltaf==1) then
!        do i=1,col_2_m
!           psi=col_2_dp*(real(i)-0.5) + col_2_pin
!           col_2_org_weight_sum(i,sp)=col_2_vol(i)*initden(psi,sp)
!        enddo
!        print *, '#1', col_2_org_weight_sum(col_2_m/2,1), dum(col_2_m/2,1) 
!        col_2_org_weight_sum(:,sp)=col_2_org_weight_sum(:,sp) + dum(:,sp)
     col_2_org_weight_sum(:,:)=col_2_w0(:,:) + dum(:,:)
  else
     col_2_org_weight_sum(:,:)=dum(:,:)
  endif
  
  call my_mpi_allreduce(col_2_org_v2_avg,dum,col_2_m*col_2_mtheta)
  col_2_org_v2_avg=dum
  call my_mpi_allreduce(col_2_org_vp_avg,dum,col_2_m*col_2_mtheta)
  col_2_org_vp_avg=dum
  ! summation to avg
  col_2_org_v2_avg=col_2_org_v2_avg/col_2_org_weight_sum ! 2 kin_en average
  col_2_org_vp_avg=col_2_org_vp_avg/col_2_org_weight_sum ! v_|| average
  call monitor_stop (CONS_COL_RED_)

  ! Apply MC collision
  do i=1, ptl%ion%num
     if(ptl%ion%gid(i)>0) then        
        psi=ptl%ion%phase(9,i)
        b=org_b(i)
        rho=ptl%ion%phase(4,i)
        mu=ptl%ion%phase(5,i)

        ! mode 2 collision -----------------------------------------------------
        !nan test
!        if(psi > 1D0 .or. psi < 2D0 )  then
!        else
!           print *,'col before scatr call', psi,b, rho,mu, i, sml_time/sml_dt
!        endif
        call scatr(psi,b,rho,mu,4,1,ptl%ion) !collision to spcies 1 only , 3: energy and pitch angle
        ptl%ion%phase(4,i)=rho
        ptl%ion%phase(5,i)=mu
        
        ! get erf(y) value and store coefficient---------------------------------------------------
        two_therm_en=2d0 * (  tempi_ev(psi)*sml_ev2j  )
        v2th=two_therm_en/ptl_mass(1)
        y=2d0*(mu*B+ptl_c2_2m(1)*(rho*B)**2)/two_therm_en ! v^2/vth^2 = En/Eth
        sqrty=sqrt(y)

        dphidy=2D0/sml_sqrtpi*exp(-y)*sqrty
        phiy=derf(sqrty)  - dphidy
        coeff_a(i)=-3D0 *sml_sqrt2pi *phiy *ptl_c_m(1)*rho*B/ (sqrty**3 *v2th) 
        coeff_b(i)=-sml_sqrt2pi * (phiy-dphidy)/(sqrty *v2th)
        coeff_d(i)=1.5D0 * sml_sqrt2pi * (phiy-dphidy)/y -1D0
     endif
  enddo
  



  !energy change compensation
  !if(k==1 .and. col_en_col_on) then

  !energy summation

  !   do i=1, ptl_num
  !      j=org_j(i)
  !      l=org_l(i)
  !      if(org_j(i)/=-1) then
  !         factor=col_2_org_v2_avg(j,l)/new_v2_avg(j,l)
  !         ptl_phase(i,5)=factor*ptl_phase(i,5)
  !         ptl_phase(i,4)=sqrt(factor)*ptl_phase(i,4)
  !      endif
  !   enddo
  !   new_v2_avg=col_2_org_v2_avg
  !endif

  ! get erf(y) value and store coefficient



  !------------------
  
  ! P, P1, P2, ..., Pn
  do k=1, 3
     
     new_weight_sum=1D-90 ! a very small number for avoiding dividing zero by zero
     new_v2_avg=0D0
     new_vp_avg=0D0
     ! get new sum
     do i=1, ptl%ion%num
        if(org_j(i)/=-1) then           
           ! sum of new weight, v^2 and v_||
           j=org_j(i)
           l=org_l(i)
           mu=ptl%ion%phase(5,i)
           rho=ptl%ion%phase(4,i)
           b=org_b(i)
           weight=ptl%ion%phase(8,i)
           if(sml_deltaf==1) then
              weight=weight*ptl%ion%phase(6,i)
           endif
           new_weight_sum(j,l)=new_weight_sum(j,l)+ weight
           new_v2_avg(j,l)=new_v2_avg(j,l)+ weight *  2D0*( mu*B+ptl_c2_2m(1)*(rho*B)**2 )/ptl_mass(1) ! 2 kin_en total sum
           new_vp_avg(j,l)=new_vp_avg(j,l)+ weight * ptl_c_m(1)*(rho*B) ! v_|| total sum
        endif
     enddo
     ! parallel sumation of weight, v2_avg, vp_avg

     call monitor_start (CONS_COL_RED2_)
     call my_mpi_allreduce(new_weight_sum,dum,col_2_m*col_2_mtheta)
     if(sml_deltaf==1) then
!           do i=1,col_2_m
!              psi=col_2_dp*(real(i)-0.5) + col_2_pin
!              new_weight_sum(i,sp)=col_2_vol(i)*initden(psi,sp)
!           enddo
!           print *, '#2', new_weight_sum(col_2_m/2,1) ,dum(col_2_m/2,1)
!           new_weight_sum(:,sp)=new_weight_sum(:,sp) + dum(:,sp)
        new_weight_sum(:,:)=col_2_w0(:,:) + dum(:,:)
     else
        new_weight_sum(:,:)=dum(:,:)
     endif
     
!     new_weight_sum=dum
     call my_mpi_allreduce(new_v2_avg,dum,col_2_m*col_2_mtheta)
     new_v2_avg=dum
     call my_mpi_allreduce(new_vp_avg,dum,col_2_m*col_2_mtheta)
     new_vp_avg=dum
     call monitor_stop (CONS_COL_RED2_)

     ! sumation to average
     new_v2_avg=new_v2_avg/new_weight_sum
     new_vp_avg=new_vp_avg/new_weight_sum
     
     delta_weight = new_weight_sum/col_2_org_weight_sum -1D0
!     delta_weight = new_weight_sum - col_2_org_weight_sum
     delta_v2 = new_v2_avg-col_2_org_v2_avg
     delta_vp = new_vp_avg-col_2_org_vp_avg


     !get momentum and energy from other species
!     if(ptl_sp_num > 1) then
!        do sp=2, ptl_sp_num
!           delta_v2(:,:,1)=delta_v2(:,:,1)-delta_v2(:,:,sp) * col_2_org_weight_sum(:,:,sp)/col_2_org_weight_sum(:,:,1)
!           delta_vp(:,:,1)=delta_vp(:,:,1)-delta_vp(:,:,sp) * col_2_org_weight_sum(:,:,sp)/col_2_org_weight_sum(:,:,1)
!        enddo
!     endif

     ! debug
!     if (sml_mype==0) then
!        print *, '*******************************************'
!        print *, delta_weight(col_2_m/2,1) , maxval(abs(delta_weight))
!        print *, '-------------------------------------------'
!        print *, delta_v2(col_2_m/2,1)/col_2_org_v2_avg(col_2_m/2,1), maxval(abs(delta_v2))
!        print *, '-------------------------------------------'
!        print *, delta_vp(col_2_m/2,1)/(col_2_org_vp_avg(col_2_m/2,1)+1d-99),maxval(abs(delta_vp))
!     endif

     
     ! Change weight
!     do sp=1, ptl_sp_num
!        if(mod(int(sml_deltaf)/int(2**(sp-1)),2)==1) then
!           wi(sp)=6
!        else
!           wi(sp)=8
!        endif
!     enddo

     do i=1, ptl%ion%num
        if(org_j(i)/=-1 ) then
           j=org_j(i)
           l=org_l(i)
           if(sml_deltaf==1) then
              ptl%ion%phase(6,i)= ptl%ion%phase(6,i)   &
                   +coeff_d(i)*delta_weight(j,l)  & 
                   +coeff_a(i)*delta_vp(j,l)       &
                   +coeff_b(i)*delta_v2(j,l)      
           else
              ptl%ion%phase(8,i) = ptl%ion%phase(8,i)* (1.D0       &
                   +coeff_d(i)*delta_weight(j,l)  & 
                   +coeff_a(i)*delta_vp(j,l)       &
                   +coeff_b(i)*delta_v2(j,l)   )
           endif
           !2003/10/16 weight >0 
           if(sml_deltaf/=1) then
              ptl%ion%phase(8,i)=max(ptl%ion%phase(8,i),1D-5)
           endif
        endif        
     enddo

  enddo


  if (sml_mype==0) then
     if( maxval(abs(delta_weight)) > 0.02 ) then
        print *, 'col3 error: weight change exceeded : time=', sml_time/sml_dt,maxval(abs(delta_weight))
        print *, 'location : ', maxloc(abs(delta_weight)),':',col_2_m,col_2_mtheta
        call err_count
     endif

!     if( maxval(delta_v2/col_2_org_v2_avg) > 0.03 ) then
!        print *, 'col3 error: energy change exceeded : time=', sml_time/sml_dt,maxval(delta_v2/col_2_org_v2_avg),col_2_org_v2_avg
!     endif
        
     !        if( max(delta_weight) > 0.2 ) then
     !           print *, 'col3 error: weight change exceeded : time=', sml_time/sml_dt,delta_weight
     !        endif

     ! store change of weight and energy due to collision routine
     ! this is only for diagnosis. no need for other processor
     
     col_2_dw_sum=col_2_dw_sum + delta_weight
     col_2_dvp_sum=col_2_dvp_sum + delta_vp
     col_2_dv2_sum=col_2_dv2_sum + delta_v2
     
  endif





  !dubug
!  if(sml_mype==0) then
!     print *, '   '
!     print *,'=========================================='           
!  endif
  
  
  ! electron + impurity collision
  do i=1, ptl%ion%num     
     if(ptl%ion%gid(i) >0) then
        ! store psi and B ------------------------------------------------------
        b=org_b(i)
        rho=ptl%ion%phase(4,i)
        mu=ptl%ion%phase(5,i)
        psi=ptl%ion%phase(9,i)
        
        ! mode 2 collision -----------------------------------------------------
        if(col_en_col_on==1) then
           call scatr(psi,b,rho,mu,4,2,ptl%ion) ! self, energy only
        endif
       
        call scatr(psi,b,rho,mu,1+2*col_imp_on,1+col_en_col_on*2,ptl%ion) ! 3 impurity and electron coll only 2002/09/18 
        ptl%ion%phase(4,i)=rho
        ptl%ion%phase(5,i)=mu
        
     endif
  enddo


end subroutine conserving_collision
        

! store real angle to ptl_angle and set sml_angle_stored=1  
subroutine store_angle(ptl)
  use ptl_module
  use sml_module
  use perf_monitor
  use eq_module
  implicit none
  type(ptl_type) :: ptl
  real (kind=8) :: r,z,theta
  integer :: i
  real (kind=8) :: t0,t1
  
  call monitor_start(ANGLE_)
  do i=1, ptl%ion%num
     if(ptl%ion%gid(i)>0) then
        r=ptl%ion%phase(1,i)
        z=ptl%ion%phase(2,i)

        theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
        if(z < 0D0) then
           theta= sml_2pi - theta
        endif
        ptl%ion%angle(i)=theta
     endif
  enddo
  sml_angle_stored=1
  call monitor_stop(ANGLE_)
  
end subroutine store_angle
