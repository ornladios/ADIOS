!*****************************************************
! neutral collision (charge exchange, ionization)
! first created 2002/09/26
! last modified
! This reoutine assumes initail maximum number and electron number is the same.
! shodul be generalized later.
!*****************************************************


subroutine neutral_col(istep,ptl)
  use ptl_module
  use neu_module
  use sml_module
  use eq_module, only : eq_x_psi
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(ptl_type) :: ptl
  type(species_type), pointer:: sp
  integer ,intent(in) :: istep
  integer :: i,k,l
  real (kind=8),external :: b_interpol, neutral_den,ranx,psi_interpol,te_neutral_col,neutral_temp
!  real (kind=8),external :: outptl
  real (kind=8) :: r,z,phi,rho,mu,b,n,f_cx,enev,pitch,Te,psi,t0
  real (kind=8) ,allocatable::  f_ion(:)
!  integer , allocatable ::pindex(:)
  integer :: dumi, old_ptl_num, num_ionized, newgid_start, ierror
  integer, allocatable :: global_ionized(:)
  real (kind=8) :: zmax, zdum, v, dum , dum1(3), dum2(3),&
       weight_sum_out,weight_sum_ionize,new_n0,weight_sum_inside
  real (kind=8) :: coll,agg,ekin
  real (kind=8) :: t00,t11
  target ptl

  ! neutral collision simulation 
  
  if(neu_col_mode==2 .AND. (mod(istep,neu_mode2_period)==0 .OR. istep==1)) then
     ! neutral diffusion calculaton
     call monitor_start (NEUTRAL_STEP_) 
     call neutral_step(ptl%ion%mass)
     call monitor_stop (NEUTRAL_STEP_)
  endif

  sp=>ptl%ion
  
  if(istep < neu_start_time) return ! no neutral effect before a couple profile evaluations

  ! algorithm 2
  !! Ionization ----------------------------------------------------------------
  if(mod(istep,neu_ion_period)==0) then
     ! print *,'ionize',istep,neu_ion_period
     ! ionization due to electron - neutral collision
     ! recycling gid==-gid ptl
     ! choose a ptl randomly
     
!     print *, 'memory_cleaning', sml_mype
!     call memory_cleaning  ! eliminate all gid(i)<0 particles and find the sum of weights
!     print *, 'memory_cleaning end', sml_mype

     if(neu_varying_mfp==1) then  ! Change mean free path 
        call change_neu_mfp(ptl) ! set new mean free path
     endif

     ! initialize
!     outptl_num=0
!     empty_num=0
     allocate(f_ion(sp%num))
     weight_sum_ionize=0D0 ! weight sum of ionize ptl (expectation value)
     f_ion=0D0     ! ionization probability
     weight_sum_inside=0D0
     neu_actual_ionized=0D0
     ! calculate ionization probability
     do i=1, sp%num
        if(sp%gid(i) >0) then
           ! weight sumation of ionize ptl (expectation value)
           r=sp%phase(1,i)
           z=sp%phase(2,i)
           psi=sp%phase(9,i)
           if(psi<neu_ionize2_psi) weight_sum_inside=weight_sum_inside+sp%phase(8,i)
           Te=Te_neutral_col(psi) ! eV 
           n=neutral_den(r,z,psi) ! normalized unit
           f_ion(i)= n * 0.8D-8*sqrt(Te)*(exp(-13.56/Te)/(1D0+0.01D0*Te))*1D-6  ! MKS !bugfix 1D-6
           !          print *, i,n,sqrt(Te), exp(-13.56/Te)
           f_ion(i)=f_ion(i)! s^-1 --> normalized unit  , unit conversion
           if(neu_ionize_mode/=2 .OR. psi < neu_ionize2_psi) then !exclude ionize_mode==2 .and. psi > neu_ionize2_psi
              weight_sum_ionize=weight_sum_ionize+f_ion(i)*sml_dt*neu_ion_period*sp%phase(8,i) 
           endif
        endif
     enddo
     ! MPI sumation
     call monitor_start (NEUTRAL_SUM_)
     if(neu_ionize_mode==1) then
        dum1(1)=neu_weight_sum_lost
        dum1(3)=dum1(1) !for diagnosis only
        neu_weight_sum_lost=0D0 !initilize for next use
     else  ! ionize mode 2
        dum1(1)=max(neu_old_inside_weight-weight_sum_inside,1D-50)
!        print *,neu_old_inside_weight, weight_sum_inside_sep !debug
        neu_old_inside_weight=weight_sum_inside
        dum1(3)=neu_weight_sum_lost !for diagnosis ---
        neu_weight_sum_lost=0D0
     endif
     dum1(2)=weight_sum_ionize
     call my_mpi_allreduce(dum1,dum2,3)
     weight_sum_out=dum2(1)  !if neu_ionize_mode /=1 weight_sum_out is 'difference of weight sum inside sep'
     weight_sum_ionize=dum2(2)
     call monitor_stop (NEUTRAL_SUM_)

     if(neu_adjust_n0==1) then
        ! Finding neutral density for given recycle rate
        new_n0=neu_n0*weight_sum_out/weight_sum_ionize*neu_recycle_rate
       if(istep>=neu_start_time .and. istep < neu_start_time +neu_ion_period) then
          new_n0=min(1.d15  , new_n0)
       endif
       !if(sml_mype==0) then
       !   write(62,*) sml_dt*istep, neu_n0/(NC_norm_r**3),new_n0/(NC_norm_r**3)&
       !        ,weight_sum_out/weight_sum_ionize*neu_recycle_rate,neu_mfp
       !endif
        
        ! set new ionization probability
        f_ion=f_ion*new_n0/neu_n0
        neu_n0=new_n0 + 1D-10 ! 1D-10 is for avoiding 'devide by zero'
     endif
     i=1 ! ptl index
     old_ptl_num=sp%num
     do while(i<=old_ptl_num )! gid=-gid ptl -> f_ion(i)=0
        if(f_ion(i)*sml_dt*neu_ion_period > ranx()) then
           if(f_ion(i)*sml_dt*neu_ion_period > 1.5D0 .AND. sml_mype==0 .AND. sml_mstep/=1) then
              print *, sml_mype, f_ion(i)*sml_dt*neu_ion_period , 'Too large dt in ionization'
           endif

           l=sp%num+1
           sp%num=sp%num+1
!              print *, ptl_num,ptl_num0, ptl_num_max

           r=sp%phase(1,i)
           z=sp%phase(2,i)
           psi=sp%phase(9,i)
           phi=sp%phase(3,i)
           rho=sp%phase(4,i)
           mu=sp%phase(5,i)
           ! set newly ionized ptl's position  & weight
           sp%phase(1,l)=r
           sp%phase(2,l)=z
           sp%phase(9,l)=psi
           sp%phase(3,l)=phi

           sp%phase(8,l)=sp%phase(8,i)
           if(neu_ionize_mode/=2 .or. psi < neu_ionize2_psi) then
              neu_actual_ionized=neu_actual_ionized+sp%phase(8,l)
           endif
           sp%phase(6:7,l)=0D0  ! delta-f scheme is not implimented for neutral collision

           ! ionize -- 2002/10/09 same energy ionization inside separatrix
           b=b_interpol(r,z,phi)
           if(psi<eq_x_psi .and. neu_col_mode==1) then
              call rho_mu_to_ev_pitch2(rho,mu,b,enev,pitch,1)
              pitch=2D0*ranx() - 1D0
              call ev_pitch_to_rho_mu2(enev,pitch,b,rho,mu,1)
           else  ! Maxwellian - outside of separatrix or mode2
              !perp velocity
              t0=neutral_temp(r,z,psi)
              zmax=1.d0 - dexp(-7.d0)
              zdum=zmax*ranx()
              v= dsqrt(-2.d0*dlog(1.d0-zdum)*t0)
              mu=0.5D0*v**2 /b
              ! parallel velocity
              zdum=zmax*ranx()
              v= dsqrt(-2.d0*dlog(1.d0-zdum)*t0)
              v= v*dcos(sml_pi*ranx())
              rho=v/b
           endif

           sp%phase(4,l)=rho
           sp%phase(5,l)=mu
           ! add electron with same energy and temp.
           ptl%elec%num=ptl%elec%num+1
           ptl%elec%phase(:,ptl%elec%num)=ptl%ion%phase(:,l)
           ptl%elec%phase(4,ptl%elec%num)=ptl%ion%phase(4,l)*sqrt(ptl_mass(2))/ptl_charge(2)
         
        endif
        i=i+1

     enddo
     if(neu_ionize_mode==2) then
        neu_old_inside_weight=neu_old_inside_weight+neu_actual_ionized
     endif
     if( sp%num> sp%maxnum) then
        print *, sml_mype , 'small memory space',real(i)/real(sp%num),neu_n0
     endif

     deallocate(f_ion)

     ! print out n0
     if(sml_mype==0 .and.sml_istep >neu_ion_period  ) then
        write(62,*) sml_time/sml_tran, neu_n0,dum2(3)*sml_e_charge/(sml_dt*neu_ion_period) 
     endif

! set global IDs for new particles and increment max global ID
     call monitor_start (NEUTRAL_GATHER_)
     num_ionized = sp%num - old_ptl_num
     if (sml_mype==0) then
        allocate(global_ionized(sml_totalpe))
     endif
! gather number of new particles on each processor
     call MPI_GATHER(num_ionized,1,MPI_INTEGER,global_ionized,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     if (sml_mype==0) then
        dumi = 0
        do i=1,sml_totalpe
           dumi = dumi + global_ionized(i)
        enddo
        ptl%ion%maxgid = ptl%ion%maxgid + dumi
! refill global array with newgid_start values and scatter
        global_ionized(sml_totalpe) = ptl%ion%maxgid - global_ionized(sml_totalpe)
        do i=sml_totalpe-1,1,-1
           global_ionized(i) = global_ionized(i+1) - global_ionized(i)
        enddo
     endif
     call MPI_SCATTER(global_ionized,1,MPI_INTEGER,newgid_start,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     if (sml_mype==0) then
        deallocate(global_ionized)
     endif
! assign global IDs to new particles
     do i=1,num_ionized
        sp%gid(old_ptl_num+i) = newgid_start + i
        ptl%elec%gid(old_ptl_num+i) = newgid_start + i
     enddo
     call monitor_stop (NEUTRAL_GATHER_)

  endif ! end of ionization

  ! Elastic Collisions -------------------

  if(mod(istep,neu_col_period)==0) then
     ! Pitch angle scattering by neutrals -> change mu and rho
     do i=1, sp%num
        if(sp%gid(i) > 0) then
           r=sp%phase(1,i)
           z=sp%phase(2,i)
           psi=sp%phase(9,i)
           rho=sp%phase(4,i)
           mu=sp%phase(5,i)
           b=b_interpol(r,z,phi)
           n=neutral_den(r,z,psi)
           call rho_mu_to_ev_pitch2(rho,mu,b,enev,pitch,1)
           ekin=enev*sml_ev2j
           !cccccccccccccccccccccccccccccccccccccc
           !ccc   pitch angle scattering for small value of ion-neutral
           !collision frequency, coll
           !ccccc Boozer and Kuo-Petravic Phys. Fluids 24, 851 (1981)
           !C============
           coll=1.55D-13*(n)*dsqrt(enev/1D3/sp%mass)
           !coll = col_col*(sml_en_order/(ekin))**1.5
           agg = ranx() - .5D0
           dum = 1.D0 - pitch**2
           dum = max(0.D0,dum)
           agg = dsign(1.D0,agg)*dsqrt(dum*coll*sml_dt*neu_col_period*0.5D0)
           pitch = pitch*(1.D0 - coll*sml_dt*neu_col_period*.5D0) + agg  !2002/09
           pitch = min( 0.9999999D0,pitch)
           pitch = max(-0.9999999D0,pitch)

           !cc   change rho ,rmu,  due to scattering

           rho = pitch*dsqrt(2.D0*ekin/ptl_mass(1))/b
           mu = ekin/b - ptl_c2_2m(1)*rho*rho*b
           sp%phase(4,i) = rho
           sp%phase(5,i) =mu
           !debug
           !if(mu<0) then
           !   print *, 'minus mu in neutral',ekin*sml_norme2ev,pitch,rho,mu
           !   stop
           !endif
        endif
     enddo

  endif

  ! Charge exchange----------------------------------------------------------------
  if(mod(istep,neu_cx_period)==0) then
     !  print *, 'cx',istep,neu_cx_period
     ! charge exchange. -> change ion mu, rho

     do i=1, sp%num
        if(sp%gid(i) > 0) then
           r=sp%phase(1,i)
           z=sp%phase(2,i)
           psi=sp%phase(9,i)
           rho=sp%phase(4,i)
           mu=sp%phase(5,i)
           b=b_interpol(r,z,phi)
           n=neutral_den(r,z,psi)
           call rho_mu_to_ev_pitch2(rho,mu,b,enev,pitch,1)
           ! probability of charge exchange per sec
           f_cx=( n * 1.1D-8* (enev)**0.3 /sqrt(sp%mass) *1D-6 )  
           ! 1.1 K_i(eV)^0.3 * 10^-8 / Sqrt(Mi(AU))  cm s^-1   -> 1D-6 m s^-1
           f_cx = f_cx  ! s^-1 --> normalized unit. 
           if(sml_dt*neu_cx_period*f_cx > ranx() ) then
              !Charge exchange -- give random change of pitch, Ti=Tn 2002/10/09 inside separatrix  -> same energy   
              ! Ti=/= Tn exactly, because of temp. dependance of f_cx
              if(psi<eq_x_psi .and. neu_col_mode==1) then
                 pitch=2D0*ranx() - 1D0
                 call ev_pitch_to_rho_mu2(enev,pitch,b,rho,mu,1)
              else  ! outside of seapratrix , Maxwellian dist
                 t0=neutral_temp(r,z,psi)
                 !perp velocity
                 zmax=1.d0 - dexp(-7.d0)
                 zdum=zmax*ranx()
                 v= dsqrt(-2.d0*dlog(1.d0-zdum)*t0)
                 mu=0.5D0*v**2 /b
                 ! parallel velocity
                 zdum=zmax*ranx()
                 v= dsqrt(-2.d0*dlog(1.d0-zdum)*t0)
                 v= v*dcos(sml_pi*ranx())
                 rho=v/b
              endif

              sp%phase(4,i)=rho
              sp%phase(5,i)=mu
           endif
        endif
     enddo
  endif

end subroutine neutral_col

real (kind=8) function neutral_den(r,z,psi)
  use neu_module
  use sml_module
  use eq_module
  implicit none
  real (kind=8) ,external :: neutral_den1
  real (kind=8), intent(in) :: r , z ,psi
  integer :: ipsi,itheta,itheta_p1
  real (kind=8) :: a1,a2,theta

  if(neu_col_mode==1) then
     neutral_den=neutral_den1(r,z)
  else     
     if(neu_grid_min_psi < psi .AND. psi < neu_grid_max_psi .and. z>eq_x_z) then
        theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
        if(z < 0D0) then
           theta= sml_2pi - theta
        endif
        ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1 
        ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
        a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi +1
        
        itheta= int(theta/neu_grid_dtheta) + 1 
        itheta= min(neu_grid_mtheta, max(1, itheta))
        a2= theta/neu_grid_dtheta - itheta +1
        itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1
        
        neutral_den=neu_grid_den(ipsi,itheta)*(1D0-a1)*(1D0-a2) + &
             neu_grid_den(ipsi+1,itheta)*a1*(1D0-a2) + &
             neu_grid_den(ipsi,itheta_p1)*(1D0-a1)*a2 + &
             neu_grid_den(ipsi+1,itheta_p1)*a1*a2 
        neutral_den=neutral_den*neu_n0  ! propotional to neu_n0
     elseif(neu_grid_min_psi >= psi ) then
        neutral_den=1D-10  ! very small number instead of zero
        
     else
        neutral_den=neutral_den1(r,z)
     endif
     neutral_den=max(neutral_den,1D-10)  !very small lower bound for safety
  endif

end function neutral_den
  

real (kind=8) function neutral_temp(r,z,psi)
  use neu_module
  use sml_module
  use eq_module
  implicit none
  real (kind=8),external :: init_tempi_ev
  real (kind=8),intent(in) :: r,z,psi
  integer :: ipsi,itheta,itheta_p1
  real (kind=8) :: a1,a2,theta

  if(neu_col_mode==1) then
     neutral_temp=init_tempi_ev(psi)*sml_ev2j*neu_temp_factor
  else ! mode 2
     if(neu_grid_min_psi < psi .AND. psi < neu_grid_max_psi .and. z>eq_x_z) then
        theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
        if(z < 0D0) then
           theta= sml_2pi - theta
        endif
        ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1 
        ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
        a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi +1
        
        itheta= int(theta/neu_grid_dtheta) + 1 
        itheta= min(neu_grid_mtheta, max(1, itheta))
        a2= theta/neu_grid_dtheta - itheta +1
        itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1
        
        neutral_temp=neu_grid_temp(ipsi,itheta)*(1D0-a1)*(1D0-a2) + &
             neu_grid_temp(ipsi+1,itheta)*a1*(1D0-a2) + &
             neu_grid_temp(ipsi,itheta_p1)*(1D0-a1)*a2 + &
             neu_grid_temp(ipsi+1,itheta_p1)*a1*a2 
!     elseif(neu_grid_min_psi >= psi ) then
!        neutral_temp=???   
     else
        neutral_temp=init_tempi_ev(psi)*sml_ev2j
     endif
  endif

end function neutral_temp

real (kind=8) function neutral_den1(r,z) !return nuetral density in normalized unit
  use neu_module
  use eq_module
  use sml_module, only : sml_2pi
  implicit none
  real (kind=8) , intent(in) :: r,z
  real (kind=8) , external :: psi_interpol
  real (kind=8) :: theta, n_theta, psi,n_r, alpha, r_s, z_s, d,del
  integer :: i

  ! n_theta
  ! finding ther value of (r,z) position
  theta=acos( (r-eq_axis_r)/dsqrt((r-eq_axis_r)**2+(z-eq_axis_z)**2) )
  if( z < 0D0 ) then
     theta=sml_2pi-theta
  endif
  ! Del value -- angle distance from x-point
  del=min(abs(theta-neu_theta_x),min( abs(theta-neu_theta_x+sml_2pi), abs(theta-neu_theta_x-sml_2pi)))
  ! theta dependance of density
  n_theta=neu_n0*(    1D0 + neu_delta_n*exp( -(del/neu_delta_theta)**2 )    )

  ! n_r
  ! radical dependance of density - distance from separatrix
  psi=psi_interpol(r,z,0,0)
  if(psi>eq_x_psi) then
     n_r=1
  else  ! distance from separatrix, ignoring non-orthogonal effect
     i=int(theta/sml_2pi*real(neu_sep_mtheta))+1
     i=min(neu_sep_mtheta-1,max(1,i))
     alpha=theta/sml_2pi*real(neu_sep_mtheta) + 1 - i
     r_s= (1D0-alpha)*neu_sep_r(i) + alpha*neu_sep_r(i+1)
     z_s= (1D0-alpha)*neu_sep_z(i) + alpha*neu_sep_z(i+1)
     d=dsqrt( (r-r_s)**2 + (z-z_s)**2 )
     n_r=dexp( - d/neu_mfp  )
     
  endif
  
  neutral_den1=n_theta*n_r
! print *, neu_n0, n_theta, n_r
  
end function neutral_den1



real (kind=8) function te_neutral_col(psi) 
! electron temperature in eV
  use eq_module, only : eq_x_psi
  implicit none
  real (kind=8), external :: tempe_ev
  real (kind=8) :: psi,xd
  
!  if(xd<1) then
!     te_neutral_col=tempe_static(1D0)*1D3
!  else
     te_neutral_col=tempe_ev(psi)
!  endif
end function te_neutral_col

     
subroutine change_neu_mfp(ptl)  !set new mean free path
  use ptl_module
  use diag_module, only : diag_flow_vol,diag_flow_npsi,diag_flow_pin,diag_flow_pout,diag_flow_dp
  use eq_module, only : eq_x_z,eq_x_psi,eq_den_edge
  use neu_module
  use sml_module, only: sml_mype
  implicit none
  type(ptl_type) :: ptl
  integer :: in, out,i,j
  real (kind=8) :: den(diag_flow_npsi),r,z,psi,dum(diag_flow_npsi),maxden
  real (kind=8) :: psi_interpol

  ! finding density as a function of psi
  den=0D0  ! den is weight sum for a while
  do i=1, ptl%ion%num
     r=ptl%ion%phase(1,i)
     z=ptl%ion%phase(2,i)
     
     if(ptl%ion%gid(i) > 0) then
        psi=psi_interpol(r,z,0,0)
        if(diag_flow_pin<psi .AND. psi<diag_flow_pout .AND. z>eq_x_z ) then
           j=int((psi-diag_flow_pin)/diag_flow_dp)+1
           j=min(diag_flow_npsi,max(1,j))

           den(j)=den(j)+ptl%ion%phase(8,i)
           
        endif
     endif
  enddo
  
  call my_mpi_allreduce(den,dum,diag_flow_npsi)
  ! den is density.
  den=dum/(diag_flow_vol) *1D-6  ! cm-3

  in=int( (0.98*eq_x_psi-diag_flow_pin)/diag_flow_dp ) + 1
  in=min(diag_flow_npsi,max(1,in))
  out=int( (eq_x_psi-diag_flow_pin)/diag_flow_dp   ) + 1
  out=min(diag_flow_npsi,max(1,out))

  maxden=maxval( den(in:out) ) ! finding maxmum density from 0.98*Psi_x ~ Psi_x
  
!  if(sml_mype==1 ) then
!     print *, neu_mfp, maxden, col_denb_edge,col_denb_edge/maxden,in,out,flow_stat_n
!  endif
  
  neu_mfp=neu_mfp0*eq_den_edge/maxden
  
  
end subroutine change_neu_mfp


!!$subroutine neutral_inside_weight
!!$  use ptl_module
!!$  use neu_module
!!$  real (kind=8) :: psi, psi_interpol
!!$  integer :: i
!!$  neu_old_inside_weight=0D0
!!$
!!$  do i=1, ptl_num
!!$     psi=psi_interpol(ptl_phase(i,1), ptl_phase(i,2),0,0)
!!$     if(psi < eq_x_psi) neu_old_inside_weight=neu_old_inside_weight+ptl_weight(i)
!!$  enddo
!!$end subroutine neutral_inside_weight

