! ***************************************************
subroutine neutral_step(mass)
!
! 2-dimensional neutral particle simulation routine
! After each evaluation of the plasma profile, stop the plasma routine
! and evaluate neutral density and temperature profiles for the same physical
! duration.
! Assume a hypothetical boundary at Psi_hypo=1.03 or z=zmin.  At Psi >= Psi_hypo
! or z<= zmin, assume Maxwellian neutral distribution function at room 
! temperature.
! **************************************************
  use sml_module
  use eq_module
  use neu_module
  use sml_module
  use ptl_module
  implicit none
  real (kind=8), intent(in) :: mass
  real (kind=8) :: rstart,zstart,slope,r0,z0,r,z,vr,vz,vangle,tn,v,pol
  integer :: ipol,k,tstep,ipart,ptest,istep
  real (kind=8) :: f_ion, f_cx, f_el,rate_ion,rate_cx,rate_el,drate
  real (kind=8) :: denp,teev,tiev,enev,ranx,xmax,sigma0,jacob,tekev
  real (kind=8) :: weight, den, del,theta,psi_interpol
!  integer , parameter :: mpol=1000
! neutral-ion Elastic scattering crosssection
!  print *, 'neutral_step'

  if(sml_mype==0) print *, 'Neutral profile renewal...'
  sigma0=5d-19/sqrt(2D0)  !m**2
  
  neu_grid_den=0D0
  neu_grid_temp=0D0
  neu_grid_wsum=0D0
  neu_grid_esum=0D0


!  slope=sml_pi*0.5D0
  r0=1d0 
  z0=0d0

  ! Make the neutral time equal to the plasma time interval
  ! Initiation position
  do ipol=1,neu_mpol
!     print *, ipol
     pol=2.d0*sml_pi*(real(ipol-1)+ real(sml_mype)/real(sml_totalpe))/real(neu_mpol)  ! mype is for phase difference
  ! start simulation at each initiation position
     do k=1, neu_istep_max
        tstep=k-1
        do ipart=1, neu_num
           
           ! neut_num is the particle number per processor
           call startrz(pol,rstart,zstart,jacob,slope)
           !Get (R,Z) of psi=1.03 at the starting position.
           r=rstart
           z=zstart
           tn=neu_temp0
           call maxwell_neut(r0,z0,r,z,slope,tn,vr,vz,vangle)    
           v=sqrt(vr**2+vz**2)

           ! get theta value and neutral density
           theta=dacos( (r-eq_axis_r)/dsqrt((r-eq_axis_z)**2+(z-eq_axis_z)**2) )
           if( z < 0D0 ) then
              theta=sml_2pi-theta
           endif
           del=min(abs(theta-neu_theta_x),min( abs(theta-neu_theta_x+sml_2pi), abs(theta-neu_theta_x-sml_2pi)))
           den=(1D0 +(neu_delta_n-1.d0)*exp( -(del/neu_delta_theta)**2 ) ) ! density normalize to neu_n0


           ! weight calculation -----------------------------------
           weight=jacob*sml_2pi/real(neu_mpol)*  &  ! area of segment
                den* &  !neutral density in norm unit
                v*abs(sin(slope-vangle))*neu_dt  & !
                /real(neu_num*sml_totalpe) !---------------------

           do istep=tstep, neu_istep_max-1
              ! move forward.
              r= r+vr*neu_dt
              z= z+vz*neu_dt

              ! If particle is outside the hypothetical boundary, stop simulation.
              ! If ptest=-2, neutral entered divertor.
              ! If ptest=-1, neutral is on or outside Psi = 1.03.
              call neut_pos_test(r,z,ptest)
!              if(ptest<0) write(143,*) r,z
              if(ptest < 0) go to 100
              call plasma(r,z,denp,teev,tiev)  !te, ti-- ev, denp -- mks
              !              print *, denp,te,ti
              !------------------------------------------------
              ! Get local collision probability of the neutral particle
		    !
              ! Probability of ionization per sec 
              
              f_ion=denp*0.8D-8*sqrt(Teev)*(exp(-13.56/Teev)/(1D0+0.01D0*Teev))*1D-6  
              ! denp is in MKS here.
              f_ion=f_ion ! s^-1 --> conversion to normalized unit

              ! probability of charge exchange per sec 
              enev=tiev
              f_cx= denp * 1.1D-8* (enev)**0.3 /sqrt(mass) *1D-6   
              ! 1.1 K_i(eV)^0.3 * 10^-8 / Sqrt(Mi(AU))  cm s^-1   -> 1D-6 m s^-1
              f_cx = f_cx  ! s^-1 --> conversion to normalized unit. 

              !probablity of elastic collision
              f_el = denp*sigma0*v*0.5D0*sml_pi  !check !!!!

              ! -----------------------------------------------------
              ! Rate parameters
              drate= neu_dt*f_ion
              rate_ion=drate-drate**2/2.d0+drate**3/6.d0-drate**4/24.d0
              drate= neu_dt*f_cx
              rate_cx=drate-drate**2/2.d0+drate**3/6.d0-drate**4/24.d0
              drate= neu_dt*f_el
              rate_el=drate-drate**2/2.d0+drate**3/6.d0-drate**4/24.d0
              
              !write(123,*) rate_ion,rate_cx,rate_el
!              rate_cx=0D0
!              rate_el=0D0
!              rate_ion=0D0

              ! Ionization of the neutral particle.  
              !  The neutral disappears upon ionization.
!             if(rate_ion .ge. ranx().and. psi_interpol(r,z,0,0) < NC_x_psi) go to 100  !2003/08/12 - ionize inside only
              if(rate_ion .ge. ranx()) go to 100  !2003/08/12 - ionize everywhere
              ! Charge Exchange of the neutral particle.
              ! Vangle is randomized and te energy jumps to local Ti.
              if(rate_cx .ge. ranx()) then
                 vangle=vangle+2.d0*sml_pi*ranx()
                 v=dsqrt(2.d0*tiev*sml_ev2j/ptl_mass(1))
                 vr=v*dcos(vangle)
                 vz=v*dsin(vangle)
              endif
              
              
              ! Elastic collision of the neutral particle.
              ! Vangle is randomized keeping the same energy.
              if(rate_el .ge. ranx()) then
                 vangle=vangle+2.d0*sml_pi*ranx()
                 vr=v*dcos(vangle)
                 vz=v*dsin(vangle)
              endif

           enddo
           call neut_fluid(r,z,vr,vz,weight)  ! Record neutral density and temperature
100        continue
        enddo
     enddo
!     write(144,*) theta,den
  enddo

  call neut_fluid_final
  if(sml_mype==0) print *, 'Neutral profile modified.'
end subroutine



Subroutine maxwell_neut(r0,z0,r,z,slope,tn,vr,vz,vangle)
  !***********************************************************
  ! SLOPE is the local slope of the imaginary wall at Psi = 1.03
  ! Created on 5/17/2003
  !**********************************************************
  use sml_module
  use ptl_module
  implicit none
  real(kind=8) , intent(in) :: r,z,r0,z0,slope, tn
  real(kind=8) , intent(out) :: vr, vz, vangle
  real (kind=8) :: xdum,xmax, ranx, v,tn2
  ! r-z velocity
  xmax=1.d0-dexp(-7.d0)      
  xdum=xmax*ranx()
  tn2=tn*sml_ev2j/ptl_mass(1)
  v=dsqrt(-2.d0*dlog(1.d0-xdum)*tn2)
  vangle=2.d0*sml_pi*ranx()
!  vangle=slope+0.5D0*sml_pi  !debug   
  vr=v*dcos(vangle)
  vz=v*dsin(vangle)
!  if( (r.gt.r0.and. (slope .ge. vangle .or. vangle .ge. (slope+sml_pi))) .or. &
!  (r.le.r0 .and. slope .le. vangle .and. vangle .le. (slope+sml_pi))) then
!     vangle=vangle+sml_pi
!     vr=v*dcos(vangle)
!     vz=v*dsin(vangle)
!  endif

end subroutine maxwell_neut

subroutine neut_fluid(r,z,vr,vz,weight)
  use neu_module
  use eq_module
  use sml_module
  implicit none
  real (kind=8),intent(in) :: r,z,vr,vz,weight
  real (kind=8), external :: psi_interpol
  real (kind=8) :: psi,theta,rs,a1,a2,energy
  integer :: ipsi,itheta,itheta_p1

  psi=psi_interpol(r,z,0,0)
  rs=(r-eq_axis_r)**2 + (z-eq_axis_z)**2
  theta=acos( (r-eq_axis_r)/sqrt(rs) )
  if(z<0)  theta = sml_2pi -theta
  
  if( psi  > neu_grid_min_psi .and. psi<= neu_grid_max_psi*1.0001D0 .and. z>eq_x_z) then
     
     ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1 
     a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi + 1
     ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
     
     itheta= int(theta/neu_grid_dtheta) + 1 
     itheta= min(neu_grid_mtheta, max(1, itheta))
     a2= theta/neu_grid_dtheta - itheta +1
     itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1

     
     neu_grid_wsum(ipsi,itheta)=neu_grid_wsum(ipsi,itheta) + weight*(1D0-a1)*(1D0-a2)
     neu_grid_wsum(ipsi+1,itheta)=neu_grid_wsum(ipsi+1,itheta) + weight*a1*(1D0-a2)
     neu_grid_wsum(ipsi,itheta_p1)=neu_grid_wsum(ipsi,itheta_p1) + weight*(1D0-a1)*a2
     neu_grid_wsum(ipsi+1,itheta_p1)=neu_grid_wsum(ipsi+1,itheta_p1) + weight*a1*a2

     energy=0.5D0*(vr**2 + vz**2)*weight
     
     neu_grid_esum(ipsi,itheta)=neu_grid_esum(ipsi,itheta) + energy*(1D0-a1)*(1D0-a2)
     neu_grid_esum(ipsi+1,itheta)=neu_grid_esum(ipsi+1,itheta) + energy*a1*(1D0-a2)
     neu_grid_esum(ipsi,itheta_p1)=neu_grid_esum(ipsi,itheta_p1) + energy*(1D0-a1)*a2
     neu_grid_esum(ipsi+1,itheta_p1)=neu_grid_esum(ipsi+1,itheta_p1) + energy*a1*a2  ! bug fix 2003/06/07

  endif
  
  
!  write(141,*) r,z,vr,vz
  
end subroutine neut_fluid


subroutine neut_fluid_final
  use neu_module
  use sml_module
  use eq_module
  implicit none
  real (kind=8) :: wsum(neu_grid_mpsi,neu_grid_mtheta), esum(neu_grid_mpsi, neu_grid_mtheta)
  integer :: i,j
  real (kind=8) , external :: tempi_ev
  save wsum, esum

  call  my_mpi_allreduce(neu_grid_wsum, wsum, neu_grid_mtheta*neu_grid_mpsi)
  call  my_mpi_allreduce(neu_grid_esum, esum, neu_grid_mtheta*neu_grid_mpsi)

  neu_grid_den=wsum/neu_grid_vol
  neu_grid_temp=esum/wsum
  do j=1, neu_grid_mtheta
     do i=1, neu_grid_mpsi
         if(neu_grid_den(i,j)==0D0) then
		neu_grid_temp(i,j)=tempi_ev(neu_grid_dpsi*real(i-1)+neu_grid_min_psi)*sml_ev2j
	 endif
     enddo
  enddo	
  if(sml_mype==0) then
     do j=1,neu_grid_mtheta 
        do i=1,neu_grid_mpsi
           write(145,*) (neu_grid_dpsi*real(i-1)+neu_grid_min_psi)/eq_x_psi,neu_grid_dtheta*real(j-1),&
                neu_grid_den(i,j),neu_grid_temp(i,j)*sml_j2ev,neu_grid_vol(i,j),wsum(i,j)
        enddo
        write(145,*) ' '
     enddo
  endif
  close(145)

end subroutine


subroutine plasma(r,z,denp,teev,tiev)
  use sml_module
  implicit none
  real (kind=8) ,intent(in):: r,z
  real (kind=8) , intent(out) :: denp,teev,tiev
  real (kind=8) , external :: psi_interpol,tempe_ev,tempi_ev,den_ion
  real (kind=8) :: psi,xd

  psi=psi_interpol(r,z,0,0)
  teev=tempe_ev(psi)  ! collision.f90 -- ev
  tiev=tempi_ev(psi)  ! collision.f90 -- ev
  denp=den_ion(psi) ! collision.f90   MKS  m^-3
  teev=max(teev,1D-9)
  tiev=max(tiev,1D-9)
  denp=max(denp,1D0) 
 
!  denp0,denp1,te0,te1,ti0,ti1,width,rs

  ! Plasma density and temperature in a test case
  ! ne=ni assumed
  
!  denp0=5.d19 !pedestal plasma density in mks
!  denp1=2.d18 !scrape-off plasma density in mks
!  te0=0.5d0 !keV
!  te1=0.01d0 !kev
!  ti0=1d0 !keV
!  ti1=0.02d0 !kev
  
!  width=0.01 !pedestal width in meter
!  rs=0.03 !separatrix distance in mks from the ficitious wall at r=0.
!  denp=denp0*(dtanh(2d0*(r-rs)/width))+denp0+denp1
!  te=te0*(dtanh(2d0*(r-rs)/width))+te0+te1
!  ti=ti0*(dtanh(2d0*(r-rs)/width))+ti0+ti1
end subroutine plasma

subroutine  neut_pos_test(r,z,ptest)
  use eq_module
  use sml_module
  implicit none
  real (kind=8) :: r,z
  integer :: ptest
  real (kind=8) :: theta, rw,zw,tmp1,tmp2,rs
  
  rs=(r-eq_axis_r)**2 + (z-eq_axis_z)**2
  theta=dacos((r-eq_axis_r)/dsqrt(rs))
  if(z < 0D0) then
     theta= sml_2pi - theta
  endif
  
  call startrz(theta,rw,zw,tmp1,tmp2)
  if( rs*0.9999D0 <= (rw-eq_axis_r) **2 +(zw-eq_axis_z)**2 ) then
     ptest=1
  else
     ptest=-1
  endif
endsubroutine neut_pos_test


subroutine startrz(theta,r,z,jacob,slope)
  use neu_module
  use sml_module
  implicit none
  real (kind=8) ,intent(in) :: theta
  real (kind=8) ,intent(out) :: r,z,jacob,slope
!  integer ,intent(out) :: flag
  integer :: i,j
  real (kind=8) :: aa,dr,dz
!  flag=1
!  ! check if theta is in valid range
!  if( neu_theta1 < theta .AND.  theta <neu_theta2 ) then
!     flag=0
!     return
!  endif

  ! find proper index of pol
  i= int(theta/neu_dtheta) + 1
  i=min(neu_mtheta,max(1,i))
  aa= theta/neu_dtheta -real(i-1)

  j=i+1
  if(i+1 > neu_mtheta)  j=1 ! cyclic
     
  r=(1D0-aa)*neu_r0(i) + aa*neu_r0(j)
  z=(1D0-aa)*neu_z0(i) + aa*neu_z0(j)
  jacob=(1D0-aa)*neu_jacob(i) + aa*neu_jacob(j)
!  jacob=1D0
  dr=neu_r0(j)-neu_r0(i)
  dz=neu_z0(j)-neu_z0(i)
  slope= acos( dr/sqrt(dz**2+dr**2) )
  if(dz<0) then
     slope= sml_2pi - slope 
  endif
  
end subroutine startrz


subroutine neutral2_setup
  use neu_module
  use eq_module
  use sml_module
  implicit none
  real (kind=8):: max_r,min_r,max_z,min_z,theta,zd,rd,right,left,mid,r_tmp,z_tmp,&
       max_psi,max_psi_r,max_psi_z,factor,psi,tantheta
  real (kind=8) ,external :: psi_interpol
  integer :: i,j,ip,im,ifactor

  ! for each angle theta, find psi=neu_psi0 (1.03?) position
  ! using binary search
  ! left is magnetic axis and right is ?
  max_r=sml_bd_max_r
  min_r=sml_bd_min_r
  max_z=sml_bd_max_z
  min_z=0.5*(eq_x_z+sml_bd_min_z   )
  min_z=max(min_z, sml_bd_min_z)

  do i=1, neu_mtheta
     theta=real(i-1)*neu_dtheta
!     print *,theta/sml_2pi
     left=0D0
     if( theta <sml_pi) then
        zd=max_z
     else
        zd=-min_z
     endif
     if( theta < 0.5D0*sml_pi .or. theta > 1.5D0*sml_pi ) then
        rd=max_r-eq_axis_r
     else
        rd=eq_axis_r-min_r
     endif
     
!     right=sqrt( (rd*cos(theta))**2 + (zd*sin(theta))**2 )
     tantheta=tan(theta)
     right=sqrt( (1D0 + tantheta**2) / (1D0/rd**2 + tantheta**2/zd**2) )
!     print *,'R1',right,sml_bd_max_r-1D0,rd
     max_psi=0D0
     do ifactor=50, 100
        factor = ifactor/100D0        
        !linear search
        r_tmp=factor*right*cos(theta)+1D0
        z_tmp=factor*right*sin(theta)
!        if( z_tmp < eq_x_z ) then
!           r_tmp= (r_tmp-1D0)*eq_x_z/z_tmp +1D0
!           z_tmp=eq_x_z
!        endif
        psi=psi_interpol(r_tmp,z_tmp,0,0)
        if( z_tmp < eq_x_z ) then  ! xd_z cut-off
           psi=neu_psi_edge * (1D0 + (z_tmp-eq_x_z)**2)
        endif
        if(max_psi < psi) then
           max_psi_r=r_tmp
           max_psi_z=z_tmp
           max_psi=psi
        endif
     enddo
!     print *,'R2',r_tmp
     if(max_psi < neu_psi_edge) then
        neu_r0(i)=max_psi_r
        neu_z0(i)=max_psi_z
     else
           
        !binary search
        right=sqrt((max_psi_r-eq_axis_r)**2+(max_psi_z-eq_axis_z)**2)
        left=0.3D0*right
        
        do j=1, 30
           mid=0.5D0*(right+left)
           r_tmp=mid*cos(theta)+1D0
           z_tmp=mid*sin(theta)
           psi=psi_interpol(r_tmp,z_tmp,0,0)
           if( z_tmp < eq_x_z ) then  ! xd_z cut-off
              psi=neu_psi_edge * (1D0 + (z_tmp-eq_x_z)**2)
           endif
           if(psi < neu_psi_edge) then
              left=mid
           else
              right=mid
           endif
        enddo
        neu_r0(i)=mid*cos(theta)+1D0
        neu_z0(i)=mid*sin(theta)
        
     endif

  enddo

  do i=1, neu_mtheta
     !jacobian ? --> area (ignore psi direction)
     ip=mod(i+1  -1,neu_mtheta) +1
     im=mod(i-1  -1+neu_mtheta,neu_mtheta) +1     
     neu_jacob(i)=sml_2pi*neu_r0(i) *0.5D0/neu_dtheta*( &
!     neu_jacob(i)=sml_2pi*0.5D0/neu_dtheta*( &
          sqrt( (neu_r0(ip)-neu_r0(i))**2 + (neu_z0(ip)-neu_z0(i))**2 ) + &
          sqrt( (neu_r0(i)-neu_r0(im))**2 + (neu_z0(i)-neu_z0(im))**2 )  )
     
!     if(sml_mype==0) then
!        write(190,*) neu_r0(i),neu_z0(i),neu_jacob(i)
!     endif
  enddo
        
  call get_volume_neu   ! monte-carlo volume calculation
  
end subroutine neutral2_setup
     
  
subroutine get_volume_neu
  use ptl_module
  use sml_module
  use eq_module
  use neu_module
  implicit none
  integer valid, total,itheta,ipsi,itheta_p1
  real (kind=8) :: rdim, zdim, roffset, zoffset, tvolume, vvolume, &
       r,z,psi,phi,dum(neu_grid_mpsi,neu_grid_mtheta),dum1(1),dum2(1),theta
  real (kind=8) , external :: ranx, psi_interpol
  real (kind=8) :: a1,a2

  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z

  valid=0
  total=0
  
  neu_grid_vol=1D-50   
  do while(valid<neu_monte_num)
     r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) ) 
     z=zdim*ranx() + zoffset
!     r=rdim*ranx() + roffset  ! flat   --- debug
     psi=psi_interpol(r,z,0,0)     
     theta=dacos((r-eq_axis_r)/dsqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
     if(z < 0D0) then
        theta= sml_2pi - theta
     endif
     total=total+1
     if(neu_grid_min_psi < psi .AND. psi < neu_grid_max_psi .and. z>eq_x_z) then
        valid=valid+1

        ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1 
        ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
        a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi +1
        
        itheta= int(theta/neu_grid_dtheta) + 1 
        itheta= min(neu_grid_mtheta, max(1, itheta))
        a2= theta/neu_grid_dtheta - itheta +1
        itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1
        
        neu_grid_vol(ipsi,itheta)=neu_grid_vol(ipsi,itheta) + (1D0-a1)*(1D0-a2)
        neu_grid_vol(ipsi+1,itheta)=neu_grid_vol(ipsi+1,itheta) + a1*(1D0-a2)
        neu_grid_vol(ipsi,itheta_p1)=neu_grid_vol(ipsi,itheta_p1) + (1D0-a1)*a2
        neu_grid_vol(ipsi+1,itheta_p1)=neu_grid_vol(ipsi+1,itheta_p1) + a1*a2

     endif
  enddo

  ! neu_grid_vol sum
  call my_mpi_allreduce(neu_grid_vol,dum,neu_grid_mpsi*neu_grid_mtheta)
  neu_grid_vol=dum
  ! total sum
  dum1(1)=real(total)
  call my_mpi_allreduce(dum1,dum2,1) !dum2(1) is sum of total


!  sml_marker_den= real(ptl_num)*dum2(1)/ (  rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) * sml_monte_num)  ! normalized unit  is normalization constant
  
  ! volume of each flux shell -- shell volume = (ptl number of each shell) / (sum of total) * volume
    neu_grid_vol=neu_grid_vol/dum2(1)*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) ) ! bug fix 2002/07/04

end subroutine 



 
