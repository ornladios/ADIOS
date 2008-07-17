subroutine collision(istep,ptl)
  use col_module
  use ptl_module
  use sml_module
  implicit none
  integer, intent(in) :: istep
  type(ptl_type) :: ptl
  real (kind=8) :: r, z, phi, psi, b, rho, mu
  integer ::  i,flag,flag2
  real (kind=8), external :: b_interpol, psi_interpol

  if(mod(istep,col_period)==0) then
     flag=2* 2**2 - 1 ! all collision
     flag2= 1 + col_en_col_on*2    ! +1, pitch angle, +2 energy col

     if(col_mode==1) then
        do i=1, ptl%ion%num
           r=ptl%ion%phase(1,i)
           z=ptl%ion%phase(2,i)
           phi=ptl%ion%phase(3,i)
           b=b_interpol(r,z,phi)
           psi=ptl%ion%phase(9,i)
           rho = ptl%ion%phase(4,i)
           mu = ptl%ion%phase(5,i)
           ! to do : find correct iflag and iflag2 from input
           call scatr(psi,b,rho,mu,flag,flag2,ptl%ion)
        enddo
      else if (col_mode==2) then
        call conserving_collision(ptl,1)
     endif

     if(col_mode/=0 .and. sml_electron_on==1) then
        !for electron collision
        flag=4 ! ion collision
        flag2= 1    ! +1, pitch angle, +2 energy col
        do i=1, ptl%elec%num
           r=ptl%elec%phase(1,i)
           z=ptl%elec%phase(2,i)
           phi=ptl%elec%phase(3,i)
           b=b_interpol(r,z,phi)
           psi=ptl%elec%phase(9,i)
           rho = ptl%elec%phase(4,i)
           mu = ptl%elec%phase(5,i)
           ! to do : find correct iflag and iflag2 from input
           call scatr(psi,b,rho,mu,flag,flag2,ptl%elec)
        enddo
     endif
     
  endif
end subroutine

!cccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccc
real (kind=8) function den_ion(psi_in)
  use col_module
  implicit none
  real (kind=8) , intent(in) :: psi_in
  real (kind=8) , external :: init_den
  integer :: i
  real (kind=8) :: aa,bb,psi

  if (col_varying_bg==1) then
     psi=min(col_vb_pout, max(col_vb_pin,psi_in))
     i=int((psi-col_vb_pin)/col_vb_dp) +1
     i=min(col_vb_m-1,max(1,i))
     bb=(psi-col_vb_pin)/col_vb_dp + 1 - i
     aa= 1D0- bb
     den_ion=col_vb_den(i)*aa + col_vb_den(i+1)*bb
  else
     den_ion=init_den(psi)
  endif
end function den_ion

real (kind=8) function tempi_ev(psi_in)
  use col_module
  implicit none
  real (kind=8), intent(in) :: psi_in
  real (kind=8) , external :: init_tempi_ev
  integer :: i
  real (kind=8) :: aa,bb,psi

  if ( col_varying_bg==1) then
     psi=min(col_vb_pout, max(col_vb_pin,psi_in))
     i=int((psi-col_vb_pin)/col_vb_dp) +1
     i=min(col_vb_m-1,max(1,i))
     bb=(psi-col_vb_pin)/col_vb_dp + 1 - i
     aa= 1D0- bb
     tempi_ev=col_vb_temp(i)*aa + col_vb_temp(i+1)*bb
  else
     tempi_ev=init_tempi_ev(psi)
  endif

  tempi_ev=max(tempi_ev,1D0) !prevent minus temperature
end function tempi_ev


!cccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccc
real (kind=8) function den_imp(psi)
  use col_module 
  IMPLICIT NONE
  real (kind=8) , intent(in):: psi
  real (kind=8) :: a,b
  
  a=0.5D0 *(col_den_imp_edge-col_den_imp_out)
  b=0.5D0 *(col_den_imp_edge+col_den_imp_out)
  den_imp= a*dtanh(2D0*(col_den_imp_ped_c-psi)/col_den_imp_ped_width) +b

end function den_imp


!cccccccccccccccccccccccccccccccccccccc
! New scatr routine
!cccccccccccccccccccccccccccccccccccccc
subroutine scatr(psi,b,rho,mu,iflag,iflag2,sp)

  use sml_module
  use ptl_module
  use col_module
  IMPLICIT NONE

  real (kind=8), intent(in) :: psi,b
  real (kind=8), intent(inout) :: rho,mu
  real (kind=8), external :: tempe_ev, tempi_ev, den_ion, den_imp, ranx

  integer, intent(in) :: iflag,iflag2
  type(species_type), intent(in) :: sp

  real (kind=8) :: dt,xd,ekmin,ekin,pitch,dnb,dni,dne,ti,te,ti_code, te_code
  real (kind=8) :: dum,agg,esig,del_pitch
  real (kind=8) :: colb,colbs,fac0b,cole,coles,fac0e,coli,colis,fac0i
  real (kind=8) :: e_mass, e_charge     !electron mass, charge relative to proton
  real (kind=8) :: s_mass, s_charge
  e_mass=1.d0/1836d0   !relative mass of electron wrt proton
  e_charge=1.d0


! iflag : +2 impurity collsion, +1 electron collision, +2^(1+n) ion speices #n collision
! iflag2: +3 pitch angle + energy

#ifdef NAN_TEST
  !NaN test
  if((psi > 1D0 .or. psi < 2D0) .and. (b>1D0 .or.b<2D0)  )  then

  else
     print *,'col scatr entry', psi,b, rho,mu 
     stop
  endif
#endif

!!!-----------------------------------------------Common part for all - just one species

     s_mass=sp%mass
     s_charge=sp%charge


  dt = sml_dt * col_period

  ekmin = 1.d-3 * sml_ev2j


  call rho_mu_to_ev_pitch2(rho,mu,b,ekin,pitch,sp%type)
  ekin=ekin*sml_ev2j         ! conversion to code unit
  ekin = max(ekmin, ekin)        ! 2002/09/17 added for preventing small energy

  dni=den_imp(psi) ! impurity density m-3
  dne=dni*col_imp_on*col_imp_charge
  dnb=den_ion(psi) ! background density m-3
  dne=dne + dnb*ptl_charge_eu !!!+ col_imp*dni*col_chgi !electron density


!!!-----------------------------------------------collision with background ions
     
  if(mod(iflag/(2*2),2)==1) then
     
     call find_freq(ekin, s_mass, s_charge, dnb, tempi_ev(psi), &
             & ptl_mass_au, ptl_charge_eu, colb, colbs, fac0b) !!! energy in [eV]

     if(mod(iflag2,2)==1) then  !pitch angle scattering
        !!print*, 'ion--pitch called'
        agg = ranx() - .5d0
        dum = 1.d0 - pitch**2
        dum = max(0.d0, dum)
        del_pitch = dsign(1.d0,agg)*sqrt(dum*colb*dt*0.5d0)
        pitch = pitch*(1.d0-colb*dt*0.5d0) + del_pitch
     endif
     if(mod(iflag2/2,2)==1) then !energy collision
        !!print*, 'ion--energy called'
        agg = ranx()- .5d0
        esig = dsign(1.d0,agg)*sqrt(2*ekin*tempi_ev(psi)*sml_ev2j*colbs*dt)
        ekin = ekin - colbs*dt*fac0b + esig
        
        ekin = max(ekin,ekmin)
        if(ekin > 1D0 .or. ekin < 2D0 )  then
           
        else
           print *,'colbs',psi,rho,mu,colbs, fac0b, esig
           ekin=ekmin
           stop
        endif
     endif
     
  endif


!!!-----------------------------------------------collision with electrons
  if(mod(int(iflag),2)==1) then
     call find_freq(ekin, s_mass, s_charge, dne, tempe_ev(psi), e_mass, e_charge, &
          & cole, coles, fac0e)
     if(mod(iflag2/2,2)==1) then !energy collision
!     if(col_en_col_on ==1) then
        agg = ranx()-0.5d0
        esig = dsign(1.d0,agg)*sqrt(2.d0*ekin*tempe_ev(psi)*sml_ev2j*coles*dt)
        ekin = ekin - coles*fac0e*dt + esig
        ekin = max(ekin,ekmin)   ! truncate slowing down
     endif
     agg = ranx() - .5D0
     dum = 1.d0 - pitch**2
     dum = max(0.d0, dum)
!     if(dum > 1. .or. dum < 2.) then
!     else
!        print *, 'nan found in col ele', ekin*j2ev,pitch,b,rho,mu,iflag,iflag2
!        print *, 'back freqs', colb,colbs,fac0b
!        print *, 'elec freqs', cole,coles,fac0e,del_pitch,dum,agg
!        stop
!     endif
     del_pitch = dsign(1.D0,agg)*sqrt(dum*cole*dt*.5D0)
     pitch = pitch*(1.D0 - cole*dt*.5D0) + del_pitch

  endif
!  if(pitch > 1. .or. pitch < 2.) then
!  else
!     print *, 'nan found in col', ekin*j2ev,pitch,b,rho,mu,iflag,iflag2
!     print *, 'back freqs', colb,colbs,fac0b
!     print *, 'elec freqs', cole,coles,fac0e,del_pitch,dum,agg
!     stop
!  endif

!!!-----------------------------------------------collision with impurities
  if(mod(int(iflag/2),2)==1) then
     call find_freq(ekin, s_mass, s_charge, dni, tempi_ev(psi), &
          & col_imp_mass, col_imp_charge, coli, colis, fac0i)

     if(mod(iflag2/2,2)==1) then !energy collision
!     if(col_en_col_on==1) then
        agg = ranx()-0.5d0
        esig = dsign(1.d0,agg)*sqrt(2.d0*ekin*tempi_ev(psi)*sml_ev2j*colis*dt)
        ekin = ekin - colis*fac0i*dt + esig
        ekin = max(ekin,ekmin)   ! truncate slowing down
     endif
     agg = ranx() - .5D0
     dum = 1.d0 - pitch**2
     dum = max(0.d0, dum)
     del_pitch = dsign(1.D0,agg)*sqrt(dum*coli*dt*.5D0)
     pitch = pitch*(1.D0 - coli*dt*.5D0) + del_pitch
  endif

!!!----------------------------------------------final rho and mu
  
  pitch = min(1.D0,pitch)
  pitch = max(-1.D0,pitch)

  call ev_pitch_to_rho_mu2(ekin*sml_j2ev, pitch, b, rho, mu,sp%type)
  
!  if(rho > 1. .or. rho < 2.) then
!  else
!     print *, 'nan found in col', ekin*sml_j2ev,pitch,b,rho,mu,iflag,iflag2
!     print *, 'back freqs', colb,colbs,fac0b
!     print *, 'elec freqs', cole,coles,fac0e
!     print *, 'imp  freqs', coli, colis, fac0i,del_pitch,dum,agg
!     stop
!  endif
  
end subroutine scatr



! This routine uses cgs unit -- be careful
subroutine find_freq(en_a, mass, charge, dn_b, en_b_ev, mass_b, charge_b, freq_scat, freq_slow, freq_fac0)
  ! calculate collision frequencies and convert into code unit
  ! approx to psi function
  ! pitch angle scattering, small value of col!  (alpha)
  ! Boozer and Kuo-Petravic Phys. Fluids 24, 851 (1981)
  ! profiles-density (cm-3) dnb (background), dni (impurity), temp (kev)
  ! psi(x) = (2/sqrt(pi)Int[dt*sqrt(t)*exp(-t)]
  ! psi(x) = 1 - 1/(1 + p),
  ! p = x**1.5[a0 + a1*x +a2*x**2 +a3*x**3 +a4*x**4]
  ! R.White,M.Redi Jan (2000) error in psi(x) < 1.D-3, asymptotically correct
  ! relative error dpsi/psi < 1.D-2

!  use NC_module, only : nc_norm_r, nc_norm_t
  use sml_module !, only : sml_prot, sml_zprt, sml_en_order_kev, sml_en_order
!  use col_module, only : col_massb, col_massi, col_accel
  implicit NONE

  real(kind=8), intent(in) :: en_a, mass, charge, mass_b, charge_b, dn_b, en_b_ev !, mass_b, charge_b
  real(kind=8), intent(out) :: freq_scat, freq_slow, freq_fac0
  real(kind=8) :: dn, dumb, dd,d3, dum, ap0, ap1, ap2, ap3, ap4, ee_ev, cnst, vprt, &
       & vt, ap_psi, f, g, gp, bmax, bmin1, massab, bmin2, bmin, clog, dnu_b
  data ap0 /.75225/,ap1 /-.238/,ap2 /.311/,ap3 /-.0956/,ap4 /.0156/

  vprt = sqrt(2D0*en_a/mass)
!!  vt = dsqrt(2d0*(en_b_ev/1.d3*sml_en_order/sml_en_order_kev)*mass/mass_b)
  vt = sqrt(2.d0*(en_b_ev*sml_ev2j)/mass_b)

!!  norm_r_cgs = nc_norm_r * 100.d0 !!! conversion MKS to cgs
!!  cnst = 2.41D11*(charge/mass)**2/(norm_r_cgs/nc_norm_t*vprt)**3  !100 is conversion MKS to CGS
!!  cnst = cnst*col_accel

  cnst = 2.41D11*(charge/mass)**2/(100d0*vprt)**3  !100 is conversion MKS to CGS
  dn = dn_b/1d6

!!  cnst = cnst*col_accel

  ee_ev = en_a*sml_j2ev   !!! test ptl's kinetic energy in [eV]

  ! calculate psi(x), f, g, gp
  dumb = vprt/vt
  dd = dumb**2   ! dd = x, dumb = sqrt(x)
  d3 = dumb**3
  dum = d3*(ap0 +ap1*dd +ap2*dd**2 +ap3*dd**3 +ap4*dd**4)
  ap_psi = 1.D0 - 1.D0/(1.D0 + dum) !psi(x)

  f = (2.D0 - 1.D0/dd)*ap_psi + 2.257D0*dumb*dexp(-dd) 
  g = 2D0*mass*ap_psi/mass_b
  gp = 2.257D0*dexp(-dd)*(mass/mass_b)*dumb

  ! find Coulomb logarithm
!!  bmax  = 7.4d-3*dsqrt(en_b_ev/1000d0*1.d13/(charge**2*dn_b))
!!  bmin1 = 1.4d-10*charge*charge_b/(ee_ev + en_b_ev)*1000d0
!!  massab = mass*mass_b/(mass+mass_b)
!!  bmin2 = 6.2d-10/(dsqrt(massab*1836.d0*en_b_ev/1000d0))


  bmax  = 7.4d-3*sqrt(en_b_ev*1.d10/(charge**2*dn))
  bmin1 = 1.4d-7*charge*charge_b/(ee_ev + en_b_ev)
  massab = mass*mass_b/(mass+mass_b)
  bmin2 = 6.2d-10/(sqrt(massab*1836.d0*en_b_ev/1.d3))
  bmin = max(bmin1,bmin2)
  clog = dlog(bmax/bmin)

  ! collision frequencies - scattering, slowing down
  dnu_b = cnst*dn*charge_b**2
  freq_scat = dabs(clog*dnu_b*f)
  freq_slow = dabs(clog*dnu_b*g)
  freq_fac0 = en_a*(1.d0 - mass_b*gp/(g*mass)) 

  return

end subroutine find_freq

!get ion density and temperature for collision routine
subroutine col_snapshot(ptl)
  use eq_module, only : eq_x_psi, eq_x_z, eq_den_out
  use ptl_module
  use col_module
  use sml_module,only : sml_mype, sml_j2ev,sml_deltaf
  use perf_monitor
  implicit none
  type(ptl_type) :: ptl
  integer :: i,j
  real (kind=8) :: dum(col_vb_m)
  real (kind=8) :: B,br,bz,bphi,e
  real (kind=8) :: psi,r,z,weight,phi,aa,bb
  real (kind=8),external :: psi_interpol

  col_vb_temp=0D0
  col_vb_den=1D-99  ! to avoid divide by zero
  do i=1, ptl%ion%num
     r=ptl%ion%phase(1,i)
     z=ptl%ion%phase(2,i)
     phi=ptl%ion%phase(3,i)

     if(ptl%ion%gid(i) > 0) then
        psi=psi_interpol(r,z,0,0)
        if(col_vb_pin<psi .AND. psi<col_vb_pout .AND. z>eq_x_z ) then ! 2002/10/10 , region 2 condition modified
           j=int((psi-col_vb_pin)/col_vb_dp) +1
           j=min(col_vb_m-1,max(1,j))
           bb=(psi-col_vb_pin)/col_vb_dp + 1 - j
           aa= 1D0- bb
           !energy calculation
           call bvec_interpol(r,z,phi,br,bz,bphi)
           B=dsqrt(br**2+bz**2+bphi**2)
           E=0.5D0*(ptl%ion%phase(4,i)*B)**2 + ptl%ion%phase(5,i)*B
           if(sml_deltaf==1 .and. col_varying_bg==1) then 
              weight=ptl%ion%phase(6,i)*ptl%ion%phase(8,i)
           else
              weight=ptl%ion%phase(8,i)
           endif
           col_vb_temp(j)=col_vb_temp(j)+ weight*E*aa  ! temp is energy sum
           col_vb_den(j)=col_vb_den(j)+weight*aa    
           col_vb_temp(j+1)=col_vb_temp(j+1)+ weight*E*bb  ! temp is energy sum
           col_vb_den(j+1)=col_vb_den(j+1)+weight*bb
        endif
     endif
  enddo
  
  call monitor_start (COL_SNAP_RED_) 
  call my_mpi_allreduce(col_vb_temp,dum,col_vb_m)
  col_vb_temp=dum
  call my_mpi_allreduce(col_vb_den,dum,col_vb_m)
  col_vb_den=dum
  call monitor_stop (COL_SNAP_RED_)


  if(sml_deltaf==1 .and. col_varying_bg==1) then
     print *, 'not implimented yet'
     stop
  endif

  col_vb_temp=col_vb_temp/dum ! temp is average energy per ptl
  col_vb_temp=col_vb_temp * sml_j2ev *2D0/3D0
  col_vb_temp=max(col_vb_temp,1D-3) !10/21/03
  col_vb_den=dum/col_vb_vol ! m-3
  col_vb_den=max(col_vb_den,0.5D0*eq_den_out) !05/05/05

  if (sml_mype==0 ) then
     do i=1,col_vb_m
        psi=col_vb_pin+(real(i)-0.5D0)*col_vb_dp
        write (34,2000) psi/eq_x_psi,col_vb_den(i),col_vb_temp(i),col_vb_vol(i)
     enddo
     write (34,*) ' '
     !close(34) !debug
  endif

2000 format(4(e19.13,' '))

end subroutine col_snapshot

