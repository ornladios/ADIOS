! last modified 2002/05/18
subroutine rf_save_position(ptl)
  use rf_module
  use ptl_module
  implicit none
  type(ptl_type) :: ptl
  rf_r_save(1:ptl%ion%num)=ptl%ion%phase(1,1:ptl%ion%num)

end subroutine rf_save_position


!last modified 2002/05/18
! not verified -- check the formula and unit
subroutine rf(ptl)
  use rf_module
  use ptl_module
  use sml_module
  implicit none
  type(ptl_type) :: ptl
  integer :: m,sp
  real (kind=8) :: r_old, r_new,z,r_mid
  real (kind=8) :: delta_en,v_r,Q_rf_dt
  real (kind=8) :: E_rf,b,mu
  real (kind=8), external :: b_interpol,ranx

  if(rf_on==1) then
     do m=1, ptl%ion%num
        !     if(mass(m)>0.45 .AND. mass(m)<0.55) then ! hydrogen only
        R_old=rf_r_save(m)
        R_new=ptl%ion%phase(1,m)
        if( (R_old - RF_R_res)*(R_new - RF_R_res) < 0 ) then
           z=ptl%ion%phase(2,m)
           mu=ptl%ion%phase(5,m)
           r_mid=0.5D0*(r_old+r_new)
           b=b_interpol(r_mid,z) 
           v_R=abs((R_new-R_old)/sml_dt)
           E_rf=rf_E0*max( (1D0-((z-rf_za)/(0.5*rf_lap))**4),0D0)
           !v_R minimum
           v_R=max(v_R,1D-4*sqrt(sml_en_order))
           Q_rf_dt=0.5D0*sml_pi*(ptl%ion%charge*E_rf)**2 * R_mid/v_R*sml_dt 
           delta_en= 0.5D0* Q_rf_dt + 0.5D0*(ranx()-0.5D0) * sqrt(48D0* mu*B*Q_rf_dt)
!           print *, ptl_phase(m,5)*B*sml_norme2ev, delta_en*sml_norme2ev
           ptl%ion%phase(5,m)=ptl%ion%phase(5,m) + delta_en/B
           ptl%ion%phase(5,m)=max(ptl_phase(5,m),0D0)
        endif
     enddo
  end if

end subroutine rf

