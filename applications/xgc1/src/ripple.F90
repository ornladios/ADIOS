subroutine get_rippl(r_in,z_in,ripp,drippdr,drippdz) 
!  use sml_module, only : sml_norm_r
  use rpl_module 
  implicit none
  real (kind=8), intent(in) :: r_in,z_in
  real (kind=8), intent(out) :: ripp,drippdr,drippdz
  real (kind=8) :: tau
  
  !real (kind=8)  :: rpl_N_coil, rpl_R0, rpl_tau0, rlp_ratio, rpl_elong
  !r0 : center of elipse
  !elong : elongation of elipse
  !tau0 : increase rate
  !ratio : relative strength to B at the center of elipse
  !ripple model
  ! ripple strenght = ratio*B*exp(tau/tau0)
  ! tau= (R-R0)**2 + elong*z**2
  ! ripp = ripple strength / B

  tau= (r_in - rpl_R0)**2 + rpl_elong*(z_in)**2
  ripp= rpl_ratio*exp(tau/rpl_tau0)
  
  drippdr=0D0 !incomplete
  drippdz=0D0 !incomplete


end subroutine get_rippl
