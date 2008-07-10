

!-----------------------------------------------------------------------
!         Neutrino energy advection controls
!
!      ivc_x: neutrino energy advection switch.
!
!          ivc_x = 0: neutrino advection bypassed.
!          ivc_x = 1: neutrino advection included.
!
!      ivc_y: neutrino energy advection due to y-hydro switch.
!
!          ivc_y = 0: y-hydro neutrino energy advection bypassed.
!          ivc_y = 1: y-hydro neutrino energy advection included.
!
!      tau_advct: optical depth above which neutrinos are advected in
!       energy like a gamma = 4/3 gass
!
!      rhomin_y_eadvect: density below which y-neutrino energy advection
!       is turned off
!
!-----------------------------------------------------------------------

ivc                          1                                          ivc_x
ivc                          1                                          ivc_y
tau                                 5.00000000e+00                      tau_advct
rhomin                              1.00000000e+12                      rhomin_y_eadvect

!-----------------------------------------------------------------------
!         Time step control criteria
!
!      t_cntrl_burn(n): n-neutrino zero-moment change time step
!       criterion due to advection, i.e., maximum permitted 
!       abs( dpsi0(j,k,n)/psi0(j,k,n) + psivmin(n) ) due to neutrino
!       energy advection.
!
!      psivmin(n): parameters used in determining the psi0 change time
!       step due to neutrino energy advection.
!
!      dpsivmx(n) = max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psivmin(n) ) ) .
!
!-----------------------------------------------------------------------

tcntrl                       1      1.00000000E-01                      t_cntrl_e_advct
tcntrl                       2      1.00000000E-01                      t_cntrl_e_advct
tcntrl                       3      1.00000000E-01                      t_cntrl_e_advct
tcntrl                       4      1.00000000E-01                      t_cntrl_e_advct
psivtl                       1      1.00000000e-01                      psivmin
psivtl                       2      1.00000000e-01                      psivmin
psivtl                       3      1.00000000e-01                      psivmin
psivtl                       4      1.00000000e-01                      psivmin
