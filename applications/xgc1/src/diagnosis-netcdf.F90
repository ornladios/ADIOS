
subroutine diag_flow_ncinit
  use netcdf
  use diag_module
  use eq_module, only : eq_x_psi
  use sml_module, only : sml_electron_on
  implicit none

  real (kind=8), dimension(diag_flow_npsi) :: psi1
  integer :: j
  ! netCDF error status return
  integer :: iret
  character (len=32) :: var_names(diag_flow_nvars2)
!!$  =(/ 'psi', 'density', 'toroidal_flow', 'poloidal_flow',&
!!$         'parallel_flow', 'poloidal_comp._of_ExB', 'Radial_flow_times_grad_psi', 'Radial_current_density', &
!!$       'v_para_times_B', 'parallel_mean_energy', 'perp._temperature', 'parallel_temperature', &
!!$       'density(df)','toroidal_flow(df)', 'poloidal_flow(df)',&
!!$       'parallel_flow(df)', 'poloidal_comp._of_ExB(df)', 'Radial_flow_times_grad_psi(df)', 'Radial_current_density(df)', &
!!$       'v_para_times_B(df)', 'parallel_mean_energy(df)', 'perp._temperature(df)', 'parallel_temperature(df)', &
!!$       'delta-f_ratio', 'flow_shell_volume' /)
  character (len=32) :: var_units(diag_flow_nvars2)
!!$  =(/ 'A.U.', 'm^-3', 'm/s', 'm/s', &
!!$       'm/s', 'm/s', 'A.U.', 'A.U.', &
!!$       'mT/s','keV', 'keV', 'keV', &
!!$       'm^-3', 'm/s', 'm/s', &
!!$       'm/s', 'm/s', 'A.U.', 'A.U.', &
!!$       'mT/s','keV', 'keV', 'keV', &
!!$       'Unitless', 'm^3' /)
  character (len=5) :: sp_name(2)=(/'ion__', 'elec_'/)
  integer :: sp_type
#endif 
  var_names(1) ='psi'
  var_names(2) ='density(df)'
  var_names(3) ='toroidal_flow(df)'
  var_names(4) ='poloidal_flow(df)'
  var_names(5) ='parallel_flow(df)'
  var_names(6) ='poloidal_comp._of_ExB(df)'
  var_names(7) ='Radial_flow_times_grad_psi(df)'
  var_names(8) ='Radial_current_density(df)'
  var_names(9) ='v_para_times_B(df)'
  var_names(10)='parallel_mean_energy(df)'
  var_names(11)='perp._temperature(df)'
  var_names(12)='parallel_temperature(df)'
  var_names(13)='density'
  var_names(14)='toroidal_flow'
  var_names(15)='poloidal_flow'
  var_names(16)='parallel_flow'
  var_names(17)='poloidal_comp._of_ExB'
  var_names(18)='Radial_flow_times_grad_psi'
  var_names(19)='Radial_current_density'
  var_names(20)='v_para_times_B'
  var_names(21)='parallel_mean_energy'
  var_names(22)='perp._temperature'
  var_names(23)='parallel_temperature'
  var_names(24)='delta-f_ratio'
  var_names(25)='flow_shell_volume'

  var_units(1 )='A.U.'
  var_units(2 )='m^-3'
  var_units(3 )='m/s'
  var_units(4 )='m/s'
  var_units(5 )='m/s'
  var_units(6 )='m/s'
  var_units(7 )='A.U.'
  var_units(8 )='A.U.'
  var_units(9 )='mT/s'
  var_units(10)='keV'
  var_units(11)='keV'
  var_units(12)='keV'
  var_units(13)='m^-3'
  var_units(14)='m/s'
  var_units(15)='m/s'
  var_units(16)='m/s'
  var_units(17)='m/s'
  var_units(18)='A.U.'
  var_units(19)='A.U.'
  var_units(20)='mT/s'
  var_units(21)='keV'
  var_units(22)='keV'
  var_units(23)='keV'
  var_units(24)='Unitless'
  var_units(25)= 'm^3'

  ! create netCDF file and enter define mode
  iret = nf90_create('xgc.flowdiag.cdf',NF90_CLOBBER, diag_flow_ncfileid)
  call check_err(iret)

  ! define netCDF dimensions
  iret = nf90_def_dim(diag_flow_ncfileid, 'samples', diag_flow_npsi, diag_psi_ncdim)
  call check_err(iret)

  iret = nf90_def_dim(diag_flow_ncfileid, 'timesteps', NF90_UNLIMITED, diag_time_ncdim)
  call check_err(iret)

  ! define netCDF variables
  diag_psi_ncdims(1) = diag_psi_ncdim
  iret = nf90_def_var(diag_flow_ncfileid, 'psi', NF90_DOUBLE, diag_psi_ncdims, diag_psi_ncvarid)
  call check_err(iret)

  diag_time_ncdims(1) = diag_time_ncdim
  iret = nf90_def_var(diag_flow_ncfileid, 'time', NF90_DOUBLE, diag_time_ncdims, diag_time_ncvarid)
  call check_err(iret)
  iret = nf90_put_att(diag_flow_ncfileid, diag_time_ncvarid, 'units', 'transit times')
  call check_err(iret)

  diag_flow_ncdims(1) = diag_psi_ncdim
  diag_flow_ncdims(2) = diag_time_ncdim

     do j=2, diag_flow_nvars2
        iret = nf90_def_var(diag_flow_ncfileid,sp_name(sp_type)// var_names(j), NF90_DOUBLE, diag_flow_ncdims, diag_flow_ncvarid(j,sp_type))
        call check_err(iret)
        iret = nf90_put_att(diag_flow_ncfileid, diag_flow_ncvarid(j,sp_type), 'units', var_units(j))
        call check_err(iret)
     enddo
  enddo
  ! set netCDF attribute for ElVis to start monitoring netCDF file
  iret = nf90_put_att(diag_flow_ncfileid, NF90_GLOBAL, 'running', 'true')
  call check_err(iret)

  ! leave netCDF define mode
  iret = nf90_enddef(diag_flow_ncfileid)
  call check_err(iret)

  ! put data into netCDF psi variable
  do j=1, diag_flow_npsi
     psi1(j) = (diag_flow_pin + diag_flow_dp * (real(j)-0.5D0)) / eq_x_psi
  enddo
  iret = nf90_put_var(diag_flow_ncfileid, diag_psi_ncvarid, psi1)
  call check_err(iret)

  ! sync output file
  iret = nf90_sync(diag_flow_ncfileid)
  call check_err(iret)

  ! initialize program variables
  diag_time_nccount(1) = 1
  diag_flow_ncstart(1) = 1
  diag_flow_nccount(1) = diag_flow_npsi
  diag_flow_nccount(2) = 1
end subroutine diag_flow_ncinit

subroutine diag_flow_ncfinal
  use netcdf
  use diag_module
  implicit none

! netCDF error status return
  integer :: iret

! sete netCDF attribute for ElVis to stop monitoring netCDF file
  iret = nf90_redef(diag_flow_ncfileid)
  call check_err(iret)
  iret = nf90_put_att(diag_flow_ncfileid, NF90_GLOBAL, 'running', 'false')
  call check_err(iret)

! leave netCDF define mode and close file
  iret = nf90_enddef(diag_flow_ncfileid)
  call check_err(iret)
  iret = nf90_close(diag_flow_ncfileid)
  call check_err(iret)

end subroutine diag_flow_ncfinal

subroutine diag_avg_ncfinal
  use netcdf
  use diag_module
  implicit none

! netCDF error status return
  integer :: iret

! sete netCDF attribute for ElVis to stop monitoring netCDF file
  iret = nf90_redef(diag_avg_ncfileid)
  call check_err(iret)
  iret = nf90_put_att(diag_avg_ncfileid, NF90_GLOBAL, 'running', 'false')
  call check_err(iret)

! leave netCDF define mode and close file
  iret = nf90_enddef(diag_avg_ncfileid)
  call check_err(iret)
  iret = nf90_close(diag_avg_ncfileid)
  call check_err(iret)
  
end subroutine diag_avg_ncfinal

subroutine diag_avg_ncinit
  use netcdf
  use diag_module
  use eq_module, only : eq_x_psi
  use sml_module, only : sml_electron_on
  implicit none  
  real (kind=8), dimension(diag_flow_npsi) :: psi1
  integer, parameter :: nvars=diag_avg_nv1+diag_avg_nv2
  integer :: j
  ! netCDF error status return
  integer :: iret
  character (len=32) :: var_names(nvars) 
!!$  =(/ 'psi', 'radial_flux', 'radial_E_para_flux', 'radial_E_perp_flux',&
!!$       'parallel_flow', 'parallel_mean_E', 'perp_temperature', 'v_r_v_para_corel', &
!!$       'average_gradpsi_sqr', 'full_weight', 'density', 'total_E_flux', &
!!$       'grad_T','Thermal_conductivity(psi)', 'Thermal_conductivity(m)' /)
  character (len=32) :: var_units(nvars)
!!$  =(/ 'A.U.', 'psi/s', 'J psi/s', 'J psi/s', &
!!$       'm/s', 'eV', 'eV', 'psi m /s^2', &
!!$       'psi^2/m^2','Number', 'm^-3', 'J psi/s', &
!!$       'J/psi', 'psi^2/s', 'm^2/s' /)
  character (len=5) :: sp_name(2)=(/'ion__', 'elec_'/)
  integer :: sp_type

  var_names(1 )='psi'
  var_names(2 )='radial_flux(avg)'
  var_names(3 )='radial_E_para_flux(avg)'
  var_names(4 )='radial_E_perp_flux(avg)'
  var_names(5 )='parallel_flow(avg)'
  var_names(6 )='parallel_mean_E(avg)'
  var_names(7 )='perp_temperature(avg)'
  var_names(8 )='v_r_v_para_corel(avg)'
  var_names(9 )='average_gradpsi_sqr(avg)'
  var_names(10)='full_weight(avg)'
  var_names(11)='density(avg)'
  var_names(12)='total_E_flux(avg)'
  var_names(13)='grad_T(avg)'
  var_names(14)='Thermal_conductivity(psi:avg)'
  var_names(15)='Thermal_conductivity(m:avg)'

  var_units(1)= 'A.U.'
  var_units(2)= 'psi/s'
  var_units(3)= 'J psi/s'
  var_units(4)= 'J psi/s'
  var_units(5)= 'm/s'
  var_units(6)= 'eV'
  var_units(7)= 'eV'
  var_units(8)= 'psi m /s^2'
  var_units(9)= 'psi^2/m^2'
  var_units(10)= 'Number'
  var_units(11)= 'm^-3'
  var_units(12)= 'J psi/s'
  var_units(13)= 'J/psi'
  var_units(14)= 'psi^2/s'
  var_units(15)= 'm^2/s'
  

! create netCDF file and enter define mode
  iret = nf90_create('xgc.fluxdiag.cdf',NF90_CLOBBER, diag_avg_ncfileid)
  call check_err(iret)

! define netCDF dimensions
  iret = nf90_def_dim(diag_avg_ncfileid, 'samples', diag_avg_npsi, diag_avg_psi_ncdim)
  call check_err(iret)

  iret = nf90_def_dim(diag_avg_ncfileid, 'timesteps', NF90_UNLIMITED, diag_avg_time_ncdim)
  call check_err(iret)

! define netCDF variables
  diag_avg_psi_ncdims(1) = diag_avg_psi_ncdim
  iret = nf90_def_var(diag_avg_ncfileid, 'psi', NF90_DOUBLE, diag_avg_psi_ncdims, diag_avg_psi_ncvarid)
  call check_err(iret)

  diag_avg_time_ncdims(1) = diag_avg_time_ncdim
  iret = nf90_def_var(diag_avg_ncfileid, 'time', NF90_DOUBLE, diag_avg_time_ncdims, diag_avg_time_ncvarid)
  call check_err(iret)
  iret = nf90_put_att(diag_avg_ncfileid, diag_avg_time_ncvarid, 'units', 'transit times')
  call check_err(iret)

  diag_avg_ncdims(1) = diag_avg_psi_ncdim
  diag_avg_ncdims(2) = diag_avg_time_ncdim

  do sp_type=1, 1+sml_electron_on
     do j=2, nvars
        iret = nf90_def_var(diag_avg_ncfileid,sp_name(sp_type)// var_names(j), NF90_DOUBLE, diag_avg_ncdims, diag_avg_ncvarid(j,sp_type))
        call check_err(iret)
        iret = nf90_put_att(diag_avg_ncfileid, diag_avg_ncvarid(j,sp_type), 'units', var_units(j))
        call check_err(iret)
     enddo
  enddo
! set netCDF attribute for ElVis to start monitoring netCDF file
  iret = nf90_put_att(diag_avg_ncfileid, NF90_GLOBAL, 'running', 'true')
  call check_err(iret)

! leave netCDF define mode
  iret = nf90_enddef(diag_avg_ncfileid)
  call check_err(iret)

! put data into netCDF psi variable
  do j=1, diag_avg_npsi
     psi1(j) = (diag_avg_pin + diag_avg_dp * (real(j)-0.5D0)) / eq_x_psi
  enddo
  iret = nf90_put_var(diag_avg_ncfileid, diag_avg_psi_ncvarid, psi1)
  call check_err(iret)

! sync output file
  iret = nf90_sync(diag_avg_ncfileid)
  call check_err(iret)

! initialize program variables
  diag_avg_time_nccount(1) = 1
  diag_avg_ncstart(1) = 1
  diag_avg_nccount(1) = diag_avg_npsi
  diag_avg_nccount(2) = 1

end subroutine diag_avg_ncinit
