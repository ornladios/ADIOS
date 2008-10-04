!=========================================================================================
  module reference_m
!=========================================================================================
! module for reference variables

  implicit none
!-----------------------------------------------------------------------------------------
! reals

  real a_ref        !reference velocity (speed of sound) (m/s)
  real l_ref        !reference length (m)
  real rho_ref      !reference density (kg/m^3)
  real lambda_ref   !reference conductivity (W/m-K)
  real mu_ref       !reference viscosity (Pa-s)
  real t_ref        !reference temperature (K)
  real t_o          !freestream temperature (K)
  real p_ref        !reference pressure (Pa)
  real time_ref     !reference time (s)
  real cp_ref       !reference specific heat at constant pres (J/kg-K)
  real univ_gascon  !universal gas constant (J/mol-K)
  real g_ref        !reference ratio specific heats
  real pres_atm     !standard atmospheric pressure from Chemkin (Pa)

  real re         !acoustic Reynolds number = rho_ref * a_ref * l_ref / mu_ref
  real re_real    !convective Reynolds number = rho_ref * u_ref * l_ref / mu_ref
  real mach_no    !Mach number = u_ref / a_ref (note that u_ref is never calculated)

!-----------------------------------------------------------------------------------------
  contains
!=========================================================================================
  subroutine initialize_reference(io,pr)
!=========================================================================================
! set reference values
!-----------------------------------------------------------------------------------------
  use topology_m, only : myid
  use chemkin_m, only : ickwrk, rckwrk

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed

  integer io
  real pr

! local declarations

  real dummy
!----------------------------------------------------------------------------------------
! calculate reference quantities

! set reynolds number
  re = re_real / mach_no

  univ_gascon = 8.314                     !J/mol-K
  t_ref = t_o * ( g_ref - 1 )             !K
  cp_ref = (a_ref**2) / t_ref             !J/kg-K
!jcs  mu_ref = pr * lambda_ref / cp_ref       !kg/m-s
!jcs  l_ref = re * mu_ref / rho_ref / a_ref   !m
!jcs 
  l_ref = pr * lambda_ref / cp_ref * re / rho_ref / a_ref   !m
  mu_ref = rho_ref*l_ref*a_ref       !kg/m-s

  p_ref = (a_ref**2)*rho_ref              !Pa
  time_ref=l_ref/a_ref                    !sec

  
!----------------------------------------------------------------------------------------
! get atmospheric pressure from chemkin

  call ckrp(ickwrk,rckwrk,dummy,dummy,pres_atm)   !dynes
  pres_atm=pres_atm/10.0                          !Pa
!----------------------------------------------------------------------------------------
! write all reference quantities to file and screen

  if(myid.eq.0) then

    write(io,*) 'initializing reference module...'
    write(io,*)

    write(io,*) 'the various reference values are as follows:'
    write(io,*)

    write(io,1) 'universal gas constant (J/mol-K)   = ',univ_gascon
    write(io,1) 'freestream temperature (K)         = ',t_o
    write(io,1) 'reference ratio of specifice heats = ',g_ref
    write(io,1) 'reference speed of sound (m/s)     = ',a_ref
    write(io,1) 'reference density (kg/m^3)         = ',rho_ref
    write(io,1) 'reference conductivity (W/m-K)     = ',lambda_ref
    write(io,1) 'reference temperature (K)          = ',t_ref
    write(io,1) 'reference pressure (Pa)            = ',p_ref
    write(io,1) 'standard atmospheric pressure (Pa) = ',pres_atm
    write(io,1) 'reference time (s)                 = ',time_ref
    write(io,1) 'reference specific heat (J/kg-K)   = ',cp_ref
    write(io,1) 'reference length (m)               = ',l_ref
    write(io,1) 'reference viscosity (Pa-s)         = ',mu_ref
    write(io,1) 'acoustic Reynolds number           = ',re
    write(io,1) 'Mach number                        = ',mach_no
    write(io,1) 'convective Reynolds number         = ',re*mach_no
    write(io,*)

  endif

  1 format(a40,1pe12.5)
  2 format(60('-'))
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine initialize_reference
!-----------------------------------------------------------------------------------------
  end module reference_m
