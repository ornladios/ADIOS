SUBROUTINE nse_test( jr_min, jr_max, rho, t, ye, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         nse_test
!    Module:       nse_test
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/28/00
!
!    Purpose:
!      To test the first mass zone with nse(j) = 0 to determine whether
!       it should be switched to nse. If so, the new temperature of the mass zone
!       is computed such that the internal energy of the zone after being switched
!       to nse is the same as before. The index nse(j) is then switched from 0 to 1.
!
!        xn(jz,1)     : Carbon-12    mass fraction
!        xn(jz,2)     : Oxygen-16    mass fraction
!        xn(jz,3)     : Neon-20      mass fraction
!        xn(jz,4)     : Magnesium-24 mass fraction
!        xn(jz,5)     : Silicon-28   mass fraction
!        xn(jz,6)     : Nickel-56    mass fraction
!        xn(jz,7)     : Alpha        mass fraction
!        xn(jz,8)     : Free neutron mass fraction
!        xn(jz,9)     : Free proton  mass fraction
!        xn(jz,10)    : Heavy nucleons (given by EOS at the
!                        termination of NSE) mass fraction
!
!        nse(j)       : 0, nse is not assumed for mass zone j
!                     : 1, nse is assumed for mass zone j
!
!    Subprograms called:
!  flash_x    : flashes material to NSE
!  deflash_x  : deflashes material out of NSE
!  set_cube_j : resets EOS, nuclear rates, and neutrino opacities
!
!    Input arguments:
!  jr_min     : minimum radial zone for which thermodynamic
!                variables are to be evaluated
!  jr_max     : maximum radial zone for which thermodynamic
!                variables are to be evaluated
!  rho        : shifted matter density array (g/cm**3).
!  t          : shifted matter matter temperature array (K).
!  ye         : shifted matter matter electron fraction array.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module,
!  cycle_module, edit_module, eos_snc_x_module, nucbrn_module,
!  parallel_module, prb_cntl_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx
USE numerical_module, ONLY: zero, one, epsilon
USE physcnst_module, ONLY: rmu, cm3fm3

USE cycle_module, ONLY: ncycle
USE edit_module, ONLY: nprint, nlog
USE eos_snc_x_module, ONLY: nse, eosrho, aesv, nuc_number
USE nucbrn_module, ONLY: xn, a_nuc_rep, z_nuc_rep, a_name
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY: tnse, tdnse
USE t_cntrl_module, ONLY: dtnmh, dtnph

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min           ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max           ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray           ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho  ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t    ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ye   ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1)                :: eos_reset        ! EOS reset flag
CHARACTER (len=1)                :: v_const          ! parameter that is kept constant during flash

LOGICAL                          :: first = .true.

INTEGER                          :: j                ! radial zone on a processor index
INTEGER                          :: j_test           ! first nonnse zone
INTEGER                          :: j_tstm1          ! last nse boundary
INTEGER                          :: i_radial         ! radial zone index

INTEGER                          :: i                ! composition index
INTEGER                          :: n_nucp1          ! nuc_number + 1
INTEGER                          :: i_He             ! helium index
INTEGER                          :: i_neut           ! neutron index
INTEGER                          :: i_prot           ! proton index
INTEGER                          :: i_Ni             ! 56Ni index

REAL(KIND=double)                :: kfm              ! ( # nucleons/gram )( cm3/fm3 )
REAL(KIND=double)                :: brydns           ! nucleons/fm^3
REAL(KIND=double)                :: t_nonNSE         ! value of t before flashing
REAL(KIND=double)                :: s_nonNSE         ! value of s before flashing
REAL(KIND=double)                :: s_NSE            ! value of s after flashing\

REAL(KIND=double)                :: x_n              ! neutron mass fraction after deflashing
REAL(KIND=double)                :: x_p              ! proton mass fraction after deflashing
REAL(KIND=double)                :: x_He             ! heilum mass fraction after deflashing
REAL(KIND=double)                :: x_Ni             ! nickel mass fraction after deflashing
REAL(KIND=double)                :: x_A              ! auxiliary heavy nucleus mass fraction after deflashing
!REAL(KIND=double), PARAMETER     :: x_Fe = 0.85d0   ! maxximum allowed mass fraction of 56Fe on deflashing
REAL(KIND=double), PARAMETER     :: x_Fe = 0.50d0    ! maxximum allowed mass fraction of 56Fe on deflashing
REAL(KIND=double), PARAMETER     :: A_min = 4.0d1    ! minimum heavy A for deflashing
REAL(KIND=double)                :: ZA_min           ! minimum ratio of Z to A for deflashing
REAL(KIND=double)                :: ZA_Aux           ! charge to mass ratio of auxiliary nucleus

REAL(KIND=double), PARAMETER     :: dtreset = 1.d-06 ! time step after flashing
REAL(KIND=double), PARAMETER     :: rho_set = 1.d+10 ! used for determining time step after flashing

REAL(KIND=double), PARAMETER     :: rho_h = 2.d+08   ! upper density for interpolating tnse
REAL(KIND=double), PARAMETER     :: rho_l = 5.d+07   ! lower density for interpolating tnse
REAL(KIND=double), PARAMETER     :: t_h   = 9.5d+9   ! upper temperature for interpolating tnse
REAL(KIND=double), PARAMETER     :: t_l   = 8.7d+9   ! lower temperature for interpolating tnse
REAL(KIND=double)                :: slope            ! lower temperature for interpolating tnse

REAL(KIND=double), PARAMETER     :: delta_t = 2.d+09 ! deflashing temperature 2.d+09

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Radial zone ',i3,' radial ray',i3,' has been switched to nse at cycle', &
& i7,' told=',1pe10.3,' tnew=',1pe10.3,' s_nonNSE=',es11.3,' s_NSE=',es11.3)
  201 FORMAT (' Radial zone ',i3,' radial ray',i3,' has been deflashed at cycle',i7, &
& ' xn =',1pe10.3,' xp=',1pe10.3,' xalpha =',1pe10.3,' x_Ni=',1pe10.3,' x_A=',1pe10.3, &
& ' a=',1pe10.3,' z=',1pe10.3,' s_NSE=',es11.3,' s_nonNSE=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constant
!-----------------------------------------------------------------------

IF ( first ) THEN
  first              = .false.

  kfm                = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
  slope              = ( t_h - t_l )/( rho_h - rho_l )
  ZA_min             = ( 26.d0/56.d0 ) * x_Fe + 0.5d0 * ( 1.d0 - x_Fe )

!........indexing

  i_neut           = 0
  i_prot           = 0
  i_He             = 0
  i_Ni             = 0
  n_nucp1          = nuc_number + 1

  DO i = 1,nuc_number
    IF (      a_name(i) == '  n  ' ) THEN
      i_neut       = i
    ELSE IF ( a_name(i) == '  p  ' ) THEN
      i_prot       = i
    ELSE IF ( a_name(i) == '  4He' ) THEN
      i_He         = i
    ELSE IF ( a_name(i) == ' 56Ni' ) THEN
      i_Ni         = i
    END IF
  END DO ! i

END IF ! first

!-----------------------------------------------------------------------
!
!             \\\\\ FIND NSE - NONNSE BOUNDARY /////
!
!  In the direction of increaseing zone indeces, j_test is the first
!   nonnse zone index. Return if there is no nse-nonnse interface.
!
!-----------------------------------------------------------------------

j_test               = 0

DO j = jr_min,jr_max
  IF ( nse(j,ij_ray,ik_ray) == 0 ) THEN
    j_test           = j
    EXIT
  END IF
END DO

IF ( j_test == 0 ) RETURN

!-----------------------------------------------------------------------
!
!             \\\\\ TEST ZONE FOR FLASHING /////
!
!  Test zone witn index j_test for switching to nse.
!
!-----------------------------------------------------------------------

tnse                 = 1.d+20
IF ( t(j_test) >= tnse  .and.  brydns >= eosrho ) THEN

!-----------------------------------------------------------------------
!  Store non-nse temperature and entropy
!-----------------------------------------------------------------------

  t_nonNSE           = t(j_test)
  s_nonNSE           = aesv(j_test,3,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Flash zone
!-----------------------------------------------------------------------

  v_const            = 'p'
  CALL flash_x( j_test, ij_ray, ik_ray, v_const )

!-----------------------------------------------------------------------
!  Reset thermodynamic quantities and neutrino rates
!-----------------------------------------------------------------------

  nse(j_test,ij_ray,ik_ray)  = 1
  CALL set_cube_j( j_test, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Reset the time step
!-----------------------------------------------------------------------

  dtnmh = DMIN1( dtreset, dtreset * rho_set/rho(j_test) )
  dtnph = DMIN1( dtreset, dtreset * rho_set/rho(j_test) )

!-----------------------------------------------------------------------
!  Record the event
!-----------------------------------------------------------------------

  s_NSE              = aesv(j_test,3,ij_ray,ik_ray)
  i_radial           = myid * i_ray_dim + i_ray
  WRITE (nprint,101) j_test,i_radial,ncycle,t_nonNSE,t(j_test),s_nonNSE,s_NSE
  WRITE (nlog,101) j_test,i_radial,ncycle,t_nonNSE,t(j_test),s_nonNSE,s_NSE
  RETURN

END IF

!-----------------------------------------------------------------------
!
!             \\\\\ TEST ZONE FOR DEFLASHING /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Test zone with index j_test - 1 for switching to nonnse.
!  Return if j_test = 1 since all zones are in non-NSE
!-----------------------------------------------------------------------

j_tstm1              = j_test - 1
IF ( j_tstm1 == 1 ) RETURN

!-----------------------------------------------------------------------
!  Deflash at t_flash - delta_t if the heavy nucleus, split inti 56Ni
!   and 56Fe, will result in less than a mass fraction x_Fe if 56Fe;
!   otherwise deflash when t goes below tdnse
!-----------------------------------------------------------------------

ZA_Aux               = aesv(j_tstm1,11,ij_ray,ik_ray)                            &
&                    / ( aesv(j_tstm1,10,ij_ray,ik_ray) + epsilon )
IF (      ( t(j_tstm1) < tdnse )                              .or.               &
&       ( ( ZA_Aux > ZA_min )  .and.  ZA_Aux <= 0.5d0                            &
&                              .and.  aesv(j_tstm1,10,ij_ray,ik_ray) >= A_min    &
&                              .and.  ( t(j_tstm1) < tnse - delta_t )            &
&                              .and.  ( brydns < 1.1d0 * eosrho     ) )         ) THEN

!-----------------------------------------------------------------------
!  Deflash zone
!-----------------------------------------------------------------------

  s_NSE              = aesv(j_tstm1,3,ij_ray,ik_ray)

  eos_reset          = 'y'
  CALL deflash_x( j_tstm1, rho, t, ye, ij_ray, ik_ray, eos_reset )

!-----------------------------------------------------------------------
!  Reset thermodynamic quantities and neutrino rates
!-----------------------------------------------------------------------

  nse(j_tstm1,ij_ray,ik_ray) = 0
  CALL set_cube_j( j_tstm1, rho, t, ye, ij_ray, ik_ray )
  s_nonNSE           = aesv(j_tstm1,3,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Record the event
!-----------------------------------------------------------------------

  x_n                = zero
  x_p                = zero
  x_He               = zero
  x_Ni               = zero
  x_A                = zero
  
  IF ( i_neut /= 0 ) x_n  = xn(j_tstm1,i_neut)
  IF ( i_prot /= 0 ) x_p  = xn(j_tstm1,i_prot)
  IF ( i_He   /= 0 ) x_He = xn(j_tstm1,i_He  )
  IF ( i_Ni   /= 0 ) x_Ni = xn(j_tstm1,i_Ni  )
  x_A                = xn(j_tstm1,n_nucp1)

  i_radial           = myid * i_ray_dim + i_ray
  WRITE (nprint,201) j_tstm1,i_radial,ncycle,x_n,x_p,x_He,x_Ni,x_A, &
& a_nuc_rep(j_tstm1),z_nuc_rep(j_tstm1),s_NSE,s_nonNSE
  WRITE (nlog     ,201) j_tstm1,i_radial,ncycle,x_n,x_p,x_He,x_Ni,x_A, &
& a_nuc_rep(j_tstm1),z_nuc_rep(j_tstm1),s_NSE,s_nonNSE

  RETURN

END IF ! flash_x or deflash criteria

RETURN
END SUBROUTINE nse_test
