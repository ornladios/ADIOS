SUBROUTINE sctnncal( enuin, enuout, enubk, enubkp1, enubkp, enubkpp1, rho, &
& t, xn, xp, s_d )
!-----------------------------------------------------------------------
!
!    File:         sctnncal
!    Module:       sctnncal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/26/02
!
!    Purpose:
!      To compute the neutrino-nucleon inelastic scattering by the
!       analytic fitting formula given by Mannestad & Raffelt 1998, 
!       Apj, 507, 339.
!
!    Subprograms called:
!  g_sctnn
!  s_sctnn
!
!    Input arguments:
!  enuin     : incoming neutrino energy (zone centered) (MeV)
!  enuout    : outgoing neutrino energy (zone centered) (MeV)
!  enubk     : lower edge of incoming neutrino energy zone (MeV)
!  enubkp1   : upper edge of incoming neutrino energy zone (MeV)
!  enubkp    : lower edge of outgoing neutrino energy zone (MeV)
!  enubkpp1  : upper edge of outgoing neutrino energy zone (MeV)
!  rho       : matter density (g/cm3)
!  t         : matter temperature (K)
!  xn        : free neutron mass fraction
!  xp        : free proton mass fraction
!
!    Output arguments:
!  s_d       : down and isoenergetic neutrino nucleon inelastic
!                 scattering kernel
!
!    Include files: 
!        kind_module, numerical_module, physcnst_module
!
!----------------------------------------------------------------------c

USE kind_module
USE numerical_module
USE physcnst_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: enuin         ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: enuout        ! zone centered outgoing neutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: enubk         ! inner zone edge of incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: enubkp1       ! outer zone edge of incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: enubkp        ! inner zone edge of outgoing neutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: enubkpp1      ! outer zone edge of outgoing neutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: rho           ! density (g/cm3)
REAL(KIND=double), INTENT(in)  :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)  :: xn            ! free neutron mass fraction
REAL(KIND=double), INTENT(in)  :: xp            ! free proton mass fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out) :: s_d           ! down and isoenergetic scattering kernal

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                        :: first = .true.

INTEGER, PARAMETER             :: nleg_24 = 24  ! number of integration points
INTEGER                        :: i_eot         ! summation index for outgoing nu integration
INTEGER                        :: i_ein         ! summation index for incoming nu integration

REAL(KIND=double), PARAMETER   :: tthird = 2./3.

REAL(KIND=double)              :: C_AG_FnB      ! Raffelt coefficient
REAL(KIND=double)              :: conv
REAL(KIND=double)              :: coef
REAL(KIND=double)              :: tmev          ! temperature (meV)
REAL(KIND=double)              :: rho_14        ! rho/1e+14
REAL(KIND=double)              :: t_10          ! t/1e+10
REAL(KIND=double)              :: eta_star      ! degeneracy parameter
REAL(KIND=double)              :: gamma         ! spin fluctuation rate
REAL(KIND=double)              :: y             ! pion mass parameter
REAL(KIND=double)              :: sb            ! dimensionless fitting parameter
REAL(KIND=double)              :: gb            ! dimensionless fitting parameter
REAL(KIND=double)              :: xe_inu        ! upper limit of incoming nu integration
REAL(KIND=double)              :: xe_inl        ! lower limit of incoming nu integration
REAL(KIND=double)              :: a_ein         ! midpoint of incoming nu integration
REAL(KIND=double)              :: b_ein         ! halfwidth of incoming nu integration
REAL(KIND=double)              :: c_ein         ! halfwidth of incoming nu integration
REAL(KIND=double)              :: e_in          ! value of incoming nu integration variable
REAL(KIND=double)              :: xe_otu        ! upper limit of outgoing nu integration
REAL(KIND=double)              :: xe_otl        ! lower limit of outgoing nu integration
REAL(KIND=double)              :: a_eot         ! midpoint of outgoing nu integration
REAL(KIND=double)              :: b_eot         ! halfwidth of outgoing nu integration
REAL(KIND=double)              :: c_eot         ! halfwidth of outgoing nu integration
REAL(KIND=double)              :: e_ot          ! value of outgoing nu integration variable
REAL(KIND=double)              :: x_in_ot       ! value of outgoing nu integration variable
REAL(KIND=double)              :: fexp          ! declaration of function
REAL(KIND=double)              :: s_d_in        ! partial sum of scattering kernal
REAL(KIND=double)              :: s_d_ot        ! partial sum of scattering kernal

REAL(KIND=double), DIMENSION(24) :: xe_24
REAL(KIND=double), DIMENSION(24) :: wte_24

      xe_24 =    (/ 0.995187219997021d0  ,   0.974728555971309d0  ,   &
&                   0.938274552002733d0  ,   0.886415527004401d0  ,   &
&                   0.820001985973903d0  ,   0.740124191578554d0  ,   &
&                   0.648093651936976d0  ,   0.545421471388840d0  ,   &
&                   0.433793507626045d0  ,   0.315042679696163d0  ,   &
&                   0.191118867473616d0  ,   0.640568928626056d-1 ,   &
&                  -0.640568928626056d-1 ,  -0.191118867473616d0  ,   &
&                  -0.315042679696163d0  ,  -0.433793507626045d0  ,   &
&                  -0.545421471388840d0  ,  -0.648093651936976d0  ,   &
&                  -0.740124191578554d0  ,  -0.820001985973903d0  ,   &
&                  -0.886415527004401d0  ,  -0.938274552002733d0  ,   &
&                  -0.974728555971309d0  ,  -0.995187219997021d0  /)

      wte_24 =   (/ 0.123412297999872d-1 ,   0.285313886289337d-1 ,   &
&                   0.442774388174198d-1 ,   0.592985849154368d-1 ,   &
&                   0.733464814110803d-1 ,   0.861901615319533d-1 ,   &
&                   0.976186521041139d-1 ,   0.107444270115966d0  ,   &
&                   0.115505668053726d0  ,   0.121670472927803d0  ,   &
&                   0.125837456346828d0  ,   0.127938195346752d0  ,   &
&                   0.127938195346752d0  ,   0.125837456346828d0  ,   &
&                   0.121670472927803d0  ,   0.115505668053726d0  ,   &
&                   0.107444270115966d0  ,   0.976186521041139d-1 ,   &
&                   0.861901615319533d-1 ,   0.733464814110803d-1 ,   &
&                   0.592985849154368d-1 ,   0.442774388174198d-1 ,   &
&                   0.285313886289337d-1 ,   0.123412297999872d-1 /)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Initialize....................................................

IF (first) THEN
  C_AG_FnB     = 3.79d+4/cvel
  conv         = 3.d0 / ( 2.d0 * pi**2 )
  coef         = conv * C_AG_FnB
END IF

rho_14         = rho/1.e+14
tmev           = kmev * t
t_10           = t/1.e+10

!........Compute eta_star, the neutron effective degeneracy parameter..

eta_star       = 3.04d0 * rho_14**tthird / t_10

!........Compute gamma, the spin fluctuation rate......................

gamma          = 8.6d0 * rho_14 / sqrt(t_10)

!........Compute y, the pion mass parameter.

y              = 1.94 / t_10

!........Compute the dimensionless fitting parameter, gb...............

call g_sctnn(y,eta_star,gb)

!........Case 1: enuin = enuout........................................

IF ( enuin == enuout ) THEN

!.............Integrate from enubk to enuin

  xe_otu       = enuin
  xe_otl       = enubk
  s_d_ot       = zero
  a_eot        = half * ( xe_otu + xe_otl )
  b_eot        = half * ( xe_otu - xe_otl )

  DO i_eot = 1,nleg_24
    c_eot      = xe_24(i_eot) * b_eot
    e_ot       = a_eot + c_eot
    x_in_ot    = ( enuin - e_ot )/tmev
    CALL s_sctnn(x_in_ot,y,eta_star,sb)
    s_d_ot     = s_d_ot                                           &
&              + gamma/( x_in_ot**2 + ( half * gamma * gb )**2 )  &
&              * sb/tmev * e_ot**2 * wte_24(i_eot)

  END DO

  s_d          = one/( enuin**2 * ( enubkp1 - enubk ) )           &
&              * coef * rho_14 * ( xn + xp ) * s_d_ot * b_eot

!.............Integrate from enuin to enubkp1, compute upscattering from
!             downscattering * fexp( - x_in_ot )

  xe_otu       = enubkp1
  xe_otl       = enuin
  a_eot        = half * ( xe_otu + xe_otl )
  b_eot        = half * ( xe_otu - xe_otl )

  DO i_eot = 1,nleg_24
    c_eot      = xe_24(i_eot) * b_eot
    e_ot       = a_eot + c_eot
    x_in_ot    = - ( enuin - e_ot )/tmev
    CALL s_sctnn(x_in_ot,y,eta_star,sb)
    s_d_ot     = s_d_ot                                           &
&              + gamma/( x_in_ot**2 + ( half * gamma * gb )**2 )  &
&              * sb/tmev * e_ot**2 * wte_24(i_eot)                &
&              * fexp( - x_in_ot )
  END DO

   s_d         = s_d                                              &
&              + one/( enuin**2 * ( enubkp1 - enubk ) )           &
&              * coef * rho_14 * ( xn + xp ) * s_d_ot * b_eot

  RETURN

END IF ! enuin == enuout

!........Case 2: enubk = enubkpp1 (contiguous energy zones)............

IF ( enubk == enubkpp1 ) THEN

!.............Outer integral - integrate incoming neutrino energy from enubk to enubkp1

  s_d_in      = zero
  xe_inu      = enubkp1
  xe_inl      = enubk
  a_ein       = half * ( xe_inu + xe_inl )
  b_ein       = half * ( xe_inu - xe_inl )
  outer: DO i_ein = 1,nleg_24
    c_ein      = xe_24(i_ein) * b_ein
    e_in       = a_ein + c_ein

!.............Inner integral - integrate outgoing neutrino energy from enubkp to enubkpp1

    s_d_ot     = zero
    xe_otu     = enubkpp1
    xe_otl     = enubkp
    a_eot      = half * ( xe_otu + xe_otl )
    b_eot      = half * ( xe_otu - xe_otl )
    inner: DO i_eot = 1,nleg_24
      c_eot     = xe_24(i_eot) * b_eot
      e_ot      = a_eot + c_eot

!.............x: the dimensionless neutrino energy sum

      x_in_ot   = ( e_in - e_ot )/tmev

!.............Compute the dimensionless fitting parameter, sb

      CALL s_sctnn(x_in_ot,y,eta_star,sb)

!.............Differential kernel

      s_d_ot    = s_d_ot                                           &
&               + gamma/( x_in_ot**2 + ( half * gamma * gb )**2 )  &
&               * sb/tmev * e_ot**2 * wte_24(i_eot)

    END DO inner

    s_d_in      = s_d_in + s_d_ot * wte_24(i_ein) * e_in**3

  END DO outer

!.............FInal value

  s_d         = s_d_in * b_ein * b_eot * coef * rho_14 * ( xn + xp )  &
&             / ( enuin**3  * ( enubkp1  - enubk  )                   &
&             *   enuout**2 * ( enubkpp1 - enubkp ) )


  RETURN

END IF ! enubk = enubkpp1

!........Case 3: enubk > enubkpp1 (noncontiguous energy zones).........

x_in_ot      = ( enuin - enuout )/tmev
CALL s_sctnn(x_in_ot,y,eta_star,sb)
s_d          = gamma/( x_in_ot**2 + ( half * gamma * gb )**2 )  &
&            * sb/tmev * coef * rho_14 * ( xn + xp )

RETURN
END SUBROUTINE sctnncal
