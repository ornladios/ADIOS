SUBROUTINE w_gt_wp_e_lt_wp( xl, xu, eta, h0i, h0ii, h1i, h1ii )
!-----------------------------------------------------------------------
!
!    File:         w_gt_wp_e_lt_wp
!    Module:       w_gt_wp_e_lt_wp
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/26/03
!
!    Purpose:
!      To integrate for the case w > w', e < wp, the quantities
!
!          f(e)*(1.-f(e+w-wp))*hl  (l = i,ii)
!
!      and
!
!          f(e+w-wp)*(1.-f(e))*hl  (l = i,ii)
!
!  e        : (electron energy)/kt    (integration variable)
!  w        : (in beam neutrino energy)/kt
!  wp       : (out beam neutrino energy)/kt
!  eta      : (electron chemical potential - mc2)/kt
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  nleg     : number of Gauss-Legendre integration points
!  xl       : lower limit of integration
!  xu       : upper limit of integration
!  eta      : (electron chemical potential - mc2)/kt
!
!    Output arguments:
!  h0i      : zero moment of the "i" in neutrino scattering function
!  h0ii     : zero moment of the "ii" in neutrino scattering function
!  h1i      : first moment of the "i" in neutrino scattering function
!  h1ii     : first moment of the "ii" in neutrino scattering function
!
!    Variables that must be passed through common:
!      none
!
!    Modules:
!  kind_module, numerical_module
!  nes_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero, half, one

USE nes_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)      :: xl            ! lower limit of integration
REAL(KIND=double), INTENT(in)      :: xu            ! upper limit of integration
REAL(KIND=double), INTENT(in)      :: eta           ! (electron chemical potential - mc2)/kt

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout)   :: h0i           ! partial integral contributed by this subroutine
REAL(KIND=double), INTENT(inout)   :: h0ii          ! partial integral contributed by this subroutine
REAL(KIND=double), INTENT(inout)   :: h1i           ! partial integral contributed by this subroutine
REAL(KIND=double), INTENT(inout)   :: h1ii          ! partial integral contributed by this subroutine

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                            :: first = .true.

INTEGER, PARAMETER                 :: nleg = 24     ! number of points of Gauss-Lagendre quadrature
INTEGER                            :: i             ! summation index

REAL(KIND=double), DIMENSION(nleg) :: xe          ! points of Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nleg) :: wte         ! weights of Gauss-Lagendre quadrature

REAL(KIND=double)                  :: su0i          ! partial integral
REAL(KIND=double)                  :: su0ii         ! partial integral
REAL(KIND=double)                  :: su1i          ! partial integral
REAL(KIND=double)                  :: su1ii         ! partial integral

REAL(KIND=double)                  :: e_mid         ! midpoint energy
REAL(KIND=double)                  :: e_del         ! half the energy width
REAL(KIND=double)                  :: e_var         ! scaled integration point

REAL(KIND=double)                  :: h0            ! integrand variable
REAL(KIND=double)                  :: hp            ! integrand variable
REAL(KIND=double)                  :: res           ! integrand variable
REAL(KIND=double)                  :: rep           ! integrand variable
REAL(KIND=double)                  :: ffc           ! integrand variable
REAL(KIND=double)                  :: ep            ! e + w - wp
REAL(KIND=double)                  :: xme           ! e - eta
REAL(KIND=double)                  :: emxd          ! eta - ep

REAL(KIND=double)                  :: fexp          ! exponential
REAL(KIND=double)                  :: ff            ! product of fermi distributions

EXTERNAL fexp
EXTERNAL ff

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Get quadrature points and weights
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  CALL gquad( nleg, xe, wte, nleg )
END IF

!-----------------------------------------------------------------------
!  Return if xl = xu
!-----------------------------------------------------------------------

IF ( xl == xu ) RETURN

!-----------------------------------------------------------------------
!  Initialize for integration
!-----------------------------------------------------------------------

su0i               = zero
su0ii              = zero
su1i               = zero
su1ii              = zero

e_mid              = half * ( xu + xl )
e_del              = half * ( xu - xl )

!-----------------------------------------------------------------------
!  Integrate
!-----------------------------------------------------------------------

DO i = 1,nleg

  e_var            = xe(i) * e_del
  e                = e_mid + e_var

  e2               = e * e
  e3               = e * e2
  e4               = e * e3
  e5               = e * e4
  e6               = e * e5
  e7               = e * e6

  h0               = r4_15 * e5 + r4_3 * e4 * w  + r8_3 * e3 * w2
  hp               = r4_15 * e5 - r4_3 * e4 * wp + r8_3 * e3 * wp2

  res              = w2pwp2 * h0 - r16_35 * e7 - r4_5 * e6 * ( 3.d0 * w - wp ) &
&      - r2_15 * e5 * ( 37.d0 * w2 - 26.d0 * wwp + wp2 ) - r2_3 * e4 * w * w_wp * ( 7.d0 * w - wp ) &
&      - r4_3 * e3 * w2 * w_wp2
  rep              =  w2pwp2 * hp - r16_35 * e7 + r4_5 * e6 * ( 3.d0 * wp - w ) &
&      - r2_15 * e5 * ( 37.d0 * wp2 - 26.d0 * wwp + w2 ) - r2_3 * e4 * wp * w_wp * ( 7.d0 * wp - w ) &
&      - r4_3 * e3 * wp2 * w_wp2

  res              = res/wwp
  rep              = rep/wwp

  ep               = e + w_wp
  xme              = e - eta
  emxd             = eta - ep
  ffc              = ff(xme,emxd)

  su0i             = su0i  + h0  * ffc * wte(i)
  su0ii            = su0ii + hp  * ffc * wte(i)
  su1i             = su1i  + res * ffc * wte(i)
  su1ii            = su1ii + rep * ffc * wte(i)

END DO

su0i               = e_del * su0i
su0ii              = e_del * su0ii
su1i               = e_del * su1i
su1ii              = e_del * su1ii

h0i                = h0i  + su0i
h0ii               = h0ii + su0ii
h1i                = h1i  + su1i
h1ii               = h1ii + su1ii

RETURN

END SUBROUTINE w_gt_wp_e_lt_wp
