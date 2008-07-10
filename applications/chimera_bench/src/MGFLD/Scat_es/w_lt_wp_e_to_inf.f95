SUBROUTINE w_lt_wp_e_to_inf( xl, eta, h0i, h0ii, h1i, h1ii )
!-----------------------------------------------------------------------
!
!    File:         w_lt_wp_e_to_inf
!    Module:       w_lt_wp_e_to_inf
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
!  nlag     : number of Gauss-Laguerre integration points
!  xl       : lower limit of integration
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

INTEGER                            :: i             ! summation index
INTEGER, PARAMETER                 :: nlag = 24     ! number of points of Gauss-Laguerre quadrature

REAL(KIND=double), DIMENSION(nlag) :: xa          ! points of Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nlag) :: wta         ! weights of Gauss-Lagendre quadrature

REAL(KIND=double)                  :: su0i          ! partial integral
REAL(KIND=double)                  :: su0ii         ! partial integral
REAL(KIND=double)                  :: su1i          ! partial integral
REAL(KIND=double)                  :: su1ii         ! partial integral

REAL(KIND=double)                  :: a0            ! integrand factor for computing the zero moment
REAL(KIND=double)                  :: b0            ! integrand factor for computing the zero moment
REAL(KIND=double)                  :: c0            ! integrand factor for computing the zero moment
REAL(KIND=double)                  :: u0            ! integrand factor for computing the zero moment
REAL(KIND=double)                  :: v0            ! integrand factor for computing the zero moment
REAL(KIND=double)                  :: a1            ! integrand factor for computing the first moment
REAL(KIND=double)                  :: b1            ! integrand factor for computing the first moment
REAL(KIND=double)                  :: c1            ! integrand factor for computing the first moment
REAL(KIND=double)                  :: u1            ! integrand factor for computing the first moment
REAL(KIND=double)                  :: v1            ! integrand factor for computing the first moment
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
  CALL glaquad( nlag, xa, wta, nlag )
  DO i = 1,nlag
    wta(i) = fexp(xa(i)) * wta(i)
  END DO
END IF

!-----------------------------------------------------------------------
!  Initialize for integration
!-----------------------------------------------------------------------

su0i               = zero
su0ii              = zero
su1i               = zero
su1ii              = zero

!-----------------------------------------------------------------------
!  Integrate
!-----------------------------------------------------------------------

DO i = 1,nlag

  e                = xa(i) + xl

  a0               =  r4_15 * w5
  b0               =  r4_3 * w4
  c0               =  r8_3 * w3
  u0               =  r8_3 * wp2 * w3 - 4.d0 * wp * w4 + r8_5 * w5
  v0               = -r16_3 * wp * w3 + 4.d0 * w4

  a1               =  r36_105 * w6/wp - r16_30 * w5
  b1               =  r8_5 * w5/wp    - r28_15 * w4
  c1               =  r12_5 * w4/wp   - r4_3 * w3
  u1               = -r4_3 * wp2 * w3 + r16_5 * wp * w4 - r16_5 * w5 + r8_7 * w6/wp
  v1               =  r8_3 * wp * w3  - r28_5 * w4 + r16_5 * w5/wp

  e2               = e * e
  h0               = a0 + b0 * e + c0 * e2
  hp               = u0 + v0 * e + c0 * e2
  res              = a1 + b1 * e + c1 * e2
  rep              = u1 + v1 * e + c1 * e2

  ep               = e + w_wp
  xme              = eta - e
  emxd             = ep - eta
  ffc              = ff(xme,emxd)

  su0i             = su0i  + h0  * ffc * wta(i)
  su0ii            = su0ii + hp  * ffc * wta(i)
  su1i             = su1i  + res * ffc * wta(i)
  su1ii            = su1ii + rep * ffc * wta(i)

END DO

h0i                = h0i  + su0i
h0ii               = h0ii + su0ii
h1i                = h1i  + su1i
h1ii               = h1ii + su1ii

END SUBROUTINE w_lt_wp_e_to_inf
