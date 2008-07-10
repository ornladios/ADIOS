SUBROUTINE paircal( enu, enubar, cmpe, t, j0i, j0ii, j1i, j1ii )
!-----------------------------------------------------------------------
!
!    File:         paircal
!    Module:       paircal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/5/03
!
!    Purpose:
!      To integrates the quantities
!
!          Fe(Ee)*Febar(Enu + Enubar - Ee)*Phi(Enu,Enubar,Ee)
!
!      and
!
!      [1-Fe(Ee)]*[1-Febar(Enu + Enubar - Ee)]*Phi(Enu,Enubar,Ee)
!
!  Fe                   : electron occupation number
!  Febar                : positron occupation number
!  Fe                   : electron energy
!  Enu                  : neutrino energy
!  Enubar               : antineutrino energy
!
!      Phi(Enu,Enubar,Ee) = cpair1*j0i + cpair2*j0ii
!
!    Subprograms called:
!  pr_w_gt_wp_e_lt_wp   : to integrate from 0 to wp for the case w > wp
!  pr_w_gt_wp_lt_e_lt_w : to integrate from wp to w for the case w > wp
!  pr_w_gt_wp_w_lt_e    : to integrate from w to w+wp for the case w > wp
!  pr_w_lt_wp_e_lt_w    : to integrate from 0 to w for the case w < wp
!  pr_w_lt_wp_e_lt_wp   : to integrate from w to wp for the case w < wp
!  pr_w_lt_wp_wp_lt_e   : to integrate from wp to w+wp for the case w < wp
!
!    Input arguments:
!  enu                  : neutrino energy (MeV)
!  enubar               : antineutrino energy (MeV)
!  cmpe                 : electron chemical potential (MeV)
!  t                    : matter temperature (K)
!
!    Output arguments:
!  j0i                  : zero moment of the "i" pair annihilation kernal
!  j0ii                 : zero moment of the "ii" pair annihilation kernal
!  j1i                  : first moment of the "i" pair annihilation kernal
!  j1ii                 : first moment of the "ii" pair annihilation kernal
!
!    Variables that must be passed through common:
!      none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  nes_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, one
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, pi, kmev

USE nes_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: enu           ! zone centered neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: enubar        ! zone centered antineutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: cmpe          ! electron chemical potential (MeV)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: j0i           ! zero moment of the "i" pair annihilation kernal
REAL(KIND=double), INTENT(out)   :: j0ii          ! zero moment of the "ii" pair annihilation kernal
REAL(KIND=double), INTENT(out)   :: j1i           ! first moment of the "i" pair annihilation kernal
REAL(KIND=double), INTENT(out)   :: j1ii          ! first moment of the "ii" pair annihilation kernal

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

REAL(KIND=double)                :: g2            ! combination of physical constants
REAL(KIND=double)                :: coef          ! combination of physical constants

REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: eta           ! cmpe/tmev
REAL(KIND=double)                :: beta          ! 1/tmev

REAL(KIND=double)                :: fexp          ! exponential

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize.
!
!  coef = (1/(2*pi)**2)*(1/!)*(1/(hb*!)**3)*g**2/pi
!
!  g2 has units of cm3/MeV2 s.
!  coef has units of 1/[ MeV5 cm ].
!  j0i, etc. have units of 1/[ MeV3 cm ].
!-----------------------------------------------------------------------

IF ( first ) THEN

  first            = .false.
  g2               = ( Gw/mp**2 )**2 * hbar**2 * cvel**3 ! cm3/MeV2 s
  coef             = 2.d0 * pi * ( 1.d0/cvel ) * ( 1.d0/( 2.d0 * pi * hbar * cvel ) )**3 * g2 / pi

END IF

w                  = enu
w2                 = w * w
w3                 = w2 * w
w4                 = w3 * w
w5                 = w4 * w
wp                 = enubar
wp2                = wp * wp
wp3                = wp2 * wp
wp4                = wp3 * wp
wp5                = wp4 * wp
wwp                = w * wp
w_p_wp             = w + wp
w_p_wp2            = w_p_wp * w_p_wp
w_p_wp3            = w_p_wp2 * w_p_wp
tmev               = t * kmev
beta               = one/tmev
eta                = beta * cmpe

j0i                = zero
j0ii               = zero
j1i                = zero
j1ii               = zero

!-----------------------------------------------------------------------
!  w > wp
!-----------------------------------------------------------------------

IF ( w >= wp ) THEN

  CALL pr_w_gt_wp_e_lt_wp( zero, wp, eta, beta, j0i, j0ii, j1i, j1ii )
  CALL pr_w_gt_wp_lt_e_lt_w( wp, w, eta, beta, j0i, j0ii, j1i, j1ii )
  CALL pr_w_gt_wp_w_lt_e( w, w_p_wp, eta, beta, j0i, j0ii, j1i, j1ii )
  
ELSE

!-----------------------------------------------------------------------
!  w < wp
!-----------------------------------------------------------------------

  CALL pr_w_lt_wp_e_lt_w( zero, w, eta, beta, j0i, j0ii, j1i, j1ii )
  CALL pr_w_lt_wp_e_lt_wp( w, wp, eta, beta, j0i, j0ii, j1i, j1ii )
  CALL pr_w_lt_wp_wp_lt_e( wp, w_p_wp, eta, beta, j0i, j0ii, j1i, j1ii )
  
END IF

!-----------------------------------------------------------------------
!  .Final touches
!-----------------------------------------------------------------------

j0i                = j0i /( enu * enu * enubar * enubar )
j0ii               = j0ii/( enu * enu * enubar * enubar )
j1i                = j1i /( enu * enu * enubar * enubar )
j1ii               = j1ii/( enu * enu * enubar * enubar )

j0i                = coef * j0i
j0ii               = coef * j0ii
j1i                = coef * j1i
j1ii               = coef * j1ii

RETURN
END SUBROUTINE paircal
