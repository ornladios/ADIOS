SUBROUTINE scatiicr( rho, t, e_nu, xion, aion, zion, ciicr )
!-----------------------------------------------------------------------
!
!    File:         scatiicr_H
!    Module:       scatiicr_H
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To compute correction factors for the isoenergetic
!         scattering rates for ion-ion correlation effects
!         using a prescription given by C. Horowitz.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!   rho       : matter density (g/cm**3)
!   t         : matter temperature (K)
!   e_nu         : neutrino energy
!   xion      : mass fraction of heavy ions
!   aion      : mass number of heavy ions
!   zion      : charge number of heavy ions
!
!    Output arguments:
!
!   ccicr     : correction factor for ion-ion correlation effects
!
!    Input arguments (common):
!
!   iicor     : 0, ion-ion correlation effects  omitted
!               not 0, ion-ion correlation effects  computed
!
!    Output arguments (common):
!        none
!
!    Include files:
!  numerical_module, physcnst_module
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, one, third
USE physcnst_module, ONLY : pi, kmev, hbarc, rmu

USE prb_cntl_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: e_nu          ! neutrino energy
REAL(KIND=double), INTENT(in)    :: xion          ! heavy nucleus mass fraction
REAL(KIND=double), INTENT(in)    :: aion          ! heavy nucleus mass number
REAL(KIND=double), INTENT(in)    :: zion          ! heavy nucleus charge number

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT)   :: ciicr         ! correlation correction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: i             ! do index

REAL(KIND=double), PARAMETER     :: cgamma = 1.96d-5 ! e**2/4*pi*rmu (e**2 = 4*pi*alpha (natural units) 
!                                                    ! = 4*pi*alpha*hbar*c (cgs units) = 4*pi*(4.8e-10)**2
REAL(KIND=double), PARAMETER     :: cx = 9.44d+9
REAL(KIND=double), PARAMETER     :: amin = 1.d-10

REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: gamma         ! zion**2*e**2/4*pi*a*kt = cgamma*zion**2/tmev*(rho*xion/aion)**1/3
REAL(KIND=double)                :: gamma12       ! gamma^{1/2}
REAL(KIND=double)                :: nion          ! number density of ions
REAL(KIND=double)                :: a             ! (3/4*pi*n)**1/3 = (3*rmu*aion/4*pi*rho*xion)**1/3
REAL(KIND=double)                :: hbarc_cm      ! hbar (Mev cm)
REAL(KIND=double)                :: aE            ! a/hbarc_cm
REAL(KIND=double)                :: ebar          ! e_nu * aE - Wigner cell size/neutrino wavelength
REAL(KIND=double)                :: estar         ! 3.d0 + 4.d0/gamma12
REAL(KIND=double)                :: arg           ! argument of exponential

REAL(KIND=double)                :: beta0
REAL(KIND=double), DIMENSION(6)  :: beta
REAL(KIND=double), DIMENSION(4,4) :: beta3

REAL(KIND=double)                :: fexp          ! exponential function

EXTERNAL fexp


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN

  first            = .false.

  beta3(1,1)       = -7.362056d+0
  beta3(1,2)       =  0.5371365d0
  beta3(1,3)       = -0.1078845d0
  beta3(1,4)       =  4.189612d-3

  beta3(2,1)       =  3.4489581d0
  beta3(2,2)       = -0.40251656d0
  beta3(2,3)       =  9.0877878d-2
  beta3(2,4)       = -3.4353581d-3

  beta3(3,1)       = -0.74128645d0
  beta3(3,2)       =  0.11019855d0
  beta3(3,3)       = -2.5359361d-2
  beta3(3,4)       =  9.0487744d-4

  beta3(4,1)       =  5.9573285d-2
  beta3(4,2)       = -1.0186552d-2
  beta3(4,3)       =  2.2791369d-3
  beta3(4,4)       = -7.4614597d-5
  
END IF


!-----------------------------------------------------------------------
!  Return if xion < amin  or  aion < amin  or  zion < amin or iicor = 0.
!-----------------------------------------------------------------------

IF ( DMIN1( xion, aion, zion ) < amin  .or.  iicor == 0 ) THEN
  ciicr            = one
  RETURN
END IF

!-----------------------------------------------------------------------
!  Compute gamma, a, and ebar
!-----------------------------------------------------------------------

tmev               = kmev * t
gamma              = cgamma * ( zion**2/tmev ) * ( rho * xion/aion )**third
gamma              = DMAX1( gamma, one )
gamma              = DMIN1( gamma, 150.d0 )
gamma12            = DSQRT(gamma)
nion               = rho * xion/( rmu * aion )
a                  = ( 3.d0/( 4.d0 * pi * nion ) )**third
hbarc_cm           = 1.d-13 * hbarc
aE                 = a/hbarc_cm
ebar               = e_nu * aE

!-----------------------------------------------------------------------
!  ciicr = 1 if ebar > estar
!-----------------------------------------------------------------------

estar              = 3.d0 + 4.d0/gamma12
IF ( ebar >= estar ) THEN
  ciicr            = one
  RETURN
END IF

!-----------------------------------------------------------------------
!  Compute ciicr if ebar > estar
!-----------------------------------------------------------------------

beta0              = DLOG(0.3d0/( 0.3d0 + 3.d0 * gamma ) )
beta(1)            = zero
beta(2)            = 6.667d0
DO i = 3,6
  beta(i)          = beta3(i-2,1) + beta3(i-2,2) * gamma12 + beta3(i-2,3) * gamma + beta3(i-2,4) &
&                  * gamma * gamma12
END DO

arg                = beta0
DO i = 1,6
  arg              = arg + beta(i) * ebar**(DBLE(i))
END DO

ciicr              = one/( one + fexp( -arg ) )

RETURN
END SUBROUTINE scatiicr
