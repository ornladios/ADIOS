!***********************************************************************
!
! physcnst_module
!
!***********************************************************************

MODULE physcnst_module

USE kind_module
SAVE


REAL(KIND=double), PARAMETER  :: Gw        = 1.027d-5      ! weak interaction constant (proton mass).

REAL(KIND=double), PARAMETER  :: Gw_MeV    = 1.16639d-11   ! weak interaction constant (MeV).

REAL(KIND=double), PARAMETER  :: sin2W     = 0.2315d+00    ! weak mixing angle.

REAL(KIND=double), PARAMETER  :: gv        = 1.0d+00       ! Fermi beta decay constant. This is unity as
!                                                            the coupling for Fermi decay is absorbed 
!                                                            by the weak interaction constant.

REAL(KIND=double), PARAMETER  :: ga        = 1.21d+00      ! Gamow-Teller beta decay constant. This is 
!                                                            not unity as the oupling for Gamow-Teller
!                                                            decay is renormalized. 

REAL(KIND=double), PARAMETER  :: cv        = 0.96d+00      ! weak interaction constant, which in the 
!                                                            standard model is given by
!                                                                    1       2
!                                                               cv = _ + 2sin (W)
!                                                                    2
!                                                            where W is the 'Weinberg angle'.  

REAL(KIND=double), PARAMETER  :: ca        = 0.5d+00       ! weak interaction constant, which in the 
!                                                            standard model is given by
!                                                               cv = 1/2 .

REAL(KIND=double), PARAMETER  :: mu_p      = 1.793d0       ! proton magnetic moment

REAL(KIND=double), PARAMETER  :: mu_n      = -1.913d0      ! neutron magnetic moment

REAL(KIND=double), PARAMETER  :: cvel      = 2.9979d+10    ! velocity of light (cm/s)

REAL(KIND=double), PARAMETER  :: pi        = 3.1415926535897932385d0 ! pi

REAL(KIND=double), PARAMETER  :: h         = 4.13567d-21   ! Planck's constant (MeV s)

REAL(KIND=double), PARAMETER  :: hbar      = h/( 2.d0*pi ) ! Planck's constant divided by 2*pi (MeV s)

REAL(KIND=double), PARAMETER  :: hbarc     = 197.33d+00    ! hbar * c (GeV fm)

REAL(KIND=double), PARAMETER  :: g         = 6.672d-08     ! Gravitational constant (dynes cm**2/gm**2)

REAL(KIND=double), PARAMETER  :: ergmev    = 1.602d-6      ! ergs per MeV

REAL(KIND=double), PARAMETER  :: ergfoe    = 1.d-51        ! ergs/s per foe

REAL(KIND=double), PARAMETER  :: kmev      = 8.61739d-11   ! Boltzmann's constant (MeV/K)

REAL(KIND=double), PARAMETER  :: kmev_inv  = 1.d0/kmev     ! Boltzmann's constant inverse (K/MeV)

REAL(KIND=double), PARAMETER  :: cm3fm3    = 1.d-39        ! cm**3/fm**3

REAL(KIND=double), PARAMETER  :: rmu       = 1.674d-24     ! mean baryon mass (gm)

REAL(KIND=double), PARAMETER  :: me        = 0.511d+00     ! electron mass (MeV)

REAL(KIND=double), PARAMETER  :: mp        = 938.27d+00    ! proton mass (MeV)

REAL(KIND=double), PARAMETER  :: mn        = 939.55d+00    ! neutron mass (MeV)

REAL(KIND=double), PARAMETER  :: mb        = 938.91d+00    ! baryon mass (MeV)

REAL(KIND=double), PARAMETER  :: mp_ex     = 7.2889705d+00 ! proton mass excess (MeV)

REAL(KIND=double), PARAMETER  :: mn_ex     = 8.0713171d+00 ! neutron mass excess (MeV)

REAL(KIND=double), PARAMETER  :: dmnp      = 1.2936d+00    ! mn - mp

REAL(KIND=double), PARAMETER  :: msolar    = 1.989d+33     ! mass of the sun (gm)

REAL(KIND=double), PARAMETER  :: gamhv     = 2.5d+00       ! mass of the sun (gm)

REAL(KIND=double), PARAMETER  :: wnm       = - 16.0d+00    ! mass of the sun (gm)

REAL(KIND=double), PARAMETER  :: ws        = 31.5d+00      ! mass of the sun (gm)

REAL(KIND=double), PARAMETER  :: xk0       = 180.0d+00     ! mass of the sun (gm)

REAL(KIND=double), PARAMETER  :: xkzafac   = 2.0d+00       ! mass of the sun (gm)

END MODULE physcnst_module

