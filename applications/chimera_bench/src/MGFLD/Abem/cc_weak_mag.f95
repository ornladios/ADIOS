SUBROUTINE cc_weak_mag( e_nu, xi_n_wm, xib_p_wm )
!-----------------------------------------------------------------------
!
!    File:         cc_weak_mag
!    Module:       cc_weak_mag
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/05/07
!
!    Purpose:
!      To compute the weak magnetism correction for neutrino absorption
!       on neutrons and antineutrino absorption on protons.
!
!    Subprograms called:
!       none
!
!    Input arguments:
!  e_nu           : neutrino energy (MeV)
!
!    Output arguments:
!  xi_n_wm        : weak magnetism correction for neutrino absorption on neutrons
!  xib_p_wm       : weak magnetism correction for antineutrino absorption on protons
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module
!  prb_cntl_module
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: half, one, third
USE physcnst_module, ONLY: sin2W, mu_p, mu_n, mp, mn, ga

USE prb_cntl_module, ONLY: icc_wm

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)        :: e_nu          ! neutrino energy (MeV)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)       :: xi_n_wm       ! weak magnetism correction for neutrino absorption on neutrons
REAL(KIND=double), INTENT(out)       :: xib_p_wm      ! weak magnetism correction for antineutrino absorption on protons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                              :: first = .true.

REAL(KIND=double)                    :: CV            ! charged vector current coupling constant
REAL(KIND=double)                    :: CA            ! charged axial vector current coupling constant
REAL(KIND=double)                    :: F2            ! neutron-proton form factor
REAL(KIND=double)                    :: CV_2          ! CV * CV
REAL(KIND=double)                    :: CA_2          ! CA * CA
REAL(KIND=double)                    :: CV_F2CA       ! (CV + F2)*CA
REAL(KIND=double)                    :: CVF2          ! CV * F2
REAL(KIND=double)                    :: F22           ! F2 * F2
REAL(KIND=double)                    :: CV3CA         ! ( CV_2 + 3.d0 * CA_2 )
REAL(KIND=double)                    :: m             ! 1/2 * ( mp + mn )

REAL(KIND=double)                    :: e             ! e_nu/mc2
REAL(KIND=double)                    :: e2            ! e**2
REAL(KIND=double)                    :: zeta          ! 1 + 2e
REAL(KIND=double)                    :: zeta3         ! zeta**3
REAL(KIND=double)                    :: ln_zeta       ! ln(zeta)

REAL(KIND=double)                    :: chi_wm_rec    ! cross section correction due to weak magnetism and recoil
REAL(KIND=double)                    :: chi_wm_rec1   ! cross section correction 1 due to weak magnetism and recoil
REAL(KIND=double)                    :: chi_wm_rec2   ! cross section correction 2 due to weak magnetism and recoil
REAL(KIND=double)                    :: chi_wm_rec3   ! cross section correction 3 due to weak magnetism and recoil
REAL(KIND=double)                    :: chi_wm_rec4   ! cross section correction 4 due to weak magnetism and recoil
REAL(KIND=double)                    :: chi_wm_rec5   ! cross section correction 5 due to weak magnetism and recoil
REAL(KIND=double)                    :: chi_rec       ! cross section correction due to recoil
REAL(KIND=double)                    :: chi_rec1      ! cross section correction 1 due to recoil
REAL(KIND=double)                    :: chi_rec2      ! cross section correction 2 due to recoil


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set weak magnetism corrections to zero if inc_wm = 0
!-----------------------------------------------------------------------

IF ( icc_wm == 0 ) THEN
  xi_n_wm          = one
  xib_p_wm         = one
  RETURN
END IF ! icc_wm == 0

!-----------------------------------------------------------------------
!  Initialize coupling constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  CV               = 1.d0
  CA               = 1.26
  CV_2             = CV * CV
  CA_2             = CA * CA
  F2               = 3.706d0
  CV_F2CA          = ( CV + F2 ) * CA
  CVF2             = CV * F2
  F22              = F2 * F2
  CV3CA            = ( CV_2 + 3.d0 * CA_2 )
  m                = 0.5d0 * ( mn + mp )
  first            = .false.
END IF ! first

!-----------------------------------------------------------------------
!
!         \\\\\ GENERAL NEUTRINO-NUCLEON WEAK MAGNETISM /////
!         \\\\\       AND RECOIL CORRECTION TEMRS       /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Dimensionless energies
!-----------------------------------------------------------------------

e                  = e_nu/m
e2                 = e * e
zeta               = 1.d0 + 2.d0 * e
zeta3              = zeta * zeta * zeta
ln_zeta            = DLOG( zeta )

!-----------------------------------------------------------------------
!  chi_wm_rec terms
!-----------------------------------------------------------------------

chi_wm_rec1        = 1.d0 + 4.d0 * e + 16.d0 * e2 * third
chi_wm_rec2        = 3.d0 * ( 1.d0 + 4.d0 * e * third ) * ( 1.d0 + 4.d0 * e * third )
chi_wm_rec3        = 4.d0 * e * ( 1.d0 + 4.d0 * e * third )
chi_wm_rec4        = 8.d0 * third * e2
chi_wm_rec5        = 5.d0 * third * e2 * ( 1.d0 + 2.d0 * e/5.d0 )

!-----------------------------------------------------------------------
!
!     \\\\\ GENERAL NEUTRINO-NUCLEON RECOIL CORRECTION TEMRS /////
!
!-----------------------------------------------------------------------

chi_rec1           = 1.d0/e - ( 1.d0/( 2.d0 * e2 ) ) * ln_zeta
chi_rec2           = ( zeta * ln_zeta - 2.d0 * e + 4.d0 * e2 )/( 2.d0 * e2 * zeta )

!-----------------------------------------------------------------------
!
!       \\\\\ NEUTRINO ABSORPTION WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

chi_wm_rec         = ( CV_2 * chi_wm_rec1 + CA_2 * chi_wm_rec2        &
&                  + CV_F2CA * chi_wm_rec3 + CVF2 * chi_wm_rec4     &
&                  + F22 * chi_wm_rec5 )/( zeta3 * CV3CA )
chi_rec            = ( CV_2 * chi_rec1 + CA_2 * chi_rec2 )/CV3CA
xi_n_wm            = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!     \\\\\ ANTINEUTRINO ABSORPTION WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

chi_wm_rec         = ( CV_2 * chi_wm_rec1 + CA_2 * chi_wm_rec2        &
&                  - CV_F2CA * chi_wm_rec3 + CVF2 * chi_wm_rec4     &
&                  + F22 * chi_wm_rec5 )/( zeta3 * CV3CA )
chi_rec            = ( CV_2 * chi_rec1 + CA_2 * chi_rec2 )/CV3CA
xib_p_wm           = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE cc_weak_mag
