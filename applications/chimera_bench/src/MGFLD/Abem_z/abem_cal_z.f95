SUBROUTINE abem_cal_z( n, e_in, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
& cmpe, absornp, emitnp, ye_cube )
!-----------------------------------------------------------------------
!
!    File:         abem_cal_z
!    Module:       abem_cal_z
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/07
!
!    Purpose:
!      To call calculate the neutrino absorption and emission rates.
!      (1) No recoil or thermal motions: Subroutine abemfrnpetanp is called if
!       ye is <= 0.5. If ye > 0.5, the free neutron and proton mass fractions
!       are adjusted to reflect ye > 0.5 (the LS EOS cannot handle ye > 0.5),
!       and subroutine abemfrnp is called to compute the absorption and
!       emission rates.
!      (2) Recoil and thermal motions Included: Subroutine nu_N_absr_momts is
!       called. If ye > 0.5, the free neutron and proton mass fractions
!       are adjusted to reflect ye > 0.5, and neutron and proton chemical
!       potentials are computed.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  abemfrnpetanp_y : calculates rates of neutrino absorption and emission
!                     when chemical potential data are available
!  abemfrnp_y      : calculates rates of neutrino absorption and emission
!                     when chemical potential data are not available
!
!    Input arguments:
!  n               : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  e_in            : neutrino energy (MeV)!       
!  rho             : matter density (g/cm**3)
!  t               : matter temperature (K)
!  xneut           : free neutron mass fraction
!  xprot           : free proton mass fraction
!  xh              : heavy nucleus mass fraction
!  ah              : heavy nucleus mass number
!  zh              : heavy nucleus charge number
!  cmpn            : free neutron chemical potential (excluding rest mass) (MeV)
!  cmpp            : free proton chemical potential (excluding rest mass) (MeV)
!  cmpe            : electron chemical potential (including rest mass) (MeV)
!  iaefnp          : 0; inverse mean free paths set to zero
!                    1; inverse mean free paths computed
!  ye_cube         : electron fraction at the cube corner
!
!    Output arguments:
!  absor           : absorption inverse mean free path (/cm)
!  emitnp          : emission inverse mean free path (/cm)!  
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module,
!  abem_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, epsilon
USE physcnst_module, ONLY : cvel, hbar, gv, ga, mn, mp, me, Gw, pi, kmev, &
& dmnp, rmu

USE abem_z_module
USE prb_cntl_module, ONLY : iaefnp, i_aeps, rhoaefnp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)              :: n             ! neutrino flavor index

REAL(KIND=double), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye_cube       ! electron fraction at the cube corner
REAL(KIND=double), INTENT(in)    :: xneut         ! free neutron mass fraction
REAL(KIND=double), INTENT(in)    :: xprot         ! free proton mass fraction
REAL(KIND=double), INTENT(in)    :: xh            ! heavy nuclei mass fraction
REAL(KIND=double), INTENT(in)    :: ah            ! heavy nuclei mass number
REAL(KIND=double), INTENT(in)    :: zh            ! heavy nuclei charge number
REAL(KIND=double), INTENT(in)    :: cmpn          ! neutron chemical porential
REAL(KIND=double), INTENT(in)    :: cmpp          ! proton chemical porential
REAL(KIND=double), INTENT(in)    :: cmpe          ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout) :: absornp       ! inverse mean free path for absorption on free nucleons
REAL(KIND=double), INTENT(inout) :: emitnp        ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: i_abemetanp

REAL(KIND=double)                :: xneutp        ! free neutron mass fraction
REAL(KIND=double)                :: xprotp        ! free proton mass fraction
REAL(KIND=double)                :: xhe           ! helium mass fraction
REAL(KIND=double)                :: yep           ! electron fraction

REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: m_trgt_i      ! mass of the initial target particle (MeV)
REAL(KIND=double)                :: m_trgt_f      ! mass of the final target particle (MeV)
REAL(KIND=double)                :: m_lep         ! mass of the final lepton (MeV)
REAL(KIND=double)                :: cmp_trgt_i    ! chemical potential of the initial target particle (MeV)
REAL(KIND=double)                :: cmp_trgt_f    ! chemical potential of the transformed target particle (MeV)
REAL(KIND=double)                :: cmp_lep       ! chemcal potential of the secondary lepton (MeV)
REAL(KIND=double)                :: ab_r0_nu      ! zero moment of he inverse mean free path for neutrino absorption on free neutrons
REAL(KIND=double)                :: ab_r1_nu      ! first moment of he inverse mean free path for neutrino absorption on free neutrons
REAL(KIND=double)                :: e_out_e       ! mean energy of the emitted electron
REAL(KIND=double)                :: ab_r0_nub     ! zero moment of he inverse mean free path for antineutrino absorption on free protons
REAL(KIND=double)                :: ab_r1_nub     ! first moment of he inverse mean free path for antineutrino absorption on free protons
REAL(KIND=double)                :: e_out_p       ! mean energy of the emitted positron

REAL(KIND=double)                :: xi_n_wm       ! weak magnetism correction for antineutrino absorption on neutrons
REAL(KIND=double)                :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering

REAL(KIND=double), PARAMETER     :: x_min = 1.d-30 ! minimum mass fraction fraction

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ADJUST XNEUT AND XPROT IF NECESSARY ////
!
!  Adjust xneut and xprot if ye > yep or if ye < 0.03 to cover regions
!   out of range of the EOS.
!
!-----------------------------------------------------------------------

xneutp             = xneut
xprotp             = xprot
i_abemetanp        = .true.

IF ( ye_cube > 0.5d0 ) THEN

  yep              = DMIN1( ye_cube, 0.99d0 )
  xhe              = DMAX1( one - xneut - xprot - xh, zero )
  i_abemetanp      = .false.

!-----------------------------------------------------------------------
!
!  Are there enough nucleons to shuffle between neutrons and protons
!   to reproduce yep?
!  If so, adjust xneut and xprot to reproduce yep.
!
!-----------------------------------------------------------------------

  IF ( xneut + xprot + 0.5d0 * xhe + xh * zh/( ah + epsilon ) > yep ) THEN

    xprotp         = DMAX1( yep - 0.5d0 * xhe - xh * zh/( ah + epsilon ), x_min )
    xneutp         = DMAX1( xneut + xprot - xprotp, x_min )

  ELSE ! if not, do the best possible

    xprotp         = xneut + xprot
    xneutp         = 1.d-30

  END IF ! xneut + xprot + 0.5d0 * xhe + xh * zh/( ah + epsilon ) > yep
  
END IF ! ye_cube > 0.5d0


IF ( ye_cube < 0.03d0 ) THEN
  xprotp           = ye_cube
  xneutp           = one - ye_cube
  i_abemetanp      = .false.
END IF ! ye_cube < 0.03d0

!-----------------------------------------------------------------------
!
!             \\\\\ COMPUTE EMISSION AND ABSORPTION RATES ////
!
!-----------------------------------------------------------------------

IF ( i_aeps == 0 ) THEN
  IF ( i_abemetanp ) THEN

!-----------------------------------------------------------------------
!  Neutrino and antineutrino absorption on neutrons and protons with
!   approximate nucleon blocking but no recoil or thermal motions
!-----------------------------------------------------------------------

    CALL abemfrnpetanp( n, e_in, rho, t, xneutp, xprotp, cmpn, cmpp, cmpe, &
&    absornp, emitnp )
  ELSE

!-----------------------------------------------------------------------
!  Neutrino and antineutrino absorption on neutrons and protons with
!   no recoil, thermal motions, or nucleon blocking
!-----------------------------------------------------------------------

    CALL abemfrnp( n, e_in, rho, t, xneutp, xprotp, cmpe, absornp, emitnp )
  END IF ! i_abemetanp
ELSE ! i_aeps /= 0

!-----------------------------------------------------------------------
!  Neutrino absorption on neutrons with recoil, thermal motions, and
!   nucleon blocking
!-----------------------------------------------------------------------

  IF ( n == 1 ) THEN
    tmev           = kmev * t
    m_trgt_i       = mn
    m_trgt_f       = mp
    m_lep          = me
    cmp_trgt_i     = cmpn + dmnp + mn
    cmp_trgt_f     = cmpp + dmnp + mp
    cmp_lep        = cmpe
    CALL nu_N_absr_momts( n, e_in, tmev, m_trgt_i, m_trgt_f, m_lep, &
&    cmp_trgt_i, cmp_trgt_f, cmp_lep, ab_r0_nu, ab_r1_nu, e_out_e )
    absornp        = ab_r0_nu
!    WRITE (6,3001) n,m_trgt_i,m_trgt_f,m_lep,cmp_trgt_i,cmp_trgt_f,cmp_lep,e_in,e_out_e
! 3001 FORMAT (' n=',i4,' m_trgt_i=',es11.3,' m_trgt_f=',es11.3,' m_lep=',es11.3, &
!& ' cmp_trgt_i=',es11.3,' cmp_trgt_f=',es11.3,' cmp_lep=',es11.3,' e_in=',es11.3,' e_out_e=',es11.3)
  
!-----------------------------------------------------------------------
!  Antieutrino absorption on protons with recoil, thermal motions, and
!   nucleon blocking
!-----------------------------------------------------------------------

  ELSE IF ( n == 2 ) THEN
    tmev           = kmev * t
    m_trgt_i       = mp
    m_trgt_f       = mn
    m_lep          = me
    cmp_trgt_i     = cmpp + dmnp + mp
    cmp_trgt_f     = cmpn + dmnp + mn
    cmp_lep        = - cmpe
    CALL nu_N_absr_momts( n, e_in, tmev, m_trgt_i, m_trgt_f, m_lep, &
&    cmp_trgt_i, cmp_trgt_f, cmp_lep, ab_r0_nub, ab_r1_nub, e_out_p )
    absornp        = ab_r0_nub
!    WRITE (6,3001) n,m_trgt_i,m_trgt_f,m_lep,cmp_trgt_i,cmp_trgt_f,cmp_lep,e_in,e_out_p

!-----------------------------------------------------------------------
!  Absorption and emission on free nucleons is zero for mu and tau
!   neutrinos and antineutrinos
!-----------------------------------------------------------------------

  ELSE
    absornp        = zero
    emitnp         = zero
  END IF ! n == 1

END IF ! i_aeps == 0

!-----------------------------------------------------------------------
!  Weak magnetism corrections for neutrino and antineutrino absorptions
!   on nucleons.
!-----------------------------------------------------------------------

CALL cc_weak_mag( e_in, xi_n_wm, xib_p_wm )

IF ( n == 1 ) THEN
  absornp          = xi_n_wm * absornp
  emitnp           = xi_n_wm * emitnp
ELSE IF ( n == 2 ) THEN
  absornp          = xib_p_wm * absornp
  emitnp           = xib_p_wm * emitnp
END IF ! n == 1

RETURN
END SUBROUTINE abem_cal_z
