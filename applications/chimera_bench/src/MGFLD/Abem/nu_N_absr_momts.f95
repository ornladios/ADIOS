SUBROUTINE nu_N_absr_momts( n, enu_in, tmev, m_trgt_i, m_trgt_f, m_lep, &
& cmp_trgt_i, cmp_trgt_f, cmp_lep, ab_r0, ab_r1, e_out )
!-----------------------------------------------------------------------
!
!    File:         nu_N_absr_momts
!    Module:       nu_N_absr_momts
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/25/02
!
!    Purpose:
!      To compute moments of the neutrino emission and absorption
!       inverse mean free paths.
!
!    Input arguments:
!   n          : neutrino type (1 = electron neutrino, 2 = electron
!                 antineutrino, 3 = muon or tau neutrino , 4 = muon or
!                 tau antineutrino
!   enu_in     : absorbed or emitted neutrino (MeV)
!   tmev       : temperature (MeV)
!   m_trgt_i   : mass of the target particle (MeV)
!   m_trgt_f   : mass of the transformed target particle (MeV)
!   m_lep      : mass of the created lepton (MeV)
!   cmp_trgt_i : neutron chemical potential (including rest mass) (MeV)
!   cmp_trgt_f : proton chemical potential (including rest mass) (MeV)
!   cmp_lep    : electron chemical potential (including rest mass) (MeV)
!
!    Output arguments:
!   ab_r0      : zero angular moment of the inverse absorption mean free
!                 path (cm^{-1})
!   ab_r1      : first angular moment of the inverse absorption mean
!                 free path (cm^{-1})
!   e_out      : mean energy of the emitted lepton (MeV)
!
!    The variables that must be passed through common are
!   cv        : weak interaction coefficient
!   ca        : weak interaction coefficient
!
!    Subprograms called:
!      cc_difcs
!
!    Include files:
!      dmn_parm.inc, physcnst.cmn
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : cvel, gv, ga

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: n             ! neutrino flavor index

REAL(KIND=double), INTENT(in)    :: enu_in        ! incoming neutrino energy
REAL(KIND=double), INTENT(in)    :: tmev          ! temperature (MeV)
REAL(KIND=double), INTENT(in)    :: m_trgt_i      ! mass of the initial target particle (MeV)
REAL(KIND=double), INTENT(in)    :: m_trgt_f      ! mass of the final target particle (MeV)
REAL(KIND=double), INTENT(in)    :: m_lep         ! mass of the final lepton (MeV)
REAL(KIND=double), INTENT(in)    :: cmp_trgt_i    ! chemical potential of the initial target particle (MeV)
REAL(KIND=double), INTENT(in)    :: cmp_trgt_f    ! chemical potential of the transformed target particle (MeV)
REAL(KIND=double), INTENT(in)    :: cmp_lep       ! chemcal potential of the secondary lepton (MeV)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: ab_r0         ! zero angular moment of the inverse absorption mean free path (cm^{-1})
REAL(KIND=double), INTENT(out)   :: ab_r1         ! first angular moment of the inverse absorption mean free path (cm^{-1})
REAL(KIND=double), INTENT(out)   :: e_out         ! mean energy of the emitted lepton

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

REAL(KIND=double)                :: m_cm          ! meters/cm
REAL(KIND=double)                :: m_trgt_i2     ! m_trgt_i^{2}
REAL(KIND=double)                :: m_trgt_f2     ! m_trgt_f^{2}
REAL(KIND=double)                :: m_lep2        ! m_lep/cm^{2}
REAL(KIND=double)                :: dm_trgt       ! m_trgt_i - m_trgt_f
REAL(KIND=double)                :: dm_trgtp      ! - dm_trgt

INTEGER                          :: i_a           ! summation index of angular Gauss-Lagendre quadrature
REAL(KIND=double)                :: xu_a          ! upper limit of lepton angular quadrature
REAL(KIND=double)                :: xl_a          ! lower limit of lepton angular quadrature
REAL(KIND=double)                :: mid_a         ! midpoint of lepton angular quadrature
REAL(KIND=double)                :: width_a       ! half-width of lepton angular quadrature
REAL(KIND=double)                :: c_a           ! scaled points of angular quadrature
REAL(KIND=double)                :: costh         ! points of angular quadrature

REAL(KIND=double)                :: ab_r0_e       ! zero moment of the neutrino absorption inverse mean free path per angle
REAL(KIND=double)                :: ab_r1_e       ! first moment of the neutrino absorption inverse mean free path per angle
REAL(KIND=double)                :: e_out_e       ! parameter to estimate energy quadrature widths per angle

REAL(KIND=double), PARAMETER     :: mult = 5.d0     ! boundary multiplier
REAL(KIND=double)                :: t_m             ! paremeter of the energy limits
REAL(KIND=double)                :: prin            ! paremeter of the energy limits
REAL(KIND=double)                :: radical         ! paremeter of the energy limits

INTEGER                          :: i_e           ! summation index of ouitgoing lepton energy Gauss-Lagendre quadrature
REAL(KIND=double)                :: xu_e          ! upper limit of lepton energy quadrature
REAL(KIND=double)                :: xl_e          ! lower limit of lepton energy quadrature
REAL(KIND=double)                :: mid_e         ! midpoint of lepton energy quadrature
REAL(KIND=double)                :: width_e       ! half-width of lepton energy quadrature
REAL(KIND=double)                :: c_e           ! scaled points of energy quadrature
REAL(KIND=double)                :: e_lep_f       ! energy quadrature points

REAL(KIND=double)                :: absr_eomega   ! absorption cross section per unit energy and solid angle

INTEGER, PARAMETER               :: nleg_a = 64   ! number of points of ouitgoing lepton angular Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION (nleg_a) :: x_a      ! scaled points of angular quadrature
REAL(KIND=double), DIMENSION (nleg_a) :: wt_a     ! scaled weights of angular quadrature

INTEGER, PARAMETER               :: nleg_e = 64   ! number of points of ouitgoing lepton energy Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION (nleg_e) :: x_e      ! scaled points of energy quadrature
REAL(KIND=double), DIMENSION (nleg_e) :: wt_e     ! scaled weights of energy quadrature

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN

!........Get integration points and weights

  CALL gquad(nleg_a,x_a,wt_a,nleg_a)
  CALL gquad(nleg_e,x_e,wt_e,nleg_e)

!........Once only

  first            = .false.

END IF ! first

!........Initialize constants

m_cm               = 1.d-2
m_trgt_i2          = m_trgt_i * m_trgt_i
m_trgt_f2          = m_trgt_f * m_trgt_f
m_lep2             = m_lep * m_lep
dm_trgt            = m_trgt_i - m_trgt_f
dm_trgtp           = - dm_trgt

!-----------------------------------------------------------------------
!  Initialize for double (angle and energy) integration
!-----------------------------------------------------------------------

ab_r0              = zero
ab_r1              = zero
e_out              = zero

!-----------------------------------------------------------------------
!
!           \\\\\ INTEGRATION OVER FINAL LEPTON ANGLE /////
!
!-----------------------------------------------------------------------

xu_a               = one
xl_a               = - one
mid_a              = half * ( xu_a + xl_a )
width_a            = half * ( xu_a - xl_a)

DO i_a = 1,nleg_a
  c_a              = x_a(i_a) * width_a
  costh            = mid_a + c_a

!-----------------------------------------------------------------------
!  Initialize for final lepton energy integration
!-----------------------------------------------------------------------

  ab_r0_e          = zero
  ab_r1_e          = zero
  e_out_e          = zero

!-----------------------------------------------------------------------
!
!          \\\\\ INTEGRATION OVER GFINAL LEPTON ENERGY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine the energy quadrature limits
!-----------------------------------------------------------------------

  t_m             = 2.d0 * mult * tmev/m_trgt_i

  prin            = one - t_m * costh
  radical         = DSQRT( 2.d0 * t_m * ( one - costh ) - t_m * t_m * ( one - costh * costh ) )
  xu_e            = enu_in * ( prin + radical )/( one - t_m ) + dm_trgt
  xl_e            = enu_in * ( prin - radical )/( one - t_m ) + dm_trgt
  xl_e            = DMAX1( xl_e, zero   )

!      IF ( i_a == 2   .or.  i_a == 32  .or.  i_a == 64  .or.  i_a == 16  .or.  i_a == 48 ) &
!& WRITE (6,3005) t_m,prin,radical,xu_e,xl_e
! 3005 format (' t_m=',1pe10.3,' prin=',1pe10.3,' radical=',1pe10.3,' xu_e=',1pe10.3,' xl_e=',1pe10.3)

!-----------------------------------------------------------------------
!  Integrate over final lepton energy
!-----------------------------------------------------------------------

  mid_e            = half * ( xu_e + xl_e )
  width_e          = half * ( xu_e - xl_e )

  DO i_e = 1,nleg_e
    c_e            = x_e(i_e) * width_e
    e_lep_f        = mid_e + c_e
    IF ( e_lep_f <= zero ) CYCLE
    CALL cc_difcs( enu_in, e_lep_f, costh, tmev, m_trgt_i, dm_trgtp, &
&    cmp_trgt_i, cmp_trgt_f, cmp_lep, gv, ga, absr_eomega )

!    IF ( i_a == 2  .or.  i_a == 32  .or.  i_a == 64  .or. i_a == 16  .or.  i_a == 48 ) &
!& WRITE (6,3001) enu_in,e_lep_f,costh,absr_eomega
! 3001 FORMAT (' enu_in=',1pe10.3,' e_lep_f=',1pe10.3,' costh=',1pe10.3,' absr_eomega=',1pe10.3)

    absr_eomega    = DMAX1( absr_eomega, zero )
    ab_r0_e        = ab_r0_e + absr_eomega * wt_e(i_e)
    ab_r1_e        = ab_r1_e + absr_eomega * wt_e(i_e) * costh
    e_out_e        = e_out_e + absr_eomega * wt_e(i_e) * e_lep_f
  END DO ! i_e = 1,nleg_e
 
  ab_r0_e          = ab_r0_e * width_e
  ab_r1_e          = ab_r1_e * width_e
  e_out_e          = e_out_e * width_e

!-----------------------------------------------------------------------
!  Sum final result over angle and average the final electron energy
!-----------------------------------------------------------------------

!  IF ( i_a == 2   .or.  i_a == 32  .or.  i_a == 64  .or. i_a == 16  .or.  i_a == 48 ) &
!& WRITE (6,3002) enu_in,e_lep_f,ab_r0_e,ab_r1_e
! 3002 format (' enu_in=',1pe10.3,' e_lep_f=',1pe10.3,' ab_r0_e=',1pe10.3,' ab_r1_e=',1pe10.3)

  ab_r0            = ab_r0 + ab_r0_e * wt_a(i_a)
  ab_r1            = ab_r1 + ab_r1_e * wt_a(i_a)
  e_out            = e_out + e_out_e * wt_a(i_a) 

END DO ! i_a = 1,nleg_a

e_out              = e_out/( ab_r0 + epsilon )
ab_r0              = m_cm * ab_r0 * width_a
ab_r1              = m_cm * ab_r1 * width_a

RETURN
END SUBROUTINE nu_N_absr_momts
