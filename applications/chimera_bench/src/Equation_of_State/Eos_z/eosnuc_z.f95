SUBROUTINE eosnuc_z
!-----------------------------------------------------------------------
!
!    File:         eosnuc_z
!    Module:       eosnuc_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To compute the nuclear component of the equation of state
!       for material not in nuclear statistical equilibrium. The
!       composition of the material is given by the array xn(k,i).
!
!      The nucleons and nuclei are assumed to comprise an ideal
!       Boltzmann gas.
!
!    Variables that must be passed through common:
!  xn(k,i) : mass fraction of the ith nucleus in radial zone k.
!  Cooperstein-BCK eos variables.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, eos_bck_module, eos_snc_x_module, eos_snc_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnc
USE numerical_module, ONLY : zero, third, one, epsilon
USE physcnst_module, ONLY : mb, pi, hbarc, wnm, ws, dmnp

USE edit_module, ONLY : nlog
USE eos_bck_module, ONLY : upack, dbck, u1pack, zabck, jshel, d0, un, uhat, &
& dtran, xnbck, xpbck, xabck, xhbck, ahbck, b, eh, ph, sh, tbck, ed, pd, &
& therm, sd
USE eos_snc_x_module, ONLY : a_name
USE eos_snc_z_module, ONLY : xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, nuc_number, &
& a_nuc, z_nuc, be_nuc

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                           :: first = .true.

INTEGER, PARAMETER                :: ncd = 300      ! composition dimension
INTEGER                           :: i              ! composition index
INTEGER                           :: ii             ! composition index
INTEGER                           :: n_nucp1        ! nuc_number + 1
INTEGER                           :: n_hvy          ! number of heavy nuclei (not counting representative heavy nucleus)
INTEGER                           :: i_He           ! helium index
INTEGER                           :: i_neut         ! neutron index
INTEGER                           :: i_prot         ! proton index
INTEGER, DIMENSION(ncd)           :: i_hvy          ! heavy nucleus index

REAL(KIND=double), DIMENSION(ncd) :: thetab         ! mass number of nucleus i
REAL(KIND=double), DIMENSION(ncd) :: phib           ! mass number of nucleus i
REAL(KIND=double), DIMENSION(ncd) :: bcompb         ! mass number of nucleus i
REAL(KIND=double), DIMENSION(ncd) :: xmstarb        ! mass number of nucleus i
REAL(KIND=double), DIMENSION(ncd) :: afb            ! mass number of nucleus i

REAL(KIND=double)                 :: a_harmonic     ! harmonic mean heavy nuclear mass number
REAL(KIND=double)                 :: a_harmonic_inv ! inverse of the harmonic mean heavy nuclear mass number
REAL(KIND=double)                 :: zhbck          ! neutrino absorption and emission time step
REAL(KIND=double)                 :: av_inv         ! inverse of the harmonic mean nuclear mass number
REAL(KIND=double)                 :: zv             ! mean nuclear charge number

REAL(KIND=double), PARAMETER      :: gprot = 2.d0   ! proton spin degeneracy
REAL(KIND=double), PARAMETER      :: gnewt = 2.d0   ! neutron spin degeneracy
REAL(KIND=double), PARAMETER      :: galph = 1.d0   ! alpha particle spin degeneracy
REAL(KIND=double), PARAMETER      :: ghvy  = 1.d0   ! heavy nucleus spin degeneracy

REAL(KIND=double)                 :: c0             ! constant for computing the chemical potential
REAL(KIND=double)                 :: a0             ! constant for computing the nuclear level density
REAL(KIND=double)                 :: xms            ! effective nuclear mass at high density
REAL(KIND=double)                 :: xm0            ! effective nuclear mass at low density
REAL(KIND=double)                 :: xm0ms          ! xm0 - xms
REAL(KIND=double)                 :: tc             ! nuclear level density parameter
REAL(KIND=double)                 :: a13            ! cube root of the nuclear mass number
REAL(KIND=double)                 :: piby2          ! constant

REAL(KIND=double)                 :: x_pheavy       ! mass fraction of bound protons
REAL(KIND=double)                 :: x_nuclei       ! mass fraction of nuclei
REAL(KIND=double)                 :: x_nucleon      ! mass fraction of free nucleons
REAL(KIND=double)                 :: xalph          ! mass fraction of helium
REAL(KIND=double)                 :: x_heavy        ! mass fraction of heavy nuclei

REAL(KIND=double)                 :: av             ! harmonic mean nuclear mass number
REAL(KIND=double)                 :: rn             ! number of nuclei per baryon
REAL(KIND=double)                 :: sn             ! entropy of free neutrons
REAL(KIND=double)                 :: sp             ! entropy of free protons
REAL(KIND=double)                 :: sa             ! entropy of alphas
REAL(KIND=double)                 :: strh           ! entropy of heavy nuclei
REAL(KIND=double)                 :: up             ! proton chemical potential
REAL(KIND=double)                 :: ua             ! alpha chemical potential
REAL(KIND=double)                 :: x_nuc_i        ! mass fraction of nucleus i
REAL(KIND=double)                 :: a_i            ! mass number of nucleus i
REAL(KIND=double)                 :: utrhvy_i       ! chemical potential of nucleus i
REAL(KIND=double)                 :: strh_i         ! translational entropy of nucleus i

!-----------------------------------------------------------------------
!        Formats.
!-----------------------------------------------------------------------

  101 FORMAT (' nnc =',i4,' > ncd=',i4,' in subroutine eosnuc_z.')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.

!........dimension consistency

  IF ( nnc > ncd ) THEN
    WRITE (nlog,101) nnc,ncd
    STOP
  END IF

!........Constants

  c0               = ( mb/( 2.d0 * pi * hbarc**2 ) )**( 3.d0/2.d0 )
  a0               = .067d0
  piby2            = 1.570796327d0
  xms              = 0.7d0
  xm0              = 2.0d0
  xm0ms            = xm0 - xms
  d0               = 0.16d0
  upack            = dbck/d0
  u1pack           = 1.0d0 - upack
  tc               = 12.d0
  n_nucp1          = nuc_number + 1

!........indexing

  n_hvy            = 0
  i_He             = 0
  i_neut           = 0
  i_prot           = 0
  i_hvy            = 0
  ii               = 0

  DO i = 1,nuc_number
    IF (      a_name(i) == '  n  ' ) THEN
      i_neut       = i
    ELSE IF ( a_name(i) == '  p  ' ) THEN
      i_prot       = i
    ELSE IF ( a_name(i) == '  4He' ) THEN
      i_He         = i
    ELSE
      ii           = ii + 1
      n_hvy        = n_hvy + 1
      i_hvy(ii)    = i
    END IF
  END DO

!........Nuclear parameters

  DO i = 1,n_hvy
    zabck          = DMIN1( z_nuc(i_hvy(i))/a_nuc(i_hvy(i)), 1.d0 )
    phib(i)        = fphi(zabck)
    thetab(i)      = 1.0d0
    a13            = f3(a_nuc(i_hvy(i)))
    bcompb(i)      = a_nuc(i_hvy(i)) * ( fbulk(zabck) + fsurf(zabck)/a13 + fcoul(zabck) * a13 * a13 )
    afb(i)         = a0 * faf(phib(i) * thetab(i))
  END DO

END IF ! first

!-----------------------------------------------------------------------
!  Load heavy nucleus parameter
!-----------------------------------------------------------------------

be_nuc(n_nucp1)    = be_nuc_rep(jshel)
a_nuc(n_nucp1)     = a_nuc_rep(jshel)
z_nuc(n_nucp1)     = z_nuc_rep(jshel)
zabck              = DMIN1( z_nuc(n_nucp1)/( a_nuc(n_nucp1) + epsilon ), 1.d0 )
phib(n_nucp1)      = fphi(zabck)
thetab(n_nucp1)    = 1.0d0
afb(n_nucp1)       = a0 * faf(phib(n_nucp1) * thetab(n_nucp1))
un                 = zero
uhat               = zero
dtran              = 0.16d0

!-----------------------------------------------------------------------
!             Composition.
!
!              Compute
!  xnbck      : free neutron mass fraction
!  xpbck      : free proton mass fraction
!  xabck      : helium mass fraction
!  xhbck      : heavy nuclei mass fraction excluding helium
!  zabck      : average charge to mass ratio for heavy nuclei
!  ahbck      : mass averaged heavy nucleus mass number
!  a_harmonic : harmonic average og heavy nucleus mass number
!
!  xnbck and xpbck are computed from xnbck + xbck and ye by the
!   requirement that xnbck + xpbck be conserved and that ye be the 
!   proton fraction.
!
!  x_pheavy    : mass fraction of bound protons including helium
!  x_nucleon   : mass fraction of free nucleons
!  x_nuclei    : mass fraction of nuclei plus helium
!  x_heavy     : mass fraction of heavy nuclei excluding helium and
!                 the auxiliary nucleus
!-----------------------------------------------------------------------

!........load proton, neutrn and helium mass fractions

xnbck              = zero
xpbck              = zero
xabck              = zero
IF ( i_neut /= 0 ) xnbck = xn(jshel,i_neut)
IF ( i_prot /= 0 ) xpbck = xn(jshel,i_prot)
IF ( i_He /= 0 )   xabck = xn(jshel,i_He)
xhbck              = zero
zabck              = zero
a_harmonic         = zero

!........load xhbck, x_pheavy, x_nuclei

x_pheavy           = zero
x_nuclei           = zero
DO i = 1,n_hvy
  xhbck            = xhbck    + xn(jshel,i_hvy(i))
  x_pheavy         = x_pheavy + z_nuc(i_hvy(i)) * xn(jshel,i_hvy(i))/a_nuc(i_hvy(i))
  x_nuclei         = x_nuclei + xn(jshel,i_hvy(i))
END DO
xhbck              = xhbck + xn(jshel,n_nucp1)

x_pheavy           = x_pheavy + z_nuc(i_He) * xn(jshel,i_He)/( a_nuc(i_He) + epsilon )
x_pheavy           = x_pheavy + z_nuc(n_nucp1) * xn(jshel,n_nucp1)/( a_nuc(n_nucp1) + epsilon )
x_nuclei           = x_nuclei + xn(jshel,i_He) + xn(jshel,n_nucp1)
x_nucleon          = DMAX1( one - x_nuclei, zero )

ahbck              = zero
x_heavy            = zero
DO i = 1,n_hvy
  x_heavy          = x_heavy + xn(jshel,i_hvy(i))
  ahbck            = ahbck + xn(jshel,i_hvy(i)) * a_nuc(i_hvy(i))
END DO
ahbck              = ahbck/( x_heavy + epsilon )

a_harmonic_inv     = zero
zhbck              = zero
DO i = 1,n_hvy
  a_harmonic_inv   = a_harmonic_inv +              xn(jshel,i_hvy(i))/( a_nuc(i_hvy(i)) + epsilon )
  zhbck            = zhbck     + z_nuc(i_hvy(i)) * xn(jshel,i_hvy(i))/( a_nuc(i_hvy(i)) + epsilon )
END DO
a_harmonic_inv     = a_harmonic_inv +                xn(jshel,n_nucp1)/( a_nuc_rep(jshel) + epsilon )
zhbck              = zhbck     + z_nuc_rep(jshel)  * xn(jshel,n_nucp1)/( a_nuc_rep(jshel) + epsilon )

a_harmonic         = xhbck/( a_harmonic_inv + epsilon )
zhbck              = zhbck * a_harmonic/( xhbck + epsilon )
zabck              = DMIN1( zhbck/( a_harmonic + epsilon ), 1.d0 )

!-----------------------------------------------------------------------
!  Update xpbck and xnbck from updated ye
!-----------------------------------------------------------------------

IF ( xhbck <= zero ) THEN
  zabck            = 0.5d0
  a_harmonic       = 56.d0
END IF

xalph              = 1.d0 - xnbck - xpbck - xhbck

!-----------------------------------------------------------------------
!        Compute
!  b  : binding energy/particle in unburnt shell
!  av : harmonic average mass number of all particles (appropriate
!        for computing the kinetic energy)
!  rn : number of nuclei per baryon
!-----------------------------------------------------------------------

av_inv             = zero
zv                 = zero
b                  = zero
DO i = 1,n_nucp1
  av_inv           = av_inv +             xn(jshel,i)/( a_nuc(i) + epsilon )
  zv               = zv     + z_nuc(i)  * xn(jshel,i)/( a_nuc(i) + epsilon )
  b                = b      - be_nuc(i) * xn(jshel,i)/( a_nuc(i) + epsilon )
END DO
av                 = one/( av_inv + epsilon )
zv                 = zv * av
rn                 = one/( av + epsilon )

!-----------------------------------------------------------------------
!        Compute
!  eh : energy of nuclei w/o translation
!  ph : pressure of nuclei w/o translation
!  sh : entropy of heavy nuclei (intrinsic)
!  ed : energy of drip (including heavy translation)
!  pd : pressure of drip (including heavy translation)
!  sd : entropy of drip (including heavy translation)
!
!  e0 = ( egy0 - ( ye - yeFe )*dmnp ) is the energy per baryon of an
!   isolated iron nucleus plus the rest mass difference per baryon
!   between matter at ye and that at the ye of iron. The value of eh is
!   relative to that energy.!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  eh, ph, and sh
!-----------------------------------------------------------------------

eh              = zero
ph              = zero
sh              = zero

DO i = 1,n_hvy
  xmstarb(i)    = xms + xm0ms /( one + tbck/tc )**2
  sh            = sh + xn(jshel,i_hvy(i)) * afb(i) * tbck * ( 2.d0 * xmstarb(i) - 2.d0 * tbck/tc &
&               / ( one + tbck/tc ) * ( xmstarb(i) - xms ) )
  eh            = eh + tbck * sh
END DO

xmstarb(n_nucp1) = xms + xm0ms /( one + tbck/tc )**2
sh              = sh + xn(jshel,n_nucp1) * afb(n_nucp1) * tbck * ( 2.d0 * xmstarb(n_nucp1) - 2.d0 * tbck/tc &
&               / ( one + tbck/tc ) * ( xmstarb(n_nucp1) - xms ) )
eh              = eh + tbck * sh

eh              = eh + b

!-----------------------------------------------------------------------
!  Drip: ed and pd
!-----------------------------------------------------------------------

ed              = 1.5d0 * rn * tbck
pd              = rn * dbck * tbck

therm           = c0 * tbck * dsqrt(tbck) / dbck

!-----------------------------------------------------------------------
!  Neutrons: un and sn (dmnp is subracted from un to be bring it in
!   conformaty with the LS and BCK EOS
!-----------------------------------------------------------------------

IF ( xnbck <= epsilon ) THEN
  un            = -100.d0
  sn            = zero
ELSE
  un            = tbck * DLOG( DABS( xnbck/( therm * gnewt ) ) )
  sn            = ( xnbck   ) * ( 2.5d0 - un/tbck )
  un            = un - dmnp
END IF ! xnbck < 0

!-----------------------------------------------------------------------
!  Protons: up and sp (dmnp is subracted from up to be bring it in
!   conformaty with the LS and BCK EOS.
!-----------------------------------------------------------------------

IF ( xpbck <= epsilon ) THEN
  up            = -100.d0
  sp            = zero
ELSE
  up            = tbck * DLOG( DABS( xpbck/( therm * gprot ) ) )
  sp            = ( xpbck    ) * ( 2.5d0 - up/tbck )
  up            = up - dmnp
END IF ! xpbck < 0

!-----------------------------------------------------------------------
!  Helium: ua and sa
!-----------------------------------------------------------------------

IF ( xabck <= epsilon ) THEN
  ua            = -100.d0
  sa            = zero
ELSE
  ua            = tbck * DLOG( DABS( xabck/( 32.d0 * therm * galph ) ) )
  sa            = ( xabck/4.d0 ) * ( 2.5d0 - ua/tbck )
END IF ! xabck < 0

!-----------------------------------------------------------------------
!  Heavies: utrhvy_i and sh
!-----------------------------------------------------------------------

strh            = zero
DO i = 1,n_hvy
  x_nuc_i       = xn(jshel,i_hvy(i))
  IF ( x_nuc_i <= epsilon ) THEN
    utrhvy_i    = -100.d0
    strh_i      = zero
  ELSE
    a_i         = a_nuc(i_hvy(i))
    utrhvy_i    = tbck * DLOG( DABS( x_nuc_i/( a_i * a_i * DSQRT(a_i) * therm * ghvy ) ) )
    strh_i      = ( x_nuc_i/a_i ) * ( 2.5d0 - utrhvy_i/tbck )
  END IF ! xhbck < 0
  strh          = strh + strh_i
END DO ! i

x_nuc_i         = xn(jshel,n_nucp1)
IF ( x_nuc_i <= epsilon ) THEN
  utrhvy_i      = -100.d0
  strh_i        = zero
ELSE
  a_i           = a_nuc_rep(jshel)
  utrhvy_i      = tbck * DLOG( DABS( x_nuc_i/( a_i * a_i * DSQRT(a_i) * therm * ghvy ) ) )
  strh_i        = ( x_nuc_i/a_i ) * ( 2.5d0 - utrhvy_i/tbck )
END IF ! xhbck < 0
strh            = strh + strh_i

!-----------------------------------------------------------------------
!  sd and uhat
!-----------------------------------------------------------------------

sd              = sn + sp + sa + strh + sh
uhat            = un - up

!***********************************************************************
! what about uhat and un etc?  this is a problem of first order?!! 

RETURN

CONTAINS

REAL (KIND=double) FUNCTION f3(z)
REAL (KIND=double) :: z
f3              = SIGN( ( DABS(z) )**third , z )
END FUNCTION f3

REAL (KIND=double) FUNCTION fphi(x1)
REAL (KIND=double) :: x1
fphi            = one - 3.0d0 * ( 0.5d0 - x1 )**2
END FUNCTION fphi

REAL (KIND=double) FUNCTION fbulk(x2)
REAL (KIND=double) :: x2
fbulk           = wnm + ws * ( one - 2.0d0 * x2 ) * ( one - 2.0d0 * x2 )
END FUNCTION fbulk

REAL (KIND=double) FUNCTION fsurf(x3)
REAL (KIND=double) :: x3
fsurf           = 290.d0 * x3 * x3 * ( one - x3 ) * ( one - x3 )
END FUNCTION fsurf

REAL (KIND=double) FUNCTION fcoul(x4)
REAL (KIND=double) :: x4
fcoul           = 0.75d0 * x4 * x4
END FUNCTION fcoul

REAL (KIND=double) FUNCTION faf(x5)
REAL (KIND=double) :: x5
faf             = DSIN( piby2 * DMIN1( one, x5 ) )/x5**(2.d0 * third)
END FUNCTION faf

END SUBROUTINE eosnuc_z
