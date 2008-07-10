SUBROUTINE eosnuc_e( dbck, tbck, xn, nnc, a_nuc_rep, z_nuc_rep, be_nuc_rep )
!-----------------------------------------------------------------------
!
!    File:         eosnuc_e
!    Module:       eosnuc_e
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/15/04
!
!    Purpose:
!      To compute the nuclear component of the equation of state
!       for material not in nuclear statistical equilibrium. The
!       composition of the material is given by the array xn(i).
!
!      The nucleons and nuclei are assumed to comprise an ideal
!       Boltzmann gas.
!
!    Variables that must be passed through common:
!      xn(i): mass fraction of the ith nucleus in radial zone j.
!      Cooperstein-BCK eos variables.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nnc         : composition array dimension
!  dbck        : density (baryons fm^{-3})
!  tbck        : temperature (MeV)
!  xn          : nuclaer mass fractions
!  a_nuc_rep   : mass number of nucleus representative nucleus
!  z_nuc_rep   : charge number of nucleus representative nucleus
!  be_nuc_rep  : binding energy of nucleus representative nucleus
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  eos_bck_module, eos_bck_module, eos_snc_x_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, third, one, epsilon
USE physcnst_module, ONLY : mb, pi, hbarc, wnm, ws

USE edit_module, ONLY : nlog
USE eos_bck_module, ONLY : zabck, b, eh, sh, ed
USE eos_snc_x_module, ONLY : nuc_number, a_nuc, z_nuc, be_nuc, a_name

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nnc             ! composition array dimension

REAL(KIND=double), INTENT(in)     :: dbck            ! density (baryons fm^{-3})
REAL(KIND=double), INTENT(in)     :: tbck            ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION (nnc) :: xn ! nuclaer mass fractions
REAL(KIND=double), INTENT(in)     :: a_nuc_rep       ! mass number of nucleus representative nucleus
REAL(KIND=double), INTENT(in)     :: z_nuc_rep       ! charge number of nucleus representative nucleus
REAL(KIND=double), INTENT(in)     :: be_nuc_rep      ! binding energy of nucleus representative nucleus

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

REAL(KIND=double)                 :: av_inv         ! neutrino absorption and emission time step
REAL(KIND=double)                 :: zv             ! neutrino absorption and emission time step

REAL(KIND=double), PARAMETER      :: gprot = 2.d0   ! neutrino absorption and emission time step
REAL(KIND=double), PARAMETER      :: gnewt = 2.d0   ! neutrino absorption and emission time step
REAL(KIND=double), PARAMETER      :: galph = 1.d0   ! neutrino absorption and emission time step
REAL(KIND=double), PARAMETER      :: ghvy  = 1.d0   ! neutrino absorption and emission time step
REAL(KIND=double), PARAMETER      :: x_min = 1.d-30 ! neutrino absorption and emission time step

REAL(KIND=double)                 :: c0             ! constant for computing the chemical potential
REAL(KIND=double)                 :: a0             ! constant for computing the nuclear level density
REAL(KIND=double)                 :: xms            ! effective nuclear mass at high density
REAL(KIND=double)                 :: xm0            ! effective nuclear mass at low density
REAL(KIND=double)                 :: xm0ms          ! xm0 - xms
REAL(KIND=double)                 :: tc             ! nuclear level density parameter
REAL(KIND=double)                 :: a13            ! cube root of the nuclear mass number
REAL(KIND=double)                 :: piby2          ! constant

REAL(KIND=double)                 :: av             ! harmonic mean nuclear mass number
REAL(KIND=double)                 :: rn             ! number of nuclei per baryon

!-----------------------------------------------------------------------
!        Formats.
!-----------------------------------------------------------------------

  101 FORMAT (' nnc =',i4,' > ncd=',i4,' in subroutine eosnuc_e.')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( nnc > ncd ) THEN
  WRITE (nlog,101) nnc,ncd
  STOP
END IF

IF ( first ) THEN
  first            = .false.

!........Constants

  c0               = ( mb/( 2.d0 * pi * hbarc**2 ) )**( 3.d0/2.d0 )
  a0               = .067d0
  piby2            = 1.570796327d0
  xms              = 0.7d0
  xm0              = 2.0d0
  xm0ms            = xm0 - xms
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
  END DO ! i

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

be_nuc(n_nucp1)    = be_nuc_rep
a_nuc(n_nucp1)     = a_nuc_rep
z_nuc(n_nucp1)     = z_nuc_rep
zabck              = DMIN1( z_nuc(n_nucp1)/( a_nuc(n_nucp1) + epsilon ), 1.d0 )
phib(n_nucp1)      = fphi(zabck)
thetab(n_nucp1)    = 1.0d0
afb(n_nucp1)       = a0 * faf(phib(n_nucp1) * thetab(n_nucp1))

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
  av_inv           = av_inv +             xn(i)/( a_nuc(i) + epsilon )
  zv               = zv     + z_nuc(i)  * xn(i)/( a_nuc(i) + epsilon )
  b                = b      - be_nuc(i) * xn(i)/( a_nuc(i) + epsilon )
END DO
av                 = one/( av_inv + epsilon )
zv                 = zv * av
rn                 = one/( av + epsilon )

!-----------------------------------------------------------------------
!        Compute
!  eh : energy of nuclei w/o translation
!  ed : energy of drip (including heavy translation)
!
!  e0 = ( egy0 - ( ye - yeFe )*dmnp ) is the energy per baryon of an
!   isolated iron nucleus plus the rest mass difference per baryon
!   between matter at ye and that at the ye of iron. The value of eh is
!   relative to that energy.!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  eh
!-----------------------------------------------------------------------

eh              = zero
sh              = zero

DO i = 1,n_hvy
  xmstarb(i)    = xms + xm0ms /( one + tbck/tc )**2
  sh            = sh + xn(i_hvy(i)) * afb(i) * tbck * ( 2.d0 * xmstarb(i) - 2.d0 * tbck/tc &
&               / ( one + tbck/tc ) * ( xmstarb(i) - xms ) )
  eh            = eh + tbck * sh
END DO

xmstarb(n_nucp1) = xms + xm0ms /( one + tbck/tc )**2
sh              = sh + xn(n_nucp1) * afb(n_nucp1) * tbck * ( 2.d0 * xmstarb(n_nucp1) - 2.d0 * tbck/tc &
&               / ( one + tbck/tc ) * ( xmstarb(n_nucp1) - xms ) )
eh              = eh + tbck * sh

!-----------------------------------------------------------------------
!  ed
!-----------------------------------------------------------------------

ed              = 1.5d0 * rn * tbck

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

END SUBROUTINE eosnuc_e
