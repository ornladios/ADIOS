SUBROUTINE scatical( rho, t, enu, xn, xp, xhe, xh, ah, zh, rmdnns, rmdnps, &
& rmdnbns, rmdnbps, rmdnhes, rmdnhs, coh, cohb )
!-----------------------------------------------------------------------
!
!    File:         scatical
!    Module:       scatical
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To calculate the zero and first legendre coefs for the n-type isoenergetic
!       scattering functions.
!
!      These are included in the multi-group diffusion equations, which
!       have the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2 = xw + yw*psi0 + zw*dpsi0/dr
!
!       where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!      Neutrino-nucleus (nucleon) isoenergetic scattering contributes to the
!       terms a0w,a1w, b0w, b1w, c0w, c1w as follows:
!
!            b1w  =  K w2 [ sct1(w) - sct0(w) ]
!
!       where
!
!            K    = 2*pi/! * 1/(hc)**3
!
!       and where sct0(w) and sct1(w) are the zero and first Legendre moments
!        for neutrino-nucleus (nucleon) iso-energetic scattering.
!
!    Subprograms called:
!      etaxx, scatiicr
!
!    Input arguments:
!
!   rho       : matter density (g/cm**3)
!   t         : matter temperature (K)
!   enu       : neutrino energy (MeV)
!   xn        : mass fraction of free neutrons
!   xp        : mass fraction of free protons
!   xhe       : mass fraction of helium nuclei
!   xh        : mass fraction of heavy nuclei
!   ah        : mean heavy nucleus mass number
!   zh        : mean heavy nucleus charge number
!   in        : 0, neutrino-neutron isoenergetic scattering omitted
!               1, neutrino-neutron isoenergetic scattering included
!   ip        : 0, neutrino-proton isoenergetic scattering omitted
!               1, neutrino-proton isoenergetic scattering included
!   ihe       : 0, neutrino-helium isoenergetic scattering omitted
!               1, neutrino-helium isoenergetic scattering included
!   iheavy    : 0, neutrino-heavy nucleus isoenergetic scattering omitted
!               1, neutrino-heavy nucleus isoenergetic scattering included
!   iscat     : 0, all neutrino scattering processes omitted
!               1, neutrino scattering processes not necessarily omitted
!
!    Output arguments:
!
!   rmdnns    : neutrino-neutron scattering function
!   rmdnps    : neutrino-proton scattering function
!   rmdnbns   : neutrino-neutron scattering function
!   rmdnbps   : neutrino-proton scattering function
!   rmdnhes   : neutrino-helium scattering function
!   rmdnhs    : neutrino-nucleus scattering function
!   coh       : b1w  =  K w2 [ sct1(w) - sct0(w) ] = rmdnns + rmdnps + rmdnhes + rmdnhs
!   cohb      : b1w  =  K w2 [ sct1(w) - sct0(w) ] = rmdnbns + rmdnbps + rmdnhes + rmdnhs
!
!    Input arguments (common):
!
!   cv        : weak interaction constant
!   ca        : weak interaction constant
!   ga        : weak interaction constant
!
!    Output arguments (common):
!
!   rmdnpsi   : ith moment of the neutrino-proton scattering function       
!   rmdnnsi   : ith moment of the neutrino-neutron scattering function      
!   rmdnnsi   : ith moment of the neutrino-helium scattering function      
!   rmdnnsi   : ith moment of the neutrino-nucleus scattering function       
!
!    Include files:
!      kind_module, numerical_module, physcnst_module
!      prb_cntl_module, scat_i_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : Gw, mp, hbar, cvel, pi, cv, ga, rmu

USE prb_cntl_module, ONLY : iscat, in, ip, ihe, iheavy, isctn, rhosctnemn, rhosctntmn
USE scat_i_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)     :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)     :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)     :: enu           ! neutrino energy (MeV)
REAL(KIND=double), INTENT(in)     :: xn            ! neutron mass fraction
REAL(KIND=double), INTENT(in)     :: xp            ! proton mass fraction
REAL(KIND=double), INTENT(in)     :: xhe           ! helium mass fraction
REAL(KIND=double), INTENT(in)     :: xh            ! heavy nucleus mass fraction
REAL(KIND=double), INTENT(in)     :: ah            ! heavy nucleus mass number
REAL(KIND=double), INTENT(in)     :: zh            ! heavy nucleus charge number

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)    :: rmdnns        ! mfp^-1 for neutrino-neutron scattering
REAL(KIND=double), INTENT(out)    :: rmdnps        ! mfp^-1 for neutrino-proton scattering
REAL(KIND=double), INTENT(out)    :: rmdnbns       ! mfp^-1 for antineutrino-neutron scattering
REAL(KIND=double), INTENT(out)    :: rmdnbps       ! mfp^-1 for antineutrino-proton scattering
REAL(KIND=double), INTENT(out)    :: rmdnhes       ! mfp^-1 for neutrino-helium scattering
REAL(KIND=double), INTENT(out)    :: rmdnhs        ! mfp^-1 for neutrino-heavy nucleus scattering
REAL(KIND=double), INTENT(out)    :: coh           ! sum of above for neutrinos
REAL(KIND=double), INTENT(out)    :: cohb          ! sum of above for antineutrinos

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                           :: first = .true.

REAL(KIND=double)                 :: tthird        ! 2/3
REAL(KIND=double)                 :: g2            ! constant
REAL(KIND=double)                 :: cc            ! constant

REAL(KIND=double)                 :: hvp           ! proton vector coupling constant
REAL(KIND=double)                 :: hap           ! proton axial vector coupling constant
REAL(KIND=double)                 :: hvn           ! neutron vector coupling constant
REAL(KIND=double)                 :: han           ! neutron axial vector coupling constant
REAL(KIND=double)                 :: cv0           ! coupling constant
REAL(KIND=double)                 :: cv1           ! coupling constant

REAL(KIND=double)                 :: ap0           ! proton zero moment coupling constant
REAL(KIND=double)                 :: ap1           ! proton first moment coupling constant
REAL(KIND=double)                 :: an0           ! neutron zero moment coupling constant
REAL(KIND=double)                 :: an1           ! neutron first moment coupling constant

REAL(KIND=double)                 :: eeche         ! helium opacity calculation
REAL(KIND=double)                 :: eche          ! helium opacity calculation
REAL(KIND=double)                 :: eche2         ! helium opacity calculation
REAL(KIND=double)                 :: eche3         ! helium opacity calculation
REAL(KIND=double)                 :: b0he          ! helium opacity calculation
REAL(KIND=double)                 :: saghe         ! helium opacity calculation
REAL(KIND=double)                 :: sbghe         ! helium opacity calculation

REAL(KIND=double)                 :: eec           ! heavy nucleus opacity calculation
REAL(KIND=double)                 :: ec            ! heavy nucleus opacity calculation
REAL(KIND=double)                 :: ec2           ! heavy nucleus opacity calculation
REAL(KIND=double)                 :: ec3           ! heavy nucleus opacity calculation
REAL(KIND=double)                 :: b0h           ! heavy nucleus opacity calculation
REAL(KIND=double)                 :: sag           ! heavy nucleus opacity calculation
REAL(KIND=double)                 :: sbg           ! heavy nucleus opacity calculation
REAL(KIND=double)                 :: xnh           ! heavy nucleus opacity calculation

REAL(KIND=double)                 :: aisov         ! vector coupling constant
REAL(KIND=double)                 :: aisosc        ! vector coupling constant
REAL(KIND=double)                 :: a01           ! vector coupling constant
REAL(KIND=double)                 :: a02           ! coupling constant
REAL(KIND=double)                 :: e2            ! neutrino energy squared
REAL(KIND=double)                 :: etann         ! neutron number corrected for blocking
REAL(KIND=double)                 :: xnn           ! neutron number corrected for blocking
REAL(KIND=double)                 :: etapp         ! proton number corrected for blocking
REAL(KIND=double)                 :: xnp           ! proton number corrected for blocking
REAL(KIND=double)                 :: xnhe          ! helium number
REAL(KIND=double)                 :: xheaa         ! heavy nucleus number * a_he^2
REAL(KIND=double)                 :: xhaa          ! heavy nucleus number * ah^2

REAL(KIND=double)                 :: ciicr         ! ion-ion correlation correction

REAL(KIND=double)                 :: xi_p_wm       ! weak magnetism correction for neutrino-proton scattering
REAL(KIND=double)                 :: xi_n_wm       ! weak magnetism correction for neutrino-neutron scattering
REAL(KIND=double)                 :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering
REAL(KIND=double)                 :: xib_n_wm      ! weak magnetism correction for antineutrino-neutron scattering

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN

  first            = .false.
  tthird           = 2.d0/3.d0
  g2               = ( Gw/mp**2 )**2 * hbar**5 * cvel**6
  cc               = ( 2.d0 * pi )/( cvel * ( 2.d0 * pi * hbar * cvel )**3 ) * ( 4.d0 * pi ) * g2

  hvp              =  one - cv
  hap              =  0.5d0 * ga
  hvn              = -0.5d0
  han              = -0.5d0 * ga
  ap0              =  hvp**2 + 3.d0 * hap**2
  ap1              =  hvp**2 - hap**2
  an0              =  hvn**2 + 3.d0 * han**2
  an1              =  hvn**2 - han**2
  cv0              =  half * ( hvp + hvn )
  cv1              =  hvp - hvn
  aisosc           =  cv0
  a01              =  cv0 * cv0

END IF

!-----------------------------------------------------------------------
!  Initialize  mfp^-1's
!-----------------------------------------------------------------------

rmdnps             = zero
rmdnns             = zero
rmdnbps            = zero
rmdnbns            = zero
rmdnhes            = zero
rmdnhs             = zero
coh                = zero
cohb               = zero
rmdnps0            = zero
rmdnns0            = zero
rmdnbps0           = zero
rmdnbns0           = zero
rmdnhes0           = zero
rmdnhs0            = zero
rmdnps1            = zero
rmdnns1            = zero
rmdnbps1           = zero
rmdnbns1           = zero
rmdnhes1           = zero
rmdnhs1            = zero

!-----------------------------------------------------------------------
!  No neutrino nuclei (nucleon) isoenergetic scattering if iscat = 0.
!-----------------------------------------------------------------------

IF ( iscat == 0 ) RETURN

!-----------------------------------------------------------------------
!  Zero neutrino-nucleus scattering rate if ah and zh below 1.d-10
!-----------------------------------------------------------------------

IF ( ah < 1.d-10  .or.  zh < 1.d-10 ) THEN
  a02              = zero
ELSE
  aisov            = half * cv1 * ( 2.d0 * zh - ah )/ah
  a02              = ( aisosc + aisov ) * ( aisosc + aisov )
END IF

!-----------------------------------------------------------------------
!  Quatities needed to compute the scattering functions.
!-----------------------------------------------------------------------

e2                 = enu * enu
CALL etaxx( rho, t, xn, xp, etann, etapp )
xnn                = etann
xnp                = etapp
xnhe               = ( xhe/4.d0             ) * rho/rmu
xnh                = ( xh /( ah + epsilon ) ) * rho/rmu
xheaa              = xnhe * 16.d0
xhaa               = xnh * ah * ah
b0he               = 1.21d-5
b0h                = 4.80d-6 * ( ah )**tthird

!-----------------------------------------------------------------------
!  Quantities needed to compute the coherent scattering on helium.
!-----------------------------------------------------------------------
 
eche               = 4.d0 * b0he * e2
eche2              = eche * eche
eche3              = eche * eche2

IF ( eche < 0.1d0 ) THEN
  saghe            = 1.d0 - tthird * eche
  sbghe            = ( 5.d0 - eche2 )/15.d0
ELSE
  eeche            = zero
IF ( eche < 15.d0 ) eeche = DEXP( -2.d0 * eche )
  saghe            = ( eche - 0.5d0 + 0.5d0 * eeche )/eche2
  sbghe            = ( eche2 - 1.5d0 * eche + 1.d0 - ( 0.5d0 * eche + 1.d0 ) * eeche )/ eche3
END IF ! eche < 0.1

!-----------------------------------------------------------------------
!  Quantities needed to compute the coherent scattering on heavy nuclei.
!-----------------------------------------------------------------------

ec                 = 4.d0 * b0h * e2
ec2                = ec * ec
ec3                = ec * ec2

IF ( ec < 0.1d0 ) THEN
  sag              = 1.d0 - tthird * ec
  sbg              = ( 5.d0 - ec2 )/15.d0
ELSE
  eec              = 0.d0
  IF ( ec < 15.d0 ) eec = DEXP( -2.d0 * ec )
  sag              = ( ec - 0.5d0 + 0.5d0 * eec )/ec2
  sbg              = ( ec2 - 1.5d0 * ec + 1.d0 - ( 0.5d0 * ec + 1.d0 ) * eec )/ec3
END IF

!-----------------------------------------------------------------------
!  Inverse mean free paths for coherent scattering.
!-----------------------------------------------------------------------

rmdnps0            = cc * e2 * xnp * ap0
rmdnns0            = cc * e2 * xnn * an0
rmdnbps0           = rmdnps0
rmdnbns0           = rmdnns0
rmdnhes0           = cc * e2 * xheaa * a01 * saghe
rmdnhs0            = cc * e2 * xhaa * a02 * sag
rmdnps1            = cc * e2 * xnp * ( -ap1/3.d0 )
rmdnns1            = cc * e2 * xnn * ( -an1/3.d0 )
rmdnbps1           = rmdnps1
rmdnbns1           = rmdnns1
rmdnhes1           = cc * e2 * xheaa * a01 * ( -sbghe )
rmdnhs1            = cc * e2 * xhaa * a02 * ( -sbg )
rmdnps             = rmdnps0  + rmdnps1
rmdnns             = rmdnns0  + rmdnns1
rmdnbps            = rmdnps
rmdnbns            = rmdnns
rmdnhes            = rmdnhes0 + rmdnhes1
rmdnhs             = rmdnhs0  + rmdnhs1

!-----------------------------------------------------------------------
!  Ion-ion correlation correction for coherent scattering.
!-----------------------------------------------------------------------

CALL scatiicr( rho, t, enu, xh, ah, zh, ciicr )

rmdnhs             = rmdnhs  * ciicr
rmdnhs0            = rmdnhs0 * ciicr
rmdnhs1            = rmdnhs1 * ciicr

!-----------------------------------------------------------------------
!  Weak magnetism corrections for neutrino and antineutrino neutron and
!   proton scattering.
!-----------------------------------------------------------------------

CALL nc_weak_mag( enu, xi_p_wm, xi_n_wm, xib_p_wm, xib_n_wm )

rmdnps0            = rmdnps0  * xi_p_wm
rmdnns0            = rmdnns0  * xi_n_wm
rmdnbps0           = rmdnbps0 * xib_p_wm
rmdnbns0           = rmdnbns0 * xib_n_wm
rmdnps1            = rmdnps1  * xi_p_wm
rmdnns1            = rmdnns1  * xi_n_wm
rmdnbps1           = rmdnbps1 * xib_p_wm
rmdnbns1           = rmdnbns1 * xib_n_wm
rmdnps             = rmdnps   * xi_p_wm
rmdnns             = rmdnns   * xi_n_wm
rmdnbps            = rmdnbps  * xib_p_wm
rmdnbns            = rmdnbns  * xib_n_wm

!-----------------------------------------------------------------------
!  Incorporate scattering keys.
!-----------------------------------------------------------------------

IF ( in      == 0 ) rmdnns  = zero
IF ( ip      == 0 ) rmdnps  = zero
IF ( ihe     == 0 ) rmdnhes = zero
IF ( iheavy  == 0 ) rmdnhs  = zero
IF ( isctn   /= 0  .and.  rho > DMAX1( rhosctnemn, rhosctntmn ) ) THEN
  rmdnps0          = zero
  rmdnns0          = zero
  rmdnbps0         = zero
  rmdnbns0         = zero
  rmdnps1          = zero
  rmdnns1          = zero
  rmdnbps1         = zero
  rmdnbns1         = zero
  rmdnps           = zero
  rmdnns           = zero
  rmdnbps          = zero
  rmdnbns          = zero
END IF

!-----------------------------------------------------------------------
!  Coherent (net) scattering inverse mean free path
!-----------------------------------------------------------------------

coh                = rmdnps  + rmdnns  + rmdnhes + rmdnhs
cohb               = rmdnbps + rmdnbns + rmdnhes + rmdnhs

RETURN
END SUBROUTINE scatical
