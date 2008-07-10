SUBROUTINE sctednns( n, j, ij_ray, ik_ray, k, rho, t, ye, a0w, b0w, c0w, a1w, &
& b1w, c1w, rmdnns0, rmdnns1, rmdnns, arnns, brnns, crnns, aenns, benns, cenns, &
& arscte, brscte, crscte, arsctu, brsctu, crsctu, arscti, brscti, crscti, arsctd, &
& brsctd, crsctd, aesctu, besctu, cesctu, aescti, bescti, cescti, aesctd, besctd, &
& cesctd )
!-----------------------------------------------------------------------
!
!    File:         sctednns
!    Module:       sctednns
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/09/03
!
!    Purpose:
!      To compute quantities which are needed for editing neutrino-nucleon inelastic
!       scattering. In particular the quantities arns, brns, crns, aens, bens, and cens
!       are computed here. These are used for computing the total rate, rns, of
!       neutrino-nucleon inelastic scattering per incident neutrino of energy unu(k),
!       and the ratio, qns, of the mean energy transferred from the neutrino to the
!       nucleons to the incident neutrino energy. The quantities arscte, brscte, crscte,
!       arsctu, brsctu, crsctu, arscti, brscti, crscti, arsctd, brsctd, crsctd,
!       aesctu, besctu, cesctu, aescti, bescti, cescti, aesctd, besctd, cesctd are also computed
!       here. These are used to compute the rate of down-, isoenergetic-, and up-scattering,
!       and the mean change in neutrino energy for down-, and up-scattering. These
!       coefficients are used in subroutine editng. Finally, the quantities rmdns0, rmdns1,
!       rmdns, a0w, b0w, c0w, a1w, b1w,c1w are computed. The quantities rmdns0, rmdns1,
!       and rmdns are respectively the zero moment, first moment and total inverse mean
!       free path for neutrino-nucleon inelastic scattering. The quantities a0w, b0w, c0w, a1w,
!       b1w, and c1w are terms giving the in and out neutrino-nucleon inelastic scattering rates.
!
!      The quantity rnns is given by
!
!         rnns = tnes/psi0 = !*[ arnns*psi0 + brnns*psi1 + crnns ]/psi0
!
!       where
!
!         arnns  =  K   Int{ w2'dw'[ phi0out(w,w')( 1 - psi0(w') ) ] }
!         brnns  = -K/3 Int{ w2'dw'[ phi1out(w,w') ]psi1(w') }
!         crnns  =  0.0
!
!      (Note that tnes*( 4*pi/( hc)**3 )*w**2*dw is the net scattering rate for all
!       incident neutrinos in dw about w, so that
!
!   tnnns*( 4*pi/( hc)**3 )*w**2*dw/[( 4*pi/( hc)**3 )*w**2*dw*psi0] = tnnns/psi0
!       is the scattering rate per incident neutrino of energy w.
!
!      The quantity qnns is given by
!
!         qnns = - (efinal - einitial)/einitial
!
!       where efinal is the sum of the final energies of all the neutrinos scattered in a
!       unit time per unit incident neutrino energy, and einitial is the sum of the initial
!       energies of all the neutrinos scattered in a unit time per incident neutrino energy.
!       The quantities efinal and einitial are given by
!
!         efinal = K'*!*[ aenns*psi0 + benns*psi1 + cenns ]
!
!       where
!
!         aenns  =  K   Int{ w3'dw'[ phi0out(w,w')( 1 - psi0(w') ) ] }
!         benns  = -K/3 Int{ w3'dw'[ phi1out(w,w') ]psi1(w') }
!         cenns  =  0.0
!
!       and
!
!         einitial = K'*!*[ arnns*psi0 + brnns*psi1 + crnns ]*w
!
!    Variables that must be passed through common:
!
!  nnugp(n)    : number of energy zones for neutrinos of type n
!  unu(k)      : energy of energy zone k (MeV)
!  dunu(k)     : energy width of energy zone k (MeV)
!  psi0(j,k,n) : zeroth moment of of the neutrino occupation probability for neutrinos of
!                type n, energy zone k, radial zone j
!  psi1(j,k,n) : first moment of of the neutrino occupation probability for neutrinos of
!                type n, energy zone k, radial zone j
!
!    Input arguments:
!
!  n           : neutrino type
!  j           : radial zone number
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!  k           : energy zone number
!  rho         : matter density (g/cm**3)
!  t           : matter temperature (K)
!  cmpe        : electron chemical potential (MeV)
!
!    Input arguments (common):
!
!  isctnn      : 0 - neutrino-nucleon inelastic scattering turned off; isctn subroutines
!                 are bypassed; isctn scattering function arrays, if used, must be
!                 zeroed elsewhere
!              : 1 - neutrino-nucleon inelastic scattering turned on
!  iscat       : 0 - all scattering processes turned off.
!                1 - any scattering process is allowed to proceed
!                 if its key is turned on.
!  rhosctnnemn : density below which e-neutrino-nucleon inelastic scattering is turned off,
!                 function arrays are zeroed.
!  rhosctnnemx : density above which e-neutrino-nucleon inelastic scattering is turned off,
!                 function arrays are zeroed.
!  rhosctnntmn : density below which t-neutrino-nucleon inelastic scattering is turned off,
!                 function arrays are zeroed.
!  rhosctnntmx : density above which t-neutrino-nucleon inelastic scattering is turned off,
!                 function arrays are zeroed.
!
!    Output arguments:
!
!  a0w         : coefficient for computing in and out scattering rates
!  b0w         : coefficient for computing in and out scattering rates
!  c0w         : coefficient for computing in and out scattering rates
!  a1w         : coefficient for computing in and out scattering rates
!  b1w         : coefficient for computing in and out scattering rates
!  c1w         : coefficient for computing in and out scattering rates
!
!  arnns,brnns,crnns,aenns,benns,cenns
!              : coefficients from which the rate of neutrino-nucleon scattering per
!          n-neutrino of energy group k, and the fraction of n-neutrino energy
!          transferred to electrons per scattering event can be computed.
!          These coefficients are used in subroutine editn.
!
!  arscte,brscte,crscte,arsctu,brsctu,crsctu,arscti,brscti,crscti,
!  arsctd,brsctd,crsctd,aesctu,besctu,cesctu,aescti,bescti,cescti,
!  aesctd,besctd,cesctd
!              : coefficients from which the rate of down-, isoenergetic-, and up-scattering,
!          and the mean change in neutrino energy for down-, and up-scattering can be
!          computed for n-neutrino-electron scattering. These coefficients are used in
!          subroutine editng.
!
!  rmdnns0    : zero angular moment of the n-neutrino-nucleon inelastic
!          scattering inverse mean free path (returned from subroutine
!          sctednns for a given matter composition and initial and final
!          neutrino energy).
!
!  rmdnns1    : first angular moment of the n-neutrino-nucleon inelastic 
!          scattering inverse mean free path (returned from subroutine
!          sctednns for a given matter composition and initial and final
!          neutrino energy).
!
!  rmdnns     : n-neutrino-nucleon inelastic scattering inverse mean free
!          path (returned from subroutine sctednns for a given matter
!          composition and initial and final neutrino energy).
!
!    Subprograms called:
!      sctnnkrnl
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nez
USE numerical_module, ONLY : zero, one
USE physcnst_module, ONLY : kmev, pi

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY :  iscat, isctnn, rhosctnnemn, rhosctnnemx, rhosctnntmn, rhosctnntmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: j             ! radial zone index
INTEGER, INTENT(in)               :: k             ! incoming neutrino energy zone index
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)     :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)     :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)     :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)    :: rmdnns        ! neutrino inverse mean free path
REAL(KIND=double), INTENT(out)    :: rmdnns0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double), INTENT(out)    :: rmdnns1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double), INTENT(out)    :: a0w           ! coefficient for computing change in psi0
REAL(KIND=double), INTENT(out)    :: b0w           ! coefficient for computing change in psi0
REAL(KIND=double), INTENT(out)    :: c0w           ! coefficient for computing change in psi0
REAL(KIND=double), INTENT(out)    :: a1w           ! coefficient for computing change in psi1
REAL(KIND=double), INTENT(out)    :: b1w           ! coefficient for computing change in psi1
REAL(KIND=double), INTENT(out)    :: c1w           ! coefficient for computing change in psi1

REAL(KIND=double), INTENT(out)    :: arscte        ! coefficient for computing NNS scattering rate
REAL(KIND=double), INTENT(out)    :: brscte        ! coefficient for computing NNS scattering rate
REAL(KIND=double), INTENT(out)    :: crscte        ! coefficient for computing NNS scattering rate
REAL(KIND=double), INTENT(out)    :: arsctu        ! coefficient for computing NNS upscattering rate
REAL(KIND=double), INTENT(out)    :: brsctu        ! coefficient for computing NNS upscattering rate
REAL(KIND=double), INTENT(out)    :: crsctu        ! coefficient for computing NNS upscattering rate
REAL(KIND=double), INTENT(out)    :: arscti        ! coefficient for computing NNS isoscattering rate
REAL(KIND=double), INTENT(out)    :: brscti        ! coefficient for computing NNS isoscattering rate
REAL(KIND=double), INTENT(out)    :: crscti        ! coefficient for computing NNS isoscattering rate
REAL(KIND=double), INTENT(out)    :: arsctd        ! coefficient for computing NNS downscattering rate
REAL(KIND=double), INTENT(out)    :: brsctd        ! coefficient for computing NNS downscattering rate
REAL(KIND=double), INTENT(out)    :: crsctd        ! coefficient for computing NNS downscattering rate
REAL(KIND=double), INTENT(out)    :: aesctu        ! coefficient for computing NNS upscattering energy
REAL(KIND=double), INTENT(out)    :: besctu        ! coefficient for computing NNS upscattering energy
REAL(KIND=double), INTENT(out)    :: cesctu        ! coefficient for computing NNS upscattering energy
REAL(KIND=double), INTENT(out)    :: aescti        ! coefficient for computing NNS isoscattering energy
REAL(KIND=double), INTENT(out)    :: bescti        ! coefficient for computing NNS isoscattering energy
REAL(KIND=double), INTENT(out)    :: cescti        ! coefficient for computing NNS isoscattering energy
REAL(KIND=double), INTENT(out)    :: aesctd        ! coefficient for computing NNS downscattering energy
REAL(KIND=double), INTENT(out)    :: besctd        ! coefficient for computing NNS downscattering energy
REAL(KIND=double), INTENT(out)    :: cesctd        ! coefficient for computing NNS downscattering energy
REAL(KIND=double), INTENT(out)    :: arnns         ! coefficient for computing NNS rate
REAL(KIND=double), INTENT(out)    :: brnns         ! coefficient for computing NNS rate
REAL(KIND=double), INTENT(out)    :: crnns         ! coefficient for computing NNS rate
REAL(KIND=double), INTENT(out)    :: aenns         ! coefficient for computing NNS fractional eneergy transfer
REAL(KIND=double), INTENT(out)    :: benns         ! coefficient for computing NNS fractional eneergy transfer
REAL(KIND=double), INTENT(out)    :: cenns         ! coefficient for computing NNS fractional eneergy transfer

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

LOGICAL                           :: sctn_off

INTEGER                           :: kp            ! outcoming neutrino energy zone index
INTEGER                           :: istat         ! allocation status

REAL(KIND=double)                 :: enuin         ! incoming neutrino energy/kt
REAL(KIND=double)                 :: enuout        ! outgoing neutrino energy/kt
REAL(KIND=double)                 :: tmev          ! temperature (MeV)

REAL(KIND=double)                 :: w2dw          ! w2 * neutrino energy bin width
REAL(KIND=double)                 :: w3dw          ! w3 * neutrino energy bin width
REAL(KIND=double)                 :: x0            ! psi0
REAL(KIND=double)                 :: x1            ! psi1

REAL(KIND=double)                 :: c1            ! integrated scattering function
REAL(KIND=double)                 :: zt            ! integrated scattering function
REAL(KIND=double)                 :: asp           ! integrated scattering function
REAL(KIND=double)                 :: gp            ! integrated scattering function

REAL(KIND=double)                 :: f0in          ! zero moment of the NNNS in scattering function
REAL(KIND=double)                 :: f0ot          ! zero moment of the NNNS out scattering function
REAL(KIND=double)                 :: f1in          ! first moment of the NNNS in scattering function
REAL(KIND=double)                 :: f1ot          ! first moment of the NNNS out scattering function

REAL(KIND=double)                 :: art           ! coefficient for computing the mean NNNS rate
REAL(KIND=double)                 :: brt           ! coefficient for computing the mean NNNS rate
REAL(KIND=double)                 :: aqrt          ! coefficient for computing arnes
REAL(KIND=double)                 :: bqrt          ! coefficient for computing brnes
REAL(KIND=double)                 :: cqrt          ! coefficient for computing brnes
REAL(KIND=double)                 :: aqet          ! coefficient for computing aenes
REAL(KIND=double)                 :: bqet          ! coefficient for computing benes
REAL(KIND=double)                 :: cqet          ! coefficient for computing benes

REAL(KIND=double)                 :: artu          ! coefficient for computing the mean up-scattering NNNS rate
REAL(KIND=double)                 :: brtu          ! coefficient for computing the mean up-scattering NNNS rate
REAL(KIND=double)                 :: aetu          ! coefficient for computing the mean up-scattering NNNS energy change
REAL(KIND=double)                 :: betu          ! coefficient for computing the mean up-scattering NNNS energy change
REAL(KIND=double)                 :: cetu          ! coefficient for computing the mean up-scattering NNNS energy change

REAL(KIND=double)                 :: arti          ! coefficient for computing the mean iso-scattering NNNS rate
REAL(KIND=double)                 :: brti          ! coefficient for computing the mean iso-scattering NNNS rate
REAL(KIND=double)                 :: aeti          ! coefficient for computing the mean iso-scattering NNNS energy change
REAL(KIND=double)                 :: beti          ! coefficient for computing the mean iso-scattering NNNS energy change

REAL(KIND=double)                 :: artd          ! coefficient for computing the mean dow-nscattering NNNS rate
REAL(KIND=double)                 :: brtd          ! coefficient for computing the mean down-scattering NNNS rate
REAL(KIND=double)                 :: aetd          ! coefficient for computing the mean down-scattering NNNS energy change
REAL(KIND=double)                 :: betd          ! coefficient for computing the mean down-scattering NNNS energy change
REAL(KIND=double)                 :: cetd          ! coefficient for computing the mean down-scattering NNNS energy change

REAL(KIND=double)                 :: fexp          ! exponential function

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat0        ! zero moment of the NNNS function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat1        ! first moment of the NNNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat0d       ! d(scat0)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat1d       ! d(scat1)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat0t       ! d(scat0)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat1t       ! d(scat1)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat0y       ! d(scat0)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scat1y       ! d(scat1)/d(ye)

EXTERNAL fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for sctednns arrays')
 2001 FORMAT (' Deallocation problem for sctednns arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  No scattering if nnugp(n) = 0, iscat = 0  or isctnn = 0
!  No scattering if rho outside prescribed boundaries
!-----------------------------------------------------------------------

sctn_off           = .false.
IF ( nnugp(n) == 0  .or.  iscat == 0  .or.  isctnn == 0 ) sctn_off = .true.
IF ( n < 3   .and.  rho < rhosctnnemn ) sctn_off = .true.
IF ( n < 3   .and.  rho > rhosctnnemx ) sctn_off = .true.
IF ( n == 3  .and.  rho < rhosctnntmn ) sctn_off = .true.
IF ( n == 3  .and.  rho > rhosctnntmx ) sctn_off = .true.

IF ( sctn_off ) THEN

  rmdnns           = zero
  rmdnns0          = zero
  rmdnns1          = zero
  a0w              = zero
  b0w              = zero
  c0w              = zero
  a1w              = zero
  b1w              = zero
  c1w              = zero

  arscte           = zero
  brscte           = zero
  crscte           = zero
  arsctu           = zero
  brsctu           = zero
  crsctu           = zero
  arscti           = zero
  brscti           = zero
  crscti           = zero
  arsctd           = zero
  brsctd           = zero
  crsctd           = zero
  aesctu           = zero
  besctu           = zero
  cesctu           = zero
  aescti           = zero
  bescti           = zero
  cescti           = zero
  aesctd           = zero
  besctd           = zero
  cesctd           = zero

  arnns            = zero
  brnns            = zero
  crnns            = zero
  aenns            = zero
  benns            = zero
  cenns            = zero

  RETURN

END IF

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (scat0(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scat1(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scat0d(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0d    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scat1d(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1d    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scat0t(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0t    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scat1t(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1t    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scat0y(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0y    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scat1y(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1y    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize scattering functions
!-----------------------------------------------------------------------

scat0              = zero
scat1              = zero
scat0d             = zero
scat1d             = zero
scat0t             = zero
scat1t             = zero
scat0y             = zero
scat1y             = zero

!-----------------------------------------------------------------------
!  Compute variables related to n-type neutrinos.
!
!  tmev       : t (MeV)
!  enuin      : (in beam neutrino energy)/kT
!  enuout     : (out beam neutrino energy)/kT
!-----------------------------------------------------------------------

tmev               = t * kmev
enuin              = unu(j,k)/tmev

!-----------------------------------------------------------------------
!  Initialize coefficients for computing scattering rates
!-----------------------------------------------------------------------

c1                 = zero
zt                 = zero
asp                = zero
gp                 = zero

art                = zero
brt                = zero
artu               = zero
brtu               = zero
arti               = zero
brti               = zero
artd               = zero
brtd               = zero
aetu               = zero
betu               = zero
cetu               = zero
aeti               = zero
beti               = zero
aetd               = zero
betd               = zero
cetd               = zero
aqrt               = zero
bqrt               = zero
cqrt               = zero
aqet               = zero
bqet               = zero
cqet               = zero

!-----------------------------------------------------------------------
!  Get scattering kernels
!-----------------------------------------------------------------------

CALL sctnnkrnl(n,j,ij_ray, ik_ray, k, rho, t, ye, scat0, scat1, scat0d, &
& scat1d, scat0t, scat1t, scat0y, scat1y, nez )

!-----------------------------------------------------------------------
!                   ||||| Begin loop over kp |||||
!-----------------------------------------------------------------------

DO kp = 1,nnugp(n)

  enuout           = unu(j,kp)/tmev

  w2dw             = 2.d0 * pi * unu(j,kp) * unu(j,kp) * dunu(j,kp)
  w3dw             = unu(j,kp) * w2dw
  x0               = psi0(j,kp,n)
  x1               = psi1(j,kp,n)

  IF ( kp < k ) THEN

    f0ot           = scat0(kp)
    f0in           = fexp( enuout - enuin ) * scat0(kp)
    f1in           = zero
    f1ot           = zero

  ELSE IF ( kp == k ) THEN

    f0ot           = scat0(kp)
    f0in           = scat0(kp)
    f1in           = zero
    f1ot           = zero

  ELSE IF ( kp > k ) THEN

    f0in           = scat0(kp)
    f0ot           = fexp( enuin - enuout ) * scat0(kp)
    f1in           = zero
    f1ot           = zero

  END IF

  zt               = zt   + ( f0in * x0 + f0ot * ( one - x0 ) ) * w2dw
  asp              = asp  + x1 * w2dw * ( f1in - f1ot )
  c1               = c1   + f0in * x0 * w2dw
  gp               = gp   + x1 * w2dw * f1in

  art              = art  + f0in * x0 * w2dw
  brt              = brt  + f1in * x1 * w2dw
  aqrt             = aqrt + f0ot * ( one - x0 ) * w2dw
  bqrt             = bqrt - f1ot * x1 * w2dw
  aqet             = aqet + f0ot * ( one - x0 ) * w3dw
  bqet             = bqet - f1ot * x1 * w3dw

  IF ( kp < k ) THEN

    artu           = artu + f0in * x0 * w2dw
    brtu           = brtu + f1in * x1 * w2dw
    aetu           = aetu - f0in * x0 * w2dw
    betu           = betu - f1in * x1 * w2dw
    cetu           = cetu + f0in * x0 * w2dw
    aetd           = aetd - f0ot * ( one - x0 ) * w2dw
    betd           = betd + f1ot * x1 * w2dw

  ELSE IF ( kp == k ) THEN

    arti           = arti + f0in * x0 * w2dw
    brti           = brti + f1in * x1 * w2dw
    aeti           = aeti + ( f0in - f0ot ) * ( one - x0 ) * w2dw
    beti           = beti + ( f1ot - f1in ) * x1 * w2dw

  ELSE IF ( kp > k ) THEN

    artd           = artd + f0in * x0 * w2dw
    brtd           = brtd + f1in * x1 * w2dw
    aetd           = aetd - f0in * x0 * w2dw 
    betd           = betd - f1in * x1 * w2dw
    cetd           = cetd + f0in * x0 * w2dw
    aetu           = aetu - f0ot * ( one - x0 ) * w2dw
    betu           = betu + f1ot * x1 * w2dw

  END IF

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

END DO ! kp

rmdnns             =  zt
rmdnns0            =  zt
a0w                = -zt
b0w                = -asp/3.d0
c0w                =  c1
a1w                = -asp
b1w                = -zt
c1w                =  gp

arscte             = -art
brscte             = -brt/3.d0
crscte             =  art
arsctu             = -artu
brsctu             = -brtu/3.d0
crsctu             =  artu
arscti             = -arti
brscti             = -brti/3.d0
crscti             =  arti
arsctd             = -artd
brsctd             = -brtd/3.d0
crsctd             =  artd
aesctu             =  aetu
besctu             =  betu/3.d0
cesctu             =  cetu
aescti             =  aeti
bescti             =  beti/3.d0
cescti             =  zero
aesctd             =  aetd
besctd             =  betd/3.d0
cesctd             =  cetd

arnns              =  aqrt
brnns              =  bqrt
crnns              =  zero                                                       
aenns              =  aqet
benns              =  bqet
cenns              =  zero                                                       

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (scat0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scat1, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scat0d, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0d    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scat1d, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1d    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scat0t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0t    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scat1t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1t    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scat0y, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat0y    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scat1y, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scat1y    '; WRITE (nlog,2001) var_name; END IF

RETURN 
END SUBROUTINE sctednns
