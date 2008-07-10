SUBROUTINE sctednAs( n, j, ij_ray, ik_ray, k, rho, t, ye, a0w, b0w, c0w, a1w, &
& b1w, c1w, rmdnAs0, rmdnAs1, rmdnAs, arnAs, brnAs, crnAs, aenAs, benAs, cenAs, &
& arscte, brscte, crscte, arsctu, brsctu, crsctu, arscti, brscti, crscti, arsctd, &
& brsctd, crsctd, aesctu, besctu, cesctu, aescti, bescti, cescti, aesctd, besctd, &
& cesctd )
!-----------------------------------------------------------------------
!
!    File:         sctednAs
!    Module:       sctednAs
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/09/03
!
!    Purpose:
!      To compute quantities which are needed for editing neutrino-nucleus inelastic
!       scattering. In particular the quantities aras, bras, cras, aeas, beas, and ceas
!       are computed here. These are used for computing the total rate, rns, of
!       neutrino-nucleus inelastic scattering per incident neutrino of energy unu(k),
!       and the ratio, qns, of the mean energy transferred from the neutrino to the
!       nucleis to the incident neutrino energy. The quantities arscte, brscte, crscte,
!       arsctu, brsctu, crsctu, arscti, brscti, crscti, arsctd, brsctd, crsctd,
!       aesctu, besctu, cesctu, aescti, bescti, cescti, aesctd, besctd, cesctd are also computed
!       here. These are used to compute the rate of down-, isoenergetic-, and up-scattering,
!       and the mean change in neutrino energy for down-, and up-scattering. These
!       coefficients are used in subroutine editng. Finally, the quantities rmdas0, rmdas1,
!       rmdas, a0w, b0w, c0w, a1w, b1w,c1w are computed. The quantities rmdas0, rmdas1,
!       and rmdas are respectively the zero moment, first moment and total inverse mean
!       free path for neutrino-nucleus inelastic scattering. The quantities a0w, b0w, c0w, a1w,
!       b1w, and c1w are terms giving the in and out neutrino-nucleus inelastic scattering rates.
!
!      The quantity rns is given by
!
!         rns = tnes/psi0 = !*[ aras*psi0 + bras*psi1 + cras ]/psi0
!
!       where
!
!         arnAs  =  K   Int{ w2'dw'[ phi0out(w,w')( 1 - psi0(w') ) ] }
!         brnAs  = -K/3 Int{ w2'dw'[ phi1out(w,w') ]psi1(w') }
!         crnAs  =  0.0
!
!      (Note that tnes*( 4*pi/( hc)**3 )*w**2*dw is the net scattering rate for all
!       incident neutrinos in dw about w, so that
!
!   tnAns*( 4*pi/( hc)**3 )*w**2*dw/[( 4*pi/( hc)**3 )*w**2*dw*psi0] = tnAns/psi0
!       is the scattering rate per incident neutrino of energy w.
!
!      The quantity qnAs is given by
!
!         qnAs = - (efinal - einitial)/einitial
!
!       where efinal is the sum of the final energies of all the neutrinos scattered in a
!       unit time per unit incident neutrino energy, and einitial is the sum of the initial
!       energies of all the neutrinos scattered in a unit time per incident neutrino energy.
!       The quantities efinal and einitial are given by
!
!         efinal = K'*!*[ aenAs*psi0 + benAs*psi1 + cenAs ]
!
!       where
!
!         aenAs  =  K   Int{ w3'dw'[ phi0out(w,w')( 1 - psi0(w') ) ] }
!         benAs  = -K/3 Int{ w3'dw'[ phi1out(w,w') ]psi1(w') }
!         cenAs  =  0.0
!
!       and
!
!         einitial = K'*!*[ arnAs*psi0 + brnAs*psi1 + crnAs ]*w
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
!  isctnA      : 0 - neutrino-nuclei inelastic scattering turned off; isctn subroutines
!                 are bypassed; isctn scattering function arrays, if used, must be
!                 zeroed elsewhere
!              : 1 - neutrino-nucleus inelastic scattering turned on
!  iscat       : 0 - all scattering processes turned off.
!                1 - any scattering process is allowed to proceed
!                 if its key is turned on.
!  rhosctnAemn : density below which e-neutrino-nucleus inelastic scattering is turned off,
!                 function arrays are zeroed.
!  rhosctnAemx : density above which e-neutrino-nucleus inelastic scattering is turned off,
!                 function arrays are zeroed.
!  rhosctnAtmn : density below which t-neutrino-nucleus inelastic scattering is turned off,
!                 function arrays are zeroed.
!  rhosctnAtmx : density above which t-neutrino-nucleus inelastic scattering is turned off,
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
!  arnAs,brnAs,crnAs,aenAs,benAs,cenAs
!              : coefficients from which the rate of neutrino-nucleus scattering per
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
!  rmdnAs0    : zero angular moment of the n-neutrino-nucleus inelastic
!          scattering inverse mean free path (returned from subroutine
!          sctednAs for a given matter composition and initial and final
!          neutrino energy).
!
!  rmdnAs1    : first angular moment of the n-neutrino-nucleus inelastic 
!          scattering inverse mean free path (returned from subroutine
!          sctednAs for a given matter composition and initial and final
!          neutrino energy).
!
!  rmdnAs     : n-neutrino-nucleus inelastic scattering inverse mean free
!          path (returned from subroutine sctednAs for a given matter
!          composition and initial and final neutrino energy).
!
!    Subprograms called:
!      sctnAkrnl
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
USE prb_cntl_module, ONLY :  iscat, isctnA, rhosctnAemn, rhosctnAemx, rhosctnAtmn, rhosctnAtmx

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

REAL(KIND=double), INTENT(out)    :: rmdnAs        ! neutrino inverse mean free path
REAL(KIND=double), INTENT(out)    :: rmdnAs0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double), INTENT(out)    :: rmdnAs1       ! first moment of the neutrino inverse mean free path

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
REAL(KIND=double), INTENT(out)    :: arnAs         ! coefficient for computing NNS rate
REAL(KIND=double), INTENT(out)    :: brnAs         ! coefficient for computing NNS rate
REAL(KIND=double), INTENT(out)    :: crnAs         ! coefficient for computing NNS rate
REAL(KIND=double), INTENT(out)    :: aenAs         ! coefficient for computing NNS fractional eneergy transfer
REAL(KIND=double), INTENT(out)    :: benAs         ! coefficient for computing NNS fractional eneergy transfer
REAL(KIND=double), INTENT(out)    :: cenAs         ! coefficient for computing NNS fractional eneergy transfer

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

LOGICAL                           :: sctn_off

INTEGER                           :: kp            ! outcoming neutrino energy zone index
INTEGER                           :: istat         ! allocation status

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

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0_ot     ! zero moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1_ot     ! first moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0d_ot    ! d(scant0)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1d_ot    ! d(scant1)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0t_ot    ! d(scant0)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1t_ot    ! d(scant1)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0y_ot    ! d(scant0)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1y_ot    ! d(scant1)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0_in     ! zero moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1_in     ! first moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0d_in    ! d(scant0)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1d_in    ! d(scant1)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0t_in    ! d(scant0)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1t_in    ! d(scant1)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0y_in    ! d(scant0)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1y_in    ! d(scant1)/d(ye)

EXTERNAL fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for sctednAs arrays')
 2001 FORMAT (' Deallocation problem for sctednAs arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  No scattering if nnugp(n) = 0, iscat = 0  or isctnA = 0
!  No scattering if rho outside prescribed boundaries
!-----------------------------------------------------------------------

sctn_off           = .false.
IF ( nnugp(n) == 0  .or.  iscat == 0  .or.  isctnA == 0 ) sctn_off = .true.
IF ( n < 3   .and.  rho < rhosctnAemn ) sctn_off = .true.
IF ( n < 3   .and.  rho > rhosctnAemx ) sctn_off = .true.
IF ( n == 3  .and.  rho < rhosctnAtmn ) sctn_off = .true.
IF ( n == 3  .and.  rho > rhosctnAtmx ) sctn_off = .true.

IF ( sctn_off ) THEN

  rmdnAs           = zero
  rmdnAs0          = zero
  rmdnAs1          = zero
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

  arnAs            = zero
  brnAs            = zero
  crnAs            = zero
  aenAs            = zero
  benAs            = zero
  cenAs            = zero

  RETURN

END IF

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (scatn0_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_ot '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_ot '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0d_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1d_ot(nez), STAT = istat) 
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0t_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1t_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0y_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1y_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_ot'; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scatn0_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_in '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_in '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0d_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1d_in(nez), STAT = istat) 
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0t_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1t_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0y_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1y_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_in'; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize scattering functions
!-----------------------------------------------------------------------

scatn0_ot          = zero
scatn1_ot          = zero
scatn0d_ot         = zero
scatn1d_ot         = zero
scatn0t_ot         = zero
scatn1t_ot         = zero
scatn0y_ot         = zero
scatn1y_ot         = zero

scatn0_in          = zero
scatn1_in          = zero
scatn0d_in         = zero
scatn1d_in         = zero
scatn0t_in         = zero
scatn1t_in         = zero
scatn0y_in         = zero
scatn1y_in         = zero

!-----------------------------------------------------------------------
!  Compute variables related to n-type neutrinos.
!
!  tmev       : t (MeV)
!  enuin      : (in beam neutrino energy)/kT
!  enuout     : (out beam neutrino energy)/kT
!-----------------------------------------------------------------------

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

CALL sctnAkrnl( n, j, ij_ray, ik_ray, k, rho, t, ye, scatn0_ot, scatn1_ot, &
& scatn0d_ot, scatn1d_ot, scatn0t_ot, scatn1t_ot, scatn0y_ot, scatn1y_ot, &
& scatn0_in, scatn1_in, scatn0d_in, scatn1d_in, scatn0t_in, scatn1t_in, &
& scatn0y_in, scatn1y_in, nez )

!-----------------------------------------------------------------------
!                   ||||| Begin loop over kp |||||
!-----------------------------------------------------------------------

DO kp = 1,nnugp(n)

  w2dw             = 2.d0 * pi * unu(j,kp) * unu(j,kp) * dunu(j,kp)
  w3dw             = unu(j,kp) * w2dw
  x0               = psi0(j,kp,n)
  x1               = psi1(j,kp,n)

  f0ot             = scatn0_ot(kp)
  f0in             = scatn0_in(kp)
  f1in             = zero
  f1ot             = zero

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

END DO

rmdnAs             =  zt
rmdnAs0            =  zt
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

arnAs              =  aqrt
brnAs              =  bqrt
crnAs              =  zero                                                       
aenAs              =  aqet
benAs              =  bqet
cenAs              =  zero                                                       

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (scatn0_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_ot '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_ot '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0d_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1d_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0t_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1t_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0y_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1y_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_ot'; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (scatn0_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_in '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_in '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0d_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1d_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0t_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1t_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0y_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1y_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_in'; WRITE (nlog,2001) var_name; END IF

RETURN 
END SUBROUTINE sctednAs
