SUBROUTINE editma( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editma
!    Module:       editma
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/22/04
!
!    Purpose:
!      To edit masses and mass averaged data.
!
!    Subprograms call:
!  date_and_time_print : prints date and time
!
!    Input arguments:
!  jr_min     : inner radial zone of region for which configuration edit is to be made.
!  jr_max     : outer radial zone of region for which configuration edit is to be made.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!  prnttest   : true  - test to see IF printing criteria is satisfied.
!               false - bypass printing criteria test.
!  nprint     : unit number of print file.
!  nedma(i)   : editc counter for data set i.
!  intdma(i)  : number of cycles between edits of data set i.
!  idxdma(i)  : edit jr_min, jr_max, and every idxedc(i) radial zone between them for data set i.
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module, nu_dist_module, shock_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nez, nez2, nnu
USE numerical_module, ONLY : zero, half, epsilon, ncoef, third, frpi, csqinv
USE physcnst_module, ONLY : pi, h, cvel, g, rmu, msolar

USE edit_module, ONLY : prnttest, nedma, intdma, idxema, nprint, nlog, head
USE eos_snc_x_module, ONLY : aesv
USE mdl_cnfg_module, ONLY : r, rho, t, ye, dmrst, rstmss, dmgrv, grvmss, &
& p0, rho0, gamgr, wgr
USE nu_dist_module, ONLY : unu, dunu, psi0, apnu, stress_x
USE nu_energy_grid_module, ONLY : nnugp
USE shock_module, ONLY : pq_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

LOGICAL                          :: first = .true.
LOGICAL                          :: lprint1
LOGICAL                          :: lprint2
LOGICAL                          :: lprint3

INTEGER                          :: j             ! radial zone index
INTEGER                          :: jm1           ! j - 1
INTEGER                          :: jp1           ! j + 1
INTEGER                          :: jv            ! do index
INTEGER                          :: k             ! do index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: n_time        ! used for data-and-time
INTEGER                          :: istat         ! allocation status

REAL(KIND=double), PARAMETER     :: fthird = 4.d0/3.d0
REAL(KIND=double)                :: tpigc2        ! 2*pi*g/c**2
REAL(KIND=double)                :: fcoef         ! coefficient for computing neutrino number flux

REAL(KIND=double)                :: yenu          ! electron neutrino fraction
REAL(KIND=double)                :: yanu          ! electron antineutrino fraction
REAL(KIND=double)                :: yxnu          ! x-neutrino fraction
REAL(KIND=double)                :: yl            ! lepton fraaction

REAL(KIND=double)                :: xn            ! neutron mass fraciton
REAL(KIND=double)                :: xp            ! proton mass fraciton
REAL(KIND=double)                :: xhe           ! helium mass fraction
REAL(KIND=double)                :: xnuc          ! heavy nuclei mass fraction

REAL(KIND=double)                :: rj2           ! r(j) * r(j)

REAL(KIND=double)                :: dratio        ! current to initial density ratio
REAL(KIND=double)                :: pneut         ! initial pressure scaled to current density by 4/3 law
REAL(KIND=double)                :: pdfct         ! p/pneut
REAL(KIND=double)                :: pndfct        ! ( p + p_neutrino )/pneut
REAL(KIND=double)                :: prhoc2        ! p0/rho0^4/3/c^{2}
REAL(KIND=double)                :: expa          ! argument of exponential
REAL(KIND=double)                :: pneutr        ! extrapolated pressure
REAL(KIND=double)                :: pdfcr         ! pressure ratio
REAL(KIND=double)                :: pndfcr        ! pressure ratio

REAL(KIND=double)                :: strsnn        ! neutrino stress (nonrelativistic)
REAL(KIND=double)                :: strsnr        ! neutrino stress (relativistic)
REAL(KIND=double)                :: strspn        ! pressure stress (nonrelativistic)
REAL(KIND=double)                :: strspr        ! pressure stress (relativistic)
REAL(KIND=double)                :: strsgn        ! gravitational stress (nonrelativistic)
REAL(KIND=double)                :: strsgr        ! gravitational stress (relativistic)

REAL(KIND=double)                :: pmsgs         ! ratio of pressure to gravitaional stress
REAL(KIND=double)                :: pmnsgs        ! ratio of pressure + neutrino to gravitaional stress
REAL(KIND=double)                :: pmsgsr        ! ratio of pressure to gravitaional stress (relativistic)
REAL(KIND=double)                :: ptsgsr        ! ratio of pressure + neutrino to gravitaional stress (relativistic)

REAL(KIND=double)                :: rhorm         ! enclosed mass weighted density
REAL(KIND=double)                :: trm           ! enclosed mass weighted temperature
REAL(KIND=double)                :: xnrm          ! enclosed mass weighted neutron fraction
REAL(KIND=double)                :: xprm          ! enclosed mass weighted proton fraction
REAL(KIND=double)                :: xherm         ! enclosed mass weighted helium fraction
REAL(KIND=double)                :: xarm          ! enclosed mass weighted heavy nucleus fraction
REAL(KIND=double)                :: yerm          ! enclosed mass weighted electron fraction
REAL(KIND=double)                :: ylrm          ! enclosed mass weighted lepton fraction
REAL(KIND=double)                :: yenurm        ! enclosed mass weighted electron neutrino fraction
REAL(KIND=double)                :: yanurm        ! enclosed mass weighted electron antineutrino fraction
REAL(KIND=double)                :: yxnurm        ! enclosed mass weighted x neutrino fraction

REAL(KIND=double)                :: pdfctm        ! enclosed mass weighted pressure deficit
REAL(KIND=double)                :: pndfcm        ! enclosed mass weighted material + neutrino pressure deficit
REAL(KIND=double)                :: pdfctr        ! enclosed mass weighted pressure deficit (relativistic)
REAL(KIND=double)                :: pndctr        ! enclosed mass weighted material + neutrino pressure deficit (relativistic)

REAL(KIND=double)                :: pmgsa         ! enclosed mass weighted ratio of pressure to gravitaional stress
REAL(KIND=double)                :: ptgsa         ! enclosed mass weighted ratio of pressure + neutrino to gravitaional stress
REAL(KIND=double)                :: pmgsra        ! enclosed mass weighted ratio of pressure to gravitaional stress (rel)
REAL(KIND=double)                :: ptgsra        ! enclosed mass weighted ratio of pressure + neutrino to gravitaional stress (rel)

REAL(KIND=double)                :: srm           ! enclosed mass weighted entropy

REAL(KIND=double)                :: dmrsts        ! zone rest mass (solar masses)
REAL(KIND=double)                :: rstmsss       ! enclosed rest mass (solar masses)
REAL(KIND=double)                :: dmgrvp        ! zone gravitational mass (solar masses)
REAL(KIND=double)                :: grvmssp       ! enclosed gravitational mass (solar masses)

REAL(KIND=double)                :: fexp          ! exponential function

EXTERNAL fexp

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rnnu

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: ac
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: dk

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1a)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
  101 FORMAT (46x,'Masses and Mass averaged data')
  103 FORMAT (45x,31('-')/)
  105 FORMAT ('   j     rho         T        xn          xp         xhe&
&       xa         ye         yl        yenu       yanu       yxnu'/)
  107 FORMAT (1x,i4,11(1pe11.3))
  201 FORMAT ('   j  p-deficit  p+n dfct   p gr dfct p+n gr dfct pms/gsnr &
& pmns/gs nr  pms/gs gr pmns/gs gr    s ave'/)
  203 FORMAT (1x,i4,9(1pe11.3))
  301 FORMAT ('   j    rst mss       grv mss    d rst mss    d grv mss'/)
  303 FORMAT (1x,i4,2(1pe14.6),2(1pe12.4))
 1001 FORMAT (' Allocation problem for array ',a10,' in editma')
 2001 FORMAT (' Deallocation problem for array ',a10,' in editma')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rnnu(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnnu      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ac(20,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ac        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dk(nx,20), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dk        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rnnu               = zero
ac                 = zero
dk                 = zero

!-----------------------------------------------------------------------
!  Calculate quantities needed for this edit
!-----------------------------------------------------------------------

IF ( .not. prnttest         .or.  &
&    nedma(1) >= intdma(1)  .or.  &
&    nedma(2) >= intdma(2)  .or.  &
&    nedma(3) >= intdma(3) ) THEN
  CALL w_cal ( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, wgr, nx )
  wgr(jr_max+1)      = wgr(jr_max)
  CALL nu_U( jr_min, jr_max, rho, nx, nnu )
  DO n = 1,nnu
    CALL flux( jr_min, jr_max, n )
  END DO
END IF

!**********************************************************************!
!                                                                      !
!                  Mass averaged data                       part 1     !
!                                                                      !
!----------------------------------------------------------------------!
! ac(1,j)     'rho'        Mass averaged value of the density          !
!                           interior to zone j (g/cm**3)               !
! ac(2,j)     'T'          Mass averaged value of the temperature      !
!                           interior to zone j (k)                     !
! ac(3,j)     'xn'         Mass averaged value of the neutron          !
!                           fraction interior to zone j                !
! ac(4,j)     'xp'         Mass averaged value of the proton           !
!                           fraction interior to zone j                !
! ac(5,j)     'xhe'        Mass averaged value of the helium           !
!                           fraction interior to zone j                !
! ac(6,j)     'xa'         Mass averaged value of the heavy nucleus    !
!                           fraction interior to zone j                !
! ac(7,j)     'ye'         Mass averaged value of the electron         !
!                           fraction interior to zone j                !
! ac(8,j)     'yl'         Mass averaged value of the lepton           !
!                           fraction interior to zone j                !
! ac(9,j)     'yenu'       Mass averaged value of the e-neutrino       !
!                           fraction interior to zone j                !
! ac(10,j)    'yanu'       Mass averaged value of the e-antineutrino   !
!                           fraction interior to zone j                !
! ac(11,j)    'yxnu'       Mass averaged value of the t-neutrino       !
!                           fraction interior to zone j                !
!----------------------------------------------------------------------!
!        Print header                                                  !
!----------------------------------------------------------------------!

lprint1            = .true.
n_time             = nprint

IF ( prnttest ) THEN
  nedma(1)         = nedma(1) + 1
  IF ( nedma(1) < intdma(1) ) lprint1 = .false.
  IF ( lprint1 ) nedma(1) = 0
END IF ! prnttest

IF ( lprint1 ) THEN
  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,101)
  WRITE (nprint,103)
  WRITE (nprint,105)
END IF ! lprint1

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  tpigc2           = 2.d+00 * pi * g /( cvel**2 )
  fcoef            = 4.d+00 * pi * ( cvel/3.d+00 )/( h * cvel )**3
  first            = .false.
END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rhorm              = zero
trm                = zero
xnrm               = zero
xprm               = zero
xherm              = zero
xarm               = zero
yerm               = zero
ylrm               = zero
yenurm             = zero
yanurm             = zero
yxnurm             = zero

pdfctm             = zero
pndfcm             = zero
pdfctr             = zero
pndctr             = zero
pmgsa              = zero
ptgsa              = zero
pmgsra             = zero
ptgsra             = zero
srm                = zero

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!        Begin radial loop.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Neutrino and lepton fractions
!-----------------------------------------------------------------------

  DO n = 1,nnu
    rnnu(n)        = zero
    IF ( nnugp(n) /= 0 ) THEN
    rnnu(n)        = ( ncoef/rho(j) ) * SUM( unu(j,:) * unu(j,:) * dunu(j,:) * psi0(j,:,n) )
    END IF ! nnugp(n) /= 0
  END DO
  yenu             = rnnu(1) * rmu
  yanu             = rnnu(2) * rmu
  yxnu             = rnnu(3) * rmu
  IF ( yenu < epsilon ) yenu = epsilon
  IF ( yanu < epsilon ) yanu = epsilon
  IF ( yxnu < epsilon ) yxnu = epsilon
  yl               = ye(j) + yenu - yanu

!-----------------------------------------------------------------------
!  Pressure deficits
!-----------------------------------------------------------------------

  dratio           = rho(j)/( rho0(j) + rho(jr_max) * 1.d-6 )
  pneut            = p0(j) * dratio**fthird
  pdfct            = aesv(j,1,ij_ray,ik_ray)/( pneut + 1.d-10 * aesv(j,1,ij_ray,ik_ray) )
  pndfct           = ( aesv(j,1,ij_ray,ik_ray) + apnu(j) )/( pneut + 1.d-10 * aesv(j,1,ij_ray,ik_ray) )
  prhoc2           = ( p0(j)/( rho0(j)**fthird ) ) * csqinv
  expa             = 2.78d0 * prhoc2 * 3.d0 * ( rho(j)**third - rho0(j)**third )
  pneutr           = p0(j) * dratio**(fthird) * fexp(expa)
  pdfcr            = aesv(j,1,ij_ray,ik_ray)/( pneutr + 1.d-10 * aesv(j,1,ij_ray,ik_ray) )
  pndfcr           = ( aesv(j,1,ij_ray,ik_ray) + apnu(j) )/( pneutr + 1.d-10*aesv(j,1,ij_ray,ik_ray) )

  jp1              = j + 1
  jm1              = j - 1
  rj2              = r(j) * r(j)

!-----------------------------------------------------------------------
!  Neutrino stress (nonrelativistic)
!-----------------------------------------------------------------------

strsnn             = zero
DO n = 1,nnu
  strsnn           = strsnn + stress_x(j,n,ij_ray,ik_ray)
END DO

!-----------------------------------------------------------------------
!  Neutrino stress (relativistic)
!-----------------------------------------------------------------------

  strsnr           = strsnn * gamgr(j)/( half * ( wgr(jp1) + wgr(j) ) )

!-----------------------------------------------------------------------
!  Pressure stress (nonrelativistic)
!-----------------------------------------------------------------------

  strspn           = -frpi * rj2 * ( ( ( aesv(jp1,1,ij_ray,ik_ray) + pq_x(jp1,ij_ray,ik_ray) ) &
&                  - ( aesv(j,1,ij_ray,ik_ray) + pq_x(j,ij_ray,ik_ray) ) )/( half * ( dmrst(jp1) + dmrst(j) ) ) )

!-----------------------------------------------------------------------
!  Pressure stress (relativistic)
!-----------------------------------------------------------------------

  strspr           = -frpi * r(j) * r(j) * ( ( ( aesv(jp1,1,ij_ray,ik_ray) + pq_x(jp1,ij_ray,ik_ray) ) &
&                  - ( aesv(j,1,ij_ray,ik_ray) + pq_x(j,ij_ray,ik_ray) ) ) * gamgr(j) &
&                  /( 0.25d0*( wgr(jp1) + wgr(j) ) * ( dmrst(jp1) + dmrst(j) ) ) )

!-----------------------------------------------------------------------
!  Gravitational stress (nonrelativistic)
!-----------------------------------------------------------------------

  strsgn           = rstmss(j) * g/rj2

!-----------------------------------------------------------------------
!  Gravitational stress (relativistic)
!-----------------------------------------------------------------------

  strsgr           = g * grvmss(j)/rj2 + tpigc2 * ( aesv(jp1,1,ij_ray,ik_ray) &
&                  + pq_x(jp1,ij_ray,ik_ray) + apnu(jp1)                      &
&                  + aesv(j,1,ij_ray,ik_ray) + pq_x(j,ij_ray,ik_ray) + apnu(j) ) * r(j)

!-----------------------------------------------------------------------
!  Stress ratios
!-----------------------------------------------------------------------

  pmsgs            = strspn/strsgn
  pmnsgs           = ( strspn + strsnn )/strsgn
  pmsgsr           = strspr/strsgr
  ptsgsr           = ( strspr + strsnr )/strsgr

!-----------------------------------------------------------------------
!  Mass weighting
!-----------------------------------------------------------------------

  rhorm            = rhorm + rho(j) * dmrst(j)
  trm              = trm   + t(j)   * dmrst(j)
  xn               = aesv(j,7,ij_ray,ik_ray)
  xnrm             = xnrm  + xn     * dmrst(j)
  xp               = aesv(j,8,ij_ray,ik_ray)
  xprm             = xprm  + xp     * dmrst(j)
  xnuc             = aesv(j,9,ij_ray,ik_ray)
  xarm             = xarm  + xnuc   * dmrst(j)
  xhe              = dmax1( 1.0d+00 - xn - xp - xnuc, zero )
  xherm            = xherm + xhe    * dmrst(j)

  yerm             = yerm   + ye(j)  * dmrst(j)
  ylrm             = ylrm   + yl     * dmrst(j)
  yenurm           = yenurm + yenu   * dmrst(j)
  yanurm           = yanurm + yanu   * dmrst(j)
  yxnurm           = yxnurm + yxnu   * dmrst(j)

  pdfctm           = pdfctm + pdfct  * dmrst(j)
  pndfcm           = pndfcm + pndfct * dmrst(j)
  pdfctr           = pdfctr + pdfcr  * dmrst(j)
  pndctr           = pndctr + pndfcr * dmrst(j)
  pmgsa            = pmgsa  + pmsgs  * dmrst(j)
  ptgsa            = ptgsa  + pmnsgs * dmrst(j)
  pmgsra           = pmgsra + pmsgsr * dmrst(j)
  ptgsra           = ptgsra + ptsgsr * dmrst(j)

  srm              = srm    + aesv(j,3,ij_ray,ik_ray)  * dmrst(j)

!-----------------------------------------------------------------------
!  Mass averaging
!-----------------------------------------------------------------------

  ac(1,j)          = rhorm  /  rstmss(j)
  ac(2,j)          = trm    /  rstmss(j)
  ac(3,j)          = xnrm   /  rstmss(j)
  ac(4,j)          = xprm   /  rstmss(j)
  ac(5,j)          = xherm  /  rstmss(j)
  ac(6,j)          = xarm   /  rstmss(j)
  ac(7,j)          = yerm   /  rstmss(j)
  ac(8,j)          = ylrm   /  rstmss(j)
  ac(9,j)          = yenurm /  rstmss(j)
  ac(10,j)         = yanurm /  rstmss(j)
  ac(11,j)         = yxnurm /  rstmss(j)

  dk(j,1)          = pdfctm /  rstmss(j)
  dk(j,2)          = pndfcm /  rstmss(j)
  dk(j,3)          = pdfctr /  rstmss(j)
  dk(j,4)          = pndctr /  rstmss(j)
  dk(j,5)          = pmgsa  /  rstmss(j)
  dk(j,6)          = ptgsa  /  rstmss(j)
  dk(j,7)          = pmgsra /  rstmss(j)
  dk(j,8)          = ptgsra /  rstmss(j)
  dk(j,9)          = srm    /  rstmss(j)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!        End radial loop.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END DO

!-----------------------------------------------------------------------
!  Print mass averages
!-----------------------------------------------------------------------

IF ( lprint1 ) THEN
  DO jv = jr_min,jr_max,idxema(1)
    j              = jr_max - jv + jr_min
    WRITE (nprint,107) j,(ac(k,j),k=1,11)
  END DO
END IF ! lprint1

!**********************************************************************!
!                                                                      !
!                  Mass averaged data                       part 2     !
!                                                                      !
!----------------------------------------------------------------------!
! dk(j,1)     'p deficit'  Mass averaged value of the pressure         !
!                           deficit, i.e., the ratio of the actual     !
!                           pressure to the pressure required for      !
!                           neutral stability, interior to zone j      !
! dk(j,2)     'p+n dfct'   Mass averaged value of the pressure         !
!                           deficit, including the neutrino isotropic  !
!                           pressure, interior to zone j               !
! dk(j,3)     'p gr dfct'  Mmass averaged value of the pressure        !
!                           deficit interior to zone j with a general  !
!                           relativistic correction to the neutral     !
!                           pressure support added                     !
! dk(j,4)    'p+n gr dfct' Mass averaged value of the pressure         !
!                           deficit, including the neutrino isotropic  !
!                           pressure, interior to zone j with a gener- !
!                           al relativistic correction to the neutral  !
!                           pressure support added                     !
! pmsgs       'pms/gs nr'  Mass averaged value interior to j of the    !
!                           ratio of the material pressure stress to   !
!                           the gravitational stress, computed         !
!                           nonrelativistically                        !
! pmnsgs      'pmns/gs nr' Mass averaged value interior to j of the    !
!                           ratio of the material plus neutrino        !
!                           pressure stress to the gravitational       !
!                           stress, computed nonrelativistically       !
! pmsgsr      'pms/gs gr'  Mass averaged value interior to j of the    !
!                           ratio of the material pressure stress to   !
!                           the gravitational stress, computed relati- !
!                           vistically                                 !
! ptsgsr      'pmns/gs gr' Mass averaged value interior to j of the    !
!                           ratio of the material plus neutrino        !
!                           pressure stress to the gravitational       !
!                           stress, computed relativistically          !
! dk(j,9)     's ave'      Mass averaged value interior to j of the    !
!                           matter entropy                             !
!----------------------------------------------------------------------!
!        Print header                                                  !
!----------------------------------------------------------------------!

lprint2            = .true.
 
IF ( prnttest ) THEN
  nedma(2)         = nedma(2) + 1
  IF ( nedma(2) < intdma(2) ) lprint2 = .false.
  IF ( lprint2 ) nedma(2) = 0
END IF ! prnttest

IF ( lprint2 ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,101)
  WRITE (nprint,103)
  WRITE (nprint,201)

!-----------------------------------------------------------------------
!  Print mass averages
!-----------------------------------------------------------------------

  DO jv = jr_min,jr_max,idxema(2)
    j                = jr_max - jv + jr_min
    WRITE (nprint,203) j,(dk(j,k),k=1,9)
  END DO

END IF

!**********************************************************************!
!                                                                      !
!                  Mass averaged data                       part 3     !
!                                                                      !
!----------------------------------------------------------------------!
! rstmsss     'rst mss'    Rest mass enclosed by zone j (solar masses) !
! grvmssp     'grv mss'    Gravitational mass enclosed by zone j       !
!                           (solar masses)                             !
! dmrsts      'd rst mss'  Rest mass of enclosed bu zone boundaries    !
!                           j and j-1 (solar masses)                   !
! dmgrvp      'd grv mss'  Gravitational mass of enclosed bu zone      !
!                           boundaries j and j-1 (solar masses)        !
!----------------------------------------------------------------------!
!        Print header                                                  !
!----------------------------------------------------------------------!

lprint3            = .true.

IF ( prnttest ) THEN
  nedma(3)         = nedma(3) + 1
  IF ( nedma(3) < intdma(3) ) lprint3 = .false.
  IF ( lprint3 ) nedma(3) = 0
END IF ! prnttest

IF ( lprint3 ) THEN 
  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,101)
  WRITE (nprint,103)
  WRITE (nprint,301)

!-----------------------------------------------------------------------
!  Print masses
!-----------------------------------------------------------------------

  DO jv = jr_min,jr_max,idxema(2)
    j                = jr_max - jv + jr_min

    rstmsss          = rstmss(j)/msolar
    grvmssp          = grvmss(j)/msolar
    dmrsts           = dmrst(j)/msolar
    dmgrvp           = dmgrv(j)/msolar

    WRITE (nprint,303) j,rstmsss,grvmssp,dmrsts,dmgrvp

  END DO
END IF ! lprint3

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (rnnu, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnnu      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ac, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ac        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dk, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dk        '; WRITE (nlog,2001) var_name; END IF

RETURN
END
