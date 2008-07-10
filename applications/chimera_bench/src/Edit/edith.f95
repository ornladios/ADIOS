SUBROUTINE edith( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         edith
!    Module:       edith
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/18/96
!
!    Purpose:
!      To edit hydrodynamic data.
!
!    Subprograms called:
!      date_and_time_print
!
!    Input arguments:
!  jr_min     : inner radial zone of region for which configuration edit
!                is to be made.
!  jr_max     : outer radial zone of region for which configuration edit
!                is to be made.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!  prnttest   : true  - test to see IF printing criteria is satisfied.
!               false - bypass printing criteria test.
!  iprint     : = 0   - do not print to print file.
!               \= 0 - print to print file.
!  nprint     : unit number of print file.
!  iplot      : = 0  - do not print to plot file.
!               \= 0 - print to plot file.
!  nplot      : unit number of plot file.
!  nedh(i)    : editc counter for data set i.
!  intedh(i)  : number of cycles between edits of data set i.
!  idxedh(i)  : edit jr_min, jr_max, and every idxedc(i) radial zone between
!                them for data set i.
!
!    Include files:
!  array_module, array_module, numerical_module, physcnst_module
!  convect_module, cycle_module, edit_module, el_eos_module,
!  eos_m4c_module, eos_snc_x_module, incrmnt_module, mdl_cnfg_module,
!  nu_dist_module, nu_energy_grid_module, shock_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nez2, nnu
USE numerical_module, ONLY : zero, half, frpi, epsilon, csqinv, ncoef, &
& third, pi2, sxtnpi
USE physcnst_module, ONLY : pi, cvel, rmu, g, cm3fm3, msolar, kmev, ergmev

USE convect_module, ONLY : pcnvct
USE cycle_module, ONLY : ncycle, nrst, intnu_trns
USE edit_module, ONLY : nprint, nlog, prnttest, nedh, intedh, idxedh, &
& head, nustrss, pstrss, gstrss, rstrss
USE el_eos_module, ONLY : PPRESS, PU
USE eos_m4c_module, ONLY: BPRESS, BPROUT, BPRNUC, BPRALF, BUOUT, BUNUC, &
& BU, BUALFA
USE eos_bck_module, ONLY: pe, pd, ph, prad, ee, ed, eh, erad, d0, dtran, &
& bind=>b, theta, upack, xmstar
USE eos_snc_x_module, ONLY: nse, aesv, aesvd, aesvt, aesvy, gam1, gam2, gam3
USE incrmnt_module, ONLY: dtmpmn
USE mdl_cnfg_module, ONLY: u, r, rhor, rho, t, ye, dmgrv, grvmss, gamgr, &
& agr, wgr, dmrst, rstmss, rho0, p0
USE nu_dist_module, ONLY: rhonun, rhonur, apmnun, apmnur, apnun, apnur, &
& aunu, fluxnu, unu, dunu, psi0, apnu, stwt, stress_nu=>stress_x
USE nu_energy_grid_module, ONLY : nnugp
USE shock_module, ONLY : pq_x
USE t_cntrl_module, ONLY : dtnmh, dth, dtj, dtjr

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
CHARACTER (len=1)                :: sch, ldx, sfr

LOGICAL                          :: first = .true.

INTEGER                          :: istat         ! allocation status
INTEGER                          :: jp1           ! j + 1
INTEGER                          :: jm1           ! j - 1
INTEGER                          :: jv            ! do index
INTEGER                          :: j             ! radial do index
INTEGER                          :: k             ! do index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: n_time        ! used for data-and-time

REAL(KIND=double), PARAMETER     :: fthird = 4.d0/3.d0
REAL(KIND=double)                :: c4            ! c^{4}
REAL(KIND=double)                :: tpigc2        ! 2*pi*g/c^{2}
REAL(KIND=double)                :: kdcgs         ! ( fm3/cm3 )( gram/# nucleons )

REAL(KIND=double)                :: c_sound       ! speed of sound (cm s^{-1})
REAL(KIND=double)                :: rho10         ! rho * 10^{-10}
REAL(KIND=double)                :: r7            ! r * 10^{-7}
REAL(KIND=double)                :: u9            ! u * 10^{-9}
REAL(KIND=double)                :: rhor7         ! rho10 * r7^{3}
REAL(KIND=double)                :: u9dr7         ! u9/r7
REAL(KIND=double)                :: u92r7         ! u9^{2} * r7
REAL(KIND=double)                :: csjp1         ! c_sound for j + 1
REAL(KIND=double)                :: csnd2         ! c_sound^{2}
REAL(KIND=double)                :: udc           ! u/c_sound
REAL(KIND=double)                :: uff2          ! free fall speed squared
REAL(KIND=double)                :: u2uff2        ! u^{2}/uff2
REAL(KIND=double)                :: uduff         ! u/uff

REAL(KIND=double)                :: dmdt          ! mass accretion rate (msolar^{-1} s^{-1})

REAL(KIND=double)                :: pml           ! log of the current matter pressure
REAL(KIND=double)                :: prml          ! log of the matter pressure at the preceeding timestep
REAL(KIND=double)                :: pmnl          ! log of the current matter + neutrino pressure
REAL(KIND=double)                :: prmnl         ! log of the matter + neutrino pressure at the preceeding timestep
REAL(KIND=double)                :: rhol          ! log of the current density
REAL(KIND=double)                :: rhorl         ! log of the density at the preceeding timestep
REAL(KIND=double)                :: gammae        ! effective matter adiabatic index
REAL(KIND=double)                :: gnmmae        ! effective matter + neutrino adiabatic index

REAL(KIND=double)                :: edc2          ! ratio of material internal energy to rest mass energy
REAL(KIND=double)                :: unudc2        ! ratio of the neutrino (all types) energy to the rest mass energy
REAL(KIND=double)                :: ufnu          ! ratio of the neutrino (all types) energy to the rest mass energy
REAL(KIND=double)                :: rschdr        ! measure of relativistic effects
REAL(KIND=double)                :: u2dc2         ! measure of relativistic effects
REAL(KIND=double)                :: r3pmc         ! measure of relativistic effects
REAL(KIND=double)                :: pdrhoc        ! measure of relativistic effects
REAL(KIND=double)                :: gmdrm         ! ratio of the gravitational mass to the rest mass

REAL(KIND=double)                :: yl            ! lepton fraction
REAL(KIND=double)                :: yerm          ! mass weighted ye
REAL(KIND=double)                :: ye2rm         ! mass weighted ye^{2}
REAL(KIND=double)                :: ylrm          ! mass weighted yl
REAL(KIND=double)                :: yl2rm         ! mass weighted yl^{2}
REAL(KIND=double)                :: gam1rm        ! mass weighted adiabatic exponent
REAL(KIND=double)                :: gamerm        ! mass weighted effective adiabatic exponent
REAL(KIND=double)                :: gnmerm        ! mass weighted adiabatic exponent including neutrinos

REAL(KIND=double)                :: dratio        ! ratio of density to initial density
REAL(KIND=double)                :: pneut         ! scaled initial pressure
REAL(KIND=double)                :: prhoc2        ! p0/rho0^{1/3}/c^{2}
REAL(KIND=double)                :: expa          ! argument of exponential
REAL(KIND=double)                :: pneutr        ! scaled pressure

REAL(KIND=double)                :: dylch         ! mass weighted lepton fraction function
REAL(KIND=double)                :: ylch          ! enclosed mass weighted lepton fraction function
REAL(KIND=double)                :: ylchav        ! mass aveerage of ylch
REAL(KIND=double)                :: rmcht         ! Chandrasekhar mass
REAL(KIND=double)                :: rmch          ! Chandrasekhar mass

REAL(KIND=double)                :: cmpe          ! electron chemical potential
REAL(KIND=double)                :: eta           ! electron chemical potential/kT
REAL(KIND=double)                :: ylcht         ! mass weighted lepton fraction function
REAL(KIND=double)                :: ylchvt        ! mass aveerage of ylcht
REAL(KIND=double)                :: rmchs         ! Chandrasekhar mass
REAL(KIND=double)                :: rmchss        ! Chandrasekhar mass
REAL(KIND=double)                :: ro10          ! rho * 10^{-10}
REAL(KIND=double)                :: ror10         ! rhor * 10^{-10}
REAL(KIND=double)                :: dro10t        ! d(rho10)/dt

REAL(KIND=double)                :: dmj           ! percentage enclosed mass

REAL(KIND=double)                :: v             ! EOS variable
REAL(KIND=double)                :: vd            ! d(EOS variable)/drho
REAL(KIND=double)                :: vt            ! d(EOS variable)/dT
REAL(KIND=double)                :: vy            ! d(EOS variable)/dye

REAL(KIND=double)                :: kp            ! ( erg/cm3 ) / ( mev/fm3 )
REAL(KIND=double)                :: ku            ! ( # nucleons/gram )( erg/MeV )
REAL(KIND=double)                :: p_total       ! total pressure
REAL(KIND=double)                :: u_total       ! total energy
REAL(KIND=double)                :: pedp          ! electron pressure/total pressure
REAL(KIND=double)                :: pbdp          ! drip pressure/total pressure
REAL(KIND=double)                :: padp          ! heavy nucleus pressure/total pressure
REAL(KIND=double)                :: praddp        ! neutrino pressure/total pressure
REAL(KIND=double)                :: pcnvctdp      ! convective pressure/total pressure
REAL(KIND=double)                :: uedu          ! electron energy/total energy
REAL(KIND=double)                :: ubdu          ! drip energy/total energy
REAL(KIND=double)                :: uadu          ! heavy nucleus energy/total energy
REAL(KIND=double)                :: uraddu        ! neutrino energy/total energy

REAL(KIND=double)                :: ro0           ! EOS parameter
REAL(KIND=double)                :: rotran        ! EOS parameter

REAL(KIND=double)                :: rj2           ! r^{2}
REAL(KIND=double)                :: strsnn        ! neutrino stress (nonrelativistic)
REAL(KIND=double)                :: strssn        ! neutrino stress (relativistic)
REAL(KIND=double)                :: strspn        ! pressure stress (nonrelativistic)
REAL(KIND=double)                :: strspr        ! pressure stress (relativistic)
REAL(KIND=double)                :: strssp        ! isotropic neutrino stress (relativistic)
REAL(KIND=double)                :: strsgn        ! gravitational stress (nonrelativistic)
REAL(KIND=double)                :: strsgr        ! gravitational stress (relativistic)
REAL(KIND=double)                :: pmsgs         ! strspn/strsgn
REAL(KIND=double)                :: pmnsgs        ! ( strspn + strsnn )/strsgn
REAL(KIND=double)                :: pmsgsr        ! - strspr/strsgr
REAL(KIND=double)                :: ptsgsr        ! - ( strspr + strssn )/strsgr
REAL(KIND=double)                :: rtsgsr        ! rstrss(j)/gstrss(j)

REAL(KIND=double)                :: vs2           ! zone-centered sound speed squared
REAL(KIND=double)                :: dmass         ! zone-edged mass
REAL(KIND=double)                :: drho          ! density decrement
REAL(KIND=double)                :: dp            ! pressure decrement
REAL(KIND=double)                :: stbltyt       ! stability parameter
REAL(KIND=double)                :: stbltyl       ! stability parameter
REAL(KIND=double)                :: fg2           ! stability parameter
REAL(KIND=double)                :: fg2p          ! stability parameter
REAL(KIND=double)                :: stbltys       ! stability parameter
REAL(KIND=double)                :: stbltysf      ! stability parameter

REAL(KIND=double)                :: dtmev         ! temperature change due to convection
REAL(KIND=double)                :: denom         ! denominator
REAL(KIND=double)                :: dddsyp        ! convective parameter
REAL(KIND=double)                :: dddysp        ! convective parameter
REAL(KIND=double)                :: dlnddlnsyp    ! convective parameter
REAL(KIND=double)                :: dlnddlnysp    ! convective parameter
REAL(KIND=double)                :: drhobuoy      ! convective parameter

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rnnu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: y_nu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: yla     ! lepton fraction array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: wBV     ! Brunt-Vaisala frequency (positive if unstable) [s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: twBV    ! Brunt-Vaisala growth time (positive if unstable) [s]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: wBV_s   ! Brunt-Vaisala frequency due to entropy (positive if unstable) [s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: twBV_s  ! Brunt-Vaisala growth time due to entropy (positive if unstable) [s]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: wBV_yl  ! Brunt-Vaisala frequency due to lepton fraction (positive if unstable) [s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: twBV_yl ! Brunt-Vaisala growth time due to lepton fraction (positive if unstable) [s]

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: ac
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: dz

REAL(KIND=double)                :: fexp          ! exponential function

EXTERNAL fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 format (1x)
    3 format (1x,a128)
    5 format (1x,128('*'))
  101 format (51x,'Hydrodynamic data')
  103 format (50x,19('-')/)
  105 format ('   j    u9/r7   rho10*r7**3 u9**2*r7   u/vsound     u/uff&
&      dm/dt     gamma1     gamma2     gamma3     gammae     gnmmae'/)
  107 format (1x,i4,11(1pe11.3))
  201 format (38x,'Hydrodynamic data - Relativistic parameters')
  203 format (37x,45('-')/)
  205 format ('   j   emat/c2    enu/c2   u*fnu/d*c4   rsch/r      u2/c2&
&   4pir3p/mc2   p/rhoc2     gamma        a        gm/rm        w'/)
  207 format (1x,i4,11(1pe11.3))
  301 format (38x,'Hydrodynamic data - Masses and mass averages')
  303 format (37x,46('-')/)
  305 format ('   j    yeave     ye2ave      ylave      yl2ave   gamma1ave&
&  gammaeave  gnmmaeave   dm grv     dm rst      %dm         chi'/)
  307 format (1x,i4,11(1pe11.3))
  309 format (1x,'Chandrasekhar mass =',1pe10.3)
  311 format (1x,'Chandrasekhar mass =',1pe10.3,' t effects included')
  401 format (33x,'Hydrodynamic data - Pressure, energy, and derivatives')
  403 format (32x,55('-')/)
  405 format ('   j      p        dp/dd      dp/dt     dp/dye        u &
&      du/dd      du/dt     du/dye'/)
  407 format (1x,i4,10(1pe11.3))
  501 format (34x,'Hydrodynamic data - Partial pressures and energies')
  503 format (33x,52('-')/)
  505 format ('   j    pe/p     pb-trns/p  pa-cprs/p   prad/p    pcnvct/p&
&     ue/u     ub-trns/u  ua-cprs/u   urad/u'/)
  507 format (1x,i4,9(1pe11.3))
  601 format (39x,'Hydrodynamic data - Nuclear matter properties')
  603 format (38x,47('-')/)
  605 format ('   j    be/b       d-sat     d-tran    d-a/d-sat  packfrac    mstar/m'/)
  607 format (1x,i4,6(1pe11.3))
  701 format (38x,'Hydrodynamic data - Pressure and stress ratios')
  703 format (37x,48('-')/)
  705 format ('   j  nu strss  nu isostrss p deficit  p+n dfct   p gr dfct&
& p+n gr dfct pms/gs nr pmns/gs nr  pms/gs gr pmns/gs gr  rs/gs gr'/)
  707 format (1x,i4,11(1pe11.3))
  801 format (25x,'Hydrodynamic data - Hydro timesteps, Convective Stability and stress ratios')
  803 format (24x,77('-')/)
  805 format ('   j     dth(j)      dtj(j)      dtjr(j)    dT (hydro)    Thorn &
&    Ledoux     Schwarzs   neut finger   snu/snunr    sp/spnr     sg/sgnr  &
&    wBV         twBV       twBV_s     twBV_yl '/)
  807 format (1x,i4,5(1pe12.3),3(1pe11.3,a1),7(1pe12.3))
  901 format (40x,'Hydrodynamic data - Ledoux Convection')
  903 format (39x,39('-')/)
  905 format ('   j   aledoux    psclht      lmix      usound     ulmixl &
&     ulcnvct      dt         dye      pcnvct     scnvct   droye/dros'/)
  907 format (1x,i4,11(1pe11.3))
 1001 FORMAT (' Allocation problem for array ',a10,' in edith')
 2001 FORMAT (' Deallocation problem for array ',a10,' in edith')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rnnu(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnnu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (y_nu(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'y_nu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (yla(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yla       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wBV(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wBV       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (twBV(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'twBV      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wBV_s(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wBV_s     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (twBV_s(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'twBV_s    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wBV_yl(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wBV_yl    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (twBV_yl(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'twBV_yl   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ac(20,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ac        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dz(20,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dz        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rnnu               = zero
y_nu               = zero
ac                 = zero
dz                 = zero

IF ( first ) THEN
  first            = .false.
  c4               = cvel**4               ! c**4
  tpigc2           = 2.d0 * pi * g/cvel**2 ! 2.*pi*g/!**2
  kdcgs            = rmu/cm3fm3            ! ( fm3/cm3 )( gram/# nucleons )
  kp               = ergmev/cm3fm3         ! ( erg/cm3 ) / ( mev/fm3 )
  ku               = ergmev/rmu            ! ( # nucleons/gram )( erg/mev )
END IF ! first
n_time             = nprint

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                       part 1     !
!                                                                      !
!----------------------------------------------------------------------!
! u9dr7       'u9/r7'      (u(j)/1.e+9)/(r(j)/1.e+7); this quantity    !
!                           is constant for homologous collapse        !
! rhor7      'rh010*r7**3' (rho(j)/1.e+10) * (r(j)/1.e+7)**3; this     !
!                           quantity is constant when the density pro- !
!                           file is proportional to r**(-3)            !
! u92r7       'u9**2*r7'   (u(j)/1.e+9)**2 * (r(j)/1.e+7); this        !
!                           quantity is constant when the infall velo- !
!                           city is proportional to r**(-1/2), as it   !
!                           is approximately in the outer core during  !
!                           infall                                     !
! udc         'u/vsound'   Mach number of the flow                     !
! uduff       'u/uff'      Ratio of the absolute magnitude of the      !
!                           material velocity over the free-fall       !
!                           velocity                                   !
! dmdt        'dm/dt'      Mass accretion rate; 4*pi*r**2*rho*u        !
! gamma1      'gamma1'     First adiabatic exponent                    !
! gamma2      'gamma2'     Second adiabatic exponent                   !
! gamma3      'gamma3'     Third adiabatic exponent                    !
! gammae      'gammae'     Effective adiabatic exponent, i.e., the     !
!                           value of gamma1 as seen by a comoving      !
!                           fluid element                              !
! gnmmae      'gnmmae'     Effective adiabatic exponent including      !
!                           isotropic neutrino pressure                !
!----------------------------------------------------------------------!

IF ( prnttest ) THEN
  nedh(1)          = nedh(1) + 1
  IF ( nedh(1) < intedh(1)  .or.  intedh(1) < 0 ) GO TO 1100
  nedh(1)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,101)
WRITE (nprint,103)
WRITE (nprint,105)

c_sound            = zero

DO jv = jr_min,jr_max,idxedh(1)
  j                = jr_max - jv + jr_min

  jp1              = j + 1

  rho10            = 1.d-10 * rho(j)
  r7               = 1.d-7 * r(j)
  u9               = 1.d-9 * u(j)
  rhor7            = rho10 * r7 * r7 * r7
  u9dr7            = u9/r7
  u92r7            = u9 * u9 * r7

!-----------------------------------------------------------------------
!  Sound speed, free fall speed
!-----------------------------------------------------------------------

  dz(5,j)          = gam1(j,ij_ray,ik_ray)
  csjp1            = c_sound
  csnd2            = gam1(j,ij_ray,ik_ray) * aesv(j,1,ij_ray,ik_ray)/rho(j)
  IF ( csnd2 > zero ) THEN
    c_sound        = DSQRT(csnd2)
    udc            = DABS( u(j)/c_sound )
  ELSE
    c_sound        = zero
    udc            = zero
  END IF
  uff2             = 2.d0 * g * grvmss(j)/r(j)
  u2uff2           = u(j) * u(j)/uff2
  IF ( u2uff2 > zero ) THEN
    uduff          = DSQRT(u2uff2)
  ELSE
    uduff          = zero
  END IF

!-----------------------------------------------------------------------
!  Mass accretion rate
!-----------------------------------------------------------------------

  dmdt             = frpi * r(j)**2 * u(j) * half * ( rho(j) + rho(j+1) )/msolar

!-----------------------------------------------------------------------
!  Effective adiabatic exponent
!-----------------------------------------------------------------------

  IF ( rhonun(j) /= rhonur(j)  .and.  ncycle > intnu_trns + nrst + 1 ) THEN
    pml            = DLOG( DABS(apmnun(j))             + epsilon )
    prml           = DLOG( DABS(apmnur(j))             + epsilon )
    pmnl           = DLOG( DABS(apmnun(j)  + apnun(j)) + epsilon )
    prmnl          = DLOG( DABS(apmnur(j)  + apnur(j)) + epsilon )
    rhol           = DLOG( DABS(rhonun(j))             + epsilon )
    rhorl          = DLOG( DABS(rhonur(j))             + epsilon )
    IF ( rhol .ne. rhorl ) THEN
      gammae       = ( pml  - prml  )/( rhol - rhorl )
      gnmmae       = ( pmnl - prmnl )/( rhol - rhorl )
      dz(6 ,j)     = gammae
      dz(11,j)     = gnmmae
    END IF ! rhol ne rhorl
  ELSE
    dz(6 ,j)       = zero
    dz(11,j)       = zero
    gammae         = zero
    gnmmae         = zero
  END IF ! rhonun(j) ne rhonur(j)  and  ncycle < intnu_trns + nrst + 1

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

WRITE (nprint,107) j, u9dr7, rhor7, u92r7, udc, uduff, dmdt, &
& gam1(j,ij_ray,ik_ray), gam2(j,ij_ray,ik_ray), gam3(j,ij_ray,ik_ray), &
& gammae, gnmmae

END DO

!**********************************************************************!
!                  Hydrodynamic data.                      part 2      !
!                                                                      !
!        Relativistic parameters.                                      !
!----------------------------------------------------------------------!
! edc2        'emat/c2'    Ratio of material internal energy to rest   !
!                           mass energy in zone j - 1/2 : relative     !
!                           contribution of the mat internal energy to !
!                           the gravitational energy                   !
! unudc2      'enu/c2'     Ratio of the neutrino (all types) energy    !
!                           to the rest mass energy in zone j - 1/2 :  !
!                           relative contribution of the neutrino      !
!                           energy to the gravitational energy         !
! ufnu        'u*fnu/d*c4' (mat velocity)*(nu flux)/(rho*!**4)   :     !
!                           relative contribution of the neutrino flux !
!                           to the gravitational mass                  !
! rschdr      'rsch/r'     Ratio of the schwarzchild radius of the     !
!                           mass enclosed by zone j to the radius of   !
!                           zone j - a measure of relativistic         !
!                           effects                                    !
! u2dc2       'u2/c2'      (mat velocity)**2/!**2 ) - a measure of     !
!                           relativistic effects                       !
! r3pmc       '4pir3p/mc2' (4*pi*radius**3*pressure/(mass*!**2) - a    !
!                           measure of relativistic effects            !
! pdrhoc      'p/rhoc2'    pressure/(density*!**2) - a measure of      !
!                           relativistic effects                       !
! gamgr(j)    'gamma'      d(coordinate radius)/d(proper radius) of    !
!                           zone j                                     !
! agr(j)      'a'          d(proper time)/d(coordinate time) of zone   !
!                           j                                          !
! gmdrm       'gm/rm'      Ratio of the gravitational mass to the      !
!                           rest mass of zone j                        !
! wgr(j)      'w'          Relativistic enthalpy of zone j             !
!----------------------------------------------------------------------!

 1100 CONTINUE

IF ( prnttest ) THEN
  nedh(2)          = nedh(2) + 1
  IF ( nedh(2) .lt. intedh(2)  .or.  intedh(2) .le. 0 ) GO TO 2100
  nedh(2)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,201)
WRITE (nprint,203)
WRITE (nprint,205)

DO jv = jr_min,jr_max,idxedh(2)
  j                = jr_max - jv + jr_min

  jm1              = j - 1
  edc2             = aesv(j,2,ij_ray,ik_ray) * csqinv
  unudc2           = csqinv * aunu(j)
  ufnu             = half * ( u(j) + u(jm1) ) * SUM(fluxnu(j,:))/( c4 * rho(j) )
  rschdr           = 2.d0 * g * ( grvmss(j) * csqinv )/r(j)
  u2dc2            = u(j) * u(j) * csqinv
  r3pmc            = frpi * ( aesv(j,1,ij_ray,ik_ray)/grvmss(j) ) * ( r(j) * r(j) * csqinv * r(j) )
  pdrhoc           = ( aesv(j,1,ij_ray,ik_ray) * csqinv )/rho(j)
  gmdrm            = dmgrv(j)/dmrst(j)

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

  WRITE (nprint,207) j, edc2, unudc2, ufnu, rschdr, u2dc2, r3pmc, pdrhoc, &
&  gamgr(j), agr(j), gmdrm, wgr(j)

END DO

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                      part 3      !
!                                                                      !
!        Masses and mass averages.                                     !
!----------------------------------------------------------------------!
! dz(1,j)     'yeave'      Mass average of electron fraction interior  !
!                           to zone j                                  !
! dz(2,j)     'ye2ave'     Mass average of (electron fraction)**2      !
!                           interior to zone j                         !
! dz(3,j)     'ylave'      (Mass average of lepton fraction interior   !
!                           to zone j                                  !
! dz(4,j)     'yl2ave'     Mass average of (lepton fraction)**2        !
!                           interior to zone j                         !
! dz(5,j)     'gamma1ave'  Mass average of first adiabatic exponent    !
!                           interior to zone j                         !
! dz(6,j)     'gammaeave'  Mass average of effective adiabatic expo-   !
!                           nent interior to zone j                    !
! dz(j,7)     'gnmmaeave'  Mass average of effective adiabatic expo-   !
!                           nent (including the isotropic neutrino     !
!                           pressure) interior to zone j               !
! dmgrv(j)    'dm grv'     Gravitational mass of zone j (g)            !
! dmrst(j)    'dm rst'     Rest mass of zone j (g)                     !
! dmj         '%dm'        Percent of total core rest mass contained   !
!                           in zone j - 1/2                            !
! dz(10,j)    'chi'        Ratio of collapse rate to free fall rate    !
!----------------------------------------------------------------------!

 2100 CONTINUE

IF ( prnttest ) THEN
  nedh(3)          = nedh(3) + 1
  IF ( nedh(3) < intedh(3)  .or.  intedh(3) <= 0 ) GO TO 4100
  nedh(3)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,301)
WRITE (nprint,303)
WRITE (nprint,305)

rmch               = zero
rmchss             = zero
yerm               = zero
ye2rm              = zero
ylrm               = zero
yl2rm              = zero
gam1rm             = zero
gamerm             = zero
gnmerm             = zero
ylch               = zero
ylcht              = zero

DO j = jr_min,jr_max

  DO n = 1,nnu
    rnnu(n)        = zero

    IF ( nnugp(n) /= 0 ) THEN
      DO k = 1,nnugp(n)
        rnnu(n)    = rnnu(n) + ( ncoef/rho(j) ) * unu(j,k)**2 * dunu(j,k) * psi0(j,k,n)
      END DO
    END IF ! nnugp(n) ne 0
    y_nu(n)        = DMAX1( rnnu(n) * rmu, epsilon )

  END DO

  yl               = ye(j) + y_nu(1) - y_nu(2)

!-----------------------------------------------------------------------
!  Mass weighting
!-----------------------------------------------------------------------

  yerm             = yerm + ye(j)      * dmrst(j)
  dz(1,j)          = yerm/rstmss(j)

  ye2rm            = ye2rm + ye(j)**2  * dmrst(j)
  dz(2,j)          = ye2rm/rstmss(j)

  ylrm             = ylrm + yl         * dmrst(j)
  dz(3,j)          = ylrm/rstmss(j)

  yl2rm            = yl2rm + yl**2     * dmrst(j)
  dz(4,j)          = yl2rm/rstmss(j)

  gam1rm           = gam1rm + dz(5,j)  * dmrst(j)
  dz(7,j)          = gam1rm/rstmss(j)

  gamerm           = gamerm + dz(6,j)  * dmrst(j)
  dz(8,j)          = gamerm/rstmss(j)

  gnmerm           = gnmerm + dz(11,j) * dmrst(j)
  dz(9,j)          = gnmerm/rstmss(j)

!-----------------------------------------------------------------------
!  Pressure deficits
!-----------------------------------------------------------------------

  dratio           = rho(j)/( rho0(j) + rho(jr_max) * 1.d-6 )
  pneut            = p0(j) * dratio**fthird
  ac(1,j)          = aesv(j,1,ij_ray,ik_ray)/( pneut + 1.d-10 * aesv(j,1,ij_ray,ik_ray) )
  ac(2,j)          = ( aesv(j,1,ij_ray,ik_ray) + apnu(j) )/( pneut + 1.d-10 * aesv(j,1,ij_ray,ik_ray) )
  prhoc2           = ( p0(j)/( rho0(j)**fthird ) ) * csqinv
  expa             = 2.78d0 * prhoc2 * 3.d0 * ( rho(j)**third - rho0(j)**third )
  pneutr           = p0(j) * dratio**fthird * fexp(expa)
  ac(3,j)          = aesv(j,1,ij_ray,ik_ray)/( pneutr + 1.d-10 * aesv(j,1,ij_ray,ik_ray) )
  ac(4,j)          = ( aesv(j,1,ij_ray,ik_ray) + apnu(j) )/( pneutr + 1.d-10 * aesv(j,1,ij_ray,ik_ray) )

!-----------------------------------------------------------------------
!  Chandrasekhar masses
!-----------------------------------------------------------------------

  dylch            = ( (   ye(j)**fthird + 1.2599d0 * ( stwt(1) * y_nu(1)**fthird &
&                  + stwt(2) * y_nu(2)**fthird + stwt(3) * y_nu(3)**fthird ) )**1.5 )*dmrst(j)
  ylch             = ylch + dylch
  ylchav           = ylch/rstmss(j)
  rmcht            = 5.76d0 * ylchav
  IF ( rmcht < rstmss(j)/msolar ) rmch = rmcht
  cmpe             = aesv(j,6,ij_ray,ik_ray)
  eta              = cmpe/( kmev * t(j) )
  ylcht            = ylcht + dylch + ( pi2/( eta * eta ) ) * dmrst(j)
  ylchvt           = ylcht/rstmss(j)
  rmchs            = 5.76d0 * ylchvt
  IF ( rmchs < rstmss(j)/2.e+33 ) rmchss = rmchs
  ro10             = rho(j)/1.d+10
  ror10            = rhor(j)/1.d+10
  dro10t           = ( ro10 - ror10 )/( dtnmh + 1.d-24 )
  dz(10,j)         = 224.d0 * ( (ro10)**1.5 )/( dro10t + 1.d-10 )

END DO

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

DO jv = jr_min,jr_max,idxedh(3)
  j                = jr_max - jv + jr_min

  dmj              = 1.e+02 * dmrst(j)/rstmss(jr_max)
  WRITE (nprint,307) j, dz(1,j), dz(2,j), dz(3,j), dz(4,j), dz(7,j), dz(8,j), &
&  dz(9,j), dmgrv(j), dmrst(j), dmj, dz(10,j)

END DO

WRITE (nprint,309) rmch
WRITE (nprint,311) rmchss

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                      part 4      !
!                                                                      !
!	       Pressure, energy, and derivatives.                          !
!----------------------------------------------------------------------!
! aesv(j,1,ij_ray,ik_ray)   'p'      Mat pressure (dynes)              !
! aesvd(j,1,ij_ray,ik_ray)  'dp/dd'  d(mat pressure)/d(density         !
!                                     (dynes/g/cm**3)                  !
! aesvt(j,1,ij_ray,ik_ray)  'dp/dt'  d(mat pressure)/d(temperature     !
!                                     (dynes/k)                        !
! aesvy(j,1,ij_ray,ik_ray)  'dp/dr'  d(mat pressure)/d(ye) (dynes)     !
! aesv(j,2,ij_ray,ik_ray)   'u'      Mat energy (ergs/g)               !
! aesvd(j,2,ij_ray,ik_ray)  'du/dd'  d(mat energy)/d(density           !
!                                     (ergs/g/g/cm**3)                 !
! aesvt(j,2,ij_ray,ik_ray)  'du/dt'  d(mat energy)/d(temperature)      !
!                                     (ergs/g/k)                       !
! aesvy(j,2,ij_ray,ik_ray)  'du/dr'  d(mat energy)/d(ye) (ergs/g)      !
!----------------------------------------------------------------------!

 4100 CONTINUE

IF ( prnttest ) THEN
  nedh(4)          = nedh(4) + 1
  IF ( nedh(4) .lt. intedh(4)  .or.  intedh(4) .le. 0 ) GO TO 5100
  nedh(4)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,401)
WRITE (nprint,403)
WRITE (nprint,405)

DO jv = jr_min,jr_max,idxedh(4)
  j                = jr_max - jv + jr_min

  WRITE (nprint,407) j, aesv(j,1,ij_ray,ik_ray), aesvd(j,1,ij_ray,ik_ray), &
&  aesvt(j,1,ij_ray,ik_ray), aesvy(j,1,ij_ray,ik_ray), aesv(j,2,ij_ray,ik_ray), &
&  aesvd(j,2,ij_ray,ik_ray), aesvt(j,2,ij_ray,ik_ray), aesvy(j,2,ij_ray,ik_ray)

END DO

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                      part 5      !
!                                                                      !
!        Partial pressures and energies.                               !
!----------------------------------------------------------------------!
! pedp        'pe/p'       Ratio of electron-positron pressure to the  !
!                           total material pressure                   !
! pbdp        'pb-trns/p'  Ratio of the translation pressure of        !
!                           nuclei and dripped nucleons to the total   !
!                           material pressure                          !
! padp	      'pa-cprs/p'  Ratio of the pressure of nuclei other than  !
!                           translational to the total material        !
!                           pressure                                   !
! praddp      'prad/p'     Ratio of the photon pressure to the total   !
!                           material pressure                          !
! pcnvctdp    'pcnvct/p'   Ratio of the turbulent pressure to the      !
!                           total meterial pressure                    !
! uedu	      'ue/u'       Ratio of electron-positron energy to the    !
!                           total material energy                      !
! ubdu	      'ub-trns/u'  Ratio of the translational energy of        !
!                           nuclei and dripped nucleons to the total   !
!                           material energy                            !
! uadu        'ua-cprs/u'  Ratio of the energy of nuclei other than    !
!                           translational to the total meterial        !
!                           energy                                     !
! uraddu      'urad/u'     Ratio of the photon energy to the total     !
!                           material energy                            !
!----------------------------------------------------------------------!

5100 CONTINUE

IF ( prnttest ) THEN
nedh(5)     = nedh(5) + 1
IF ( nedh(5) < intedh(5)  .or.  intedh(5) <= 0 ) GO TO 6100
nedh(5)     = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,501)
WRITE (nprint,503)
WRITE (nprint,505)

DO jv = jr_min,jr_max,idxedh(5)
j                = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Equation of state components
!-----------------------------------------------------------------------

pd               = zero
ph               = zero
prad             = zero
BPROUT           = zero
BPRALF           = zero
BPRNUC           = zero
PPRESS           = zero
ed               = zero
eh               = zero
erad             = zero
BUOUT            = zero
BUALFA           = zero
BUNUC            = zero
PU               = zero

CALL eqstta_x( 1, j, ij_ray,ik_ray, rho(j), t(j), ye(j), v, vd, vt, vy )
p_total          = v
CALL eqstta_x( 2, j, ij_ray,ik_ray, rho(j), t(j), ye(j), v, vd, vt, vy )
u_total          = v

pedp             = kp * pe                      /( p_total + epsilon )
pbdp             = kp * ( pd + BPROUT + BPRALF )/( p_total + epsilon )
padp             = kp * ( ph + BPRNUC          )/( p_total + epsilon )
praddp           = kp * ( prad + PPRESS        )/( p_total + epsilon )
pcnvctdp         =        pcnvct(j)             /( p_total + epsilon )
uedu             = ku * ee       /( u_total + epsilon )
ubdu             = ku * ( ed + BUOUT + BUALFA  )/( u_total + epsilon )
uadu             = ku * ( eh + BUNUC           )/( u_total + epsilon ) 
uraddu           = ku * ( erad + PU            )/( u_total + epsilon )

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

WRITE (nprint,507) j, pedp, pbdp, padp, praddp, pcnvctdp, uedu, ubdu, &
& uadu, uraddu

END DO

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                      part 6      !
!                                                                      !
!        Nuclear matter properties.                                    !
!----------------------------------------------------------------------!
! bind        'be/b'       Binding energy of heavy nucleus (mev)       !
! d0          'd-sat'      Saturation density of dense phase, or       !
!                           density of nuclear matter (g/cm3)          !
! dtran	      'd-tran'     Transition density to nuclear matter        !
!                           (g/cm3)                                    !
! theta       'd-a/d-sat'  Ratio of the density of dense phase to the  !
!                           saturation density of the dense phase      !
! upack       'packfrac'   Packing fraction of dense phase             !
! xmstar      'mstar/m'    Ratio of effective mass to rest mass for    !
!                           heavy nuclei or nuclear matter             !
!----------------------------------------------------------------------!

6100 CONTINUE

IF ( prnttest ) THEN
nedh(6)          = nedh(6) + 1
IF ( nedh(6) < intedh(6)  .or.  intedh(6) <= 0 ) GO TO 7100
nedh(6)     = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,601)
WRITE (nprint,603)
WRITE (nprint,605)

DO jv = jr_min,jr_max,idxedh(6)
j                = jr_max - jv + jr_min

d0               = zero
dtran            = zero
bind             = zero
ro0              = zero
theta            = zero
upack            = zero
xmstar           = zero
IF ( nse(j,ij_ray,ik_ray) == 0 ) CALL eqstta_x( 12, j, ij_ray,ik_ray, &
& rho(j), t(j), ye(j), v, vd, vt, vy )
ro0              = kdcgs * d0
rotran           = kdcgs * dtran

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

WRITE (nprint,607) j, bind, ro0, rotran, theta, upack, xmstar

END DO

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                       part 7     !
!                                                                      !
!        Pressure and stress ratios.                                   !
!----------------------------------------------------------------------!
! strssn      'nu strss'   Force on matter due to neutrino flow        !
!                           (dynes/gm)                                 !
! strssp     'nu isostrss' Force on matter due to neutrinos, computed  !
!                           assuming neutrino isotropy in each zone    !
!                           (dynes/gm)                                 !
! ac(1,j)     'p deficit'  Pressure deficit; the ratio of the          !
!                           actual pressure to the pressure required   !
!                           for neutral stability; i.e.,               !
!                           current pressure                           !
!                           /( (initial pressure)*(density/initial     !
!                           density)**(4/3)                            !
! ac(2,j)     'p+n dfct'   Pressure deficit with the neutrino          !
!                           isotropic pressure included                !
! ac(3,j)     'p gr dfct'  Pressure deficit with general relativistic  !
!                           correction to pressure for neutral support !
!                           included                                   !
! ac(4,j)    'p+n gr dfct' Pressure deficit with general relativistic  !
!                           correction and neutrino isotropic pressure !
!                           included                                   !
! pmsgs       'pms/gs nr'  Ratio of the material pressure stress to    !
!                           the gravitational stress, computed         !
!                           nonrelativistically                        !
! pmnsgs      'pmns/gs nr' Ratio of the material plus neutrino         !
!                           pressure stress to the gravitational       !
!                           stress, computed nonrelativistically       !
! pmsgsr      'pms/gs gr'  Ratio of the material pressure stress to    !
!                           the gravitational stress, computed relati- !
!                           vistically                                 !
! ptsgsr      'pmns/gs gr' Ratio of the material plus neutrino         !
!                           pressure stress to the gravitational       !
!                           stress, computed relativistically          !
! rtsgsr      'rs/gs gr'   Ratio of the pressure term in the gravita-  !
!                           tional stress to tje total gravitational   !
!                           stress                                     !
!----------------------------------------------------------------------!

7100 CONTINUE

IF ( prnttest ) THEN
nedh(7)           = nedh(7) + 1
IF ( nedh(7) < intedh(7)  .or.  intedh(7) <= 0 ) GO TO 8100
nedh(7)           = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,701)
WRITE (nprint,703)
WRITE (nprint,705)

CALL stress( jr_min, jr_max, ij_ray, ik_ray )

DO jv = jr_min,jr_max,idxedh(7)
j                = jr_max - jv + jr_min

jp1              = j + 1
jm1              = j - 1
rj2              = r(j) * r(j)

!-----------------------------------------------------------------------
!  Neutrino stress (nonrelativistic)
!-----------------------------------------------------------------------

strsnn           = SUM( stress_nu(j,:,ij_ray,ik_ray) )                          

!-----------------------------------------------------------------------
!  Neutrino stress (relativistic)
!-----------------------------------------------------------------------

strssn           = agr(j) * nustrss(j)

!-----------------------------------------------------------------------
!  Pressure stress (nonrelativistic)
!-----------------------------------------------------------------------

IF ( j /= jr_max ) THEN
  strspn         = -frpi * rj2 * ( ( ( aesv(jp1,1,ij_ray,ik_ray) + pq_x(jp1,ij_ray,ik_ray)) &
&                -  aesv(j,1,ij_ray,ik_ray) - pq_x(j,ij_ray,ik_ray) )                       &
&                / ( half * ( dmrst(jp1) + dmrst(j) ) ) )
ELSE
  strspn         = -frpi * rj2 * ( aesv(jp1,1,ij_ray,ik_ray) - aesv(j,1,ij_ray,ik_ray)      &
&                - pq_x(j,ij_ray,ik_ray) )/dmrst(j)
END IF ! j ne jr_max

!-----------------------------------------------------------------------
!  Pressure stress (relativistic)
!-----------------------------------------------------------------------

strspr           = agr(j) * pstrss(j)

!-----------------------------------------------------------------------
!  Neutrino stress (neutrino distribution assumed locally isotropic)
!-----------------------------------------------------------------------

IF ( j /= jr_max ) THEN
  strssp         = agr(j) * ( -sxtnpi * rj2 * gamgr(j) * ( ( apnu(jp1) - apnu(j) ) &
&                  /( ( wgr(jp1) + wgr(j) ) * ( rstmss(jp1) - rstmss(jm1) ) ) ) )
ELSE
  strssp         = zero
END IF ! j ne jr_max

!-----------------------------------------------------------------------
!  Gravitational stress (nonrelativistic)
!-----------------------------------------------------------------------

strsgn           = rstmss(j) * g/rj2

!-----------------------------------------------------------------------
!  Gravitational stress (relativistic)
!-----------------------------------------------------------------------

strsgr           = agr(j) * ( gstrss(j) + rstrss(j) )

!-----------------------------------------------------------------------
!  Stress ratios
!-----------------------------------------------------------------------

pmsgs            = strspn/( strsgn + epsilon )
pmnsgs           = ( strspn + strsnn )/( strsgn + epsilon )
pmsgsr           = - strspr/( strsgr + epsilon )
ptsgsr           = - ( strspr + strssn )/( strsgr + epsilon )
rtsgsr           = rstrss(j)/( gstrss(j) + rstrss(j) + epsilon )
ac(5,j)          = strssn/( strsnn + epsilon )
ac(6,j)          = strspr/( strspn + epsilon )
ac(7,j)          = strsgr/( strsgn + epsilon )

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

WRITE (nprint,707) j, strssn, strssp, ac(1,j), ac(2,j), ac(3,j), ac(4,j), &
& pmsgs, pmnsgs, pmsgsr, ptsgsr, rtsgsr

END DO

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                      part 8      !
!                                                                      !
!        Hydro timesteps, convective stability, and stress ratios.     !
!----------------------------------------------------------------------!        
! dth(j)      'dth(j)'     Minimum hydro timestep criteria for zone    !        
!                           j (sec)                                    !        
! dtj(j)      'dtj(j)'     Hydro timestep of zone j (sec)              !        
! dtjr(j)     'dtjr(j)'    Hydro timestep of zone j for preceding      !        
!                           cycle (sec)                                !        
! dtmpmn(j,1,ij_ray,ik_ray)                                            !
!             'dt (hydro)' Matter temperature change in the last time  !
!                           cycle due to hydrodynamics (k)             !
! stbltyt     'Thorn'      "Thorn" criterion for convection            !        
!                           c_sound(j)*c_sound(j+1)*(rho(j+1) - rho(j))!        
!                           /( 0.5*(rstmss(j+1) - rstmss(j-1)) )       !        
!                           - (p(j+1) - p(j)                           !        
!                           /( 0.5*(rstmss(j+1) - rstmss(j-1))         !
!                           Unstable IF stbltyt > 0.                   !
! stbltyl     'Ledoux'     Ledoux criterion for convection.            !        
!                           Unstable IF stbltyl > 0.                   !
! stbltys     'Schwarz...' Schwarzschild criterion for convection.     !        
!                           Unstable IF stbltys > 0.                   !
! ac(5,j)     'snu/snunr'  Ratio of the relativistic neutrino stress   !
!                           to the nonrelativistic neutrino stress     !
! ac(6,j)     'sp/spnr'    Ratio of the relativistic pressure stress   !
!                           to the nonrelativistic pressure stress     !
! ac(7,j)     'sg/sgnr'    Ratio of the relativistic gravitational     !
!                           stress to the nonrelativistic gravitation- !
!                           al stress                                  !
! wBV(j)      'wBV_yl'     Brunt-Vaisala frequency at high densities   !
!                           (positive if unstable) [s^{-1}]            !
! wBV(j)      'twBV_yl'    Brunt-Vaisala growth time  at high          !
!                           densities (positive if unstable) [s]       !
!----------------------------------------------------------------------!        

8100       CONTINUE

IF ( prnttest ) THEN
nedh(8)          = nedh(8) + 1
IF ( nedh(8) < intedh(8)  .or.  intedh(8) <= 0 ) GO TO 9100
nedh(8)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,801)
WRITE (nprint,803)
WRITE (nprint,805)

csjp1            = zero

!-----------------------------------------------------------------------
!
!                   \\\\\ CONVECTIVE STABILITY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Ledoux stability at high densities
!-----------------------------------------------------------------------

wBV              = zero
twBV             = zero
wBV_s            = zero
twBV_s           = zero
wBV_yl           = zero
twBV_yl          = zero

CALL esrgneldnu( jr_min, jr_max, nx, rho, t, ye, ij_ray, ik_ray, yla )
CALL eqstld( jr_min, jr_max, rho, t, yla, ij_ray, ik_ray )

DO j = jr_min, jr_max
  CALL Brunt_Vaisala_yl( j, ij_ray, ik_ray, nx, r, rho, yla, gstrss, wBV(j), &
&  twBV(j), wBV_s(j), twBV_s(j),  wBV_yl(j), twBV_yl(j) )
END DO ! j = jr_min, jr_max

!-----------------------------------------------------------------------
!  Ledoux stability at high densities
!-----------------------------------------------------------------------

DO jv = jr_min,jr_max,idxedh(8)
  j              = jr_max - jv + jr_min

  IF ( j /= jr_max ) THEN
    csjp1        = c_sound
    csnd2        = gam1(j,ij_ray,ik_ray) * aesv(j,1,ij_ray,ik_ray)/rho(j)
    IF ( csnd2 > zero ) THEN
      c_sound    = DSQRT(csnd2)
      udc        = DABS( u(j)/c_sound )
    ELSE
      c_sound    = zero
      udc        = zero
    END IF ! csnd2 > zero
    vs2          = c_sound * csjp1
    dmass        = half * ( dmrst(j) + dmrst(jp1) )
    drho         = rho(jp1) * ( 1.d0 + csqinv * aesv(jp1,2,ij_ray,ik_ray) )  &
&                - rho(j) * ( 1.d0 + csqinv * aesv(j,2,ij_ray,ik_ray) )
    dp           = aesv(jp1,1,ij_ray,ik_ray) - aesv(j,1,ij_ray,ik_ray)
    stbltyt      = vs2 * ( drho/dmass ) - dp/dmass

    stbltyl      = DLOG( rho(j+1)/rho(j) ) - ( 2.d0/ ( gam1(j,ij_ray,ik_ray) &
&                + gam1(j+1,ij_ray,ik_ray) ) ) &
&                * DLOG( aesv(j+1,1,ij_ray,ik_ray)/aesv(j,1,ij_ray,ik_ray) )

    fg2          = ( gam2(j  ,ij_ray,ik_ray) - 1. )/gam2(j  ,ij_ray,ik_ray)
    fg2p         = ( gam2(j+1,ij_ray,ik_ray) - 1. )/gam2(j+1,ij_ray,ik_ray)
    stbltys      = half * ( fg2 + fg2p ) * DLOG( aesv(j+1,1,ij_ray,ik_ray)   &
&                / aesv(j,1,ij_ray,ik_ray) ) - DLOG( t(j+1)/t(j) )

    stbltysf     = - ( aesvy(j,1,ij_ray,ik_ray)/aesvd(j,1,ij_ray,ik_ray) )   &
&                * DLOG( ye(j+1)/ye(j) ) * ye(j)/rho(j)

  ELSE ! j = jr_max

    stbltyt      = zero
    stbltyl      = zero
    stbltys      = zero
    stbltysf     = zero

  END IF ! j ne jr_max
  sch            = 'n'
  ldx            = 'n'
  sfr            = 'n'
!       CALL cvstblty( j, jr_min, jr_max, sch, ldx, sfr )

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

WRITE (nprint,807) j, dth(j), dtj(j), dtjr(j), dtmpmn(j,1,ij_ray,ik_ray), &
& stbltyt, stbltyl, ldx, stbltys, sch, stbltysf, sfr, ac(5,j), ac(6,j),   &
& ac(7,j), wBV(j), twBV(j), twBV_s(j), twBV_yl(j)

END DO

!**********************************************************************!
!                                                                      !
!                  Hydrodynamic data.                      part 9      !
!                                                                      !
!        Ledoux Convection.                                            !
!----------------------------------------------------------------------!        
!                           1     drho            drho                 !
! aledoux(j)  'aledoux'    --- [( ---- )      - ( ---- )     ].        !
!                          rho     dr   star       dr   blob           !
! psclht(j)   'psclht'     Pressure scale-height at the interface      !
!                           between radial zone j and j+1 (cm).        !
! lmix(j)     'lmix'       Mixing length at the interface between      !
!                           radial zone j and j+1 (cm).                !
! usound(j)   'usound'     Sound velocity at the interface between     !
!                           radial zone j and j+1 (cm).                !
! ulmixl(j)   'ulmixl'     'Mixing length' velocity at the interface   !
!                           between radial zone j and j+1 (cm).        !
! ulcnvct(j)  'ulcnvct'    Convective velocity at the interface        !
!                           between radial zone j and j+1 (cm).        !
! dtmev       'dt'         Change in the temperature of radial zone j  !
!                           due to convection (MeV).                   !
! dyecnvt(j,ij_ray,ik_ray) 'dye'   Change in the electron fraction of radial   !
!                           zone j due to convection.                  !
! pcnvct(j)   'pcnvct'     Turbulent pressure due to convection        !
!                           (dynes/cm3).                               !
! scnvct(j)   'scnvct'     Rate at which energy is transferred to      !
!                           matter due to dissipation of turbulent     !
!                           eddies (ergs/g).                           !
! drhobuoy    'droye/dros' Ratio of (dlnro/dlnye)p,s to                !
!                           (dlnro/dlns)p,ye.                          !
!----------------------------------------------------------------------!        

9100 CONTINUE

IF ( prnttest ) THEN
nedh(9)     = nedh(9) + 1
IF ( nedh(9) < intedh(9)  .or.  intedh(9) <= 0 ) GO TO 9900
nedh(9)     = 0
END IF ! prnttest


WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,901)
WRITE (nprint,903)
WRITE (nprint,905)

DO jv = jr_min,jr_max,idxedh(9)
j                 = jr_max - jv + jr_min

dtmev             = dtmpmn(j,9,ij_ray,ik_ray) * kmev
denom             = aesvd(j,3,ij_ray,ik_ray) * aesvt(j,1,ij_ray,ik_ray)   &
&                 - aesvt(j,3,ij_ray,ik_ray) * aesvd(j,1,ij_ray,ik_ray)
dddsyp            = aesvt(j,1,ij_ray,ik_ray)/( denom + epsilon )
dddysp            = ( aesvy(j,1,ij_ray,ik_ray) * aesvt(j,3,ij_ray,ik_ray) &
&                 - aesvy(j,3,ij_ray,ik_ray) * aesvt(j,1,ij_ray,ik_ray) ) &
&                 / ( denom + epsilon )
dlnddlnsyp        = aesv(j,3,ij_ray,ik_ray) * dddsyp/rho(j)
dlnddlnysp        = ye(j) * dddysp/rho(j)
drhobuoy          = dlnddlnysp/( dlnddlnsyp + epsilon )

!----------------------------------------------------------------------!
!        Print.                                                        !
!----------------------------------------------------------------------!
!       IF ( iprint .ne. 0 )
!     *  WRITE (nprint,907) j,aledoux(j),psclht(j),lmix(j),usound(j),
!     *  ulmixl(j),ulcnvct(j),dtmev,dyecnvt(j,ij_ray,ik_ray),pcnvct(j),
!     *  scnvct(j),drhobuoy
!                                                                      !
!       IF ( iplot .ne. 0 )
!     *  WRITE (nprint,107) j,aledoux(j),psclht(j),lmix(j),usound(j),
!     *  ulmixl(j),ulcnvct(j),dtmev,dyecnvt(j,ij_ray,ik_ray),pcnvct(j),
!     *  scnvct(j),drhobuoy
!                                                                      !

END DO

9900 CONTINUE

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (rnnu, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'rnnu      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (y_nu, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'y_nu      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (yla, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'yla       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (wBV, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'wBV       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (twBV, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'twBV      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (wBV_s, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'wBV_s     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (twBV_s, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'twBV_s    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (wBV_yl, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'wBV_yl    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (twBV_yl, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'twBV_yl   '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ac, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'ac        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dz, STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'dz        '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE edith
