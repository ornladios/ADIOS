SUBROUTINE editu( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editu
!    Module:       editu
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/18/96
!
!    Purpose:
!      To edit energy data.
!
!    Subprograms called:
!  date_and_time_print : prints date and time
!  gammaj_x            : computes the EOS gammas
!  nu_number           : computes neutrino numbers and energies in zones
!  nu_U                : computes the neutrino energy per volume and per unit mass
!  w_cal               : computes the relativistic enthalpies
!  flux                : computes neutrino fluxes
!  gammaj_x            : computes EOS gammas
!  eqstt_x             : interpolates selected EOS quantities
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
!  prnttest   : true  - test to see if printing criteria is satisfied.
!               false - bypass printing criteria test.
!  iprint     : = 0   - do not print to print file.
!               \= 0 - print to print file.
!  nprint     : unit number of print file.
!  iplot      : = 0  - do not print to plot file.
!               \= 0 - print to plot file.
!  nplot      : unit number of plot file.
!  nedu(i)    : editc counter for data set i.
!  intedu(i)  : number of cycles between edits of data set i.
!  idxedu(i)  : edit jr_min, jr_max, and every idxedc(i) radial zone between them for data set i.
!
!    Include files:
!  array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, eos_snc_x_module,
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module,
!  shock_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nx, nez, nnu
USE numerical_module, ONLY : zero, third, half, one, epsilon, csqinv, frpi, &
& frpith
USE physcnst_module, ONLY : cvel, g, ergfoe, ergmev, msolar, rmu, kmev

USE cycle_module
USE e_advct_module, ONLY : dunujvdt, unujv
USE edit_module, ONLY : head, nprint, nlog, nedu, intedu, idxedu, prnttest, &
& pdv, twrk, g_pot
USE eos_snc_x_module, ONLY : aesv, duesrc
USE mdl_cnfg_module, ONLY : grvmss, dmgrv, rstmss, dmrst, u, r, uad, &
& gamgr, wgr, rho, rhor, t, tr, ye, yer
USE nu_dist_module, ONLY : unujinfty, aunu, unujrad, unu, dunu, dunujeadt, &
& dunujnisdt, dunujpadt, dunujtdt, unurad, fluxnu, psi0, psi1, unujea, &
& unujnis, unujpa, unujt
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY :  irelhy 
USE shock_module, ONLY : pq_x
USE t_cntrl_module, ONLY : time

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
LOGICAL                          :: print_sub = .false.
LOGICAL                          :: l_alloc
LOGICAL                          :: lshock

INTEGER                          :: i             ! do index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: jv            ! do index
INTEGER                          :: jminm         ! jr_min - 1
INTEGER                          :: jmaxp         ! jr_max + 1
INTEGER                          :: jmaxm         ! jr_max - 1
INTEGER                          :: k             ! do index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: n_time        ! used for data-and-time

INTEGER                          :: j0            ! zone of outer edge of shock
INTEGER                          :: j0p1          ! j0 + 1
INTEGER                          :: j0p2          ! j0 + 2
INTEGER                          :: j0p3          ! j0 + 3
INTEGER                          :: j1            ! zone of inner edge of shock
INTEGER                          :: j1m1          ! j1 - 1
INTEGER                          :: j1m2          ! j1 - 2
INTEGER                          :: j1m3          ! j1 - 3
INTEGER                          :: js1           ! zone index immediately behind shock middle 
INTEGER                          :: js0           ! js1 + 1, zone index immediately ahead of shock middle 
INTEGER                          :: jd            ! special index
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: csq           ! square of the velocity of light
REAL(KIND=double)                :: uke           ! kinetic energy per unit mass of mass zone j-1/2 (g/cm3)
REAL(KIND=double)                :: ge            ! GR gravitational energy per unit gravitational mass of zone j-1/2 (ergs/g)

REAL(KIND=double)                :: unulsj        ! total neutrino energy radiated out of zone j
REAL(KIND=double)                :: radtot        ! total neutrino energy radiated out of grid
REAL(KIND=double)                :: duesrt        ! total glltch neutrino
REAL(KIND=double)                :: pdvjp1        ! work done by material in grid
REAL(KIND=double)                :: utotgr        ! total GR energy
REAL(KIND=double)                :: utotgres      ! total GR energy minus glitch energy
REAL(KIND=double)                :: utotnr        ! total Newtonian energy
REAL(KIND=double)                :: utotnres      ! total Newtonian energy minus glitch energy

REAL(KIND=double)                :: gemev         ! gravitational binding energy per baryon for matter in zone j-1/2 (MeV)
REAL(KIND=double)                :: ejhmev        ! internal energy per baryon rel fe for matter in zone j-1/2 (MeV)
REAL(KIND=double)                :: ukemev        ! kinetic energy per baryon for matter in zone j-1/2 (MeV)
REAL(KIND=double)                :: pde           ! (3.*aesv(j,1,ij_ray,ik_ray)/rho(j)-(aesv(j,2,ij_ray,ik_ray) rel fe) for matter  in zone j-1/2 (MeV)
REAL(KIND=double)                :: pdez          ! pde * dmrst(j) * ergfoe
REAL(KIND=double)                :: pdemev        ! (3.*aesv(j,1,ij_ray,ik_ray)/rho(j)-(aesv(j,2,ij_ray,ik_ray) rel fe) for matter  in zone j-1/2 (MeV)
REAL(KIND=double)                :: fpirme        ! (4.*pi*aesv(j,1,ij_ray,ik_ray)*r(j)**3-(e sum) = "net ram pressure" (foes)
REAL(KIND=double)                :: ezsum         ! Total GR energy enclosed by mass shell j (foes)
REAL(KIND=double)                :: pdesum        ! (4.*pi*aesv(j,1,ij_ray,ik_ray)*r(j)**3-(e sum) = "net ram pressure" (foes)
REAL(KIND=double)                :: totwrk        ! (4.*pi*aesv(j,1,ij_ray,ik_ray)*r(j)**3-(e sum) = "net ram pressure" (foes)

REAL(KIND=double)                :: ezjmh         ! total energy of zone j-1/2
REAL(KIND=double), PARAMETER     :: befegm = 1.561d+18 ! binding energy of iron (ergs/gm)
REAL(KIND=double)                :: ejmh          ! internal energy - befegm - uad(j)
REAL(KIND=double)                :: shkrmp        ! 4 pi r^{3} p ergmev

REAL(KIND=double)                :: s0            ! entropy ahead of shock
REAL(KIND=double)                :: s1            ! entropy behind shock
REAL(KIND=double)                :: s01           ! (s0 + s1)/2
REAL(KIND=double)                :: msht          ! enclosed mass at shock location (gm)
REAL(KIND=double)                :: msh           ! enclosed mass at shock location (solar masses)
REAL(KIND=double)                :: dush          ! velocity change across shock
REAL(KIND=double)                :: lsh           ! rate of work done by shock on preshocked matter (foes/sec)
REAL(KIND=double)                :: rsp           ! radius at shock
REAL(KIND=double)                :: rhosh         ! density at shock
REAL(KIND=double)                :: rj1           ! radius behind shock
REAL(KIND=double)                :: pj1           ! pressure behind shock
REAL(KIND=double)                :: rho1          ! density behind shock
REAL(KIND=double)                :: u0            ! velocity ahead of shock
REAL(KIND=double)                :: u1            ! velocity behind of shock
REAL(KIND=double)                :: rho00         ! density ahead of shock
REAL(KIND=double)                :: rj0           ! radius ahead of shock
REAL(KIND=double)                :: uesc          ! escape velocity 
REAL(KIND=double)                :: udif          ! uesc - u0
REAL(KIND=double)                :: ramp0         ! 4 * pi * rho0 * rj0^3  udif^2 just ahead of shock
REAL(KIND=double)                :: rampes        ! 4 * pi * rho0 * rj0^3  udif^2 at shock
REAL(KIND=double)                :: rho1d0        ! rho1/rh00
REAL(KIND=double)                :: pl            ! log(pressure) (current)
REAL(KIND=double)                :: ppl           ! log(pressure) (previous)
REAL(KIND=double)                :: rhol          ! log(density) (current)
REAL(KIND=double)                :: rhorl         ! log(density) (previous)
REAL(KIND=double)                :: rmassjd       ! mass enclosed just ahead of shock
REAL(KIND=double)                :: rmassjdm1     ! mass enclosed just behind shock

REAL(KIND=double)                :: gamma1        ! EOS gamma1
REAL(KIND=double)                :: gamma2        ! EOS gamma2
REAL(KIND=double)                :: gamma3        ! EOS gamma3
REAL(KIND=double)                :: gammae        ! effective gamma
REAL(KIND=double)                :: pp            ! pressure
REAL(KIND=double)                :: dpddp         ! d(pressure)/d(density)
REAL(KIND=double)                :: dpdtp         ! d(pressure)/d(temperature)
REAL(KIND=double)                :: dpdyp         ! d(pressure)/d(ye)

REAL(KIND=double)                :: unetaf        ! Net energy given to n-neutrinos by abs and emisin zone j (foes)
REAL(KIND=double)                :: unetsf        ! Net energy given to n-neutrinos by non-iso scattering in zone j (foes)
REAL(KIND=double)                :: unetpf        ! Net energy given to n-neutrinos by pair production in zone j (foes)
REAL(KIND=double)                :: unetvf        ! Net energy given to n-neutrinos by energy advection in zone j (foes)
REAL(KIND=double)                :: unettf        ! Net energy given to n-neutrinos by all processes in zone j (foes)
REAL(KIND=double)                :: unettr        ! Net energy given to n-neutrinos by transport in zone j (foes)
REAL(KIND=double)                :: unersd        ! total energy of n-neutrinos in zone j (foes)

REAL(KIND=double)                :: uneta         ! Energy transfer rate to n-neutrinos by absorption and emission in zone j
REAL(KIND=double)                :: unets         ! Energy transfer rate to n-neutrinos by non-iso scattering in zone j
REAL(KIND=double)                :: unetp         ! Energy transfer rate to n-neutrinos by pair production in zone j
REAL(KIND=double)                :: unett         ! Energy transfer rate to n-neutrinos by by all interactions in zone j
REAL(KIND=double)                :: unetv         ! Energy transfer rate to n-neutrinos by energy advection in zone j
REAL(KIND=double)                :: unetr         ! Energy transfer rate to n-neutrinos by transport in zone j

REAL(KIND=double)                :: fluxc         ! coefficient for computing the neutrino luminosity
REAL(KIND=double)                :: r_lume         ! neutrino luminosity
REAL(KIND=double)                :: denn          ! factor for computing mean energies
REAL(KIND=double)                :: dene          ! factor for computing mean energies
REAL(KIND=double)                :: dene2         ! factor for computing rms energies
REAL(KIND=double)                :: dfnn          ! factor for computing flux weighted rms energies
REAL(KIND=double)                :: dfne2         ! factor for computing flux weighted rms energies
REAL(KIND=double)                :: rnavee        ! neutrino mean energy
REAL(KIND=double)                :: rnrmse        ! neutrino rms energy
REAL(KIND=double)                :: rfrmse        ! neutrino rms energy (flux weighted)
REAL(KIND=double)                :: coefr         ! 1/dmgrv
REAL(KIND=double)                :: fluxcj        ! 4 * pi * r(j)^2
REAL(KIND=double)                :: fluxcjm1      ! 4 * pi * r(j-1)^2
REAL(KIND=double)                :: coefn         ! unu(j)^2 dunu(j)
REAL(KIND=double)                :: coefe         ! unu(j)^3 dunu(j)
REAL(KIND=double)                :: coefe2        ! unu(j)^4 dunu(j)
REAL(KIND=double)                :: coefnm        ! unu(j-1)^2 dunu(j-1)
REAL(KIND=double)                :: coefem        ! unu(j-1)^3 dunu(j-1)
REAL(KIND=double)                :: coefe2m       ! unu(j-1)^4 dunu(j-1)

REAL(KIND=double)                :: dennjd        ! factor for computing mean energies
REAL(KIND=double)                :: denejd        ! factor for computing mean energies
REAL(KIND=double)                :: dene2jd       ! factor for computing mean energies
REAL(KIND=double)                :: dennjdm       ! factor for computing mean energies
REAL(KIND=double)                :: denejdm       ! factor for computing mean energies
REAL(KIND=double)                :: dene2jdm      ! factor for computing mean energies
REAL(KIND=double)                :: enuvjd        ! factor for computing mean energies
REAL(KIND=double)                :: enuvjdm       ! factor for computing mean energies
REAL(KIND=double)                :: enuv2jdm      ! factor for computing mean energies
REAL(KIND=double)                :: enuv2         ! factor for computing mean energies
REAL(KIND=double)                :: enuv2jd       ! factor for computing mean energies
REAL(KIND=double)                :: fluxcjd       ! 4 * pi * r^2 just ahead of shock
REAL(KIND=double)                :: fluxcjdm      ! 4 * pi * r^2 just behind shock
REAL(KIND=double)                :: r_lumjd        ! luminosity just ahead of shock
REAL(KIND=double)                :: r_lumjdm       ! luminosity just behind shock

REAL(KIND=double), PARAMETER     :: pqcrit = 0.5d+00
REAL(KIND=double), PARAMETER, DIMENSION(7) :: rhonud = (/1.d+12, 1.d+11, 1.d+10, 1.d+9, 1.d+8, 1.d+7, 1.d+6/)
REAL(KIND=double), PARAMETER, DIMENSION(7) :: rnud   = (/1.d+6, 3.d+6, 1.d+7, 3.d+7, 1.d+8, 3.d+8, 1.d+9/)
REAL(KIND=double), DIMENSION(14) :: rhonu         ! radius at shock
REAL(KIND=double), DIMENSION(14) :: rnu           ! radius at shock
REAL(KIND=double), DIMENSION(14) :: rmass         ! enclosed mass at shock
REAL(KIND=double), DIMENSION(14) :: tnu           ! temperature at shock (MeV)
REAL(KIND=double), DIMENSION(14) :: snu           ! entropy at shock (MeV)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: u_gr      ! zone GR energiesezjmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: u_nt      ! zone Newtonian energies
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: u_grt     ! GR energies summed to j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: u_ntt     ! Newtonian energies summed to j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: a1        ! working array

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: lum       ! luminosity
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: enuv      ! mean energy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: enuvrms   ! rms energy

REAL(KIND=double)                :: tnujd         ! temperature just ahead of shock (MeV)
REAL(KIND=double)                :: tnujdm1       ! temperature just behind shock (MeV)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1x)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (1x,'cycle =',i6,';',1x,'date:',i2,'/',i2,'/',i2)
  101 FORMAT (55x,'Energy data')
  103 FORMAT (54x,13('-')/)
  105 FORMAT ('   j  ge/zone    ie/zone    ke/zone    ne/zone   iwrk/zone  twrk/zone &
& ke+ge+twk  ie+ne+iwk  ge+ie+ke  ge+ie+ke+ne tote+nulss'/)
  107 FORMAT (1x,i4,11(es11.4))
  201 FORMAT (44x,'Energy data (summed from jr_min to j)')
  203 FORMAT (43x,37('-')/)
  205 FORMAT ('   j     ge         ie         ke         ne      i work     t work &
&   ke+ge+twk  ie+ne+iwk  ge+ie+ke  ge+ie+ke+ne tote+nulss'/)
  209 FORMAT ('0totne rad (energy radiated as neutrinos) =  ',es11.4,' foes')
  211 FORMAT (' pdv (work done by core on surroundings) =   ',es11.4,' foes')
  213 FORMAT (' totge + totie + totke + totne =             ',es11.4,' foes')
  215 FORMAT ('0totge + totie + totke + totne + rad + pdv = ',es11.4,' foes')
  217 FORMAT ('0totge + totie + totke + totne + rad + pdv - duesrc = ',es11.4,' foes')
  219 FORMAT ('0totge + totie + totke + totne + rad + pdv = ',es11.4,' foes, computed nonrelativisticly')
  221 FORMAT ('0totge + totie + totke + totne + rad + pdv - duesrc = ',es11.4,' foes, computed nonrelativisticly')
  301 FORMAT (44x,'Energy data (summed from jr_max to j)')
  303 FORMAT (43x,37('-')/)
  305 FORMAT ('   j    ge         ie         ke         ne       i work     t work    ke+ge+twk &
& ie+ne+iwk  ge+ie+ke  ge+ie+ke+ne nu e loss'/)
  401 FORMAT (43x,'Energy and shock ram pressure data')
  403 FORMAT (42x,36('-')/)
  405 FORMAT ('   j    ge/b       ie/b       ke/b      3p/d-e        p        e sum   3p/d-e sum  &
& 4pipr3-e   nu loss     pdv'/)
  407 FORMAT (1x,i4,10(es11.3))
  501 FORMAT (55x,'Shock data')
  503 FORMAT (54x,12('-')///)
  505 FORMAT (' No strong outward propagating shock wave is present')
  507 FORMAT (' j1 cannot be found')
  509 FORMAT (' Shock center cannot be found')
  511 FORMAT (' Shock is located at r=',es10.3,' cm from center')
  513 FORMAT (' The shock is at a density of ',es10.3,' g/cm3')
  515 FORMAT (' The shock encloses a rest mass of',es10.3,' msolar')
  517 FORMAT (' The shock is encountering a mass luminosity of',es10.3,' 1.e+50 ergs/sec')
  519 FORMAT (' The shock ram pressure is',es10.3,' foe')
  521 FORMAT (' The shock ram pressure necessary for u1=0. is',es10.3,' foe; u0=',es10.3, ' u1=',es10.3)
  523 FORMAT (' The shock ram pressure necessary for u1=uesc is',es10.3,' foe; uesc=',es10.3)
  525 FORMAT (' rho0=',es10.3,' rho1=',es10.3,' rho1/rho0=',es10.3)
  527 FORMAT (' Ahead of shock,            gamma1=',es10.3,' gamma2=',es10.3,' gamma3=',es10.3,' gammae=',es10.3)
  529 FORMAT (' Ahead of shock middle,     gamma1=',es10.3,' gamma2=',es10.3,' gamma3=',es10.3,' gammae=',es10.3)
  531 FORMAT (' Behind shock middle,       gamma1=',es10.3,' gamma2=',es10.3,' gamma3=',es10.3,' gammae=',es10.3)
  533 FORMAT (' Behind shock,              gamma1=',es10.3,' gamma2=',es10.3,' gamma3=',es10.3,' gammae=',es10.3)
  601 FORMAT (40x,'Summary of electron neutrino energy data')
  603 FORMAT (39x,42('-')/)
  605 FORMAT ('   j  enu a e    enu nis    enu p-a     enu v     enu tot   enu trsp   enu resd &
&   enu lum   e nu av n  e nu rms n e nu flx n'/)
  607 FORMAT (1x,i4,11(es11.4))
  611 FORMAT (38x,'Summary of electron antineutrino energy data')
  613 FORMAT (37x,46('-')/)
  615 FORMAT ('   j  anu a e    anu nis    anu p-a     anu v     anu tot   anu trsp   anu resd &
&   anu lum   a nu av n  a nu rms n a nu flx n'/)
  617 FORMAT (1x,i4,11(es11.4))
  621 FORMAT (42x,'Summary of x-neutrino energy data')
  623 FORMAT (41x,35('-')/)
  625 FORMAT ('   j  tnu a e    tnu nis    tnu p-a     tnu v     tnu tot   tnu trsp   tnu resd &
&   tnu lum   t nu av n  t nu rms n t nu flx n'/)
  627 FORMAT (1x,i4,11(es11.4))
  631 FORMAT (42x,'Summary of x-antineutrino energy data')
  633 FORMAT (41x,35('-')/)
  635 FORMAT ('   j tbnu a e   tbnu nis   tbnu p-a    tbnu v    tbnu tot  tbnu trsp  tbnu resd &
&  tbnu lum  tb nu av n tb nu rms n tb nu flx n'/)
  637 FORMAT (1x,i4,11(es11.4))
  701 FORMAT (36x,'Matter-electron neutrino energy transfer rates')
  703 FORMAT (35x,48('-')/)
  705 FORMAT ('   j  enur a e   enur nes   enur p-a   enur tot    enur v   enur trsp'/)
  707 FORMAT (1x,i4,6(es11.4))
  711 FORMAT (34x,'Matter-electron antineutrino energy transfer rates')
  713 FORMAT (33x,52('-')/)
  715 FORMAT ('   j  anur a e   anur nes   anur p-a   anur tot    anur v   anur trsp'/)
  717 FORMAT (1x,i4,6(es11.4))
  721 FORMAT (40x,'Matter-x-neutrino energy transfer rates')
  723 FORMAT (39x,41('-')/)
  725 FORMAT ('   j  tnur a e   tnur nes   tnur p-a   tnur tot    tnur v   tnur trsp'/)
  727 FORMAT (1x,i4,6(es11.4))
  731 FORMAT (40x,'Matter-x-antineutrino energy transfer rates')
  733 FORMAT (39x,41('-')/)
  735 FORMAT ('   j tbnur a e  tbnur nes  tbnur p-a  tbnur tot   tbnur v  tbnur trsp'/)
  737 FORMAT (1x,i4,6(es11.4))
  801 FORMAT (35x,'Summary of neutrino luminosity and mean energy data')
  803 FORMAT (34x,53('-')/)
  805 FORMAT (9x,'Ellapsed time=',es14.7/)
  807 FORMAT ('   i      r         rho      enu lum    anu lum    tnu lum   tbnu lum'/)
  809 FORMAT (1x,i4,6(es11.3))
  811 FORMAT (///)
  821 FORMAT ('   i      r         rho     enu eav1   anu eav1   tnu eav1  tbnu eav1'/)
  823 FORMAT (1x,i4,6(es11.3))
  831 FORMAT ('   i      r         rho     enu eav2   anu eav2   tnu eav2  tbnu eav2'/)
  833 FORMAT (1x,i4,6(es11.3))
  851 FORMAT ('   i      r         rho     M (solar)   T (MeV)       s'/)
  853 FORMAT (1x,i4,5(es11.3))
 1001 FORMAT (' Allocation problem for array ',a10,' in editu')
 2001 FORMAT (' Deallocation problem for array ',a10,' in editu')
 8001 FORMAT (' jd cannot be found in subroutine editu for rhonu(',i1,')=',es10.3)
 8101 FORMAT (' jd cannot be found in subroutine editu for rnu(',i1,')=',es10.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  first             = .false.
  csq               = cvel**2            ! !**2
END IF ! first
n_time              = nprint

!-----------------------------------------------------------------------
!  Allocate arrays if any of the subedits are executed
!-----------------------------------------------------------------------

print_sub           = .false.
l_alloc             = .false.

IF ( prnttest ) THEN
  DO i = 1,12
    IF ( nedu(i) + 1 >= intedu(i) ) print_sub = .true.
  END DO
ELSE
  print_sub         = .true.
END IF

IF ( print_sub ) THEN

  ALLOCATE (u_gr(nx,20), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_gr      '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (u_nt(nx,20), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_nt      '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (u_grt(nx,20), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_grt     '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (u_ntt(nx,20), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_ntt     '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (a1(nx,20), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'a1        '; WRITE (nlog,1001) var_name; END IF

  ALLOCATE (lum(14,nnu), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'lum       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (enuv(14,nnu), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'enuv      '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (enuvrms(14,nnu), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'enuvrms   '; WRITE (nlog,1001) var_name; END IF

  l_alloc           = .true.

!-----------------------------------------------------------------------
!  Initialize arrays
!-----------------------------------------------------------------------

  u_gr              = zero
  u_nt              = zero
  u_grt             = zero
  u_ntt             = zero
  a1                = zero

  lum               = zero
  enuv              = zero
  enuvrms           = zero

  CALL w_cal( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, wgr, nx )
  DO n = 1,nnu
    CALL nu_number( jr_min, jr_max, n, ij_ray, ik_ray, nx, nez, nnu, &
&    r, u, psi0, psi1 )
  END DO
  CALL nu_U( jr_min, jr_max, rho, nx, nnu )
  DO n = 1,nnu
    CALL flux( jr_min, jr_max, n )
  END DO

END IF

!**********************************************************************!
!                                                                      !
!                  Energy data                             part 1      !
!                                                                      !
!      u_gr: general relativistiv energies                             !
!      u_nt: Newtonian energies                                        !
!----------------------------------------------------------------------!
! u_(j,1):    'ge/zone'    Gravitational energy of zone j (*foes,      !
!                           i.e., 1.e+51 ergs)                         !
! u_(j,2):    'ie/zone'    Internal energy of zone j (foes)            !
! u_(j,3):    'ke/zone'    Kinetic energy of zone j (foes)             !
! u_(j,4):    'ne/zone'    Net energy of all types of neutrinos in     !
!                           zone j (foes)                              !
! u_(j,5):    'iwrk/zone'  Net work done by zone j at the expense of   !
!                           its internal energy (foes)                 !
! u_(j,6):    'twrk/zone'  Net work done by zone j at the expense of   !
!                           its gravitational and kinetic energy       !
!                           (foes)                                     !
! u_(j,7):    'ke+ge+twk'  Sum of ke/zone, ge/zone, and twrk/zone      !
!                           (foes)                                     !
!                           (should be constant during a newtonian     !
!                           calculation)                               !
! u_(j,8):    'ie+ne+iwk'  Sum of ie/zone, ne/zone, and iwrk/zone      !
!                           (foes)                                     !
!                           (should be constant during a newtonian     !
!                           calculaiton)                               !
! u_(j,9):    'ge+ie+ke'   Sum of ge/zone, ie/zone, and ke/zone        !
!                           (foes)                                     !
! u_(j,10):  'ge+ie+ke+ne' Sum of ge/zone, ie/zone, ke/zone, and       !
!                           ne/zone (foes))                            !
! totee:      'tote+nulss' Sum of ge/zone, ie/zone, ke/zone, ne/zone   !
!                           and energy lost by all neutrino types      !
!                           (foes)                                     !
!----------------------------------------------------------------------!

!-----------------------------------------------------------------------
!  The following must be executed if any of the first four subedits
!   are executed
!-----------------------------------------------------------------------

print_sub           = .false.

IF ( prnttest ) THEN
  DO i = 1,4
    IF ( nedu(i) + 1 >= intedu(i) ) print_sub = .true.
  END DO
ELSE
  print_sub         = .true.
END IF

IF ( print_sub ) THEN

  DO j = 2,jr_max

!-----------------------------------------------------------------------
!  u_ny(j,1):  Newtonian gravitational energy of mass zone j-1/2 (foes)
!-----------------------------------------------------------------------

    u_nt(j,1)       = half * ( - rstmss(j-1) * g/( r(j-1) + epsilon )      &
&                              - rstmss(j)   * g/( r(j  ) + epsilon ) )    &
&                   * dmrst(j) * ergfoe

!-----------------------------------------------------------------------
!  u_gr(j,1):  GR gravitational energy of mass zone j-1/2 (foes)
!-----------------------------------------------------------------------

    IF ( irelhy == 0 ) THEN
      u_gr(j,1)     = u_nt(j,1)
    ELSE IF ( irelhy == 1 ) THEN
      u_gr(J,1)     = g_pot(j)
    ELSE
      u_gr(j,1)     = half * ( - grvmss(j-1) * g/( r(j-1) + epsilon )      &
&                   - grvmss(j  ) * g/( r(j  ) + epsilon ) )               &
&                   * dmrst(j) * ergfoe * ( one + ( aesv(j,2,ij_ray,ik_ray) + aunu(j) ) * csqinv )
    END IF

!-----------------------------------------------------------------------
!  u_gr(j,2):  GR internal energy of mass zone j-1/2 (foes)
!-----------------------------------------------------------------------

    u_gr(j,2)       = aesv(j,2,ij_ray,ik_ray) * dmrst(j) * ergfoe

!-----------------------------------------------------------------------
!  u_nt(j,2):  Newtonian internal energy of mass zone j-1/2 (foes)
!-----------------------------------------------------------------------

    u_nt(j,2)       = u_gr(j,2)

!-----------------------------------------------------------------------
!  Kinetic energy
!-----------------------------------------------------------------------

    uke             = half * u(j) * u(j)

!-----------------------------------------------------------------------
!  u_gr(j,3):  GR kinetic energy of mass zone j-1/2 (foes)
!-----------------------------------------------------------------------

    u_gr(j,3)       = uke * dmrst(j) * ergfoe                              &
&                   * ( one + ( aesv(j,2,ij_ray,ik_ray) + aunu(j) ) * csqinv )

!-----------------------------------------------------------------------
!  u_nt(j,3):  Newtonian kinetic energy of mass zone j-1/2 (foes)
!-----------------------------------------------------------------------

    u_nt(j,3)       = uke * dmrst(j) * ergfoe

!-----------------------------------------------------------------------
!  u_gr(j,4):  GR energy of all types of neutrinos in zone j-1/2 (foes)
!-----------------------------------------------------------------------

    u_gr(j,4)       = SUM ( unujinfty(j,:) ) * ergfoe

!-----------------------------------------------------------------------
!  u_nt(j,4):  Newtonian energy of all types of neutrinos in zone
!   j-1/2 (foes)
!-----------------------------------------------------------------------

    u_nt(j,4)       = u_gr(j,4)

!-----------------------------------------------------------------------
!  u_nt(j,5):  Net work done by zone j at the expense of its internal
!   energy (foes)
!-----------------------------------------------------------------------

    u_gr(j,5)       = pdv(j) * ergfoe
    u_nt(j,5)       = u_gr(j,5)

!-----------------------------------------------------------------------
!  u_gr(j,6):  Net work done by zone j at the expense of its
!   gravitational and kinetic energy energy (foes)
!-----------------------------------------------------------------------

    u_gr(j,6)       = twrk(j) * ergfoe
    u_nt(j,6)       = u_gr(j,6)

!-----------------------------------------------------------------------
!  u_gr(j,7):  Sum of ke/zone, ge/zone, and twrk/zone (foes)
!   (Should be constant during a Newtonian calculation)
!-----------------------------------------------------------------------

    u_gr(j,7)       = u_gr(j,3) + u_gr(j,1) + u_gr(j,6)
    u_nt(j,7)       = u_nt(j,3) + u_nt(j,1) + u_nt(j,6)

!-----------------------------------------------------------------------
!  u_gr(j,8):  Sum of ie/zone, ne/zone, and iwrk/zone (foes)
!   (Should be constant during a Newtonian calculation)
!-----------------------------------------------------------------------

    u_gr(j,8)       = u_gr(j,2) + u_gr(j,4) + u_gr(j,5)
    u_nt(j,8)       = u_nt(j,2) + u_nt(j,4) + u_nt(j,5)

!-----------------------------------------------------------------------
!  u_gr(j,9):  Total GR energy of mass zone j-1/2 excluding
!   neutrinos (foes)
!-----------------------------------------------------------------------

    ge              = u_gr(j,1)/ergfoe/( dmrst(j)                          &
&                   * ( one + ( aesv(j,2,ij_ray,ik_ray) + aunu(j) ) * csqinv ) )

    u_gr(j,9)       = ( csq * dmrst(j) * ( one + aesv(j,2,ij_ray,ik_ray) * csqinv ) &
&                   * DSQRT( dabs( one + 2.d+00 * csqinv * ( uke + ge ) ) )         &
&                   - csq * dmrst(j) ) * ergfoe
    u_nt(j,9)       = u_nt(j,1) + u_nt(j,2) + u_nt(j,3)

!-----------------------------------------------------------------------
!  u_gr(j,10): Total GR energy of mass zone j-1/2 including neutrinos
!   (foes)
!-----------------------------------------------------------------------

       u_gr(j,10)   = ( csq * dmgrv(j) - csq * dmrst(j) ) * ergfoe
       u_nt(j,10)   = u_nt(j,9) + u_nt(j,4)

!-----------------------------------------------------------------------
!  u_gr(j,11): Total GR energy of mass zone j-1/2 including neutrinos
!   plus the energy of neutrinos that left the zone minus the energy
!   of neutrinos that entered the zone (foes)
!-----------------------------------------------------------------------

    u_gr(j,11)      = u_gr(j,10) + ( unujrad(j  ,1,ij_ray,ik_ray)                   &
&                   + unujrad(j  ,2,ij_ray,ik_ray) + unujrad(j  ,3,ij_ray,ik_ray)   &
&                   -   unujrad(j-1,1,ij_ray,ik_ray) + unujrad(j-1,2,ij_ray,ik_ray) &
&                   + unujrad(j-1,3,ij_ray,ik_ray) ) * ergfoe
    u_nt(j,11)      = u_nt(j,10) + ( unujrad(j  ,1,ij_ray,ik_ray)                   &
&                   + unujrad(j  ,2,ij_ray,ik_ray) + unujrad(j  ,3,ij_ray,ik_ray)   &
&                   - unujrad(j-1,1,ij_ray,ik_ray) + unujrad(j-1,2,ij_ray,ik_ray)   &
&                   + unujrad(j-1,3,ij_ray,ik_ray) ) * ergfoe

  END DO

END IF ! print_sub

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(1)           = nedu(1) + 1
  IF ( nedu(1) >= intedu(1) ) THEN
    print_sub       = .true.
    nedu(1)         = 0
  END IF ! nedu(1) >= intedu(1)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,101)
  WRITE (nprint,103)
  WRITE (nprint,105)

  DO jv = jr_min,jr_max,idxedu(1)
    j               = jr_max - jv + jr_min
    IF ( irelhy == 0 ) THEN
      WRITE (nprint,107) j,(u_nt(j,k),k=1,11)
    ELSE IF ( irelhy == 1 ) THEN
      WRITE (nprint,107) j,(u_gr(j,k),k=1,11)
    ELSE
      WRITE (nprint,107) j,(u_gr(j,k),k=1,11)
    END IF ! irelhy = 0
  END DO

END IF ! print_sub

!-----------------------------------------------------------------------
!
!           ||||| Energy data (summed from jr_min to j) |||||
!
!-----------------------------------------------------------------------

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(2)           = nedu(2) + 1
  IF ( nedu(2) >= intedu(2) ) THEN
    print_sub       = .true.
    nedu(2)         = 0
  END IF ! nedu(2) >= intedu(2)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,201)
  WRITE (nprint,203)
  WRITE (nprint,205)

  jminm             = jr_min - 1

!-----------------------------------------------------------------------
!  Initialize u_grt and u_ntt
!-----------------------------------------------------------------------

  DO j = jminm,jr_max
    DO k = 1,10
      u_grt(j,k)    = zero
      u_ntt(j,k)    = zero
    END DO
  END DO

!-----------------------------------------------------------------------
!  Sum zone energies
!-----------------------------------------------------------------------

  DO j = jr_min,jr_max

    DO k = 1,10
      u_grt(j,k)    = u_grt(j-1,k) + u_gr(j,k)
      u_ntt(j,k)    = u_ntt(j-1,k) + u_nt(j,k)
    END DO

    unulsj          = zero
    DO n = 1,nnu
      unulsj        = unulsj + unujrad(j,n,ij_ray,ik_ray)
    END DO
    unulsj          = unulsj * ergfoe
    u_grt(j,11)     = u_grt(j,10) + unulsj
    u_ntt(j,11)     = u_ntt(j,10) + unulsj

  END DO

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

  DO jv = jr_min,jr_max,idxedu(2)
    j               = jr_max - jv + jr_min

    IF ( irelhy == 0 ) THEN
      WRITE (nprint,107) j,(u_ntt(j,k),k=1,11)
    ELSE IF ( irelhy == 1 ) THEN
      WRITE (nprint,107) j,(u_grt(j,k),k=1,11)
    ELSE
      WRITE (nprint,107) j,(u_grt(j,k),k=1,11)
    END IF ! irelhy = 0

  END DO

!-----------------------------------------------------------------------
!  Total energies
!-----------------------------------------------------------------------

  radtot            = zero
  DO n = 1,nnu
    radtot          = radtot + unurad(n,ij_ray,ik_ray)
  END DO
  radtot            = radtot * ergfoe

  duesrt            = zero
  DO j = 2,jr_max
    duesrt          = duesrt + duesrc(j,ij_ray,ik_ray) * dmrst(j)
  END DO
  duesrt            = duesrt * ergfoe

  pdvjp1            = pdv(jr_max+1) * ergfoe
  utotgr            = u_grt(jr_max,10) + radtot + pdvjp1
  utotgres          = utotgr - duesrt
  utotnr            = u_ntt(jr_max,10) + radtot + pdvjp1
  utotnres          = utotnr - duesrt

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

  WRITE (nprint,209) radtot
  WRITE (nprint,211) pdvjp1
  WRITE (nprint,213) u_grt(jr_max,10)
  WRITE (nprint,215) utotgr
  WRITE (nprint,217) utotgres
  WRITE (nprint,219) utotnr
  WRITE (nprint,221) utotnres

END IF ! print_sub

!-----------------------------------------------------------------------
!
!           ||||| Energy data (summed from jr_max to j) |||||
!
!-----------------------------------------------------------------------

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(3)           = nedu(3) + 1
  IF ( nedu(3) >= intedu(3) ) THEN
    print_sub       = .true.
    nedu(3)         = 0
  END IF ! nedu(3) >= intedu(3)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,301)
  WRITE (nprint,303)
  WRITE (nprint,305)

  jmaxp             = jr_max + 1

!-----------------------------------------------------------------------
!  Initialize u_grt and u_ntt
!-----------------------------------------------------------------------

  DO j = jr_min,jmaxp
    DO k = 1,10
      u_grt(j,k)    = zero
      u_ntt(j,k)    = zero
    END DO
  END DO

!-----------------------------------------------------------------------
!  Sum energies from jr_max
!-----------------------------------------------------------------------

  DO jv = jr_min,jr_max
  j           = jr_max - jv + jr_min
    DO k = 1,10
      u_grt(j,k)    = u_grt(j+1,k) + u_gr(j,k)
      u_ntt(j,k)    = u_ntt(j+1,k) + u_nt(j,k)
    END DO
  END DO

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

  DO jv = jr_min,jr_max,idxedu(3)

    j               = jr_max - jv + jr_min
    unulsj          = ( unujrad(j,1,ij_ray,ik_ray) + unujrad(j,2,ij_ray,ik_ray) &
&                   +   unujrad(j,3,ij_ray,ik_ray) ) * ergfoe

    IF ( irelhy == 0 ) THEN
      WRITE (nprint,107) j,(u_ntt(j,k),k=1,10),unulsj
    ELSE IF ( irelhy == 1 ) THEN
      WRITE (nprint,107) j,(u_grt(j,k),k=1,10),unulsj
    ELSE
      WRITE (nprint,107) j,(u_grt(j,k),k=1,10),unulsj
    END IF ! irelhy = 0

  END DO
  
END IF ! print_sub

!**********************************************************************!
!                                                                      !
!                  Energy and shock ram pressure data                  !
!                                                                      !
!----------------------------------------------------------------------!
! gemev       'ge/b'       Gravitational binding energy per baryon     !
!                           for matter in zone j-1/2 (MeV)             !
! ejhmev      'ie/b'       Internal energy per baryon rel fe for       !
!                           matter in zone j-1/2 (MeV)                 !
! ukemev      'ke/b'       Kinetic energy per baryon for matter in     !
!                           zone j-1/2 (MeV)                           !
! pdemev      '3p/d-e'     (3.*aesv(j,1,ij_ray,ik_ray)/rho(j)          !
!                           - (aesv(j,2,ij_ray,ik_ray) rel fe)         !
!                           for matter  in zone j-1/2 (MeV)            !
! aesv(j,1,ij_ray,ik_ray)   'p'                                        !
!                          Pressure of matter in zone j-1/2            !
!                           (dynes/ cm2)                               !
! a1(j,1)     'e sum'      Total gravitational energy (rel fe) - rest  !
!                           mass energy interior to zone j             !
!                           (foes)                                     !
! a1(j,2)     '3p/d-e sum' (3.*aesv(j,1,ij_ray,ik_ray)/rho(j)          !
!                           - (aesv(j,2,ij_ray,ik_ray) rel fe)         !
!                           summed from jr_min to j (foes)             !
! fpirme      '4pipr3-e'   (4.*pi*aesv(j,1,ij_ray,ik_ray)*r(j)**3      !
!                           -(e sum) = "net ram pressure" (foes)       !
! unulsj      'nu loss'    Total neutrino energy loss from zone        !
!                           boundary j (foes)                          !
! a1(j,3)     'pdv'        Total work done on material inside zone j   !
!                           (foes)                                     !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(4)           = nedu(4) + 1
  IF ( nedu(4) >= intedu(4) ) THEN
    print_sub       = .true.
    nedu(4)         = 0
  END IF ! nedu(4) >= intedu(4)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,401)
  WRITE (nprint,403)
  WRITE (nprint,405)

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

  ezsum             = zero
  pdesum            = zero
  totwrk            = zero

  DO j = 2,jr_max

!-----------------------------------------------------------------------
!  u_gr(j,1): GR gravitational energy per unit mass of zone j-1/2 (ergs)
!-----------------------------------------------------------------------

    u_gr(j,1)       = half * ( - grvmss(j-1) * g/( r(j-1) + epsilon ) - grvmss(j  ) * g/( r(j  ) + epsilon ) )

!-----------------------------------------------------------------------
!  ezjmh: Total GR energy of mass zone j-1/2 including neutrinos and
!   nuclear binding energy (foes)
!-----------------------------------------------------------------------

    ezjmh           = u_gr(j,10) - befegm * dmrst(j) * ergfoe

!-----------------------------------------------------------------------
!  a1(j,1): Total GR energy including neutrinos and nuclear binding
!   energy enclosed by mass shell j (foes)
!-----------------------------------------------------------------------

    ezsum           = ezsum + ezjmh
    a1(j,1)         = ezsum

!-----------------------------------------------------------------------
!  a1(j,2): 3*p/rho - e enclosed by mass shell j (foes)
!-----------------------------------------------------------------------

    ejmh            = aesv(j,2,ij_ray,ik_ray) - befegm - uad(j)
    pde             = 3.d+00 * aesv(j,1,ij_ray,ik_ray)/rho(j) - ejmh
    pdez            = pde * dmrst(j) * ergfoe
    pdesum          = pdesum + pdez
    a1(j,2)         = pdesum

!-----------------------------------------------------------------------
!  a1(j,3): Total work done on material inside zone j (foes)
!-----------------------------------------------------------------------

    totwrk          = totwrk - pdv(j) - twrk(j)
    a1(j,3)         = totwrk * ergfoe

!-----------------------------------------------------------------------
!  a1(j,4): Escape velocity at mass shell j (cm/s)
!-----------------------------------------------------------------------

    a1(j,4)         = DSQRT( 2.d+00 * g * grvmss(j)/r(j) )

  END DO

  DO jv = jr_min,jr_max,idxedu(4)
    j               = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  uke: Kinetic energy per unit gravitational mass of mass zone j-1/2
!   (ergs/g)
!-----------------------------------------------------------------------

    uke             = half * ( half * u(j) * u(j) + half * u(j-1) * u(j-1) ) &
&                   * ( one + ( aesv(j,2,ij_ray,ik_ray) + aunu(j) ) * csqinv )

!-----------------------------------------------------------------------
!  ukemev: Kinetic energy per baryon in mass zone j-1/2 (MeV)
!-----------------------------------------------------------------------

    ukemev          = uke * rmu/ergmev

!-----------------------------------------------------------------------
!  ejhmev: Internal energy per baryon rel fe for matter in zone j-1/2
!   (MeV)
!-----------------------------------------------------------------------

    ejmh            = aesv(j,2,ij_ray,ik_ray) - befegm - uad(j)
    ejhmev          = ejmh * rmu/ergmev

!-----------------------------------------------------------------------
!  gemev: Gravitational binding energy per baryon for matter in zone
!   j-1/2 (MeV)
!-----------------------------------------------------------------------

    gemev           = u_gr(j,1) * rmu/ergmev

!-----------------------------------------------------------------------
!  pdemev: (3.*aesv(j,1,ij_ray,ik_ray)/rho(j)-(aesv(j,2,ij_ray,ik_ray)
!   rel fe) for matter in zone j-1/2 MeV)
!-----------------------------------------------------------------------

    pde             = 3.d0 * aesv(j,1,ij_ray,ik_ray)/rho(j) - ejmh
    pdemev          = pde * rmu/ergmev

!-----------------------------------------------------------------------
!  fpirme: (4.*pi*aesv(j,1,ij_ray,ik_ray)*r(j)**3-(e sum) = "net ram pressure"
!   (foes)
!-----------------------------------------------------------------------

    fpirme          = ( frpi * aesv(j,1,ij_ray,ik_ray) * r(j) * r(j) * r(j) )   &
&                   * ergfoe - a1(j,1)

!-----------------------------------------------------------------------
!  unulsj: Total neutrino energy loss from zone boundary j (foes)
!-----------------------------------------------------------------------

    unulsj          = ( unujrad(j,1,ij_ray,ik_ray) + unujrad(j,2,ij_ray,ik_ray) &
&                   + unujrad(j,3,ij_ray,ik_ray) ) * ergfoe

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

    WRITE (nprint,407) j, gemev, ejhmev, ukemev, pdemev, aesv(j,1,ij_ray,ik_ray), &
&    a1(j,1), a1(j,2), fpirme, unulsj, a1(j,3)

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!                  Shock data                                          !
!                                                                      !
!----------------------------------------------------------------------!
! j0:         s(j0-1) > 1.4*s(j0) and s(j-1) < s(j) for j > j0         !
! j1:         s(j1) < 1.3*s(j1+1) and s(j) > 1.3*s(j+1) for            !
!               j1 < j < j0                                            !
! s0:         ( s(j0+1) + s(j0+2) + s(j0+3) )/3.                       !
! s1:         ( s(j1-1) + s(j1-2) + s(j1-3) )/3.                       !
! s01:        ( s0 + s1 )/2.                                           !
! js0:        Zone index immediately ahead of shock middle             !
! js1:        Zone index immediately behind shock middle               !
! rsp:        Interpolated radius at which s = s01 (cm)                !
! rhosh:      Interpolated density at which s = s01 (g/cm3)            !
! msh:        Interpolated enclosed rest mass at which s = s01         !
!               (solar masses)                                         !
! u0:         Matter velocity immediately ahead of shock; u0=u(j0)     !
!               (cm/sec)                                               !
! u1:         Matter velocity immediately behind shock; u1=u(j1)       !
!               (cm/sec)                                               !
! uesc:       Escape velocity of matter at zone j1 (cm/sec)            !
! lsh:        Rate of work done by shock on preshocked matter          !
!               (foes/sec)                                             !
! shkrmp:     Shock ram pressure ( 4.*pi*r(j1)**3*p(j1) )              !
!               (foes)                                                 !
! ramp0:      Shock ram pressure necessary for u1=0. (foes)            !
! rampes:     Shock ram pressure necessary for u1=uesc (foes)          !
!----------------------------------------------------------------------!        

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(5)           = nedu(5) + 1
  IF ( nedu(5) >= intedu(5) ) THEN
    print_sub       = .true.
    nedu(5)         = 0
  END IF ! nedu(5) >= intedu(5)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( (jr_max - jr_min) < 10 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,501)
  WRITE (nprint,503)

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

  lshock           = .false.

!-----------------------------------------------------------------------
!  j0
!-----------------------------------------------------------------------

  jmaxm            = jr_max - 1
  DO j = jmaxm,jr_min,-1
    j0             = j + 1
    IF ( pq_x(j,ij_ray,ik_ray)/aesv(j,1,ij_ray,ik_ray) >= pqcrit ) THEN
      lshock = .true.
      EXIT
    END IF
  END DO
  
  IF ( .not. lshock ) WRITE (nprint,505)

!-----------------------------------------------------------------------
!  s0
!-----------------------------------------------------------------------

  IF ( lshock ) THEN
    j0p1           = MIN( j0 + 1, jr_max )
    j0p2           = MIN( j0 + 2, jr_max )
    j0p3           = MIN( j0 + 3, jr_max )
    s0             = third * ( aesv(j0p1,3,ij_ray,ik_ray)                    &
&                  + aesv(j0p2,3,ij_ray,ik_ray) + aesv(j0p3,3,ij_ray,ik_ray) )

!-----------------------------------------------------------------------
!  j1
!-----------------------------------------------------------------------

    lshock         = .false.
    DO j = jr_min+3,jr_max
      j1           = j
      IF ( pq_x(j,ij_ray,ik_ray)/aesv(j,1,ij_ray,ik_ray) >= pqcrit ) THEN
        lshock     = .true.
        EXIT
      END IF
    END DO
   IF ( .not. lshock ) WRITE (nprint,507)
  END IF ! lshock

  IF ( lshock ) THEN
    j1m1           = max( j1 - 1, jr_min )
    j1m2           = max( j1 - 2, jr_min )
    j1m3           = max( j1 - 3, jr_min )
    s1             = third * ( aesv(j1m1,3,ij_ray,ik_ray)                    &
&                  + aesv(j1m2,3,ij_ray,ik_ray) + aesv(j1m3,3,ij_ray,ik_ray) )

!-----------------------------------------------------------------------
!  s01
!-----------------------------------------------------------------------

    s01            = half * ( s0 + s1 )

!-----------------------------------------------------------------------
!  js0 and js1
!-----------------------------------------------------------------------

    lshock         = .false.
    DO j = j0,j1,-1
      js1          = j
      IF ( aesv(j,3,ij_ray,ik_ray) > s01 ) THEN
        lshock = .true.
        EXIT
      END IF ! aesv(j,3,ij_ray,ik_ray) > s01
    END DO
    IF ( .not. lshock ) WRITE (nprint,509)
  END IF ! lshock

  js0              = js1 + 1

!-----------------------------------------------------------------------
!
!   ||||| Calculate and print shock parameters of lshock = true |||||
!
!-----------------------------------------------------------------------

  IF ( lshock ) THEN

!-----------------------------------------------------------------------
!  Shock location (radius)
!-----------------------------------------------------------------------

    rsp            = rinterp( r(js0), r(js1), aesv(js0,3,ij_ray,ik_ray), s01, &
&                    aesv(js1,3,ij_ray,ik_ray) )
    WRITE (nprint,511) rsp

!-----------------------------------------------------------------------
!  Shock location (density)
!-----------------------------------------------------------------------

    rhosh          = rinterp( rho(js0), rho(js1), aesv(js0,3,ij_ray,ik_ray),  &
&                    s01, aesv(js1,3,ij_ray,ik_ray) )
    WRITE (nprint,513) rhosh

!-----------------------------------------------------------------------
!  Shock location (mass)
!-----------------------------------------------------------------------

    msht           = rstmss(js1) + frpith * ( rsp * rsp * rsp - r(js1) * r(js1) * r(js1) ) &
&                  * rho(js0)/( half * ( gamgr(js1) + gamgr(js0) ) )
    msh            = msht/msolar
    WRITE (nprint,515) msh

!-----------------------------------------------------------------------
!  Matter luminosity encountered by shock
!-----------------------------------------------------------------------

    dush           = u(j1) - u(j0)
    lsh            = ( frpi * r(j0)**2 * dush**3 * rho(j0) ) * ergfoe
    WRITE (nprint,517) lsh

!-----------------------------------------------------------------------
!  Shock ram pressure
!-----------------------------------------------------------------------

    rj1            = r(j1)
    pj1            = half * ( aesv(j1,1,ij_ray,ik_ray) + aesv(j1+1,1,ij_ray,ik_ray) )
    shkrmp         = frpi * rj1**3 * pj1 * ergfoe
    WRITE (nprint,519) shkrmp

!-----------------------------------------------------------------------
!  Shock ram pressure necessary for u1 = 0
!-----------------------------------------------------------------------

    u0             = u(j0)
    rho00          = half * ( rho(j0) + rho(j0+1) )
    rj0            = r(j0)
    ramp0          = frpi * rho00 * rj0**3 * u0**2 * ergfoe
    u1             = u(j1)
    WRITE (nprint,521) ramp0,u0,u1

!-----------------------------------------------------------------------
!  Shock ram pressure necessary for u1 = uescape
!-----------------------------------------------------------------------

    uesc           = a1(j1,4)
    udif           = uesc - u0
    rampes         = frpi * rho00 * rj0**3 * udif**2
    WRITE (nprint,523) rampes,uesc

!-----------------------------------------------------------------------
!  Pre and postshock density
!-----------------------------------------------------------------------

    rho00          = rho(j0)
    rho1           = rho(j1)
    rho1d0         = rho1/rho00
    WRITE (nprint,525) rho00,rho1,rho1d0

!----------------------------------------------------------------------
!  Pre, middle, and postshock values of gamma1, gamma2, gamma3,
!   and gammae
!----------------------------------------------------------------------!

    DO i = 1,4

      IF (i == 1) j = j0
      IF (i == 2) j = js0
      IF (i == 3) j = js1
      IF (i == 4) j = j1
      
      j = MIN( j, jr_max )

      CALL gammaj_x( j, rho(j), t(j), ij_ray, ik_ray, gamma1, gamma2, &
&      gamma3 )
      gammae       = zero

      IF (rho(j) /= rhor(j)) THEN
        CALL eqstt_x( 1, j, ij_ray, ik_ray, rho(j), tr(j), yer(j), pp, &
&        dpddp, dpdtp, dpdyp )
        pl         = LOG(aesv(j,1,ij_ray,ik_ray))
        ppl        = LOG(pp)
        rhol       = LOG(rho(j))
        rhorl      = LOG(rhor(j))
        gammae     = ( pl - ppl )/( rhol - rhorl + epsilon )
      END IF

      IF (i == 1) WRITE (nprint,527) gamma1,gamma2,gamma3,gammae
      IF (i == 2) WRITE (nprint,529) gamma1,gamma2,gamma3,gammae
      IF (i == 3) WRITE (nprint,531) gamma1,gamma2,gamma3,gammae
      IF (i == 4) WRITE (nprint,533) gamma1,gamma2,gamma3,gammae

    END DO

  END IF ! lshock

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!                  Summary of e-neutrino energy data                   !
!                                                                      !
!----------------------------------------------------------------------!
! unetaf      'enu a e'    Net energy given to e-neutrinos by          !
!                           absorption and emission in zone j (foes)   !
! unetsf      'enu nis'    Net energy given to e-neutrinos by          !
!                           nonisoenergetic scattering in zone j       !
!                           (foes)                                     !
! unetpf      'enu p-a'    Net energy given to e-neutrinos by          !
!                           electron-positron pair annihilation in     !
!                           zone j (foes)                              !
! unetvf      'enu v'      Net energy given to e-neutrinos by          !
!                           advection in  zone j (foes)                !
! unettf      'enu tot'    Net energy given to e-neutrinos by all      !
!                           microphysical processes in zone j (foes)   !
! unettr      'enu trsp'   Net energy given to e-neutrinos by          !
!                           transport into zone j (foes)               !
! unersd      'enu resd'   Total energy of e-neutrinos currently       !
!                           residing in zone j (foes)                  !
! r_lume:      'enu lum'    e-neutrino luminosity across j (foes/sec)   !
! rnavee:     'e nu av n'  e-neutrino number averaged energy in zone   !
!                           j (MeV)                                    !
! rnrmse:     'e nu rms n' e-neutrino number rms averaged energy in    !
!                           zone j (MeV)                               !
! rfrmse:     'e nu rms n' e-neutrino flux rms averaged energy in      !
!                           zone j (MeV)                               !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(6)           = nedu(6) + 1
  IF ( nedu(6) >= intedu(6) ) THEN
    print_sub       = .true.
    nedu(6)         = 0
  END IF ! nedu(6) >= intedu(6)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(1) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,601)
  WRITE (nprint,603)
  WRITE (nprint,605)

  DO jv = jr_min,jr_max,idxedu(6)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of e-neutrino energy losses and gains, luminosities and
!   mean energies
!-----------------------------------------------------------------------

    fluxc          = frpi * r(j) * r(j)
    r_lume          = fluxnu(j,1) * fluxc * ergfoe

    unetaf         = unujea (j,1,ij_ray,ik_ray) * ergfoe
    unetsf         = unujnis(j,1,ij_ray,ik_ray) * ergfoe
    unetpf         = unujpa (j,1,ij_ray,ik_ray) * ergfoe
    unetvf         = unujv  (j,1,ij_ray,ik_ray) * ergfoe
    unettf         = unujt  (j,1,ij_ray,ik_ray) * ergfoe
    unettr         = ( unujrad(j-1,1,ij_ray,ik_ray) - unujrad(j,1,ij_ray,ik_ray) ) * ergfoe
    unersd         = unujinfty(j,1) * ergfoe

    denn           = zero
    dene           = zero
    dene2          = zero
    dfnn           = zero
    dfne2          = zero

    DO k = 1,nnugp(1)
      denn         = denn  + unu(j,k)**2 * dunu(j,k) * psi0(j,k,1)
      dene         = dene  + unu(j,k)**3 * dunu(j,k) * psi0(j,k,1)
      dene2        = dene2 + unu(j,k)**4 * dunu(j,k) * psi0(j,k,1)
      dfnn         = dfnn  + unu(j,k)**2 * dunu(j,k) * psi1(j,k,1)
      dfne2        = dfne2 + unu(j,k)**4 * dunu(j,k) * psi1(j,k,1)
    END DO

    rnavee         = dene/( denn + epsilon )
    rnrmse         = DSQRT( dene2/( denn + epsilon ) + epsilon )
    rfrmse         = DSQRT( DABS( dfne2/( dfnn + epsilon ) ) + epsilon )

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

    WRITE (nprint,607) j,unetaf,unetsf,unetpf,unetvf,unettf,unettr,unersd,r_lume,rnavee,rnrmse,rfrmse

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!                  Summary of e-antineutrino energy data               !
!                                                                      !
!----------------------------------------------------------------------!
! unetaf      'anu a e'    Net energy given to e-antineutrinos by      !
!                           absorption and emission in zone j (foes)   !
! unetsf      'anu nis'    Net energy given to e-antineutrinos by      !
!                           nonisoenergetic scattering in zone j       !
!                           (foes)                                     !
! unetpf      'anu p-a'    Net energy given to e-antineutrinos by      !
!                           electron-positron pair annihilation in     !
!                           zone j (foes)                              !
! unetvf      'anu v'      Net energy given to e-antineutrinos by      !
!                           convecting with material in  zone j        !
!                          (foes)                                      !
! unettf      'anu tot'    Net energy given to e-antineutrinos by all  !
!                           microphysical processes in zone j (foes)   !
! unettr      'anu trsp'   Net energy given to e-antineutrinos by      !
!                           transport into zone j (foes)               !
! unersd      'anu resd'   Total energy of e-antineutrinos currently   !
!                           residing in zone j (foes)                  !
! r_lume:      'anu lum'    e-antineutrino luminosity across j          !
!                           (foes/sec)                                 !
! rnevee:     'e nu av n'  e-antineutrino number averaged energy in    !
!                           zone j (MeV)                               !
! rnrmse:     'e nu rms n' e-neutrino number rms averaged energy in    !
!                           zone j (MeV)                               !
! rfrmse:     'e nu rms n' e-antineutrino flux rms averaged energy in  !
!                           zone j (MeV)                               !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(7)           = nedu(7) + 1
  IF ( nedu(7) >= intedu(7) ) THEN
    print_sub       = .true.
    nedu(7)         = 0
  END IF ! nedu(7) >= intedu(7)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(2) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,611)
  WRITE (nprint,613)
  WRITE (nprint,615)

  DO jv = jr_min,jr_max,idxedu(7)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of e-antineutrino energy losses and gains, luminosities and
!   mean energies
!-----------------------------------------------------------------------

    fluxc          = frpi * r(j) * r(j)
    r_lume          = fluxnu(j,2) * fluxc * ergfoe

    unetaf         = unujea (j,2,ij_ray,ik_ray) * ergfoe
    unetsf         = unujnis(j,2,ij_ray,ik_ray) * ergfoe
    unetpf         = unujpa (j,2,ij_ray,ik_ray) * ergfoe
    unetvf         = unujv  (j,2,ij_ray,ik_ray) * ergfoe
    unettf         = unujt  (j,2,ij_ray,ik_ray) * ergfoe
    unettr         = ( unujrad(j-1,2,ij_ray,ik_ray) - unujrad(j,2,ij_ray,ik_ray) ) * ergfoe
    unersd         = unujinfty(j,2) * ergfoe

    denn           = zero
    dene           = zero
    dene2          = zero
    dfnn           = zero
    dfne2          = zero

    DO k = 1,nnugp(2)
      denn         = denn  + unu(j,k)**2 * dunu(j,k) * psi0(j,k,2)
      dene         = dene  + unu(j,k)**3 * dunu(j,k) * psi0(j,k,2)
      dene2        = dene2 + unu(j,k)**4 * dunu(j,k) * psi0(j,k,2)
      dfnn         = dfnn  + unu(j,k)**2 * dunu(j,k) * psi1(j,k,2)
      dfne2        = dfne2 + unu(j,k)**4 * dunu(j,k) * psi1(j,k,2)
    END DO

    rnavee         = dene/( denn + epsilon )
    rnrmse         = DSQRT( dene2/( denn + epsilon ) + epsilon )
    rfrmse         = DSQRT( DABS( dfne2/( dfnn + epsilon ) ) + epsilon )

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

    WRITE (nprint,617) j,unetaf,unetsf,unetpf,unetvf,unettf,unettr,unersd,r_lume,rnavee,rnrmse,rfrmse

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!                  Summary of t-neutrino energy data                   !
!                                                                      !
!----------------------------------------------------------------------!
! unetaf      'tnu a e'    Net energy given to t-neutrinos by          !
!                           absorption and emission in zone j (foes)   !
! unetsf      'tnu nis'    Net energy given to t-neutrinos by          !
!                           nonisoenergetic scattering in zone j       !
!                           (foes)                                     !
! unetpf      'tnu p-a'    Net energy given to t-neutrinos by          !
!                           electron-positron pair annihilation in     !
!                           zone j (foes)                              !
! unetvf      'tnu v'      Net energy given to t-neutrinos by          !
!                           convecting with material in  zone j        !
!                          (foes)                                      !
! unettf      'tnu tot'    Net energy given to t-neutrinos by all      !
!                           microphysical processes in zone j (foes))  !
! unettr      'tnu trsp'   Net energy given to t-neutrinos by          !
!                           transport into zone j (foes)               !
! unersd      'tnu resd'   Total energy of t-neutrinos currently       !
!                           residing in zone j (foes))                 !
! r_lumt:      'tnu lum'    t-neutrino luminosity across j (foes/sec)   !
! rnavee:     'e nu av n'  t-neutrino number averaged energy in zone   !
!                           j (MeV)                                    !
! rnrmse:     'e nu rms n' t-neutrino number rms averaged energy in    !
!                           zone j (MeV)                               !
! rfrmse:     'e nu rms n' t-neutrino flux rms averaged energy in      !
!                           zone j (MeV)                               !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(8)           = nedu(8) + 1
  IF ( nedu(8) >= intedu(8) ) THEN
    print_sub       = .true.
    nedu(8)         = 0
  END IF ! nedu(8) >= intedu(8)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(3) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,621)
  WRITE (nprint,623)
  WRITE (nprint,625)

  DO jv = jr_min,jr_max,idxedu(8)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of x-neutrino energy losses and gains, luminosities and
!   mean energies
!-----------------------------------------------------------------------

    fluxc          = frpi * r(j) * r(j)
    r_lume          = fluxnu(j,3) * fluxc * ergfoe

    unetaf         = unujea (j,3,ij_ray,ik_ray) * ergfoe
    unetsf         = unujnis(j,3,ij_ray,ik_ray) * ergfoe
    unetpf         = unujpa (j,3,ij_ray,ik_ray) * ergfoe
    unetvf         = unujv  (j,3,ij_ray,ik_ray) * ergfoe
    unettf         = unujt  (j,3,ij_ray,ik_ray) * ergfoe
    unettr         = ( unujrad(j-1,3,ij_ray,ik_ray) - unujrad(j,3,ij_ray,ik_ray) ) * ergfoe
    unersd         = unujinfty(j,3) * ergfoe

    denn           = zero
    dene           = zero
    dene2          = zero
    dfnn           = zero
    dfne2          = zero

    DO k = 1,nnugp(3)
      denn         = denn  + unu(j,k)**2 * dunu(j,k) * psi0(j,k,3)
      dene         = dene  + unu(j,k)**3 * dunu(j,k) * psi0(j,k,3)
      dene2        = dene2 + unu(j,k)**4 * dunu(j,k) * psi0(j,k,3)
      dfnn         = dfnn  + unu(j,k)**2 * dunu(j,k) * psi1(j,k,3)
      dfne2        = dfne2 + unu(j,k)**4 * dunu(j,k) * psi1(j,k,3)
    END DO

    rnavee         = dene/( denn + epsilon )
    rnrmse         = DSQRT( dene2/( denn + epsilon ) + epsilon )
    rfrmse         = DSQRT( DABS( dfne2/( dfnn + epsilon ) ) + epsilon )

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

    WRITE (nprint,627) j,unetaf,unetsf,unetpf,unetvf,unettf,unettr,unersd,r_lume,rnavee,rnrmse,rfrmse

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!                  Summary of t-antineutrino energy data               !
!                                                                      !
!----------------------------------------------------------------------!
! unetaf      'tnu a e'    Net energy given to t-antineutrinos by      !
!                           absorption and emission in zone j (foes)   !
! unetsf      'tnu nis'    Net energy given to t-antineutrinos by      !
!                           nonisoenergetic scattering in zone j       !
!                           (foes)                                     !
! unetpf      'tnu p-a'    Net energy given to t-antineutrinos by      !
!                           electron-positron pair annihilation in     !
!                           zone j (foes)                              !
! unetvf      'tnu v'      Net energy given to t-antineutrinos by      !
!                           convecting with material in  zone j        !
!                          (foes)                                      !
! unettf      'tnu tot'    Net energy given to t-antineutrinos by all  !
!                           microphysical processes in zone j (foes))  !
! unettr      'tnu trsp'   Net energy given to t-antineutrinos by      !
!                           transport into zone j (foes)               !
! unersd      'tnu resd'   Total energy of t-antineutrinos currently   !
!                           residing in zone j (foes))                 !
! r_lumt:      'tnu lum'    t-antineutrino luminosity across j         !
!                          (foes/sec)                                  !
! rnavee:     'e nu av n'  t-antineutrino number averaged energy in    !
!                           zone j (MeV)                               !
! rnrmse:     'e nu rms n' t-antineutrino number rms averaged energy   !
!                           in zone j (MeV)                            !
! rfrmse:     'e nu rms n' t-antineutrino flux rms averaged energy in  !
!                           zone j (MeV)                               !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(9)           = nedu(9) + 1
  IF ( nedu(9) >= intedu(9) ) THEN
    print_sub       = .true.
    nedu(9)         = 0
  END IF ! nedu(9) >= intedu(9)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(4) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,631)
  WRITE (nprint,633)
  WRITE (nprint,635)

  DO jv = jr_min,jr_max,idxedu(9)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of x-antineutrino energy losses and gains, luminosities and
!   mean mean energies
!-----------------------------------------------------------------------

    fluxc          = frpi * r(j) * r(j)
    r_lume          = fluxnu(j,4) * fluxc * ergfoe

    unetaf         = unujea (j,4,ij_ray,ik_ray) * ergfoe
    unetsf         = unujnis(j,4,ij_ray,ik_ray) * ergfoe
    unetpf         = unujpa (j,4,ij_ray,ik_ray) * ergfoe
    unetvf         = unujv  (j,4,ij_ray,ik_ray) * ergfoe
    unettf         = unujt  (j,4,ij_ray,ik_ray) * ergfoe
    unettr         = ( unujrad(j-1,4,ij_ray,ik_ray) - unujrad(j,4,ij_ray,ik_ray) ) * ergfoe
    unersd         = unujinfty(j,4) * ergfoe

    denn           = zero
    dene           = zero
    dene2          = zero
    dfnn           = zero
    dfne2          = zero

    DO k = 1,nnugp(4)
      denn         = denn  + unu(j,k)**2 * dunu(j,k) * psi0(j,k,4)
      dene         = dene  + unu(j,k)**3 * dunu(j,k) * psi0(j,k,4)
      dene2        = dene2 + unu(j,k)**4 * dunu(j,k) * psi0(j,k,4)
      dfnn         = dfnn  + unu(j,k)**2 * dunu(j,k) * psi1(j,k,4)
      dfne2        = dfne2 + unu(j,k)**4 * dunu(j,k) * psi1(j,k,4)
    END DO

    rnavee         = dene/( denn + epsilon )
    rnrmse         = DSQRT( dene2/( denn + epsilon ) + epsilon )
    rfrmse         = DSQRT( DABS( dfne2/( dfnn + epsilon ) ) + epsilon )

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

    WRITE (nprint,637) j,unetaf,unetsf,unetpf,unetvf,unettf,unettr,unersd,r_lume,rnavee,rnrmse,rfrmse

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!                  Matter-e-neutrino energy transfer rates             !
!                                                                      !
!----------------------------------------------------------------------!
! uneta:      'enur a e'   Energy transfer rate to e-neutrinos         !
!                           by absorption and emission in zone j       !
!                           (ergs/gm*sec)                              !
! unets:      'enur nes'   Energy transfer rate to e-neutrinos         !
!                           by nonconservative scattering in zone j    !
!                           (ergs/gm*sec)                              !
! unetp:      'enur p-a'   Energy transfer rate to e-neutrinos         !
!                           by electron-positron pair annihilation     !
!                           in zone j (ergs/gm*sec)                    !
! unett:      'enur tot'   Energy transfer rate to e-neutrinos         !
!                           by all matter-neutrino interactions        !
!                           (ergs/gm*sec)                              !
! unetv:      'enur v'     Energy transfer rate to e-neutrinos         !
!                           by advection in zone j (ergs/gm*sec)       !
! unetr:      'enur trsp'  Energy transfer rate to e-neutrinos         !
!                           by transport in zone j (ergs/gm*sec)       !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(10)          = nedu(10) + 1
  IF ( nedu(10) >= intedu(10) ) THEN
    print_sub       = .true.
    nedu(10)         = 0
  END IF ! nedu(10) >= intedu(10)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(1) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,701)
  WRITE (nprint,703)
  WRITE (nprint,705)

  DO jv = jr_min,jr_max,idxedu(10)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of e-neutrino energy transfer rates
!-----------------------------------------------------------------------

    coefr          = one/dmgrv(j)
    uneta          = coefr * dunujeadt (j,1,ij_ray,ik_ray)
    unets          = coefr * dunujnisdt(j,1,ij_ray,ik_ray)
    unetp          = coefr * dunujpadt (j,1,ij_ray,ik_ray)
    unett          =         dunujtdt  (j,1)
    unetv          = coefr * dunujvdt  (j,1,ij_ray,ik_ray)
    fluxcj         = frpi * r(j  ) * r(j  )
    fluxcjm1       = frpi * r(j-1) * r(j-1)
    unetr          = ( fluxcjm1 * fluxnu(j-1,1) - fluxcj * fluxnu(j,1) )/dmgrv(j)

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

    WRITE (nprint,707) j,uneta,unets,unetp,unett,unetv,unetr

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!        Matter-e-antineutrino energy transfer rates                   !
!                                                                      !
!----------------------------------------------------------------------!
! uneta:      'anur a e'   Energy transfer rate to e-antineutrinos     !
!                           by absorption and emission in zone j       !
!                           (ergs/gm*sec)                              !
! unets:      'anur nes'   Energy transfer rate to e-antineutrinos     !
!                           by nonconservative scattering in zone j    !
!                           (ergs/gm*sec)                              !
! unetp:      'anur p-a'   Energy transfer rate to e-antineutrinos     !
!                           by electron-positron pair annihilation     !
!                           in zone j (ergs/gm*sec)                    !
! unett:      'anur tot'   Energy transfer rate to e-antineutrinos     !
!                           by all matter-neutrino interactions        !
!                           (ergs/gm*sec)                              !
! unetv:      'anur v'     Energy transfer rate to e-antineutrinos     !
!                           by advection in zone j (ergs/gm*sec)       !
! unetr:      'anur trsp'  Energy transfer rate to e-antineutrinos     !
!                           by transport in zone j (ergs/gm*sec)       !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(11)          = nedu(11) + 1
  IF ( nedu(11) >= intedu(11) ) THEN
    print_sub       = .true.
    nedu(11)        = 0
  END IF ! nedu(11) >= intedu(11)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(2) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,711)
  WRITE (nprint,713)
  WRITE (nprint,715)

  DO jv = jr_min,jr_max,idxedu(11)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of e-antineutrino energy transfer rates
!-----------------------------------------------------------------------

    coefr          = one/dmgrv(j)
    uneta          = coefr * dunujeadt (j,2,ij_ray,ik_ray)
    unets          = coefr * dunujnisdt(j,2,ij_ray,ik_ray)
    unetp          = coefr * dunujpadt (j,2,ij_ray,ik_ray)
    unett          =         dunujtdt  (j,2)
    unetv          = coefr * dunujvdt  (j,2,ij_ray,ik_ray)
    fluxcj         = frpi * r(j  ) * r(j  )
    fluxcjm1       = frpi * r(j-1) * r(j-1)
    unetr          = ( fluxcjm1 * fluxnu(j-1,2) - fluxcj * fluxnu(j,2) )/dmgrv(j)

    WRITE (nprint,717) j,uneta,unets,unetp,unett,unetv,unetr

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!        Matter-t-neutrino energy transfer rates                       !
!                                                                      !
!----------------------------------------------------------------------!
! uneta:      'tnur a e'   Energy transfer rate to t-neutrinos         !
!                           by absorption and emission in zone j       !
!                           (ergs/gm*sec)                              !
! unets:      'tnur nes'   Energy transfer rate to t-neutrinos         !
!                           by nonconservative scattering in zone j    !
!                           (ergs/gm*sec)                              !
! unetp:      'tnur p-a'   Energy transfer rate to t-neutrinos         !
!                           by electron-positron pair annihilation     !
!                           in zone j (ergs/gm*sec)                    !
! unett:      'tnur tot'   Energy transfer rate to t-neutrinos         !
!                           by all matter-neutrino interactions        !
!                           (ergs/gm*sec)                              !
! unetv:      'tnur v'     Energy transfer rate to t-neutrinos         !
!                           by advection in zone j (ergs/gm*sec)       !
! unetr:      'tnur trsp'  Energy transfer rate to t-neutrinos         !
!                           by transport in zone j (ergs/gm*sec)       !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(12)          = nedu(12) + 1
  IF ( nedu(12) >= intedu(12) ) THEN
    print_sub       = .true.
    nedu(12)        = 0
  END IF ! nedu(12) >= intedu(12)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(3) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,721)
  WRITE (nprint,723)
  WRITE (nprint,725)

  DO jv = jr_min,jr_max,idxedu(12)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of x-neutrino energy transfer rates
!-----------------------------------------------------------------------

    coefr          = one/dmgrv(j)
    uneta          = coefr * dunujeadt (j,3,ij_ray,ik_ray)
    unets          = coefr * dunujnisdt(j,3,ij_ray,ik_ray)
    unetp          = coefr * dunujpadt (j,3,ij_ray,ik_ray)
    unett          =         dunujtdt  (j,3)
    unetv          = coefr * dunujvdt  (j,3,ij_ray,ik_ray)
    fluxcj         = frpi * r(j  ) * r(j  )
    fluxcjm1       = frpi * r(j-1) * r(j-1)
    unetr          = ( fluxcjm1 * fluxnu(j-1,3) - fluxcj * fluxnu(j,3) )/dmgrv(j)

    WRITE (nprint,727) j,uneta,unets,unetp,unett,unetv,unetr

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!        Matter-t-antineutrino energy transfer rates                   !
!                                                                      !
!----------------------------------------------------------------------!
! uneta:      'tnur a e'   Energy transfer rate to t-antineutrino      !
!                           by absorption and emission in zone j       !
!                           (ergs/gm*sec)                              !
! unets:      'tnur nes'   Energy transfer rate to t-antineutrinos     !
!                           by nonconservative scattering in zone j    !
!                           (ergs/gm*sec)                              !
! unetp:      'tnur p-a'   Energy transfer rate to t-antineutrinos     !
!                           by electron-positron pair annihilation     !
!                           in zone j (ergs/gm*sec)                    !
! unett:      'tnur tot'   Energy transfer rate to t-antineutrinos     !
!                           by all matter-neutrino interactions        !
!                           (ergs/gm*sec)                              !
! unetv:      'tnur v'     Energy transfer rate to t-antineutrinos     !
!                           by advection in zone j (ergs/gm*sec)       !
! unetr:      'tnur trsp'  Energy transfer rate to t-antineutrinos     !
!                           by transport in zone j (ergs/gm*sec)       !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(13)          = nedu(13) + 1
  IF ( nedu(13) >= intedu(13) ) THEN
    print_sub       = .true.
    nedu(13)        = 0
  END IF ! nedu(13) >= intedu(13)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugp(4) == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,731)
  WRITE (nprint,733)
  WRITE (nprint,735)

  DO jv = jr_min,jr_max,idxedu(13)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Summary of x-antineutrino energy transfer rates
!-----------------------------------------------------------------------

    coefr          = one/dmgrv(j)
    uneta          = coefr * dunujeadt (j,4,ij_ray,ik_ray)
    unets          = coefr * dunujnisdt(j,4,ij_ray,ik_ray)
    unetp          = coefr * dunujpadt (j,4,ij_ray,ik_ray)
    unett          =         dunujtdt  (j,4)
    unetv          = coefr * dunujvdt  (j,4,ij_ray,ik_ray)
    fluxcj         = frpi * r(j  ) * r(j  )
    fluxcjm1       = frpi * r(j-1) * r(j-1)
    unetr          = ( fluxcjm1 * fluxnu(j-1,4) - fluxcj * fluxnu(j,4) )/dmgrv(j)

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

    WRITE (nprint,737) j,uneta,unets,unetp,unett,unetv,unetr

  END DO

END IF ! print_sub

!**********************************************************************!
!                                                                      !
!        Summary of Neutrino Luminosity data                           !
!                                                                      !
!----------------------------------------------------------------------!
! rnu(i):     'r'          Radius at which neutrino luminosity data    !
!                           is evaluated (cm)                          !
! rhonu(i):   'rho'        Density at which neutrino luminosity data   !
!                           is evaluated (gm/cm**3)                    !
! lum(i,1)    'enu lum'    e-neutrino luminosity at rnu(i) (foes/sec)  !
! lum(i,2)    'anu lum'    e-antineutrino luminosity at rnu(i)         !
!                           (foes/sec)                                 !
! lum(i,3)    'tnu lum'    t-neutrino luminosity at rnu(i) (foes/sec)  !
! lum(i,4)    'tbnu lum'   t-antineutrino luminosity at rnu(i)         !
!                           (foes/sec)                                 !
! enuv(i,1):  'enu eav1'   e-neutrino number averaged energy at rnu(i) !
!                           (MeV)                                      !
! enuv(i,2):  'anu eav1'   e-antineutrino number averaged energy at    !
!                           rnu(i) (MeV)                               !
! enuv(i,3):  'tnu eav1'   t-neutrino number averaged energy at rnu(i) !
!                           (MeV)                                      !
! enuv(i,4):  'tbnu eav1'  t-antineutrino number averaged energy at    !
!                            rnu(i)(MeV)                               !
! enuvrms(i,1)'enu eav2'   e-neutrino number rms averaged energy at    !
!                           rnu(i) (MeV)                               !
! enuvrms(i,2)'anu eav2'   e-antineutrino number rms averaged energy   !
!                           at rnu(i) (MeV)                            !
! enuvrms(i,3)'tnu eav2'   t-neutrino number rms averaged energy at    !
!                           rnu(i) (MeV)                               !
! enuvrms(i,4)'tbnu eav2'  t-antineutrino number rms averaged energy   !
!                           at rnu(i) (MeV)                            !
! rnu(i):     'r'          Radius at which neutrino luminosity data    !
!                           is evaluated (cm)                          !
! rhonu(i):   'rho'        Density at which neutrino luminosity data   !
!                           is evaluated (gm/cm**3)                    !
! rmass(i):   'M (solar)'  Rest mass enclosed by rnu(i) (solar masses) !
! tnu(i):     'T (MeV)'    Temperature at rnu(i) (MeV)                 !
! snu(i):     's'          Entropy at rnu(i)                           !
!----------------------------------------------------------------------!

print_sub           = .false.

IF ( prnttest ) THEN
  nedu(14)          = nedu(14) + 1
  IF ( nedu(14) >= intedu(14) ) THEN
    print_sub       = .true.
    nedu(14)        = 0
  END IF ! nedu(14) >= intedu(14)
ELSE
  print_sub         = .true.
END IF ! prnttest

IF ( nnugpmx == 0 ) print_sub = .false.

IF ( print_sub ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,801)
  WRITE (nprint,803)
  WRITE (nprint,805) time
  WRITE (nprint,807)


!-----------------------------------------------------------------------
!
!         ||||| Neutrino luminosity and mean energy data |||||
!
!-----------------------------------------------------------------------

outer1: DO i = 1,7


!-----------------------------------------------------------------------
!  Find jd such that rho(jd) < rhonud(i) and rho(jd-1) > rhonud(i)
!  IF jd cannot be found, try the next index i
!-----------------------------------------------------------------------

    IF ( rho(2) < rhonud(i) ) THEN
      WRITE (nprint,8001) i,rhonud(i)
      EXIT
    END IF ! rho(2) < rhonud(i)

    DO j = 2,jr_max
      IF ( rho(j) <= rhonud(i) ) THEN
        jd         = j
        EXIT
      END IF ! rho(j) le rhonud(i)
      IF ( j == jr_max ) THEN
        WRITE (nprint,8001) i,rhonud(i)
        CYCLE outer1
      END IF
    END DO

!-----------------------------------------------------------------------
!  Density
!-----------------------------------------------------------------------

    rhonu(i)       = rhonud(i)

!-----------------------------------------------------------------------
!  Radius
!-----------------------------------------------------------------------

    rnu(i)         = rinterp( r(jd), r(jd-1), rho(jd), rhonud(i), rho(jd-1) )

!-----------------------------------------------------------------------
!  Mass
!-----------------------------------------------------------------------

    rmassjd        = rstmss(jd)/msolar
    rmassjdm1      = rstmss(jd-1)/msolar
    rmass(i)       = rinterp( rmassjd, rmassjdm1, rho(jd), rhonud(i), rho(jd-1) )

!-----------------------------------------------------------------------
!  Temperature
!-----------------------------------------------------------------------

    tnujd          = t(jd) * kmev
    tnujdm1        = t(jd-1) * kmev
    tnu(i)         = rinterp( tnujd, tnujdm1, rho(jd), rhonud(i), rho(jd-1) )

!-----------------------------------------------------------------------
!  Entropy
!-----------------------------------------------------------------------

    snu(i)         = rinterp( aesv(jd,3,ij_ray,ik_ray), aesv(jd-1,3,ij_ray,ik_ray), &
 &                   rho(jd), rhonud(i), rho(jd-1) )

!-----------------------------------------------------------------------
!  Neutrino luminosities
!-----------------------------------------------------------------------

    fluxcjd        = frpi * r(jd  ) * r(jd  )
    fluxcjdm       = frpi * r(jd-1) * r(jd-1)

    DO n = 1,nnu

!-----------------------------------------------------------------------
!  Set luminosities to zero if nnugp(n) = 0
!-----------------------------------------------------------------------

      IF ( nnugp(n) == 0 ) THEN
        lum(i,n)   = zero
        CYCLE
      END IF

!-----------------------------------------------------------------------
!  Calculate the luminosities
!-----------------------------------------------------------------------

      r_lumjd       = fluxnu(jd  ,n) * fluxcjd  * ergfoe
      r_lumjdm      = fluxnu(jd-1,n) * fluxcjdm * ergfoe
      lum(i,n)     = rinterp( r_lumjd, r_lumjdm, rho(jd), rhonud(i), rho(jd-1) )

    END DO

!-----------------------------------------------------------------------
!  Neutrino mean energies
!-----------------------------------------------------------------------

inner1:  DO n = 1,nnu

      dennjd       = zero
      denejd       = zero
      dene2jd      = zero
      dennjdm      = zero
      denejdm      = zero
      dene2jdm     = zero

!-----------------------------------------------------------------------
!  Set the mean energies to zero if nnugp(n) = 0
!-----------------------------------------------------------------------

      IF ( nnugp(n) == 0 ) THEN
        enuv(i,n)    = zero
        enuvrms(i,n) = zero
        CYCLE inner1
      END IF

!-----------------------------------------------------------------------
!  Calculate the mean energies
!-----------------------------------------------------------------------
    
      DO k = 1,nnugp(n)
        coefn      = unu(jd,k) * unu(jd,k) * dunu(jd,k)
        coefe      = coefn * unu(jd,k)
        coefe2     = coefe * unu(jd,k)
        coefnm     = unu(jd-1,k) * unu(jd-1,k) * dunu(jd-1,k)
        coefem     = coefnm * unu(jd-1,k)
        coefe2m    = coefem * unu(jd-1,k)
        dennjd     = dennjd    + coefn   * psi0(jd  ,k,n)
        denejd     = denejd    + coefe   * psi0(jd  ,k,n)
        dene2jd    = dene2jd   + coefe2  * psi0(jd  ,k,n)
        dennjdm    = dennjdm   + coefnm  * psi0(jd-1,k,n)
        denejdm    = denejdm   + coefem  * psi0(jd-1,k,n)
        dene2jdm   = dene2jdm  + coefe2m * psi0(jd-1,k,n)
      END DO

      enuvjd       = denejd /( dennjd   + epsilon )
      enuvjdm      = denejdm/( dennjdm  + epsilon )
      enuv2jd      = dene2jd /( dennjd   + epsilon )
      enuv2jdm     = dene2jdm/( dennjdm  + epsilon )
      enuv(i,n)    = rinterp( enuvjd, enuvjdm, rho(jd), rhonud(i), rho(jd-1) )
      enuv2        = rinterp( enuv2jd, enuv2jdm, rho(jd), rhonud(i), rho(jd-1) )
      enuvrms(i,n) = DSQRT( abs(enuv2) + epsilon )

    END DO inner1

  END DO outer1

outer2: DO i = 1,7

!-----------------------------------------------------------------------
!  Find jd such that r(jd) > rnud(i) and r(jd-1) < rnud(i)
!  IF jd cannot be found, try the next index i
!-----------------------------------------------------------------------

    DO j = 2,jr_max
      IF ( r(j) >= rnud(i) ) THEN
        jd         = j
        EXIT
       END IF !  r(j) ge rnud(i)
       IF ( j == jr_max ) THEN
        WRITE (nprint,8101) i,rnud(i)
        CYCLE outer2
      END IF
    END DO

!-----------------------------------------------------------------------
!  Radius
!-----------------------------------------------------------------------

    rnu(i+7)       = rnud(i)

!-----------------------------------------------------------------------
!  Density
!-----------------------------------------------------------------------

    rhonu(i+7)     = rinterp( rho(jd), rho(jd-1), r(jd), rnud(i), r(jd-1) )

!-----------------------------------------------------------------------
!  Mass
!-----------------------------------------------------------------------

    rmassjd        = rstmss(jd)/msolar
    rmassjdm1      = rstmss(jd-1)/msolar
    rmass(i+7)     = rinterp( rmassjd, rmassjdm1, r(jd), rnud(i), r(jd-1) )

!-----------------------------------------------------------------------
!  Temperature
!-----------------------------------------------------------------------

    tnujd          = t(jd) * kmev
    tnujdm1        = t(jd-1) * kmev
    tnu(i+7)       = rinterp( tnujd, tnujdm1, r(jd), rnud(i), r(jd-1) )

!-----------------------------------------------------------------------
!  Entropy
!-----------------------------------------------------------------------

    snu(i+7)       = rinterp( aesv(jd,3,ij_ray,ik_ray), aesv(jd-1,3,ij_ray,ik_ray), &
&                    r(jd), rnud(i), r(jd-1) )

!-----------------------------------------------------------------------
!  Neutrino luminosities
!-----------------------------------------------------------------------

    fluxcjd        = frpi * r(jd  ) * r(jd  )
    fluxcjdm       = frpi * r(jd-1) * r(jd-1)

    DO n = 1,nnu

!-----------------------------------------------------------------------
!  Set luminosities to zero if nnugp(n) = 0
!-----------------------------------------------------------------------

      IF ( nnugp(n) == 0 ) THEN
        lum(i+7,n) = zero
        CYCLE
      END IF

!-----------------------------------------------------------------------
!  Calculate the luminosities
!-----------------------------------------------------------------------

      r_lumjd       = fluxnu(jd  ,n) * fluxcjd  * ergfoe
      r_lumjdm      = fluxnu(jd-1,n) * fluxcjdm * ergfoe
      lum(i+7,n) = rinterp( r_lumjd, r_lumjdm, r(jd), rnud(i), r(jd-1) )

      END DO

!-----------------------------------------------------------------------
!  Neutrino mean energies
!-----------------------------------------------------------------------

inner2:  DO n = 1,nnu

        dennjd     = zero
        denejd     = zero
        dene2jd    = zero
        dennjdm    = zero
        denejdm    = zero
        dene2jdm   = zero

!-----------------------------------------------------------------------
!  Set the mean energies to zero if nnugp(n) = 0
!-----------------------------------------------------------------------

        IF ( nnugp(n) == 0 ) THEN
          enuv(i+7,n)    = zero
          enuvrms(i+7,n) = zero
          CYCLE inner2
        END IF

!-----------------------------------------------------------------------
!  Calculate the mean energies
!-----------------------------------------------------------------------
    
        DO k = 1,nnugp(n)
          coefn     = unu(jd,k) * unu(jd,k) * dunu(jd,k)
          coefe     = coefn * unu(jd,k)
          coefe2    = coefe * unu(jd,k)
          coefnm    = unu(jd-1,k) * unu(jd-1,k) * dunu(jd-1,k)
          coefem    = coefnm * unu(jd-1,k)
          coefe2m   = coefem * unu(jd-1,k)
          dennjd    = dennjd    + coefn   * psi0(jd  ,k,n)
          denejd    = denejd    + coefe   * psi0(jd  ,k,n)
          dene2jd   = dene2jd   + coefe2  * psi0(jd  ,k,n)
          dennjdm   = dennjdm   + coefnm  * psi0(jd-1,k,n)
          denejdm   = denejdm   + coefem  * psi0(jd-1,k,n)
          dene2jdm  = dene2jdm  + coefe2m * psi0(jd-1,k,n)
        END DO

        enuvjd      = denejd /( dennjd   + epsilon )
        enuvjdm     = denejdm/( dennjdm  + epsilon )
        enuv2jd     = dene2jd /( dennjd   + epsilon )
        enuv2jdm    = dene2jdm/( dennjdm  + epsilon )
        enuv(i+7,n) = rinterp( enuvjd, enuvjdm, r(jd), rnud(i), r(jd-1) )
        enuv2       = rinterp( enuv2jd, enuv2jdm, r(jd), rnud(i), r(jd-1) )
        enuvrms(i+7,n)= DSQRT( ABS(enuv2) + epsilon )

      END DO inner2

  END DO outer2

!-----------------------------------------------------------------------
!  Print
!-----------------------------------------------------------------------

  DO i = 1,14
    WRITE (nprint,809) i,rnu(i),rhonu(i),lum(i,1),lum(i,2),lum(i,3),lum(i,4)
  END DO

  WRITE (nprint,811)
  WRITE (nprint,821)
  DO i = 1,14
    WRITE (nprint,823) i,rnu(i),rhonu(i),enuv(i,1),enuv(i,2),enuv(i,3),enuv(i,4)
  END DO

  WRITE (nprint,811)
  WRITE (nprint,831)
  DO i = 1,14
    WRITE (nprint,833) i,rnu(i),rhonu(i),enuvrms(i,1),enuvrms(i,2),enuvrms(i,3),enuvrms(i,4)
  END DO

  WRITE (nprint,811)
  WRITE (nprint,851)
  DO i = 1,14
    WRITE (nprint,853) i,rnu(i),rhonu(i),rmass(i),tnu(i),snu(i)
  END DO

END IF ! print_sub

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

IF ( l_alloc ) THEN

  DEALLOCATE (u_gr, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_gr      '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (u_nt, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_nt      '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (u_grt, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_grt     '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (u_ntt, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'u_ntt     '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (a1, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'a1        '; WRITE (nlog,2001) var_name; END IF

  DEALLOCATE (lum, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'lum       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (enuv, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'enuv      '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (enuvrms, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'enuvrms   '; WRITE (nlog,2001) var_name; END IF

END IF ! print_sub

RETURN

CONTAINS
REAL (KIND=double) FUNCTION rinterp(avar,bvar,x,y,z)

REAL (KIND=double) :: avar
REAL (KIND=double) :: bvar
REAL (KIND=double) :: x
REAL (KIND=double) :: y
REAL (KIND=double) :: z

rinterp             = bvar + ( avar - bvar ) * ( y - z )/( x - z )
END FUNCTION rinterp

END SUBROUTINE editu
