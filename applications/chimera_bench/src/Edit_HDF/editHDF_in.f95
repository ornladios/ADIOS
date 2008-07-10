FUNCTION editMD_in( imin, imax, nx, js, je, ny, i_ray, i_ray_dim, nez, nnu, &
& nnc, rhop, tp, yep, rhobarp, rp, rpph, thetap, thetapph, up, vp, wp, psi0p, &
& psi1p, dtime, time_elapsed, cycle_number, xnp, be_nuc_repp, a_nuc_repp, &
& z_nuc_repp, nsep, first, e_rms_stat, e_rms_trns, e_rms_r_stat, e_rms_r_trns, &
& e_rms_d_stat, e_rms_d_trns, rsphere_mean, dsphere_mean, tsphere_mean, &
& msphere_mean, esphere_mean, r_sphere, d_sphere, t_sphere, m_sphere, lum, &
& lum_r, lum_rho, inv_fluxfact, stat_pinch_r, trns_pinch_r, stat_pinch_d, &
& trns_pinch_d, vel_x, vel_y, vel_z, xi, xiph, thetai, thetaiph, rho_ed, &
& t_ed, ye_ed, s_ed, yl_ed, dm_ed, yenu, yenubar, yxnu, yxnubar, e_nu_center, &
& e_nu_edge )

!-----------------------------------------------------------------------
!
!    File:         editMD_in
!    Module:       editMD_in
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To calculate quantities for multi-D edit dumps.
!
!    Subprograms called:
!
!  eqstz_x        : updates equation of state variables
!  gammaz_x       : updates EOS gammas
!  yl_cal         : computes the lepton and neutrino fractions
!  psi1_cal       : computes the first moment of the neutrino distribution
!                    function
!  e_rms          : computes the RMS j,n neutrino energies
!  e_rms_surface  : computes the RMS n neutrino energies at radius r_e_rms
!                    and density d_e_rms
!  nu_sphere      : computes the k,n neutrinospheres
!  nu_sphere_mean : computes the mean n neutrinospheres at radius r_lum and
!                    density d_lum
!  luminosity     : computes the j,n neutrino luminosities
!  lum_surface    : computes the luminosities at radius r_lum and density d_lum
!  inv_fluxfactor : computes the j,n neutrino inverse flux factors
!  pinch_surface  : computes the n neutrino pinch parameters at ratius r_pinch
!                    and density d_pinch
!
!    Input arguments:
!
!  imin           : inner x-array index
!  imax           : outer x-array index
!  nx             : x_array extent
!  js             : inner y-array index
!  je             : outer y-array index
!  ny             : y_array extent
!  i_ray          : index denoting a specific radial ray
!  i_ray_dim      : number of rays assigned to a processor
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!  nnc            : neutrino abundance array extent
!  rhop           : density (cm^{-3})
!  tp             : temperature (K)
!  yep            : electron fraction
!  rhobarp        : mean density (cm^{-3})
!  rp             : zone-edged radial coordinate (cm)
!  rpph           : zone-centered radial coordinate (cm)
!  thetap         : zone-edged y (angular) coordinate
!  thetapph       : zone-centered y (angular) coordinate
!  up             : x-component of velocity of zone (cm s^{-1})
!  vp             : y-component of velocity of zone (cm s^{-1})
!  wp             : z-component of velocity of zone (cm s^{-1})
!  psi0p          : zeroth angular moment of the NDF
!  psi1p          : first angular moment of the NDF
!  dtime          : time step used at current cycle
!  time_ellapsed  : ellapsed time
!  cycle_number   : cycle number
!  xnp            : initial mass fractions
!  be_nuc_repp    : binding energy of mean heavy nucleus (MeV)
!  a_nuc_repp     : mass number of mean heavy nucleus
!  z_nuc_repp     : charge number of mean heavy nucleus
!  nesp           : nuclear statistical equilibrium flag
!  first          : initial call flag
!
!    Output arguments:
!
!  e_rms_stat     : sqrt[ SUM psi0 * w5dw/SUM w3ww ] (MeV) (RMS static neutrino energy)
!  e_rms_trns     : sqrt[ SUM psi1 * w5dw/SUM w3ww ] (MeV) (RMS transport neutrino energy)
!  e_rms_r_stat   : rms static neutrino energy at radius r_e_rms (MeV)
!  e_rms_r_trns   : rms transport neutrino energy at radius r_e_rms (MeV)
!  e_rms_d_stat   : rms static neutrino energy at density d_e_rms (MeV)
!  e_rms_d_trns   : rms transport neutrino energy at density d_e_rms (MeV)
!  r_sphere       : k,n neutrinosphere radius (cm)
!  d_sphere       : k,n neutrinosphere density (g cm^{-3})
!  r_sphere       : k,n neutrinosphere temperature (MeV)
!  m_sphere       : k,n neutrinosphere enclosed mass (g)
!  rsphere_mean   : mean n-neutrinosphere radius (cm)
!  dsphere_mean   : mean n-neutrinosphere density (g cm^{-3})
!  tsphere_mean   : mean n-neutrinosphere temperature (MeV)
!  msphere_mean   : mean n-neutrinosphere enclosed mass (g)
!  esphere_mean   : mean n-neutrinosphere energy at the mean nu_sphere (MeV)
!  lum            : n-luminosity at radius r_lum (foes)
!  lum_r          : n-luminosity at radius r_lum (foes)
!  lum_rho        : n-luminosity at density d_lum(foes)
!  inv_fluxfact   : j,n neutrino inverse flux factor
!  stat_pinch_r   : static spectral pinch factor at radius r_pinch
!  trns_pinch_r   : transport spectral pinch factor at radius r_pinch
!  stat_pinch_d   : static spectral pinch factor at density d_pinch
!  trns_pinch_d   : transport spectral pinch factor at density d_pinch
!  vel_x          : x-component of velocity (cm s^{-1})
!  vel_y          : y-component of velocity (cm s^{-1})
!  vel_z          : z-component of velocity (cm s^{-1})
!  xi             : zone-edged x-coordinate (cm)
!  xiph           : zone-edged x-coordinate (cm)
!  thetai         : zone-edged y-coordinate
!  thetaiph       : zone-centered y-coordinate
!  rho_ed         : density (g cm^{-3})
!  t_ed           : temperature (MeV)
!  ye_ed          : electron fraction
!  s_ed           : entropy
!  yl_ed          : lepton fraction
!  dm_ed          : zone mass (solar masses)
!  yenu           : electron neutrino fraction
!  yenubar        : electron antineutrino fraction
!  yxnu           : muon or tau neutrino fraction
!  yxnubar        : muon or tau antineutrino fraction
!
!    Include files:
!      kind_module, numerical_module, physcnst_module
!      cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!      nucbrn_module, nu_dist_module, nu_energy_grid_module,
!      t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : kmev, msolar

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : r_lum, d_lum, r_e_rms, d_e_rms, i_MDedit, dt_MDedit1, &
& dt_MDedit2, r_pinch, d_pinch
USE eos_snc_x_module, ONLY : xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, nse, aesv
USE mdl_cnfg_module, ONLY : rho, t, ye, dr, r, u, dmrst, rstmss, jr_min, jr_max
USE nucbrn_module, ONLY : xn_n =>xn, be_nuc_rep_n =>be_nuc_rep, a_nuc_rep_n =>a_nuc_rep, &
& z_nuc_rep_n =>z_nuc_rep, fescrn, fascrn, uburn_n =>uburn, nse_n =>nse
USE nu_dist_module, ONLY : psi0, psi1, vol, rjmh, j_sphere, r_spherep  => r_sphere, &
& d_spherep  => d_sphere, t_spherep  => t_sphere, m_spherep  => m_sphere
USE nu_energy_grid_module, ONLY : nnugp, unui, unubi
USE t_cntrl_module, ONLY: dtime_hydro, dtnph, time, t_bounce

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

LOGICAL, INTENT(in)                             :: first         ! initial model print flag

INTEGER, INTENT(in)                             :: imin          ! inner x-array index
INTEGER, INTENT(in)                             :: imax          ! outer x-array index
INTEGER, INTENT(in)                             :: nx            ! x-array extent

INTEGER, INTENT(in)                             :: js            ! inner y-array index
INTEGER, INTENT(in)                             :: je            ! outer y-array index
INTEGER, INTENT(in)                             :: ny            ! y-array extent

INTEGER, INTENT(in)                             :: i_ray_dim     ! number of rays assigned to a processor
INTEGER, INTENT(in)                             :: i_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)                             :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)                             :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)                             :: nnc           ! composition array extent

INTEGER, INTENT(in)                             :: cycle_number  ! cycle number

INTEGER, INTENT(in), DIMENSION(nx,i_ray_dim)    :: nsep          ! nuclear sttistical equilibrium flag

REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: rhop         ! density (cm^{-3})
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: tp           ! temperature (K)
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: yep          ! electron fraction
REAL(KIND =double), INTENT(in), DIMENSION(nx)                   :: rhobarp      ! mean density (cm^{-3})
REAL(KIND =double), INTENT(in), DIMENSION(nx+1)                 :: rp           ! zone-edged radial coordinate (cm)
REAL(KIND =double), INTENT(in), DIMENSION(nx)                   :: rpph         ! zone-centered radial coordinate (cm)
REAL(KIND =double), INTENT(in), DIMENSION(ny+1)                 :: thetap       ! zone-edged y (angular) coordinate
REAL(KIND =double), INTENT(in), DIMENSION(ny)                   :: thetapph     ! zone-centered y (angular) coordinate
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: up           ! x-component of zone-centered velocity (cm s^{-1})
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: vp           ! y-component of zone-centered velocity (cm s^{-1})
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: wp           ! z-component of zone-centered velocity (cm s^{-1})
REAL(KIND =double), INTENT(in), DIMENSION(nx,nez,nnu,i_ray_dim) :: psi0p        ! zeroth angular moment of the NDF
REAL(KIND =double), INTENT(in), DIMENSION(nx,nez,nnu,i_ray_dim) :: psi1p        ! first angular moment of the NDF

REAL(KIND =double), INTENT(in), DIMENSION(nx,nnc,i_ray_dim)     :: xnp          ! Composition mass fractions
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: be_nuc_repp  ! binding energy of mean heavy nucleus (MeV)
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: a_nuc_repp   ! mass number of mean heavy nucleus
REAL(KIND =double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: z_nuc_repp   ! charge number of mean heavy nucleus

REAL(KIND =double), INTENT(in)                                  :: dtime        ! time step used at current cycle
REAL(KIND =double), INTENT(in)                                  :: time_elapsed ! elapsed time

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND =double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim)    :: e_rms_stat   ! sqrt[ SUM psi0 * w5dw/SUM w3ww ] (MeV) (RMS static neutrino energy)
REAL(KIND =double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim)    :: e_rms_trns   ! sqrt[ SUM psi1 * w5dw/SUM w3ww ] (MeV) (RMS transport neutrino energy)

REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: e_rms_r_stat ! rms static neutrino energy at radius r_e_rms (MeV)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: e_rms_r_trns ! rms transport neutrino energy at radius r_e_rms (MeV)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: e_rms_d_stat ! rms static neutrino energy at density d_e_rms (MeV)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: e_rms_d_trns ! rms transport neutrino energy at density d_e_rms (MeV)

REAL(KIND =double), INTENT(out), DIMENSION(nez,nnu,i_ray_dim)   :: r_sphere     ! k,n neutrinosphere radius (cm)
REAL(KIND =double), INTENT(out), DIMENSION(nez,nnu,i_ray_dim)   :: d_sphere     ! k,n neutrinosphere density (g cm^{-3})
REAL(KIND =double), INTENT(out), DIMENSION(nez,nnu,i_ray_dim)   :: t_sphere     ! k,n neutrinosphere temperature (MeV)
REAL(KIND =double), INTENT(out), DIMENSION(nez,nnu,i_ray_dim)   :: m_sphere     ! k,n neutrinosphere enclosed mass (g)

REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: rsphere_mean ! mean n-neutrinosphere radius (cm)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: dsphere_mean ! mean n-neutrinosphere density (g cm^{-3})
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: tsphere_mean ! mean n-neutrinosphere temperature (MeV)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: msphere_mean ! mean n-neutrinosphere enclosed mass (g)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: esphere_mean ! mean n-neutrinosphere energy at the mean nu_sphere (MeV)

REAL(KIND =double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim)    :: lum          ! n-luminosity (foes)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: lum_r        ! n-luminosity at radius r_lum (foes)
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: lum_rho      ! n-luminosity at density d_lum (foes)

REAL(KIND =double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim)    :: inv_fluxfact ! inverse flux factor

REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: stat_pinch_r ! static spectral pinch factor at radius r_pinch
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: trns_pinch_r ! transport spectral pinch factor at radius r_pinch
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: stat_pinch_d ! static spectral pinch factor at density d_pinch
REAL(KIND =double), INTENT(out), DIMENSION(nnu,i_ray_dim)       :: trns_pinch_d ! transport spectral pinch factor at density d_pinch

REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: vel_x        ! x-component of velocity (cm s^{-1})
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: vel_y        ! y-component of velocity (cm s^{-1})
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: vel_z        ! z-component of velocity (cm s^{-1})
REAL(KIND =double), INTENT(out), DIMENSION(nx+1)                :: xi           ! zone-edged x-coordinate (cm)
REAL(KIND =double), INTENT(out), DIMENSION(nx)                  :: xiph         ! zone-centered x-coordinate (cm)
REAL(KIND =double), INTENT(out), DIMENSION(ny+1)                :: thetai       ! zone-edged y-coordinate
REAL(KIND =double), INTENT(out), DIMENSION(ny)                  :: thetaiph     ! zone-centered y-coordinate
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: rho_ed       ! density (g cm^{-3})
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: t_ed         ! temperature (MeV)
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: ye_ed        ! electron fraction
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: s_ed         ! entropy
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: yl_ed        ! lepton fraction
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: yenu         ! electron neutrino fraction
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: yenubar      ! electron antineutrino fraction
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: yxnu         ! muon or tau neutrino fraction
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: yxnubar      ! muon or tau antineutrino fraction
REAL(KIND =double), INTENT(out), DIMENSION(nx,i_ray_dim)        :: dm_ed        ! zone mass (solar masses)

REAL(KIND =double), INTENT(out), DIMENSION(nez)                 :: e_nu_center  ! zone-centered neutrino energy
REAL(KIND =double), INTENT(out), DIMENSION(nez+1)               :: e_nu_edge    ! zone-edged neutrino energy

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                                         :: l_MDedit

INTEGER                                         :: jr_maxp       ! jr_max + 1
INTEGER                                         :: i             ! do index
INTEGER                                         :: j             ! radial zone index
INTEGER                                         :: k             ! neutrino energy index
INTEGER                                         :: n             ! neutrino flavor index
INTEGER                                         :: nc            ! compoisition index
INTEGER                                         :: itime         ! used to determine when to generate a HDF edit
INTEGER                                         :: itimeprev     ! used to determine when to generate a HDF edit
INTEGER                                         :: editMD_in     ! return code

REAL(KIND =double)                               :: t_tb          ! time from bounce
REAL(KIND =double)                               :: tmult         ! used to determine when to generate a HDF edit

REAL(KIND =double), DIMENSION(nx)                :: yl            ! lepton fraction
REAL(KIND =double), DIMENSION(nx)                :: yenu_ray      ! electron neutrino fraction
REAL(KIND =double), DIMENSION(nx)                :: yenubar_ray   ! electron antineutrino fraction
REAL(KIND =double), DIMENSION(nx)                :: yxnu_ray      ! muon or tau neutrino fraction
REAL(KIND =double), DIMENSION(nx)                :: yxnubar_ray   ! muon or tau antineutrino fraction
REAL(KIND =double), DIMENSION(nx)                :: rhobar        ! mean density (MGFLD indexed( (cm^{-3})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!        \\\\\ EXAMINE CRITERIA FOR GENERATING AN HDF EDIT /////
!
!-----------------------------------------------------------------------

!........Initialize.....................................................

dtnph                     = dtime
time                      = time_elapsed
editMD_in                 = 4
dm_ed                     = zero

!........Return if i_editMG  = 0.........................................

IF ( i_MDedit == 0 ) RETURN 

!........Initialize.....................................................

l_MDedit                  = .false.

IF ( t_bounce <= zero ) THEN

!........Generate an HDF edit every multiple of dt_MDedit1 ms............
!........up to bounce bounce............................................

  tmult                   = 1.d+3/dt_MDedit1
  itime                   = int( time * tmult )
  itimeprev               = int( ( time - dtnph ) * tmult )
  IF ( itime /= itimeprev ) l_MDedit  = .true.

ELSE

  t_tb                    = time - t_bounce

!........Generate an HDF edit at bounce..................................

  IF ( t_tb - dtnph <= zero  .and.  t_tb > zero ) l_MDedit  = .true.

!........Generate an HDF edit 1 ms after bounce..........................

  IF ( t_tb - dtnph <= 1.d-3  .and.  t_tb > 1.d-3 ) l_MDedit  = .true.

!........Generate an HDF edit 10 ms after bounce.........................

  IF ( t_tb - dtnph <= 1.d-2  .and.  t_tb > 1.d-2 ) l_MDedit  = .true.

!........Generate an HDF edit every multiple of dt_MDedit2 ms............
!....... after bounce...................................................

  tmult                   = 1.d+3/dt_MDedit2
  itime                   = int( t_tb * tmult )
  itimeprev               = int( ( t_tb - dtnph ) * tmult )
  IF ( itime /= itimeprev ) l_MDedit  = .true.

END IF

!........Return if l_MDedit  = .false...................................

IF ( .not. l_MDedit ) RETURN

!-----------------------------------------------------------------------
!
!            \\\\\ LOAD VARIABLES INTO MGFLD ARRAYS /////
!
!-----------------------------------------------------------------------

!........Transfer integer scalar parameters.............................

jr_min                    = imin + 1
jr_max                    = imax + 1
jr_maxp                   = jr_max + 1
ncycle                    = cycle_number

!........Transfer zone-centered independent variables to mgfld arrays...

DO i = imin,imax
  j                       = i - imin + 2
  rho (j)                 = rhop (i,i_ray)
  t   (j)                 = tp   (i,i_ray)
  ye  (j)                 = yep  (i,i_ray)
  u   (j)                 = up   (i,i_ray)
  rjmh(j)                 = rpph(i)
  rhobar(j)               = rhobarp(i)
END DO

DO n  = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k  = 1,nnugp(n)
    DO i  = imin,imax
      psi0(i-imin+2,k,n)  = psi0p(i,k,n,i_ray)
      psi1(i-imin+1,k,n)  = psi1p(i,k,n,i_ray)
    END DO
    psi1(imax-imin+2,k,n) = psi1p(imax+1,k,n,i_ray)
  END DO
END DO

!........Transfer zone-edgeed independent variables to mgfld arrays.....

DO i = imin,imax+1
  r(i-imin+1)             = rp(i)
END DO

!........Transfer zone-edgeed composition variables.....................

DO i = imin,imax
  j                       = i-imin+2
  nse(j,i_ray)            = nsep   (i,i_ray)
  be_nuc_rep(j)           = be_nuc_repp(i,i_ray)
  a_nuc_rep (j)           = a_nuc_repp (i,i_ray)
  z_nuc_rep (j)           = z_nuc_repp (i,i_ray)
END DO

DO nc = 1,nnc
  DO i = imin,imax+1
    xn(i-imin+2,nc)       = xnp(i,nc,i_ray)
  END DO
END DO

DO i = imin,imax
  j                       = i-imin+2
  nse_n   (j)             = nsep   (i,i_ray)
  be_nuc_rep_n(j)         = be_nuc_repp(i,i_ray)
  a_nuc_rep_n (j)         = a_nuc_repp(i,i_ray)
  z_nuc_rep_n (j)         = z_nuc_repp(i,i_ray)
END DO

DO nc = 1,nnc
  DO i = imin,imax
    xn_n(i-imin+2,nc)     = xnp(i,nc,i_ray)
  END DO
END DO

!........Derived zone-centered and zone-edged variables.................

!........radial zone thickness

DO j = jr_min,jr_max
  dr(j)                   = r(j) - r(j-1)
END DO

!........Newtonian rest masses

IF ( r(1) == zero ) THEN
  dmrst(1)                = zero
  dm_ed(1,i_ray)          = zero
ELSE
  dmrst(1)                = frpith * r(1)**3 * rhobar(2)
  dm_ed(1,i_ray)          = dmrst(1)/msolar
END IF

rstmss(1)                 = dmrst(1)

DO j = jr_min,jr_max
  vol(j)                  = frpith * ( r(j) * ( r(j) + r(j-1) ) + r(j-1) * r(j-1) ) * dr(j)
  dmrst(j)                = vol(j) * rhobar(j)
  dm_ed(j,i_ray)          = dmrst(j)/msolar
  rstmss(j)               = rstmss(j-1) + dmrst(j)
END DO

!........Update EOS variables...........................................

CALL eqstz_x( jr_min, jr_maxp, rho, t, ye, i_ray )
CALL gammaz_x( jr_min, jr_max, rho, t, i_ray )

!-----------------------------------------------------------------------
!
!              \\\\\ COMPUTE QUANTITIES TO EDIT /////
!
!-----------------------------------------------------------------------

!........Model configuration

CALL yl_cal ( jr_min, jr_max, rho, ye, psi0, yenu_ray, yenubar_ray, yxnu_ray, &
& yxnubar_ray, yl, nx, nez, nnu )

DO i = imin,imax
  j                       = i-imin+2
  vel_x(i,i_ray)          = up(i,i_ray)
  vel_y(i,i_ray)          = vp(i,i_ray)
  vel_z(i,i_ray)          = wp(i,i_ray)
  xi(i)                   = rp(i)
  xiph(i)                 = rpph(i)
  rho_ed(i,i_ray)         = rhop(i,i_ray)
  t_ed(i,i_ray)           = tp(i,i_ray) * kmev
  ye_ed(i,i_ray)          = yep(i,i_ray)
  s_ed(i,i_ray)           = aesv(j,3,i_ray)
  yl_ed(i,i_ray)          = yl(j)
  yenu(i,i_ray)           = yenu_ray(j)
  yenubar(i,i_ray)        = yenubar_ray(j)
  yxnu(i,i_ray)           = yxnu_ray(j)
  yxnubar(i,i_ray)        = yxnubar_ray(j)
END DO

xi(imax+1)                = rp(imax+1)

DO j = js,je
  thetai(j)               = thetap(j)
  thetaiph(j)             = thetapph(j)
END DO

thetai(je+1)              = thetap(je+1)

!........Compute the neutrinospheres

CALL nu_sphere( jr_min, jr_maxp, i_ray, i_ray_dim, r, rho, t, rstmss, nx, nez, nnu, &
& j_sphere, r_sphere, d_sphere, t_sphere, m_sphere )

DO n = 1,nnu
  DO k = 1,nez
    r_spherep(k,n,i_ray)  = r_sphere(k,n,i_ray)
    d_spherep(k,n,i_ray)  = d_sphere(k,n,i_ray)
    t_spherep(k,n,i_ray)  = t_sphere(k,n,i_ray)
    m_spherep(k,n,i_ray)  = m_sphere(k,n,i_ray)
  END DO
END DO

!........Compute psi1 if first

IF ( first ) CALL psi1_cal( jr_min, jr_max, i_ray, i_ray_dim, rho, t, ye, &
& r, u, rstmss, psi0, psi1, nx, nez, nnu, 1 )

!........Compute neutrino rms energies

CALL e_rms( jr_min, jr_maxp, i_ray, i_ray_dim, nx, nnu, e_rms_stat, e_rms_trns )

!........Compute the neutrino rms energies at r_e_rms and d_e_rms

CALL e_rms_surface( jr_min, jr_max, i_ray, i_ray_dim, nnu, r_e_rms, d_e_rms, &
& e_rms_r_stat, e_rms_r_trns, e_rms_d_stat, e_rms_d_trns )

!........Compute the mean neutrinospheres

CALL nu_sphere_mean( jr_min, jr_maxp, i_ray, i_ray_dim, e_rms_trns, j_sphere, &
& r_sphere, d_sphere, t_sphere, m_sphere, nx, nez, nnu, rsphere_mean, dsphere_mean, &
& tsphere_mean, msphere_mean, esphere_mean )

!........Compute the neutrino luminosities

CALL luminosity( jr_min, jr_max, i_ray, i_ray_dim, nx, nnu, lum )

!........Compute the neutrino luminosities at r_lum and d_lum

CALL lum_surface( jr_min, jr_max, i_ray, i_ray_dim, nnu, r_lum, d_lum, lum_r, &
& lum_rho )

!........Compute neutrino inverse flux factors

CALL inv_fluxfactor( jr_min, jr_max, i_ray, i_ray_dim, nx, nnu, inv_fluxfact )

!........Compute the pinching parameters at r_pibch and d_pinch

CALL pinch_surface( jr_min, jr_max, i_ray, i_ray_dim, nnu, r_pinch, d_pinch, &
& stat_pinch_r, trns_pinch_r, stat_pinch_d, trns_pinch_d )

!........Neutrino zone energies

DO k = 1,nez
  e_nu_center(k)          = unui(k)
  e_nu_edge(k)            = unubi(k)
END DO
e_nu_edge(nez+1)          = unubi(nez+1)

editMD_in  = 0
END FUNCTION editMD_in
