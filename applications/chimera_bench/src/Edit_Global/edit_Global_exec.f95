SUBROUTINE edit_Global_exec( imin, imax, nx, nnu, time, t_tb, x_c_in, &
& mass_a_ray, rho_min, rho_max, rho_mean, sigma_rho, t_min, t_max, t_mean, &
& sigma_t, u_min, u_max, u_mean, sigma_u, v_min, v_max, v_mean, sigma_v, &
& w_min, w_max, w_mean, sigma_w, s_min, s_max, s_mean, sigma_s, ye_min, &
& ye_max, ye_mean, sigma_ye, dudt_nuc_min, dudt_nuc_max, dudt_nuc_mean, &
& dudt_nuc_zone, dudt_nu_min, dudt_nu_max, dudt_nu_mean, dudt_nu_zone, &
& mach_min, mach_max, mach_mean, sigma_mach, lum_min, lum_max, lum_mean, &
& sigma_lum, e_rms_min, e_rms_max, e_rms_mean, sigma_e_rms, grav_pot_a_ray, &
& ke_a_ray, e_no_bind_a_ray, e_bind_a_ray, e_bind_fnl_a_ray, xn_a_ray, &
& xn_mean, a_mean, nnc, nse_min, nse_max, c_shock, p_SASI, l_global_ned, &
& l_global_ted )
!-----------------------------------------------------------------------
!
!    File:         edit_Global_exec
!    Module:       edit_Global_exec
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/22/05
!
!    Purpose:
!      To perform 2D edits of the problem at preselected intervals
!
!    Subprograms called:
!  edit2Du          : performs an edit of the x-velocities
!  edit2Dv          : performs an edit of the y-velocities
!
!    Input arguments:
!  imin             : minimum x-array index for the edit
!  imax             : maximum x-array index for the edit
!  nx               : x-array extent
!  nnu              : neutrino flavor array extent
!  time             : elapsed time
!  t_tb             : time from core bounce
!  x_c_in           : radial midpoint of zone (cm)
!  mass_a_ray       : mass/angular ray (g)
!  rho_min          : minimum density/(radial zone) (g cm^{-3})
!  rho_max          : maximum density/(radial zone) (g cm^{-3})
!  rho_mean         : mass averaged density/radius (g cm^{-3})
!  sigma_rho        : RMS density/(radial zone) (g cm^{-3})
!  t_min            : minimum temperature/(radial zone) (MeV)
!  t_max            : maximum temperature/(radial zone) (MeV)
!  t_mean           : mass averaged temperature/radius (MeV)
!  sigma_t          : RMS temperature/(radial zone) (MeV)
!  u_min            : minimum radial velocity/(radial zone) (cm s^{-1})
!  u_max            : maximum radial velocity/(radial zone) (cm s^{-1})
!  u_mean           : mass averaged radial velocity/radius (cm s^{-1})
!  sigma_u          : RMS radial velocity/(radial zone) (MeV)
!  v_min            : minimum y-velocity/(radial zone) (cm s^{-1})
!  v_max            : maximum y-velocity/(radial zone) (cm s^{-1})
!  v_mean           : mass averaged y-velocity/radius (cm s^{-1})
!  sigma_v          : RMS y-velocity velocity/(radial zone) (MeV)
!  w_min            : minimum z-velocity/(radial zone) (cm s^{-1})
!  w_max            : maximum z-velocity/(radial zone) (cm s^{-1})
!  w_mean           : mass averaged z-velocity/radius (cm s^{-1})
!  sigma_w          : RMS z-velocity velocity/(radial zone) (MeV)
!  s_min            : minimum entropy(radial zone)
!  s_max            : maximum entropy/(radial zone)
!  s_mean           : mass averaged entropy/radius
!  sigma_s          : RMS entropy velocity/(radial zone) (MeV)
!  ye_min           : minimum electron fraction(radial zone)
!  ye_max           : maximum electron fraction/(radial zone)
!  ye_mean          : mass averaged electron fraction/radius
!  sigma_ye         : RMS electron fraction velocity/(radial zone) (MeV)
!  dudt_nuc_min     : minimum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
!  dudt_nuc_max     : maximum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
!  dudt_nuc_mean    : mass averaged energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
!  dudt_nuc_zone    : energy generation rate by nuclear reactions (B s^{-1} radial-zone^{-1})
!  dudt_nu_min      : minimum neutrino energy deposition rate (ergs g^{-1} s^{1})
!  dudt_nu_max      : maximum neutrino energy deposition rate (ergs g^{-1} s^{1})
!  dudt_nu_mean     : mass averaged neutrino energy deposition rate (ergs g^{-1} s^{1})
!  dudt_nu_zone     : neutrino energy deposition rate  (B s^{-1} radial-zone^{-1})
!  mach_min         : minimum mach number/radius
!  mach_max         : maximum mach number/radius
!  mach_mean        : mass averaged mach number/radius
!  sigma_mach       : mass averaged mach number/radius
!  lum_mean         : mass averaged n-neutrino luminosity (foes s^{1})
!  lum_min          : minimum n-neutrino luminosity (foes s^{1}
!  lum_max          : maximum  n-neutrino luminosity (foes s^{1})
!  sigma_lum        : RMS neutrino n-neutrino luminosity (foes s^{1})
!  e_rms_mean       : mass averaged neutrino n-neutrino rms energy (MeV)
!  e_rms_min        : minimum neutrino n-neutrino rms energy (MeV)
!  e_rms_max        : maximum neutrino n-neutrino rms energy (MeV)
!  sigma_e_rms      : RMS neutrino n-neutrino rms energy (MeV))
!  grav_pot_a_ray   : angular ray gravitational potential (ergs)
!  ke_a_ray         : angular ray kinetic energy (ergs)
!  e_no_bind_a_ray  : angular ray internal minus binding energy (ergs)
!  e_bind_a_ray     : angular ray current binding energy (ergs)
!  xn_a_ray         : composition nasses (solar masses)
!  xn_mean          : ngular average of the composition mass fractions
!  a_mean           : angular average of the mean nuclear mass number along an angular ray
!  e_bind_fnl_a_ray : angular ray final binding energy (ergs)
!  nnc              : composition array extent
!  nse_min          : minimum value of NSE-nonNSE boundary
!  nse_max          : maximum value of NSE-nonNSE boundary
!  c_shock          : shock location
!  p_SASI           : power in the convection or SASI as a function of the Legendre mode
!  l_global_ned     : cycle criterion edit flag
!  l_global_ted     : time criterion edit flag
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  edit_module
!      
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE edit_module, ONLY: nned_global, nted_global, data_path, nprint, nlog

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER                               :: nleg = 10        ! highest Legendre polynomial to be used 

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

LOGICAL, INTENT(in)                              :: l_global_ned     ! cycle criterion edit flag
LOGICAL, INTENT(in)                              :: l_global_ted     ! time criterion edit flag

INTEGER, INTENT(in)                              :: imin             ! minimum x-array index for the edit
INTEGER, INTENT(in)                              :: imax             ! maximum x-array index for the edit
INTEGER, INTENT(in)                              :: nx               ! x-array extent

INTEGER, INTENT(in)                              :: nnu              ! neutrino flavor array extent
INTEGER, INTENT(in)                              :: nnc              ! composition array extent
INTEGER, INTENT(in)                              :: nse_min          ! minimum value of NSE-nonNSE boundary
INTEGER, INTENT(in)                              :: nse_max          ! maximum value of NSE-nonNSE boundary

CHARACTER (len=2), INTENT(in), DIMENSION(nx)     :: c_shock          ! shock location

REAL(KIND=double), INTENT(in)                    :: time             ! elapsed time
REAL(KIND=double), INTENT(in)                    :: t_tb             ! time from core bounce
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: x_c_in           ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mass_a_ray       ! mass/(angular ray) (g)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: rho_mean         ! mass averaged density/radius (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: rho_min          ! minimum density/(radial zone) (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: rho_max          ! maximum density/(radial zone) (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_rho        ! RMS density/(radial zone) (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: t_mean           ! mass averaged temperature/radius (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: t_min            ! minimum temperature/(radial zone) (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: t_max            ! maximum temperature/(radial zone) (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_t          ! RMS temperature/(radial zone) (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: u_mean           ! mass averaged radial velocity/radius (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: u_min            ! minimum radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: u_max            ! maximum radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_u          ! RMS radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: v_mean           ! mass averaged y-velocity/radius (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: v_min            ! minimum y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: v_max            ! maximum y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_v          ! RMS y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: w_mean           ! mass averaged z-velocity/radius (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: w_min            ! minimum z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: w_max            ! maximum z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_w          ! RMS z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: s_mean           ! mass averaged entropy/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: s_min            ! minimum entropy/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: s_max            ! maximum entropy/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_s          ! RMS entropy/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ye_mean          ! mass averaged electron fraction/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ye_min           ! minimum electron fraction/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ye_max           ! maximum electron fraction/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_ye         ! RMS electron fraction/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_mean    ! mass averaged energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_min     ! minimum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_max     ! maximum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_zone    ! energy generation rate by nuclear reactions (B s^{-1} radial-zone^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_mean     ! mass averaged neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_min      ! minimum neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_max      ! maximum neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_zone     ! neutrino energy deposition rate  (B s^{-1} radial-zone^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mach_min         ! minimum mach number/radius
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mach_max         ! maximum mach number/radius
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mach_mean        ! mass averaged mach number/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_mach       ! RMS mach number/radius
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: lum_mean         ! mass averaged n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: lum_min          ! minimum n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: lum_max          ! maximum  n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: sigma_lum        ! RMS neutrino n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: e_rms_mean       ! mass averaged neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: e_rms_min        ! minimum neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: e_rms_max        ! maximum neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: sigma_e_rms      ! RMS neutrino n-neutrino rms energy (MeV))
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: grav_pot_a_ray   ! angular ray gravitational potential (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ke_a_ray         ! angular ray kinetic energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: e_no_bind_a_ray  ! angular ray internal minus binding energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: e_bind_a_ray     ! angular ray current binding energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: e_bind_fnl_a_ray ! angular ray final binding energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc) :: xn_a_ray         ! composition nasses (solar masses)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc) :: xn_mean          ! angular average of the composition mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: a_mean           ! angular average of the mean nuclear mass number along an angular ray
REAL(KIND=double), INTENT(in), DIMENSION(nleg)   :: p_SASI           ! power in the convection or SASI as a function of the Legendre mode

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)                      :: outfile_n     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_t     ! character string containing name of model files

INTEGER                                  :: istat         ! open file flag
INTEGER                                  :: istat_n       ! open new file flag (passed to subroutine)
INTEGER, PARAMETER                       :: iunit = 20    ! unit number to print models
INTEGER                                  :: iunitp        ! unit number to print models

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 8001 format (' File ',a128,/' cannot be opened in subroutime edit_Global_exec')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ PRINT IF CYCLE CRITERIA SATISFIED /////
!
!-----------------------------------------------------------------------

!........min, max, and mass averaged data

IF ( l_global_ned ) THEN

  nned_global      = nned_global + 1

  WRITE (outfile_n,'(a27,i5.5,a2)') '/Models_Global_n/model_mean', &
&  nned_global,'.d'
  outfile_n        = TRIM(data_path)//TRIM(outfile_n)
  OPEN (UNIT=iunit, FILE=TRIM(outfile_n), RECL=512, STATUS='new', IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=iunit, FILE=TRIM(outfile_n), RECL=512, STATUS='old')

  iunitp           = iunit
  CALL edit_Global_mean( imin, imax, nx, time, t_tb, x_c_in, mass_a_ray, &
&  rho_min, rho_max, rho_mean, sigma_rho, t_min, t_max, t_mean, sigma_t, &
&  u_min, u_max, u_mean, sigma_u, v_min, v_max, v_mean, sigma_v, w_min, &
&  w_max, w_mean, sigma_w, s_min, s_max, s_mean, sigma_s, ye_min, ye_max, &
&  ye_mean, sigma_ye, mach_min, mach_max, mach_mean, sigma_mach, c_shock, &
&  iunitp )

  CALL edit_Global_energetics( imin, imax, nx, time, t_tb, x_c_in,  &
&  grav_pot_a_ray, ke_a_ray, e_no_bind_a_ray, e_bind_a_ray, e_bind_fnl_a_ray, &
&  dudt_nuc_min, dudt_nuc_max, dudt_nuc_mean, dudt_nuc_zone, dudt_nu_min, &
&  dudt_nu_max, dudt_nu_mean, dudt_nu_zone, c_shock, iunitp )

  CALL edit_Global_composition( imin, imax, nx, nnc, time, t_tb, x_c_in,  &
&  xn_a_ray, xn_mean, a_mean, nse_min, nse_max, c_shock, iunitp )

  CALL edit_Global_neutrinos( imin, imax, nx, nnu, time, t_tb, x_c_in,  &
&  nse_min, nse_max, lum_mean, lum_min, lum_max, sigma_lum, e_rms_mean, &
&  e_rms_min, e_rms_max, sigma_e_rms, c_shock, iunitp )

  CLOSE (unit=iunit,status='keep')

  WRITE (outfile_n,'(a24)') '/Plot_Files/Power_SASI.d'
  outfile_n        = TRIM(data_path)//TRIM(outfile_n)
  OPEN (UNIT=iunit, FILE=TRIM(outfile_n), RECL=512, STATUS='new', IOSTAT=istat )
  istat_n          = istat
  IF ( istat /= 0 ) OPEN (UNIT=iunit, FILE=TRIM(outfile_n), RECL=512, &
&  STATUS='old', POSITION='append', IOSTAT=istat )
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8001) outfile_n
    WRITE (nlog,8001) outfile_n
    STOP
  END IF ! istat /= 0
  
  iunitp           = iunit
  CALL power_SASI_plot( time, t_tb, p_SASI, istat_n, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_global_ned

!-----------------------------------------------------------------------
!
!             \\\\\ PRINT IF TIME CRITERIA SATISFIED /////
!
!-----------------------------------------------------------------------

IF ( l_global_ted ) THEN

  nted_global      = nted_global + 1

!........min, max, and mass averaged data

  WRITE (outfile_t,'(a27,i5.5,a2)') '/Models_Global_t/model_mean', &
&  nted_global,'.d'
  outfile_t        = TRIM(data_path)//TRIM(outfile_t)
  OPEN (UNIT=iunit, FILE=TRIM(outfile_t), RECL=512, STATUS='new', IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=iunit, FILE=TRIM(outfile_t), RECL=512, STATUS='old')

  iunitp           = iunit
  CALL edit_Global_mean( imin, imax, nx, time, t_tb, x_c_in, mass_a_ray, &
&  rho_min, rho_max, rho_mean, sigma_rho, t_min, t_max, t_mean, sigma_t, &
&  u_min, u_max, u_mean, sigma_u, v_min, v_max, v_mean, sigma_v, w_min, &
&  w_max, w_mean, sigma_w, s_min, s_max, s_mean, sigma_s, ye_min, ye_max, &
&  ye_mean, sigma_ye, mach_min, mach_max, mach_mean, sigma_mach, c_shock, &
&  iunitp )

  CALL edit_Global_energetics( imin, imax, nx, time, t_tb, x_c_in,  &
&  grav_pot_a_ray, ke_a_ray, e_no_bind_a_ray, e_bind_a_ray, e_bind_fnl_a_ray, &
&  dudt_nuc_min, dudt_nuc_max, dudt_nuc_mean, dudt_nuc_zone, dudt_nu_min, &
&  dudt_nu_max, dudt_nu_mean, dudt_nu_zone, c_shock, iunitp )

  CALL edit_Global_composition( imin, imax, nx, nnc, time, t_tb, x_c_in,  &
&  xn_a_ray, xn_mean, a_mean, nse_min, nse_max, c_shock, iunitp )

  CALL edit_Global_neutrinos( imin, imax, nx, nnu, time, t_tb, x_c_in,  &
&  nse_min, nse_max, lum_mean, lum_min, lum_max, sigma_lum, e_rms_mean, &
&  e_rms_min, e_rms_max, sigma_e_rms, c_shock, iunitp )

  CLOSE (unit=iunit,status='keep')

  WRITE (outfile_n,'(a24)') '/Plot_Files/Power_SASI.d'
  outfile_n        = TRIM(data_path)//TRIM(outfile_n)
  OPEN (UNIT=iunit, FILE=TRIM(outfile_n), RECL=512, STATUS='new', IOSTAT=istat )
  istat_n          = istat
  IF ( istat /= 0 ) OPEN (UNIT=iunit, FILE=TRIM(outfile_n), RECL=512, &
&  STATUS='old', POSITION='append', IOSTAT=istat )
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8001) outfile_n
    WRITE (nlog,8001) outfile_n
    STOP
  END IF ! istat /= 0
  
  iunitp           = iunit
!  CALL power_SASI_plot( time, t_tb, p_SASI, istat_n, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_global_ted


RETURN
END SUBROUTINE edit_Global_exec
