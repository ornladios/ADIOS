SUBROUTINE edit_MD_exec( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz,         &
& nnu, time, t_tb, x_in, y_in, z_in, rhoMD, uMD, vMD, wMD, pMD, eMD, sMD, aMD,   &
& yeMD, v_csoundMD, e_nu_MD, f_nu_MD, lumMD, e_rmsMD, r_nse, r_O1, r_xO,         &
& r_shock, r_shock_mn, r_shock_mx, rsphere_mean, dsphere_mean, tsphere_mean,     &
& msphere_mean, esphere_mean, r_gain, tau_advMD, tau_heat_nuMD, tau_heat_nucMD,  &
& n_growMD, wbvMD, ylMD, i_comp, dudt_nucMD, dudt_nuMD, grav_x_cMD, grav_y_cMD,  &
& grav_z_cMD, l_editMDu, l_editMDv, l_editMDw, l_editMDs, l_editMDd, l_editMDe,  &
& l_editMDp, l_editMDenu, l_editMDfnu, l_editMDa, l_editMDx, l_editMDye,         &
& l_editMDcm, l_editMDnu, l_editMDnc, l_editMDnl, l_editMDne, l_editMDgx,        &
& l_editMDgy, l_editMDgz, l_editMDBVw, l_editMDyl, l_MDedit )
!-----------------------------------------------------------------------
!
!    File:         edit_MD_exec
!    Module:       edit_MD_exec
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/22/05
!
!    Purpose:
!      To perform MD edits of the problem at preselected intervals
!
!    Subprograms called:
!  editMDu        : performs an edit of the x-velocities
!  editMDv        : performs an edit of the y-velocities
!
!    Input arguments:
!  imin           : minimum x-array index for the edit
!  imax           : maximum x-array index for the edit
!  nx             : x-array extent
!  jmin           : minimum y-array index for the edit
!  jmax           : maximum y-array index for the edit
!  ny             : y_array extent
!  kmin           : minimum z-array index for the edit
!  kmax           : maximum z-array index for the edit
!  nz             : z dimension
!  nnu            : neutrino flavor array extent
!  time           : elapsed time
!  t_tb           : time from core bounce
!  x_in           : radial midpoint of zone [cm]
!  y_in           : y (angular) midpoint of zone
!  z_in           : z (azimuthal) midpoint of zone
!  rhoMD          : density [g s^{-3}]
!  uMD            : x (radial) velocity of zone [cm s^{-1}]
!  vMD            : y (angular) velocity of zone [cm s^{-1}]
!  wMD            : z (azimuthal) velocity of zone [cm s^{-1}]
!  pMD            : pressure of zone [ergs cm^{-3}]
!  eMD            : internal energy of zone [ergs g^{-1}]
!  sMD            : entropy of zone
!  aMD            : mean nuclear mass number
!  yeMD           : electron fraction of zone
!  v_csoundMD     : mach number of zone
!  e_nu_MD        : neutrino energy density [ergs g^{-1}]
!  f_nu_MD        : neutrino energy flux [ergs cm^{-2} s^{-1}]
!  lumMD          : neutrino luminosity [foes]
!  e_rmsMD        : neutrino rms energy ([MeV])
!  r_nse          : radius of NSE-nonNSE boundary
!  r_O1           : radius of x(16O)=0.1 boundary
!  r_xO           : radius of x(16O)=0.5 boundary
!  r_shock        : radius of shock maximum [cm]
!  r_shock_mn     : minimum estimateed shock radius [cm]
!  r_shock_mx     : maximum estimateed shock radius [cm]
!  rsphere_mean   : mean neutrinosphere radius [cm]
!  dsphere_mean   : mean neutrinosphere density [g cm^{-3)]
!  tsphere_mean   : mean neutrinosphere temperature (K)
!  msphere_mean   : mean neutrinosphere enclosed mass (g)
!  esphere_mean   : mean neutrinosphere rms neutrino energy ([MeV])
!  r_gain         : gain radius
!  tau_advMD      : advection time scale [s]
!  tau_heat_nuMD  : neutrino heating time scale [s]
!  tau_heat_nucMD : nuclear heating time scale [s] 
!  n_growMD       : number of convective e-foldings in fluid from shock to gain radius
!  wBVMD          : Brunt Vaisala frequency (positive if unstable) [s^{-1}]
!  ylMD           : lepton fraction
!  i_comp         : integer denoting the dominant element
!  dudt_nucMD     : energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
!  dudt_nuMD      : energy deposition rate by neutrinos (ergs g^{-1} s^{1})
!  grav_x_cMD     : zone-centered x-component of gravitational acceleration [cm s^{-2} g^{-1}]
!  grav_y_cMD     : zone-centered y-component of gravitational acceleration [cm s^{-2} g^{-1}]
!  grav_z_cMD     : zone-centered z-component of gravitational acceleration [cm s^{-2} g^{-1}]
!  l_editMDu      : editMDu flag
!  l_editMDv      : editMDv flag
!  l_editMDw      : editMDw flag
!  l_editMDs      : editMDs flag
!  l_editMDd      : editMDd flag
!  l_editMDe      : editMDe flag
!  l_editMDo      : editMDp flag
!  l_editMDenu    : editMDenu flag
!  l_editMDfnu    : editMDfnu flag
!  l_editMDa      : editMDa flag
!  l_editMDx      : editMDx flag
!  l_editMDye     : editMDye flag
!  l_editMDcm     : editMDcm flag
!  l_editMDnu     : editMDnu flag
!  l_editMDnc     : editMDnc flag
!  l_editMDnl     : editMDnl flag
!  l_editMDne     : editMDne flag
!  l_editMDgx     : editMDgx flag
!  l_editMDgy     : editMDgy flag
!  l_editMDgz     : editMDgz flag
!  l_editMDBVw    : editMDBVw flag
!  l_editMDyl     : editMDyl flag
!  l_MDedit       : time flag
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

USE edit_module, ONLY: n_editMDu, n_editMDv, n_editMDw, n_editMDs, n_editMDd, &
& n_editMDe, n_editMDp, n_editMDenu, n_editMDfnu, n_editMDa, n_editMDx,       &
& n_editMDye, n_editMDcm, n_editMDnu, n_editMDnc, n_editMDnl, n_editMDne,     &
& n_editMDgx, n_editMDgy, n_editMDgz, n_editMDBVw, n_editMDyl, n_editMD,      &
& data_path, nlog

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

LOGICAL, INTENT(in)                                    :: l_editMDu      ! editMDu flag
LOGICAL, INTENT(in)                                    :: l_editMDv      ! editMDv flag
LOGICAL, INTENT(in)                                    :: l_editMDw      ! editMDw flag
LOGICAL, INTENT(in)                                    :: l_editMDs      ! editMDs flag
LOGICAL, INTENT(in)                                    :: l_editMDd      ! editMDd flag
LOGICAL, INTENT(in)                                    :: l_editMDe      ! editMDe flag
LOGICAL, INTENT(in)                                    :: l_editMDp      ! editMDp flag
LOGICAL, INTENT(in)                                    :: l_editMDenu    ! editMDenu flag
LOGICAL, INTENT(in)                                    :: l_editMDfnu    ! editMDfnu flag
LOGICAL, INTENT(in)                                    :: l_editMDa      ! editMDa flag
LOGICAL, INTENT(in)                                    :: l_editMDx      ! editMDx flag
LOGICAL, INTENT(in)                                    :: l_editMDye     ! editMDye flag
LOGICAL, INTENT(in)                                    :: l_editMDcm     ! editMDcm flag
LOGICAL, INTENT(in)                                    :: l_editMDnu     ! editMDnu flag
LOGICAL, INTENT(in)                                    :: l_editMDnc     ! editMDnc flag
LOGICAL, INTENT(in)                                    :: l_editMDnl     ! editMDnl flag
LOGICAL, INTENT(in)                                    :: l_editMDne     ! editMDne flag
LOGICAL, INTENT(in)                                    :: l_editMDgx     ! editMDgx flag
LOGICAL, INTENT(in)                                    :: l_editMDgy     ! editMDgy flag
LOGICAL, INTENT(in)                                    :: l_editMDgz     ! editMDgz flag
LOGICAL, INTENT(in)                                    :: l_editMDBVw    ! editMDBVw flag
LOGICAL, INTENT(in)                                    :: l_editMDyl     ! editMDyl flag
LOGICAL, INTENT(in)                                    :: l_MDedit       ! time flag

INTEGER, INTENT(in)                                    :: imin           ! minimum x-array index for the edit
INTEGER, INTENT(in)                                    :: imax           ! maximum x-array index for the edit
INTEGER, INTENT(in)                                    :: nx             ! x-array extent

INTEGER, INTENT(in)                                    :: jmin           ! minimum y-array index for the edit
INTEGER, INTENT(in)                                    :: jmax           ! maximum y-array index for the edit
INTEGER, INTENT(in)                                    :: ny             ! y_array extent

INTEGER, INTENT(in)                                    :: kmin           ! minimum z-array index for the edit
INTEGER, INTENT(in)                                    :: kmax           ! maximum z-array index for the edit
INTEGER, INTENT(in)                                    :: nz             ! z_array extent

INTEGER, INTENT(in)                                    :: nnu            ! neutrino flavor array extent

INTEGER, INTENT(in), DIMENSION(nx,ny)                  :: i_comp         ! integer denoting the dominant element

REAL(KIND=double), INTENT(in)                          :: time           ! elapsed time
REAL(KIND=double), INTENT(in)                          :: t_tb           ! time from core bounce
REAL(KIND=double), INTENT(in), DIMENSION(nx)           :: x_in           ! radial midpoint of zone [cm]
REAL(KIND=double), INTENT(in), DIMENSION(ny)           :: y_in           ! y (angular) midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nz)           :: z_in           ! z (azimuthal) midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: rhoMD          ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: uMD            ! x (radial) velocity of zone [cm s^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: vMD            ! y (angular) velocity of zone [cm s^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: wMD            ! z (azimuthal) velocity of zone [cm s^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: pMD            ! pressure of zone [ergs cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: eMD            ! internal energy of zone [ergs g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: sMD            ! entropy of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: aMD            ! mean nuclear mass number
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: yeMD           ! electron fraction of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: v_csoundMD     ! mach number of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,ny,nz) :: e_nu_MD        ! neutrino energy density [ergs g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,ny,nz) :: f_nu_MD        ! neutrino energy flux [ergs cm^{-2} s^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,ny,nz) :: lumMD          ! neutrino luminosity [foes]
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,ny,nz) :: e_rmsMD        ! SQRT( SUM psi0 * w5dw/SUM w3ww ) ([MeV])

REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: r_nse          ! radius of NSE-nonNSE boundary
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: r_O1           ! radius of x(16O)=0.1 boundary
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: r_xO           ! radius of x(16O)=0.5 boundary
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: r_shock        ! radius of shock maximum
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: r_shock_mn     ! minimum estimateed shock radius
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: r_shock_mx     ! maximum estimateed shock radius
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)    :: rsphere_mean   ! mean neutrinosphere radius
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)    :: dsphere_mean   ! mean neutrinosphere density
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)    :: tsphere_mean   ! mean neutrinosphere temperature
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)    :: msphere_mean   ! mean neutrinosphere enclosed mass
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)    :: esphere_mean   ! mean neutrinosphere energy
REAL(KIND=double), INTENT(in), DIMENSION(nnu+1,ny,nz)  :: r_gain         ! gain radius
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: tau_advMD      ! advection time scale [s]
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: tau_heat_nuMD  ! neutrino heating time scale [s]
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: tau_heat_nucMD ! nuclear heating time scale [s] 
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)        :: n_growMD       ! number of convective e-foldings in fluid from shock to gain radius
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: wBVMD          ! Brunt Vaisala frequency (positive if unstable) [s^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: ylMD           ! lepton fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: dudt_nucMD     ! energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: dudt_nuMD      ! energy deposition rate by neutrinos (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: grav_x_cMD     ! zone-centered x-component of gravitational acceleration [cm s^{-2} g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: grav_y_cMD     ! zone-centered y-component of gravitational acceleration [cm s^{-2} g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz)     :: grav_z_cMD     ! zone-centered z-component of gravitational acceleration [cm s^{-2} g^{-1}]

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)                      :: outfile_u     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_v     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_w     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_p     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_e     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_s     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_d     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_a     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_ye    ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_mach  ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_nl    ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_ne    ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_gx    ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_gy    ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_gz    ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_x     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_comp  ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_nuc   ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_nu    ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_enu   ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_fnu   ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_BVw   ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_yl    ! character string containing name of model files

INTEGER                                  :: istat         ! open file flag
INTEGER, PARAMETER                       :: iunit = 20    ! unit number to print models
INTEGER                                  :: iunitp        ! unit number to print models
INTEGER                                  :: n_editMDmach  ! model number of mach number edits
INTEGER                                  :: n             ! neutrino flavor array index

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Unable to open',a128)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ PRINT IF PRINT CRITERIA SATISFIED /////
!
!-----------------------------------------------------------------------

!........x-velocities

IF ( l_editMDu ) THEN

  n_editMDu        = n_editMDu + 1

  WRITE (outfile_u,'(a21,i5.5,a2)') '/Models_MD_n/u/modelu', &
&  n_editMDu,'.d'
  outfile_u        = TRIM(data_path)//TRIM(outfile_u)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_u),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_u)

  iunitp           = iunit
  CALL editMDu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, uMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDu

!........y-velocities

IF ( l_editMDv ) THEN

  n_editMDv        = n_editMDv + 1

  WRITE (outfile_v,'(a21,i5.5,a2)') '/Models_MD_n/v/modelv', &
&  n_editMDv,'.d'
  outfile_v        = TRIM(data_path)//TRIM(outfile_v)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_v),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_v)

  iunitp           = iunit
  CALL editMDv( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, vMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDv

!........z-velocities

IF ( l_editMDw ) THEN

  n_editMDw        = n_editMDw + 1

  WRITE (outfile_w,'(a21,i5.5,a2)') '/Models_MD_n/w/modelw', &
&  n_editMDw,'.d'
  outfile_w        = TRIM(data_path)//TRIM(outfile_w)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_w),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_w)

  iunitp           = iunit
  CALL editMDw( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, wMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDw

!........entropy

IF ( l_editMDs ) THEN

  n_editMDs        = n_editMDs + 1

  WRITE (outfile_s,'(a21,i5.5,a2)') '/Models_MD_n/s/models', &
&  n_editMDs,'.d'
  outfile_s        = TRIM(data_path)//TRIM(outfile_s)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_s),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_s)

  iunitp           = iunit
  CALL editMDs( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, sMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDs

!........density

IF ( l_editMDd ) THEN

  n_editMDd        = n_editMDd + 1

  WRITE (outfile_d,'(a21,i5.5,a2)') '/Models_MD_n/d/modeld', &
&  n_editMDd,'.d'
  outfile_d        = TRIM(data_path)//TRIM(outfile_d)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_d),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_d)

  iunitp           = iunit
  CALL editMDd( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, rhoMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDd

!........pressure

IF ( l_editMDp ) THEN

  n_editMDp        = n_editMDp + 1

  WRITE (outfile_p,'(a21,i5.5,a2)') '/Models_MD_n/p/modelp', &
&  n_editMDp,'.d'
  outfile_p        = TRIM(data_path)//TRIM(outfile_p)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_p),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_p)

  iunitp           = iunit
  CALL editMDp( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, pMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDp

!........energy

IF ( l_editMDe ) THEN

  n_editMDe        = n_editMDe + 1

  WRITE (outfile_e,'(a21,i5.5,a2)') '/Models_MD_n/e/modele', &
&  n_editMDe,'.d'
  outfile_e        = TRIM(data_path)//TRIM(outfile_e)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_e),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_e)

  iunitp           = iunit
  CALL editMDe( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, eMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDe

!........neutrino energy density

IF ( l_editMDenu ) THEN

  n_editMDenu      = n_editMDenu + 1

  WRITE (outfile_enu,'(a25,i5.5,a2)') '/Models_MD_n/enu/modelenu', &
&  n_editMDenu,'.d'
  outfile_enu        = TRIM(data_path)//TRIM(outfile_enu)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_enu),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_enu)

  iunitp           = iunit
  CALL editMDenu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, e_nu_MD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDenu

!........neutrino flux

IF ( l_editMDfnu ) THEN

  n_editMDfnu      = n_editMDfnu + 1

  WRITE (outfile_fnu,'(a25,i5.5,a2)') '/Models_MD_n/fnu/modelfnu', &
&  n_editMDfnu,'.d'
  outfile_fnu        = TRIM(data_path)//TRIM(outfile_fnu)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_fnu),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_fnu)

  iunitp           = iunit
  CALL editMDfnu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, f_nu_MD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDfnu

!........mean mass number

IF ( l_editMDa ) THEN

  n_editMDa        = n_editMDa + 1

  WRITE (outfile_a,'(a21,i5.5,a2)') '/Models_MD_n/a/modela', &
&  n_editMDa,'.d'
  outfile_a        = TRIM(data_path)//TRIM(outfile_a)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_a),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_a)

  iunitp           = iunit
  CALL editMDa( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, aMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDa

!........critical bpundaries

IF ( l_editMDx ) THEN

  n_editMDx        = n_editMDx + 1

  WRITE (outfile_x,'(a21,i5.5,a2)') '/Models_MD_n/x/modelx', &
&  n_editMDx,'.d'
  outfile_x        = TRIM(data_path)//TRIM(outfile_x)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_x),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_x)

  iunitp           = iunit
  CALL editMDx( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, nnu, time, &
&  t_tb, y_in, z_in, r_nse, r_O1, r_xO, r_shock, r_shock_mn, r_shock_mx, &
&  rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, esphere_mean, &
&  tau_advMD, tau_heat_nuMD, tau_heat_nucMD, n_growMD, r_gain, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDx

!........electron fraction

IF ( l_editMDye ) THEN

  n_editMDye       = n_editMDye + 1

  WRITE (outfile_ye,'(a23,i5.5,a2)') '/Models_MD_n/ye/modelye', &
&  n_editMDye,'.d'
  outfile_ye       = TRIM(data_path)//TRIM(outfile_ye)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_ye),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_ye)

  iunitp           = iunit
  CALL editMDye( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, yeMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDye

!........mach number

IF ( l_editMDu  .and.  l_editMDv  .and.  l_editMDw ) THEN

  n_editMDmach     = MIN( n_editMDu, n_editMDv )
  WRITE (outfile_mach,'(a27,i5.5,a2)') '/Models_MD_n/mach/modelmach', &
&  n_editMDu,'.d'
  outfile_mach     = TRIM(data_path)//TRIM(outfile_mach)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_mach),STATUS='replace',IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_mach)

  iunitp           = iunit
  CALL editMDmach( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, v_csoundMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDu and l_editMDv and l_editMDw

!........neutrino luminosity

IF ( l_editMDnl ) THEN

  n_editMDnl       = n_editMDnl + 1

  DO n = 1,nnu
    WRITE (outfile_nl,'(a16,i1,a8,i1,a1,i5.5,a2)') '/Models_MD_n/nl_',n,'/modelnl', &
&    n,'_',n_editMDnl,'.d'
    outfile_nl     = TRIM(data_path)//TRIM(outfile_nl)
    OPEN (UNIT=iunit,FILE=TRIM(outfile_nl),STATUS='replace',IOSTAT=istat)
    IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_nl)
  END DO

  iunitp           = iunit
  CALL editMDnl( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, n, nnu, x_in, &
& y_in, z_in, time, t_tb, lumMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDnl

!........neutrino rms energy

IF ( l_editMDne ) THEN

  n_editMDne       = n_editMDne + 1

  DO n = 1,nnu
    WRITE (outfile_ne,'(a16,i1,a8,i1,a1,i5.5,a2)') '/Models_MD_n/ne_',n,'/modelne', &
&    n,'_',n_editMDne,'.d'
    outfile_ne     = TRIM(data_path)//TRIM(outfile_ne)
    OPEN (UNIT=iunit,FILE=TRIM(outfile_ne),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
    IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_ne)
  END DO

  iunitp           = iunit
  CALL editMDne( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, n, nnu, x_in, &
&  y_in, z_in, time, t_tb, e_rmsMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDne

!........x-gravitational acceleration

IF ( l_editMDgx ) THEN

  n_editMDgx       = n_editMDgx + 1

  WRITE (outfile_gx,'(a23,i5.5,a2)') '/Models_MD_n/gx/modelgx', &
&  n_editMDgx,'.d'
  outfile_gx       = TRIM(data_path)//TRIM(outfile_gx)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_gx),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_gx)

  iunitp           = iunit
  CALL editMDgx( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, grav_x_cMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDgx

!........y-gravitational acceleration

IF ( l_editMDgy ) THEN

  n_editMDgy       = n_editMDgy + 1

  WRITE (outfile_gy,'(a23,i5.5,a2)') '/Models_MD_n/gy/modelgy', &
&  n_editMDgy,'.d'
  outfile_gy       = TRIM(data_path)//TRIM(outfile_gy)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_gy),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_gy)

  iunitp           = iunit
  CALL editMDgy( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, grav_y_cMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDgy

!........z-gravitational acceleration

IF ( l_editMDgz ) THEN

  n_editMDgz       = n_editMDgz + 1

  WRITE (outfile_gz,'(a23,i5.5,a2)') '/Models_MD_n/gz/modelgz', &
&  n_editMDgz,'.d'
  outfile_gz       = TRIM(data_path)//TRIM(outfile_gz)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_gz),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_gz)

  iunitp           = iunit
  CALL editMDgz( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, grav_z_cMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDgz

!........dominant composition

IF ( l_editMDcm ) THEN

  n_editMDcm       = n_editMDcm + 1

  WRITE (outfile_comp,'(a27,i5.5,a2)') '/Models_MD_n/comp/modelcomp', &
&  n_editMDcm,'.d'
  outfile_comp   = TRIM(data_path)//TRIM(outfile_comp)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_comp),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_comp)

  iunitp         = iunit
  CALL editMDcomp( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, i_comp, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDcm

!........Energy generation rate by nuclear reactions

IF ( l_editMDnc ) THEN

  n_editMDnc       = n_editMDnc + 1

  WRITE (outfile_nuc,'(a35,i5.5,a2)') '/Models_MD_n/dudt_nuc/model_dudt_nuc', &
&  n_editMDnc,'.d'
  outfile_nuc   = TRIM(data_path)//TRIM(outfile_nuc)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_nuc),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_nuc)

  iunitp         = iunit
  CALL editMDdudt_nuc( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, dudt_nucMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDnc

!........Energy deposition rate by neutrinos

IF ( l_editMDnu ) THEN

  n_editMDnu       = n_editMDnu + 1

  WRITE (outfile_nu,'(a33,i5.5,a2)') '/Models_MD_n/dudt_nu/model_dudt_nu', &
&  n_editMDnu,'.d'
  outfile_nu    = TRIM(data_path)//TRIM(outfile_nu)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_nu),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_nu)

  iunitp         = iunit
  CALL editMDdudt_nu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, dudt_nuMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........Brunt-Vaisala_ye frequency

IF ( l_editMDBVw ) THEN

  n_editMDBVw      = n_editMDBVw + 1

  WRITE (outfile_BVw,'(a26,i5.5,a2)') '/Models_MD_t/BVw/model_BVw', &
&  n_editMDBVw,'.d'
  outfile_BVw   = TRIM(data_path)//TRIM(outfile_BVw)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_BVw),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_BVw)

  iunitp         = iunit
  CALL editMDBVw( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, wBVMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDBVw

!........lepton fraction

IF ( l_editMDyl ) THEN

  n_editMDyl       = n_editMDyl + 1

  WRITE (outfile_yl,'(a23,i5.5,a2)') '/Models_MD_n/yl/modelyl', &
&  n_editMDyl,'.d'
  outfile_yl       = TRIM(data_path)//TRIM(outfile_yl)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_yl),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_yl)

  iunitp           = iunit
  CALL editMDyl( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, ylMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_editMDyl


END IF ! l_editMDnu

!-----------------------------------------------------------------------
!
!             \\\\\ PRINT IF TIME CRITERIA SATISFIED /////
!
!-----------------------------------------------------------------------

IF ( l_MDedit ) THEN

  n_editMD         = n_editMD + 1

!........x-velocities

  WRITE (outfile_u,'(a21,i5.5,a2)') '/Models_MD_t/u/modelu', &
&  n_editMD,'.d'
  outfile_u        = TRIM(data_path)//TRIM(outfile_u)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_u),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_u)

  iunitp           = iunit
  CALL editMDu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, uMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........y-velocities

  WRITE (outfile_v,'(a21,i5.5,a2)') '/Models_MD_t/v/modelv', &
&  n_editMD,'.d'
  outfile_v        = TRIM(data_path)//TRIM(outfile_v)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_v),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_v)

  iunitp           = iunit
  CALL editMDv( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, vMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........z-velocities

  WRITE (outfile_w,'(a21,i5.5,a2)') '/Models_MD_t/w/modelw', &
&  n_editMD,'.d'
  outfile_w        = TRIM(data_path)//TRIM(outfile_w)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_w),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_w)

  iunitp           = iunit
  CALL editMDw( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, wMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........entropy

  WRITE (outfile_s,'(a21,i5.5,a2)') '/Models_MD_t/s/models', &
&  n_editMD,'.d'
  outfile_s        = TRIM(data_path)//TRIM(outfile_s)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_s),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_s)

  iunitp           = iunit
  CALL editMDs( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, sMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........density

  WRITE (outfile_d,'(a21,i5.5,a2)') '/Models_MD_t/d/modeld', &
&  n_editMD,'.d'
  outfile_d        = TRIM(data_path)//TRIM(outfile_d)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_d),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_d)

  iunitp           = iunit
  CALL editMDd( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, rhoMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........pressure

  WRITE (outfile_p,'(a21,i5.5,a2)') '/Models_MD_t/p/modelp', &
&  n_editMD,'.d'
  outfile_p        = TRIM(data_path)//TRIM(outfile_p)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_p),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_p)

  iunitp           = iunit
  CALL editMDp( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, pMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........energy

  WRITE (outfile_e,'(a21,i5.5,a2)') '/Models_MD_t/e/modele', &
&  n_editMD,'.d'
  outfile_e        = TRIM(data_path)//TRIM(outfile_e)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_e),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_e)

  iunitp           = iunit
  CALL editMDe( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, eMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........neutrino energy density

  WRITE (outfile_enu,'(a25,i5.5,a2)') '/Models_MD_t/enu/modelenu', &
&  n_editMD,'.d'
  outfile_enu        = TRIM(data_path)//TRIM(outfile_enu)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_enu),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_enu)

  iunitp           = iunit
  CALL editMDenu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, e_nu_MD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........neutrino flux

  WRITE (outfile_fnu,'(a25,i5.5,a2)') '/Models_MD_t/fnu/modelfnu', &
&  n_editMD,'.d'
  outfile_fnu        = TRIM(data_path)//TRIM(outfile_fnu)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_fnu),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_fnu)

  iunitp           = iunit
  CALL editMDfnu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, f_nu_MD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........mean mass number

  WRITE (outfile_a,'(a21,i5.5,a2)') '/Models_MD_t/a/modela', &
&  n_editMD,'.d'
  outfile_a        = TRIM(data_path)//TRIM(outfile_a)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_a),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_a)

  iunitp           = iunit
  CALL editMDa( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, aMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........critical boundaries

  WRITE (outfile_x,'(a21,i5.5,a2)') '/Models_MD_t/x/modelx', &
&  n_editMD,'.d'
  outfile_x        = TRIM(data_path)//TRIM(outfile_x)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_x),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_x)

  iunitp           = iunit
  CALL editMDx( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, nnu, time, &
&  t_tb, y_in, z_in, r_nse, r_O1, r_xO, r_shock, r_shock_mn, r_shock_mx, &
&  rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, esphere_mean, &
&  tau_advMD, tau_heat_nuMD, tau_heat_nucMD, n_growMD, r_gain, iunitp )

  CLOSE (unit=iunit,status='keep')

!........electron fraction

  WRITE (outfile_ye,'(a23,i5.5,a2)') '/Models_MD_t/ye/modelye', &
&  n_editMD,'.d'
  outfile_ye       = TRIM(data_path)//TRIM(outfile_ye)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_ye),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_ye)

  iunitp           = iunit
  CALL editMDye( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, yeMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........mach number

  WRITE (outfile_mach,'(a27,i5.5,a2)') '/Models_MD_t/mach/modelmach', &
&  n_editMD,'.d'
  outfile_mach     = TRIM(data_path)//TRIM(outfile_mach)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_mach),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_mach)

  iunitp           = iunit
  CALL editMDmach( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
& z_in, time, t_tb, v_csoundMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........neutrino luminosity

  DO n = 1,nnu
    WRITE (outfile_nl,'(a16,i1,a8,i1,a1,i5.5,a2)') '/Models_MD_t/nl_',n,'/modelnl', &
&    n,'_',n_editMD,'.d'
    outfile_nl     = TRIM(data_path)//TRIM(outfile_nl)
    OPEN (UNIT=iunit,FILE=TRIM(outfile_nl),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_nl)

    iunitp         = iunit
    CALL editMDnl( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, n, nnu, x_in, &
&    y_in, z_in, time, t_tb, lumMD, iunitp )
  END DO

  CLOSE (unit=iunit,status='keep')

!........neutrino rms energy

  DO n = 1,nnu
    WRITE (outfile_ne,'(a16,i1,a8,i1,a1,i5.5,a2)') '/Models_MD_t/ne_',n,'/modelne', &
&    n,'_',n_editMD,'.d'
    outfile_ne     = TRIM(data_path)//TRIM(outfile_ne)
    OPEN (UNIT=iunit,FILE=TRIM(outfile_ne),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_ne)

    iunitp         = iunit
    CALL editMDne( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, n, nnu, x_in, &
&    y_in, z_in, time, t_tb, e_rmsMD, iunitp )
  END DO

  CLOSE (unit=iunit,status='keep')

!........x-gravitational acceleration

  WRITE (outfile_gx,'(a23,i5.5,a2)') '/Models_MD_t/gx/modelgx', &
&  n_editMD,'.d'
  outfile_gx       = TRIM(data_path)//TRIM(outfile_gx)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_gx),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_gx)

  iunitp           = iunit
  CALL editMDgx( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, grav_x_cMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........y-gravitational acceleration

  WRITE (outfile_gy,'(a23,i5.5,a2)') '/Models_MD_t/gy/modelgy', &
&  n_editMD,'.d'
  outfile_gy       = TRIM(data_path)//TRIM(outfile_gy)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_gy),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_gy)

  iunitp           = iunit
  CALL editMDgy( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, grav_y_cMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........z-gravitational acceleration

  WRITE (outfile_gz,'(a23,i5.5,a2)') '/Models_MD_t/gz/modelgz', &
&  n_editMD,'.d'
  outfile_gz       = TRIM(data_path)//TRIM(outfile_gz)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_gz),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_gz)

  iunitp           = iunit
  CALL editMDgz( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, grav_z_cMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........dominant composition

  WRITE (outfile_comp,'(a27,i5.5,a2)') '/Models_MD_t/comp/modelcomp', &
&  n_editMD,'.d'
  outfile_comp   = TRIM(data_path)//TRIM(outfile_comp)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_comp),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_comp)

  iunitp         = iunit
  CALL editMDcomp( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, i_comp, iunitp )

  CLOSE (unit=iunit,status='keep')

!........Energy generation rate by nuclear reactions

  WRITE (outfile_nuc,'(a35,i5.5,a2)') '/Models_MD_t/dudt_nuc/model_dudt_nuc', &
&  n_editMD,'.d'
  outfile_nuc   = TRIM(data_path)//TRIM(outfile_nuc)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_nuc),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_nuc)

  iunitp         = iunit
  CALL editMDdudt_nuc( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, dudt_nucMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........Energy deposition rate by neutrinos

  WRITE (outfile_nu,'(a33,i5.5,a2)') '/Models_MD_t/dudt_nu/model_dudt_nu', &
&  n_editMD,'.d'
  outfile_nu    = TRIM(data_path)//TRIM(outfile_nu)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_nu),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_nu)

  iunitp         = iunit
  CALL editMDdudt_nu( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, dudt_nuMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........Brunt-Vaisala_ye frequency

  WRITE (outfile_BVw,'(a26,i5.5,a2)') '/Models_MD_t/BVw/model_BVw', &
&  n_editMD,'.d'
  outfile_BVw   = TRIM(data_path)//TRIM(outfile_BVw)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_BVw),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_BVw)

  iunitp         = iunit
  CALL editMDBVw( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, wBVMD, iunitp )

  CLOSE (unit=iunit,status='keep')

!........lepton fraction

  WRITE (outfile_yl,'(a23,i5.5,a2)') '/Models_MD_t/yl/modelyl', &
&  n_editMD,'.d'
  outfile_yl       = TRIM(data_path)//TRIM(outfile_yl)
  OPEN (UNIT=iunit,FILE=TRIM(outfile_yl),STATUS='replace', FORM='unformatted', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001) TRIM(outfile_yl)

  iunitp           = iunit
  CALL editMDyl( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, y_in, &
&  z_in, time, t_tb, ylMD, iunitp )

  CLOSE (unit=iunit,status='keep')

END IF ! l_MDedit

RETURN
END SUBROUTINE edit_MD_exec
