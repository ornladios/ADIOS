SUBROUTINE HDFeditHDF_out( is, ie, nx, js, je, ny, i_ray, i_ray_dim, nez, nnu, nnc, rhop,&
& tp, yep, rp, rpph, thetap, thetapph, up, vp, wp, psi0p, psi1p, dtime, time_ellapsed, &
& cycle_number, xnp, be_nucp, a_nucp, z_nucp, nsep, numcall, e_rms_stat, &
& e_rms_trns, e_rms_r_stat, e_rms_r_trns, e_rms_d_stat, e_rms_d_trns, &
& rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, esphere_mean, &
& r_sphere, d_sphere, t_sphere, m_sphere, lum, lum_r, lum_rho, inv_fluxfact, & 
& stat_pinch_r, trns_pinch_r, stat_pinch_d, trns_pinch_d, &
& vel_x, vel_y, vel_z, xi, xiph, thetai, thetaiph, rho_ed, t_ed, ye_ed, &
& s_ed, yl_ed, dm_ed, yenu, yenubar, yxnu, yxnubar, e_nu_center, e_nu_edge )
!-----------------------------------------------------------------------
!
!    File:         HDFeditHDF_out
!    Module:       HDFeditHDF_out 
!    Type:         Subprogram
!    Author:       Charlotte Dirk, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         03/12/05
!
!    Purpose:
!        To WRITE HDF file for quantities for multi-D edit dumps.
!
!    Subprograms called:
!   
!    Input arguments:
!
!  is             : inner x-array index
!  ie             : outer x-array index
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
!  be_nucp        : binding energy of mean heavy nucleus (MeV)
!  a_nucp         : mass number of mean heavy nucleus
!  z_nucp         : charge number of mean heavy nucleus
!  nesp           : nuclear statistical equilibrium flag
!  first          : initial call flag  (not used)
!
!	Array arguments:
!
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
!  xi             : zone-edged    x-coordinate (cm)
!  xiph           : zone-centered x-coordinate (cm)
!  thetai         : zone-edged    y-coordinate
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
!    Output arguments:
!
!    Include files:
!      kind_module, numerical_module, physcnst_module
!      cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!      nucbrn_module, nu_dist_module, nu_energy_grid_module,
!      parallel_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module , ONLY : double, single
!USE param
!USE root
!use grid
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : kmev, msolar

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : r_lum, d_lum, r_e_rms, d_e_rms, i_HDFedit, dt_HDFedit1, dt_HDFedit2, &
& r_pinch, d_pinch
USE eos_snc_x_module, ONLY : xn, be_nuc, a_nuc, z_nuc, nse, aesv
USE mdl_cnfg_module, ONLY : rho, t, ye, dr, r, u, dmrst, rstmss, jr_max
USE nucbrn_module, ONLY : xn_n=>xn, be_nuc_n=>be_nuc, a_nuc_n=>a_nuc, z_nuc_n=>z_nuc, &
& fescrn, fascrn, uburn_n=>uburn, nse_n=>nse
USE nu_dist_module, ONLY : psi0, psi1, vol, rjmh, j_sphere, r_spherep => r_sphere, &
& d_spherep => d_sphere, t_spherep => t_sphere, m_spherep => m_sphere
USE parallel_module, ONLY : myid
USE t_cntrl_module, ONLY: dtime_hydro, dtnph, time, t_bounce

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                             :: numcall       ! inner x-array index

INTEGER, INTENT(in)                             :: is            ! inner x-array index
INTEGER, INTENT(in)                             :: ie            ! outer x-array index
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

REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: rhop         ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: tp           ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: yep          ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                 :: rp           ! zone-edged radial coordinate (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)                 :: thetap       ! zone-edged y (angular) coordinate
REAL(KIND=double), INTENT(in), DIMENSION(nx)                   :: rpph         ! zone-centered radial coordinate (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny)                   :: thetapph     ! zone-centered y (angular) coordinate
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: up           ! x-component of zone-centered velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: vp           ! y-component of zone-centered velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: wp           ! z-component of zone-centered velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,i_ray_dim) :: psi0p        ! zeroth angular moment of the NDF
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,i_ray_dim) :: psi1p        ! first angular moment of the NDF

REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc,i_ray_dim)     :: xnp          ! Composition mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: be_nucp      ! binding energy of mean heavy nucleus (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: a_nucp       ! mass number of mean heavy nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)         :: z_nucp       ! charge number of mean heavy nucleus

REAL(KIND=double), INTENT(in)                                  :: dtime         ! time step used at current cycle
REAL(KIND=double), INTENT(in)                                  :: time_ellapsed ! ellapsed time
!
!----------------------------------------------------------------------
!
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,i_ray_dim)    :: e_rms_stat   ! sqrt[ SUM psi0 * w5dw/SUM w3ww ] (MeV) (RMS static neutrino energy)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,i_ray_dim)    :: e_rms_trns   ! sqrt[ SUM psi1 * w5dw/SUM w3ww ] (MeV) (RMS transport neutrino energy)

REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: e_rms_r_stat ! rms static neutrino energy at radius r_e_rms (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: e_rms_r_trns ! rms transport neutrino energy at radius r_e_rms (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: e_rms_d_stat ! rms static neutrino energy at density d_e_rms (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: e_rms_d_trns ! rms transport neutrino energy at density d_e_rms (MeV)

REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,i_ray_dim)   :: r_sphere     ! k,n neutrinosphere radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,i_ray_dim)   :: d_sphere     ! k,n neutrinosphere density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,i_ray_dim)   :: t_sphere     ! k,n neutrinosphere temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,i_ray_dim)   :: m_sphere     ! k,n neutrinosphere enclosed mass (g)

REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: rsphere_mean ! mean n-neutrinosphere radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: dsphere_mean ! mean n-neutrinosphere density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: tsphere_mean ! mean n-neutrinosphere temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: msphere_mean ! mean n-neutrinosphere enclosed mass (g)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: esphere_mean ! mean n-neutrinosphere energy at the mean nu_sphere (MeV)

REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,i_ray_dim)    :: lum          ! n-luminosity (foes)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: lum_r        ! n-luminosity at radius r_lum (foes)
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: lum_rho      ! n-luminosity at density d_lum (foes)

REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,i_ray_dim)    :: inv_fluxfact ! inverse flux factor

REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: stat_pinch_r ! static spectral pinch factor at radius r_pinch
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: trns_pinch_r ! transport spectral pinch factor at radius r_pinch
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: stat_pinch_d ! static spectral pinch factor at density d_pinch
REAL(KIND=double), INTENT(in), DIMENSION(nnu,i_ray_dim)       :: trns_pinch_d ! transport spectral pinch factor at density d_pinch

REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: vel_x        ! x-component of velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: vel_y        ! y-component of velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: vel_z        ! z-component of velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                  :: xi           ! zone-edged x-coordinate (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                  :: xiph         ! zone-centered x-coordinate (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                  :: thetai       ! zone-edged y-coordinate
REAL(KIND=double), INTENT(in), DIMENSION(nx)                  :: thetaiph     ! zone-centered y-coordinate
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: rho_ed       ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: t_ed         ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: ye_ed        ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: s_ed         ! entropy
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: yl_ed        ! lepton fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: yenu         ! electron neutrino fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: yenubar      ! electron antineutrino fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: yxnu         ! muon or tau neutrino fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: yxnubar      ! muon or tau antineutrino fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,i_ray_dim)        :: dm_ed        ! zone mass (solar masses)

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nez)                 :: e_nu_center  ! zone-centered neutrino energy
REAL(KIND=double), INTENT(out), DIMENSION(nez+1)               :: e_nu_edge    ! zone-edged neutrino energy

      CHARACTER*10  :: digits  = '0123456789'       
      CHARACTER*21  :: filename1 = 'Data3/HDF/HDF_radhyd_'  ! constant
      CHARACTER*3   :: filename2 = '000'               ! ray number
      CHARACTER*3   :: filename3 = '000'               ! time sequence
      CHARACTER*28  :: filename                        ! full filename
      CHARACTER*16  :: coordsys
      CHARACTER*80  :: string
      CHARACTER*9   :: form = '(12e11.4)'
      CHARACTER*12  :: units
      
      INTEGER       :: i, j, k, n1, n2, n3, indx
      INTEGER       :: rank,shape(3),ret
      REAL(KIND=single) :: data(nx*ny), xscale(nx), yscale(ny), nuscale(nnu), ezscale(nez)
      REAL(KIND=single) :: xscale1(nx+1), yscale1(ny+1), ezscale1(nez+1)

      INTEGER       :: dssdims,dssdast,dssdisc,dsadata,dspdata
      external         dssdims,dssdast,dssdisc,dsadata,dspdata

      ! **********************
      Print *,'start - HDFeditHDF_out '

      ! *************************
      coordsys = 'spherical polar' // char(0)
      
      n1 = numcall / 100
      n2 = numcall / 10 - n1 * 10
      n3 = numcall - n1 * 100 - n2 * 10
      filename3 = digits(n1+1:n1+1)//digits(n2+1:n2+1)//digits(n3+1:n3+1)
      
      indx = i_ray_dim * myid + i_ray ! zero to i_ray;  i_ray+1 to 2*i_ray; etc
      n1 = indx / 100
      n2 = indx / 10 - n1 * 10
      n3 = indx - n1 * 100 - n2 * 10
      filename2 = digits(n1+1:n1+1)//digits(n2+1:n2+1)//digits(n3+1:n3+1)
      
      filename = filename1//filename2//'_'//filename3

      DO  i=is,ie
         xscale(i-is+1) = real( rpph(i))   ! zone centered
      ENDDO
                     ! single ray
         yscale(1) = real( thetapph(indx))   
      DO  i=is,ie+1
         xscale1(i-is+1) = real( rp(i))    ! zone edged (future use)
      ENDDO
      DO  j=js,je+1
         yscale1(j ) = real( thetap(indx))   
      ENDDO
      j = 1
      DO  k=1,nnu
         nuscale(k) = real(k)
      ENDDO
      DO  k=1,nez
         ezscale(k) = real(e_nu_center(k))
      ENDDO
      DO  k=1,nez+1
         ezscale1(k) = real(e_nu_edge(k))
      ENDDO
          
      !shape  of zone-edge matrix
      rank     = 2
      shape(1) = ie-is+2
      shape(2) = je-js+2

      ret = dssdims(rank,shape)
      ret = dssdisc(1,shape(1),xscale1)
      ret = dssdisc(2,shape(2),yscale1)
      
      !shape of zone-center matrix 
      rank     = 2
      shape(1) = ie-is+1
      shape(2) = je-js+1
      !shape of zone-center matrix with single ray  (the shape usually used)
      shape(2) = 1
      
      ret = dssdims(rank,shape)
      ret = dssdisc(1,shape(1),xscale)
      ret = dssdisc(2,shape(2),yscale)
!*************************************************************
!******************   write single ray   *********************
!*************************************************************
!  1-velocity - x or r -component of velocity (cm s^{-1})
!
      units = 'cm / s'
         DO  i=is,ie
            data(i) = vel_x (i,i_ray) 
         ENDDO
      WRITE(string,"('1-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dspdata(filename,rank,shape,data)
!
!  2-velocity
!
         DO  i=is,ie
            data(i) = vel_y (i,i_ray)
         ENDDO
      WRITE(string,"('2-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)
!
!  3-velocity
!
         DO  i=is,ie
            data(i) = vel_z (i,i_ray) 
         ENDDO
      WRITE(string,"('3-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  density (g cm^{-3})
!
      units = 'g / (cm^3)'
         DO  i=is,ie
            data(i) =   rho_ed (i,i_ray) 
         ENDDO
      WRITE(string,"('DENSITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)
!
!   temperature (MeV)
!
      units = 'MeV'
         DO  i=is,ie
            data(i) =   t_ed (i,i_ray)
         ENDDO
      WRITE(string,"('TEMPERATURE AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  electron fraction
!
      units = ' - '
         DO  i=is,ie
            data(i) =   ye_ed  (i,i_ray) 
         ENDDO
      WRITE(string,"('ELECTRON FRACTION AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
! entropy
!
         DO  i=is,ie
            data(i) = s_ed (i,i_ray) 
         ENDDO
      WRITE(string,"('ENTROPY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  lepton fraction
!
         DO  i=is,ie
            data(i) =   yl_ed (i,i_ray) 
         ENDDO
      WRITE(string,"('LEPTON FRACTION AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  electron neutrino fraction
!
         DO  i=is,ie
            data(i) =   yenu (i,i_ray) 
         ENDDO
      WRITE(string,"('ELECTRON NEUTRINO FRACTION AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  electron antineutrino fraction
!
         DO  i=is,ie
            data(i) =  yenubar (i,i_ray)
         ENDDO
      WRITE(string,"('ELECTRON ANTINEUTRINO FRACTION AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  muon or tau neutrino fraction
!
         DO  i=is,ie
            data(i) =   yxnu (i,i_ray) 
         ENDDO
      WRITE(string,"('MUON OR TAU NEUTRINO FRACTION AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!   muon or tau antineutrino fraction
!
         DO  i=is,ie
            data(i) =   yxnubar  (i,i_ray) 
         ENDDO
      WRITE(string,"('MUON OR TAU ANTINEUTRINO FRACTION AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  zone mass (g)
!
      units = 'g'
         DO  i=is,ie
            data(i) =   dm_ed (i,i_ray) 
         ENDDO
      WRITE(string,"('zone mass (g) AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!--------
      rank     = 3
      shape(1) = ie-is+1
      shape(2) = nnu
      shape(3) = 1
      
      ret = dssdims(rank,shape)
      ret = dssdisc(1,shape(1),xscale)
      ret = dssdisc(2,shape(2),nuscale)
      ret = dssdisc(3,shape(3),yscale)
!--------
!  RMS static neutrino energy (MEV)
!
      units = 'MeV'
         DO  k=1,nnu
            DO  i=is,ie
               indx =  + (k)*shape(1) + (i-is) + 1
               data(indx) =   e_rms_stat (i,k,i_ray) 
            ENDDO
         ENDDO
      WRITE(string,"('RMS STATIC NEUTRINO ENERGY (MEV) AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  RMS transport neutrino energy (MEV)
!
         DO  k=1,nnu
            DO  i=is,ie
               indx =   + (k)*shape(1) + (i-is) + 1
               data(indx) =   e_rms_trns (i,k,i_ray)  
            ENDDO
         ENDDO
      WRITE(string,"('RMS TRANSPORT NEUTRINO ENERGY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!--------
      rank     = 2
      shape(1) = nnu
      shape(2) = 1

      ret = dssdims(rank,shape)
      ret = dssdisc(1,shape(1),nuscale)
      ret = dssdisc(2,shape(2),yscale)
!--------
!
! rms static neutrino energy at radius r_e_rms (MeV)
!
         DO  k=1,nnu
!         DO  i=is,ie
            indx =  k
            data(indx) =  e_rms_r_stat (k,i_ray)
         ENDDO
      WRITE(string,"('RMS STATIC NEUTRINO ENERGY AT RADIUS R_E_RMS AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  rms transport neutrino energy at radius r_e_rms (MeV)
!
         DO  k=1,nnu
 !        DO  i=is,ie
            indx =  k
            data(indx) =  e_rms_r_trns  (k,i_ray) 
         ENDDO
      WRITE(string,"('RMS TRANSPORT NEUTRINO ENERGY AT RADIUS R_E_RMS AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
! rms static neutrino energy at density d_e_rms (MeV)
!
         DO  k=1,nnu
!         DO 270 i=is,ie
            indx =  k  
            data(indx) =  e_rms_r_stat (k,i_ray) 
         ENDDO

      WRITE(string,"('RMS STATIC NEUTRINO ENERGY AT DENSITY D_E_RMS AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  rms transport neutrino energy at density d_e_rms (MeV) 
!
         DO  k=1,nnu
            indx =  k 
            data(indx) =   e_rms_r_trns  (k,i_ray) 
         ENDDO

      WRITE(string,"('RMS TRANSPORT NEUTRINO ENERGY AT DENSITY D_E_RMS AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!------------
!
!  static spectral pinch factor at radius r_pinch
!
      units = ' - '
         DO  k=1,nnu
            indx =   k 
             data(indx) = stat_pinch_r (k,i_ray) 
         ENDDO
      WRITE(string, "('STATIC SPECTRAL PINCH FACTOR AT RADIUS R_PINCH AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
! transport spectral pinch factor at radius r_pinch
!
         DO  k=1,nnu
            indx = k 
            data(indx) =  trns_pinch_r (k,i_ray) 
         ENDDO
      WRITE(string, "('TRANSPORT SPECTRAL PINCH FACTOR AT RADIUS R_PINCH AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  static spectral pinch factor at density d_pinch
!
         DO  k=1,nnu
            indx =  k 
            data(indx) =   stat_pinch_d  (k,i_ray) 
         ENDDO
      WRITE(string, "('STATIC SPECTRAL PINCH FACTOR AT DENSITY D_PINCH AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  transport spectral pinch factor at density d_pinch
!
         DO  k=1,nnu
            indx =  k 
            data(indx) =  trns_pinch_d (k,i_ray) 
         ENDDO
      WRITE(string,"('TRANSPORT SPECTRAL PINCH FACTOR AT DENSITY D_PINCH AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!--------
      rank     = 3
      shape(1) = nez
      shape(2) = nnu
      shape(3) = 1

      ret = dssdims(rank,shape)
      ret = dssdisc(1,shape(1),ezscale)
      ret = dssdisc(2,shape(2),nuscale)
      ret = dssdisc(3,shape(3),yscale)
!--------
!
!  k,n neutrinosphere radius (cm)
!
      units = 'cm'
         DO  k=1,nnu
            DO  i=1,nez
               indx =   (k-1)*nez + i
               data(indx) =   r_sphere  (i,k,i_ray) 
            ENDDO
         ENDDO
      WRITE(string,"('K,N NEUTRINOSPHERE RADIUS AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  k,n neutrinosphere density (g cm^{-3})
!
      units = 'g / (cm^3)'
         DO  k=1,nnu
            DO  i=1,nez
               indx = (k-1)*nez + i
               data(indx) =   d_sphere  (i,k,i_ray) 
            ENDDO
         ENDDO
      WRITE(string,"('K,N NEUTRINOSPHERE DENSITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!  mean n-neutrinosphere temperature (MeV)
!
      units = 'MeV'
         DO  k=1,nnu
            DO  i=1,nez
               indx = (k-1)*nez + i
               data(indx) =   t_sphere  (i,k,i_ray) 
            ENDDO
         ENDDO
      WRITE(string, "(' MEAN N-NEUTRINOSPHERE TEMPERATURE AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

!
!   k,n neutrinosphere enclosed mass (g)
!
      units = 'g'
         DO  k=1,nnu
            DO  i=1,nez
               indx =( k-1)*nez + i
               data(indx) =  m_sphere (i,k,i_ray)/msolar
            ENDDO
      ENDDO
      WRITE(string, "('K,N NEUTRINOSPHERE  ENCLOSED MASS AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,units,form,coordsys)
      ret = dsadata(filename,rank,shape,data)

      RETURN
      END SUBROUTINE HDFeditHDF_out