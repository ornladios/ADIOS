SUBROUTINE radhyd_to_edit_HDF( i_ray_min, i_ray_max, i_ray_dim, nx, ny, nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_edit_HDF
!    Module:       radhyd_to_edit_HDF
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To load variables into 3-d arrays for porting into MGFLD multi-D edit
!         via subroutine editMD_in.
!
!    Subprograms called:
!  editMD_in   : calculates quantities for the multi-D edit dumps
!  editMD_out  : writes multi-D data to HDF files for postprocessing
!
!    Input arguments:
!  i_ray_min   : minimum ray index
!  i_ray_max   : maximum ray index
!  i_ray_dim   : number of rays to put on a processor
!  nx          : x_array extent
!  nx          : y_array extent
!  nez         : neutrino energy array extent
!  nnu         : neutrino flavor array extent
!  nnc         : neutrino abundance array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  radhyd_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, rho_c, t_c, ye_c, &
& rhobar, u_c, v_c, w_c, x_e=>x_ef, x_c=>x_cf, y_e=>y_ef, y_c=>y_cf, &
& psi0_c, psi1_e, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, nse_c, &
& dtnph, time, ncycle
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: i_ray_min        ! minimum index denoting a specific radial ray
INTEGER, INTENT(in)              :: i_ray_max        ! maximum index denoting a specific radial ray
INTEGER, INTENT(in)              :: i_ray_dim        ! number of rays to put on a processor
INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: ny               ! y-array extent
INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc              ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.
INTEGER                          :: numcall = 0      ! minimum x-array index
INTEGER                          :: ret              !
INTEGER                          :: is               ! minimum x-array index
INTEGER                          :: ie               ! maximum x-array index
INTEGER                          :: js               ! minimum y-array index
INTEGER                          :: je               ! maximum y-array index
INTEGER                          :: i_ray            ! ndex denoting a specific radial ray
INTEGER                          :: editMD_in        ! return code
REAL(KIND=double), DIMENSION(nx,nnu,i_ray_dim)      :: e_rms_stat   ! sqrt[ SUM psi0 * w5dw/SUM w3ww ] (MeV)
REAL(KIND=double), DIMENSION(nx,nnu,i_ray_dim)      :: e_rms_trns   ! sqrt[ SUM psi1 * w5dw/SUM w3ww ] (MeV)

REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: e_rms_r_stat ! rms static neutrino energy at radius r_e_rms (MeV)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: e_rms_r_trns ! rms transport neutrino energy at radius r_e_rms (MeV)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: e_rms_d_stat ! rms static neutrino energy at density d_e_rms (MeV)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: e_rms_d_trns ! rms transport neutrino energy at density d_e_rms (MeV)

REAL(KIND=double), DIMENSION(nez,nnu,i_ray_dim)     :: r_sphere     ! k,n neutrinosphere radius (cm)
REAL(KIND=double), DIMENSION(nez,nnu,i_ray_dim)     :: d_sphere     ! k,n  neutrinosphere density (g cm^{-3})
REAL(KIND=double), DIMENSION(nez,nnu,i_ray_dim)     :: t_sphere     ! k,n  neutrinosphere temperature (MeV)
REAL(KIND=double), DIMENSION(nez,nnu,i_ray_dim)     :: m_sphere     ! k,n  neutrinosphere enclosed mass (g)

REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: rsphere_mean ! mean neutrinosphere radius (cm)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: dsphere_mean ! mean neutrinosphere density (g cm^{-3})
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: tsphere_mean ! mean neutrinosphere temperature (MeV)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: msphere_mean ! mean neutrinosphere enclosed mass (g)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: esphere_mean ! mean neutrinosphere energy (MeV)


REAL(KIND=double), DIMENSION(nx,nnu,i_ray_dim)      :: lum          ! n-neutrino luminosity (foes)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: lum_r        ! n-neutrino luminosity at radius          (foes)
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: lum_rho      ! n-neutrino luminosity at density d_lum (foes)

REAL(KIND=double), DIMENSION(nx,nnu,i_ray_dim)      :: inv_fluxfact ! inverse flux factor

REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: stat_pinch_r ! static spectral pinch factor at radius r_pinch
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: trns_pinch_r ! transport spectral pinch factor at radius r_pinch
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: stat_pinch_d ! static spectral pinch factor at density d_pinch
REAL(KIND=double), DIMENSION(nnu,i_ray_dim)         :: trns_pinch_d ! transport spectral pinch factor at density d_pinch

REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: vel_x        ! x-component of velocity (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: vel_y        ! y-component of velocity (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: vel_z        ! z-component of velocity (cm s^{-1})
REAL(KIND=double), DIMENSION(nx+1)                  :: xi           ! zone-edged x-coordinate (cm)
REAL(KIND=double), DIMENSION(nx)                    :: xiph         ! zone-centered x-coordinate (cm)
REAL(KIND=double), DIMENSION(ny+1)                  :: thetai       ! zone-edged y-coordinate
REAL(KIND=double), DIMENSION(ny)                    :: thetaiph     ! zone-centered y-coordinate
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: rho_ed       ! density (g cm^{-3})
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: t_ed         ! temperature (MeV)
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: ye_ed        ! electron fraction
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: s_ed         ! entropy
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: yl_ed        ! lepton fraction
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: yenu         ! electron neutrino fraction
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: yenubar      ! electron antineutrino fraction
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: yxnu         ! muon or tau neutrino fraction
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: yxnubar      ! muon or tau antineutrino fraction
REAL(KIND=double), DIMENSION(nx,i_ray_dim)          :: dm_ed        ! zone mass (solar masses)
REAL(KIND=double), DIMENSION(nez)                   :: e_nu_center  ! zone-centered neutrino energy
REAL(KIND=double), DIMENSION(nez+1)                 :: e_nu_edge    ! zone-edged neutrino energy


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Load integer scalars...........................................

is                    = imin
ie                    = imax
js                    = jmin
je                    = jmax

!-----------------------------------------------------------------------
!
!            \\\\\ TRANSFER VARIABLES TO EDIT_HDF /////
!
!-----------------------------------------------------------------------

!........Begin loop over radial rays....................................
!.......................................................................

DO i_ray = i_ray_min,i_ray_max

!........Transfer variables to MGFLD edit

  ret = editMD_in( is, ie, nx, js, je, ny, i_ray, i_ray_dim, nez, nnu, nnc, &
&  rho_c, t_c, ye_c, rhobar, x_e, x_c, y_e, y_c, u_c, v_c, w_c, psi0_c, &
&  psi1_e, dtnph, time, ncycle, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, &
&  nse_c, first, e_rms_stat, e_rms_trns, e_rms_r_stat, e_rms_r_trns, e_rms_d_stat, &
&  e_rms_d_trns, rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, &
&  esphere_mean, r_sphere, d_sphere, t_sphere, m_sphere, lum, lum_r, lum_rho, &
&  inv_fluxfact, stat_pinch_r, trns_pinch_r, stat_pinch_d, trns_pinch_d, &
&  vel_x, vel_y, vel_z, xi, xiph, thetai, thetaiph, rho_ed, t_ed, ye_ed, &
&  s_ed, yl_ed, dm_ed, yenu, yenubar, yxnu, yxnubar, e_nu_center, e_nu_edge)

!-----------------------------------------------------------------------
!
!          \\\\\ TRANSFER VARIABLES FROM MGFLD EDIT_HDF /////
!
!-----------------------------------------------------------------------
 
  IF (ret == 0)  THEN
    CALL HDFeditMD_out( is, ie, nx, js, je, ny, i_ray, i_ray_dim, nez, nnu, nnc, rho_c, &
&    t_c, ye_c,  x_e, x_c, y_e, y_c, u_c, v_c, w_c, psi0_c, psi1_e, dtnph, time, &
&    ncycle, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, nse_c, numcall, e_rms_stat, &
&    e_rms_trns, e_rms_r_stat, e_rms_r_trns, e_rms_d_stat, e_rms_d_trns, &
&    rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, esphere_mean, &
&    r_sphere, d_sphere, t_sphere, m_sphere, lum, lum_r, lum_rho, inv_fluxfact, &
&    stat_pinch_r, trns_pinch_r, stat_pinch_d, trns_pinch_d, &
&    vel_x, vel_y, vel_z, xi, xiph, thetai, thetaiph, rho_ed, t_ed, ye_ed, &
&    s_ed, yl_ed, dm_ed, yenu, yenubar, yxnu, yxnubar, e_nu_center, e_nu_edge )

    IF ( i_ray == i_ray_max ) numcall = numcall + 1
  END IF

  ret                  = 8
!........End loop over radial rays......................................
!.......................................................................

END DO

first                  = .false.

RETURN
END SUBROUTINE radhyd_to_edit_HDF
