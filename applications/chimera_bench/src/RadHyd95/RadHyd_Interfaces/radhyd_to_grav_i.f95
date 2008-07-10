SUBROUTINE radhyd_to_grav_i( nx, ny, nz, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_grav_i
!    Module:       radhyd_to_grav_i
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To compute the gravitational potentials and accelerations at the
!        beginning of a cycle and return and store them in radial_ray_module
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x_array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  ny         : y_array extent
!  nz         : z_array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, mdl_cnfg_module, nu_dist_module, prb_cntl_module,
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, frpi, third

USE edit_module, ONLY : nlog
USE mdl_cnfg_module, ONLY : grvmss, rstmss, dmgrv, dmrst, agr, agrh
USE nu_dist_module, ONLY : unue, dunue, unube, unu, dunu, unub
USE prb_cntl_module, ONLY : irelhy, ilapsehy, i_grav
USE radial_ray_module, ONLY : imin, imax, x_ei, x_ci, dx_ci, y_ei, y_ci, &
& dy_ci, z_ei, dz_ci, rho_c, vx_bar, v2_bar, rhobar, e_bar, e_nu_c_bar,  &
& f_nu_e_bar, p_bar, agr_e, agr_c, agr_e_r, agr_c_r, grav_x_c, grav_y_c, &
& grav_z_c, grav_pot_c, gtot_pot_c, grav_pot_c_i, grav_x_e, grav_y_e,    &
& grav_z_e, grav_pot_e, gtot_pot_e, unu_c, unub_c, dunu_c, unue_e,       &
& unube_e, dunue_e

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                      :: nx          ! x-array extent
INTEGER, INTENT(in)                      :: ij_ray_dim  ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                      :: ik_ray_dim  ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                      :: ny          ! y-array extent
INTEGER, INTENT(in)                      :: nz          ! z-array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                  :: jr_min      ! inner shifted x-array index
INTEGER                                  :: jr_max      ! outer shifted x-array index
INTEGER                                  :: ij_ray      ! j-index of a radial ray
INTEGER                                  :: ik_ray      ! k-index of a radial ray
INTEGER                                  :: mode        ! initial evaluation of Legendre polynomials switch

REAL(KIND=double), DIMENSION(nx)         :: rho         ! density (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)         :: p           ! angular averaged pressure (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx+1)       :: r           ! radius (cm)
REAL(KIND=double), DIMENSION(nx)         :: gpot_e      ! gravitational potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gpot_c      ! grav. force defined at zone-centered (dynes g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gforce_e    ! gravitational potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gforce_c    ! grav. force defined at zone-centered (dynes g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gpot_dis_e  ! gravitational disassembly potential defined at zone-edge (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gpot_dis_c  ! gravitational disassembly potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gpotr_e     ! PN gravitational potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gpotr_c     ! PN grav. force defined at zone-centered (dynes g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gforcer_e   ! PN gravitational potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gforcer_c   ! PN grav. force defined at zone-centered (dynes g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gpotn_e     ! Newt gravitational potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gpotn_c     ! Newt grav. force defined at zone-centered (dynes g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gforcen_e   ! Newt gravitational potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)         :: gforcen_c   ! Newt grav. force defined at zone-centered (dynes g^{-1})

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' ilapsehy = ',i4,' is not supported at this time')
 1003 FORMAT (' irelhy = ',i4,' is not supported at this time')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Initialize
!-----------------------------------------------------------------------

jr_min                          = imin + 1
jr_max                          = imax + 1
grav_x_c                        = zero
grav_y_c                        = zero
grav_pot_c                      = zero
grav_pot_c_i                    = zero
grav_x_e                        = zero
grav_y_e                        = zero
grav_pot_e                      = zero
gtot_pot_c                      = zero
gtot_pot_e                      = zero

!-----------------------------------------------------------------------
!  Store previous values of the lapse function
!-----------------------------------------------------------------------

agr_e_r                           = agr_e
agr_c_r                           = agr_c

!-----------------------------------------------------------------------
!
!       \\\\\ TRANSFER VARIABLES TO SHIFTED MGFLD ARRAYS /////
!
!-----------------------------------------------------------------------

rho(jr_min:jr_max)              = rhobar (imin:imax)
p  (jr_min:jr_max)              = p_bar  (imin:imax)
r  (imin:imax+1)                = x_ei   (imin:imax+1)

dmrst(jr_min:jr_max)            = frpi * dx_ci(imin:imax) * ( r(imin:imax) * r(imin+1:imax+1) &
&                               + dx_ci(imin:imax) * dx_ci(imin:imax) * third ) * rhobar(imin:imax)
dmgrv(jr_min:jr_max)            = dmrst(jr_min:jr_max)

rstmss(1)                       = zero
grvmss(1)                       = zero
rstmss(jr_min:jr_max)           = rstmss(jr_min-1:jr_max-1) + dmrst(jr_min:jr_max)
grvmss(jr_min:jr_max)           = rstmss(jr_min:jr_max)

!-----------------------------------------------------------------------
!
!           \\\\\ COMPUTE NEWTONIAN SPHERICAL GRAVITY /////
!
!-----------------------------------------------------------------------

IF ( i_grav == 1 ) THEN

  gpotr_e                       = zero
  gpotr_c                       = zero
  gforcer_e                     = zero
  gforcer_c                     = zero
  gpot_dis_e                    = zero
  gpot_dis_c                    = zero

  IF ( irelhy == 0 ) THEN

    CALL grav_Newton( imin, imax, x_ei, x_ci, dx_ci, rhobar, gpot_e, gpot_c, &
&    gforce_e, gforce_c, gpot_dis_e, gpot_dis_c )
    agr                         = one
    agrh                        = one

!-----------------------------------------------------------------------
!
!             \\\\\ COMPUTE PN GR SPHERICAL GRAVITY /////
!
!-----------------------------------------------------------------------

  ELSE IF ( irelhy == 1 ) THEN

    CALL greff( imin, imax, x_ei, x_ci, dx_ci, rhobar, e_bar, e_nu_c_bar, &
&    f_nu_e_bar, p_bar, vx_bar, v2_bar, gpot_e, gpot_c, gforce_e, gforce_c, &
&    gpot_dis_e, gpot_dis_c )

!-----------------------------------------------------------------------
!  Reset neutrino energies for change in lapse, store in nu_dist_module
!  All radial rays have the same neutrino energies and lapse, so
!   calculate it once for j_ray = 1, ik_ray = 1 and use for all rays
!-----------------------------------------------------------------------

    ij_ray                      = 1
    ik_ray                      = 1
    IF ( ilapsehy == 0 ) THEN
      agr                       = one
      agrh                      = one
    ELSE IF ( ilapsehy == 1 ) THEN
      CALL agr_cal( jr_min, jr_max, ij_ray, ik_ray, rho, p, r, agr, agrh, nx )
      CALL agr_nu_cal( jr_min, jr_max )
      CALL enu_cal( jr_min, jr_max )
    ELSE
      WRITE (nlog,1001) ilapsehy
    END IF ! ilapsehy == 0

  ELSE
    WRITE (nlog,1003) irelhy
  END IF ! irelhy == 0

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    grav_x_c  (imin:imax  ,  ij_ray,ik_ray) = gforce_c  (imin:imax  )
    grav_x_e  (imin:imax+1,  ij_ray,ik_ray) = gforce_e  (imin:imax+1)
    grav_pot_c(imin:imax  ,  ij_ray,ik_ray) = gpot_c    (imin:imax  )
    grav_pot_e(imin:imax+1,  ij_ray,ik_ray) = gpot_e    (imin:imax+1)
    gtot_pot_c(imin:imax  ,  ij_ray,ik_ray) = gpot_dis_c(imin:imax  )
    gtot_pot_e(imin:imax+1,  ij_ray,ik_ray) = gpot_dis_e(imin:imax+1)
    agr_c     (imin:imax  ,  ij_ray,ik_ray) = agrh      (jr_min:jr_max)
    agr_e     (imin:imax+1,  ij_ray,ik_ray) = agr       (imin:imax+1)
    unu_c     (imin:imax  ,:,ij_ray,ik_ray) = unu       (jr_min:jr_max,:)
    unub_c    (imin:imax  ,:,ij_ray,ik_ray) = unub      (jr_min:jr_max,:)
    dunu_c    (imin:imax  ,:,ij_ray,ik_ray) = dunu      (jr_min:jr_max,:)
    unue_e    (imin:imax+1,:,ij_ray,ik_ray) = unue      (imin:imax+1,:)
    unube_e   (imin:imax+1,:,ij_ray,ik_ray) = unube     (imin:imax+1,:)
    dunue_e   (imin:imax+1,:,ij_ray,ik_ray) = dunue     (imin:imax+1,:)
  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

grav_pot_c_i                    = grav_pot_c

RETURN

ELSE IF ( i_grav == 2 ) THEN

  gforcer_e                     = zero
  gforcer_c                     = zero
  gforcen_e                     = zero
  gforcen_c                     = zero
  gpotr_e                       = zero
  gpotr_c                       = zero
  gpotn_e                       = zero
  gpotn_c                       = zero

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE NEWTONIAN NONSPHERICAL GRAVITY /////
!
!-----------------------------------------------------------------------

  IF ( irelhy == 0 ) THEN

    CALL grav_Newton( imin, imax, x_ei, x_ci, dx_ci, rhobar, gpot_e, &
&    gpot_c, gforce_e, gforce_c, gpot_dis_e, gpot_dis_c )

    mode                        = 0
    
    CALL poisson( mode, imin, imax, nx, ij_ray_dim, ik_ray_dim, ny, nz,  & 
&    x_ei, x_ci, dx_ci, y_ei, y_ci, dy_ci, z_ei, dz_ci, rho_c, grav_x_c, &
&    grav_y_c, grav_z_c, grav_pot_c, grav_x_e, grav_y_e, grav_z_e,       &
&    grav_pot_e )

    agr                         = one
    agrh                        = one

!-----------------------------------------------------------------------
!
!             \\\\\ COMPUTE PN NONSPHERICAL GRAVITY /////
!
!-----------------------------------------------------------------------

  ELSE IF ( irelhy == 1 ) THEN

    CALL grav_Newton( imin, imax, x_ei, x_ci, dx_ci, rhobar, gpotn_e, &
&    gpotn_c, gforcen_e, gforcen_c, gpot_dis_e, gpot_dis_c )

    CALL greff( imin, imax, x_ei, x_ci, dx_ci, rhobar, e_bar, e_nu_c_bar, &
&    f_nu_e_bar, p_bar, vx_bar, v2_bar, gpotr_e, gpotr_c, gforcer_e,      &
&    gforcer_c, gpot_dis_e, gpot_dis_c )

    mode                        = 0
    
    CALL poisson( mode, imin, imax, nx, ij_ray_dim, ik_ray_dim, ny, nz,  & 
&    x_ei, x_ci, dx_ci, y_ei, y_ci, dy_ci, z_ei, dz_ci, rho_c, grav_x_c, &
&    grav_y_c, grav_z_c, grav_pot_c, grav_x_e, grav_y_e, grav_z_e,       &
&    grav_pot_e )

!-----------------------------------------------------------------------
!  Reset neutrino energies for change in lapse, store in nu_dist_module
!  All radial rays have the same neutrino energies and lapse, so
!   calculate it once for j_ray = 1, ik_ray = 1 and use for all rays
!-----------------------------------------------------------------------

    ij_ray                      = 1
    ik_ray                      = 1
    IF ( ilapsehy == 0 ) THEN
      agr                       = one
      agrh                      = one
    ELSE IF ( ilapsehy == 1 ) THEN
      CALL agr_cal( jr_min, jr_max, ij_ray, ik_ray, rho, p, r, agr, agrh, nx )
      CALL agr_nu_cal( jr_min, jr_max )
      CALL enu_cal( jr_min, jr_max )
    ELSE
      WRITE (nlog,1001) ilapsehy
    END IF ! ilapsehy == 0

  ELSE
    WRITE (nlog,1001) irelhy
  END IF ! irelhy == 0

END IF ! i_grav == 1

!-----------------------------------------------------------------------
!
!         \\\\\ TRANSFER VARIABLES TO RADIAL_RAY_MODULE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             ||||| BEGIN LOOP OVER RADIAL ARRAYS |||||
!
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    grav_x_c  (imin:imax  ,  ij_ray,ik_ray) = grav_x_c  (imin:imax  ,ij_ray,ik_ray) &
&                                           + gforcer_c  (imin:imax  ) - gforcen_c  (imin:imax  )
    grav_x_e  (imin:imax+1,  ij_ray,ik_ray) = grav_x_e  (imin:imax+1,ij_ray,ik_ray) &
&                                           + gforcer_e  (imin:imax+1) - gforcen_e  (imin:imax+1)
    grav_pot_c(imin:imax  ,  ij_ray,ik_ray) = grav_pot_c(imin:imax  ,ij_ray,ik_ray) &
&                                           + gpotr_c    (imin:imax  ) - gpotn_c    (imin:imax  )
    grav_pot_e(imin:imax+1,  ij_ray,ik_ray) = grav_pot_e(imin:imax+1,ij_ray,ik_ray) &
&                                           + gpotr_e    (imin:imax+1) - gpotn_e    (imin:imax+1)
    gtot_pot_c(imin:imax  ,  ij_ray,ik_ray) = gpot_dis_c(imin:imax  )
    gtot_pot_e(imin:imax+1,  ij_ray,ik_ray) = gpot_dis_e(imin:imax+1)
    agr_c     (imin:imax  ,  ij_ray,ik_ray) = agrh      (jr_min:jr_max)
    agr_e     (imin:imax+1,  ij_ray,ik_ray) = agr       (imin:imax+1)
    unu_c     (imin:imax  ,:,ij_ray,ik_ray) = unu       (jr_min:jr_max,:)
    unub_c    (imin:imax  ,:,ij_ray,ik_ray) = unub      (jr_min:jr_max,:)
    dunu_c    (imin:imax  ,:,ij_ray,ik_ray) = dunu      (jr_min:jr_max,:)
    unue_e    (imin:imax+1,:,ij_ray,ik_ray) = unue      (imin:imax+1,:)
    unube_e   (imin:imax+1,:,ij_ray,ik_ray) = unube     (imin:imax+1,:)
    dunue_e   (imin:imax+1,:,ij_ray,ik_ray) = dunue     (imin:imax+1,:)
  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

grav_pot_c_i                    = grav_pot_c

!-----------------------------------------------------------------------
!
!              ||||| END LOOP OVER RADIAL ARRAYS |||||
!
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE radhyd_to_grav_i
