SUBROUTINE setup_restart_GR( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nx, ny, nz, nez, nnu, rhop, rp, drp, agr_e, agr_c )
!-----------------------------------------------------------------------
!
!    File:         setup_restart_GR
!    Module:       setup_restart_GR
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/06
!
!    Purpose:
!      To compute, from the initial problem configuration, the GR quantities 
!       needed to begin the execution of the simulation.
!
!    Subprograms called:
!  radhyd_to_grav_i : computes the GR variables in the post-Newtinian approximation
!  setup_rel        : computes the GR variables in the fully relativistic case
!    
!    Input arguments:
!  imin             : inner x-zone array index
!  imax             : outer x-zone array index
!  ij_ray           : j-index of a radial ray
!  ik_ray           : k-index of a radial ray
!  ij_ray_dim       : number of y-zones on a processor before swapping
!  ik_ray_dim       : number of z-zones on a processor before swapping
!  nx               : x_array extent
!  ny               : y_array extent
!  nz               : z_array extent
!  nez              : neutrino energy array extent
!  nnu              : neutrino flavor array extent
!  rhop             : density (cm^{-3})
!  rp               : radial zone radii (cm)
!  drp              : radial zone thickness (cm)
!
!    Output arguments:
!  agr_e            : unsifted zone-edged lapse function
!  agr_c            : unsifted zone-centered lapse function
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, mdl_cnfg_module, prb_cntl_module
!     
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, frpith

USE edit_module, ONLY : nprint, nlog
USE mdl_cnfg_module, ONLY : jr_min, jr_max, r, gamgr, gamgrr, grvmss, &
& rstmss, agr, dmgrv, dmrst, agrr, agrh, rho, dr, wgr, wgrr
USE prb_cntl_module, ONLY : irelhy

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin            ! inner x_array index
INTEGER, INTENT(in)              :: imax            ! outer x_array index

INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z-zones on a processor before swapping with z

INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: ny              ! y-array extent
INTEGER, INTENT(in)              :: nz              ! s-array extent
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)             :: rp           ! radial zone radii (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)               :: drp          ! radial zone thickness (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: rhop ! density (g cm^{-3})

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)  :: agr_e       ! unsifted zone-edged lapse function
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: agr_c       ! unsifted zone-centered lapse function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! iteration index
INTEGER, PARAMETER               :: it_max = 5      ! number of iterations to get PN gravity
INTEGER                          :: j               ! radial zone index
INTEGER                          :: jr_maxp         ! jr_max + 1

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Variable initialization complete in setup_restart_GR')
  103 FORMAT (' Newtonian rest masses computed in setup_restart_GR')
  105 FORMAT (' GR variables computed in the Newtonian approximation')
  107 FORMAT (' GR variables computed in the Post-Newtonian approximation')
  109 FORMAT (' GR variables computed in the fully relativistic case')

  201 FORMAT (' irelhy =',i5,' not supported in setup_restart_GR')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize shifted radial array index boundaries
!-----------------------------------------------------------------------

jr_min                           = imin + 1
jr_max                           = imax + 1
jr_maxp                          = jr_max + 1

!-----------------------------------------------------------------------
!  Set quatities at inner edge of configuration
!-----------------------------------------------------------------------

r(1)                             = zero
gamgr(1)                         = one
gamgrr(1)                        = one
grvmss(1)                        = zero
rstmss(1)                        = zero

IF ( agr(1) == zero ) agr(1) = agr(2)
IF ( agr_e(1,ij_ray,ik_ray) == zero ) agr_e(1,ij_ray,ik_ray) = agr_e(2,ij_ray,ik_ray)

IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,101)

!-----------------------------------------------------------------------
!
!                  \\\\\ TRANSFER VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to shifted mgfld arrays
!-----------------------------------------------------------------------

rho   (imin+1:imax+1)            = rhop   (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer zone-edgeed independent variables to shifted mgfld arrays
!-----------------------------------------------------------------------

r(imin:imax+1)                   = rp(imin:imax+1)
r(imax-imin+3)                   = r(imax-imin+2) + ( r(imax-imin+2) - r(imax-imin+1) )

!-----------------------------------------------------------------------
!  Derived zone-centered and zone-edged variables
!-----------------------------------------------------------------------

!........radial zone thickness

dr(imin+1:imax+1)                = drp(imin:imax)

!........Newtonian rest masses

IF ( r(1) == zero ) THEN
  dmrst(1)                       = zero
ELSE
  dmrst(1)                       = frpith * r(1)**3 * rho(2)
END IF

rstmss(1)                        = dmrst(1)

DO j = jr_min,jr_max
  dmrst(j)                       = frpith * ( r(j) * ( r(j) + r(j-1) ) + r(j-1) * r(j-1) ) * dr(j) * rho(j)
  rstmss(j)                      = rstmss(j-1) + dmrst(j)
END DO

IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,103)

!-----------------------------------------------------------------------
!
!                     \\\\\ GR VARIABLES ////
!
!-----------------------------------------------------------------------

GR: SELECT CASE (irelhy)

!-----------------------------------------------------------------------
!  irelhy = 0 : Nonrelativistic hydro and transport
!-----------------------------------------------------------------------

 CASE(0) GR

  wgr   (1:jr_max)               = one
  wgrr  (1:jr_max)               = one
  gamgr (1:jr_max)               = one
  gamgrr(1:jr_max)               = one
  agr   (1:jr_max)               = one
  agrr  (1:jr_max)               = one
  agr_e (1:imax+1,ij_ray,ik_ray) = one
  agr_c (1:imax  ,ij_ray,ik_ray) = one
  dmgrv (1:jr_max)               = dmrst (1:jr_max)
  grvmss(1:jr_max)               = rstmss(1:jr_max)

  agr  (jr_max+1)                = agr(jr_max)
  agr_e(imax+2,ij_ray,ik_ray)    = agr_e(imax+1,ij_ray,ik_ray)
  agrr (jr_max+1)                = agrr(jr_max)
  agrh (imax+1)                  = agrh(imax)
  agr_c(imax+1,ij_ray,ik_ray)    = agr_c(imax,ij_ray,ik_ray)

  CALL radhyd_to_grav_i( nx, ny, nz, ij_ray_dim, ik_ray_dim )

  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim )  WRITE (nlog,105)

  RETURN

 CASE(1) GR

!-----------------------------------------------------------------------
!  irelhy = 1 : post-newtonian
!-----------------------------------------------------------------------

  CALL agr_nu_cal( jr_min, jr_max )
  CALL enu_cal( jr_min, jr_max )
  DO i = 1,it_max
    CALL radhyd_to_nu_energy_flux( ij_ray_dim, ik_ray_dim, nx, nez, nnu )
    CALL angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )
    CALL radhyd_to_grav_i( nx, ny, nz, ij_ray_dim, ik_ray_dim )
  END DO ! i = 1,it_max

  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim )  WRITE (nlog,107)

  RETURN

!-----------------------------------------------------------------------
!  irelhy = 2: Relativistic hydro and transport
!-----------------------------------------------------------------------

 CASE(2) GR

  CALL setup_rel( jr_min, jr_max, ij_ray, ik_ray )

  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim )  WRITE (nlog,109)

  RETURN

!-----------------------------------------------------------------------
!  Stop if irelhy /= 0 and irelhy /= 1
!-----------------------------------------------------------------------

 CASE DEFAULT GR

  WRITE (nprint, 201) irelhy
  WRITE (nlog, 201) irelhy
  STOP

END SELECT GR

RETURN
END SUBROUTINE setup_restart_GR
