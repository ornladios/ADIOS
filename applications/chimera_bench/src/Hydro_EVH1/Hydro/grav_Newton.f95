SUBROUTINE grav_Newton( imin, imax, r_e, r_c, dr_c, rhobar, gpot_e, gpot_c, &
& gforce_e, gforce_c, gpot_dis_e, gpot_dis_c )
!-----------------------------------------------------------------------
!
!    File:         grav_Newton
!    Module:       grav_Newton
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To calculate the Mewtonian gravitational potentials and forces at
!       the zone edges and centers assuming spherical symmetry.
!
!    Input arguments:
!
!  imin         : minimum x-array index
!  imax         : maximum x-array index
!  nx           : x-array extent
!  jmin         : minimum yx-array index
!  jmax         : maximum y-array index
!  nx           : y-array extent
!  x_c          : radial midpoint of zone (cm)
!  mass_a_ray   : mass/angular ray (g)
!
!    Output arguments:
!
!  grav_pot     : gravitational potential
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: nx
USE numerical_module, ONLY : zero, half, frpi, third, frpith
USE physcnst_module, ONLY : g
     
IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER                                         :: imin        ! minimum unshifted radial index
INTEGER                                         :: imax        ! maximum unshifted radial index

REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: r_c         ! unshifted zone-centered radius (cm)
REAL(KIND=double), DIMENSION(nx+1), INTENT(in)  :: r_e         ! unshifted zone-edged radius (cm)
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: dr_c        ! unshifted zone-thicknesses (cm)
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: rhobar      ! unshifted mean density (g cm^{-3})

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: gpot_e      ! gravitational potential defined at zone-edge (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: gpot_c      ! gravitational potential defined at zone-center (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: gforce_e    ! grav. force defined at zone-centered (dynes g^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: gforce_c    ! grav. force defined at zone-centered (dynes g^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: gpot_dis_e  ! gravitational disassembly potential defined at zone-edge (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: gpot_dis_c  ! gravitational disassembly potential defined at zone-center (ergs g^{-1})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                         :: ii

REAL(KIND=double), DIMENSION(nx)                :: m_neut_e    ! unshifted zone_edged Newtonian mass (g)
REAL(KIND=double), DIMENSION(nx)                :: m_neut_c    ! unshifted zone_centered Newtonian mass (g)
REAL(KIND=double), DIMENSION(nx)                :: dvol        ! radial zone volume (cm^{3})
REAL(KIND=double), DIMENSION(nx)                :: g_acc_c     ! zone-centered gravitational acceleration due to TOV mass
REAL(KIND=double), DIMENSION(nx)                :: g_acc_e     ! zone-edged gravitational acceleration due to TOV mass
REAL(KIND=double)                               :: mas0        ! mass at the outermost zone edge point
REAL(KIND=double)                               :: phi0        ! effective potential at the outer edge defined as -G*mas0/R_out

!-----------------------------------------------------------------------
!  Radial zone volumes
!-----------------------------------------------------------------------

dvol (1:imax)             = frpi * dr_c(1:imax) * ( r_e(1:imax) * ( r_e(1:imax) + dr_c(1:imax) ) &
&                         + dr_c(1:imax) * dr_c(1:imax) * third ) 

!-----------------------------------------------------------------------
!  Radial zone masses
!-----------------------------------------------------------------------

ii                        = imin
IF ( r_e(ii) == zero ) THEN
  m_neut_e(ii)            = zero
ELSE
  m_neut_e(ii)            = rhobar(ii) * frpith * r_e(ii)**3
END IF

DO ii = imin,imax
  m_neut_e(ii+1)          = m_neut_e(ii) + rhobar(ii) * dvol(ii)
END  DO

m_neut_c(imin:imax)       = half * ( m_neut_e(imin:imax ) + m_neut_e(imin+1:imax+1) )

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE POTENTIAL AT THE ZONE EDGE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Boundary condition of potential
!-----------------------------------------------------------------------

mas0                      = m_neut_e(imax+1)
phi0                      = -g * mas0/r_e(imax+1)

!-----------------------------------------------------------------------
!  Compute potential at the zone edge
!-----------------------------------------------------------------------

g_acc_c(imin+1:imax)      = half * ( m_neut_e(imin+1:imax  )/( r_e(imin+1:imax  ) * r_e(imin+1:imax  ) )      &
&                         +          m_neut_e(imin+2:imax+1)/( r_e(imin+2:imax+1) * r_e(imin+2:imax+1) ) )
g_acc_c(imin)             = half * m_neut_e(imin+1)/( r_e(imin+1) * r_e(imin+1) )

gpot_e(imax+1)            = phi0

DO ii = imax,imin,-1
  gpot_e(ii)              = gpot_e(ii+1) - g * g_acc_c(ii) * dr_c(ii) 
END DO

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE POTENTIAL AT THE ZONE CENTERS /////
!
!-----------------------------------------------------------------------

g_acc_e(imin)             = zero
g_acc_e(imin+1:imax+1)    = m_neut_e(imin+1:imax+1)/( r_e(imin+1:imax+1) * r_e(imin+1:imax+1) )

!-----------------------------------------------------------------------
!  Integrate to compute zone center potential
!-----------------------------------------------------------------------

DO ii = imax-1,imin,-1

!...>>First half

  gpot_c(ii)              = gpot_c(ii+1) - g * half * g_acc_e(ii+1) * dr_c(ii+1) 

!....>>Second half

  gpot_c(ii)              = gpot_c(ii  ) - g * half * g_acc_e(ii+1) * dr_c(ii  )

END DO

!-----------------------------------------------------------------------
!
!           \\\\\ COMPUTE FINAL POTENTIALS AND FORCES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute force at the zone edge from potential at zone centers
!-----------------------------------------------------------------------

gforce_e(imin)            = zero
gforce_e(imin+1:imax)     = - ( gpot_c(imin+1:imax) - gpot_c(imin:imax-1))/( r_c(imin+1:imax) - r_c(imin:imax-1) )
gforce_e(imax+1)          = - g * mas0/r_e(imax+1)**2

!-----------------------------------------------------------------------
!  Compute force at the zone centers from potential at zone edges
!-----------------------------------------------------------------------

gforce_c(imin:imax)       = - ( gpot_e(imin+1:imax+1) - gpot_e(imin:imax) )/( r_e(imin+1:imax+1) - r_e(imin:imax) )

!-----------------------------------------------------------------------
!  Compute potential at the center from zone-edge average for better
!   consistency
!-----------------------------------------------------------------------

gpot_c(imin:imax)         = half * ( gpot_e(imin:imax) + gpot_e(imin+1:imax+1) )

!-----------------------------------------------------------------------
!
!      \\\\\ COMPUTE DISASSEMBLY POTENTIAL AT THE ZONE EDGE /////
!
!-----------------------------------------------------------------------

IF ( r_e(imin) == zero ) THEN
  gpot_dis_e(imin)        = zero
ELSE
  gpot_dis_e(imin)        = - g * m_neut_e(imin)/r_e(imin)
END IF

gpot_dis_e(imin+1:imax+1) = - g * m_neut_e(imin+1:imax+1)/r_e(imin+1:imax+1)

gpot_dis_c(imin:imax)     = half * ( gpot_dis_e(imin:imax) + gpot_dis_e(imin+1:imax+1) )

RETURN
END SUBROUTINE grav_Newton
