SUBROUTINE greff( imin, imax, r_e, r_c, dr_c, rhobar, e_bar, ei_nu, &
& flux_nu, p_gas, vx, v_square, gpot_e, gpot_c, gforce_e, gforce_c, &
& gpot_dis_e, gpot_dis_c )
!-----------------------------------------------------------------------
! 	gravitational Force/mass on the zone edge : gforce_e
!  	potential is computed on the zone edge points : gpot_e
!
!       calling subprogram mass_tov_cent to compute TOV mass
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: nx
USE numerical_module, ONLY : zero, half, frpi, frpith, third
USE physcnst_module, ONLY : g, cvel, ergmev, rmu

USE tov_potential_module, ONLY : effpot_e, effpot_c

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
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: e_bar       ! unshifted mean gas internal energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: ei_nu       ! unshifted mean neutrino internal energy (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: flux_nu     ! unshifted mean neutrino internal energy flux (ergs cm^{-2} s^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: p_gas       ! unshifted mean gas pressure (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: vx          ! unshifted mean gas x-velocity (cm s^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(in)  :: v_square    ! unshifted mean gas velocity squared (cm^{2} s^{-2})

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
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                         :: ii

REAL(KIND=double), PARAMETER                    :: c_sq     = cvel**2
REAL(KIND=double), PARAMETER                    :: c_frth   = cvel**4
REAL(KIND=double), PARAMETER                    :: c_sq_inv = 1.d0/cvel**2
REAL(KIND=double), PARAMETER                    :: tol      = 1.0d-6
REAL(KIND=double), PARAMETER                    :: frpi_inv = 1.0d0/frpi
REAL(KIND=double), PARAMETER                    :: UTOT0    = 8.9d0      ! change in the zero of energy (MeV)
REAL(KIND=double), PARAMETER                    :: ku       = ergmev/rmu ! ( # nucleons/gram )( erg/mev )
REAL(KIND=double), PARAMETER                    :: e_bar0   = ku * UTOT0

REAL(KIND=double), DIMENSION(nx)                :: e_gas       ! unshifted mean gas energy (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx)                :: p_nu        ! unshifted mean neutrino pressure (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx)                :: m_tov_e     ! unshifted zone_edged tov mass (g)
REAL(KIND=double), DIMENSION(nx)                :: m_tov_c     ! unshifted zone_centered tov mass (g)
REAL(KIND=double), DIMENSION(nx)                :: hen1        ! gas enthalphy+1
REAL(KIND=double), DIMENSION(nx)                :: dvol        ! radial zone volume (cm^{3})
REAL(KIND=double), DIMENSION(nx)                :: g_acc_c     ! zone-centered gravitational acceleration due to TOV mass
REAL(KIND=double), DIMENSION(nx)                :: g_acc_e     ! zone-edged gravitational acceleration due to TOV mass
REAL(KIND=double), DIMENSION(nx)                :: Gam2_inv    ! gamma factor^(-2)
REAL(KIND=double), DIMENSION(nx)                :: m_grav      ! effective gravitational mass
REAL(KIND=double)                               :: mas0        ! mass at the outermost zone edge point
REAL(KIND=double)                               :: phi0        ! effective potential at the outer edge defined as -G*mas0/R_out

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

e_gas                     = zero
p_nu                      = zero
m_tov_e                   = zero
m_tov_c                   = zero
hen1                      = zero
dvol                      = zero
g_acc_c                   = zero
g_acc_e                   = zero
Gam2_inv                  = zero
m_grav                    = zero

!-----------------------------------------------------------------------
!  Preprocessing inputs: gas energy & neutrino pressure
!-----------------------------------------------------------------------

e_gas(1:imax)             = ( e_bar(1:imax) - e_bar0 ) * rhobar(1:imax)
p_nu (1:imax)             = third * ei_nu(1:imax)
dvol (1:imax)             = frpi * dr_c(1:imax) * ( r_e(1:imax) * ( r_e(1:imax) + dr_c(1:imax) ) &
&                         + dr_c(1:imax) * dr_c(1:imax) * third ) 

DO ii = imin,imax
  IF ( rhobar(ii)/rhobar(imin) < tol ) THEN
    hen1(ii)              = 1.0d0
  ELSE
    hen1(ii)              = 1.0d0 + ( e_gas(ii) + p_gas(ii) )/( rhobar(ii) * c_sq )
  END IF
END DO

!-----------------------------------------------------------------------
!  Compute TOV mass at the zone-edges and zone centers
!-----------------------------------------------------------------------

CALL mass_tov( imin, imax, r_c, r_e, dr_c, rhobar, e_gas, ei_nu, flux_nu, &
& vx, v_square, m_tov_e, m_tov_c )

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE POTENTIAL AT THE ZONE EDGE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Boundary condition of potential
!-----------------------------------------------------------------------

mas0                      = m_tov_e(imax+1)
phi0                      = -g * mas0/r_e(imax+1)

!-----------------------------------------------------------------------
!  Compute gamma inverse
!-----------------------------------------------------------------------

Gam2_inv(imin:imax)       = 1.0d0/( 1.0d0 + v_square(imin:imax) * c_sq_inv &
&                         - 2.0d0 * g * m_tov_c(imin:imax)/( c_sq * r_c(imin:imax) ) )

!-----------------------------------------------------------------------
!  Compute potential at the zone edge
!-----------------------------------------------------------------------

g_acc_c(imin+1:imax)      = half * ( m_tov_e(imin+1:imax  )/( r_e(imin+1:imax  ) * r_e(imin+1:imax  ) )      &
&                         +          m_tov_e(imin+2:imax+1)/( r_e(imin+2:imax+1) * r_e(imin+2:imax+1) ) )
g_acc_c(imin)             = half * m_tov_e(imin+1)/( r_e(imin+1) * r_e(imin+1) )

gpot_e(imax+1)            = phi0

DO ii = imax,imin,-1

  gpot_e(ii)              = gpot_e(ii+1) - g * ( g_acc_c(ii) * dr_c(ii) + dvol(ii)                           &
&                         * ( ( p_gas(ii) + p_nu(ii) )/( c_sq * r_c(ii) ) ) ) * Gam2_inv(ii) * hen1(ii)

END DO

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE POTENTIAL AT THE ZONE CENTERS /////
!
!-----------------------------------------------------------------------

g_acc_e(imin)             = zero
g_acc_e(imin+1:imax+1)    = m_tov_e(imin+1:imax+1)/( r_e(imin+1:imax+1) * r_e(imin+1:imax+1) )


!-----------------------------------------------------------------------
!  Compute zone center potential at imax
!-----------------------------------------------------------------------

gpot_c(imax)              = phi0 - g * 5.0d-1 * ( g_acc_e(imax) * dr_c(imax) + dvol(imax)                     &
&                         * ( ( p_gas(imax) + p_nu(imax) )/( c_sq * r_c(imax) ) ) ) * Gam2_inv(imax) * hen1(imax) 

!-----------------------------------------------------------------------
!  Integrate to compute zone center potential
!-----------------------------------------------------------------------

DO ii = imax-1,imin,-1

!...>>First half

  gpot_c(ii)              = gpot_c(ii+1) - g * 5.0d-1 * ( g_acc_e(ii+1) * dr_c(ii+1) + dvol(ii+1)             &
&                         * ( ( p_gas(ii+1) + p_nu(ii+1) )/( c_sq * half * ( r_e(ii+1) + r_c(ii+1) ) ) ) )    &
&                         * Gam2_inv(ii+1) * hen1(ii+1)

!....>>Second half

  gpot_c(ii)              = gpot_c(ii  ) - g * 5.0d-1 * ( g_acc_e(ii+1) * dr_c(ii  ) + dvol(ii  )             &
&                         * ( ( p_gas(ii  ) + p_nu(ii  ) )/( c_sq * half * ( r_e(ii+1) + r_c(ii)   ) ) ) )    &
&                         * Gam2_inv(ii  ) * hen1(ii  )

END DO

!-----------------------------------------------------------------------
!
!           \\\\\ COMPUTE FINAL POTENTIALS AND FORCES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute force at the zone edge from potential at zone centers
!-----------------------------------------------------------------------

gforce_e(imin)            = 0.0d0
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

gpot_c(imin:imax)         = 0.5d0 * ( gpot_e(imin:imax) + gpot_e(imin+1:imax+1) )

!-----------------------------------------------------------------------
!
!      \\\\\ COMPUTE DISASSEMBLY POTENTIAL AT THE ZONE EDGE /////
!
!-----------------------------------------------------------------------  

IF ( r_e(imin) == zero ) THEN
  gpot_dis_e(imin)        = zero
ELSE
  gpot_dis_e(imin)        = - g * m_tov_e(imin)/r_e(imin)
END IF

m_grav(imin)              = m_tov_c(imin)
DO ii = imin,imax
  m_grav(ii+1)            = m_grav(ii) + m_tov_c(ii+1) - m_tov_c(ii) + dvol(ii)                           &
&                         * ( ( p_gas(ii) + p_nu(ii) )/c_sq ) * Gam2_inv(ii) * hen1(ii)
END DO


gpot_dis_e(imin+1:imax+1) = - g * m_grav(imin+1:imax+1)/r_e(imin+1:imax+1)

gpot_dis_c(imin:imax)     = half * ( gpot_dis_e(imin:imax) + gpot_dis_e(imin+1:imax+1) )

!-----------------------------------------------------------------------
!  Load tov_potential_module
!-----------------------------------------------------------------------

effpot_e(imin:imax+1)     = gpot_e(imin:imax+1)
effpot_c(imin:imax)       = gpot_c(imin:imax)

END Subroutine greff
