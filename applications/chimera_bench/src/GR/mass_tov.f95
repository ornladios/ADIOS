
SUBROUTINE mass_tov( imin, imax, r_c, r_e, dr_c, rhobar, e_gas, ei_nu, & 
& flux_nu_e, vx, v_square, m_tov_e, m_tov_c )
!-----------------------------------------------------------------------
! Subroutine mass_tov_cent
!    When provided energy density, flux and v_square, etc., 
!    this computes 'TOV mass' enclosed inside each zone-centered radius
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: nx
USE numerical_module, ONLY : zero, half, frpi, frpith, third

USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid
USE physcnst_module, ONLY : g, cvel

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER                                        :: imin         ! minimum unshifted radial index
INTEGER                                        :: imax         ! maximum unshifted radial index

REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: r_c          ! unshifted zone-centered radius (cm)
REAL(KIND=double), DIMENSION(nx+1), INTENT(in) :: r_e          ! unshifted zone-edged radius (cm)
REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: dr_c         ! unshifted zone-thicknesses (cm)
REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: rhobar       ! unshifted mean density (g cm^{-3})
REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: e_gas        ! unshifted mean gas energy (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: ei_nu        ! unshifted mean neutrino internal energy (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: flux_nu_e    ! unshifted mean zone-edged neutrino internal energy flux (ergs cm^{-2} s^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: vx           ! unshifted mean gas x-velocity (cm s^{-1})
REAL(KIND=double), DIMENSION(nx),   INTENT(in) :: v_square     ! unshifted mean gas velocity squared (cm^{2} s^{-2})

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: m_tov_e     ! TOV mass enclosed in each zone-edge radius (g)
REAL(KIND=double), DIMENSION(nx),   INTENT(out) :: m_tov_c     ! TOV mass enclosed in each zone-centered radius (g)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, SAVE                                  :: first = .true.

INTEGER                                        :: ii

REAL(KIND=double), PARAMETER                   :: c_sq     = cvel**2
REAL(KIND=double), PARAMETER                   :: cvel_inv = 1.d0/cvel
REAL(KIND=double), PARAMETER                   :: tol      = 1.0d-6

REAL(KIND=double), DIMENSION(nx)               :: flux_nu_c    ! unshifted mean zone-centered neutrino internal energy flux (ergs cm^{-2} s^{-1})
REAL(KIND=double), DIMENSION(nx)               :: dvol         ! radial zone volume (cm^{3})
REAL(KIND=double), DIMENSION(nx)               :: mas_tmp_e    ! temporary value of zone-edged TOV mass (g)
REAL(KIND=double), DIMENSION(2000), SAVE       :: mas_tov_prev ! previous value of zone-edged TOV mass (g)
REAL(KIND=double)                              :: Gam          ! zone-centered gamma factor * cvel
REAL(KIND=double)                              :: relax        ! relaxation factor
REAL(KIND=double)                              :: error        ! relative error

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT ('!! Too Relativistic in m_tov_c !!')
 1003 FORMAT ('!! Too Relativistic in m_tov_c !!, myid=',i6)

!-----------------------------------------------------------------------
!  0) pre-process input data  : nu-flux has to be averaged at the zone-center
!	Simple geometrical linear interpolation in r
!    ## this averaging scheme may be rather crude
!-----------------------------------------------------------------------

DO ii = imin,imax
  flux_nu_c(ii)           = ( r_e(ii+1) - r_c(ii) )/dr_c(ii) * flux_nu_e(ii)   &
&                         + ( r_c(ii)   - r_e(ii) )/dr_c(ii) * flux_nu_e(ii+1)
END DO
dvol (1:imax)             = frpi * dr_c(1:imax) * ( r_e(1:imax) * ( r_e(1:imax) + dr_c(1:imax) ) &
&                         + dr_c(1:imax) * dr_c(1:imax) * third ) 

!-----------------------------------------------------------------------
! 1) give an initial guess of the zone-edged mass distribution
!-----------------------------------------------------------------------

IF ( first ) THEN
  first                   = .false.
  ii                      = imin
  IF ( r_e(ii) == zero ) THEN
    m_tov_e(ii)           = zero
  ELSE
    m_tov_e(ii)           = rhobar(ii) * frpith * r_e(ii)**3
  END IF

  DO ii = 1,imax
    m_tov_e(ii+1)         = m_tov_e(ii) + rhobar(ii) * dvol(ii)
  END  DO
ELSE ! first = false
  m_tov_e(imin:imax+1)    = mas_tov_prev(imin:imax+1)
END IF ! first

m_tov_c(imin:imax)        = half * ( m_tov_e(imin:imax ) + m_tov_e(imin+1:imax+1) )

!-----------------------------------------------------------------------
!                  ||||| Iteration starts here |||||
!-----------------------------------------------------------------------

error                     = 1.0d0

iteration: DO 

  error                   = 0.0d0 ! reset error

!------------------------------------------------------------------------------------
! 2) integrate density integrand
!    Gam is computed with the preceding value of m_tov_c and is used to compute the
!     new value of m_tov_c
!------------------------------------------------------------------------------------

  ii = imin
  IF ( r_e(ii) == zero ) THEN
    mas_tmp_e(ii)         = zero
  ELSE
    mas_tmp_e(ii)         = rhobar(ii) * frpith * r_e(ii)**3
  END IF

  DO ii = imin,imax
    Gam                   = DSQRT( c_sq + v_square(ii) - 2.0d0 * g * m_tov_c(ii)/r_c(ii) ) * cvel_inv

    mas_tmp_e(ii+1)       = mas_tmp_e(ii)                                            &
&                         + Gam * ( rhobar(ii) + ( e_gas(ii) + ei_nu(ii)             &
&                         + vx(ii) * flux_nu_c(ii)/( c_sq * Gam ) )/c_sq )           &
&                         * dvol(ii)

    error                 = error + DABS( mas_tmp_e(ii+1) - m_tov_e(ii+1) ) /DABS(m_tov_e(ii+1) )

    IF ( ( 2.0d0 * g * mas_tmp_e(ii+1)/( c_sq * r_e(ii+1)) - 1.0d0 ) > - tol ) THEN
      WRITE (nlog,1001)
      WRITE (nprint,1003) myid
      STOP    
    END IF ! 2.0d0 * g * mas_tmp_e(ii+1)/( c_sq * r_e(ii+1)) - 1.0d0 ) > - tol

  END DO

!   --- relaxation factor for iteration

  IF ( error > 1.0d0 ) THEN
    relax                 = 0.1d0
  ELSE IF ( error <= 1.0d0 ) THEN
    relax                 = 1.0d0
  END IF


! -------------------   revised  mass distribution

  m_tov_e(imin+1:imax+1)  = relax * mas_tmp_e(imin+1:imax+1) + ( 1.0d0 - relax) * m_tov_e(imin+1:imax+1)
  m_tov_c(imin:imax)      = half * ( m_tov_e(imin:imax ) + m_tov_e(imin+1:imax+1) )
		
  If ( error < tol ) EXIT  ! if error is small enough, exit iteration do loop

!-----------------------------------------------------------------------
!                   ||||| Iteration ends here |||||
!-----------------------------------------------------------------------

END DO iteration

mas_tov_prev(imin:imax+1) = m_tov_e(imin:imax+1)

RETURN
END SUBROUTINE mass_tov
