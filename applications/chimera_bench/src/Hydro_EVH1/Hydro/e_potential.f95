SUBROUTINE  e_potential( xf, r )
!=======================================================================  
!  Calculate gravitational potential energy.
!=======================================================================

USE kind_module, ONLY: double
USE array_module, ONLY: max_12
USE numerical_module, ONLY : zero, half, frpi, third, frpith
USE physcnst_module, ONLY: g

USE edit_module, ONLY: nlog
USE evh1_zone, ONLY: imax
USE evh1_sweep, ONLY: egrav
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in), DIMENSION(max_12) :: xf       ! zone positions (cm)
REAL(KIND=double), INTENT(in), DIMENSION(max_12) :: r        ! zone densities (g cm^{-3})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                :: n           ! padded zone index
INTEGER                                :: nmin        ! padded zone index minimum
INTEGER                                :: nmax        ! padded zone index maximum

REAL(KIND=double), DIMENSION(max_12)   :: dvol        ! padded radial zone volume (cm^{3})
REAL(KIND=double), DIMENSION(max_12)   :: dxf         ! padded radial zone thickness (cm)
REAL(KIND=double), DIMENSION(max_12)   :: m_neut_e    ! padded zone_edged Newtonian mass (g)
REAL(KIND=double), DIMENSION(max_12)   :: egrav_e     ! padded zone_edged gravitational potential (ergs g^{-1})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nmin                      = 7
nmax                      = imax + 6

!-----------------------------------------------------------------------
!  Radial zone thickness and volumes
!-----------------------------------------------------------------------

dxf(nmin:nmax)            = xf(nmin+1:nmax+1) - xf(nmin:nmax)
dvol(nmin:nmax)           = frpi * dxf(nmin:nmax) * ( xf(nmin:nmax) * ( xf(nmin+1:nmax+1) ) &
&                         + dxf(nmin:nmax) * dxf(nmin:nmax) * third ) 

!-----------------------------------------------------------------------
!  Radial zone masses
!-----------------------------------------------------------------------

IF ( xf(nmin) == zero ) THEN
  m_neut_e(nmin)          = zero
ELSE
  m_neut_e(nmin)          = r(nmin) * frpith * xf(nmin)**3
END IF

DO n = nmin,nmax
  m_neut_e(n+1)           = m_neut_e(n) + r(n) * dvol(n)
END  DO

!-----------------------------------------------------------------------
!  Gravitational potential
!-----------------------------------------------------------------------

IF ( xf(nmin) == zero ) THEN
  egrav_e(nmin)           = zero
ELSE
  egrav_e(nmin)           = - g * m_neut_e(nmin)/xf(nmin)
END IF

egrav_e(nmin+1:nmax+1)    = - g * m_neut_e(nmin+1:nmax+1)/xf(nmin+1:nmax+1)

egrav(nmin:nmax)          = half * ( egrav_e(nmin:nmax) + egrav_e(nmin+1:nmax+1) )

RETURN
END SUBROUTINE e_potential
