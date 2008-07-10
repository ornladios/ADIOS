SUBROUTINE g_brem( y_p, eta_star_p, g_fit )
!-----------------------------------------------------------------------
!
!    File:         g_brem
!    Module:       g_brem
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/23/02
!
!    Purpose:
!      To compute g_fit-component of the kernal for neutrino bremss-
!       trahlung and inelastic neutrino scattering in a nucleon
!       field.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  y           : m_{pi}^{2}/M_{N}c^{2}kT
!  eta_star    : neutron degeneracy parameter
!
!    Output arguments:
!  g           : dimensionless quantity related to S_{sigma}
!
!    Modules:
!  kind_module, numerical_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : one
      
IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(IN)  :: y_p                  ! pion mass parameter
REAL(KIND=double), INTENT(IN)  :: eta_star_p           ! degeneracy parameter

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT) :: g_fit                ! dimensionless fitting parameter

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)              :: y                    ! use in place of y_p
REAL(KIND=double)              :: eta_star             ! use in place of eta_star_p

REAL(KIND=double), PARAMETER   :: y_min = 1.d-10
REAL(KIND=double), PARAMETER   :: eta_min = 1.d-10
REAL(KIND=double)              :: y2                   ! y*y
REAL(KIND=double)              :: eta_star_1           ! 1/eta_star
REAL(KIND=double)              :: alpha1_denom
REAL(KIND=double)              :: alpha1,alpha2,alpha3
REAL(KIND=double)              :: p1,p2
REAL(KIND=double)              :: fexp                 ! declare function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Replace input variables for modification
!-----------------------------------------------------------------------

y                   = y_p
eta_star            = eta_star_p

!-----------------------------------------------------------------------
!  Prevent singular behavior of terms
!-----------------------------------------------------------------------

y                  = MAX( y       , y_min   )
eta_star           = MAX( eta_star, eta_min )

!-----------------------------------------------------------------------
!  Compute alpha1
!-----------------------------------------------------------------------

y2                 = y * y
eta_star_1         = one/eta_star
alpha1_denom       = 25.d0 * y2 + one
alpha1             = ( 0.5d0 + eta_star_1 )/( one + eta_star_1 )   &
&                  * one/alpha1_denom                              &
&                  + ( 0.5d0 + eta_star/15.6d0 ) * 25.d0 * y2      &
&                 / alpha1_denom

!-----------------------------------------------------------------------
!  Compute alpha2
!-----------------------------------------------------------------------

alpha2             = ( 0.63d0 + 0.04d0 * eta_star**1.45d0 )        &
&                  / ( one + 0.02d0 * eta_star**2.5d0 )

!-----------------------------------------------------------------------
!  Compute alpha3
!-----------------------------------------------------------------------

alpha3             = 1.2d0 * fexp( 0.6d0 * eta_star                &
&                  - 0.4d0 * eta_star**1.5d0 )

!-----------------------------------------------------------------------
!  Compute p1
!-----------------------------------------------------------------------

p1                 = ( 1.8d0 + 0.45d0 * eta_star )                 &
&                  / ( one + 0.15d0 * eta_star**1.5d0 )

!-----------------------------------------------------------------------
!  Compute p2
!-----------------------------------------------------------------------

p2                 = 2.3d0 - 0.05 * eta_star                       &
&                  / ( one + 0.025d0 * eta_star )

!-----------------------------------------------------------------------
!  Compute g_fit
!-----------------------------------------------------------------------

g_fit              = ( alpha1 + alpha2 * y**p1 )                   &
&                  / ( one + alpha3 * y**p2                        &
&                  + alpha2 * y**( p1 + 2.d0 )/13.75d0 )


RETURN
END SUBROUTINE g_brem
