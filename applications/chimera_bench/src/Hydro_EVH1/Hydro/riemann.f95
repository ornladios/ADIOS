SUBROUTINE riemann ( lmin, lmax, game, gamc, prgh, urgh, vrgh, gergh, &
& gcrgh,  plft, ulft, vlft, gelft, gclft, pmid, umid, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         riemann
!    Module:       riemann
!    Type:         Subprogram
!    Author:       The Virginia Numerical Bull Session ideal hydrodynamics PPMLR
!    Modifier:     S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To solve the Riemann shock tube problem for the left and right
!       input states, using the Newton interation procedure described
!       in van Leer (1979).
!
!    Input arguments:
!
!  lmin  : minimum padded array index
!  lmax  : maximum padded array index
!  game  : 1 + p/(rho*ei)
!  gamc  : adiabatic gamma_1
!  prgh  : mean pressure over the right causal domain
!  urgh  : mean velocity over the right causal domain
!  vrgh  : mean pecific volume (initially density) over the right causal domain
!  gergh : mean 1 + p/(rho*ei) over the right causal domain
!  gcrgh : mean adiabatic gamma_1 over the right causal domain
!  plft  : mean pressure over the right causal domain
!  ulft  : mean velocity over the right causal domain
!  vlft  : mean pecific volume (initially density) over the right causal domain
!  gelft : mean 1 + p/(rho*ei) over the right causal domain
!  gclft : mean adiabatic gamma_1 over the right causal domain
!  j_ray : index denoting a specific radial or angular ray
!
!    Output arguments:
!
!  pmid  : time averaged pressure at the zone interface
!  umid  : time averaged velocity at the zone interface
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: max_12

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                  :: lmin     ! zone number of first physical zone
INTEGER, INTENT(in)                  :: lmax     ! zone number of first ghost zone on right (lmax=nmax+1)
INTEGER, INTENT(in)                  :: ij_ray   ! y-index of an x ray, x-index of an y ray, x index of a z array
INTEGER, INTENT(in)                  :: ik_ray   ! z-index of an x ray, z-index of an y ray, y index of a z array

REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: game     ! 1 + p/(rho*ei)
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: gamc     ! adiabatic gamma_1
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: plft     ! mean pressure over the left causal domain
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: prgh     ! mean pressure over the right causal domain
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: ulft     ! mean velocity over the left causal domain
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: urgh     ! mean velocity over the right causal domain
REAL(KIND=double), DIMENSION(max_12)              :: vlft     ! mean pecific volume (initially density) over the left causal domain
REAL(KIND=double), DIMENSION(max_12)              :: vrgh     ! mean pecific volume (initially density) over the right causal domain
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: gclft    ! mean adiabatic gamma_1 over the left causal domain
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: gcrgh    ! mean adiabatic gamma_1 over the right causal domain
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: gelft    ! mean 1 + p/(rho*ei over the left causal domain
REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: gergh    ! mean 1 + p/(rho*ei over the right causal domain

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(max_12), INTENT(out) :: pmid     ! time averaged pressure at the zone interface
REAL(KIND=double), DIMENSION(max_12), INTENT(out) :: umid     ! time averaged velocity at the zone interface

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                              :: l        ! zone index
INTEGER                              :: n        ! iteration index

REAL(KIND=double), DIMENSION(max_12) :: gamfl    ! 0.5*(1+gclft)/gclft
REAL(KIND=double), DIMENSION(max_12) :: gamfr    ! 0.5*(1+gcrgh)/gcrgh
REAL(KIND=double), DIMENSION(max_12) :: clft     ! mean sound speed over the left causal domain
REAL(KIND=double), DIMENSION(max_12) :: crgh     ! mean sound speed over the right causal domain
REAL(KIND=double), DIMENSION(max_12) :: wlft     ! mean wave over the left causal domain
REAL(KIND=double), DIMENSION(max_12) :: wrgh     ! mean wave over the right causal domain
REAL(KIND=double), DIMENSION(max_12) :: zlft     ! time averaged pressure at the zone interface
REAL(KIND=double), DIMENSION(max_12) :: zrgh     ! time averaged velocity at the zone interface
REAL(KIND=double), DIMENSION(max_12) :: umidl    ! "old" value of umid on left during iteration
REAL(KIND=double), DIMENSION(max_12) :: umidr    ! "old" value of umid on right during iteration
REAL(KIND=double), DIMENSION(max_12) :: pmold    ! "old" value of pmid during iteration

REAL(KIND=double), PARAMETER         :: tol = 1.0d-3
REAL(KIND=double), PARAMETER         :: smallp = 1.d-15

!-----------------------------------------------------------------------
!  Obtain first guess for pmid by assuming wlft, wrgh = clft, crgh
!
!  pmid and umid are given by
!
!     pmid - p_left
!     -------------  + (umid - u_left ) = 0
!        w_left
!
!     pmid - p_right
!     -------------- + (umid - u_right) = 0
!        w_right
!-----------------------------------------------------------------------

DO l = lmin, lmax
  clft(l)            = DSQRT( gclft(l) * plft(l) * vlft(l) )
  crgh(l)            = DSQRT( gcrgh(l) * prgh(l) * vrgh(l) )
  vlft(l)            = 1.d0/vlft(l)
  vrgh(l)            = 1.d0/vrgh(l)
  pmid(l)            = prgh(l) - plft(l) - crgh(l) * ( urgh(l) - ulft(l) )     
  pmid(l)            = plft(l) + pmid(l) * clft(l)/( clft(l) + crgh(l) )      
  pmid(l)            = DMAX1( smallp, pmid(l) ) 
END DO

!-----------------------------------------------------------------------
! Iterate 5 times using the Newton method to converge on correct Pmid
!     -use previous guess for pmid to get wavespeeds: wlft, wrgh
!     -find the slope in the u-P plane for each state: zlft, zrgh
!     -use the wavespeeds and pmid to guess umid on each side: umidl, umidr
!     -project tangents from (pmid,umidl) and (pmid,umidr) to get new pmid
!     -make sure pmid does not fall below floor value for pressure
!-----------------------------------------------------------------------

DO l = lmin, lmax
  DO n = 1, 8
    pmold(l)         = pmid(l)
    gamfl(l)         = 0.5d0 * ( gelft(l) + 1.0d0 )/gelft(l)
    gamfr(l)         = 0.5d0 * ( gergh(l) + 1.0d0 )/gergh(l)
    wlft (l)         = 1.0d0 + gamfl(l) * ( pmid(l) - plft(l) )/plft(l)   
    wrgh (l)         = 1.0d0 + gamfr(l) * ( pmid(l) - prgh(l) )/prgh(l)
    wlft (l)         = clft(l) * DSQRT( wlft(l) )      
    wrgh (l)         = crgh(l) * DSQRT( wrgh(l) )   
    zlft (l)         = 4.0d0 * vlft(l) * wlft(l) * wlft(l)      
    zrgh (l)         = 4.0d0 * vrgh(l) * wrgh(l) * wrgh(l)   
    zlft (l)         = -zlft(l) * wlft(l)/( zlft(l) - ( 1.d0 + gelft(l) ) * ( pmid(l) - plft(l) ) )
    zrgh (l)         =  zrgh(l) * wrgh(l)/( zrgh(l) - ( 1.d0 + gergh(l) ) * ( pmid(l) - prgh(l) ) )
    umidl(l)         = ulft(l) - ( pmid(l) - plft(l) )/wlft(l)      
    umidr(l)         = urgh(l) + ( pmid(l) - prgh(l) )/wrgh(l)
    pmid (l)         = pmid(l) + ( umidr(l) - umidl(l) ) * ( zlft (l) * zrgh (l) )/( zrgh(l) - zlft(l) )   
    pmid (l)         = DMAX1( smallp, pmid(l) )
    IF ( DABS( pmid(l) - pmold(l) )/pmid(l) < tol ) EXIT
  END DO

!-----------------------------------------------------------------------
!  Calculate umid by averaging umidl, umidr based on new pmid
!-----------------------------------------------------------------------

  umidl(l)           = ulft(l) - ( pmid(l) - plft(l) )/wlft(l)      
  umidr(l)           = urgh(l) + ( pmid(l) - prgh(l) )/wrgh(l)   
  umid (l)           = 0.5d0 * ( umidl(l) + umidr(l) )
END DO

RETURN
END SUBROUTINE riemann
