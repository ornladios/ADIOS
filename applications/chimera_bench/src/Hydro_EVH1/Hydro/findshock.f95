SUBROUTINE findshock( jr_min, jr_max, ij_ray, ik_ray, pqmin, jjshock, &
& jjshockmx, jjshockmn, nx, x_c, m_c, r_shock, m_shock )
!-----------------------------------------------------------------------
!
!    File:         findshock
!    Module:       findshock
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/09/03
!
!    Purpose:
!      To find the outer and inner zones of a shock.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_max      : outer shifted radial zone to start the search
!  jr_min      : inner shifted radial zone to end search 
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!  pqmin       : minimum ratio of pseudoviscosity to pressure used
!                 to infer the presence of a shock
!  x_c         : x (radial) zone-centered coordinate (cm)
!  m_c         : enclosed mass to radial zone-center (g)
!  nx          : x-array extent
!
!    Output arguments:
!  jjshock     : radial zone at shock maximum
!  jjshockmx   : outer radial zone of shock
!  jjshockmn   : inner radial zone of shock
!  r_shock     : shock radius (cm)
!  m_shock     : rest mass enclosed by shock (solar masses)
!
!    Include files:
!  kind_module, numerical_module:
!  eos_snc_x_module, shock_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon
USE physcnst_module, ONLY : msolar

USE eos_snc_x_module, ONLY : aesv
USE shock_module, ONLY : pq_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min          ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max          ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray
INTEGER, INTENT(in)              :: nx              ! x-array extent

REAL(KIND=double), INTENT(in)    :: pqmin           ! minimum ratio of pseudoviscosity to pressure used
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: x_c ! x (radial) unshifted zone-centered coordinate (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: m_c ! enclosed mass to radial zone-center (g)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out)             :: jjshock         ! radial zone at shock maximum
INTEGER, INTENT(out)             :: jjshockmn       ! minimum radial zone index
INTEGER, INTENT(out)             :: jjshockmx       ! maximum radial zone index

REAL(KIND=double), INTENT(out)   :: r_shock         ! shock radius (g)
REAL(KIND=double), INTENT(out)   :: m_shock         ! rest mass enclosed by shock (solar masses)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j               ! radial zone index
INTEGER, DIMENSION(1)            :: i_shk           ! radial zone index of shock maximum

REAL(KIND=double)                :: pq_max          ! pq_x(j,i_ray)/( aesv(j,1,i_ray) at jjshock (i.e., its maximum value)
REAL(KIND=double)                :: pq_m            ! pq_x(j,i_ray)/( aesv(j,1,i_ray) at jjshock - 1
REAL(KIND=double)                :: pq_p            ! pq_x(j,i_ray)/( aesv(j,1,i_ray) at jjshock + 1
REAL(KIND=double), DIMENSION(nx) :: r_c             ! x (radial) shifted zone-centered coordinate (cm)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                 \\\\\ RADIAL INDEX OF SHOCK ///// 
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find jjshockmx, the maximum radial index for which
!     pq_x(j,ij_ray,ik_ray)/( aesv(j,1,ij_ray,ik_ray) > pqmin
!-----------------------------------------------------------------------

DO j = jr_max,jr_min,-1
  jjshockmx         = j 
  IF ( pq_x(j,ij_ray,ik_ray)/( aesv(j,1,ij_ray,ik_ray) + epsilon ) >= pqmin ) EXIT
END DO

IF ( jjshockmx == jr_min ) THEN
  jjshockmn         = jr_min
  jjshock           = jr_min
  RETURN
END IF ! jjshockmx = jr_min

!-----------------------------------------------------------------------
!  Find jjshockmn, the minimum radial index for which
!     pq_x(j,ij_ray,ik_ray)/( aesv(j,1,ij_ray,ik_ray) > pqmin
!-----------------------------------------------------------------------

DO j = jjshockmx-1,jr_min,-1
  jjshockmn         = j + 1
  IF ( pq_x(j,ij_ray,ik_ray)/( aesv(j,1,ij_ray,ik_ray) + epsilon ) < pqmin ) EXIT
END DO


!-----------------------------------------------------------------------
!  jjshock is the index giving the maximum pq_x(:,i_ray)/aesv(:,1,i_ray)
!-----------------------------------------------------------------------

i_shk               = MAXLOC( pq_x(:,ij_ray,ik_ray)                     &
&                   / ( aesv(:,1,ij_ray,ik_ray) + epsilon ),            &
&                   pq_x(:,ij_ray,ik_ray)/( aesv(:,1,ij_ray,ik_ray) + epsilon ) >= pqmin )
jjshock             = MAX( i_shk(1), jr_min )

!-----------------------------------------------------------------------
!
!                      \\\\\ SHOCK RADIUS ///// 
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Load r_c
!-----------------------------------------------------------------------

r_c(jr_min:jr_max)  = x_c(jr_min-1:jr_max-1)

!-----------------------------------------------------------------------
!  Set r_shock to zero if jjshock = jr_min
!-----------------------------------------------------------------------

IF ( jjshock == jr_min ) THEN
  r_shock           = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Set r_shock to x_c(jr_max) if jjshock = jr_max
!-----------------------------------------------------------------------

IF ( jjshock == jr_max ) THEN
  r_shock           = x_c(jr_max)
  RETURN
END IF

!-----------------------------------------------------------------------
!  Quadratically fit pq_x(:,i_ray)/aesv(:,1,i_ray) and set r_shock to
!   the maximum
!-----------------------------------------------------------------------

pq_max              = pq_x(jjshock  ,ij_ray,ik_ray)/( aesv(jjshock  ,1,ij_ray,ik_ray) + epsilon )
pq_m                = pq_x(jjshock-1,ij_ray,ik_ray)/( aesv(jjshock-1,1,ij_ray,ik_ray) + epsilon )
pq_p                = pq_x(jjshock+1,ij_ray,ik_ray)/( aesv(jjshock+1,1,ij_ray,ik_ray) + epsilon )

r_shock             = r_max( x_c(jjshock-1), x_c(jjshock), x_c(jjshock+1), &
& pq_m, pq_max,  pq_p )

m_shock             = r_max( m_c(jjshock-1), m_c(jjshock), m_c(jjshock+1), &
& pq_m, pq_max,  pq_p )
m_shock             = m_shock/msolar

RETURN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS
REAL (KIND=double) FUNCTION r_max( x1, x2, x3, y1, y2, y3 )

REAL (KIND=double) :: x1
REAL (KIND=double) :: x2
REAL (KIND=double) :: x3
REAL (KIND=double) :: y1
REAL (KIND=double) :: y2
REAL (KIND=double) :: y3
REAL (KIND=double) :: a
REAL (KIND=double) :: b
REAL (KIND=double) :: denom

denom               = 1.d0/( ( x2 - x1 ) * ( x3 - x2 ) * ( x1 - x3 ) )

a                   = ( ( y2 - y1 ) * ( x3 - x1 ) - ( y3 - y1 ) * ( x2 - x1 ) ) * denom
b                   = ( ( y3 - y1 ) * ( x2 - x1 ) * ( x2 + x1 )           &
&                   -   ( y2 - y1 ) * ( x3 - x1 ) * ( x3 + x1 ) ) * denom
r_max               = - b/( 2.d0 * a + epsilon )
END FUNCTION r_max

END SUBROUTINE findshock
