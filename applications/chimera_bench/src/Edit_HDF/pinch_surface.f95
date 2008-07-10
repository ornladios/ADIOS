SUBROUTINE pinch_surface( jr_min, jr_max, i_ray, i_ray_dim, nnu, r_pinch, &
& d_pinch, stat_pinch_r, trns_pinch_r, stat_pinch_d, trns_pinch_d )
!-----------------------------------------------------------------------
!
!    File:         pinch_surface
!    Module:       pinch_surface
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To calculate the neutrino spectral pinching at the radius r_pinch and
!       the density d_pinch.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min          : inner radial zone number
!  jr_max          : outer radial zone number
!  i_ray           : index denoting a specific radial ray
!  i_ray_dim       : number radial rays assigned to a processor
!  nnu             : neutrino flavor array extent
!  r_pinch         : radius at which to compute the special pinching
!  d_pinch         : density at which to compute the special pinching
!
!    Output arguments:
!  stat_pinch_r    : static spectral pinch factor at radius r_pinch
!  trns_pinch_r    : transport spectral pinch factor at radius r_pinch
!  stat_pinch_d    : static spectral pinch factor at density d_pinch
!  trns_pinch_d    : transport spectral pinch factor at density d_pinch
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon

USE edit_module, ONLY : nprint
USE mdl_cnfg_module, ONLY : r, rho
USE nu_dist_module, ONLY : psi0, psi1, unue, dunue, rjmh
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: jr_min        ! minimum radial index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial index
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)               :: i_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)               :: i_ray_dim     ! number radial rays assigned to a processor

REAL(KIND=double), INTENT(in)     :: r_pinch       ! radius at which to compute the special pinching
REAL(KIND=double), INTENT(in)     :: d_pinch       ! density at which to compute the special pinching

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: stat_pinch_r ! static spectral pinch factor at radius r_pinch
REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: trns_pinch_r ! transport spectral pinch factor at radius r_pinch
REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: stat_pinch_d ! static spectral pinch factor at density d_pinch
REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: trns_pinch_d ! transport spectral pinch factor at density d_pinch

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                           :: l_r_pinch
LOGICAL                           :: l_d_pinch

INTEGER                           :: j               ! radial zone index
INTEGER                           :: jd              ! particular radial zone index
INTEGER                           :: k               ! neutrino energy index
INTEGER                           :: n               ! neutrino flavor index

REAL(KIND=double)                 :: w2dw            ! unue(j)**2 * dunue(j)
REAL(KIND=double)                 :: w3dw            ! unue(j)**3 * dunue(j)
REAL(KIND=double)                 :: w4dw            ! unue(j)**3 * dunue(j)

REAL(KIND=double)                 :: psi0j           ! psi0(j,k,n) interpolated to the zone edge
REAL(KIND=double)                 :: psi0j2          ! psi0(j,k,n) * unue(j)**2 * dunue(j)
REAL(KIND=double)                 :: psi0j3          ! psi0(j,k,n) * unue(j)**3 * dunue(j)
REAL(KIND=double)                 :: psi0j4          ! psi0(j,k,n) * unue(j)**4 * dunue(j)
REAL(KIND=double)                 :: psi1j2          ! psi1(j,k,n) * unue(j)**2 * dunue(j)
REAL(KIND=double)                 :: psi1j3          ! psi1(j,k,n) * unue(j)**3 * dunue(j)
REAL(KIND=double)                 :: psi1j4          ! psi1(j,k,n) * unue(j)**4 * dunue(j)

REAL(KIND=double)                 :: e_density       ! mean neutrin energy density
REAL(KIND=double)                 :: e2_density      ! mean squaare neutrin energy density

REAL(KIND=double)                 :: stat_pinch_jd   ! static spectral pinch factor at jd
REAL(KIND=double)                 :: trns_pinch_jd   ! transport spectral pinch factor at jd
REAL(KIND=double)                 :: stat_pinch_jdm1 ! static spectral pinch factor at jd-1
REAL(KIND=double)                 :: trns_pinch_jdm1 ! transport spectral pinch factor at jd-1

 1001 format (' jd cannot be found in subroutine pinch_surface for r_pinch=',es11.3)
 4001 format (' jd cannot be found in subroutine pinch_surface for d_pinch=',es11.3)

!-----------------------------------------------------------------------
!        Initialize.
!-----------------------------------------------------------------------

DO n = 1,nnu
  stat_pinch_r(n,i_ray) = zero
  trns_pinch_r(n,i_ray) = zero
  stat_pinch_d(n,i_ray) = zero
  trns_pinch_d(n,i_ray) = zero
END DO

!-----------------------------------------------------------------------
!        Return if no neutrinos.
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!          \\\\\ SPECTRAL PINCHING AT RADIUS R_PINCH /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Find jd such that r(jd) > r_pinch and r(jd-1) < r_pinch.
!-----------------------------------------------------------------------

l_r_pinch               = .false.
DO j = jr_min,jr_max
  IF ( r(j) >= r_pinch ) THEN
    jd                  = j
    l_r_pinch           = .true.
    EXIT
  END IF ! r(j) > r_pinch
END DO

IF ( .not. l_r_pinch ) WRITE (nprint,1001) r_pinch

!-----------------------------------------------------------------------
!        Compute the pinching factors at radius r_pinch.
!-----------------------------------------------------------------------

IF ( l_r_pinch ) THEN

  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN

      psi0j2          = zero
      psi0j3          = zero
      psi0j4          = zero
      psi1j2          = zero
      psi1j3          = zero
      psi1j4          = zero
      DO k = 1,nnugp(n)
        w2dw          = unue(jd,k)**2 * dunue(jd,k)
        w3dw          = unue(jd,k)**3 * dunue(jd,k)
        w4dw          = unue(jd,k)**4 * dunue(jd,k)
        psi0j         = ( ( rjmh(jd+1)**2 - r(jd)**2 ) * psi0(jd,k,n) &
&                     + ( r(jd)**2 - rjmh(jd)**2 ) * psi0(jd+1,k,n) ) &
&                     / ( rjmh(jd+1)**2 - rjmh(jd)**2 )
        psi0j2        = psi0j2 + w2dw  * psi0j
        psi0j3        = psi0j3 + w3dw  * psi0j
        psi0j4        = psi0j4 + w4dw  * psi0j
        psi1j2        = psi1j2 + w2dw  * psi1(jd,k,n)
        psi1j3        = psi1j3 + w3dw  * psi1(jd,k,n)
        psi1j4        = psi1j4 + w4dw  * psi1(jd,k,n)
      END DO ! k
      e_density       = psi0j3/( psi0j2 + epsilon )
      e2_density      = psi0j4/( psi0j2 + epsilon )
      stat_pinch_jd   = 0.75d0 * e2_density/( e_density * e_density + epsilon )
      e_density       = psi1j3/( psi1j2 + epsilon )
      e2_density      = psi1j4/( psi1j2 + epsilon )
      trns_pinch_jd   = 0.75d0 * e2_density/( e_density * e_density + epsilon )

      psi0j2          = zero
      psi0j3          = zero
      psi0j4          = zero
      psi1j2          = zero
      psi1j3          = zero
      psi1j4          = zero
      DO k = 1,nnugp(n)
        w2dw          = unue(jd-1,k)**2 * dunue(jd-1,k)
        w3dw          = unue(jd-1,k)**3 * dunue(jd-1,k)
        w4dw          = unue(jd-1,k)**4 * dunue(jd-1,k)
        psi0j         = ( ( rjmh(jd)**2 - r(jd-1)**2 ) * psi0(jd-1,k,n) + ( r(jd-1)**2 - rjmh(jd-1)**2 ) * psi0(jd,k,n) ) &
&                     / ( rjmh(jd)**2 - rjmh(jd-1)**2 )
        psi0j2        = psi0j2 + w2dw  * psi0j
        psi0j3        = psi0j3 + w3dw  * psi0j
        psi0j4        = psi0j4 + w4dw  * psi0j
        psi1j2        = psi1j2 + w2dw  * psi1(jd-1,k,n)
        psi1j3        = psi1j3 + w3dw  * psi1(jd-1,k,n)
        psi1j4        = psi1j4 + w4dw  * psi1(jd-1,k,n)
      END DO ! k
      e_density       = psi0j3/( psi0j2 + epsilon )
      e2_density      = psi0j4/( psi0j2 + epsilon )
      stat_pinch_jdm1 = 0.75d0 * e2_density/( e_density * e_density + epsilon )
      e_density       = psi1j3/( psi1j2 + epsilon )
      e2_density      = psi1j4/( psi1j2 + epsilon )
      trns_pinch_jdm1 = 0.75d0 * e2_density/( e_density * e_density + epsilon )

    END IF ! nnugp(n) /= 0

    stat_pinch_r(n,i_ray)= rinterp( stat_pinch_jd, stat_pinch_jdm1, r(jd), r_pinch, r(jd-1) )
    trns_pinch_r(n,i_ray)= rinterp( trns_pinch_jd, trns_pinch_jdm1, r(jd), r_pinch, r(jd-1) )

  END DO

!-----------------------------------------------------------------------
!        End computation of pinch factors at radius r_pinch.
!-----------------------------------------------------------------------

END IF ! l_r_pinch

!-----------------------------------------------------------------------
!
!          \\\\\ SPECTRAL PINCHING AT DENSITY D_PINCH /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Find jd such that rho(jd) < d_pinch and rho(jd-1) > d_pinch.
!-----------------------------------------------------------------------

l_d_pinch             = .false.
DO j = jr_min,jr_max
  IF ( rho(j) <= d_pinch ) THEN
    jd                = j
    l_d_pinch         = .true.
    EXIT
  END IF ! rho(j) <= d_pinch
END DO

IF ( .not. l_d_pinch ) WRITE (nprint,4001) d_pinch

!-----------------------------------------------------------------------
!        Compute the luminosities at density d_lum.
!-----------------------------------------------------------------------

IF ( l_d_pinch ) THEN

  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN

      psi0j2          = zero
      psi0j3          = zero
      psi0j4          = zero
      psi1j2          = zero
      psi1j3          = zero
      psi1j4          = zero
      DO k = 1,nnugp(n)
        w2dw          = unue(jd,k)**2 * dunue(jd,k)
        w3dw          = unue(jd,k)**3 * dunue(jd,k)
        w4dw          = unue(jd,k)**4 * dunue(jd,k)
        psi0j         = ( ( rjmh(jd+1)**2 - r(jd)**2 ) * psi0(jd,k,n) &
&                     + ( r(jd)**2 - rjmh(jd)**2 ) * psi0(jd+1,k,n) ) &
&                     / ( rjmh(jd+1)**2 - rjmh(jd)**2 )
        psi0j2        = psi0j2 + w2dw  * psi0j
        psi0j3        = psi0j3 + w3dw  * psi0j
        psi0j4        = psi0j4 + w4dw  * psi0j
        psi1j2        = psi1j2 + w2dw  * psi1(jd,k,n)
        psi1j3        = psi1j3 + w3dw  * psi1(jd,k,n)
        psi1j4        = psi1j4 + w4dw  * psi1(jd,k,n)
      END DO ! k
      e_density       = psi0j3/( psi0j2 + epsilon )
      e2_density      = psi0j4/( psi0j2 + epsilon )
      stat_pinch_jd   = 0.75d0 * e2_density/( e_density * e_density + epsilon )
      e_density       = psi1j3/( psi1j2 + epsilon )
      e2_density      = psi1j4/( psi1j2 + epsilon )
      trns_pinch_jd   = 0.75d0 * e2_density/( e_density * e_density + epsilon )

      psi0j2          = zero
      psi0j3          = zero
      psi0j4          = zero
      psi1j2          = zero
      psi1j3          = zero
      psi1j4          = zero
      DO k = 1,nnugp(n)
        w2dw          = unue(jd-1,k)**2 * dunue(jd-1,k)
        w3dw          = unue(jd-1,k)**3 * dunue(jd-1,k)
        w4dw          = unue(jd-1,k)**4 * dunue(jd-1,k)
        psi0j         = ( ( rjmh(jd)**2 - r(jd-1)**2 ) * psi0(jd-1,k,n) + ( r(jd-1)**2 - rjmh(jd-1)**2 ) * psi0(jd,k,n) ) &
&                     / ( rjmh(jd)**2 - rjmh(jd-1)**2 )
        psi0j2        = psi0j2 + w2dw  * psi0j
        psi0j3        = psi0j3 + w3dw  * psi0j
        psi0j4        = psi0j4 + w4dw  * psi0j
        psi1j2        = psi1j2 + w2dw  * psi1(jd-1,k,n)
        psi1j3        = psi1j3 + w3dw  * psi1(jd-1,k,n)
        psi1j4        = psi1j4 + w4dw  * psi1(jd-1,k,n)
      END DO ! k
      e_density       = psi0j3/( psi0j2 + epsilon )
      e2_density      = psi0j4/( psi0j2 + epsilon )
      stat_pinch_jdm1 = 0.75d0 * e2_density/( e_density * e_density + epsilon )
      e_density       = psi1j3/( psi1j2 + epsilon )
      e2_density      = psi1j4/( psi1j2 + epsilon )
      trns_pinch_jdm1 = 0.75d0 * e2_density/( e_density * e_density + epsilon )

    END IF ! nnugp(n) /= 0

    stat_pinch_d(n,i_ray)= rinterp( stat_pinch_jd, stat_pinch_jdm1, rho(jd), d_pinch, rho(jd-1) )
    trns_pinch_d(n,i_ray)= rinterp( trns_pinch_jd, trns_pinch_jdm1, rho(jd), d_pinch, rho(jd-1) )

  END DO

!-----------------------------------------------------------------------
!        End computation of pinch factors at density d_pinch.
!-----------------------------------------------------------------------

END IF ! l_d_pinch

!-----------------------------------------------------------------------
!        Return.
!-----------------------------------------------------------------------

RETURN

CONTAINS
  REAL (KIND=double) FUNCTION rinterp(a,b,x,y,z)

  REAL (KIND=double) :: a
  REAL (KIND=double) :: b
  REAL (KIND=double) :: x
  REAL (KIND=double) :: y
  REAL (KIND=double) :: z

  rinterp      = b + ( a - b ) * ( y - z )/( x - z )

END FUNCTION rinterp

END SUBROUTINE pinch_surface
