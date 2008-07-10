SUBROUTINE e_rms_surface( jr_min, jr_max, i_ray, i_ray_dim, nnu, r_e_rms, &
& d_e_rms, e_rms_r_stat, e_rms_r_trns, e_rms_d_stat, e_rms_d_trns )
!-----------------------------------------------------------------------
!
!    File:         e_rms_surface
!    Module:       e_rms_surface
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To calculate the rms neutrino energies at the given radius r_e_rms, and
!       at the given density d_e_rms.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!   jr_min          : inner radial zone number
!   jr_max          : outer radial zone number
!   i_ray           : index denoting a specific radial ray
!   i_ray_dim       : number radial rays assigned to a processor
!   nnu             : neutrino flavor array extent
!   r_e_rms         : radius at which to compute rms neutrino energies
!   d_e_rms         : density at which to compute rms neutrino energies
!
!    Output arguments:
!   e_rms_stat      : rms static neutrino energy at radius r_e_rms
!   e_rms_r_trns    : rms transport neutrino energy at radius r_e_rms
!   e_rms_d_stat    : rms static neutrino energy at density d_e_rms
!   e_rms_d_trns    : rms transport neutrino energy at density d_e_rms
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon
USE physcnst_module, ONLY : ergfoe

USE edit_module, ONLY : nprint
USE mdl_cnfg_module, ONLY : r, rho
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: jr_min            ! minimum radial index
INTEGER, INTENT(in)               :: jr_max            ! maximum radial index
INTEGER, INTENT(in)               :: nnu               ! neutrino flavor array extent
INTEGER, INTENT(in)               :: i_ray             ! index denoting a specific radial ray
INTEGER, INTENT(in)               :: i_ray_dim         ! number radial rays assigned to a processor

REAL(KIND=double), INTENT(in)     :: r_e_rms           ! radius at which to compute rms neutrino energies
REAL(KIND=double), INTENT(in)     :: d_e_rms           ! density at which to compute rms neutrino energies

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: e_rms_r_stat ! rms static neutrino energy at radius r_e_rms
REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: e_rms_r_trns ! rms transport neutrino energy at radius r_e_rms
REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: e_rms_d_stat ! rms static neutrino energy at density d_e_rms
REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: e_rms_d_trns ! rms transport neutrino energy at density d_e_rms

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                           :: l_r_e_rms
LOGICAL                           :: l_d_e_rms

INTEGER                           :: j                 ! radial zone index
INTEGER                           :: jd                ! particular radial zone index
INTEGER                           :: k                 ! neutrino energy index
INTEGER                           :: n                 ! neutrino flavor index

REAL(KIND=double)                 :: w3dw              ! unu(j)**3 * dunu(j)
REAL(KIND=double)                 :: w5dw              ! unu(j)**5 * dunu(j)

REAL(KIND=double)                 :: psi0j3            ! psi0(j,k,n) * unu(j)**3 * dunu(j)
REAL(KIND=double)                 :: psi0j5            ! psi0(j,k,n) * unu(j)**5 * dunu(j)
REAL(KIND=double)                 :: psi1j3            ! psi1(j,k,n) * unu(j)**3 * dunu(j)
REAL(KIND=double)                 :: psi1j5            ! psi1(j,k,n) * unu(j)**5 * dunu(j)

REAL(KIND=double)                 :: e_rms_r_stat_jd   ! rms static neutrino energy at zone jd
REAL(KIND=double)                 :: e_rms_r_trns_jd   ! rms transport neutrino energy at zone j
REAL(KIND=double)                 :: e_rms_r_stat_jdm1 ! rms static neutrino energy at zone jd-1
REAL(KIND=double)                 :: e_rms_r_trns_jdm1 ! rms transport neutrino energy at zone jd-1
REAL(KIND=double)                 :: e_rms_d_stat_jd   ! rms static neutrino energy at zone jd
REAL(KIND=double)                 :: e_rms_d_trns_jd   ! rms transport neutrino energy at zone j
REAL(KIND=double)                 :: e_rms_d_stat_jdm1 ! rms static neutrino energy at zone jd-1
REAL(KIND=double)                 :: e_rms_d_trns_jdm1 ! rms transport neutrino energy at zone jd-1

 1001 format (' jd cannot be found in subroutine e_rms_surface for r_e_rms=',es11.3)
 4001 format (' jd cannot be found in subroutine e_rms_surface for d_e_rms',es11.3)

!-----------------------------------------------------------------------
!        Initialize.
!-----------------------------------------------------------------------

DO n = 1,nnu
  e_rms_r_stat(n,i_ray) = zero
  e_rms_r_trns(n,i_ray) = zero
  e_rms_d_stat(n,i_ray) = zero
  e_rms_d_trns(n,i_ray) = zero
END DO

!-----------------------------------------------------------------------
!        Return if no neutrinos.
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!          \\\\\ RMS NEUTRINO ENERGIES AT RADIUS RENU /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Find jd such that r(jd) > r_e_rms and r(jd-1) < r_e_rms.
!-----------------------------------------------------------------------

l_r_e_rms                = .false.
DO j = jr_min,jr_max
  IF ( r(j) >= r_e_rms ) THEN
    jd                   = j
    l_r_e_rms            = .true.
    EXIT
  END IF ! r(j) > r_e_rms
END DO

IF ( .not. l_r_e_rms ) WRITE (nprint,1001) r_e_rms

!-----------------------------------------------------------------------
!        Compute the luminosities at radius r_e_rms.
!-----------------------------------------------------------------------

IF ( l_r_e_rms ) THEN

  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN

      psi0j3             = zero
      psi0j5             = zero
      psi1j3             = zero
      psi1j5             = zero
      DO k = 1,nnugp(n)
        w3dw             = unu(jd,k)**3 * dunu(jd,k)
        w5dw             = unu(jd,k)**5 * dunu(jd,k)
        psi0j3           = psi0j3 + w3dw  * psi0(jd,k,n)
        psi0j5           = psi0j5 + w5dw  * psi0(jd,k,n)
        psi1j3           = psi1j3 + w3dw  * psi1(jd,k,n)
        psi1j5           = psi1j5 + w5dw  * psi1(jd,k,n)
      END DO ! k
      e_rms_r_stat_jd    = DSQRT( DABS( psi0j5 /( psi0j3 + epsilon ) ) + epsilon )
      e_rms_r_trns_jd    = DSQRT( DABS( psi1j5 /( psi1j3 + epsilon ) ) + epsilon )

      psi0j3             = zero
      psi0j5             = zero
      psi1j3             = zero
      psi1j5             = zero
      DO k = 1,nnugp(n)
        w3dw             = unu(jd-1,k)**3 * dunu(jd-1,k)
        w5dw             = unu(jd-1,k)**5 * dunu(jd-1,k)
        psi0j3           = psi0j3 + w3dw  * psi0(jd-1,k,n)
        psi0j5           = psi0j5 + w5dw  * psi0(jd-1,k,n)
        psi1j3           = psi1j3 + w3dw  * psi1(jd-1,k,n)
        psi1j5           = psi1j5 + w5dw  * psi1(jd-1,k,n)
      END DO ! k
      e_rms_r_stat_jdm1  = DSQRT( DABS( psi0j5 /( psi0j3 + epsilon ) ) + epsilon )
      e_rms_r_trns_jdm1  = DSQRT( DABS( psi1j5 /( psi1j3 + epsilon ) ) + epsilon )

    END IF ! nnugp(n) /= 0

    e_rms_r_stat(n,i_ray)= rinterp( e_rms_r_stat_jd, e_rms_r_stat_jdm1, r(jd), r_e_rms, r(jd-1) )
    e_rms_r_trns(n,i_ray)= rinterp( e_rms_r_trns_jd, e_rms_r_trns_jdm1, r(jd), r_e_rms, r(jd-1) )

  END DO

!-----------------------------------------------------------------------
!        End computation of luminosities at r_e_rms.
!-----------------------------------------------------------------------

END IF ! l_r_e_rms

!-----------------------------------------------------------------------
!
!      \\\\\ MEAN NEUTRINO LUMINOSITIES AT DENSITY RHOLUM /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Find jd such that rho(jd) < d_e_rms and rho(jd-1) > d_e_rms.
!-----------------------------------------------------------------------

l_d_e_rms                = .false.
DO j = jr_min,jr_max
  IF ( rho(j) <= d_e_rms ) THEN
    jd                   = j
    l_d_e_rms            = .true.
    EXIT
  END IF ! rho(j) <= d_e_rms
END DO

IF ( .not. l_d_e_rms ) WRITE (nprint,4001) d_e_rms

!-----------------------------------------------------------------------
!        Compute the luminosities at density d_e_rms.
!-----------------------------------------------------------------------

IF ( l_d_e_rms ) THEN

  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN

      psi0j3             = zero
      psi0j5             = zero
      psi1j3             = zero
      psi1j5             = zero
      DO k = 1,nnugp(n)
        w3dw             = unu(jd,k)**3 * dunu(jd,k)
        w5dw             = unu(jd,k)**5 * dunu(jd,k)
        psi0j3           = psi0j3 + w3dw  * psi0(jd,k,n)
        psi0j5           = psi0j5 + w5dw  * psi0(jd,k,n)
        psi1j3           = psi1j3 + w3dw  * psi1(jd,k,n)
        psi1j5           = psi1j5 + w5dw  * psi1(jd,k,n)
      END DO ! k
      e_rms_d_stat_jd = DSQRT( DABS( psi0j5 /( psi0j3 + epsilon ) ) + epsilon )
      e_rms_d_trns_jd = DSQRT( DABS( psi1j5 /( psi1j3 + epsilon ) ) + epsilon )

      psi0j3             = zero
      psi0j5             = zero
      psi1j3             = zero
      psi1j5             = zero
      DO k = 1,nnugp(n)
        w3dw             = unu(jd-1,k)**3 * dunu(jd-1,k)
        w5dw             = unu(jd-1,k)**5 * dunu(jd-1,k)
        psi0j3           = psi0j3 + w3dw  * psi0(jd-1,k,n)
        psi0j5           = psi0j5 + w5dw  * psi0(jd-1,k,n)
        psi1j3           = psi1j3 + w3dw  * psi1(jd-1,k,n)
        psi1j5           = psi1j5 + w5dw  * psi1(jd-1,k,n)
      END DO ! k
      e_rms_d_stat_jdm1  = DSQRT( DABS( psi0j5 /( psi0j3 + epsilon ) ) + epsilon )
      e_rms_d_trns_jdm1  = DSQRT( DABS( psi1j5 /( psi1j3 + epsilon ) ) + epsilon )

    END IF ! nnugp(n) /= 0

    e_rms_d_stat(n,i_ray)= rinterp( e_rms_d_stat_jd, e_rms_d_stat_jdm1, rho(jd), d_e_rms, rho(jd-1) )
    e_rms_d_trns(n,i_ray)= rinterp( e_rms_d_trns_jd, e_rms_d_trns_jdm1, rho(jd), d_e_rms, rho(jd-1) )

  END DO

!-----------------------------------------------------------------------
!        End computation of luminosities at d_e_rms.
!-----------------------------------------------------------------------

END IF ! l_d_e_rms

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

END SUBROUTINE e_rms_surface
