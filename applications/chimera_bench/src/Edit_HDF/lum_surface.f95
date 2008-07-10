SUBROUTINE lum_surface( jr_min, jr_max, i_ray, i_ray_dim, nnu, r_lum, d_lum, &
& lum_r, lum_rho )
!-----------------------------------------------------------------------
!
!    File:         lum_surface
!    Module:       lum_surface
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To calculate the neutrino luminosities at the given radius r_lum, and
!       at the given density d_lum.
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
!  r_lum           : radius at which to compute luminosities
!  d_lum           : density at which to compute luminosities
!
!    Output arguments:
!  lum_r           : luminosity at radius r_lum
!  lum_rho         : luminosity at density d_lum
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpi
USE physcnst_module, ONLY : ergfoe, cvel

USE edit_module, ONLY : nprint
USE mdl_cnfg_module, ONLY : r, rho
USE nu_dist_module, ONLY : psi1, ecoefae, stwt
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

REAL(KIND=double), INTENT(in)     :: r_lum         ! radius at which to compute luminosities
REAL(KIND=double), INTENT(in)     :: d_lum         ! density at which to compute luminosities

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: lum_r   ! luminosity at radius r_lum
REAL(KIND=double), INTENT(out), DIMENSION(nnu,i_ray_dim) :: lum_rho ! luminosity at density d_lum

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                           :: l_r_lum
LOGICAL                           :: l_d_lum

INTEGER                           :: j             ! radial zone index
INTEGER                           :: jd            ! particular radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double)                 :: flux_jd       ! neutrino flux at zone jd
REAL(KIND=double)                 :: flux_jdm1     ! neutrino flux at zone jd-1
REAL(KIND=double)                 :: area_jd       ! area at jd
REAL(KIND=double)                 :: area_jdm1     ! area at zone jd-1
REAL(KIND=double)                 :: r_lumjd        ! luminosity at zone jd (determined by r_lum)
REAL(KIND=double)                 :: r_lumjdm1      ! luminosity at zone jd-1 (determined by r_lum)
REAL(KIND=double)                 :: d_lumjd      ! luminosity at zone jd (determined by d_lum)
REAL(KIND=double)                 :: d_lumjdm1    ! luminosity at zone jd-1 (determined by d_lum)

 1001 format (' jd cannot be found in subroutine lum_surface for r_lum',es11.3)
 4001 format (' jd cannot be found in subroutine lum_surface for d_lum=',es11.3)

!-----------------------------------------------------------------------
!        Initialize.
!-----------------------------------------------------------------------

DO n = 1,nnu
  lum_r(n,i_ray)       = zero
  lum_rho(n,i_ray)     = zero
END DO

!-----------------------------------------------------------------------
!        Return if no neutrinos.
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!        \\\\\ MEAN NEUTRINO LUMINOSITIES AT RADIUS RLUM /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Find jd such that r(jd) > r_lum and r(jd-1) < r_lum.
!-----------------------------------------------------------------------

l_r_lum                = .false.
DO j = jr_min,jr_max
  IF ( r(j) >= r_lum ) THEN
    jd                 = j
    l_r_lum            = .true.
    EXIT
  END IF ! r(j) > r_lum
END DO

IF ( .not. l_r_lum ) WRITE (nprint,1001) r_lum

!-----------------------------------------------------------------------
!        Compute the luminosities at radius r_lum.
!-----------------------------------------------------------------------

IF ( l_r_lum ) THEN

  area_jd              = frpi * r(jd  ) * r(jd  )
  area_jdm1            = frpi * r(jd-1) * r(jd-1)
  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN
      flux_jd          = zero
      flux_jdm1        = zero
      DO k = 1,nnugp(n)
        flux_jd        = flux_jd   + cvel * psi1(jd,  k,n) * ecoefae(jd,  k) * stwt(n)
        flux_jdm1      = flux_jdm1 + cvel * psi1(jd-1,k,n) * ecoefae(jd-1,k) * stwt(n)
      END DO
    END IF ! nnugp(n) /= 0
    r_lumjd            = flux_jd   * area_jd   * ergfoe
    r_lumjdm1          = flux_jdm1 * area_jdm1 * ergfoe
    lum_r(n,i_ray)     = rinterp( r_lumjd, r_lumjdm1, r(jd), r_lum, r(jd-1) )
  END DO

!-----------------------------------------------------------------------
!        End computation of luminosities at r_lum.
!-----------------------------------------------------------------------

END IF ! l_r_lum

!-----------------------------------------------------------------------
!
!      \\\\\ MEAN NEUTRINO LUMINOSITIES AT DENSITY RHOLUM /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Find jd such that rho(jd) < d_lum and rho(jd-1) > d_lum.
!-----------------------------------------------------------------------

l_d_lum                = .false.
DO j = jr_min,jr_max
  IF ( rho(j) <= d_lum ) THEN
    jd                 = j
    l_d_lum            = .true.
    EXIT
  END IF ! rho(j) <= d_lum
END DO

IF ( .not. l_d_lum ) WRITE (nprint,4001) d_lum

!-----------------------------------------------------------------------
!        Compute the luminosities at density d_lum.
!-----------------------------------------------------------------------

IF ( l_d_lum ) THEN

  area_jd              = frpi * r(jd  ) * r(jd  )
  area_jdm1            = frpi * r(jd-1) * r(jd-1)
  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN
      flux_jd          = zero
      flux_jdm1        = zero
      DO k = 1,nnugp(n)
        flux_jd        = flux_jd   + cvel * psi1(jd,  k,n) * ecoefae(jd,  k) * stwt(n)
        flux_jdm1      = flux_jdm1 + cvel * psi1(jd-1,k,n) * ecoefae(jd-1,k) * stwt(n)
      END DO
    END IF ! nnugp(n) /= 0
    d_lumjd            = flux_jd   * area_jd   * ergfoe
    d_lumjdm1          = flux_jdm1 * area_jdm1 * ergfoe
    lum_rho(n,i_ray)= rinterp( d_lumjd, d_lumjdm1, rho(jd), d_lum, rho(jd-1) )
  END DO

!-----------------------------------------------------------------------
!        End computation of luminosities at d_lum.
!-----------------------------------------------------------------------

END IF ! l_d_lum

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

  rinterp              = b + ( a - b ) * ( y - z )/( x - z )

END FUNCTION rinterp

END SUBROUTINE lum_surface
