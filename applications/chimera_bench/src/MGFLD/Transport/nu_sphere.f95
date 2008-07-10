SUBROUTINE nu_sphere( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& r, rho, t, rstmss, nx, nez, nnu, j_sphere, r_sphere, d_sphere, t_sphere, &
& m_sphere )
!-----------------------------------------------------------------------
!
!    File:         nu_sphere
!    Module:       nu_sphere
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/22/00
!
!    Purpose:
!      To calculate the radius, density, temperature, and enclosed
!       mass of the neutrino-spheres.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min              : inner radial zone index
!  jr_max              : outer radial zone index
!  ij_ray              : j-index of a radial ray
!  ik_ray              : k-index of a radial ray
!  ij_ray_dim          : number of y-zones on a processor before swapping
!  ik_ray_dim          : number of z-zones on a processor before swapping
!  r                   : radius (cm)
!  rho                 : density (g cm^{-3}
!  t                   : temperature (K)
!  rstmss              : rest mass (g)
!  nx                  : x_array extent
!  nez                 : neutrino energy array extent
!  nnu                 : neutrino flavor array extent
!
!    Output arguments:
!  j_sphere(k,n,ij_ray,ik_ray) : radial zone outwardly contiguous to the n-neutrinosphere for energy zone k
!  r_sphere(k,n,ij_ray,ik_ray) : radius of the n-neutrinosphere for energy zone k (cm)
!  d_sphere(k,n,ij_ray,ik_ray) : density of the n-neutrinosphere for energy zone k (g cm^{-3})
!  t_sphere(k,n,ij_ray,ik_ray) : temperature of the n-neutrinosphere for energy zone k (MeV)
!  m_sphere(k,n,ij_ray,ik_ray) : enclosed rest mass of the n-neutrinosphere for energy zone k (g)
!
!    Output arguments (common):
!        none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : kmev

USE nu_dist_module, ONLY : rhs1, gamgr_nu
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nx            ! radial array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent


REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: r        ! radial zone radii (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: rho      ! shifted density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: t        ! shifted temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: rstmss   ! enclosed rest mass (g)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim)           :: j_sphere   ! zone index of the neutrinospheres

REAL(KIND=double), INTENT(out), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: r_sphere   ! neutrinosphere radius
REAL(KIND=double), INTENT(out), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: d_sphere   ! neutrinosphere density
REAL(KIND=double), INTENT(out), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: t_sphere   ! neutrinosphere temperature
REAL(KIND=double), INTENT(out), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: m_sphere   ! neutrinosphere enclosed mass


!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index

REAL(KIND=double)                :: taut          ! transport optical depth to bottom of zone j
REAL(KIND=double)                :: tautp         ! transport optical depth to bottom of zone j+1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ N-NEUTRINO EMISSION SUEFACES /////
!
!-----------------------------------------------------------------------

CALL gamgr_nu_cal(jr_min,jr_max)

DO n = 1,nnu

!-----------------------------------------------------------------------
!  Set nutrinospheres to zero if nnugp(n) == 0
!-----------------------------------------------------------------------

  IF ( nnugp(n) == 0 ) THEN
    j_sphere(:,n,ij_ray,ik_ray) = 1
    r_sphere(:,n,ij_ray,ik_ray) = zero
    d_sphere(:,n,ij_ray,ik_ray) = zero
    t_sphere(:,n,ij_ray,ik_ray) = zero
    m_sphere(:,n,ij_ray,ik_ray) = zero
    CYCLE
  END IF

!-----------------------------------------------------------------------
!  Compute the nutrinospheres for each neutrino energy and flavor
!-----------------------------------------------------------------------

  DO k = 1,nnugp(n)

    taut                    = zero

    DO j = jr_max,jr_min,-1
      tautp                 = taut
      taut                  = taut - rhs1(j,k,n) * ( r(j) - r(j-1) )/( half * ( gamgr_nu(j) + gamgr_nu(j-1) ) )
      IF ( taut >= one ) THEN

!-----------------------------------------------------------------------
!  Interpolate r to taut = 1
!-----------------------------------------------------------------------

        j_sphere(k,n,ij_ray,ik_ray) = j
        r_sphere(k,n,ij_ray,ik_ray) = r(j)      + ( r  (j-1)     - r  (j) )     * ( one - tautp )/( taut - tautp + epsilon )
        d_sphere(k,n,ij_ray,ik_ray) = rho(j)    + ( rho(j-1)     - rho(j) )     * ( one - tautp )/( taut - tautp + epsilon )
        t_sphere(k,n,ij_ray,ik_ray) = t(j)      + ( t  (j-1)     - t  (j) )     * ( one - tautp )/( taut - tautp + epsilon )
        m_sphere(k,n,ij_ray,ik_ray) = rstmss(j) + ( rstmss (j-1) - rstmss (j) ) * ( one - tautp )/( taut - tautp + epsilon )
        EXIT

      END IF ! taut > 1

!-----------------------------------------------------------------------
!  IF taut < 1 at center, use central values
!-----------------------------------------------------------------------

      IF ( j == jr_min ) THEN
        j_sphere(k,n,ij_ray,ik_ray) = jr_min
        r_sphere(k,n,ij_ray,ik_ray) = r(jr_min-1)
        d_sphere(k,n,ij_ray,ik_ray) = rho(jr_min)
        t_sphere(k,n,ij_ray,ik_ray) = t(jr_min)
        m_sphere(k,n,ij_ray,ik_ray) = rstmss(jr_min-1)
      END IF ! j = jr_min

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Convert temperatures to MeV
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  t_sphere(:,n,ij_ray,ik_ray) = t_sphere(:,n,ij_ray,ik_ray) * kmev
END DO

RETURN
END SUBROUTINE nu_sphere
