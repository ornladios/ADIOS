SUBROUTINE nu_stress_y( jmin, jmax, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
& i_radial, j_radial, nx, ny, nez, nnu, rho, r, rhobar, y, psi0, nu_strs_cy, &
& nu_strs_ey )
!-----------------------------------------------------------------------
!
!    File:         nu_stress_y
!    Module:       nu_stress_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/29/99
!
!    Purpose:
!      To compute the neutrino stress (force per unit mass)
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jmin                : minimum angular zone to calculate neutrino stress
!  jmax                : maximum angular zone to calculate neutrino stress
!  ji_ray              : i (radial) index of a specific angular ray
!  jk_ray              : k (azimuthal) index of a specific angular ray
!  j_ray_dim           : the number of radial zones on a processor after swapping with y
!  ik_ray_dim          : the number of z-zones on a processor before swapping with z
!  i_radial            : the unshifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!  j_radial            : the shifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!  nx                  : x-array extent
!  ny                  : y-array extent
!  nez                 : neutrino energy array extent
!  nnu                 : neutrino flavor array extent
!  rho                 : density (g cm^{-3})
!  r                   : unshifted radius (cm)
!  y                   : y-coordinate
!  rhobar              : angularly averaged density [g cm^{-3}]
!  psi0                : zero angular moments of the neutrino occupation number
!
!    Output arguments:
!  nu_strs_cy          : y-component of zone-center neutrino stress
!  nu_strs_ey          : y-component of zone-edge neutrino stress
!
!    Output arguments (common):
!  strsnu_y(j,k,n)     : force per unit mass in radial zone j - 1/2 due to
!                         n-type neutrinos of group k (dynes/gram)
!  stress_y(j,n,ji_ray,jk_ray) : force per unit mass in radial zone j - 1/2 due to
!                         n-type neutrinos (dynes/gram)
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  e_advct_module, edit_module, nu_dist_module, nu_energy_grid_module,
!  parallel_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, third, half, frpi, frpith
USE physcnst_module, ONLY : pi

USE e_advct_module, ONLY : rhomin_y_eadvect
USE edit_module, ONLY : nprint, nlog
USE nu_dist_module, ONLY : stwt, ecoefa, strsnu_y, stress_y
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY : inutrn

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jmin          ! minimum radial zone index
INTEGER, INTENT(in)              :: jmax          ! maximum radial zone index
INTEGER, INTENT(in)              :: ji_ray        ! i (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray        ! k (azimuthal) index of a specific angular ray
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ji_ray, jk_ray
INTEGER, INTENT(in)              :: j_radial      ! the shifted radial zone corresponding to ji_ray, jk_ray
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(ny)         :: rho    ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)       :: r      ! unshifted radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: rhobar ! angularly averaged density [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)       :: y      ! y coordinate
REAL(KIND=double), INTENT(in), DIMENSION(ny,nez,nnu) :: psi0   ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(ny)   :: nu_strs_cy  ! y-component of neutrino stress
REAL(KIND=double), INTENT(out), DIMENSION(ny+1) :: nu_strs_ey  ! y-component of neutrino stress

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: rimh          ! radius at interface i
REAL(KIND=double)                :: riv2          ! volume centered radius at zone i
REAL(KIND=double)                :: ri            ! mean radius radius at zone i
REAL(KIND=double)                :: riph          ! radius at interface i+1
REAL(KIND=double)                :: strskj        ! stress_x at j, energy k, flavor n
REAL(KIND=double)                :: p_difference  ! pressure difference

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: u_nu    ! neutrino energy density
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: u_nupad ! padded neutrino energy density
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dvol    ! zone volume (cm^{3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dmrst   ! zone mass (g)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: area    ! zone interface area (cm^{2})

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' jmax=',i4,' > ny_max=',i4,' in subroutine nu_stress_y')
 1001 FORMAT (' Allocation problem for array ',a10,' in nu_stress_y')
 2001 FORMAT (' Deallocation problem for array ',a10,' in nu_stress_y')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Zero stresses, return if nnugpmx = 0 or inutrn = 0
!-----------------------------------------------------------------------

stress_y(:,:,ji_ray,jk_ray) = zero
strsnu_y                  = zero
nu_strs_cy                = zero
nu_strs_ey                = zero

IF ( nnugpmx == 0  .or.  inutrn == 0 ) RETURN

!-----------------------------------------------------------------------
!  Return if rhobar < rhomin_y_eadvect
!-----------------------------------------------------------------------

IF ( rhobar(i_radial) < rhomin_y_eadvect ) RETURN

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (u_nu(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u_nupad(ny+2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nupad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmrst(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmrst     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (area(ny+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'area      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

u_nu                      = zero
u_nupad                   = zero
dvol                      = zero
dmrst                     = zero
area                      = zero

!-----------------------------------------------------------------------
!
!                  \\\\\ COMPUTE THE STRESS /////
!
!
!  Zone Centered Stress, stress:
!
!                             1
!               (P    - P   ) - (A      + A     )
!                 j+1    j-1  2   j-1/2    j+1/2
!     stress  = ---------------------------------
!           j                dM
!                              j
!
!  Zone Edged Stress:
!
!                              
!                    (P    - P   ) A
!                      j+1    j-1   j+1/2
!     stress     = - --------------------
!           j+1/2     1/2 ( dM  + dM   )
!                           J     j+1
!
!  Zone Centered Stress, stress:
!
!              1
!              - (stress      A     + stress      A     )
!              2        j+1/2  j+1/2        j-1/2  j-1/2
!     stress  = -----------------------------------------
!           j                   dM
!                                 j
!  where
!
!           2
!     dV = r  ( r     - r     ) (cos theta      - cos theta      ) d phi 
!       j   i    i+1/2   i-1/2             j+1/2            j-1/2       k
!
!     A      = r  sin theta      ( r     - r     ) 
!      j+1/2    i          j+1/2    i+1/2   i-1/2
!
!     dM = dV  rho
!       j    j    ijk
!
!                               sin theta
!                                        j+1/2
!     ==> A   /dM = ------------------------------------------
!          j+1   j   r (cos theta      - cos theta      ) rho
!                     i          j+1/2            j-1/2      ijk
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Radii
!-----------------------------------------------------------------------

rimh                      = r(i_radial)
riph                      = r(i_radial+1)
riv2                      = third * ( riph * ( riph + rimh ) + rimh * rimh )
ri                        = half * ( riph + rimh )

!-----------------------------------------------------------------------
!  Areas, volumes, and mass
!-----------------------------------------------------------------------

area(jmin:jmax+1)         = 2.d0 * ri * DSIN(y(jmin:jmax+1))
dvol(jmin:jmax)           = 4.d0 * riv2 * ( DCOS(y(jmin:jmax)) - DCOS(y(jmin+1:jmax+1)) )
dmrst(jmin:jmax)          = dvol(jmin:jmax) * rho(jmin:jmax)

!-----------------------------------------------------------------------
!  Neutrino energy densities and zone-edge stress
!-----------------------------------------------------------------------

DO j = jmin,jmax
  DO k = 1,nnugpmx
    DO n = 1,nnu
      if ( nnugp(n) == 0 ) CYCLE
      u_nu(j,k)           = ecoefa(j_radial,k) * psi0(j,k,n) * stwt(n)
      u_nupad(j+1)        = u_nupad(j+1) + u_nu(j,k)
    END DO ! n = 1,nnu
  END DO ! k = 1,nnugpmx
END DO ! j = jmin,jmax

DO j = jmin+1,jmax
  DO k = 1,nnugpmx
    DO n = 1,nnu
      if ( nnugp(n) == 0 ) CYCLE
      p_difference        = third * ( u_nu(j,k) - u_nu(j-1,k) )
      strskj              = - area(j) * p_difference &
&                         / ( half * ( dmrst(j) + dmrst(j-1) ) )
      strsnu_y(j,k,n)     = strskj
      stress_y(j,n,ji_ray,jk_ray) = stress_y(j,n,ji_ray,jk_ray) + strskj
    END DO !  = 1,nnu
  END DO ! k = 1,nnugpmx
END DO ! j = jmin+1,jmax

!DO j = jmin+1,jmax
!  DO n = 1,nnu
!    nu_strs_ey(j)         = nu_strs_ey(j) + stress_y(j,n,ji_ray,jk_ray)
!  END DO !  = 1,nnu
!END DO ! j = jmin+1,jmax

nu_strs_ey                = zero

!-----------------------------------------------------------------------
!  Zone-centered stress
!-----------------------------------------------------------------------

u_nupad(jmin)             = u_nupad(jmin+1)
u_nupad(jmax+2)           = u_nupad(jmax+1)
!DO j = jmin,jmax
!  p_difference            = third * ( u_nupad(j+2) - u_nupad(j) )
!  nu_strs_cy(j)           = - p_difference * half * ( area(j) + area(j+1) )/dmrst(j)
!END DO ! j
nu_strs_cy                = zero

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (u_nu, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nu      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u_nupad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nupad   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dmrst, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmrst     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (area, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'area      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE nu_stress_y
