SUBROUTINE nu_stress_z( kmin, kmax, ki_ray, kj_ray, ij_ray_dim, k_ray_dim, &
& i_radial, j_radial, nx, ny, nz, nez, nnu, rho, r, rhobar, y, z, psi0,    &
& nu_strs_cz, nu_strs_ez )
!-----------------------------------------------------------------------
!
!    File:         nu_stress_z
!    Module:       nu_stress_z
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
!  kmin                : minimum z (azimuthal) zone to calculate neutrino stress
!  kmax                : maximum z (azimuthal) zone to calculate neutrino stress
!  ki_ray              : x (radial) index of a specific z (azimuthal) ray
!  kj_ray              : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim          : the number of y-zones on a processor before swapping with y
!  k_ray_dim           : the number of radial zones on a processor after swapping with z
!  i_radial            : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  j_radial            : the shifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nx                  : x-array extent
!  ny                  : y-array extent
!  nz                  : z-array extent
!  nez                 : neutrino energy array extent
!  nnu                 : neutrino flavor array extent
!  rho                 : density (g cm^{-3})
!  r                   : unshifted radius (cm)
!  y                   : y (angular) coordinate
!  z                   : z (azimuthal) coordinate
!  rhobar              : angularly averaged density [g cm^{-3}]
!  psi0                : zero angular moments of the neutrino occupation number
!
!    Output arguments:
!  nu_strs_cz          : z-component of zone-center neutrino stress
!  nu_strs_ez          : z-component of zone-edge neutrino stress
!
!    Output arguments (common):
!  strsnu_z(iz,k,n)     : force per unit mass in radial zone iz - 1/2 due to
!                         n-type neutrinos of group k (dynes/gram)
!  stress_z(iz,n,kj_ray,ki_ray) : force per unit mass in radial zone iz - 1/2 due to
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

USE e_advct_module, ONLY : rhomin_z_eadvect
USE edit_module, ONLY : nprint, nlog
USE nu_dist_module, ONLY : stwt, ecoefa, strsnu_z, stress_z
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY : inutrn

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: kmin          ! minimum radial zone index
INTEGER, INTENT(in)              :: kmax          ! maximum radial zone index
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! the number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! the number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ki_ray, kj_ray
INTEGER, INTENT(in)              :: j_radial      ! the shifted radial zone corresponding to ki_ray, kj_ray
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nz)         :: rho    ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)       :: r      ! unshifted radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: rhobar ! angularly averaged density [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)       :: y      ! y (angular) coordinate
REAL(KIND=double), INTENT(in), DIMENSION(nz+1)       :: z      ! z (azimuthal) coordinate
REAL(KIND=double), INTENT(in), DIMENSION(nz,nez,nnu) :: psi0   ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nz)   :: nu_strs_cz  ! z-component of neutrino stress
REAL(KIND=double), INTENT(out), DIMENSION(nz+1) :: nu_strs_ez  ! z-component of neutrino stress

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

CHARACTER (len=10)               :: var_name

INTEGER                          :: iz            ! z-zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: rimh          ! radius at interface i
REAL(KIND=double)                :: riv2          ! volume centered radius at zone i
REAL(KIND=double)                :: ri            ! mean radius radius at zone i
REAL(KIND=double)                :: riph          ! radius at interface i+1
REAL(KIND=double)                :: strskj        ! stress_x at iz, energy k, flavor n
REAL(KIND=double)                :: p_difference  ! pressure difference

REAL(KIND=double), DIMENSION(300)              :: costheta ! zone interface area (cm^{2})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: u_nu     ! neutrino energy density
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: u_nupad  ! padded neutrino energy density
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dvol     ! zone volume (cm^{3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dmrst    ! zone mass (g)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: area     ! zone interface area (cm^{2})

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' kmax=',i4,' > ny_max=',i4,' in subroutine nu_stress_z')
 1001 FORMAT (' Allocation problem for array ',a10,' in nu_stress_z')
 2001 FORMAT (' Deallocation problem for array ',a10,' in nu_stress_z')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Zero stresses, return if nnugpmx = 0 or inutrn = 0
!-----------------------------------------------------------------------

stress_z(:,:,kj_ray,ki_ray) = zero
strsnu_z                  = zero
nu_strs_cz                = zero
nu_strs_ez                = zero

IF ( nnugpmx == 0  .or.  inutrn == 0 ) RETURN

!-----------------------------------------------------------------------
!  Return if rhobar < rhomin_z_eadvect
!-----------------------------------------------------------------------

IF ( rhobar(i_radial) < rhomin_z_eadvect ) RETURN

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (u_nu(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u_nupad(nz+2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nupad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmrst(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmrst     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (area(nz+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'area      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

u_nu                      = zero
u_nupad                   = zero
dvol                      = zero
dmrst                     = zero
area                      = zero

IF ( first ) THEN
  costheta(1:ny+1)        = DCOS(y(1:ny+1))
  first                   = .false.
END IF ! first

!-----------------------------------------------------------------------
!
!                  \\\\\ COMPUTE THE STRESS /////
!
!
!  Zone Centered Stress, stress:
!
!                             1
!               (P    - P   ) - (A      + A     )
!                 iz+1    iz-1  2   iz-1/2    iz+1/2
!     stress  = ---------------------------------
!           iz                dM
!                              iz
!
!  Zone Edged Stress:
!
!                              
!                    (P    - P   ) A
!                      iz+1    iz-1   iz+1/2
!     stress     = - --------------------
!           iz+1/2     1/2 ( dM  + dM   )
!                           J     iz+1
!
!  Zone Centered Stress, stress:
!
!              1
!              - (stress      A     + stress      A     )
!              2        iz+1/2  iz+1/2        iz-1/2  iz-1/2
!     stress  = -----------------------------------------
!           iz                   dM
!                                 iz
!  where
!
!           2
!     dV = r  ( r     - r     ) (cos theta      - cos theta      ) d phi 
!       iz   i    i+1/2   i-1/2             iz+1/2            iz-1/2       k
!
!     A      = r  sin theta      ( r     - r     ) 
!      iz+1/2    i          iz+1/2    i+1/2   i-1/2
!
!     dM = dV  rho
!       iz    iz    ijk
!
!                               sin theta
!                                        iz+1/2
!     ==> A   /dM = ------------------------------------------
!          iz+1   iz   r (cos theta      - cos theta      ) rho
!                     i          iz+1/2            iz-1/2      ijk
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

area(kmin:kmax+1)         = 4.d0 * ri * ( y(kj_ray+1) - y(kj_ray) )
dvol(kmin:kmax)           = riv2 * ( - costheta(kj_ray+1) + costheta(kj_ray) ) &
&                         * ( z(kmin+1:kmax+1) - z(kmin:kmax) )
dmrst(kmin:kmax)          = dvol(kmin:kmax) * rho(kmin:kmax)

!-----------------------------------------------------------------------
!  Neutrino energy densities and zone-edge stress
!-----------------------------------------------------------------------

DO iz = kmin,kmax
  DO k = 1,nnugpmx
    DO n = 1,nnu
      if ( nnugp(n) == 0 ) CYCLE
      u_nu(iz,k)          = ecoefa(j_radial,k) * psi0(iz,k,n) * stwt(n)
      u_nupad(iz+1)       = u_nupad(iz+1) + u_nu(iz,k)
    END DO ! n = 1,nnu
  END DO ! k = 1,nnugpmx
END DO ! iz = kmin,kmax

DO iz = kmin+1,kmax
  DO k = 1,nnugpmx
    DO n = 1,nnu
      if ( nnugp(n) == 0 ) CYCLE
      p_difference        = third * ( u_nu(iz,k) - u_nu(iz-1,k) )
      strskj              = - area(iz) * p_difference &
&                         / ( half * ( dmrst(iz) + dmrst(iz-1) ) )
      strsnu_z(iz,k,n)    = strskj
      stress_z(iz,n,kj_ray,ki_ray) = stress_z(iz,n,kj_ray,ki_ray) + strskj
    END DO ! n = 1,nnu
  END DO ! k = 1,nnugpmx
END DO ! iz = kmin+1,kmax

!DO iz = kmin+1,kmax
!  nu_strs_ez(iz)          = nu_strs_ez(iz) + stress_z(iz,n,kj_ray,ki_ray)
!END DO ! iz = kmin+1,kmax
nu_strs_ez                = zero

!-----------------------------------------------------------------------
!  Zone-centered stress
!-----------------------------------------------------------------------

u_nupad(kmin)             = u_nupad(kmin+1)
u_nupad(kmax+2)           = u_nupad(kmax+1)
!DO iz = kmin,kmax
!  p_difference            = third * ( u_nupad(iz+2) - u_nupad(iz) )
!  nu_strs_cz(iz)          = - p_difference * half * ( area(iz) + area(iz+1) )/dmrst(iz)
!END DO ! iz
nu_strs_cz                = zero

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
END SUBROUTINE nu_stress_z
