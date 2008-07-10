SUBROUTINE nu_stress_x( jr_min, jr_max, ij_ray, ik_ray, nx, nez, nnu, rho, &
& rhobar, r, rhs1, dc, psi0, nu_strs_c, nu_strs_e )
!-----------------------------------------------------------------------
!
!    File:         nu_stress_x
!    Module:       nu_stress_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/29/99
!
!    Purpose:
!      To compute the x-componsnt of neutrino stress (force per unit mass)
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min               : minimum radial zone to calculate neutrino stress
!  jr_max               : maximum radial zone to calculate neutrino stress
!  ij_ray               : j-index of a radial ray
!  ik_ray               : k-index of a radial ray
!  nx                   : x-array extent
!  nez                  : neutrino energy array extent
!  nnu                  : neutrino flavor array extent
!  rho                  : shifted density (g cm^{-3})
!  rhobar               : angularly averaged density [g cm^{-3}]
!  r                    : radius (cm)
!  rhs1                 : shifted right-hand side of transport equation, first moment
!  dc                   : the diffusion coefficient for neutrinos
!  psi0                 : shifted zero moment of the neutrino occupation probability
!
!    Output arguments:
!  nu_strs_c            : zone-centered x-component of neutrino stress
!  nu_strs_e            : zone-edge x-component of neutrino stress
!
!    Output arguments (common):
!  strsnu_x(j,k,n)      : force per unit mass in radial zone j - 1/2 due to
!                         n-type neutrinos of group k (dynes/gram)
!  stress_ex(j,n,ij_ray,ik_ray)
!                       : force per unit mass in radial zone j - 1/2 due to
!                         n-type neutrinos (dynes/gram)
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  e_advct_module, edit_module, nu_dist_module, nu_energy_grid_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, third, half, frpi
USE physcnst_module, ONLY : pi

USE e_advct_module, ONLY : rhomin_y_eadvect
USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : stwt, ecoefa, ecoefae, strsnu_x, stress_ex=>stress_x
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : idiff, inutrn

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: rhobar ! angularly averaged density [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: r      ! radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu)    :: dc     ! the diffusion coefficient for neutrinos

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx)    :: nu_strs_c     ! zone-centered x-component of neutrino stress
REAL(KIND=double), INTENT(out), DIMENSION(nx)    :: nu_strs_e     ! zone-edge x-component of neutrino stress

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu) :: rhs1   ! the right-hand side of transport equation, first moment
REAL(KIND=double), INTENT(inout), DIMENSION(nx)         :: rho    ! shifted density (g cm^{-3})
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu) :: psi0   ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                               :: var_name

LOGICAL                                          :: first = .true.

INTEGER                                          :: j             ! radial zone index
INTEGER                                          :: k             ! neutrino energy index
INTEGER                                          :: n             ! neutrino flavor index
INTEGER                                          :: istat         ! allocation status
INTEGER                                          :: j_tautrns     ! zone at which tautrns = tautrns_edge
INTEGER                                          :: j_nonu        ! the first shifted index for which rhobar < rhomin_y_eadvect

REAL(KIND=double)                                :: rjmaxmh       ! zone-centered radius inside of boundary
REAL(KIND=double)                                :: rjmaxph       ! zone-centered radius outside of boundary
REAL(KIND=double)                                :: psi0_ratio    ! ratio of psi0 across boundary

REAL(KIND=double)                                :: rhoj          ! zone-centered radius
REAL(KIND=double), SAVE                          :: pi8           ! 8*pi
REAL(KIND=double)                                :: zltotjmh      ! inverse transport mean free path at j-1/2
REAL(KIND=double)                                :: zltotjph      ! inverse transport mean free path at j+1/2
REAL(KIND=double)                                :: zltotj        ! inverse transport mean free path at j
REAL(KIND=double)                                :: psi1it        ! psi1 for computing stress_ex
REAL(KIND=double)                                :: strskj_e      ! zone edge stress_ex at j, energy k, flavor n
REAL(KIND=double)                                :: strskj_c      ! zone center stress_ex at j, energy k, flavor n
REAL(KIND=double), PARAMETER                     :: tautrns_surf = 1.d0

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: u_nu          ! shifted neutrino energy density (ergs cm^{-3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: tautrns       ! transport optical depth
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol          ! shifted zone volume (cm^{3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dmrst         ! shifted zone mass (g)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: area          ! zone interface area (cm^{2})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rjmh2         ! shifted zone centered radius squared (cm^{2})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: stress_cx     ! total shifted neutrino energy density (ergs cm^{-3})

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in nu_stress_x')
 2001 FORMAT (' Deallocation problem for array ',a10,' in nu_stress_x')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Zero stresses, return if nnugpmx = 0 or inutrn = 0
!-----------------------------------------------------------------------

stress_ex(:,:,ij_ray,ik_ray)  = zero
strsnu_x                      = zero
nu_strs_e                     = zero
nu_strs_c                     = zero

IF ( nnugpmx == 0  .or.  inutrn == 0 ) RETURN

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (u_nu(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nu      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tautrns(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tautrns   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmrst(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmrst     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (area(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'area      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rjmh2(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rjmh2     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (stress_cx(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'stress_cx '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

u_nu                          = zero
tautrns                       = zero
dvol                          = zero
dmrst                         = zero
area                          = zero
rjmh2                         = zero

IF ( first ) THEN
  pi8                         = 8.d+00 * pi
  first                       = .false.
END IF

stress_cx(:,:)                = zero
rho(jr_max+1)                 = rho(jr_max)

!-----------------------------------------------------------------------
!  Zone-centered radii, areas volumes, and masses
!-----------------------------------------------------------------------

rjmh2(jr_min:jr_max)          = third * ( r(jr_min:jr_max)                        &
&                             * ( r(jr_min:jr_max) + r(jr_min-1:jr_max-1) )       &
&                             + r(jr_min-1:jr_max-1) * r(jr_min-1:jr_max-1) )
area(jr_min-1:jr_max)         = r(jr_min-1:jr_max) * r(jr_min-1:jr_max)
dvol(jr_min:jr_max)           = rjmh2(jr_min:jr_max) * ( r(jr_min:jr_max) - r(jr_min-1:jr_max-1) )
dmrst(jr_min:jr_max)          = dvol(jr_min:jr_max) * rho(jr_min:jr_max)
dmrst(1)                      = zero
dmrst(jr_max+1)               = dmrst(jr_max)

!-----------------------------------------------------------------------
!  Boundary values of psi0
!-----------------------------------------------------------------------

rjmaxmh                       = half * ( r(jr_max) + r(jr_max-1) )
rjmaxph                       = r(jr_max) + half * ( r(jr_max) - r(jr_max-1) )
psi0_ratio                    = ( rjmaxmh * rjmaxmh )/( rjmaxph * rjmaxph )
psi0(jr_max+1,:,:)            = psi0_ratio * psi0(jr_max,:,:)


!-----------------------------------------------------------------------
!
!            \\\\\ COMPUTE THE TRANSPORT OPTICAL DEPTH /////
!
!-----------------------------------------------------------------------

tautrns                       = zero
DO j = jr_max,jr_min,-1
  tautrns(j,:,:)              = tautrns(j+1,:,:) - rhs1(j,:,:) * ( r(j) - r(j-1) )
END DO ! j

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE THE NEUTRINO ENERGY DENSITIES /////
!
!-----------------------------------------------------------------------

DO n = 1,nnu
  u_nu(jr_min:jr_max,:,n)     = ecoefa(jr_min:jr_max,:) * psi0(jr_min:jr_max,:,n) * stwt(n)
END DO
u_nu(1,:,:)                   = u_nu(2,:,:)
u_nu(jr_max+1,:,:)            = u_nu(jr_max,:,:)

!-----------------------------------------------------------------------
!
!                  \\\\\ COMPUTE THE STRESS /////
!
!-----------------------------------------------------------------------



SELECT CASE (idiff)

  CASE(:3,5)

    DO n = 1,nnu
      IF ( nnugp(n) == 0 ) CYCLE
      DO k = 1,nnugp(n)
        rhs1(jr_max+1,k,n)     = rhs1(jr_max,k,n)

!-----------------------------------------------------------------------
!  Compute j_tautrns
!-----------------------------------------------------------------------

        j_tautrns              = jr_min
        DO j = jr_max,jr_min,-1
          IF ( tautrns(j,k,n) > tautrns_surf ) THEN
            j_tautrns          = j
            EXIT
          END IF ! tautrns(j,k,n) > tautrns_surf
        END DO ! j
        j_tautrns              = MIN( j_tautrns, jr_max - 1 )

!-----------------------------------------------------------------------
!  Compute stress for large transport optical depth
!-----------------------------------------------------------------------

        DO j = jr_min,j_tautrns
          strskj_e             = - area(j) * third * ( u_nu(j+1,k,n) - u_nu(j,k,n) ) &
&                              / ( half * ( dmrst(j+1) + dmrst(j) ) )
          strskj_c             = - half * ( area(j) + area(j-1) ) * third * ( u_nu(j+1,k,n) - u_nu(j-1,k,n) ) &
&                              / ( dmrst(j) + half * ( dmrst(j+1) + dmrst(j-1) ) )
          strsnu_x(j,k,n)      = strskj_e
          stress_ex(j,n,ij_ray,ik_ray)                                             &
&                              = stress_ex(j,n,ij_ray,ik_ray) + strskj_e
          stress_cx(j,n)       = stress_cx(j,n)       + strskj_c
        END DO ! j = jr_min,j_tautrns

!-----------------------------------------------------------------------
!  Compute stress for small transport optical depth
!-----------------------------------------------------------------------

        DO j = j_tautrns+1, jr_max
          rhoj                 = ( rho(j) * dmrst(j) + rho(j+1) * dmrst(j+1) )     &
&                              / ( dmrst(j) + dmrst(j+1) )
          zltotjmh             = - rhs1(j  ,k,n)
          zltotjph             = - rhs1(j+1,k,n)
          zltotj               = DSQRT( zltotjmh * zltotjph + 1.d-100 )
          psi1it               = -dc(j,k,n) * area(j) * rhoj                                   &
&                              * ( psi0(j+1,k,n) * ecoefa(j+1,k) - psi0(j,k,n) * ecoefa(j,k) ) &
&                              / ( half * ( dmrst(j) + dmrst(j+1) ) * ecoefae(j,k) )
          strskj_e             = ( ecoefae(j,k) * stwt(n)/rhoj ) * ( zltotj * psi1it )
          strsnu_x(j,k,n)      = strskj_e
          stress_ex(j,n,ij_ray,ik_ray)                                             &
&                              = stress_ex(j,n,ij_ray,ik_ray) + strskj_e
        END DO ! j

        DO j = j_tautrns+1,jr_max
          stress_cx(j,n)       = half * ( stress_ex(j,n,ij_ray,ik_ray) + stress_ex(j-1,n,ij_ray,ik_ray) )
        END DO ! j = j_tautrns+1,jr_max

      END DO ! k = 1, nnugp(n)
    END DO ! n = 1,nnu

  CASE(4)

    DO n = 1, nnu
      if ( nnugp(n) == 0 ) CYCLE
      DO k = 1, nnugp(n)
        DO j = jr_min, jr_max
          strskj_e             = - area(j) * third * ( u_nu(j+1,k,n) - u_nu(j,k,n) ) &
&                              / ( half * ( dmrst(j+1) + dmrst(j) ) )
          strskj_c             = - half * ( area(j) + area(j-1) ) * third * ( u_nu(j+1,k,n) - u_nu(j-1,k,n) ) &
&                              / ( dmrst(j) + half * ( dmrst(j+1) + dmrst(j-1) ) )
          strsnu_x(j,k,n)      = strskj_e
          stress_ex(j,n,ij_ray,ik_ray)                                               &
&                              = stress_ex(j,n,ij_ray,ik_ray) + strskj_e
          stress_cx(j,n)       = stress_cx(j,n)       + strskj_c
        END DO ! j = jr_min, jr_max
      END DO ! k = 1, nnugp(n)
    END DO ! n = 1, nnu

END SELECT

j_nonu                         = MAXLOC( rhobar, DIM = 1, MASK = rhobar < rhomin_y_eadvect ) + 1
DO j = j_nonu,jr_max
  nu_strs_e(j)                 = SUM( stress_ex(j,:,ij_ray,ik_ray) )
  nu_strs_c(j)                 = SUM( stress_cx(j,:) )
END DO ! j = j_nonu,jr_max

nu_strs_e(jr_max)              = nu_strs_e(jr_max-1)
nu_strs_e(jr_max+1)            = nu_strs_e(jr_max-1)
nu_strs_c(jr_max)              = nu_strs_c(jr_max-1)

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (u_nu, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_nu      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (tautrns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tautrns   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dmrst, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmrst     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (area, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'area      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rjmh2, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rjmh2     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (stress_cx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'stress_cx '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE nu_stress_x
