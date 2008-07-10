SUBROUTINE hydro_z_sweep_state_time_step( kmin, kmax, ki_ray, kj_ray, &
& ij_ray_dim, k_ray_dim, rhol, rhoi, tl, ti, nz, dt_z, jdt_z )
!-----------------------------------------------------------------------
!
!    File:         hydro_z_sweep_state_time_step
!    Module:       hydro_z_sweep_state_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/18/07
!
!    Purpose:
!      To select the time step restricted by density and temperature change
!       criteria during the z-sweep.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  kmin         : minimum y-aray index
!  kmax         : maximum y-aray index
!  ki_ray       : x (radial) index of a specific z (azimuthal) ray
!  kj_ray       : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim   : number of y (angular) zones on a processor before swapping with y
!  k_ray_dim    : number of x (radial) zones on a processor after swapping with z
!  rhol         : padded density after the y-lagrangian step (g cm^{-3})
!  rhoi         : padded density before the y-lagrangian step (g cm^{-3})
!  tl           : padded temperature after the y-lagrangian step (K)
!  ti           : padded temperature before the y-lagrangian step (K)
!  nz           : y-array extent
!
!    Output arguments:
!  dt_z         : minimum hydro time step restrictions for criterion i, angular zone k (s)
!  jdt_z        : radial zone causing minimum time step restrictions for criterion i, angular zone k
!
!    Include files:
!  kind_module, array_module, numerical_module
!  convect_module, nu_dist_module, nu_energy_grid_module, prb_cntl.cmn,
!  t_cntrl.cmn
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
Use array_module, ONLY : max_12
USE numerical_module, ONLY: zero

USE convect_module, ONLY: ulcnvct
USE nu_dist_module, ONLY: psi0
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY: ihydro, ilcnvct
USE t_cntrl_module, ONLY: dtnph, dtst1, ttst1, tcntrl

IMPLICIT none
SAVE
!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                          :: kmin       ! minimum y-aray index
INTEGER, INTENT(in)                          :: kmax       ! maximum y-aray index
INTEGER, INTENT(in)                          :: ki_ray     ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                          :: kj_ray     ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                          :: ij_ray_dim ! number of y (angular) zones on a processor before swapping with y
INTEGER, INTENT(in)                          :: k_ray_dim  ! number of x (radial) zones on a processor after swapping with z
INTEGER, INTENT(in)                          :: nz         ! z-array extent

REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rhol       ! density after the z-lagrangian step (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rhoi       ! density before the z-lagrangian step (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: tl         ! temperature after the z-lagrangian step (K)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: ti         ! temperature before the z-lagrangian step (K)

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(3,ij_ray_dim,k_ray_dim)           :: jdt_z       ! zone causing dtp

REAL(KIND=double), INTENT(out), DIMENSION(3,ij_ray_dim,k_ray_dim) :: dt_z        ! minimum allowed time step

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i                 ! do index
INTEGER                          :: k                 ! radial zone index

REAL(KIND=double), PARAMETER     :: dtmax = 1.d+20    ! times step without criteria
REAL(KIND=double), PARAMETER     :: dtmaxi = 1.d-20   ! dtmax^-1
REAL(KIND=double)                :: d                 ! working variable
REAL(KIND=double)                :: dmax              ! maximum value of a set of quantities

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!        Timestep criteria.
!
!  dt_z(i,kj_ray,ki_ray)  : minimum timestep given by criterion i
!
!  jdt_z(i,kj_ray,ki_ray) : zone at which minimum timestep given by
!   criterion i occurs
!
!  dt(2,kj_ray,ki_ray)    : Density change time step; tcntrl(8) is
!   the maximum permitted abs( d(rho)/rho ) for tne lagrangian hydro
!   z-sweep.
!
!  dt(3,kj_ray,ki_ray)    : Temperature change time step due to hydro;
!   tcntrl(9) is the maximum permitted abs( dtmpmn(k,2,kj_ray,ki_ray)/t(k) )
!   for the lagrangian hydro z-sweep.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize dt_z(i,kj_ray,ki_ray), and jdt_z(i,kj_ray,ki_ray)
!-----------------------------------------------------------------------

DO i = 2,3
  dt_z(i,kj_ray,ki_ray)  = dtmax
  jdt_z(i,kj_ray,ki_ray) = 0
END DO

!-----------------------------------------------------------------------
!
!            \\\\\ DENSITY CHANGE TIME STEP CONTROL /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Density change time step restriction due to hydro : dt_z(2,kj_ray,ki_ray)
!
!  tcntrl(8) is the maximum permitted
!     abs( d(rho)/rho )
!   in the hydro y-sweep.
!  Density change time step control is used if ihydro > 0 and tcntrl(8) > 0
!-----------------------------------------------------------------------

IF ( ihydro > 0  .and.  tcntrl(8) > zero ) THEN

  dmax                 = zero

  DO k = kmin,kmax

!----------------------------------------------------------------------!
!  Use density time step control if  rhol(k) >= dtst1
!----------------------------------------------------------------------!

    IF ( rhol(k) >= dtst1 ) THEN

!----------------------------------------------------------------------!
!  Get maximum relative change of density
!----------------------------------------------------------------------!

      d                = DABS( ( rhol(k) - rhoi(k) )/rhol(k) )

      IF ( d > dmax ) THEN
        dmax           = d
        jdt_z(2,kj_ray,ki_ray) = k
      END IF ! d > dmax
    END IF ! rhol(k) > dtst1

  END DO

!-----------------------------------------------------------------------
!  Use density change time step control only if
!     dmax = max( abs( rhol(k) - rhoi(k) )/rhol(k) ) > 0.
!-----------------------------------------------------------------------

  IF ( dmax > zero ) THEN
    dt_z(2,kj_ray,ki_ray)    = tcntrl(8) * dtnph/dmax    
  END IF ! dmax > zero
  
END IF ! ihydro > 0  and  tcntrl(8) > 0

!-----------------------------------------------------------------------
!
!          \\\\\ TEMPERATURE CHANGE TIME STEP CONTROL /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step restriction due to hydro : dt_z(3,kj_ray,ki_ray)
!
!  tcntrl(9) is the maximum permitted 
!     abs( tl(k,kj_ray,ki_ray) - ti(k,i_ray) )/ti(k,kj_ray,ki_ray).
!  Temperature change time step control is bypassed if tcntrl(9) < 0.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step control is used if tcntrl(9) > 0
!-----------------------------------------------------------------------

IF ( tcntrl(9) > zero ) THEN

dmax                 = zero

  DO k = kmin,kmax

!-----------------------------------------------------------------------
!  Use temperature time step control if tl(k) >= ttst1
!-----------------------------------------------------------------------

    IF ( tl(k) >= ttst1 ) THEN

!........Get maximum relative change of temperature

      d              = DABS( ( tl(k) - ti(k) )/tl(k) )
      IF ( d > dmax ) THEN
        dmax         = d
        jdt_z(3,kj_ray,ki_ray) = k
      END IF ! d > dmax
    END IF ! t(k) ge ttst1

  END DO

!-----------------------------------------------------------------------
!  Temperature change time step control due to hydro only if
!     dmax = max( abs( tl(k,kj_ray,ki_ray) - ti(k,i_ray) )/ti(k,kj_ray,ki_ray) ) > 0.
!-----------------------------------------------------------------------

  IF ( dmax > zero ) then
    dt_z(3,kj_ray,ki_ray)    = tcntrl(9) * dtnph/dmax
  END IF ! dmax > zero

END IF ! tcntrl(9) > zero

RETURN
END SUBROUTINE hydro_z_sweep_state_time_step