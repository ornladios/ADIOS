SUBROUTINE hydro_y_sweep_state_time_step( jmin, jmax, ji_ray, jk_ray, &
& j_ray_dim,ik_ray_dim, rhol, rhoi, tl, ti, ny, dt_y, jdt_y )
!-----------------------------------------------------------------------
!
!    File:         hydro_y_sweep_state_time_step
!    Module:       hydro_y_sweep_state_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To select the time step restricted by density and temperature change
!       criteria during the y-sweep.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jmin         : minimum y-aray index
!  jmax         : maximum y-aray index
!  ji_ray       : x (radial) index of a specific y (angular) ray
!  jk_ray       : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim    : the number of radial zones on a processor after swapping with y
!  ik_ray_dim   : the number of z-zones on a processor before swapping with z
!  rhol         : padded density after the y-lagrangian step (g cm^{-3})
!  rhoi         : padded density before the y-lagrangian step (g cm^{-3})
!  tl           : padded temperature after the y-lagrangian step (K)
!  ti           : padded temperature before the y-lagrangian step (K)
!  ny           : y-array extent
!
!    Output arguments:
!  dt_y         : minimum hydro time step restrictions for criterion i, angular zone j (s)
!  jdt_y        : radial zone causing minimum time step restrictions for criterion i, angular zone j
!
!    Include files:
!  kind_module, array_module, numerical_module
!  convect_module, nu_dist_module, nu_energy_grid_module, prb_cntl.cmn,
!  t_cntrl.cmn
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

INTEGER, INTENT(in)                          :: jmin       ! minimum y-aray index
INTEGER, INTENT(in)                          :: jmax       ! maximum y-aray index
INTEGER, INTENT(in)                          :: ji_ray     ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)                          :: jk_ray     ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)                          :: j_ray_dim  ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)                          :: ik_ray_dim ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)                          :: ny         ! y-array extent

REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rhol       ! density after the y-lagrangian step (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rhoi       ! density before the y-lagrangian step (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: tl         ! temperature after the y-lagrangian step (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: ti         ! temperature before the y-lagrangian step (K)

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(3,j_ray_dim,ik_ray_dim)           :: jdt_y       ! zone causing dtp

REAL(KIND=double), INTENT(out), DIMENSION(3,j_ray_dim,ik_ray_dim) :: dt_y        ! minimum allowed time step

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i                 ! do index
INTEGER                          :: j                 ! y (angular) zone index

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
!  dt_y(i,ji_ray,jk_ray)  : minimum timestep given by criterion i
!
!  jdt_y(i,ji_ray,jk_ray) : zone at which minimum timestep given by
!   criterion i occurs
!
!  dt(2,ji_ray,jk_ray)    : Density change time step; tcntrl(5) is the
!   maximum permitted abs( d(rho)/rho ) for tne lagrangian hydro y-sweep.
!
!  dt(3,ji_ray,jk_ray)    : Temperature change time step due to hydro;
!   tcntrl(6) is the maximum permitted abs( dtmpmn(j,2,ji_ray,jk_ray)/t(j) )
!   for the lagrangian hydro y-sweep.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize dt_y(i,ji_ray,jk_ray), and jdt_y(i,ji_ray,jk_ray)
!-----------------------------------------------------------------------

DO i = 2,3
  dt_y(i,ji_ray,jk_ray)  = dtmax
  jdt_y(i,ji_ray,jk_ray) = 0
END DO

!-----------------------------------------------------------------------
!
!            \\\\\ DENSITY CHANGE TIME STEP CONTROL /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step restriction due to hydro : dt_y(2,ji_ray,jk_ray)
!
!  tcntrl(5) is the maximum permitted
!     abs( d(rho)/rho )
!   in the hydro y-sweep.
!  Density change time step control is used if ihydro > 0 and tcntrl(5) > 0
!-----------------------------------------------------------------------

IF ( ihydro > 0  .and.  tcntrl(5) > zero ) THEN

  dmax                 = zero

  DO j = jmin,jmax

!-----------------------------------------------------------------------
!  Use density time step control if  rhol(j) >= dtst1
!-----------------------------------------------------------------------

    IF ( rhol(j) >= dtst1 ) THEN

!-----------------------------------------------------------------------
!  Get maximum relative change of density
!-----------------------------------------------------------------------

      d                = DABS( ( rhol(j) - rhoi(j) )/rhol(j) )

      IF ( d > dmax ) THEN
        dmax           = d
        jdt_y(2,ji_ray,jk_ray) = j
      END IF ! d > dmax
    END IF ! rhol(j) > dtst1

  END DO

!-----------------------------------------------------------------------
!  Use density change time step control only if
!     dmax = max( abs( rhol(j) - rhoi(j) )/rhol(j) ) > 0.
!-----------------------------------------------------------------------

  IF ( dmax > zero ) THEN
    dt_y(2,ji_ray,jk_ray) = tcntrl(5) * dtnph/dmax    
  END IF ! dmax > zero
  
END IF ! ihydro > 0  and  tcntrl(5) > 0

!-----------------------------------------------------------------------
!
!          \\\\\ TEMPERATURE CHANGE TIME STEP CONTROL /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step restriction due to hydro : dt_y(3,ji_ray,jk_ray)
!
!  tcntrl(6) is the maximum permitted 
!     abs( tl(j,ji_ray,jk_ray) - ti(j,i_ray) )/ti(j,ji_ray,jk_ray).
!  Temperature change time step control is bypassed if tcntrl(6) < 0.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step control is used if tcntrl(6) > 0
!-----------------------------------------------------------------------

IF ( tcntrl(6) > zero ) THEN

dmax                 = zero

  DO j = jmin,jmax

!-----------------------------------------------------------------------
!  Use temperature time step control if tl(j) >= ttst1
!-----------------------------------------------------------------------

    IF ( tl(j) >= ttst1 ) THEN

!........Get maximum relative change of temperature

      d              = DABS( ( tl(j) - ti(j) )/tl(j) )
      IF ( d > dmax ) THEN
        dmax         = d
        jdt_y(3,ji_ray,jk_ray) = j
      END IF ! d > dmax
    END IF ! t(j) ge ttst1

  END DO

!-----------------------------------------------------------------------
!  Temperature change time step control due to hydro only if
!     dmax = max( abs( tl(j,ji_ray,jk_ray) - ti(j, tcntrl(6)) )/ti(j,ji_ray,jk_ray) ) > 0.
!-----------------------------------------------------------------------

  IF ( dmax > zero ) then
    dt_y(3,ji_ray,jk_ray) = tcntrl(6) * dtnph/dmax
  END IF ! dmax > zero

END IF ! tcntrl(6) > zero

RETURN
END SUBROUTINE hydro_y_sweep_state_time_step