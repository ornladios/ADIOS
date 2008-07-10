SUBROUTINE hydro_x_sweep_time_step( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, dr, rhol, rhoi, tl, ti, nx, dt, jdt, dtime_hydro )
!-----------------------------------------------------------------------
!
!    File:         hydro_x_sweep_time_step
!    Module:       hydro_x_sweep_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To select the time step restricted by matter processes during
!       the x-sweep.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  imin         : minimum x-aray index
!  imax         : maximum x-aray index
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  rhol         : padded density after the x-lagrangian step (g cm^{-3})
!  rhoi         : padded density before the x-lagrangian step (g cm^{-3})
!  tl           : padded temperature after the x-lagrangian step (K)
!  ti           : padded temperature before the x-lagrangian step (K)
!
!    Output arguments:
!  dt           : minimum hydro time step restrictions for criterion i, 
!                  angular zone j (s)
!  jdt          : radial zone causing minimum time step restrictions for
!                  criterion i, angular zone j
!  dtime_hydro  : minimum times step by hydro criteria
!
!    Include files:
!  kind_module, array_module, numerical_module
!  convect_module, eos_snc_x_module, it_tol_module, nu_dist_module,
!  nu_energy_grid_module, prb_cntl.cmn, t_cntrl.cmn
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
Use array_module, ONLY : max_12
USE numerical_module, ONLY: zero

USE convect_module, ONLY: ulcnvct
USE eos_snc_x_module, ONLY: gam1, aesv
USE nu_dist_module, ONLY: psi0
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY: ihydro, ilcnvct
USE t_cntrl_module, ONLY: dtnph, dtst1, ttst1, tcntrl

IMPLICIT none
SAVE
!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                          :: imin       ! minimum x-aray index
INTEGER, INTENT(in)                          :: imax       ! maximum x-aray index
INTEGER, INTENT(in)                          :: ij_ray     ! iindex denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                          :: ik_ray     ! iindex denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                          :: ij_ray_dim ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                          :: ik_ray_dim ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                          :: nx         ! x-array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: dr         ! zone width (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rhol       ! density after the x-lagrangian step (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rhoi       ! density before the x-lagrangian step (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: tl         ! temperature after the x-lagrangian step (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ti         ! temperature before the x-lagrangian step (K)

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)           :: jdt         ! zone causing dtp

REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim) :: dt          ! minimum allowed time step
REAL(KIND=double), INTENT(out)                                      :: dtime_hydro ! minimum times step by hydro criteria

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i                 ! do index
INTEGER                          :: j                 ! radial zone index
INTEGER                          :: jr_min            ! minimum radial zone index
INTEGER                          :: jr_max            ! maximum radial zone index

REAL(KIND=double), PARAMETER     :: dtmax = 1.d+20    ! times step without criteria
REAL(KIND=double), PARAMETER     :: dtmaxi = 1.d-20   ! dtmax^-1
REAL(KIND=double)                :: d                 ! working variable
REAL(KIND=double)                :: dmin              ! minimum value of a set of quantities
REAL(KIND=double)                :: dmax              ! maximum value of a set of quantities
REAL(KIND=double), PARAMETER     :: frthrd = 4.d0/3.d0

REAL(KIND=double)                :: dtnph_t           ! time step restricted by hydrodynamic processes (working variable)

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!  Initialize
!----------------------------------------------------------------------!

jr_min                  = imin + 1
jr_max                  = imax + 1
CALL hydro_time_step_initialize( jr_min, jr_max, dtime_hydro )

!-----------------------------------------------------------------------
!
!        Timestep criteria.
!
!  dt(i,ij_ray,ik_ray)  : minimum timestep given by criterion i
!
!  jdt(i,ij_ray,ik_ray) : zone at which minimum timestep given by criterion i
!                  occurs
!
!  dt(2,ij_ray,ik_ray)  : Density change time step; tcntrl(2) is the maximum
!                  permitted abs( d(rho)/rho ) for tne lagrangian hydro
!                  x-sweep..
!
!  dt(3,ij_ray,ik_ray)  : Temperature change time step due to hydro; tcntrl(3)
!                  is the maximum permitted abs( dtmpmn(j,1,j_ray)/t(j) )
!                  for tne  lagrangian hydro x-sweep.
!
!  dt(9,ij_ray,ik_ray)  : Convective time step; tcntrl(9) is the fraction of a
!                  zone width a convective blob is allowed to propagate
!                  in one time step.
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------------!
!  Initialize dt(i,ij_ray,ik_ray), and jdt(i,ij_ray,ik_ray)
!----------------------------------------------------------------------!

DO i = 2,9
  dt(i,ij_ray,ik_ray)   = dtmax
  jdt(i,ij_ray,ik_ray)  = 0
END DO

!-----------------------------------------------------------------------
!        If jdt(50,ij_ray,ik_ray) = -1, time step is set equal to tcntrl(50);
!         otherwise tcntrl(50) is the maximum allowable time step.
!-----------------------------------------------------------------------

IF ( jdt(50,ij_ray,ik_ray) == -1 ) THEN
  dtime_hydro           = tcntrl(50)
  RETURN
END IF

!-----------------------------------------------------------------------
!
!            \\\\\ DENSITY CHANGE TIME STEP CONTROL /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step restriction due to hydro : dt(2,j_ray)
!
!  tcntrl(2) is the maximum permitted
!     abs( d(rho)/rho )
!   in the hydro x-sweep.
!  Density change time step control is used if ihydro > 0 and tcntrl(2) > 0
!-----------------------------------------------------------------------

IF ( ihydro > 0  .and.  tcntrl(2) > zero ) THEN

  dmax                  = zero

  DO j = jr_min,jr_max

!----------------------------------------------------------------------!
!  Use density time step control if  rhol(j) >= dtst1
!----------------------------------------------------------------------!

    IF ( rhol(j) >= dtst1 ) THEN

!----------------------------------------------------------------------!
!  Get maximum relative change of density
!----------------------------------------------------------------------!

      d                 = DABS( ( rhol(j) - rhoi(j) )/rhol(j) )

      IF ( d > dmax ) THEN
        dmax            = d
        jdt(2,ij_ray,ik_ray) = j
      END IF ! d > dmax
    END IF ! rhol(j) > dtst1

  END DO

!-----------------------------------------------------------------------
!  Use density change time step control only if
!     dmax = max( abs( rhol(j) - rhoi(j) )/rhol(j) ) > 0.
!-----------------------------------------------------------------------

  IF ( dmax > zero ) THEN
    dt(2,ij_ray,ik_ray) = tcntrl(2) * dtnph/dmax
  END IF ! dmax > zero
  
END IF ! ihydro > 0  and  tcntrl(2) > 0

!-----------------------------------------------------------------------
!
!          \\\\\ TEMPERATURE CHANGE TIME STEP CONTROL /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step restriction due to hydro : dt(3,j_ray)
!
!  tcntrl(3) is the maximum permitted 
!     abs( tl(j,ij_ray,ik_ray) - ti(j,ij_ray,ik_ray) )/ti(j,ij_ray,ik_ray).
!  Temperature change time step control is bypassed if tcntrl(3) < 0.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature change time step control is used if tcntrl(3) > 0
!-----------------------------------------------------------------------

IF ( tcntrl(3) > zero ) THEN

dmax                    = zero

  DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Use temperature time step control if tl(j) >= ttst1
!-----------------------------------------------------------------------

    IF ( tl(j) >= ttst1 ) THEN

!........Get maximum relative change of temperature

      d                 = DABS( ( tl(j) - ti(j) )/tl(j) )
      IF ( d > dmax ) THEN
        dmax            = d
        jdt(3,ij_ray,ik_ray) = j
      END IF ! d > dmax
    END IF ! t(j) ge ttst1

  END DO

!-----------------------------------------------------------------------
!  Temperature change time step control due to hydro only if
!     dmax = max( abs( dtmpmn(j,1,ij_ray,ik_ray) )/t(j) ) > 0.
!-----------------------------------------------------------------------

  IF ( dmax > zero ) then
    dt(3,ij_ray,ik_ray) = tcntrl(3) * dtnph/dmax
  END IF ! dmax > zero

END IF ! tcntrl(3) > zero

!-----------------------------------------------------------------------
!  Convective time step ( dt(9,ij_ray,ik_ray) )
!
!  tcntrl(9) is the fraction of a zone width a convective blob is allowed
!   to propagate in one time step.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Convection crossing time step control is used if ilcnvct > 0 and
!   tcntrl(9) > 0
!-----------------------------------------------------------------------

IF ( ilcnvct > 0  .and.  tcntrl(9) > zero ) THEN
  dmin                  = dtmax

!........Get minimum convective crossing time step

  DO j = jr_min,jr_max-1
    IF ( ulcnvct(j) > zero ) THEN
      d                 = dr(j)/ulcnvct(j)
      IF ( d < dmin ) THEN
        dmin            = d
        jdt(9,ij_ray,ik_ray) = j
      END IF ! d < dmin
      d              = dr(j+1)/ulcnvct(j)
      if ( d < dmin ) THEN
        dmin         = d
        jdt(9,ij_ray,ik_ray) = j+1
      END IF ! d < dmin
    END IF ! ulcnvct(j) > 0.0
  END DO

!-----------------------------------------------------------------------
!  dmin is the minimum time for a convective blob to cross a radial zone.
!-----------------------------------------------------------------------

  dt(9,ij_ray,ik_ray)   = tcntrl(9) * dmin
  
END IF ! ilcnvct > 0  and  tcntrl(9) > 0

!-----------------------------------------------------------------------
!  Select minimum time step
!-----------------------------------------------------------------------

dtnph_t                 = tcntrl(50)

dtnph_t                 = DMIN1( dtnph_t, dt(2,ij_ray,ik_ray) )
dtnph_t                 = DMIN1( dtnph_t, dt(3,ij_ray,ik_ray) )
dtnph_t                 = DMIN1( dtnph_t, dt(9,ij_ray,ik_ray) )

dtime_hydro             = dtnph_t

RETURN
END SUBROUTINE hydro_x_sweep_time_step
