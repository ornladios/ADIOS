SUBROUTINE courant_x_time_step( imin, imax, ij_ray_dim, ik_ray_dim, nx, dr, rho, dt, jdt )
!-----------------------------------------------------------------------
!
!    File:         courant_x_time_step
!    Module:       courant_x_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/20/03
!
!    Purpose:
!      To compute the x-sweep time step restricted by the Caurant condition.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  imin         : minimum x-aray index
!  imax         : maximum x-aray index
!  ij_ray_dim, ik_ray_dim    : number of angular rays assigned to a procesor
!  nx           : x-array extent
!  dr           : final zone width (cm)
!  rho          : final density (g cm^{-3})
!
!    Output arguments:
!  dt           : time step given by process i
!  jdt          : (i,ij_ray,ik_ray) for time step given by process i
!
!    Include files:
!  kind_module, numerical_module
!  eos_snc_x_module, prb_cntl.cmn, t_cntrl.cmn
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero

USE eos_snc_x_module, ONLY: gam1, aesv
USE prb_cntl_module, ONLY: ihydro
USE t_cntrl_module, ONLY: tcntrl

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin            ! minimum x-aray index
INTEGER, INTENT(in)              :: imax            ! maximum x-aray index
INTEGER, INTENT(in)              :: ij_ray_dim, ik_ray_dim       ! number of angular rays assigned to a procesor
INTEGER, INTENT(in)              :: nx              ! x-array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)           :: dr   ! final zone width (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: rho  ! final density (g cm^{-3})

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)           :: jdt  ! (i,ij_ray,ik_ray) for time step given by process i

REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim) :: dt   ! time step given by process i

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! radial zone index
INTEGER                          :: j               ! MGFLD radial zone index
INTEGER                          :: ij_ray,ik_ray           ! index denoting a specific radial ray

REAL(KIND=double)                :: d               ! working variable
REAL(KIND=double)                :: dmin            ! minimum value of a set of quantities
REAL(KIND=double), PARAMETER     :: dtmax = 1.d+20  ! times step without criteria
REAL(KIND=double), PARAMETER     :: frthrd = 4.d0/3.d0

!-----------------------------------------------------------------------
!
!        Timestep criteria.
!
!  dt(1,ij_ray,ik_ray)  : Courant time step; tcntrl(1) is the fraction of an x-zone width
!                  a sonic disturbance is allowed to propagate in one time step.
!
!  jdt(1,ij_ray,ik_ray) : zone at which minimum timestep given by the above criterion occurs
!
!  tcntrl(1)    : the fraction of an x-zone width a sonic disturbance
!                  is allowed to propagate in one time step.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize dt(1), and jdt(1)
!-----------------------------------------------------------------------

dt(1,:,:)                 = dtmax
jdt(1,:,:)                = 0

!-----------------------------------------------------------------------
!  Only use the Courant condition if ihydro /= 0
!-----------------------------------------------------------------------

IF ( ihydro /= 0 ) THEN

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim

      dmin                = dtmax

      DO i = imin,imax
  
        j                 = i + 1

!-----------------------------------------------------------------------
!  Find minimum sound crossing time
!-----------------------------------------------------------------------

        IF ( aesv(j,1,ij_ray,ik_ray) > zero ) THEN
          d               = dr(i) * dr(i) * rho(i,ij_ray,ik_ray) &
&                         / ( DMAX1( gam1(j,ij_ray,ik_ray), frthrd ) * aesv(j,1,ij_ray,ik_ray) )
          IF ( d < dmin ) THEN
            dmin          = d
            jdt(1,ij_ray,ik_ray) = j
          END IF ! d < dmin
        END IF ! aesv(j,1,ij_ray,ik_ray) > 0

      END DO ! i = imin,imax

!-----------------------------------------------------------------------
!  Put common Courant condition in dt(1)
!-----------------------------------------------------------------------

      dt(1,ij_ray,ik_ray) = tcntrl(1) * DSQRT(dmin)

    END DO ! ij_ray = 1,ij_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

END IF ! ihydro /= 0

RETURN
END SUBROUTINE courant_x_time_step
