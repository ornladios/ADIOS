SUBROUTINE time_step_limits( dtnmh, ij_ray_dim, ik_ray_dim, dt, jdt )
!-----------------------------------------------------------------------
!
!    File:         time_step_limits
!    Module:       time_step_limits
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To set limits on the time step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  dtnmh      : time step for the cycle just elapsed
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!
!    Output arguments:
!  dt         : time step given by process i
!  jdt        : (i,ij_ray,ik_ray) for time step given by process i
!
!    Include files:
!  kind_module, numerical_module
!  t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero

USE t_cntrl_module, ONLY: tcntrl

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                     :: ij_ray_dim       ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)                                     :: ik_ray_dim       ! number of z-zones on a processor before swapping with z

REAL(KIND=double), INTENT(inout)                        :: dtnmh            ! time step for the preceding cycle

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)           :: jdt  ! (i,i_ray) for time step given by process i

REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim) :: dt   ! time step given by process i

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)                                       :: dt_10     ! time step given by maximum allowed increase

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Restrict inccrease in time step
!
!  dt(10, i_ray) : Time step given by the maximum allowed increase in 
!                     hydro time step
!  tcntrl(10)    : the maximum allowed value of dtnph/dtnmh.
!-----------------------------------------------------------------------

dt_10              = tcntrl(10) * dtnmh
dt (10,:,:)        = dt_10
jdt(10,:,:)        = 0

!-----------------------------------------------------------------------
!        Set maximum tine step
!
!  tcntrl(50) : Maximum allowed time step.
!-----------------------------------------------------------------------

dt (50,:,:)        = tcntrl(50)
jdt(50,:,:)        = 0



RETURN
END SUBROUTINE time_step_limits
