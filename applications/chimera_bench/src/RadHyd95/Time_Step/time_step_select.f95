SUBROUTINE time_step_select( ny, nz, dt_all, jdt_all, dt_process,    &
& j_radial_dt, j_angular_dt, j_azimuthal_dt, dtnph, dtnph_trans,     &
& dtime_trans_all, ncynu_trns, intnu_trns, nutrans_trns, dt_y_state, &
& dt_z_state )
!-----------------------------------------------------------------------
!
!    File:         time_step_select
!    Module:       time_step_select
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/20/05
!
!    Purpose:
!      To select the time step restricted by all processes.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ny              : y-array extent
!  nz              : z-array extent
!  dt_all          : time step restricted by all criteria on all radial rays
!  jdt_all         : zone responsible for each dt_all
!  dtime_trans_all : source and transport time step as a function of radial ray i_ray
!  intnu_trns      : maximum hydro cycles between transport steps
!
!    Output arguments:
!  dt_process      : minimum time step for process i
!  j_radial_dt     : radial ray responsible for each dt_process
!  j_angular_dt    : angular ray responsible for each dt_process
!  j_azimuthal_dt  : azimuthal ray responsible for each dt_process
!  dtnph           : time step restricted by all criteria
!  dtnph_trans     : new source and transport time step
!  nutrans_trns    : transport step key
!  ncynu_trns      : number of hydro cycles since the last transport step
!  dt_y_state      : y-hydro time step limited by state variable change criteria
!
!    Include files:
!  kind_module, numerical_module
!  t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero

USE t_cntrl_module, ONLY: dtnph_t=>dtnph, t_step_xyz

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                             :: ny              ! y-array extent
INTEGER, INTENT(in)                             :: nz              ! z-array extent
INTEGER, INTENT(in), DIMENSION(50,ny,nz)        :: jdt_all         ! jdt_all(i,ij_ray,ik_ray) is the radial zone giving the minimum dt for 
!                                                                     process i for radial ray with y-index ij_ray and z-index ik_ray
INTEGER, INTENT(in)                             :: intnu_trns      ! maximum hydro cycles between transport steps

REAL(KIND=double), INTENT(in), DIMENSION(50,ny,nz) :: dt_all       ! dt_all(i,ij_ray,ik_ray) is the minimum dt for process i for radial ray
!                                                                     with y-index ij_ray and z-index ik_ray

REAL(KIND=double), INTENT(in), DIMENSION(ny,nz) :: dtime_trans_all ! source and transport time step given by radial ray i_ray

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out)                            :: ncynu_trns      ! number of hydro cycles since the last transport step
INTEGER, INTENT(out), DIMENSION(50)             :: j_radial_dt     ! radial ray responsible for each dt_process
INTEGER, INTENT(out), DIMENSION(50)             :: j_angular_dt    ! angular ray responsible for each dt_process
INTEGER, INTENT(out), DIMENSION(50)             :: j_azimuthal_dt  ! azimuthal ray responsible for each dt_process

REAL(KIND=double), INTENT(out), DIMENSION(50)   :: dt_process      ! time step restricted by all criteria
REAL(KIND=double), INTENT(out)                  :: dtnph           ! new time step restricted by all criteria
REAL(KIND=double), INTENT(out)                  :: dtnph_trans     ! new source and transport time step
REAL(KIND=double), INTENT(out)                  :: dt_y_state      ! y-hydro time step limited by state variable change criteria
REAL(KIND=double), INTENT(out)                  :: dt_z_state      ! z-hydro time step limited by state variable change criteria

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

LOGICAL, INTENT(inout)                          :: nutrans_trns    ! transport step key

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                         :: i               ! process index
INTEGER                                         :: ij_ray          ! y-index denoting a specific radial ray
INTEGER                                         :: ik_ray          ! z-index denoting a specific radial ray

REAL(KIND=double), PARAMETER                    :: dtmax = 1.d+20  ! times step without criteria
REAL(KIND=double)                               :: dt_sum          ! accumulated time since last source and transport step

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                      \\\\\ INITIALIZE /////
!
!-----------------------------------------------------------------------

dtnph                     = dtmax
dt_process                = dtmax
dt_y_state                = dtmax
dt_z_state                = dtmax

!-----------------------------------------------------------------------
!
!      \\\\\ DT__PROCESS(I), D_RADIAL_DT(I), J_ANGULAR_DT(I) /////
!
!  dt_process(i)   : time step restricted by criterion i
!  j_radial_dt(i)  : radial zone restricting time step by criterion i
!  j_angular_dt(i) : angular zone restricting time step by criterion i
!
!-----------------------------------------------------------------------

DO i = 1,50
  DO ik_ray = 1,nz
    DO ij_ray = 1,ny
      IF ( dt_all(i,ij_ray,ik_ray) < dt_process(i) ) THEN
        dt_process(i)     = dt_all(i,ij_ray,ik_ray)
        j_radial_dt(i)    = jdt_all(i,ij_ray,ik_ray)
        j_angular_dt(i)   = ij_ray
        j_azimuthal_dt(i) = ik_ray
      END IF ! dt_all(i,i_ray) < dt_process(i)
    END DO ! ij_ray = 1,ny
  END DO ! ik_ray = 1,ny
END DO ! i = 1,50

!-----------------------------------------------------------------------
!
!                          \\\\\ DTNPH /////
!
!  dtnph with t_step_xyz = ye : minimum time step as restricted by all
!           criterion
!  dtnph with t_step_xyz = no : minimum time step as restricted by all
!           criterion excludingthe y-hydro time step restrictions.
!           The latter are handled by subcycling.
!
!-----------------------------------------------------------------------

IF ( t_step_xyz == 'ye' ) THEN
  dtnph                   = MINVAL( dt_process(1:50) )
ELSE
  dtnph                   = DMIN1( MINVAL( dt_process(1:3)), MINVAL( dt_process(10:50) ) )
  dt_y_state              = MINVAL( dt_process(5:6))
  dt_z_state              = MINVAL( dt_process(8:9))
END IF

!-----------------------------------------------------------------------
!
!                      \\\\\ DTNPHN_TRANS /////
!
!  dtnph_trans : time step for the source and transport step
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  If nutrans_trns = true, source and transport step were performed
!   during the current cycle and dtnph_trans must be determined for the
!   next cycle. Reset dt_sum, nutrans_trns, and ncynu_trns.
!
!-----------------------------------------------------------------------

IF ( nutrans_trns ) THEN
  dtnph_trans             = DMIN1( dtmax, MINVAL( dtime_trans_all(1:ny,1:nz) ) )
  dt_sum                  = zero
  nutrans_trns            = .false.
  ncynu_trns              = 0
END IF ! nutrans_trns == .true.

!-----------------------------------------------------------------------
!
!  dt_sum is the source and transport time step if executed on the next
!   cycle.
!  ncynu_trns : the number of hydro steps since the last source and
!   transport step.
!
!-----------------------------------------------------------------------

dt_sum                    = dt_sum + dtnph
ncynu_trns                = ncynu_trns + 1

!-----------------------------------------------------------------------
!
!  If ncynu_trns >= intnu_trns, execute source and transport step during
!   the next cycle with dtnph_trans equal to the accumulated time from
!   the previous implementation of source and transport.
!
!-----------------------------------------------------------------------

IF ( ncynu_trns >= intnu_trns ) THEN
  nutrans_trns            = .true.
  dtnph_trans             = dt_sum
END IF


!-----------------------------------------------------------------------
!
!  If ncynu_trns < intnu_trns, test whether skipping the next time would
!   result in dt_sum > dtnph_trans
!   
!-----------------------------------------------------------------------

IF ( ncynu_trns < intnu_trns ) THEN
  IF ( dtnph_trans < dt_sum + dtnph ) THEN
    nutrans_trns          = .true.
    dtnph_trans           = dt_sum
  END IF ! dtnph_trans < dt_sum
END IF ! ncynu_trns < intnu_trns


RETURN
END SUBROUTINE time_step_select
