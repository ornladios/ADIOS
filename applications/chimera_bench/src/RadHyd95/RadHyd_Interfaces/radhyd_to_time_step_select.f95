SUBROUTINE radyhd_to_time_step_select( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, ndim )
!-----------------------------------------------------------------------
!
!    File:         radyhd_to_time_step_select
!    Module:       radyhd_to_time_step_select
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/05/05
!
!    Purpose:
!      To select the minimum time step as given by all time step restrictions
!       on a given processor.
!
!    Input arguments:
!  nx         : x-array extent
!  ny         : y-array extent
!  nz         : z-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  j_ray_dim  : the number of radial zones on a processor after swapping with y (used in the MPI version)
!  k_ray_dim  : the number of radial zones on a processor after swapping with z (used in the MPI version)
!  ndim       : number of spatial dimensions of the simulation
!
!    Output arguments:
!        none
!
!    Subprograms called:
!      time_step_ray_select : selects the minimum time stgeo as given by all
!                              time step restrictions on a given processor.
!
!    Include files:
!  kind_module, array_module
!  angular_ray_module, azimuthal_ray_module, cycle_module, edit_module,
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc

USE angular_ray_module, ONLY : dt_y, jdt_y, dt_y_state
USE azimuthal_ray_module, ONLY : dt_z, jdt_z, dt_z_state
USE cycle_module, ONLY: ncynu_trns, intnu_trns, nutrans_trns
USE edit_module, ONLY : nlog, nprint
USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, dtnmh,     &
& dtnph, time, dx_cf, rho_c, dt, jdt, dt_process, j_radial_dt, j_angular_dt, &
& j_azimuthal_dt, dtnph_trans, dtime_trans

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: nx              ! x-array extent
INTEGER, INTENT(in)                    :: ny              ! y-array extent
INTEGER, INTENT(in)                    :: nz              ! z-array extent
INTEGER, INTENT(in)                    :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)                    :: ik_ray_dim      ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)                    :: j_ray_dim       ! the number of radial zones on a processor after swapping with y (used in the MPI version)
INTEGER, INTENT(in)                    :: k_ray_dim       ! the number of radial zones on a processor after swapping with z (used in the MPI version)
INTEGER, INTENT(in)                    :: ndim            ! number of spatial dimensions of the simulation

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: i               ! radial zone index
INTEGER                                :: i_ray           ! saved radial ray index at which minimum dt occurs
INTEGER                                :: j_ray           ! saved angular ray index at which minimum dt occurs
INTEGER                                :: k_ray           ! saved azimuthal ray index at which minimum dt occurs
INTEGER                                :: ij_ray          ! y (angular) index of a specific x (radial) ray
INTEGER                                :: ik_ray          ! z (azimuthal) index of a specific x (radial) ray
INTEGER                                :: ji_ray          ! x (radial) index of a specific y (angular) ray
INTEGER                                :: jk_ray          ! z (azimuthal) index of a specific y (angular) ray
INTEGER                                :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER                                :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray
INTEGER, DIMENSION(50,jmax,kmax)       :: jdt_all         ! jdt_all(i,ij_ray,ik_ray) is the radial zone giving the minimum dt for 
!                                                            process i for radial ray with y-index ij_ray and z-index ik_ray
INTEGER, DIMENSION(3,imax,kmax)        :: jdt_y_all       ! jdt_all(i,ji_ray,jk_ray) is the y (angular) zone giving the minimum dt for 
!                                                            process i for angular ray x-index ji_ray and z-index jk_ray
INTEGER, DIMENSION(3,jmax,imax)        :: jdt_z_all       ! jdt_all(i,kj_ray,ki_ray) is the z (azimuthal) zone giving the minimum dt for 
!                                                            process i for angular ray x-index ki_ray and y-index kj_ray

REAL(KIND=double), PARAMETER           :: dtmax = 1.d+20  ! times step without criteria
REAL(KIND=double), DIMENSION(50,jmax,kmax) :: dt_all      ! dt_all(i,ij_ray,ik_ray) is the minimum dt for process i for radial ray
!                                                             with y-index ij_ray and z-index ik_ray
REAL(KIND=double), DIMENSION(3,imax,kmax) :: dt_y_all     ! dt_y_all(i,ji_ray,jk_ray) is the minimum dt for 
!                                                            process i on y (angular) ray with x-index ji_ray and z-index jkray
REAL(KIND=double), DIMENSION(3,jmax,imax) :: dt_z_all     ! dt_y_all(i,ji_ray,jk_ray) is the minimum dt for 
!                                                            process i on y (angular) ray with x-index ji_ray and z-index jkray
REAL(KIND=double)                      :: dt_test         ! parameter to determine minimum timestep
REAL(KIND=double), DIMENSION(jmax,kmax)    :: dtime_trans_all ! source and transport time step given by radial ray i_ray

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Need MPI version of radyhd_to_time_step_select since n_proc=',i4,' > 1')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Stop if n_proc > 1.
!-----------------------------------------------------------------------

IF ( n_proc > 1 ) THEN
  WRITE (nprint,101) n_proc
  WRITE (nlog,101) n_proc
  STOP
END IF ! n_proc > 1

!-----------------------------------------------------------------------
!  Advance the time and switch dtnph to dtnmh
!-----------------------------------------------------------------------

CALL time_adv( dtnph, dtnmh, time )

!-----------------------------------------------------------------------
!  Compute Caurant time step
!-----------------------------------------------------------------------

CALL courant_x_time_step( imin, imax, ij_ray_dim, ik_ray_dim, nx, dx_cf, &
& rho_c, dt, jdt )

!-----------------------------------------------------------------------
!  Maximum time step increase restriction and maximum time step criterion
!-----------------------------------------------------------------------

CALL time_step_limits( dtnmh, ij_ray_dim, ik_ray_dim, dt, jdt )

!-----------------------------------------------------------------------
!  Select minimum time step
!-----------------------------------------------------------------------

!........Initialize

dt_all                    = dtmax
jdt_all                   = 0

!-----------------------------------------------------------------------
!  Information in dt and jdt is the minimum time step (dt) and the radial
!   zone (jdt) where the minimum time step occurs  as a function of the
!   y-index ij_ray and the z-index ik_ray of the radial ray. Here we make
!   a direct transfer to dt_all and jdt_all. In the MPI version, this
!   transfer would have to follow an MPI_Gather.
!-----------------------------------------------------------------------

DO i = 1,3
  dt_all (i,:,:)          = dt (i,:,:)
  jdt_all(i,:,:)          = jdt(i,:,:)
END DO

!-----------------------------------------------------------------------
!  Information in dt_y and jdt_y is the minimum time step (dt_y) and the
!   angular zone (jdt_y) where the minimum time step occurs as a function
!   of the x-index, ji_ray, and the z index, jk_ray, of the angular ray
!-----------------------------------------------------------------------

IF ( ndim > 1  .and.  ny > 1 ) THEN

  DO i = 1,3
    DO jk_ray = kmin,kmax
      DO ji_ray = imin,imax
        dt_y_all (i,ji_ray,jk_ray) = dt_y (i,ji_ray,jk_ray)
        jdt_y_all(i,ji_ray,jk_ray) = jdt_y(i,ji_ray,jk_ray)
      END DO ! ji_ray
    END DO ! jk_ray
  END DO ! i

!-----------------------------------------------------------------------
!  The information in dt_y_all and jdt_y_all, i = 1 - 3, is transferred
!   to dt_all and jdt_all, i = 4 - 6. Thus, dt_all(i,j_ray,k_ray) gets
!   the minimum time for process i and the angular zone, j_ray, and
!   azimuthal zone, k_ray, that it occurs; the rest of the array gets
!   set above to dtmax.
!-----------------------------------------------------------------------

  DO i = 4,6
    dt_test                 = dtmax
    i_ray                   = 1
    j_ray                   = 1
    k_ray                   = 1
    DO jk_ray = kmin,kmax
      DO ji_ray = imin,imax
        IF ( dt_test > dt_y_all(i-3,ji_ray,jk_ray) ) THEN
          dt_test           = dt_y_all (i-3,ji_ray,jk_ray)
          j_ray             = jdt_y_all(i-3,ji_ray,jk_ray)
          i_ray             = ji_ray
          k_ray             = jk_ray
        END IF ! dt_test > dt_y_all(i-3,ji_ray)
      END DO ! ji_ray = imin,imax
    END DO ! jk_ray = kmin,kmax
    IF ( j_ray /= 0 ) THEN
      dt_all (i,j_ray,k_ray) = dt_test
      jdt_all(i,j_ray,k_ray) = i_ray
    END IF ! j_ray /= 0
  END DO ! i = 4,6

END IF ! ndim > 1  .and.  ny > 1

!-----------------------------------------------------------------------
!  Information in dt_z and jdt_z is the minimum time step (dt_z) and the
!   azimuthal zone (jdt_z) where the minimum time step occurs as a
!   function of the x-index, ki_ray, and the j index, kj_ray, of the
!   azimuthal ray
!-----------------------------------------------------------------------

IF ( ndim > 1  .and.  nz > 1 ) THEN

  DO i = 1,3
    DO kj_ray = jmin,jmax
      DO ki_ray = imin,imax
        dt_z_all (i,kj_ray,ki_ray) = dt_z (i,kj_ray,ki_ray)
        jdt_z_all(i,kj_ray,ki_ray) = jdt_z(i,kj_ray,ki_ray)
      END DO ! ji_ray
    END DO ! jk_ray
  END DO ! i

!-----------------------------------------------------------------------
!  The information in dt_z_all and jdt_z_all, i = 1 - 3, is transferred
!   to dt_all and jdt_all, i = 4 - 6. Thus, dt_all(i,j_ray,k_ray) gets
!   the minimum time for process i and the angular zone, j_ray, and
!   azimuthal zone, k_ray, that it occurs; the rest of the array gets
!   set above to dtmax.
!-----------------------------------------------------------------------

  DO i = 7,9
    dt_test                 = dtmax
    i_ray                   = 1
    j_ray                   = 1
    k_ray                   = 1
    DO kj_ray = jmin,jmax
      DO ki_ray = imin,imax
        IF ( dt_test > dt_z_all(i-6,kj_ray,ki_ray) ) THEN
          dt_test           = dt_z_all (i-6,kj_ray,ki_ray)
          k_ray             = jdt_z_all(i-6,kj_ray,ki_ray)
          i_ray             = ki_ray
          j_ray             = kj_ray
        END IF
      END DO ! ji_ray = imin,imax
    END DO ! kj_ray = jmin,jmax
    IF ( k_ray /= 0 ) THEN
      dt_all (i,j_ray,k_ray)  = dt_test
      jdt_all(i,j_ray,k_ray)  = i_ray
    END IF ! j_ray /= 0
  END DO ! i = 7,9

END IF ! ndim > 1  .and.  nz > 1

!-----------------------------------------------------------------------
!  Information in dt and jdt is the minimum time step (dt) and the radial
!   zone (jdt) where the minimum time step occurs  as a function of the
!   radial ray with y-index ij_ray and z-index ik_ray. Here we make a
!   direct transfer to dt_all and jdt_all. In the MPI version, this
!   transfer would have to follow an MPI_Gather.
!-----------------------------------------------------------------------

DO i = 10,50
  dt_all (i,:,:)            = dt (i,:,:)
  jdt_all(i,:,:)            = jdt(i,:,:)
END DO

!-----------------------------------------------------------------------
!  Load dtime_trans_all
!-----------------------------------------------------------------------

dtime_trans_all             = dtmax

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    dtime_trans_all(ij_ray,ik_ray) = dtime_trans(ij_ray,ik_ray)
  END DO ! ij_ray
END DO ! ik_ray

!-----------------------------------------------------------------------
!  Select minimum time step
!-----------------------------------------------------------------------

  CALL time_step_select( ny, nz, dt_all, jdt_all, dt_process, j_radial_dt, &
&  j_angular_dt, j_azimuthal_dt, dtnph, dtnph_trans, dtime_trans_all,      &
&  ncynu_trns, intnu_trns, nutrans_trns, dt_y_state, dt_z_state )

!-----------------------------------------------------------------------
!  Select minimum time step
!-----------------------------------------------------------------------

CALL time_step_brdcst( dt_process, j_radial_dt, j_angular_dt, j_azimuthal_dt, &
& dtnmh, dtnph, dtnph_trans)

!-----------------------------------------------------------------------
!  Reinitialize dt_y and jdt_y, and dt_z and jdt_z
!-----------------------------------------------------------------------

dt_y                        = dtmax
jdt_y                       = 0
dt_z                        = dtmax
jdt_z                       = 0

RETURN
END SUBROUTINE radyhd_to_time_step_select
