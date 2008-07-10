SUBROUTINE radyhd_to_time_step_select_MPI( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, ndim )
!-----------------------------------------------------------------------
!
!    File:         radyhd_to_time_step_select_MPI
!    Module:       radyhd_to_time_step_select_MPI
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
!  j_ray_dim  : the number of radial zones on a processor after swapping with y
!  k_ray_dim  : the number of radial zones on a processor after swapping with z
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
!  angular_ray_module, azimuthal_ray_module, cycle_module, parallel_module,
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc, n_proc_y, n_proc_z

USE angular_ray_module, ONLY : dt_y, jdt_y, dt_y_state
USE azimuthal_ray_module, ONLY : dt_z, jdt_z, dt_z_state
USE cycle_module, ONLY: ncynu_trns, intnu_trns, nutrans_trns
USE parallel_module, ONLY : myid, ierr
USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, dtnmh, &
& dtnph, time, dx_cf, rho_c, dt, jdt, dt_process, j_radial_dt, j_angular_dt, &
& j_azimuthal_dt, dtnph_trans, dtime_trans

USE mpi

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
INTEGER, INTENT(in)                    :: j_ray_dim       ! the number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)                    :: k_ray_dim       ! the number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)                    :: ndim            ! number of spatial dimensions of the simulation

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: i               ! y-array index
INTEGER                                :: j               ! y-array index
INTEGER                                :: k               ! z-array index
INTEGER                                :: m               ! processor index
INTEGER                                :: mi              ! x-block index
INTEGER                                :: mj              ! y-block index
INTEGER                                :: mk              ! z-block index
INTEGER                                :: isk             ! x-array index of gathered array
INTEGER                                :: jsk             ! y-array index of gathered array
INTEGER                                :: ksk             ! z-array index of gathered array
INTEGER                                :: c_gath_send     ! gather send buffer count
INTEGER                                :: c_gath_recv     ! gather recv buffer count
INTEGER                                :: i_ray           ! saved radial ray index at which minimum dt occurs
INTEGER                                :: j_ray           ! saved angular ray index at which minimum dt occurs
INTEGER                                :: k_ray           ! saved azimuthal ray index at which minimum dt occurs
INTEGER                                :: ij_ray          ! y (angular) index of a specific x (radial) ray
INTEGER                                :: ik_ray          ! z (azimuthal) index of a specific x (radial) ray
INTEGER                                :: ji_ray          ! x (radial) index of a specific y (angular) ray
INTEGER                                :: jk_ray          ! z (azimuthal) index of a specific y (angular) ray
INTEGER                                :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER                                :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray
INTEGER                                :: i_extent        ! broadcast array extent

INTEGER, DIMENSION(50,jmax,kmax)       :: jdt_all         ! jdt_all(i,ij_ray,ik_ray) is the radial zone giving the minimum dt for 
!                                                            process i for radial ray with y-index ij_ray and z-index ik_ray
INTEGER, DIMENSION(50,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z) :: irecv_buf1 ! receive buffer for gathering data from processors
INTEGER, DIMENSION(3,j_ray_dim,ik_ray_dim,n_proc_y*n_proc_z)   :: irecv_buf3 ! receive buffer for gathering data from processors
INTEGER, DIMENSION(3,ij_ray_dim,k_ray_dim,n_proc_y*n_proc_z)   :: irecv_buf4 ! receive buffer for gathering data from processors
INTEGER, DIMENSION(3,imax,kmax)        :: jdt_y_all       ! jdt_all(i,ji_ray,jk_ray) is the y (angular) zone giving the minimum dt for 
!                                                            process i for angular ray x-index ji_ray and z-index jk_ray
INTEGER, DIMENSION(3,jmax,imax)        :: jdt_z_all       ! jdt_all(i,kj_ray,ki_ray) is the z (azimuthal) zone giving the minimum dt for 
!                                                            process i for angular ray x-index ki_ray and y-index kj_ray
INTEGER, DIMENSION(2)                  :: i_send          ! send buffer

REAL(KIND=double), PARAMETER           :: dtmax = 1.d+20  ! times step without criteria
REAL(KIND=double), DIMENSION(50,jmax,kmax) :: dt_all      ! dt_all(i,ij_ray,ik_ray) is the minimum dt for process i for radial ray
!                                                             with y-index ij_ray and z-index ik_ray
REAL(KIND=double), DIMENSION(3,imax,kmax) :: dt_y_all     ! dt_y_all(i,ji_ray,jk_ray) is the minimum dt for 
!                                                            process i on y (angular) ray with x-index ji_ray and z-index jkray
REAL(KIND=double), DIMENSION(3,jmax,imax) :: dt_z_all     ! dt_y_all(i,ji_ray,jk_ray) is the minimum dt for 
!                                                            process i on y (angular) ray with x-index ji_ray and z-index jkray
REAL(KIND=double)                      :: dt_test         ! parameter to determine minimum timestep
REAL(KIND=double), DIMENSION(jmax,kmax)    :: dtime_trans_all ! source and transport time step given by radial ray i_ray
REAL(KIND=double), DIMENSION(4)        :: d_send          ! send buffer
REAL(KIND=double), DIMENSION(50,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z) :: recv_buf1 ! receive buffer for gathering data from processors
REAL(KIND=double), DIMENSION(3,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z)  :: recv_buf2 ! receive buffer for gathering data from processors
REAL(KIND=double), DIMENSION(3,j_ray_dim,ik_ray_dim,n_proc_y*n_proc_z)   :: recv_buf3 ! receive buffer for gathering data from processors
REAL(KIND=double), DIMENSION(3,ij_ray_dim,k_ray_dim,n_proc_y*n_proc_z)   :: recv_buf4 ! receive buffer for gathering data from processors

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

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
!
!    \\\\\ GATHER TIME STEP RESTRICTIONS FOR GLOBAL SELECTION /////
!
!-----------------------------------------------------------------------

c_gath_send                   = 50 * ij_ray_dim * ik_ray_dim
c_gath_recv                   = 50 * ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( dt, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf1, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN
  DO m = 0,n_proc-1
  mj                          = MOD( m, n_proc_y )
    mk                        = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                     = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk                   = mj * ij_ray_dim + j
        dt_all(1:50,jsk,ksk)  = recv_buf1(1:50,j,k,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1
END IF ! myid == 0

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------


CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( jdt, c_gath_send, MPI_INTEGER, irecv_buf1, c_gath_recv, &
& MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
DO m = 0,n_proc-1
  mj                          = MOD( m, n_proc_y )
  mk                          = m / n_proc_y     
  DO k = 1,ik_ray_dim
    ksk                       = mk * ik_ray_dim + k
    DO j = 1,ij_ray_dim
      jsk                     = mj * ij_ray_dim + j
      jdt_all(1:50,jsk,ksk)   = irecv_buf1(1:50,j,k,m+1)
    END DO ! j = 1,ij_ray_dim
  END DO ! k = 1,ik_ray_dim
END DO ! m = 0,n_proc-1

c_gath_send                   = 3 * j_ray_dim * ik_ray_dim
c_gath_recv                   = 3 * j_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( dt_y, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf3, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN
  DO m = 0,n_proc-1
    mi                        = MOD( m, n_proc_y )
    mk                        = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                     = mk * ik_ray_dim + k
      DO i = 1,j_ray_dim
        isk                   = mi * j_ray_dim + i
        dt_y_all(1:3,isk,ksk) = recv_buf3(1:3,i,k,m+1)
      END DO ! j = 1,j_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1
END IF ! myid == 0

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------


CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( jdt_y, c_gath_send, MPI_INTEGER, irecv_buf3, c_gath_recv, &
& MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN
  DO m = 0,n_proc-1
    mi                        = MOD( m, n_proc_y )
    mk                        = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                     = mk * ik_ray_dim + k
      DO i = 1,j_ray_dim
        isk                   = mi * j_ray_dim + i
        jdt_y_all(1:3,isk,ksk) = irecv_buf3(1:3,i,k,m+1)
      END DO ! j = 1,j_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1
END IF ! myid == 0

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------

c_gath_send                   = 3 * ij_ray_dim * k_ray_dim
c_gath_recv                   = 3 * ij_ray_dim * k_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( dt_z, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf4, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN
  DO m = 0,n_proc-1
    mj                        = MOD( m, n_proc_y )
    mi                        = m / n_proc_y     
    DO j = 1,ij_ray_dim
      jsk                     = mj * ij_ray_dim + j
      DO i = 1,k_ray_dim
        isk                   = mi * k_ray_dim + i
        dt_z_all(1:3,jsk,isk) = recv_buf4(1:3,j,i,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1
END IF ! myid == 0

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( jdt_z, c_gath_send, MPI_INTEGER, irecv_buf4, c_gath_recv, &
& MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN
  DO m = 0,n_proc-1
    mj                        = MOD( m, n_proc_y )
    mi                        = m / n_proc_y     
    DO j = 1,ij_ray_dim
      jsk                     = mj * ij_ray_dim + j
      DO i = 1,k_ray_dim
        isk                   = mi * k_ray_dim + i
        jdt_z_all(1:3,jsk,isk) = irecv_buf4(1:3,j,i,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1
END IF ! myid == 0

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------

c_gath_send                   = ij_ray_dim * ik_ray_dim
c_gath_recv                   = ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( dtime_trans, c_gath_send, MPI_DOUBLE_PRECISION, dtime_trans_all, &
& c_gath_recv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!
!               ||||| BEGIN MYID = 0 EXECUTION |||||
!
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  The information in dt_y_all and jdt_y_all, i = 1 - 3, is transferred
!   to dt_all and jdt_all, i = 4 - 6. Thus, dt_all(i,j_ray,k_ray) gets
!   the minimum time for process i and the angular zone, j_ray, and
!   azimuthal zone, k_ray, that it occurs; the rest of the array gets
!   set above to dtmax.
!-----------------------------------------------------------------------

  IF ( ndim > 1 ) THEN

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
        dt_all (i,j_ray,k_ray)  = dt_test
        jdt_all(i,j_ray,k_ray)  = i_ray
      END IF ! j_ray /= 0
    END DO ! i = 4,6

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
        dt_all (i,j_ray,k_ray) = dt_test
        jdt_all(i,j_ray,k_ray) = i_ray
      END IF ! j_ray /= 0
    END DO ! i = 7,9

  END IF ! ndim > 1

!-----------------------------------------------------------------------
!  Select minimum time step if myid = 0
!-----------------------------------------------------------------------

  CALL time_step_select( ny, nz, dt_all, jdt_all, dt_process, j_radial_dt, &
&  j_angular_dt, j_azimuthal_dt, dtnph, dtnph_trans, dtime_trans_all,      &
&  ncynu_trns, intnu_trns, nutrans_trns, dt_y_state, dt_z_state )

!-----------------------------------------------------------------------
!  Laad send buffers
!-----------------------------------------------------------------------

  i_send(1)                   = intnu_trns
  i_send(2)                   = ncynu_trns

  d_send(1)                   = dtnph
  d_send(2)                   = dtnph_trans
  d_send(3)                   = dt_y_state
  d_send(4)                   = dt_z_state

!-----------------------------------------------------------------------
!
!                ||||| END MYID = 0 EXECUTION |||||
!
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!
!    \\\\\ BROADCAST TIME STEP RESTRICTIONS TO ALL PROCESSORS /////
!
!-----------------------------------------------------------------------

i_extent                      = 50

CALL MPI_BCAST( j_radial_dt   , i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST( j_angular_dt  , i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST( j_azimuthal_dt, i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST( dt_process    , i_extent, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST( d_send        , 4       , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST( i_send        , 2       , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST( nutrans_trns  , 1       , MPI_LOGICAL         , 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Unlaad send buffers
!-----------------------------------------------------------------------

intnu_trns                    = i_send(1)
ncynu_trns                    = i_send(2)

dtnph                         = d_send(1)
dtnph_trans                   = d_send(2)
dt_y_state                    = d_send(3)
dt_z_state                    = d_send(4)

!-----------------------------------------------------------------------
!  Put time step restrictions on tcntrl_module
!-----------------------------------------------------------------------

CALL time_step_brdcst( dt_process, j_radial_dt, j_angular_dt, j_azimuthal_dt, &
& dtnmh, dtnph, dtnph_trans )

!-----------------------------------------------------------------------
!  Reinitialize dt_y and jdt_y, and dt_z and jdt_z
!-----------------------------------------------------------------------

dt_y                          = dtmax
jdt_y                         = 0
dt_z                          = dtmax
jdt_z                         = 0

RETURN
END SUBROUTINE radyhd_to_time_step_select_MPI
