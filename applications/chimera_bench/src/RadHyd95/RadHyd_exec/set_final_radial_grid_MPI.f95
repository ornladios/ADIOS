 SUBROUTINE set_final_radial_grid_MPI( ij_ray_dim, ik_ray_dim, nx )
!-----------------------------------------------------------------------
!
!    File:         set_final_radial_grid_MPI
!    Module:       set_final_radial_grid_MPI
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/05/05
!
!    Purpose:
!      To set the radial grid at the end of the hydro step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  nx         : x-array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  cycle_module, edit_module, parallel_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc, n_proc_y, n_proc_z
USE numerical_module, ONLY : zero, half, frpi, epsilon

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid, ierr
USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax,  m_grid,   &
& u_e, x_ei, x_ei, dx_ci, x_ci, x_ef, dx_cf, x_cf, x_el, dx_cl, x_cl, dtnph, &
& d_omega, lagr

USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER                               :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER                               :: ik_ray_dim       ! number of z-zones on a processor before swapping
INTEGER                               :: nx               ! x-array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                               :: first_y = .true.

INTEGER                               :: i                ! x-array index
INTEGER                               :: j                ! y-array index
INTEGER                               :: c_gath_send      ! gather send buffer count
INTEGER                               :: c_gath_recv      ! gather recv buffer count
INTEGER                               :: k                ! z-array index
INTEGER                               :: m                ! processor index
INTEGER                               :: mj               ! y-block index
INTEGER                               :: mk               ! z-block index
INTEGER                               :: jsk              ! y-array index of gathered array
INTEGER                               :: ksk              ! z-array index of gathered array

REAL(KIND=double)                          :: y_shift     ! y-grid shift parameter

REAL(KIND=double), DIMENSION(nx+1)         :: u_grid      ! x-velocity of the grid
REAL(KIND=double), DIMENSION(nx,jmax,kmax) :: u_e_all     ! x-velocity of the grid
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z) :: recv_buf ! receive buffer for gathering data from processors

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ LAGRANGIAN GRID /////
!
!-----------------------------------------------------------------------

IF ( lagr == 'ye' ) THEN

  x_ef (imin:imax+1)      = x_el (imin:imax+1,1,1)
  dx_cf(imin:imax)        = dx_cl(imin:imax  ,1,1)
  x_cf (imin:imax)        = x_cl (imin:imax  ,1,1)

  RETURN

END IF ! lagr == 'ye'

!-----------------------------------------------------------------------
!
!                 \\\\\ EULERIAN, MOVING GRID OFF /////
!
!-----------------------------------------------------------------------

IF ( lagr == 'no'  .and.  m_grid == 'no' ) THEN

!-----------------------------------------------------------------------
!  Set final grid to the initial grid
!-----------------------------------------------------------------------

  x_ef                    = x_ei
  dx_cf                   = dx_ci
  x_cf                    = x_ci

  RETURN

END IF ! lagr == 'no'  .and.  m_grid == 'no'

!-----------------------------------------------------------------------
!
!                      \\\\\ MOVING GRID ON /////
!
!-----------------------------------------------------------------------

IF ( lagr == 'no'  .and.  m_grid == 'ye' ) THEN

!-----------------------------------------------------------------------
!  Gather edge velocities from all processors
!-----------------------------------------------------------------------

  c_gath_send             = nx * ij_ray_dim * ik_ray_dim
  c_gath_recv             = nx * ij_ray_dim * ik_ray_dim

  CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
  CALL MPI_GATHER( u_e, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf,    &
&  c_gath_recv, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

  IF ( myid == 0 ) THEN

    DO m = 0,n_proc-1
      mj                  = MOD( m, n_proc_y )
      mk                  = m / n_proc_y     
      DO k = 1,ik_ray_dim
        ksk               = mk * ik_ray_dim + k
        DO j = 1,ij_ray_dim
          jsk             = mj * ij_ray_dim + j
          u_e_all(1:nx,jsk,ksk)                                         &
&                         = recv_buf(1:nx,j,k,m+1)
        END DO ! j = 1,ij_ray_dim
      END DO ! k = 1,ik_ray_dim
    END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!  Compute average edge velocities if myid = 0
!-----------------------------------------------------------------------

    DO i = imin,imax+1
      u_grid(i)           = SUM( u_e_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i

  END IF ! myid == 0

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Broadcast average edge velocities to all processors
!-----------------------------------------------------------------------

  CALL MPI_BCAST( u_grid, nx+1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  x_ef (imin:imax+1)      = x_ei(imin:imax+1) + u_grid(imin:imax+1) * dtnph
  dx_cf(imin:imax)        = x_ef(imin+1:imax+1) - x_ef(imin:imax)
  x_cf (imin:imax)        = half * ( x_ef(imin+1:imax+1) + x_ef(imin:imax) )

  RETURN

END IF ! lagr == 'no'  .and.  m_grid == 'ye'

RETURN
END SUBROUTINE set_final_radial_grid_MPI
