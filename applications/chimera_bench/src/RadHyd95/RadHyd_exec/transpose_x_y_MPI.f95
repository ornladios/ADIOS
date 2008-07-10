SUBROUTINE transpose_x_y( imin, imax, nx, jmin, jmax, ny, ij_ray_dim,     &
& ik_ray_dim, j_ray_dim, nz, nez, nnu, ls, le, nnc, n_proc, n_proc_y,     &
& n_proc_z, rho_c, t_c, ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e,        &
& xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_x, agr_c,     &
& e_nu_c, grav_y_c, rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, nse_y, psi0_y, &
& psi1_y, xn_y, a_nuc_rep_y, z_nuc_rep_y, be_nuc_rep_y, flat_x_y,         &
& agr_y, grav_y_cy, e_nu_y )
!-----------------------------------------------------------------------
!
!    File:         transpose_x_y
!    Module:       transpose_x_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive certain of the x-array variables in radial_ray_module,
!       transpose them so that all values of j reside on a given processor,
!       and send them back to be loaded into angular_ray_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  imin         : lower x-array index
!  imax         : upper x-array index
!  nx           : x-array extent
!  ij_ray_dim,ik_ray_dim    : number of radial rays assigned to a processor
!  jmin         : lower y-array index
!  jmax         : upper y-array index
!  ny           : y-array extent
!  j_ray_dim    : number of angular rays assigned to a processor
!  nz           : z-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  ls           : lower composition index
!  le           : upper composition index
!  nnc          : composition array extent
!  n_proc       : number of processors assigned to the run
!  n_proc_y     : number of processors assigned to the y-indices of rays
!  n_proc_y     : number of processors assigned to the z-indices of rays
!  rho_c        : densities (g cm^{-3})
!  t_c          : temperatures (K)
!  ye_c         : electron fractions
!  ei_c         : internal energies (ergs g^{-1})
!  u_c          : velocity x-components (cm s^{-1})
!  v_c          : velocity y-components (cm s^{-1})
!  w_c          : velocity z-components (cm s^{-1})
!  psi0_c       : zero angular moments of the neutrino occupation number
!  psi0_e       : first angular moments of the neutrino occupation number
!  xn_c         : initial mass fractions
!  be_nuc_rep_c : binding energy of mean heavy nucleus
!  z_nuc_rep_c  : mass number of mean heavy nucleus
!  a_nuc_rep_c  : charge number of mean heavy nucleus
!  flat_x       : variable indicating the presence of a radial a shock
!  agr_c        : zone-centered value of the lapse function
!  e_nu_c       : neutrino energy density [ergs cm^{-3}]
!  grav_y_c     : zone-centered y-component of gravitational acceleration
!  e_nu_c       : neutrino energy density [ergs cm^{-3}]
!
!    Output arguments:
!
!  rho_y        : densities (g cm^{-3})
!  t_y          : temperatures (K)
!  ye_y         : electron fractions
!  ei_y         : internal energies (ergs g^{-1})
!  u_y          : velocity x-components (cm s^{-1})
!  v_y          : velocity y-components (cm s^{-1})
!  w_y          : velocity z-components (cm s^{-1})
!  psi0_y       : zero angular moments of the neutrino occupation number
!  psi0_y       : first angular moments of the neutrino occupation number
!  xn_y         : initial mass fractions
!  be_nuc_rep_y : binding energy of mean heavy nucleus
!  a_nuc_rep_y  : mass number of mean heavy nucleus
!  z_nuc_rep_y  : charge number of mean heavy nucleus
!  flat_x_y     : transposed variable indicating the presence of radial a shock
!  agr_y        : zone-centered value of the lapse function
!  grav_y_cy    : zone-centered y-component of gravitational acceleration
!  e_nu_y       : neutrino energy density [ergs cm^{-3}]
!
!      
!    Include files:
!  kind_module
!  nu_energy_grid_module, prb_cntl_module, shock_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : v_trans_0
USE shock_module, ONLY: j_shk_radial_p, j_shk_radial_all_p
USE edit_module, ONLY : nlog

USE parallel_module, ONLY : myid, ierr, MPI_COMM_ROW, myid_y
USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin               ! lower x-array index
INTEGER, INTENT(in)              :: imax               ! upper x-array index
INTEGER, INTENT(in)              :: nx                 ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim         ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim         ! number of z-zones on a processor before swapping with z

INTEGER, INTENT(in)              :: jmin               ! lower y-array index
INTEGER, INTENT(in)              :: jmax               ! upper y-array index
INTEGER, INTENT(in)              :: ny                 ! y-array extent
INTEGER, INTENT(in)              :: j_ray_dim          ! number of angular rays assigned to a processor

INTEGER, INTENT(in)              :: nz                 ! z-array extent
INTEGER, INTENT(in)              :: nez                ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu                ! neutrino flavor array extent

INTEGER, INTENT(in)              :: ls                 ! lower composition index
INTEGER, INTENT(in)              :: le                 ! upper composition index
INTEGER, INTENT(in)              :: nnc                ! composition array extent

INTEGER, INTENT(in)              :: n_proc             ! number of processors assigned to the run
INTEGER, INTENT(in)              :: n_proc_y           ! number of processors assigned to the y-indices of rays
INTEGER, INTENT(in)              :: n_proc_z           ! number of processors assigned to the z-indices of rays

INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)                   :: nse_c        ! nse flag

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rho_c        ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: t_c          ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: ye_c         ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: ei_c         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: u_c          ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: v_c          ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: w_c          ! z-velocity of zone (cm)

REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_c       ! zeroth angular moment of the NDS
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1_e       ! first angular moment of the NDS

REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim)     :: xn_c         ! mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: a_nuc_rep_c  ! nuclear mass number
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: z_nuc_rep_c  ! nuclear charge number
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: be_nuc_rep_c ! nuclear binding energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: flat_x       ! variable indicating the presence of radial a shock
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: agr_c        ! zone-centered value of the lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: e_nu_c       ! neutrino energy density [ergs cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: grav_y_c     ! zone-centered y-component of gravitational acceleration

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)                   :: nse_y        ! nse flag

REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: rho_y        ! density (cm^{-3})
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: t_y          ! temperature (MeV)
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: ye_y         ! temperature (MeV)
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: ei_y         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: u_y          ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: v_y          ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: w_y          ! z-velocity of zone (cm)

REAL(KIND=double), INTENT(out), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim) :: psi0_y       ! zeroth angular moment of the NDS
REAL(KIND=double), INTENT(out), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim) :: psi1_y       ! first angular moment of the NDS

REAL(KIND=double), INTENT(out), DIMENSION(ny,nnc,j_ray_dim,ik_ray_dim)     :: xn_y         ! mass fractions
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: a_nuc_rep_y  ! nuclear mass number
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: z_nuc_rep_y  ! nuclear charge number
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: be_nuc_rep_y ! nuclear binding energy
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: flat_x_y     ! variable indicating the presence of radial a shock
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: agr_y        ! zone-centered value of the lapse function
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: grav_y_cy    ! zone-centered y-component of gravitational acceleration
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: e_nu_y       ! neutrino energy density [ergs cm^{-3}]


!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                        :: i             ! x-array index
INTEGER                        :: j             ! y-array index
INTEGER                        :: k             ! z-array index
INTEGER                        :: ie            ! neutrino energy index
INTEGER                        :: l             ! composition index
INTEGER                        :: n             ! neutrino flavor index
INTEGER                        :: m             ! processor index
INTEGER                        :: mj            ! y-block index
INTEGER                        :: mk            ! z-block index
INTEGER                        :: jsk           ! y-array index of gathered array
INTEGER                        :: ksk           ! z-array index of gathered array
INTEGER                        :: npy           ! y-processor index
INTEGER                        :: ij_ray        ! j-index of a radial ray
INTEGER                        :: ik_ray        ! ie-index of a radial ray
INTEGER                        :: count_ata     ! all to all buffer count
INTEGER                        :: c_gath_send   ! gather send buffer count
INTEGER                        :: c_gath_recv   ! gather recv buffer count
INTEGER                        :: i_extent      ! broadcast array extent

INTEGER, DIMENSION(ik_ray_dim,1,ij_ray_dim,imax)                                 :: i_send ! All to all send buffer
INTEGER, DIMENSION(ik_ray_dim,1,ij_ray_dim,j_ray_dim,n_proc_y)                   :: i_recv ! All to all receive buffer
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z)                      :: i_gath ! receive buffer for gathering data from processors

REAL(KIND=double), DIMENSION(ik_ray_dim,14+nnc,ij_ray_dim,imax)                  :: send   ! All to all send buffer
REAL(KIND=double), DIMENSION(ik_ray_dim,14+nnc,ij_ray_dim,j_ray_dim,n_proc_y)    :: recv   ! All to all receive buffer
REAL(KIND=double), DIMENSION(ik_ray_dim,2,nez,nnu,ij_ray_dim,imax)               :: p_send ! All to all send buffer
REAL(KIND=double), DIMENSION(ik_ray_dim,2,nez,nnu,ij_ray_dim,j_ray_dim,n_proc_y) :: p_recv ! All to all receive buffer

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                    \\\\\ LOAD Y-ARRAYS /////
!
!-----------------------------------------------------------------------

c_gath_send                          = ij_ray_dim * ik_ray_dim
c_gath_recv                          = ij_ray_dim * ik_ray_dim

IF ( v_trans_0 == 'ye' ) THEN
  CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
  CALL MPI_GATHER( j_shk_radial_p, c_gath_send, MPI_INTEGER, i_gath, &
&  c_gath_recv, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

  IF ( myid == 0 ) THEN
    DO m = 0,n_proc-1
      mj                             = MOD( m, n_proc_y )
      mk                             = m / n_proc_y     
      DO k = 1,ik_ray_dim
        ksk                          = mk * ik_ray_dim + k
        DO j = 1,ij_ray_dim
          jsk                        = mj * ij_ray_dim + j
          j_shk_radial_all_p(jsk,ksk) = i_gath(j,k,m+1)
        END DO ! j = 1,ij_ray_dim
      END DO ! k = 1,ik_ray_dim
    END DO ! m = 0,n_proc-1
  END IF ! myid == 0

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------

  CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
   i_extent                          = ny * nz
  CALL MPI_BCAST( j_shk_radial_all_p, i_extent, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
END IF

!-----------------------------------------------------------------------
!
!            \\\\\ ALL_TO_ALL TRANSFER IF N_PROC > 1 /////
!
!-----------------------------------------------------------------------

IF ( n_proc > 1 ) THEN

  count_ata                          = ( 14+nnc ) * ik_ray_dim * j_ray_dim * ij_ray_dim

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim
      DO i = imin,imax
        send(ik_ray, 1,ij_ray,i)     = rho_c       (i,ij_ray,ik_ray)
        send(ik_ray, 2,ij_ray,i)     = t_c         (i,ij_ray,ik_ray)
        send(ik_ray, 3,ij_ray,i)     = ye_c        (i,ij_ray,ik_ray)
        send(ik_ray, 4,ij_ray,i)     = ei_c        (i,ij_ray,ik_ray)
        send(ik_ray, 5,ij_ray,i)     = u_c         (i,ij_ray,ik_ray)
        send(ik_ray, 6,ij_ray,i)     = v_c         (i,ij_ray,ik_ray)
        send(ik_ray, 7,ij_ray,i)     = w_c         (i,ij_ray,ik_ray)
        send(ik_ray, 8,ij_ray,i)     = a_nuc_rep_c (i,ij_ray,ik_ray)
        send(ik_ray, 9,ij_ray,i)     = z_nuc_rep_c (i,ij_ray,ik_ray)
        send(ik_ray,10,ij_ray,i)     = be_nuc_rep_c(i,ij_ray,ik_ray)
        send(ik_ray,11,ij_ray,i)     = flat_x      (i,ij_ray,ik_ray)
        send(ik_ray,12,ij_ray,i)     = agr_c       (i,ij_ray,ik_ray)
        send(ik_ray,13,ij_ray,i)     = grav_y_c    (i,ij_ray,ik_ray)
        send(ik_ray,14,ij_ray,i)     = e_nu_c      (i,ij_ray,ik_ray)
        DO l = ls,le
          send(ik_ray,14+l,ij_ray,i) = xn_c(i,l,ij_ray,ik_ray)
        END DO ! l = ls,le
      END DO ! i = imin,imax
    END DO ! ij_ray = 1,ij_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

  CALL MPI_ALLTOALL( send, count_ata, MPI_DOUBLE_PRECISION, recv, count_ata, &
&   MPI_DOUBLE_PRECISION, MPI_COMM_ROW, ierr )

  DO ik_ray = 1,ik_ray_dim
    DO i = 1,j_ray_dim
      DO ij_ray = 1,ij_ray_dim
        DO npy = 1,n_proc_y
          j                          = ij_ray + ij_ray_dim * ( npy - 1 )
          rho_y       (j,i,ik_ray)   = recv(ik_ray, 1,ij_ray,i,npy)
          t_y         (j,i,ik_ray)   = recv(ik_ray, 2,ij_ray,i,npy)
          ye_y        (j,i,ik_ray)   = recv(ik_ray, 3,ij_ray,i,npy)
          ei_y        (j,i,ik_ray)   = recv(ik_ray, 4,ij_ray,i,npy)
          u_y         (j,i,ik_ray)   = recv(ik_ray, 5,ij_ray,i,npy)
          v_y         (j,i,ik_ray)   = recv(ik_ray, 6,ij_ray,i,npy)
          w_y         (j,i,ik_ray)   = recv(ik_ray, 7,ij_ray,i,npy)
          a_nuc_rep_y (j,i,ik_ray)   = recv(ik_ray, 8,ij_ray,i,npy)
          z_nuc_rep_y (j,i,ik_ray)   = recv(ik_ray, 9,ij_ray,i,npy)
          be_nuc_rep_y(j,i,ik_ray)   = recv(ik_ray,10,ij_ray,i,npy)
          flat_x_y    (j,i,ik_ray)   = recv(ik_ray,11,ij_ray,i,npy)
          agr_y       (j,i,ik_ray)   = recv(ik_ray,12,ij_ray,i,npy)
          grav_y_cy   (j,i,ik_ray)   = recv(ik_ray,13,ij_ray,i,npy)
          e_nu_y      (j,i,ik_ray)   = recv(ik_ray,14,ij_ray,i,npy)
          DO l = ls,le
            xn_y      (j,l,i,ik_ray) = recv(ik_ray,14+l,ij_ray,i,npy)
          END DO ! l = ls,le
        END DO ! npy = 1,n_proc_y
      END DO ! ij_ray = 1,ij_ray_dim
    END DO ! i = 1,j_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

  count_ata                          = ik_ray_dim * j_ray_dim * ij_ray_dim

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim
      DO i = imin,imax
        i_send(ik_ray,1,ij_ray,i)    = nse_c(i,ij_ray,ik_ray)
       END DO ! i = imin,imax
    END DO ! ij_ray = 1,ij_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

  CALL MPI_ALLTOALL( i_send, count_ata, MPI_INTEGER, i_recv, count_ata, &
&   MPI_INTEGER, MPI_COMM_ROW, ierr )

  DO ik_ray = 1,ik_ray_dim
    DO i = 1,j_ray_dim
      DO ij_ray = 1,ij_ray_dim
        DO npy = 1,n_proc_y
          j                          = ij_ray + ij_ray_dim * ( npy - 1 )
          nse_y(j,i,ik_ray)          = i_recv(ik_ray, 1,ij_ray,i,npy)
        END DO ! npy = 1,n_proc_y
      END DO ! ij_ray = 1,ij_ray_dim
    END DO ! i = 1,j_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

  IF ( nnugpmx /= 0 ) THEN

    count_ata                        = 2 * nez * nnu * ik_ray_dim * j_ray_dim * ij_ray_dim

    DO ik_ray = 1,ik_ray_dim
      DO ij_ray = 1,ij_ray_dim
        DO i = imin,imax
          DO n = 1,nnu
            IF ( nnugp(n) /= 0 ) THEN
              DO ie = 1,nnugp(n)
                p_send(ik_ray,1,ie,n,ij_ray,i) = psi0_c(i,ie,n,ij_ray,ik_ray)
                p_send(ik_ray,2,ie,n,ij_ray,i) = psi1_e(i,ie,n,ij_ray,ik_ray)
              END DO ! ie = 1,nnugp(n)
            END IF ! nnugp(n) /= 0
          END DO ! n = 1,nnu
        END DO ! i = imin,imax
      END DO ! ij_ray = 1,ij_ray_dim
    END DO ! ik_ray = 1,ik_ray_dim

    CALL MPI_ALLTOALL( p_send, count_ata, MPI_DOUBLE_PRECISION, p_recv, count_ata, &
&     MPI_DOUBLE_PRECISION, MPI_COMM_ROW, ierr )

    DO ik_ray = 1,ik_ray_dim
      DO i = 1,j_ray_dim
        DO ij_ray = 1,ij_ray_dim
          DO npy = 1,n_proc_y
            j                        = ij_ray + ij_ray_dim * ( npy - 1 )
            DO n = 1,nnu
              IF ( nnugp(n) /= 0 ) THEN
                DO ie = 1,nnugp(n)
                  psi0_y(j,ie,n,i,ik_ray) = p_recv(ik_ray,1,ie,n,ij_ray,i,npy)
                  psi1_y(j,ie,n,i,ik_ray) = p_recv(ik_ray,2,ie,n,ij_ray,i,npy)
                END DO ! ie = 1,nnugp(n)
              END IF ! nnugp(n) /= 0
            END DO ! n = 1,nnu
          END DO ! npy = 1,n_proc_y
        END DO ! ij_ray = 1,ij_ray_dim
      END DO ! i = 1,j_ray_dim
    END DO ! ik_ray = 1,ik_ray_dim

  END IF ! nnugpmx /= 0

!-----------------------------------------------------------------------
!
!            \\\\\ TRANSPOSE VARIABLES IF N_PROC = 1 /////
!
!-----------------------------------------------------------------------

ELSE ! n_procs = 1

  DO i = imin,imax
    DO j = jmin,jmax
      rho_y       (j,i,:)            = rho_c       (i,j,:)
      t_y         (j,i,:)            = t_c         (i,j,:)
      ye_y        (j,i,:)            = ye_c        (i,j,:)
      ei_y        (j,i,:)            = ei_c        (i,j,:)
      u_y         (j,i,:)            = u_c         (i,j,:)
      v_y         (j,i,:)            = v_c         (i,j,:)
      w_y         (j,i,:)            = w_c         (i,j,:)
      a_nuc_rep_y (j,i,:)            = a_nuc_rep_c (i,j,:)
      z_nuc_rep_y (j,i,:)            = z_nuc_rep_c (i,j,:)
      be_nuc_rep_y(j,i,:)            = be_nuc_rep_c(i,j,:)
      nse_y       (j,i,:)            = nse_c       (i,j,:)
      flat_x_y    (j,i,:)            = flat_x      (i,j,:)
      agr_y       (j,i,:)            = agr_c       (i,j,:)
      grav_y_cy   (j,i,:)            = grav_y_c    (i,j,:)
      e_nu_y      (j,i,:)            = e_nu_c      (i,j,:)
    END DO
  END DO

  DO l = ls,le
    DO i = imin,imax
      DO j = jmin,jmax
        xn_y(j,l,i,:)                = xn_c(i,l,j,:)
      END DO
    END DO
  END DO

  DO n = 1,nnu
    IF ( nnugp(n) == 0 ) CYCLE
    DO ie = 1,nnugp(n)
      DO i = imin,imax
        DO j = jmin,jmax
          psi0_y(j,ie,n,i,:)         = psi0_c(i,ie,n,j,:) 
          psi1_y(j,ie,n,i,:)         = psi1_e(i,ie,n,j,:) 
        END DO
      END DO
    END DO
  END DO

END IF ! ! n_procs > 1

RETURN
END SUBROUTINE transpose_x_y
