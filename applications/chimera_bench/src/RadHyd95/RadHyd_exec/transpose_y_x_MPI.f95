SUBROUTINE transpose_y_x(imin, imax, nx, ij_ray_dim, ik_ray_dim, jmin, &
& jmax, ny, j_ray_dim, nez, nnu, ls, le, nnc, n_proc, n_proc_y, rho_c, &
& t_c, ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c,   &
& z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_y_x, pq_y, rho_y, t_y, ye_y,  &
& ei_y, u_y, v_y, w_y, psi0_y, psi1_y, xn_y, a_nuc_rep_y, z_nuc_rep_y, &
& be_nuc_rep_y, nse_y, flat_y, pqy_x )
!-----------------------------------------------------------------------
!
!    File:         transpose_y_x
!    Module:       transpose_y_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive certain of the y-array variables in angular_ray_module,
!       transpose them, and send them back to be loaded into raddhyd_ray_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  imin         : lower x-array index
!  imax         : upper x-array index
!  nx           : x-array extent
!  ij_ray_dim   : number of x (radial) zones on a processor before swapping with y
!  ik_ray_dim   : number of z (azimuthal) zones on a processor before swapping with z
!  jmin         : lower y-array index
!  jmax         : upper y-array index
!  ny           : y-array extent
!  j_ray_dim,   : number of radial zones on a processor after swapping with y
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  ls           : lower composition index
!  le           : upper composition index
!  nnc          : composition array extent
!  n_proc       : number of processors assigned to the run
!  n_proc_y     : number of processors assigned to the y-indices of rays
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
!  flat_y       : variable indicating the presence of an amgular shock
!  pq_y         : pseudoviscous pressure in the y-direction
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
!  flat_y_x     : variable indicating the presence of an amgular shock
!  pqy_x        : pseudoviscous pressure in the y-direction
!
!      
!    Include files:
!  kind_module
!  nu_energy_grid_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

USE parallel_module, ONLY : myid, ierr, MPI_COMM_ROW
USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin               ! inner physical x-zone
INTEGER, INTENT(in)              :: imax               ! outer physical x-zone
INTEGER, INTENT(in)              :: nx                 ! logical dimension of x-zone
INTEGER, INTENT(in)              :: ij_ray_dim         ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim         ! number of radial zones on a processor before swapping with z

INTEGER, INTENT(in)              :: jmin               ! inner physical y-zone
INTEGER, INTENT(in)              :: jmax               ! outer physical y-zone
INTEGER, INTENT(in)              :: ny                 ! logical dimension of y-zone
INTEGER, INTENT(in)              :: j_ray_dim          ! number of radial zones on a processor after swapping with y

INTEGER, INTENT(in)              :: nez                ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu                ! neutrino flavor array extent

INTEGER, INTENT(in)              :: ls                 ! lower composition index
INTEGER, INTENT(in)              :: le                 ! upper composition index
INTEGER, INTENT(in)              :: nnc                ! composition array extent

INTEGER, INTENT(in)              :: n_proc             ! number of processors assigned to the run
INTEGER, INTENT(in)              :: n_proc_y           ! number of processors assigned to the y-indices of rays

INTEGER, INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)                    :: nse_y        ! nse flag

REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: rho_y        ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: t_y          ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: ye_y         ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: ei_y         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: u_y          ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: v_y          ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: w_y          ! z-velocity of zone (cm)

REAL(KIND=double), INTENT(in), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim)  :: psi0_y       ! zeroth angular moment of the NDS
REAL(KIND=double), INTENT(in), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim)  :: psi1_y       ! first angular moment of the NDS

REAL(KIND=double), INTENT(in), DIMENSION(ny,nnc,j_ray_dim,ik_ray_dim)      :: xn_y         ! mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: a_nuc_rep_y  ! nuclear mass number
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: z_nuc_rep_y  ! nuclear charge number
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: be_nuc_rep_y ! nuclear binding energy
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: flat_y       ! variable indicating the presence of an angular shock
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: pq_y         ! pseudoviscous pressure in the y-direction

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)                   :: nse_c        ! nse flag

REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rho_c        ! density (cm^{-3})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: t_c          ! temperature (MeV)
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: ye_c         ! temperature (MeV)
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: ei_c         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: u_c          ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: v_c          ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: w_c          ! z-velocity of zone (cm)

REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_c       ! zeroth angular moment of the NDS
REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1_e       ! first angular moment of the NDS

REAL(KIND=double), INTENT(out), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim)     :: xn_c         ! mass fractions
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: a_nuc_rep_c  ! nuclear mass number
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: z_nuc_rep_c  ! nuclear charge number
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: be_nuc_rep_c ! nuclear binding energy
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: flat_y_x     ! variable indicating the presence of an angular shock
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: pqy_x        ! pseudoviscous pressure in the y-direction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                        :: i             ! x-array index
INTEGER                        :: iblk          ! x-array index after swapping with y
INTEGER                        :: npy           ! y-processor index
INTEGER                        :: j             ! y-array index
INTEGER                        :: k             ! z-array index
INTEGER                        :: l             ! composition index
INTEGER                        :: ie            ! neutrino energy index
INTEGER                        :: n             ! neutrino flavor index
INTEGER                        :: count_ata     ! all to all buffer count

INTEGER, DIMENSION(ik_ray_dim,1,j_ray_dim,jmax)                                  :: i_send ! All to all send buffer
INTEGER, DIMENSION(ik_ray_dim,1,j_ray_dim,ij_ray_dim,n_proc_y)                   :: i_recv ! All to all receive buffer

REAL(KIND=double), DIMENSION(ik_ray_dim,12+nnc,j_ray_dim,jmax)                   :: send   ! All to all send buffer
REAL(KIND=double), DIMENSION(ik_ray_dim,12+nnc,j_ray_dim,ij_ray_dim,n_proc_y)    :: recv   ! All to all receive buffer
REAL(KIND=double), DIMENSION(ik_ray_dim,2,nez,nnu,j_ray_dim,jmax)                :: p_send ! All to all send buffer
REAL(KIND=double), DIMENSION(ik_ray_dim,2,nez,nnu,j_ray_dim,ij_ray_dim,n_proc_y) :: p_recv ! All to all receive buffer

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ ALL_TO_ALL TRANSFER IF N_PROC > 1 /////
!
!-----------------------------------------------------------------------

IF ( n_proc > 1 ) THEN

  count_ata                          = ( 12+nnc ) * ik_ray_dim * j_ray_dim * ij_ray_dim

  DO k = 1,ik_ray_dim
    DO i = 1,j_ray_dim
      DO j = jmin,jmax
        send(k, 1,i,j)               = rho_y       (j,i,k)
        send(k, 2,i,j)               = t_y         (j,i,k)
        send(k, 3,i,j)               = ye_y        (j,i,k)
        send(k, 4,i,j)               = ei_y        (j,i,k)
        send(k, 5,i,j)               = u_y         (j,i,k)
        send(k, 6,i,j)               = v_y         (j,i,k)
        send(k, 7,i,j)               = w_y         (j,i,k)
        send(k, 8,i,j)               = a_nuc_rep_y (j,i,k)
        send(k, 9,i,j)               = z_nuc_rep_y (j,i,k)
        send(k,10,i,j)               = be_nuc_rep_y(j,i,k)
        send(k,11,i,j)               = flat_y      (j,i,k)
        send(k,12,i,j)               = pq_y        (j,i,k)
        send(k,12+ls:12+le,i,j)      = xn_y        (j,ls:le,i,k)
      END DO ! j = jmin,jmax
    END DO ! i = 1,j_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

  CALL MPI_ALLTOALL( send, count_ata, MPI_DOUBLE_PRECISION, recv, count_ata, &
&   MPI_DOUBLE_PRECISION, MPI_COMM_ROW, ierr )

  DO k = 1,ik_ray_dim
    DO j = 1,ij_ray_dim
      DO iblk = 1,j_ray_dim
        DO npy = 1,n_proc_y
          i                          = iblk + j_ray_dim * ( npy - 1 )
          rho_c       (i,j,k)        = recv(k, 1,iblk,j,npy)
          t_c         (i,j,k)        = recv(k, 2,iblk,j,npy)
          ye_c        (i,j,k)        = recv(k, 3,iblk,j,npy)
          ei_c        (i,j,k)        = recv(k, 4,iblk,j,npy)
          u_c         (i,j,k)        = recv(k, 5,iblk,j,npy)
          v_c         (i,j,k)        = recv(k, 6,iblk,j,npy)
          w_c         (i,j,k)        = recv(k, 7,iblk,j,npy)
          a_nuc_rep_c (i,j,k)        = recv(k, 8,iblk,j,npy)
          z_nuc_rep_c (i,j,k)        = recv(k, 9,iblk,j,npy)
          be_nuc_rep_c(i,j,k)        = recv(k,10,iblk,j,npy)
          flat_y_x    (i,j,k)        = recv(k,11,iblk,j,npy)
          pqy_x       (i,j,k)        = recv(k,12,iblk,j,npy)
          xn_c        (i,ls:le,j,k)  = recv(k,12+ls:12+le,iblk,j,npy)
        END DO ! npy = 1,n_proc_y
      END DO ! iblk = 1,j_ray_dim
    END DO ! j = 1,ij_ray_dim
  END DO ! k = 1,ik_ray_dim

  count_ata                          = ik_ray_dim * j_ray_dim * ij_ray_dim

  DO k = 1,ik_ray_dim
    DO i = 1,j_ray_dim
      DO j = jmin,jmax
        i_send(k,1,i,j)              = nse_y(j,i,k)
      END DO ! j = jmin,jmax
    END DO ! i = 1,j_ray_dim
  END DO ! k = 1,ik_ray_dim

  CALL MPI_ALLTOALL( i_send, count_ata, MPI_INTEGER, i_recv, count_ata, &
&   MPI_INTEGER, MPI_COMM_ROW, ierr )

  DO k = 1,ik_ray_dim
    DO j = 1,ij_ray_dim
      DO iblk = 1,j_ray_dim
        DO npy = 1,n_proc_y
          i                          = iblk + j_ray_dim * ( npy - 1 )
          nse_c(i,j,k)               = i_recv(k,1,iblk,j,npy)
        END DO ! npy = 1,n_proc_y
      END DO ! iblk = 1,j_ray_dim
    END DO ! j = 1,ij_ray_dim
  END DO ! k = 1,ik_ray_dim

  IF ( nnugpmx /= 0 ) THEN

    count_ata                        = 2 * nez * nnu * ik_ray_dim * j_ray_dim * ij_ray_dim

    DO k = 1,ik_ray_dim
      DO i = 1,j_ray_dim
        DO j = jmin,jmax
          DO n = 1,nnu
            IF ( nnugp(n) /= 0 ) THEN
              DO ie = 1,nnugp(n)
                p_send(k,1,ie,n,i,j) = psi0_y(j,ie,n,i,k)
                p_send(k,2,ie,n,i,j) = psi1_y(j,ie,n,i,k)
              END DO ! ie = 1,nnugp(n)
            END IF ! nnugp(n) /= 0
          END DO ! n = 1,nnu
        END DO ! j = jmin,jmax
      END DO ! i = 1,j_ray_dim
    END DO ! k = 1,ik_ray_dim

    CALL MPI_ALLTOALL( p_send, count_ata, MPI_DOUBLE_PRECISION, p_recv, count_ata, &
&     MPI_DOUBLE_PRECISION, MPI_COMM_ROW, ierr )


    DO k = 1,ik_ray_dim
      DO j = 1,ij_ray_dim
        DO iblk = 1,j_ray_dim
          DO npy = 1,n_proc_y
            i                        = iblk + j_ray_dim * ( npy - 1 )
            DO n = 1,nnu
              IF ( nnugp(n) /= 0 ) THEN
                DO ie = 1,nnugp(n)
                  psi0_c(i,ie,n,j,k) = p_recv(k,1,ie,n,iblk,j,npy)
                  psi1_e(i,ie,n,j,k) = p_recv(k,2,ie,n,iblk,j,npy)
                END DO ! ie = 1,nnugp(n)
              END IF ! nnugp(n) /= 0
            END DO ! n = 1,nnu
          END DO ! npy = 1,n_proc_y
        END DO ! iblk = 1,j_ray_dim
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim

  END IF ! nnugpmx /= 0

!-----------------------------------------------------------------------
!
!            \\\\\ TRANSPOSE VARIABLES IF N_PROC = 1 /////
!
!-----------------------------------------------------------------------

ELSE ! n_procs = 1

  DO i = imin,imax
    DO j = jmin,jmax
      rho_c       (i,j,:)            = rho_y       (j,i,:)
      t_c         (i,j,:)            = t_y         (j,i,:)
      ye_c        (i,j,:)            = ye_y        (j,i,:)
      ei_c        (i,j,:)            = ei_y        (j,i,:)
      u_c         (i,j,:)            = u_y         (j,i,:)
      v_c         (i,j,:)            = v_y         (j,i,:)
      w_c         (i,j,:)            = w_y         (j,i,:)
      a_nuc_rep_c (i,j,:)            = a_nuc_rep_y (j,i,:)
      z_nuc_rep_c (i,j,:)            = z_nuc_rep_y (j,i,:)
      be_nuc_rep_c(i,j,:)            = be_nuc_rep_y(j,i,:)
      nse_c       (i,j,:)            = nse_y       (j,i,:)
      flat_y_x    (i,j,:)            = flat_y      (j,i,:)
      pqy_x       (i,j,:)            = pq_y        (j,i,:)
    END DO
  END DO

  DO l = ls,le
    DO i = imin,imax
      DO j = jmin,jmax
        xn_c(i,l,j,:)                = xn_y(j,l,i,:)
      END DO
    END DO
  END DO

  DO n = 1,nnu
    IF ( nnugp(n) == 0 ) CYCLE
    DO ie = 1,nnugp(n)
      DO i = imin,imax
        DO j = jmin,jmax
          psi0_c(i,ie,n,j,:)         = psi0_y(j,ie,n,i,:)
          psi1_e(i,ie,n,j,:)         = psi1_y(j,ie,n,i,:)
        END DO
      END DO
    END DO
  END DO

END IF ! ! n_procs > 1

RETURN
END SUBROUTINE transpose_y_x
