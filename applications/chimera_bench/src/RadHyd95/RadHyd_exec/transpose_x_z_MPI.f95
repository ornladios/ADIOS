SUBROUTINE transpose_x_z( imin, imax, nx, kmin, kmax, nz, ij_ray_dim,     &
& ik_ray_dim, k_ray_dim, ny, nez, nnu, ls, le, nnc, n_proc, n_proc_y,     &
& n_proc_z, rho_c, t_c, ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e,        &
& xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_x, agr_c,     &
& e_nu_c, grav_z_c, rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, nse_z, psi0_z, &
& psi1_z, xn_z, a_nuc_rep_z, z_nuc_rep_z, be_nuc_rep_z, flat_x_z,         &
& agr_z, grav_z_cz, e_nu_z )
!-----------------------------------------------------------------------
!
!    File:         transpose_x_z
!    Module:       transpose_x_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive certain of the x-array variables in radial_ray_module,
!       transpose them so that all values of k reside on a given processor,
!       and send them back to be loaded into azimuthal_ray_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  imin         : lower x-array index
!  imax         : upper x-array index
!  nx           : x-array extent
!  ij_ray_dim   : number of y-zones on a processor before swapping with y
!  ik_ray_dim   : number of z-zones on a processor before swapping with z
!  kmin         : lower z-array index
!  kmax         : upper z-array index
!  ny           : y-array extent
!  k_ray_dim    : number of radial zones on a processor after swapping with y
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  ls           : lower composition index
!  le           : upper composition index
!  nnc          : composition array extent
!  n_proc       : number of processors assigned to the run
!  n_proc_y     : number of processors assigned to the y-indices of rays
!  n_proc_z     : number of processors assigned to the z-indices of rays
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
!  grav_z_c     : zone-centered z-component of gravitational acceleration
!  e_nu_c       : neutrino energy density [ergs cm^{-3}]
!
!    Output arguments:
!
!  rho_z        : densities (g cm^{-3})
!  t_z          : temperatures (K)
!  ye_z         : electron fractions
!  ei_z         : internal energies (ergs g^{-1})
!  u_z          : velocity x-components (cm s^{-1})
!  v_z          : velocity y-components (cm s^{-1})
!  w_z          : velocity z-components (cm s^{-1})
!  psi0_z       : zero angular moments of the neutrino occupation number
!  psi0_z       : first angular moments of the neutrino occupation number
!  xn_z         : initial mass fractions
!  be_nuc_rep_z : binding energy of mean heavy nucleus
!  a_nuc_rep_z  : mass number of mean heavy nucleus
!  z_nuc_rep_z  : charge number of mean heavy nucleus
!  flat_x_z     : transposed variable indicating the presence of radial a shock
!  agr_z        : zone-centered value of the lapse function
!  grav_z_cz    : zone-centered z-component of gravitational acceleration
!  e_nu_z       : neutrino energy density [ergs cm^{-3}]
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

USE parallel_module, ONLY : myid, ierr, MPI_COMM_COL
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

INTEGER, INTENT(in)              :: kmin               ! lower z-array index
INTEGER, INTENT(in)              :: kmax               ! upper z-array index
INTEGER, INTENT(in)              :: nz                 ! z-array extent
INTEGER, INTENT(in)              :: k_ray_dim          ! number of radial zones on a processor after swapping with y

INTEGER, INTENT(in)              :: ny                 ! y-array extent
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
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: grav_z_c     ! zone-centered y-component of gravitational acceleration

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)                   :: nse_z        ! nse flag

REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: rho_z        ! density (cm^{-3})
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: t_z          ! temperature (MeV)
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: ye_z         ! temperature (MeV)
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: ei_z         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: u_z          ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: v_z          ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: w_z          ! z-velocity of zone (cm)

REAL(KIND=double), INTENT(out), DIMENSION(nz,nez,nnu,ij_ray_dim,k_ray_dim) :: psi0_z       ! zeroth angular moment of the NDS
REAL(KIND=double), INTENT(out), DIMENSION(nz,nez,nnu,ij_ray_dim,k_ray_dim) :: psi1_z       ! first angular moment of the NDS

REAL(KIND=double), INTENT(out), DIMENSION(nz,nnc,ij_ray_dim,k_ray_dim)     :: xn_z         ! mass fractions
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: a_nuc_rep_z  ! nuclear mass number
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: z_nuc_rep_z  ! nuclear charge number
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: be_nuc_rep_z  ! nuclear binding energy
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: flat_x_z     ! variable indicating the presence of radial a shock
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: agr_z        ! zone-centered value of the lapse function
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: grav_z_cz    ! zone-centered y-component of gravitational acceleration
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: e_nu_z       ! neutrino energy density [ergs cm^{-3}]

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
INTEGER                        :: npz           ! z-processor index
INTEGER                        :: ij_ray        ! k-index of a radial ray
INTEGER                        :: ik_ray        ! ie-index of a radial ray
INTEGER                        :: count_ata     ! all to all buffer count
INTEGER                        :: c_gath_send   ! gather send buffer count
INTEGER                        :: c_gath_recv   ! gather recv buffer count
INTEGER                        :: i_extent      ! broadcast array extent

INTEGER, DIMENSION(ij_ray_dim,1,ik_ray_dim,imax)                                 :: i_send ! All to all send buffer
INTEGER, DIMENSION(ij_ray_dim,1,ik_ray_dim,k_ray_dim,n_proc_z)                   :: i_recv ! All to all receive buffer
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z)                      :: i_gath ! receive buffer for gathering data from processors

REAL(KIND=double), DIMENSION(ij_ray_dim,14+nnc,ik_ray_dim,imax)                  :: send   ! All to all send buffer
REAL(KIND=double), DIMENSION(ij_ray_dim,14+nnc,ik_ray_dim,k_ray_dim,n_proc_z)    :: recv   ! All to all receive buffer
REAL(KIND=double), DIMENSION(ij_ray_dim,2,nez,nnu,ik_ray_dim,imax)               :: p_send ! All to all send buffer
REAL(KIND=double), DIMENSION(ij_ray_dim,2,nez,nnu,ik_ray_dim,k_ray_dim,n_proc_z) :: p_recv ! All to all receive buffer

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                    \\\\\ LOAD Z-ARRAYS /////
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

  count_ata                          = ( 14+nnc ) * ik_ray_dim * k_ray_dim * ij_ray_dim

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim
      DO i = imin,imax
        send(ij_ray, 1,ik_ray,i)     = rho_c       (i,ij_ray,ik_ray)
        send(ij_ray, 2,ik_ray,i)     = t_c         (i,ij_ray,ik_ray)
        send(ij_ray, 3,ik_ray,i)     = ye_c        (i,ij_ray,ik_ray)
        send(ij_ray, 4,ik_ray,i)     = ei_c        (i,ij_ray,ik_ray)
        send(ij_ray, 5,ik_ray,i)     = u_c         (i,ij_ray,ik_ray)
        send(ij_ray, 6,ik_ray,i)     = v_c         (i,ij_ray,ik_ray)
        send(ij_ray, 7,ik_ray,i)     = w_c         (i,ij_ray,ik_ray)
        send(ij_ray, 8,ik_ray,i)     = a_nuc_rep_c (i,ij_ray,ik_ray)
        send(ij_ray, 9,ik_ray,i)     = z_nuc_rep_c (i,ij_ray,ik_ray)
        send(ij_ray,10,ik_ray,i)     = be_nuc_rep_c(i,ij_ray,ik_ray)
        send(ij_ray,11,ik_ray,i)     = flat_x      (i,ij_ray,ik_ray)
        send(ij_ray,12,ik_ray,i)     = agr_c       (i,ij_ray,ik_ray)
        send(ij_ray,13,ik_ray,i)     = grav_z_c    (i,ij_ray,ik_ray)
        send(ij_ray,14,ik_ray,i)     = e_nu_c      (i,ij_ray,ik_ray)
        DO l = ls,le
          send(ij_ray,14+l,ik_ray,i) = xn_c(i,l,ij_ray,ik_ray)
        END DO ! l = ls,le
      END DO ! i = imin,imax
    END DO ! ij_ray = 1,ij_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

  CALL MPI_ALLTOALL( send, count_ata, MPI_DOUBLE_PRECISION, recv, count_ata, &
&   MPI_DOUBLE_PRECISION, MPI_COMM_COL, ierr )

  DO ij_ray = 1,ij_ray_dim
    DO i = 1,k_ray_dim
      DO ik_ray = 1,ik_ray_dim
        DO npz = 1,n_proc_z
          k                          = ik_ray + ik_ray_dim * ( npz - 1 )
          rho_z       (k,ij_ray,i)   = recv(ij_ray, 1,ik_ray,i,npz)
          t_z         (k,ij_ray,i)   = recv(ij_ray, 2,ik_ray,i,npz)
          ye_z        (k,ij_ray,i)   = recv(ij_ray, 3,ik_ray,i,npz)
          ei_z        (k,ij_ray,i)   = recv(ij_ray, 4,ik_ray,i,npz)
          u_z         (k,ij_ray,i)   = recv(ij_ray, 5,ik_ray,i,npz)
          v_z         (k,ij_ray,i)   = recv(ij_ray, 6,ik_ray,i,npz)
          w_z         (k,ij_ray,i)   = recv(ij_ray, 7,ik_ray,i,npz)
          a_nuc_rep_z (k,ij_ray,i)   = recv(ij_ray, 8,ik_ray,i,npz)
          z_nuc_rep_z (k,ij_ray,i)   = recv(ij_ray, 9,ik_ray,i,npz)
          be_nuc_rep_z(k,ij_ray,i)   = recv(ij_ray,10,ik_ray,i,npz)
          flat_x_z    (k,ij_ray,i)   = recv(ij_ray,11,ik_ray,i,npz)
          agr_z       (k,ij_ray,i)   = recv(ij_ray,12,ik_ray,i,npz)
          grav_z_cz   (k,ij_ray,i)   = recv(ij_ray,13,ik_ray,i,npz)
          e_nu_z      (k,ij_ray,i)   = recv(ij_ray,14,ik_ray,i,npz)
          DO l = ls,le
            xn_z      (k,l,ij_ray,i) = recv(ij_ray,14+l,ik_ray,i,npz)
          END DO ! l = ls,le
        END DO ! npz = 1,n_proc_z
      END DO ! ik_ray = 1,ik_ray_dim
    END DO ! i = 1,k_ray_dim
  END DO ! ij_ray = 1,ij_ray_dim

  count_ata                          = ik_ray_dim * k_ray_dim * ij_ray_dim

  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim
      DO i = imin,imax
        i_send(ij_ray,1,ik_ray,i)    = nse_c(i,ij_ray,ik_ray)
       END DO ! i = imin,imax
    END DO ! ij_ray = 1,ij_ray_dim
  END DO ! ik_ray = 1,ik_ray_dim

  CALL MPI_ALLTOALL( i_send, count_ata, MPI_INTEGER, i_recv, count_ata, &
&   MPI_INTEGER, MPI_COMM_COL, ierr )

  DO ij_ray = 1,ij_ray_dim
    DO i = 1,k_ray_dim
      DO ik_ray = 1,ik_ray_dim
        DO npz = 1,n_proc_z
          k                          = ik_ray + ik_ray_dim * ( npz - 1 )
          nse_z(k,ij_ray,i)          = i_recv(ij_ray, 1,ik_ray,i,npz)
        END DO ! npz = 1,n_proc_z
      END DO ! ik_ray = 1,ik_ray_dim
    END DO ! i = 1,k_ray_dim
  END DO ! ij_ray = 1,ij_ray_dim

  IF ( nnugpmx /= 0 ) THEN

    count_ata                        = 2 * nez * nnu * ik_ray_dim * k_ray_dim * ij_ray_dim

    DO ik_ray = 1,ik_ray_dim
      DO ij_ray = 1,ij_ray_dim
        DO i = imin,imax
          DO n = 1,nnu
            IF ( nnugp(n) /= 0 ) THEN
              DO ie = 1,nnugp(n)
                p_send(ij_ray,1,ie,n,ik_ray,i) = psi0_c(i,ie,n,ij_ray,ik_ray)
                p_send(ij_ray,2,ie,n,ik_ray,i) = psi1_e(i,ie,n,ij_ray,ik_ray)
              END DO ! ie = 1,nnugp(n)
            END IF ! nnugp(n) /= 0
          END DO ! n = 1,nnu
        END DO ! i = imin,imax
      END DO ! ij_ray = 1,ij_ray_dim
    END DO ! ik_ray = 1,ik_ray_dim

    CALL MPI_ALLTOALL( p_send, count_ata, MPI_DOUBLE_PRECISION, p_recv, count_ata, &
&     MPI_DOUBLE_PRECISION, MPI_COMM_COL, ierr )

    DO ij_ray = 1,ij_ray_dim
      DO i = 1,k_ray_dim
        DO ik_ray = 1,ik_ray_dim
          DO npz = 1,n_proc_z
            k                        = ik_ray + ik_ray_dim * ( npz - 1 )
            DO n = 1,nnu
              IF ( nnugp(n) /= 0 ) THEN
                DO ie = 1,nnugp(n)
                  psi0_z(k,ie,n,ij_ray,i) = p_recv(ij_ray,1,ie,n,ik_ray,i,npz)
                  psi1_z(k,ie,n,ij_ray,i) = p_recv(ij_ray,2,ie,n,ik_ray,i,npz)
                END DO ! ie = 1,nnugp(n)
              END IF ! nnugp(n) /= 0
            END DO ! n = 1,nnu
          END DO ! npz = 1,n_proc_z
        END DO ! ik_ray = 1,ik_ray_dim
      END DO ! i = 1,k_ray_dim
    END DO ! ij_ray = 1,ij_ray_dim

  END IF ! nnugpmx /= 0

!-----------------------------------------------------------------------
!
!            \\\\\ TRANSPOSE VARIABLES IF N_PROC = 1 /////
!
!-----------------------------------------------------------------------

ELSE ! n_procs = 1

  DO i = imin,imax
    DO k = kmin,kmax
      rho_z       (k,:,i)            = rho_c       (i,:,k)
      t_z         (k,:,i)            = t_c         (i,:,k)
      ye_z        (k,:,i)            = ye_c        (i,:,k)
      ei_z        (k,:,i)            = ei_c        (i,:,k)
      u_z         (k,:,i)            = u_c         (i,:,k)
      v_z         (k,:,i)            = v_c         (i,:,k)
      w_z         (k,:,i)            = w_c         (i,:,k)
      a_nuc_rep_z (k,:,i)            = a_nuc_rep_c (i,:,k)
      z_nuc_rep_z (k,:,i)            = z_nuc_rep_c (i,:,k)
      be_nuc_rep_z(k,:,i)            = be_nuc_rep_c(i,:,k)
      nse_z       (k,:,i)            = nse_c       (i,:,k)
      flat_x_z    (k,:,i)            = flat_x      (i,:,k)
      agr_z       (k,:,i)            = agr_c       (i,:,k)
      grav_z_cz   (k,:,i)            = grav_z_c    (i,:,k)
      e_nu_z      (k,:,i)            = e_nu_c      (i,:,k)
    END DO ! k = kmin,kmax
  END DO ! i = imin,imax

  DO l = ls,le
    DO i = imin,imax
      DO k = kmin,kmax
        xn_z(k,l,:,i)                = xn_c(i,l,:,k)
      END DO ! k = kmin,kmax
    END DO ! i = imin,imax
  END DO ! l = ls,le

  DO n = 1,nnu
    IF ( nnugp(n) == 0 ) CYCLE
    DO ie = 1,nnugp(n)
      DO i = imin,imax
        DO k = kmin,kmax
          psi0_z(k,ie,n,:,i)         = psi0_c(i,ie,n,:,k) 
          psi1_z(k,ie,n,:,i)         = psi1_e(i,ie,n,:,k) 
        END DO ! k = kmin,kmax
      END DO ! i = imin,imax
    END DO ! ie = 1,nnugp(n)
  END DO ! n = 1,nnu

END IF ! ! n_procs > 1

RETURN
END SUBROUTINE transpose_x_z
