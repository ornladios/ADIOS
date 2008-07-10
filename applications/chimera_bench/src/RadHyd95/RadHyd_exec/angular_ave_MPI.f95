 SUBROUTINE angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         angular_ave_MPI
!    Module:       angular_ave
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/13/07
!
!    Purpose:
!      To calculate the angular average of selected quantities.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x-array extent
!  ny         : y-array extent
!  nz         : z-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, parallel_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc, n_proc_y, n_proc_z
USE numerical_module, ONLY : zero, half, frpi, third, epsilon

USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid, myid_y, myid_z, ierr
USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, x_ei, &
& dx_ci, y_ei, y_ci, d_omega, rho_c, rhobar, aesv_c, p_bar, e_bar, u_c, &
& v_c, w_c, vx_bar, v2_bar, e_nu_c, e_nu_c_bar, f_nu_e, f_nu_e_bar, y_shft, &
& x_vel_ns, y_vel_ns, z_vel_ns, vel_ns, d_x_vel_ns, d_y_vel_ns, d_z_vel_ns, &
& d_vel_ns, mass_ns, G_trns, cos_theta, sin_theta, cos_phi, sin_phi

USE mpi

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                            :: nx             ! x-array extent
INTEGER, INTENT(in)                            :: ny             ! y-array extent (used only in MPI version)
INTEGER, INTENT(in)                            :: nz             ! z-array extent (used only in MPI version)
INTEGER, INTENT(inout)                         :: ij_ray_dim     ! number of y-zones on a processor before swapping
INTEGER, INTENT(inout)                         :: ik_ray_dim     ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, SAVE                                  :: first = .true.

INTEGER                                        :: i              ! x-array index
INTEGER                                        :: j              ! y-array index
INTEGER                                        :: ij_ray         ! j-index of a radial ray
INTEGER                                        :: k              ! z-array index
INTEGER                                        :: ik_ray         ! k-index of a radial ray
INTEGER                                        :: m              ! processor index
INTEGER                                        :: mj             ! y-block index
INTEGER                                        :: mk             ! z-block index
INTEGER                                        :: jsk            ! y-array index of gathered array
INTEGER                                        :: ksk            ! z-array index of gathered array
INTEGER                                        :: c_gath_send    ! gather send buffer count
INTEGER                                        :: c_gath_recv    ! gather recv buffer count
INTEGER                                        :: i_extent       ! broadcast array extent

REAL(KIND=double), PARAMETER                   :: rho_ns = 1.d+11
REAL(KIND=double), PARAMETER                   :: rho_ns_max = 2.5d+14
REAL(KIND=double)                              :: d_vol          ! radial zone volume
REAL(KIND=double)                              :: d_mass         ! zome mass
REAL(KIND=double)                              :: mass_gas       ! mass of the gas exterior to the neutron star
REAL(KIND=double)                              :: x_Momentum_gas ! x-momentum of the gas exterior to the neutron star
REAL(KIND=double)                              :: y_Momentum_gas ! y-momentum of the gas exterior to the neutron star
REAL(KIND=double)                              :: z_Momentum_gas ! z-momentum of the gas exterior to the neutron star
REAL(KIND=double), SAVE                        :: x_Momentum_i   ! initial x-momentum of the gas exterior to the neutron star
REAL(KIND=double), SAVE                        :: y_Momentum_i   ! initial y-momentum of the gas exterior to the neutron star
REAL(KIND=double), SAVE                        :: z_Momentum_i   ! initial z-momentum of the gas exterior to the neutron star

REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: v2     ! density of radial zone i, angular zone j (g cm^{-3})

REAL(KIND=double), DIMENSION(nx,ny,nz)         :: rho_all        ! density of radial zone i, angular zone j (g cm^{-3})
REAL(KIND=double), DIMENSION(nx,ny,nz)         :: p_all          ! pressure of radial zone i, angular zone j (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)         :: e_all          ! internal wnergy of radial zone i, angular zone j (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)         :: vx_all         ! x (radial) velocity of zone i, j, k (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)         :: vy_all         ! y (amgular) velocity of zone i, j, k (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)         :: vz_all         ! z (azimuthal) velocity of zone i, j, k (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)         :: v2_all         ! velocity squared of radial zone i, angular zone j (cm^{2} s^{-2})
REAL(KIND=double), DIMENSION(nx,ny,nz)         :: e_nu_all       ! neutrino energy density of radial zone i, angular zone j (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx+1,ny,nz)       :: f_nu_all       ! neutrino flux at edge of radial zone i, center ofangular zone j (ergs cm^{-3})
REAL(KIND=double), DIMENSION(9*nx+1,ij_ray_dim,ik_ray_dim) :: send_buf2   ! send buffer for gathering data from processors
REAL(KIND=double), DIMENSION(9*nx+1,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z) :: recv_buf2 ! receive buffer for gathering data from processors
REAL(KIND=double), DIMENSION(7*nx+9)           :: send_buf1      ! send buffer for gathering data from processors

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                \\\\\ GATHER QUANTITIES TO AVERAGE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

v2                               = zero
v2(imin:imax,:,:)                = u_c(imin:imax,:,:) * u_c(imin:imax,:,:) &
&                                + v_c(imin:imax,:,:) * v_c(imin:imax,:,:) &
&                                + w_c(imin:imax,:,:) * w_c(imin:imax,:,:)

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

send_buf2(     1:  nx  ,:,:)     = rho_c (1:nx  ,:,:)
send_buf2(  nx+1:2*nx  ,:,:)     = aesv_c(1:nx,1,:,:)
send_buf2(2*nx+1:3*nx  ,:,:)     = aesv_c(1:nx,2,:,:)
send_buf2(3*nx+1:4*nx  ,:,:)     = u_c   (1:nx  ,:,:)
send_buf2(4*nx+1:5*nx  ,:,:)     = v_c   (1:nx  ,:,:)
send_buf2(5*nx+1:6*nx  ,:,:)     = w_c   (1:nx  ,:,:)
send_buf2(6*nx+1:7*nx  ,:,:)     = v2    (1:nx  ,:,:)
send_buf2(7*nx+1:8*nx  ,:,:)     = e_nu_c(1:nx  ,:,:)
send_buf2(8*nx+1:9*nx+1,:,:)     = f_nu_e(1:nx+1,:,:)

c_gath_send                      = ( 9 * nx + 1 ) * ij_ray_dim * ik_ray_dim
c_gath_recv                      = ( 9 * nx + 1 ) * ij_ray_dim * ik_ray_dim

!-----------------------------------------------------------------------
!  Gather quantities to average
!-----------------------------------------------------------------------

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( send_buf2, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf2, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Unpack receive buffer
!-----------------------------------------------------------------------

  DO m = 0,n_proc-1
    mj                           = MOD( m, n_proc_y )
    mk                           = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                        = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk                      = mj * ij_ray_dim + j
        rho_all (1:nx,  jsk,ksk) = recv_buf2(1:nx,         j,k,m+1)
        p_all   (1:nx,  jsk,ksk) = recv_buf2(nx+1:2*nx,    j,k,m+1)
        e_all   (1:nx,  jsk,ksk) = recv_buf2(2*nx+1:3*nx,  j,k,m+1)
        vx_all  (1:nx,  jsk,ksk) = recv_buf2(3*nx+1:4*nx,  j,k,m+1)
        vy_all  (1:nx,  jsk,ksk) = recv_buf2(4*nx+1:5*nx,  j,k,m+1)
        vz_all  (1:nx,  jsk,ksk) = recv_buf2(5*nx+1:6*nx,  j,k,m+1)
        v2_all  (1:nx,  jsk,ksk) = recv_buf2(6*nx+1:7*nx,  j,k,m+1)
        e_nu_all(1:nx,  jsk,ksk) = recv_buf2(7*nx+1:8*nx,  j,k,m+1)
        f_nu_all(1:nx+1,jsk,ksk) = recv_buf2(8*nx+1:9*nx+1,j,k,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!  Calculate rhobar
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    rhobar(imin:imax)            = rho_c(imin:imax,1,1)
  ELSE ! jmax /= jmin and/or kmax /= kmin
    DO i = imin,imax
      rhobar(i)                  = SUM( rho_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i = imin,imax
  END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate p_bar
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    p_bar(imin:imax)             = aesv_c(imin:imax,1,1,1)
  ELSE ! jmax /= jmin and/or kmax /= kmin
    DO i = imin,imax
      p_bar(i)                   = SUM( p_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i = imin,imax
  END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate e_bar
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    e_bar(imin:imax)             = aesv_c(imin:imax,2,1,1)
  ELSE ! jmax /= jmin and/or kmax /= kmin
    DO i = imin,imax
      e_bar(i)                   = SUM( e_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i = imin,imax
  END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate vx_bar
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    vx_bar(imin:imax)            = u_c(imin:imax,1,1)
  ELSE ! jmax /= jmin and/or kmax /= kmin
    DO i = imin,imax
      vx_bar(i)                  = SUM( vx_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i = imin,imax
  END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate v2_bar
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    v2_bar(imin:imax)            = v2(imin:imax,1,1)
  ELSE ! jmax /= jmin and/or kmax /= kmin
    DO i = imin,imax
      v2_bar(i)                  = SUM( v2_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i = imin,imax
  END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate e_nu_c_bar
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    e_nu_c_bar(imin:imax)        = e_nu_c(imin:imax,1,1)
  ELSE ! jmax /= jmin and/or kmax /= kmin
    DO i = imin,imax
      e_nu_c_bar(i)              = SUM( e_nu_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i = imin,imax
  END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate f_nu_e_bar
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    f_nu_e_bar(imin:imax)        = f_nu_e(imin:imax,1,1)
  ELSE ! jmax /= jmin and/or kmax /= kmin
    DO i = imin,imax+1
      f_nu_e_bar(i)              = SUM( f_nu_all(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! i = imin,imax
  END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!
!                  \\\\\ NEUTRON STAR KICK /////
!
!-----------------------------------------------------------------------

  IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    x_vel_ns                     = zero
    y_vel_ns                     = zero
    z_vel_ns                     = zero
    vel_ns                       = zero
  ELSE ! jmax /= jmin  or  kmax /= kmin

    mass_ns                      = zero
    d_x_vel_ns                   = zero
    d_y_vel_ns                   = zero
    d_z_vel_ns                   = zero
    d_vel_ns                     = zero
    x_Momentum_gas               = zero
    y_Momentum_gas               = zero
    z_Momentum_gas               = zero

!-----------------------------------------------------------------------
!  Calculate the mass of the neutron star if rhobar(imin) > rho_ns_max
!-----------------------------------------------------------------------

    IF ( rhobar(imin) > rho_ns_max ) THEN

      DO i = imin,imax
        IF ( rhobar(i) <= rho_ns ) CYCLE
        d_vol                    = frpi * dx_ci(i) * ( x_ei(i) * x_ei(i+1) + dx_ci(i) * dx_ci(i) * third )
        mass_ns                  = mass_ns + d_vol * rhobar(i)
      END DO ! i = imin,imax

!-----------------------------------------------------------------------
!  Calculate the change momentum of the material external to the
!   neutron star, and from that the change in momentum and thence the
!   change in velocity of the neutron star
!-----------------------------------------------------------------------

      mass_gas                   = zero
      DO i = imin,imax
        IF ( rhobar(i) > rho_ns ) CYCLE
        d_vol                    = frpi * dx_ci(i) * ( x_ei(i) * x_ei(i+1) + dx_ci(i) * dx_ci(i) * third )
        DO k = kmin,kmax
          DO j = jmin,jmax
            d_mass               = rho_all(i,j,k) * d_vol * d_omega(j,k)/frpi
            mass_gas             = mass_gas + d_mass
            x_Momentum_gas       = x_Momentum_gas + ( vx_all(i,j,k) * sin_theta(j) * cos_phi(k)   &
&                                               +   vy_all(i,j,k) * cos_theta(j) * cos_phi(k)     &
&                                               -   vz_all(i,j,k) * sin_phi(k)                ) * d_mass
            y_Momentum_gas       = y_Momentum_gas + ( vx_all(i,j,k) * sin_theta(j) * sin_phi(k)   &
&                                               +   vy_all(i,j,k) * cos_theta(j) * sin_phi(k)     &
&                                               -   vz_all(i,j,k) * cos_phi(k)                ) * d_mass
            z_Momentum_gas       = z_Momentum_gas + ( vx_all(i,j,k) * cos_theta(j)                &
&                                               -   vy_all(i,j,k) * sin_theta(j)              ) * d_mass
          END DO ! j = jmin,jmax
        END DO ! k = kmin,kmax
      END DO ! i = imin,imax

      IF ( first ) THEN
        x_Momentum_i             = x_Momentum_gas
        y_Momentum_i             = y_Momentum_gas
        z_Momentum_i             = z_Momentum_gas
        first                    = .false.
      END IF ! first
      d_x_vel_ns                 = - ( x_Momentum_gas - x_Momentum_i )/mass_ns
      d_y_vel_ns                 = - ( y_Momentum_gas - y_Momentum_i )/mass_ns
      d_z_vel_ns                 = - ( z_Momentum_gas - z_Momentum_i )/mass_ns
      x_Momentum_i               = x_Momentum_gas - mass_gas * d_x_vel_ns
      y_Momentum_i               = y_Momentum_gas - mass_gas * d_y_vel_ns
      z_Momentum_i               = z_Momentum_gas - mass_gas * d_z_vel_ns

    END IF ! rhobar(imin) > rho_ns

  END IF ! jmax == jmin  .and.  kmin == kmax

!-----------------------------------------------------------------------
!
!       \\\\\ BROADCAST MEAN DENSITIES TO ALL PROCESSORS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

  send_buf1(1:nx)                = rhobar    (1:nx)
  send_buf1(nx+1:2*nx)           = p_bar     (1:nx)
  send_buf1(2*nx+1:3*nx)         = e_bar     (1:nx)
  send_buf1(3*nx+1:4*nx)         = vx_bar    (1:nx)
  send_buf1(4*nx+1:5*nx)         = v2_bar    (1:nx)
  send_buf1(5*nx+1:6*nx)         = e_nu_c_bar(1:nx)
  send_buf1(6*nx+1:7*nx+1)       = f_nu_e_bar(1:nx+1)
  send_buf1(7*nx+2)              = mass_ns
  send_buf1(7*nx+3)              = d_x_vel_ns
  send_buf1(7*nx+4)              = d_y_vel_ns
  send_buf1(7*nx+5)              = d_z_vel_ns
  send_buf1(7*nx+6)              = d_vel_ns
  send_buf1(7*nx+7)              = x_Momentum_gas
  send_buf1(7*nx+8)              = y_Momentum_gas
  send_buf1(7*nx+9)              = z_Momentum_gas

!-----------------------------------------------------------------------
!               ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast angular averaged quantities
!-----------------------------------------------------------------------

i_extent                         = 7 * nx + 9
CALL MPI_BCAST( send_buf1, i_extent, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Unpack send_buf1 buffer
!-----------------------------------------------------------------------

rhobar    (1:nx)                 = send_buf1(1:nx)
p_bar     (1:nx)                 = send_buf1(nx+1:2*nx)
e_bar     (1:nx)                 = send_buf1(2*nx+1:3*nx)
vx_bar    (1:nx)                 = send_buf1(3*nx+1:4*nx)
v2_bar    (1:nx)                 = send_buf1(4*nx+1:5*nx)
e_nu_c_bar(1:nx)                 = send_buf1(5*nx+1:6*nx)
f_nu_e_bar(1:nx+1)               = send_buf1(6*nx+1:7*nx+1)
mass_ns                          = send_buf1(7*nx+2)
d_x_vel_ns                       = send_buf1(7*nx+3)
d_y_vel_ns                       = send_buf1(7*nx+4)
d_z_vel_ns                       = send_buf1(7*nx+5)
d_vel_ns                         = send_buf1(7*nx+6)
x_Momentum_gas                   = send_buf1(7*nx+7)
y_Momentum_gas                   = send_buf1(7*nx+8)
z_Momentum_gas                   = send_buf1(7*nx+9)

!-----------------------------------------------------------------------
!
!        \\\\\ PERFORM GALILEAN TRANSPORT IF G_TRNS = YE /////
!
!-----------------------------------------------------------------------

IF ( G_trns == 'no' ) THEN
  x_vel_ns                   = - x_Momentum_gas/( mass_ns + epsilon )
  y_vel_ns                   = - y_Momentum_gas/( mass_ns + epsilon )
  z_vel_ns                   = - z_Momentum_gas/( mass_ns + epsilon )
  vel_ns                     = DSQRT( x_vel_ns * x_vel_ns + y_vel_ns * y_vel_ns + z_vel_ns * z_vel_ns + epsilon )
ELSE
  x_vel_ns                   = x_vel_ns + d_x_vel_ns
  y_vel_ns                   = y_vel_ns + d_y_vel_ns
  z_vel_ns                   = z_vel_ns + d_z_vel_ns
  DO k = 1,ik_ray_dim
    ik_ray                   = ik_ray_dim * myid_z + k
    DO j = 1,ij_ray_dim
      ij_ray                 = ij_ray_dim * myid_y + j
      DO i = imin,imax
        IF ( rhobar(i) > rho_ns ) CYCLE
        u_c(i,j,k)           = u_c(i,j,k) - d_x_vel_ns * sin_theta(ij_ray) * cos_phi(ik_ray) &
&                                         - d_y_vel_ns * sin_theta(ij_ray) * sin_phi(ik_ray) &
&                                         - d_z_vel_ns * cos_theta(ij_ray)
        v_c(i,j,k)           = v_c(i,j,k) - d_x_vel_ns * cos_theta(ij_ray) * cos_phi(ik_ray) &
&                                         - d_y_vel_ns * cos_theta(ij_ray) * sin_phi(ik_ray) &
&                                         + d_z_vel_ns * sin_theta(ij_ray)
        w_c(i,j,k)           = w_c(i,j,k) + d_x_vel_ns * sin_phi(ik_ray)                     &
&                                         - d_y_vel_ns * cos_phi(ik_ray)
      END DO ! i = imin,imax
    END DO ! j = jmin,jmax
  END DO ! k = kmin,kmax
END IF ! G_trns == 'no' 

RETURN
END SUBROUTINE angular_ave
