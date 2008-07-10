 SUBROUTINE angular_ave( nx, ny, nz, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         angular_ave
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
!  ny         : y-array extent (used only in MPI version)
!  nz         : z-array extent (used only in MPI version)
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module,, numerical_module
!  edit_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero, half, frpi, third, epsilon

USE edit_module, ONLY : nprint, nlog
USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, x_ei, &
& dx_ci, y_ei, y_ci, d_omega, rho_c, rhobar, aesv_c, p_bar, e_bar, u_c, &
& v_c, w_c, vx_bar, v2_bar, e_nu_c, e_nu_c_bar, f_nu_e, f_nu_e_bar, y_shft, &
& x_vel_ns, y_vel_ns, z_vel_ns, vel_ns, d_x_vel_ns, d_y_vel_ns, d_z_vel_ns, &
& d_vel_ns, mass_ns, G_trns, cos_theta, sin_theta, cos_phi, sin_phi

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                        :: nx             ! x-array extent
INTEGER, INTENT(in)                        :: ny             ! y-array extent (used only in MPI version)
INTEGER, INTENT(in)                        :: nz             ! z-array extent (used only in MPI version)
INTEGER, INTENT(inout)                     :: ij_ray_dim     ! number of y-zones on a processor before swapping
INTEGER, INTENT(inout)                     :: ik_ray_dim     ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, SAVE                              :: first = .true.

INTEGER                                    :: i              ! x-array index
INTEGER                                    :: j              ! y-array index
INTEGER                                    :: k              ! z-array index

REAL(KIND=double), PARAMETER               :: rho_ns = 1.d+11
REAL(KIND=double), PARAMETER               :: rho_ns_max = 2.5d+14
REAL(KIND=double)                          :: d_vol          ! radial zone volume
REAL(KIND=double)                          :: d_mass         ! zome mass
REAL(KIND=double)                          :: mass_gas       ! mass of the gas exterior to the neutron star
REAL(KIND=double)                          :: x_Momentum_gas ! x-momentum of the gas exterior to the neutron star
REAL(KIND=double)                          :: y_Momentum_gas ! y-momentum of the gas exterior to the neutron star
REAL(KIND=double)                          :: z_Momentum_gas ! z-momentum of the gas exterior to the neutron star
REAL(KIND=double), SAVE                    :: x_Momentum_i   ! initial x-momentum of the gas exterior to the neutron star
REAL(KIND=double), SAVE                    :: y_Momentum_i   ! initial y-momentum of the gas exterior to the neutron star
REAL(KIND=double), SAVE                    :: z_Momentum_i   ! initial z-momentum of the gas exterior to the neutron star

REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: v2 ! KE per unit mass (ergs g^{-1})

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Need MPI version of angular_ave since n_proc=',i4,' > 1')

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
!  Initialize
!-----------------------------------------------------------------------

v2                      = zero
v2(imin:imax,:,:)       = u_c(imin:imax,:,:) * u_c(imin:imax,:,:) &
&                       + v_c(imin:imax,:,:) * v_c(imin:imax,:,:) &
&                       + w_c(imin:imax,:,:) * w_c(imin:imax,:,:)

!-----------------------------------------------------------------------
!  Calculate rhobar
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
  rhobar(imin:imax)     = rho_c(imin:imax,1,1)
ELSE ! jmax /= jmin and/or kmax /= kmin
  DO i = imin,imax
    rhobar(i)           = SUM( rho_c(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i = imin,imax
END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate p_bar
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    p_bar(imin:imax)    = aesv_c(imin:imax,1,1,1)
ELSE ! jmax /= jmin and/or kmax /= kmin
  DO i = imin,imax
    p_bar(i)            = SUM( aesv_c(i,1,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i = imin,imax
END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate e_bar
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
    e_bar(imin:imax)    = aesv_c(imin:imax,2,1,1)
ELSE ! jmax /= jmin and/or kmax /= kmin
  DO i = imin,imax
    e_bar(i)            = SUM( aesv_c(i,2,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i = imin,imax
END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate vx_bar
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
  vx_bar(imin:imax)     = u_c(imin:imax,1,1)
ELSE ! jmax /= jmin and/or kmax /= kmin
  DO i = imin,imax
    vx_bar(i)           = SUM( u_c(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i = imin,imax
END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate v2_bar
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
  v2_bar(imin:imax)     = v2(imin:imax,1,1)
ELSE ! jmax /= jmin and/or kmax /= kmin
  DO i = imin,imax
    v2_bar(i)           = SUM( v2(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i = imin,imax
END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate e_nu_c_bar
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
  e_nu_c_bar(imin:imax) = e_nu_c(imin:imax,1,1)
ELSE ! jmax /= jmin and/or kmax /= kmin
  DO i = imin,imax
    e_nu_c_bar(i)       = SUM( e_nu_c(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i = imin,imax
END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!  Calculate f_nu_e_bar
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
  f_nu_e_bar(imin:imax) = f_nu_e(imin:imax,1,1)
ELSE ! jmax /= jmin and/or kmax /= kmin
  DO i = imin,imax+1
    f_nu_e_bar(i)       = SUM( f_nu_e(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i = imin,imax
END IF ! jmax == jmin and kmin == kmax

!-----------------------------------------------------------------------
!
!                  \\\\\ NEUTRON STAR KICK /////
!
!-----------------------------------------------------------------------

IF ( jmax == jmin  .and.  kmin == kmax ) THEN
  x_vel_ns                   = zero
  y_vel_ns                   = zero
  z_vel_ns                   = zero
  vel_ns                     = zero
ELSE ! jmax /= jmin  or  kmax /= kmin

  mass_ns                    = zero
  d_x_vel_ns                 = zero
  d_y_vel_ns                 = zero
  d_z_vel_ns                 = zero
  d_vel_ns                   = zero
  x_Momentum_gas             = zero
  y_Momentum_gas             = zero
  z_Momentum_gas             = zero

!-----------------------------------------------------------------------
!  Calculate the mass of the neutron star if rhobar(imin) > rho_ns_max
!-----------------------------------------------------------------------

  IF ( rhobar(imin) > rho_ns_max ) THEN

    DO i = imin,imax
      IF ( rhobar(i) <= rho_ns ) CYCLE
      d_vol                  = frpi * dx_ci(i) * ( x_ei(i) * x_ei(i+1) + dx_ci(i) * dx_ci(i) * third )
      mass_ns                = mass_ns + d_vol * rhobar(i)
    END DO ! i = imin,imax

!-----------------------------------------------------------------------
!  Calculate the change momentum of the material external to the
!   neutron star, and from that the change in momentum and thence the
!   change in velocity of the neutron star
!-----------------------------------------------------------------------

    mass_gas                 = zero
    DO i = imin,imax
      IF ( rhobar(i) > rho_ns ) CYCLE
      d_vol                  = frpi * dx_ci(i) * ( x_ei(i) * x_ei(i+1) + dx_ci(i) * dx_ci(i) * third )
      DO k = kmin,kmax
        DO j = jmin,jmax
          d_mass             = rho_c(i,j,k) * d_vol * d_omega(j,k)/frpi
          mass_gas           = mass_gas + d_mass
          x_Momentum_gas     = x_Momentum_gas + ( u_c(i,j,k) * sin_theta(j) * cos_phi(k)   &
&                                             +   v_c(i,j,k) * cos_theta(j) * cos_phi(k)   &
&                                             -   w_c(i,j,k) * sin_phi(k)                ) * d_mass
          y_Momentum_gas     = y_Momentum_gas + ( u_c(i,j,k) * sin_theta(j) * sin_phi(k)   &
&                                             +   v_c(i,j,k) * cos_theta(j) * sin_phi(k)   &
&                                             -   w_c(i,j,k) * cos_phi(k)                ) * d_mass
          z_Momentum_gas     = z_Momentum_gas + ( u_c(i,j,k) * cos_theta(j)                &
&                                             -   v_c(i,j,k) * sin_theta(j)              ) * d_mass
        END DO ! j = jmin,jmax
      END DO ! k = kmin,kmax
    END DO ! i = imin,imax

    IF ( first ) THEN
      x_Momentum_i           = x_Momentum_gas
      y_Momentum_i           = y_Momentum_gas
      z_Momentum_i           = z_Momentum_gas
      first                  = .false.
    END IF ! first
    d_x_vel_ns               = - ( x_Momentum_gas - x_Momentum_i )/mass_ns
    d_y_vel_ns               = - ( y_Momentum_gas - y_Momentum_i )/mass_ns
    d_z_vel_ns               = - ( z_Momentum_gas - z_Momentum_i )/mass_ns
    x_Momentum_i             = x_Momentum_gas - mass_gas * d_x_vel_ns
    y_Momentum_i             = y_Momentum_gas - mass_gas * d_y_vel_ns
    z_Momentum_i             = z_Momentum_gas - mass_gas * d_z_vel_ns

  END IF ! rhobar(imin) > rho_ns

END IF ! jmax == jmin  .and.  kmin == kmax

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
  DO k = kmin,kmax
    DO j = jmin,jmax
      DO i = imin,imax
        IF ( rhobar(i) > rho_ns ) CYCLE
        u_c(i,j,k)           = u_c(i,j,k) - d_x_vel_ns * sin_theta(j) * cos_phi(k) &
&                                         - d_y_vel_ns * sin_theta(j) * sin_phi(k) &
&                                         - d_z_vel_ns * cos_theta(j)
        v_c(i,j,k)           = v_c(i,j,k) - d_x_vel_ns * cos_theta(j) * cos_phi(k) &
&                                         - d_y_vel_ns * cos_theta(j) * sin_phi(k) &
&                                         + d_z_vel_ns * sin_theta(j)
        w_c(i,j,k)           = w_c(i,j,k) + d_x_vel_ns * sin_phi(k)                &
&                                         - d_y_vel_ns * cos_phi(k)
      END DO ! i = imin,imax
    END DO ! j = jmin,jmax
  END DO ! k = kmin,kmax
END IF ! G_trns == 'no' 

RETURN
END SUBROUTINE angular_ave
