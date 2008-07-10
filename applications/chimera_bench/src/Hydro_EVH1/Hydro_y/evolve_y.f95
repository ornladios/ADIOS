SUBROUTINE evolve_y( ngeom, umid, pmid, ji_ray, jk_ray )

!-----------------------------------------------------------------------
!
!    File:         evolve_y
!    Module:       evolve_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To update the densities, velocities, and energies for the angular
!       sweeps.
!
!    Input arguments:
!
!  ngeom            : geometry key
!  umid             : time and space averaged zone edged velocities
!  pmid             : time and space averaged zone edged pressures
!  ji_ray           : x (radial) index of a specific y (angular) ray
!  jk_ray           : z (azimuthal) index of a specific y (angular) ray
!
!    Subprograms called:
!  zone_center      : computes the zone-centered volume averaged values of the coordinates
!  forces_center    : computes the forces acting on a mass shell
!  pseudo_y         : computes the artificial viscosity (used for detecting the presence of a shock)
!  tgvndeye_sweep_y : computes the t's given the e's, rho's, and ye's
!  sweepbc          : computes the ghast values of the coordinates and other dependent variables
!  paraset          : computes the coefficients for parabolic interpolations on the advanced grid
!  e_decompose_e    : computes ei and ekin from e in the width of a shock
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  angular_ray_module, cycle_module, edit_module, eos_snc_y_module,
!  evh1_global,evh1_sweep, evh1_zone, parallel_module, prb_cntl_module,
!  shock_module
!
!  Use umid and pmid from Riemann solver to update velocity, density,
!   and total energy.
!  Physical zones are from nmin to nmax.  Zone boundary numbers run from
!   nmin to nmax+1
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: max_12, j_ray_dim
USE numerical_module, ONLY : zero, third, half, epsilon 
USE physcnst_module, ONLY : pi

USE angular_ray_module, ONLY : v_e
USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nlog
USE eos_snc_y_module, ONLY : aesv
USE evh1_global, ONLY : dt, small, smallr, nlefty, nrighty
USE evh1_sweep, ONLY : u, v, w, r, xa, dvol, e, ei, dx, nmin, nmax, &
& radius, p, nshk, p_nu, e_nu
USE evh1_zone, ONLY : jmax, zparay
USE parallel_module, ONLY : myid_y
USE prb_cntl_module, ONLY : v_trans_0
USE shock_module, ONLY : pq_y, lshock, j_shk_radial_all_p

     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                   :: ngeom    ! geometry index
INTEGER, INTENT(in)                   :: ji_ray   ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)                   :: jk_ray   ! z (azimuthal) index of a specific y (angular) ray


REAL(KIND=double), DIMENSION(max_12), INTENT(in)    :: pmid ! time averaged pressure at the zone interface

!-----------------------------------------------------------------------
!        Input-output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(max_12), INTENT(inout) :: umid ! time averaged velocity at the zone interface

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                               :: n           ! padded zone index
INTEGER                               :: ntot        ! number of zones (real plus ghost)
INTEGER                               :: j_radial    ! the shifted radial zone (angular ray) corresponding to ji_ray, jk_ray

REAL(KIND=double), DIMENSION(max_12)  :: amid        ! average area
REAL(KIND=double), DIMENSION(max_12)  :: upmid       ! umid*pmid
REAL(KIND=double), DIMENSION(max_12)  :: dtbdm       ! dt/dm

REAL(KIND=double), DIMENSION(max_12)  :: uold        ! velocity before the Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: xa1         ! xa before Lagrangian update
REAL(KIND=double), DIMENSION(max_12)  :: xa2         ! zone center before Lagrangian update
REAL(KIND=double), DIMENSION(max_12)  :: xa3         ! zone center after lagrangian update
REAL(KIND=double), DIMENSION(max_12)  :: dvol1       ! zone volume before the Lagrangian step

REAL(KIND=double), DIMENSION(max_12)  :: grav0       ! gravitational force before Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: grav1       ! gravitational force after Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: fict0       ! fictitious force before Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: fict1       ! fictitious force after Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: pseudo_fict ! pseudo fictitious force
REAL(KIND=double), DIMENSION(max_12)  :: nuf0        ! neutrino stress before Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: nuf1        ! neutrino stress after Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: r0          ! density before Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: p0          ! pressure before Lagrangian step

REAL(KIND=double), DIMENSION(max_12)  :: dm          ! zone mass per steradian
REAL(KIND=double), DIMENSION(max_12)  :: ei_intl     ! internal energy before the Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: e_intl      ! total energy before the Lagrangian step
REAL(KIND=double), DIMENSION(max_12)  :: e_nu_intl   ! neutrino energy before the Lagrangian step

REAL(KIND=double)                     :: p_difference
REAL(KIND=double)                     :: dtheta
REAL(KIND=double), PARAMETER          :: pq_min = 0.1d0    ! criterion for the presence of a shock
REAL(KIND=double), PARAMETER          :: rho_max = 1.0d+15 ! density above which shock finding turned off

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

ntot                  = jmax + 12
j_radial              = j_ray_dim * myid_y + ji_ray + 1
pq_y(:,ji_ray,jk_ray) = zero

!-----------------------------------------------------------------------
!  Grid position evolution
!-----------------------------------------------------------------------

umid(nmin)            = zero
umid(nmax+1)          = zero

DO n = nmin-3, nmax + 4
  dm(n)               = r(n) * dvol(n)
  dtbdm(n)            = dt/dm(n)
  xa1(n)              = xa(n)
  dvol1(n)            = dvol(n)
  xa(n)               = xa(n) + dt * umid(n) / radius
  upmid(n)            = umid(n) * pmid(n)
END DO ! n = nmin-3, nmax + 4

DO n = nmin,nmax
  IF ( v_trans_0 == 'ye'  .and.  j_radial > j_shk_radial_all_p(n-6,jk_ray) ) THEN
    umid(n)           = zero
    upmid(n)          = zero
    xa(n)             = xa1(n)
  END IF ! v_trans_0 == 'ye'  .and.  j_radial > j_shk_radial_all_p(n-6,jk_ray)
END DO ! n = nmin,nmax

DO n = nmin, nmax+1
  v_e(n-6,ji_ray, jk_ray) = umid(n)
END DO ! n = nmin, nmax+1

xa1(nmin-4)           = xa(nmin-4)
xa1(nmax+5)           = xa(nmax+5)
xa1(nmax+6)           = xa(nmax+6)

!-----------------------------------------------------------------------
!  Fasten inner boundary
!-----------------------------------------------------------------------

xa(nmin)              = xa1(nmin)

!-----------------------------------------------------------------------
!  Calculate zone (volume) center
!-----------------------------------------------------------------------

CALL zone_center( ngeom, nmin-4, nmax+5, xa1, xa2 )
CALL zone_center( ngeom, nmin-4, nmax+5, xa , xa3 )

DO n = nmin-4, nmax+5
  dx(n)               = xa(n+1) - xa(n)
END DO ! n = nmin-4, nmax+5

!-----------------------------------------------------------------------
!  Calculate forces using coordinates at t(0) and at t+dt(1), note that
!   fictitious forces at t+dt depend on updated velocity, but we ignore
!   this
!-----------------------------------------------------------------------

CALL forces_center( ji_ray, jk_ray, xa1, xa2, u, v, w, r, umid, grav0, &
& fict0, nuf0, pseudo_fict )
CALL forces_center( ji_ray, jk_ray, xa , xa3, u, v, w, r, umid, grav1, &
& fict1, nuf1, pseudo_fict )

!-----------------------------------------------------------------------
!  Calculate dvolume and average area based on geometry of sweep
!-----------------------------------------------------------------------

SELECT CASE (ngeom)

  CASE(:-1)

    WRITE (nlog,*) 'Geometry', ngeom, ' not implemented.'
    STOP

  CASE(0)

    DO n = nmin-3, nmax+4
      dvol(n)        = dx(n)
      amid(n)        = 1.d0
    END DO ! n = nmin-3, nmax+4

  CASE(1)

    DO n = nmin-3, nmax+4
      dvol(n)        = dx(n) * ( xa(n) + 0.5d0 * dx(n) )
      amid(n)        = 0.5d0 * ( xa(n) + xa1(n) )
    END DO ! n = nmin-3, nmax+4

  CASE(2)

    DO n = nmin-3,nmax+4
      dvol(n)        = dx(n) * ( xa(n) * (xa(n) + dx(n) ) + dx(n) * dx(n) * third ) 
      amid(n)        = ( xa(n) - xa1(n)) * ( third * ( xa(n) - xa1(n) ) + xa1(n) ) + xa1(n) * xa1(n)
    END DO ! n = nmin-3, nmax+4

  CASE(3)

    DO n = nmin-3, nmax+4
      dvol(n)         = dx(n) * radius
      amid(n)         = 1.d0
    END DO ! n = nmin-3, nmax+4

  CASE(4)

    DO n = nmin-3, nmax+4
      dvol(n)         = ( DCOS(xa(n)) - DCOS(xa(n+1)) ) * radius
      dtheta          = xa(n) - xa1(n)
      IF ( dtheta == 0.0d0 ) THEN
        amid(n)       = DSIN(xa(n))
      ELSE ! dtheta /= 0.0d0
        amid(n)       = ( DCOS(xa1(n)) - DCOS(xa(n)) )/dtheta
      END IF ! dtheta == 0.0d0
    END DO ! n = nmin-3, nmax+4

  CASE(5)

    DO n = nmin-3, nmax+4
      dvol(n)         = dx(n) * radius
      amid(n)         = 1.d0
    END DO ! n = nmin-3, nmax+4

  CASE(6:)

    WRITE (nlog,*) 'Geometry', ngeom, ' not implemented.'
    STOP

END SELECT

!-----------------------------------------------------------------------
!
!      \\\\\ EVOLVE THE DENSITY, VELOCITY, AND TOTAL ENERGY /////
!
!-----------------------------------------------------------------------

DO n = nmin, nmax

!-----------------------------------------------------------------------
!  Density evolution: Lagrangian mode ==> simply due to change in
!   comoving volume
!-----------------------------------------------------------------------

  r0(n)               = r(n)
  p0(n)               = p(n)
  r(n)                = r(n) * ( dvol1(n)/dvol(n) )
  r(n)                = DMAX1( r(n), smallr )

!-----------------------------------------------------------------------
!  Velocity evolution due to pressure acceleration and forces
!-----------------------------------------------------------------------

  uold(n)             = u(n)

  p_difference        = pmid(n+1) - pmid(n)
  IF ( DABS(p_difference)/pmid(n) <= 1.d-100 ) p_difference = zero
  u(n)                = u(n) - dtbdm(n) * p_difference * 0.5 * ( amid(n+1) + amid(n) ) &
&                     + 0.5d0 * dt * ( fict0(n) + nuf0(n) + fict1(n) + nuf1(n) )
  IF ( v_trans_0 == 'ye'  .and.  j_radial > j_shk_radial_all_p(n-6,jk_ray) ) u(n) = zero

!-----------------------------------------------------------------------
!  Store initial internal and total energies in order to advance them
!   adding increment. Add neutrino energies.
!-----------------------------------------------------------------------

  ei_intl(n)          = ei(n) + e_nu(n)
  e_intl(n)           = e(n)  + e_nu(n)
  e_nu_intl(n)        = e_nu(n)

!-----------------------------------------------------------------------
!  Evolve the neutrino energy and pressure
!-----------------------------------------------------------------------

  p(n)                = p(n) - p_nu(n)
  e_nu(n)             = e_nu_intl(n) * ( r(n)/r0(n) )**third
  p_nu(n)             = third * e_nu(n) * r(n)
  p(n)                = p(n) + p_nu(n)

!-----------------------------------------------------------------------
!  Evolve internal energy if no shock is present, otherwise,
!   evolve total energy.
!  Subtract the advanced neutrino energy from the evolved energy to get
!   the matter energy.
!-----------------------------------------------------------------------

  lshock(n-6)         = .false.
  CALL pseudo_y( n-6, ji_ray, jk_ray, r0(n), v_e(n-5,ji_ray,jk_ray), v_e(n-6,ji_ray,jk_ray), r(n) )
  IF ( pq_y(n-6,ji_ray,jk_ray)/( aesv(n-6,1,ji_ray,jk_ray) + epsilon ) >= pq_min  .and. &
&      r(n) < rho_max )                                                  lshock(n-6) = .true.
  IF ( lshock(n-6) ) THEN
    e(n)              = e_intl(n) - dtbdm(n) * ( amid(n+1) * upmid(n+1) - amid(n) * upmid(n) ) &
&                     + 0.5d0 * dt * ( uold(n) * ( nuf0(n) + grav0(n) ) + u(n) * ( nuf1(n) + grav1(n) ) )
    e(n)              = e(n) - e_nu(n)
  ELSE ! not lshock(n-6)
    ei(n)             = ei_intl(n) + p(n) * ( r(n) - r0(n) )/( r(n) * r0(n) )
    ei(n)             = ei(n) - e_nu(n)
  END IF ! lshock(n-6)

END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
!  Recompute the internal energy with centered pressure and update the
!   temperature
!-----------------------------------------------------------------------

DO n = nmin, nmax

  IF ( .not. lshock(n-6) ) THEN
    CALL tgvndeye_sweep_y( n, n, ji_ray, jk_ray, r, r0 )
    ei(n)             = ei_intl(n) + half * ( p(n) + p0(n) ) * ( r(n) - r0(n) )/( r(n) * r0(n) )
    ei(n)             = ei(n) - e_nu(n)
    CALL tgvndeye_sweep_y( n, n, ji_ray, jk_ray, r, r0 )
  END IF ! .not. lshock(n-6)

END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
!  Grid change requires updated parabolic coeff
!-----------------------------------------------------------------------

CALL sweepbc( nlefty, nrighty, nmin, nmax, ji_ray, jk_ray )
CALL paraset( ntot, zparay, dx, xa, nmin-3, nmax+3, ngeom )

!-----------------------------------------------------------------------
!  Extract the ei from e if a shock is present and update the
!   temperature
!-----------------------------------------------------------------------

DO n = nmin, nmax
  IF ( lshock(n-6) ) THEN
    CALL e_decompose( n, n )
    CALL tgvndeye_sweep_y( nmin, nmax, ji_ray, jk_ray, r, r0 )
  END IF ! lshock(n-6)
END DO ! n = nmin, nmax

RETURN
END SUBROUTINE evolve_y
