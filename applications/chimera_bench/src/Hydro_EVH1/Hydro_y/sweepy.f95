SUBROUTINE sweepy( ji_ray, jk_ray )
!-----------------------------------------------------------------------
!
!    File:         sweepy
!    Module:       sweepy
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To perform a hydro sweep in the y direction
!
!    Input arguments:
!  ji_ray           : x (radial) index of a specific y (angular) ray
!  jk_ray           : z (azimuthal) index of a specific y (angular) ray
!
!    Subprograms called:
!  tgvndeye_sweep_y : computes T, p, and the EOS gammas as a function of rho, ei, and ye
!  sweepbc          : inserts the boundary values in the ghost zones
!  volume           : computes the zone volumes
!  paraset          : computes the coefficients for parabolic interpolations on the advanced grid
!  e_compose        : sums the kinetic and internal energy
!  ppm              : computes the hydrodynamics via tha ppm method
!
!    Include files:
!  kind_module
!  edit_module, evh1_global, evh1_zone, evh1_sweep,
!  mgfld_remap_module, prb_cntl_module
!-----------------------------------------------------------------------

USE kind_module

USE edit_module, ONLY: nlog
USE evh1_global, ONLY: nlefty, nrighty, ngeomy, lagrangian
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, u, v, w, e, ei, ye, temp, xa0, &
& dx0, xa, dx, radius, egrav
USE evh1_zone, ONLY: zparay
USE mgfld_remap_module, ONLY: r0i
USE prb_cntl_module, ONLY: ihydro

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: ji_ray     ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)                    :: jk_ray     ! z (azimuthal) index of a specific y (angular) ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                     :: l_write = .false.

INTEGER                                :: ntot     ! parabolic array dimension
INTEGER, PARAMETER                     :: i_dum = 1
INTEGER, PARAMETER                     :: k_dum = 1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( ihydro == 1 ) THEN

  ntot                   = nmax + 6

!-----------------------------------------------------------------------
!  The EOS takes as input rho, ei and ye and returns the needed pressure 
!  and gammas, plus temperature.
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling tgvndeye_sweep_y from sweepy, &
&  nmin, nmax, ji_ray, jk_ray, r0i, r0i=',4i4, 2es11.3 )") &
&  nmin, nmax, ji_ray, jk_ray, r0i, r0i
  CALL tgvndeye_sweep_y( nmin, nmax, ji_ray, jk_ray, r0i, r0i )

!-----------------------------------------------------------------------
!  Set boundary conditions
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling sweepbc from sweepy, &
&  nlefty, nrighty, nmin, nmax, i_dum, k_dum=',6i4)") &
&  nlefty, nrighty, nmin, nmax, i_dum, k_dum
  CALL sweepbc( nlefty, nrighty, nmin, nmax, i_dum, k_dum )

!-----------------------------------------------------------------------
!  Compute volume elements
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling volume from sweepy, ngeomy=', &
&  i4)") ngeomy
  CALL volume ( ngeomy )

!-----------------------------------------------------------------------
!  For a Lagrangian run, grid change requires updated parabolic coeff
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling paraset from sweepy, &
&  ntot, nmin-2, nmax+2, ngeomy=',4i4)") &
&  ntot, nmin-2, nmax+2, ngeomy
  CALL paraset( ntot, zparay, dx, xa, nmin-4, nmax+4, ngeomy ) ! Redundant

!-----------------------------------------------------------------------
!  Compute the total energy ( ei(n) + ekin(n) )
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling e_compose from sweepy, &
&  nmin, nmax=',2i4)") &
&  nmin, nmax
  CALL e_compose( nmin, nmax )

!-----------------------------------------------------------------------
!  Perform the y-hydrodynamics
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling ppm from sweepy, &
&  ngeomy, ntot, ji_ray, jk_ray=',4i4)") &
&  ngeomy, ntot, ji_ray, jk_ray
  CALL ppm( ngeomy, ntot, zparay, ji_ray, jk_ray )
  IF ( l_write )  WRITE (nlog,"(' Returning from ppm in sweepy')")


END IF ! ihydro == 1

IF ( l_write )  WRITE (nlog,"(' Returning from sweepy')")

RETURN
END SUBROUTINE sweepy
