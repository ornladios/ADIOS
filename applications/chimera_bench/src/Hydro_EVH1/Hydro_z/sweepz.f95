SUBROUTINE sweepz( ki_ray, kj_ray )
!-----------------------------------------------------------------------
!
!    File:         sweepz
!    Module:       sweepz
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To perform a hydro sweep in the z direction
!
!    Input arguments:
!  ki_ray       : x (radial) index of a specific z (azimuthal) ray
!  kj_ray       : y (azimuthal) index of a specific z (azimuthal) ray
!
!    Subprograms called:
!  tgvndeye_sweep_z : computes T, p, and the EOS gammas as a function of rho, ei, and ye
!  sweepbc          : inserts the boundary values in the ghost zones
!  volume           : computes the zone volumes
!  paraset          : computes the coefficients for parabolic interpolations on the advanced grid
!  e_compose        : sumes the kinetic and internal energy
!  ppm              : computes the hydrodynamics via tha ppm method
!
!    Include files:
!  kind_module
!  edit_module, evh1_global, evh1_zone, evh1_sweep,
!  mgfld_remap_module, prb_cntl_module
!-----------------------------------------------------------------------

USE kind_module

USE edit_module, ONLY: nlog
USE evh1_global, ONLY: nleftz, nrightz, ngeomz, lagrangian
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, u, v, w, e, ei, ye, temp, xa0, &
& dx0, xa, dx, radius, egrav
USE evh1_zone, ONLY: zparaz
USE mgfld_remap_module, ONLY: r0i
USE prb_cntl_module, ONLY: ihydro

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: ki_ray   ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                    :: kj_ray   ! y (azimuthal) index of a specific z (azimuthal) ray

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

  IF ( l_write )  WRITE (nlog,"(' Calling tgvndeye_sweep_z from sweepz, &
&  nmin, nmax, ki_ray, kj_ray=',4i4 )") &
&  nmin, nmax, ki_ray, kj_ray
  CALL tgvndeye_sweep_z( nmin, nmax, ki_ray, kj_ray, r0i, r0i )

!-----------------------------------------------------------------------
!  Set boundary conditions
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling sweepbc from sweepz, &
&  nleftz, nrightz, nmin, nmax, i_dum, k_dum=',6i4)") &
&  nleftz, nrightz, nmin, nmax, i_dum, k_dum
  CALL sweepbc( nleftz, nrightz, nmin, nmax, i_dum, k_dum )

!-----------------------------------------------------------------------
!  Compute volume elements
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling volume from sweepz, ngeomz=', &
&  i4)") ngeomz
  CALL volume ( ngeomz )

!-----------------------------------------------------------------------
!  For a Lagrangian run, grid change requires updated parabolic coeff
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling paraset from sweepz, &
&  ntot, nmin-2, nmax+2, ngeomz=',4i4)") &
&  ntot, nmin-2, nmax+2, ngeomz
  CALL paraset( ntot, zparaz, dx, xa, nmin-4, nmax+4, ngeomz ) ! Redundant

!-----------------------------------------------------------------------
!  Compute the total energy ( ei(n) + ekin(n) )
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling e_compose from sweepz, &
&  nmin, nmax=',2i4)") &
&  nmin, nmax
  CALL e_compose( nmin, nmax )

!-----------------------------------------------------------------------
!  Perform the y-hydrodynamics
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling ppm from sweepy, &
&  ngeomz, ntot, ki_ray, kj_ray=',4i4)") &
&  ngeomz, ntot, ki_ray, kj_ray
  CALL ppm( ngeomz, ntot, zparaz, ki_ray, kj_ray )
  IF ( l_write )  WRITE (nlog,"(' Returning from ppm in sweepz')")

END IF ! ihydro == 1

IF ( l_write )  WRITE (nlog,"(' Returning from sweepz')")

RETURN
END SUBROUTINE sweepz
