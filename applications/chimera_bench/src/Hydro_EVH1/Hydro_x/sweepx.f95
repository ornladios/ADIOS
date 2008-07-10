SUBROUTINE sweepx( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         sweepx
!    Module:       sweepx
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/07/06
!
!    Purpose:
!      To execute the hydrodynamics along the x-axis
!
!    Input arguments:
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  tgvndeye_sweep_x : computes T, p, and the EOS gammas as a function of rho, ei, and ye
!  sweepbc          : inserts the boundary values in the ghost zones
!  volume           : computes the volume of the zones
!  paraset          : computes the coefficients of a parabolic spline
!  ppm              : evolves the flow
!
!    Include files:
!  kind_module
!  edit_module, evh1_global, evh1_sweep, evh1_zone,
!  mgfld_remap_module, prb_cntl_module
!-----------------------------------------------------------------------

USE kind_module

USE edit_module, ONLY : nlog
USE evh1_global, ONLY: nleftx, nrightx, ngeomx, lagrangian
USE evh1_zone, ONLY: zparax
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, u, v, w, ei, ye, temp, xa0, dx0, &
& xa, dx, radius
USE mgfld_remap_module, ONLY: r0i
USE prb_cntl_module, ONLY: ihydro
     
IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER               :: l_write = .false.

INTEGER                          :: ntot          ! parabolic array dimension

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( ihydro == 1 ) THEN

  ntot                   = nmax + 6

!-----------------------------------------------------------------------
!  The EOS takes as input rho, ei and ye and returns the needed pressure 
!   and gammas, plus temperature.
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling tgvndeye_sweep_x from sweepx, &
&  nmin, nmax, ij_ray, ik_ray, r0i, r0i=',4i4, 2es11.3 )") &
&  nmin, nmax, ij_ray, ik_ray, r0i, r0i
  CALL tgvndeye_sweep_x( nmin, nmax, ij_ray, ik_ray, r0i, r0i )

!-----------------------------------------------------------------------
!  Set boundary conditions
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling sweepbc from sweepx, &
&  nleftx, nrightx, nmin, nmax, ij_ray, ik_ray=',6i4)") &
&  nleftx, nrightx, nmin, nmax, ij_ray, ik_ray
  CALL sweepbc( nleftx, nrightx, nmin, nmax, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Compute volume elements
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling volume from sweepx, ngeomx=', &
&  i4)") ngeomx
  CALL volume ( ngeomx )

!-----------------------------------------------------------------------
!  For a Lagrangian run, grid change requires updated parabolic coeff
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling paraset from sweepx, &
&  ntot, nmin-4, nmax+4, ngeomx=',4i4)") &
&  ntot, nmin-4, nmax+4, ngeomx
  CALL paraset( ntot, zparax, dx, xa, nmin-4, nmax+4, ngeomx ) ! Redundant

!-----------------------------------------------------------------------
!  Compute the total energy ( ei(n) + ekin(n) + egrav(n) )
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling e_compose_g from sweepx, &
&  nmin, nmax=',2i4)") &
&  nmin, nmax
  CALL e_compose_g( nmin, nmax )

!-----------------------------------------------------------------------
!  Evolve the flow
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling ppm from sweepx, &
&  ngeomx, ntot, ij_ray, ik_ray=',4i4)") &
&  ngeomx, ntot, ij_ray, ik_ray
  CALL ppm( ngeomx, ntot, zparax, ij_ray, ik_ray )
  IF ( l_write )  WRITE (nlog,"(' Returning from ppm in sweepx')")

END IF ! ihydro

IF ( l_write )  WRITE (nlog,"(' Returning from sweepx')")

RETURN
END SUBROUTINE sweepx
