SUBROUTINE ppm( ngeom, nz, para, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         ppm
!    Module:       ppm
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/07/06
!
!    Purpose:
!       Using the 1D arrays of rho, u, and P, perform a 1D lagrangian
!        hydrodynamics evolution:
!        - obtain parabolic interpolations of rho, u, P
!        - compute input states from these interpolations for Riemann problem
!        - call the Riemann solver to find the time averages umid, pmid
!        - evolve the finite difference equations to get updated
!           values of rho, u, and E
!
!    Input arguments:
!  ngeom      : geometry flag
!  nz         : parabolic array dimension
!  para       : parabolic coefficients obtained from paraset
!  ij_ray     : y-index of an x ray, x-index of an y ray, x index of a z array
!  ik_ray     : z-index of an x ray, z-index of an y ray, y index of a z array
!
!    Output arguments:
!  flatten    : sets the variable flat between 0.0 (smooth flow) and 1.0 (strong shock)
!  parabola   : calculates the interpolation parabolas themselves
!  states     : computes the left and right states (integrated over the charachteristics) for each variable
!  riemann    : solves the Riemann shock tube problem
!  evolve_i   : updates the densities, velocities, and energies along direction i
!
!    Subprograms called:
!  flatten    : determins the location of a strong shock
!  parabola   : computes piecewise parabolic fitting coefficients for quantities
!  states     : sets up the left and right states at zone interfaces for the Riemann problem
!  riemann    : solves the Riemann problem at zone interfaces
!  evolve_i   : updates the positions, veloities and energies along i-rays.
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, evh1_sweep
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: max_12
USE numerical_module, ONLY: zero

USE edit_module, ONLY : nlog
USE evh1_sweep, ONLY : nmin, nmax, gc, ge, p, r, u, radius, sweep, &
& flat_s, l_shock, xa, dx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                  :: ngeom    ! geometry index
INTEGER, INTENT(in)                  :: nz       ! parabolic array dimension
INTEGER, INTENT(in)                  :: ij_ray   ! y-index of an x ray, x-index of an y ray, x index of a z array
INTEGER, INTENT(in)                  :: ik_ray   ! z-index of an x ray, z-index of an y ray, y index of a z array

REAL(KIND=double), INTENT(in), DIMENSION(10,nz) :: para ! parabolic coefficients obtained from paraset

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                   :: l_write = .false.

INTEGER                              :: n        ! padded zone index

REAL(KIND=double), DIMENSION(max_12) :: rl       ! density at left edge of zone
REAL(KIND=double), DIMENSION(max_12) :: r6       ! density parabola coeffecient
REAL(KIND=double), DIMENSION(max_12) :: dr       ! density slope
REAL(KIND=double), DIMENSION(max_12) :: ul       ! velocity at left edge of zone
REAL(KIND=double), DIMENSION(max_12) :: u6       ! velocity parabola coeffecient
REAL(KIND=double), DIMENSION(max_12) :: du       ! velocity slope
REAL(KIND=double), DIMENSION(max_12) :: pl       ! pressure at left edge of zone
REAL(KIND=double), DIMENSION(max_12) :: p6       ! pressure parabola coeffecient
REAL(KIND=double), DIMENSION(max_12) :: dp       ! pressure slope
REAL(KIND=double), DIMENSION(max_12) :: gel      ! gamma_p at left edge of zone
REAL(KIND=double), DIMENSION(max_12) :: ge6      ! gamma_p parabola coeffecient
REAL(KIND=double), DIMENSION(max_12) :: dge      ! gamma_p slope
REAL(KIND=double), DIMENSION(max_12) :: gcl      ! gamma_1 at left edge of zone
REAL(KIND=double), DIMENSION(max_12) :: gc6      ! gamma_1 parabola coeffecient
REAL(KIND=double), DIMENSION(max_12) :: dgc      ! gamma_1 slope

REAL(KIND=double), DIMENSION(max_12) :: rlft     ! average density in the + wave state (left of interface)
REAL(KIND=double), DIMENSION(max_12) :: rrgh     ! average density in the - wave state (right of interface)
REAL(KIND=double), DIMENSION(max_12) :: ulft     ! average velocity in the + wave state (left of interface)
REAL(KIND=double), DIMENSION(max_12) :: urgh     ! average velocity in the - wave state (right of interface)
REAL(KIND=double), DIMENSION(max_12) :: plft     ! average pressure in the + wave state (left of interface)
REAL(KIND=double), DIMENSION(max_12) :: prgh     ! average pressure in the - wave state (right of interface)
REAL(KIND=double), DIMENSION(max_12) :: gelft    ! average gamma_p in the + wave state (left of interface)
REAL(KIND=double), DIMENSION(max_12) :: gergh    ! average gamma_p in the - wave state (right of interface)
REAL(KIND=double), DIMENSION(max_12) :: gclft    ! average gamma_1 in the + wave state (left of interface)
REAL(KIND=double), DIMENSION(max_12) :: gcrgh    ! average gamma_1 in the - wave state (right of interface)

REAL(KIND=double), DIMENSION(max_12) :: flat     ! flattening coefficients
REAL(KIND=double), DIMENSION(max_12) :: umid     ! time averaged velocity at the zone interface
REAL(KIND=double), DIMENSION(max_12) :: pmid     ! time averaged pressure at the zone interface

!-----------------------------------------------------------------------
!  Calculate flattening coefficients for smoothing near shocks
!-----------------------------------------------------------------------

CALL flatten( flat )
flat_s              = zero
flat_s(nmin:nmax)   = flat(nmin:nmax)

!-----------------------------------------------------------------------
!  Set logical variable to indicate shock presence
!-----------------------------------------------------------------------

l_shock             = .false.
DO n = nmin,nmax
  IF ( flat(n) > zero   ) l_shock(n) = .true.
END DO

!-----------------------------------------------------------------------
!  Interpolate parabolae for fluid variables
!-----------------------------------------------------------------------

CALL parabola( nmin-4, nmax+4, nz, para, p , dp , p6 , pl , flat, 1, 0, &
& ngeom )
CALL parabola( nmin-4, nmax+4, nz, para, r , dr , r6 , rl , flat, 1, 0, &
& ngeom )
CALL parabola( nmin-4, nmax+4, nz, para, u , du , u6 , ul , flat, 1, 0, &
& ngeom )
CALL parabola( nmin-4, nmax+4, nz, para, ge, dge, ge6, gel, flat, 1, 0, &
& ngeom ) !EVH1
CALL parabola( nmin-4, nmax+4, nz, para, gc, dgc, gc6, gcl, flat, 1, 0, &
& ngeom ) !EVH1

!-----------------------------------------------------------------------
!  Integrate parabolae to get input states for Riemann problem
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,*) ' Calling states from ppm'
CALL states( ngeom, pl, ul, rl, gel, gcl, p6, u6, r6, ge6, gc6, dp, du, &
& dr, dge, dgc, plft, ulft, rlft, gelft, gclft, prgh, urgh, rrgh, gergh, &
& gcrgh, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Call the Riemann solver to obtain the zone face averages, umid and pmid
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,*) ' Calling riemann from ppm'
CALL riemann( nmin-3, nmax+4, ge, gc, prgh, urgh, rrgh, gergh, gcrgh, &
& plft, ulft, rlft, gelft, gclft, pmid, umid, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Do lagrangian update using umid and pmid
!-----------------------------------------------------------------------

IF ( sweep == 'x' ) THEN
IF ( l_write ) WRITE (nlog,*) ' Calling evolve_x from ppm'
  CALL evolve_x( ngeom, umid, pmid, ij_ray, ik_ray )
ELSE IF ( sweep == 'y' ) THEN
IF ( l_write ) WRITE (nlog,*) ' Calling evolve_y from ppm'
  CALL evolve_y( ngeom, umid, pmid, ij_ray, ik_ray )
ELSE IF ( sweep == 'z' ) THEN
IF ( l_write ) WRITE (nlog,*) ' Calling evolve_z from ppm'
  CALL evolve_z( ngeom, umid, pmid, ij_ray, ik_ray )
END IF

RETURN
END SUBROUTINE ppm
