SUBROUTINE phi_intrp( nmin, nmax, ntot, xi, xf, egrav_i, egrav_if )
!-----------------------------------------------------------------------
!
!    File:         phi_intrp
!    Module:       phi_intrp
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/09/07
!
!    Purpose:
!      To interpolate the gravitational potential, egrav_i from the grid
!       xi to the grid xf.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nmin        : minimum padded x-array index
!  nmax        : maximum padded x-array index
!  ntot        : padded x-array extent
!  xi          : padded x grid zone midpoints at which egrav_i is defined
!  xf          : padded x grid zone midpoints to which egrav_i is to be interpolated
!  egrav_i     : gravitational potential defined on xi
!  egrav_if    : gravitational potential interpolated to xf
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!        none
!
!    Output arguments (common):
!        none
!
!    Include files:
!      kind_module
!      evh1_sweep
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE evh1_sweep, ONLY: n_min=>nmin, n_max=>nmax

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nmin          ! minimum padded x-array index
INTEGER, INTENT(in)              :: nmax          ! maximum padded x-array index
INTEGER, INTENT(in)              :: ntot          ! padded x-array extent

REAL(KIND=double), INTENT(in), DIMENSION(ntot)    :: xf       ! padded x grid zone midpoints to which egrav_i is to be interpolated

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(ntot) :: egrav_if ! gravitational potential interpolated to xf

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(ntot) :: xi       ! padded x grid zone midpoints at which egrav_i is defined
REAL(KIND=double), INTENT(inout), DIMENSION(ntot) :: egrav_i  ! gravitational potential defined on xi

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: n             ! do index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Fill ghost xones
!-----------------------------------------------------------------------

xi(n_min-1)         = - xi(n_min)
xi(n_max+1)         = xi(n_max) + ( xi(n_max) - xi(n_max-1) )
egrav_i(n_min-1)    = egrav_i(n_min)
egrav_i(n_max+1)    = egrav_i(n_max) + ( egrav_i(n_max) - egrav_i(n_max-1) )


!-----------------------------------------------------------------------
!
!                   \\\\\ INTERPOLATE EGRAV_I /////
!
!-----------------------------------------------------------------------

DO n = nmin,nmax
 egrav_if(n)        = rinterp( xi(n-1), xi(n), xi(n+1), egrav_i(n-1), &
& egrav_i(n), egrav_i(n+1), xf(n) )
END DO

RETURN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS
REAL (KIND=double) FUNCTION rinterp( x1, x2, x3, y1, y2, y3, x )

REAL (KIND=double) :: x1
REAL (KIND=double) :: x2
REAL (KIND=double) :: x3
REAL (KIND=double) :: y1
REAL (KIND=double) :: y2
REAL (KIND=double) :: y3
REAL (KIND=double) :: x
REAL (KIND=double) :: a
REAL (KIND=double) :: b
REAL (KIND=double) :: c
REAL (KIND=double) :: denom

denom               = 1.d0/( ( x2 - x1 ) * ( x3 - x2 ) * ( x1 - x3 ) )

a                   = ( ( y2 - y1 ) * ( x3 - x1 ) - ( y3 - y1 ) * ( x2 - x1 ) ) * denom
b                   = ( ( y3 - y1 ) * ( x2 - x1 ) * ( x2 + x1 )           &
&                   -   ( y2 - y1 ) * ( x3 - x1 ) * ( x3 + x1 ) ) * denom
c                   = - ( y1 * x2 * x3 * ( x3 - x2 )                      &
&                   +     x1 * y2 * x3 * ( x1 - x3 )                      &
&                   +     x1 * x2 * y3 * ( x2 - x1 ) ) * denom
rinterp             = x * ( a * x + b ) + c
END FUNCTION rinterp

END SUBROUTINE phi_intrp
