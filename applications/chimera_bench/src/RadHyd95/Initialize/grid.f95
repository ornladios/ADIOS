SUBROUTINE grid( nmin, nmax, dmin, dmax, xag, xcg, dxg, zoom, n_dim )
!-----------------------------------------------------------------------
!
!    File:         grid
!    Module:       grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/13/04
!
!    Purpose:
!      Create grid to cover physical size from dmin to dmax.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input-output arguments:
!
!  nmin  : inner zone of the grid to be covered
!  nmax  : outer zone of the grid to be covered
!  dmin  : minimum value of grid coordinate
!  dmax  : maximum value of grid coordinate
!  zoom  : increase in the grid spacing of successive zones\
!  n_dim : array dimension
!
!    Output arguments:
!
!  xag   : grid edge values
!  xcg   : grid canter values
!  dxg   : grid widths
!
!    Output arguments (common):
!        none
!
!    Include files:
!      kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: nmin     ! minimum coordinate index
INTEGER, INTENT(in)                    :: nmax     ! maximum coordinate index
INTEGER, INTENT(in)                    :: n_dim    ! array dimension

REAL(KIND=double), INTENT(in)          :: dmin     ! minimum coordinate value
REAL(KIND=double), INTENT(in)          :: dmax     ! maximum coordinate value
REAL(KIND=double), INTENT(in)          :: zoom     ! grid zoom factor

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(n_dim+1) :: xag      ! grid edges
REAL(KIND=double), INTENT(out), DIMENSION(n_dim)   :: xcg      ! grid centers
REAL(KIND=double), INTENT(out), DIMENSION(n_dim)   :: dxg      ! grid widths

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: n        ! zone index
INTEGER                                :: nzones   ! number of zones

!........Set up the grid................................................

nzones               = nmax - nmin + 1
IF ( zoom == 1.0d0 ) THEN
  dxg(1)             = ( dmax - dmin )/float(nzones)
ELSE
  dxg(1) = ( dmax - dmin ) * ( zoom - 1.0d0 )/( zoom**nzones - 1 )
END IF

xag(1)               = dmin
xcg(1)               = xag(1) + 0.5d0 * dxg(1) 
DO n = 2, nzones
  xag(n)             = xag(n-1) + dxg(n-1)
  dxg(n)             = dxg(n-1) * zoom
  xcg(n)             = xag(n) + 0.5d0 * dxg(n)
END DO
xag(nzones+1)        = dmax

RETURN
END SUBROUTINE grid
