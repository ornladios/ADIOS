SUBROUTINE coord_bc( nleft, nright, nmin, nmax, xa, dx, xa0, dx0, ntot )
!-----------------------------------------------------------------------
!
!    File:         coord_bc
!    Module:       coord_bc
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/02/04
!
!    Purpose:
!      To compute coordinates of ghost zones.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  nleft          : left-hand boundary flags
!  nright         : right-hand boundary flags
!  nmin           : minimum paddded array index
!  nmin           : maximum paddded array index
!  xa             : padded radial coordinate after Lagr update
!  dx             : padded zone thickness after Lagr update
!  xa0            : padded radial coordinate after Eul update
!  dx0            : padded zone thickness after Eul update
!  ntot           : number of zones plus ghost zones
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module
!  evh1_sweep
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE evh1_sweep, ONLY : sweep
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                 :: nleft    ! left boundary condition key
INTEGER, INTENT(in)                                 :: nright   ! right boundary condition key
INTEGER, INTENT(in)                                 :: nmin     ! minimum paddded array index
INTEGER, INTENT(in)                                 :: nmax     ! maximum paddded array index
INTEGER, INTENT(in)                                 :: ntot     ! number of zones plus ghost zones

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(ntot+1) :: xa       ! padded coordinate after Lagr update
REAL(KIND=double), INTENT(inout), DIMENSION(ntot)   :: dx       ! padded zone thickness after Lagr update
REAL(KIND=double), INTENT(inout), DIMENSION(ntot+1) :: xa0      ! padded final coordinate
REAL(KIND=double), INTENT(inout), DIMENSION(ntot)   :: dx0      ! padded final zone thickness

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                             :: n        ! zone index
INTEGER                                             :: nmax1n   ! nmax+1-n
INTEGER                                             :: nminn1   ! nmin+n-1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : external fixed (eg, uleftbc,pleftbc,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))
!------------------------------------------------------------------------

xa0 (nmax+1) = xa0(nmax) + dx0(nmax)

!------------------------------------------------------------------------
!  Left (Inner) boundary
!------------------------------------------------------------------------

SELECT CASE (nleft)  

  CASE (0)            ! Grid symmetric with left (inner) edge

    DO n = 1, 6
      nminn1           = MIN( nmin + n - 1, nmax )
      dx(nmin-n)       = dx(nminn1)
      xa(nmin-n)       = xa(nmin-n+1) - dx(nmin-n)
      dx0(nmin-n)      = dx0(nminn1)
      xa0(nmin-n)      = xa0(nmin-n+1) - dx0(nmin-n)
    END DO

  CASE (1:2)          ! Grid spacing uniform

    DO n = 1, 6
      dx(nmin-n)       = dx(nmin)
      xa(nmin-n)       = xa(nmin-n+1) - dx(nmin-n)
      dx0(nmin-n)      = dx0(nmin+n-1)
      xa0(nmin-n)      = xa0(nmin-n+1) - dx0(nmin-n)
    END DO

  CASE (3)            ! Grid symetric with right (outer) edge 

    DO n = 1, 6
      nmax1n           = MAX( nmax + 1 - n, nmin )
      dx(nmin-n)       = dx(nmax1n)
      xa(nmin-n)       = xa(nmin-n+1) - dx(nmin-n)
      dx0(nmin-n)      = dx0(nmax1n)
      xa0(nmin-n)      = xa0(nmin-n+1) - dx0(nmin-n)
    END DO

  CASE DEFAULT

  WRITE (6,*) 'sweep ',sweep,': nleft=',nleft,' not supported in coordbc'

END SELECT

!------------------------------------------------------------------------
!  Right (Outer) Boundary
!------------------------------------------------------------------------

SELECT CASE (nright) 

  CASE (0)            ! Grid symmetric with right (outer) edge    

    DO n = 1, 6
      nmax1n           = MAX( nmax + 1 - n, nmin )
      dx(nmax+n)       = dx(nmax1n)
      xa(nmax+n)       = xa(nmax+n-1) + dx(nmax+n-1)
      dx0(nmax+n)      = dx0(nmax1n)
      xa0(nmax+n)      = xa0(nmax+n-1) + dx0(nmax+n-1)
    END DO
      xa(nmax+7)       = xa(nmax+6) + dx(nmax+6)
      xa0(nmax+7)      = xa0(nmax+6) + dx0(nmax+6)

  CASE(1,2,4,5,6)     ! Grid spacing uniform     

    DO n = 1, 6
      dx(nmax+n)       = dx(nmax)
      xa(nmax+n)       = xa(nmax+n-1) + dx(nmax+n-1)
      dx0(nmax+n)      = dx0(nmax+1-n)
      xa0(nmax+n)      = xa0(nmax+n-1) + dx0(nmax+n-1)
    END DO
      xa(nmax+7)       = xa(nmax+6) + dx(nmax+6)
      xa0(nmax+7)      = xa0(nmax+6) + dx0(nmax+6)

  CASE (3)            ! Grid symetric with left (inner) edge   

    DO n = 1, 6
      nminn1           = MIN( nmin + n - 1, nmax )
      dx(nmax+n)       = dx(nminn1)
      xa(nmax+n)       = xa(nmax+n-1) + dx(nmax+n-1)
      dx0(nmax+n)      = dx0(nminn1)
      xa0(nmax+n)      = xa0(nmax+n-1) + dx0(nmax+n-1)
    END DO
      xa(nmax+7)       = xa(nmax+6) + dx(nmax+6)
      xa0(nmax+7)      = xa0(nmax+6) + dx0(nmax+6)

  CASE DEFAULT

    WRITE (6,*) 'sweep ',sweep,': nright=',nright,' not supported in coordbc'

END SELECT

RETURN
END SUBROUTINE coord_bc
