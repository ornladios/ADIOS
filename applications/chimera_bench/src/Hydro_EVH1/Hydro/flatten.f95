SUBROUTINE flatten( flat )
!-----------------------------------------------------------------------
!  Flatten looks for signs of strong shocks and sets the variable flat
!   between 0.0 (smooth flow) and 1.0 (strong shock).
!  Simplified method of C&W: eqn. A.1 and A.2.
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: max_12

USE evh1_global, ONLY : small
USE evh1_sweep, ONLY : nmin, nmax, p, u

IMPLICIT none
SAVE

INTEGER                              :: n        ! zone index

REAL(KIND=double), DIMENSION(max_12) :: flat     ! density at left edge of zone
REAL(KIND=double), DIMENSION(max_12) :: steep    ! density at left edge of zone

!REAL(KIND=double), PARAMETER         :: omega1 = 0.75d0  ! parameter as suggested in C&W
REAL(KIND=double), PARAMETER         :: omega1 = 0.5d0  ! parameter as suggested in C&W
!REAL(KIND=double), PARAMETER         :: omega2 = 5.0d0   ! parameter as suggested in C&W
REAL(KIND=double), PARAMETER         :: omega2 = 10.0d0   ! parameter as suggested in C&W
!REAL(KIND=double), PARAMETER         :: epsilon = 0.33d0 ! parameter as suggested in C&W
REAL(KIND=double), PARAMETER         :: epsilon = 1.0d0 ! parameter as suggested in C&W
REAL(KIND=double)                    :: delp1            ! 
REAL(KIND=double)                    :: delp2            ! 
REAL(KIND=double)                    :: shock            ! 
REAL(KIND=double)                    :: temp1            ! 
REAL(KIND=double)                    :: temp2            ! 

!-----------------------------------------------------------------------
!  Look for presence of a shock using pressure gradient and sign of
!   velocity jump:  shock = 1 if there is a shock in the zone, else shock = 0
!  Compute steepness parameter based on steepness of pressure jump if 
!   there is a shock.
!-----------------------------------------------------------------------

DO n = nmin-4, nmax+4
  delp1              = p(n+1) - p(n-1)
  delp2              = p(n+2) - p(n-2)
  IF ( DABS(delp2) < small ) delp2 = small
  shock              = DABS(delp1)/DMIN1( p(n+1), p(n-1) ) - epsilon
  shock              = DMAX1( 0.0d0, shock )
  IF ( shock > 0.0d0) shock = 1.0d0
  IF ( u(n-1) < u(n+1) ) shock = 0.0d0
  temp1              = ( delp1/delp2 - omega1 ) * omega2
  steep(n)           = shock * DMAX1( 0.d0, temp1 )
END DO

!-----------------------------------------------------------------------
!  Set phony boundary conditions for the steepness parameter
!-----------------------------------------------------------------------

steep(nmin-5)        = steep(nmin-4)
steep(nmax+5)        = steep(nmax+4)

!-----------------------------------------------------------------------
!  Set flatening coefficient based on the steepness in nieghboring
!   zones
!-----------------------------------------------------------------------

DO n = nmin-4, nmax+4
  temp2              = DMAX1( steep(n-1), steep(n), steep(n+1) )
  flat(n)            = DMAX1( 0.0d0, DMIN1( 0.5d0, temp2 ) )
END DO

RETURN
END SUBROUTINE flatten


