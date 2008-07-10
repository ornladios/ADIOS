SUBROUTINE t9rhofind( kstep )
!===============================================================================
!  This routine calculates t9 and rho as a function of time, either via 
!  interpolation or via from an analytic expression
!===============================================================================

USE kind_module, ONLY : double
USE conditions
USE edit_module, ONLY : nlog
USE thermo_data

LOGICAL            :: first = .true.
INTEGER            :: n, kstep
REAL (KIND=double) :: dt, rdt, dt9, drho

!     t9t=t9start
!     rhot=rhostart

!  Calculate T9 and rho by interpolation 

DO n=1,nh
  IF ( tt <= th(n) ) EXIT
END DO

IF ( n>1  .and.  n <= nh ) THEN
  rdt         = 1.0d0/( th(n) - th(n-1) )
  dt          = tt - th(n-1)
  dt9         = t9h(n) - t9h(n-1)
  drho        = rhoh(n) - rhoh(n-1)
  t9t         = dt * rdt * dt9 + t9h(n-1)
  rhot        = dt * rdt * drho + rhoh(n-1)
ELSE IF ( n == 1 ) THEN
  t9t         = t9h(1)
  rhot        = rhoh(1)
ELSE 
  t9t         = t9h(nh)
  rhot        = rhoh(nh)
  WRITE(nlog,*) 'Time beyond thermodynamic range',tt,' >',th(nh)
  IF ( first ) THEN
    WRITE(6,*) 'Time beyond thermodynamic range',tt,' >',th(nh)
    first     = .false.
  END IF ! first
END IF

!  Calculate T9 and rho by function
!     chi=1.0
!     thd=chi*446.0/sqrt(rhostart)
!     tint=(tstart-tt)/thd
!     rhot=rhostart*exp(tint)
!     t9t=t9start*exp(tint/3.)

!  Output T9 and rho
!     Write(50,"(a5,i5,3es12.4)") 'T9rho',kstep,tt,t9t,rhot

RETURN
END SUBROUTINE t9rhofind
