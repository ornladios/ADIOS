SUBROUTINE etotal(flag)
!-----------------------------------------------------------------------
! Computes the total energy and writes to device 6
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : frpi

USE evh1_sweep, ONLY: nmin, nmax, ei, e, r, dvol, dvol0
USE t_cntrl_module, ONLY: time, dtnph
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

LOGICAL                                :: flag     ! flag for using updated or initial variables

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: n        ! padded zone index

REAL(KIND=double)                      :: e_total  ! total energy
REAL(KIND=double)                      :: ei_total ! total internal energy

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 3002 FORMAT (' e_total=',1pe14.7,' ei_total=',1pe14.7,' time=',1pe10.3,' dtnph=',1pe10.3)


e_total          = 0.0d0
ei_total         = 0.0d0
 
IF (flag) THEN

  DO  n = nmin,nmax
    e_total      = e_total + e(n) * r(n) * dvol(n) * frpi
    ei_total     = ei_total + ei(n) * r(n) * dvol(n) * frpi
  END DO

ELSE
      
  DO  n = nmin,nmax
    e_total      = e_total + e(n) * r(n) * dvol0(n) * frpi
    ei_total     = ei_total + ei(n) * r(n) * dvol0(n) * frpi
  END DO

END IF
            
!WRITE (6,3002) e_total,ei_total,time,dtnph

RETURN
END SUBROUTINE etotal
