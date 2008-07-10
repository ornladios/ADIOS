SUBROUTINE en(n,y,ye,yps,y2e)  
!===============================================================================
!  This routine calculates moments of the abundance distribution needed 
!  for the EOS.
!===============================================================================

USE kind_module, ONLY : double
USE nuclear_data

INTEGER :: n

REAL (KIND=double)  :: ye,yps,y2e
REAL(KIND=double), dimension(n) :: y

ye        = SUM(zz*y)
y2e       = SUM(zz**2*y)
yps       = SUM(y) + ye 
      
yps       = yps + ye
!      ye=sum(zz*y)
!      y2e=sum(zz**2*y)
!      yps=sum(y)+ye 

!     Write(50,"(a3,3es16.8)") 'EN',ye,y2e,yps

RETURN                                                                    
END SUBROUTINE en                                                                      
