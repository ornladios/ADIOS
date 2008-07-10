SUBROUTINE norm(yy)
!---------------------------------------------------------------------------- 
!  This routine renormalizes the abundances to guarantee mass conservation.
!---------------------------------------------------------------------------- 

USE kind_module, ONLY : double

USE nuclear_data, ONLY : aa
USE nuc_number, ONLY : ny

IMPLICIT none
SAVE

REAL(KIND=double) :: xtot,rxt
REAL(KIND=double), dimension(ny) :: yy

xtot         = SUM(yy*aa)   
rxt          = 1.0d0/xtot      
yy           = yy * rxt  

RETURN                                              
END SUBROUTINE norm                                                
