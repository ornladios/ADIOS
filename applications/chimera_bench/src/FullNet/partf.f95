SUBROUTINE partf(t9) 
!===============================================================================
!  This routine calculates the nuclear partition functions as a function of T9
!===============================================================================

USE kind_module, ONLY : double

USE nuc_number, ONLY : ny
USE part_funct_data, ONLY : t9i, g, gg

IMPLICIT none
SAVE

INTEGER           :: i, ii
REAL(KIND=double) :: t9, rdt9

DO i=1,24                                                              
  IF ( t9 <= t9i(i) ) EXIT
END DO                                                                     
ii            = i                                                                      
SELECT CASE (ii)
  CASE(2:24)
    rdt9      = ( t9 - t9i(ii-1) )/( t9i(ii) - t9i(ii-1) )
    gg(1:ny)  = DEXP( rdt9 * LOG(g(ii,1:ny)) + ( 1 - rdt9 ) * LOG(g(ii-1,1:ny)) )
  CASE (1)
    gg(1:ny)  = g(1,1:ny) 
  CASE (25)
    gg(1:ny)  = g(24,1:ny)  
 END SELECT                                                               

!  gg(0) is a placeholder for non-nuclei, gamma-rays, etc. 
gg(0)         = 1.0d0
!     Write(50,"(a5,i3,es14.7)") 'PartF',ii,t9
!     Write(50,"(5(i4,es12.4))") (i, gg(i), i=1,ny) 

RETURN                                                                 
END SUBROUTINE partf                                                                      
