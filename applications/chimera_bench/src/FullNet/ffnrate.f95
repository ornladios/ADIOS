SUBROUTINE ffnrate(t9,ene) 
!-----------------------------------------------------------------------------  
!  This routine calculates the reaction rates for the weak rates of 
!  Fuller, Fowler, & Newman
!-----------------------------------------------------------------------------  

USE kind_module, ONLY : double
USE ffn_data

INTEGER i1(4)
INTEGER i,le1,lt1

REAL (KIND=double) :: t9,ene,tg(13),egl(11),enl,dt,de,ddt

DATA tg/0.01,0.1,0.2,0.4,0.7,1.0,1.5,2.,3.,5.,10.,30.,100./               
DATA egl/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11./                              

DO i=1,13
  IF ( t9 <= tg(i) ) EXIT
END DO

lt1         = i - 1
enl         = LOG10(ene)                                                           
le1         = INT(enl)
dt          = tg(lt1+1) - tg(lt1)                                                      
de          = egl(le1+1) - egl(le1)                                                    
ddt         = t9 - tg(lt1)                                                            
i1(1)       = 13 * (le1-1) + lt1                                                      
i1(2)       = i1(1) + 1                                                             
i1(3)       = 13 * le1 + lt1                                                          
i1(4)       = i1(3) + 1                                                             

DO i=1,nffn                                                              
  r1(i)     = ffnsum(i,i1(1)) + (ffnsum(i,i1(2)) - ffnsum(i,i1(1)))/dt*ddt
  r2(i)     = ffnsum(i,i1(3)) + (ffnsum(i,i1(4)) - ffnsum(i,i1(3)))/dt*ddt
  rf(i)     = r1(i) + ( r2(i) - r1(i) )/de * ( enl - egl(le1) )    
  IF ( rf(i) < -30.d0 ) THEN
    rf(i)   = 0.d0
  ELSE
    rf(i)   = 10.d0**rf(i)     
  END IF
END DO                                                                     

RETURN                                                                    
END SUBROUTINE ffnrate
