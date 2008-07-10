      subroutine yderiv
!-----------------------------------------------------------------------------  
!  This routine calculates time derivatives for each nuclear species.
!  This calculation is performed by looping over nuclei, and summing the 
!  reaction rates for each reaction which effects that nucleus.
!-----------------------------------------------------------------------------  

USE kind_module, ONLY : double

USE abundances, ONLY : yt, ydot
USE controls, ONLY : idiag
USE cross_sect_data, ONLY : n1i, n2i, n3i, csect1, csect2, csect3
USE nuclear_data, ONLY : nname
USE nuc_number, ONLY : ny
USE reac_rate_data, ONLY : mu1, mu2, mu3, n11, n21, n22, n31, n32, n33, la, le, &
& a1, a2, a3, b1, b2, b3

IMPLICIT none
SAVE

INTEGER           :: i, j, i0, i1, la1, le1, la2, le2, la3, le3
REAL(KIND=double) :: s1, s11, s2, s22, s3, s33

!-----------------------------------------------------------------------------  
!  From the cross sections and the counting array, calculate the reaction 
!  rates
!-----------------------------------------------------------------------------  

b1              = a1 * csect1(mu1)
b2              = a2 * csect2(mu2)
b3              = a3 * csect3(mu3)

!-----------------------------------------------------------------------------  
!  Calculate Ydot for each nucleus, summing over the reactions which affect it.
!-----------------------------------------------------------------------------  

DO i0=1,ny                 

  IF ( idiag >= 5 ) WRITE (50,"(a3,a6,i4)") 'NUC',nname(i0),i0

!-----------------------------------------------------------------------------  
!  Sum over the reactions with 1 reactant
!-----------------------------------------------------------------------------  

  la1           = la(1,i0)                          
  le1           = le(1,i0)                         
  s1            = 0.0d0
  DO i1=la1,le1                    
    s11         = b1(i1) * yt(n11(i1))
    s1          = s1 + s11
    IF ( idiag >= 5 ) WRITE (50,"(3a5,'  1  ',4es15.7)")  &
&      (nname(n1i(j,mu1(i1))),j=1,3),b1(i1),yt(n11(i1)),  &
&      s11,a1(i1)
  END DO       
  IF ( idiag >= 5 ) WRITE (50,*) '1->',nname(i0),la1,le1,s1

!-----------------------------------------------------------------------------  
!  Sum over the reactions with 2 reactants
!-----------------------------------------------------------------------------  

  la2           = la(2,i0)  
  le2           = le(2,i0) 
  s2            = 0.0d0      
  DO i1=la2,le2           
    s22         = b2(i1) * yt(n21(i1)) * yt(n22(i1))
    s2          = s2 + s22                                       
    IF ( idiag >= 5 ) WRITE (50,"(4a5,4es15.7)")  &
&      (nname(n2i(i,mu2(i1))),i=1,4),b2(i1),yt(n21(i1)),  &
&      yt(n22(i1)),s22
  END DO                  
  IF ( idiag >= 5 ) WRITE (50,*) '2->',nname(i0),la2,le2,s2

!-----------------------------------------------------------------------------  
!  Sum over the reactions with 3 reactants
!-----------------------------------------------------------------------------  

  la3           = la(3,i0)                          
  le3           = le(3,i0)                         
  s3            = 0.0d0                              
  DO i1=la3,le3                  
    s33         = b3(i1) * yt(n31(i1)) * yt(n32(i1)) * yt(n33(i1))
    s3          = s3 + s33
    IF ( idiag >= 5 ) WRITE (50,"(3a5,'  3  ',5es12.4)")  &
&      (nname(n3i(i,mu3(i1))),i=1,3),b3(i1),yt(n31(i1)), &
&      yt(n32(i1)),yt(n33(i1)),s33
  END DO                                       
  IF ( idiag >= 5 ) WRITE (50,*) '3->',nname(i0),la3,le3,s3

!-----------------------------------------------------------------------------  
!  Sum the 3 components of Ydot
!-----------------------------------------------------------------------------  

  ydot(i0)      = s1 + s2 + s3                         
  IF ( idiag >= 5 ) WRITE (50,"(a4,a5,2es24.16)")  &
&    'YDOT',nname(i0),yt(i0),ydot(i0)

END DO                                 

RETURN


END SUBROUTINE yderiv                                   
