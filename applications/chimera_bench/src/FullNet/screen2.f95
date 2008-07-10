SUBROUTINE screen2

USE kind_module, ONLY : double

USE controls, ONLY : iscrn
USE nuclear_data, ONLY : zz, aa
USE cross_sect_data, ONLY :  n2i, h2, nreac

IMPLICIT none
SAVE

INTEGER :: mu
REAL(KIND=double) :: z1, z2, z12, a1, a2, z, zl, dzldt, dzldr
REAL(KIND=double) :: gsi, gsiz, gsb, f, fsr, fst, gnp, gnz, em2, tau, &
& t9er, gt3, f90
REAL(KIND=double) :: hw, hi, hs, dhdz

COMMON /scrn/ gnp, gsi, zl, dzldr, dzldt, t9er                         

!     Write(20,"(a3,3es15.7)") 'EOS',GNP,GSI,ZL

DO mu=1,nreac(2)                                                      
  z1                = zz(n2i(1,mu))                                                  
  z2                = zz(n2i(2,mu))                                                  
  a1                = aa(n2i(1,mu))                                                  
  a2                = aa(n2i(2,mu))                                                  
  z12               = z1 * z2
  IF ( z12 == 0.0d0  .or.  iscrn == 0 ) THEN
    h2(mu)=0.0                                                        
  ELSE
  
!-----------------------------------------------------------------------
!  Weak and Intermediate Screening Graboske et al.         
!-----------------------------------------------------------------------

    z               = zl * z12                                                           
    gsiz            = gsi * z12                                                       
    gsb             = gsiz**0.86                                                    
    IF ( gsiz > 1.0d0 ) THEN
      f             = 0.38d0 * ( ( 1.d0 + gsiz )**1.86d0 - gsb * gsiz -1.d0 )/gsb          
    ELSE
      f             = 0.61943d0 * gsiz**0.07d0                              
    END IF
    hw              = z * (1.d0 + z * ( LOG(z) + 0.8364d0 ) )                                      
    hi              = f * z**0.86d0                                                      
    IF ( hw <= hi ) THEN
      h2(mu)        = hw         
      dhdz          = z * ( 1.d0 + 2.d0 * z * ( LOG(z) + 1.3364d0 ) )                
    ELSE
      h2(mu)        = hi                             
      dhdz          = 0.86d0 * hi                                   
    END IF
    fsr             = dhdz * dzldr                                                    
    fst             = dhdz * dzldt                                                    

!-----------------------------------------------------------------------
!  Strong screening by Itoh et al.(1990)                             
!-----------------------------------------------------------------------

    gnz             = gnp * z12 * 2.0d0/( z1**(1./3.) + z2**(1./3.) )
    IF (gnz>=0.4) THEN
      em2           = a1 * a2 * 2.0d0/( a1 + a2 )                
      tau           = 3.3722d0 * ( em2 * z12**2 * t9er )**(1./3.) 
      gt3           = 3.0d0 * gnz/tau
      f90           = ( .0455d0 * gt3 +.348d0 * gt3**3 + 9.49d0 * gt3**6 &
&      - .123d0 *gt3**12 + .101d0 * gt3**13 )/( 1.d0 + 100.d0 * gt3**4   &
&      + .267d0 * gt3**12 )
      hs            = 1.25d0 * gnz - tau * f90 
!      hs=gnz*(1.25-0.285*gt3)    
      IF ( hs < h2(mu) ) THEN
        h2(mu)      = hs               
        fsr         = gnz * ( 1.25d0 - .57d0 * gt3 )/3.d0 
        fst         = -gnz * ( 1.25d0 - .475d0 * gt3 )   
      END IF
    END IF
!         Write(20,"(3a5,i6,4es11.4)") 'SCR2',nname(n2i(1,mu)),
!    &         nname(n2i(2,mu)),dexp(h0(mu)),hw,hi,hs
  END IF
Enddo                                                             
!     Write(20,*) (nname(n2i(1,mu)),nname(n2i(2,mu)),h2(mu),mu=1,n)

RETURN
END SUBROUTINE screen2                                                              
 