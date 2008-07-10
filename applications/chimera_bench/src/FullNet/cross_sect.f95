SUBROUTINE cross_sect
!===============================================================================
!  This routine calculates the cross section for each reaction.
!===============================================================================

USE kind_module, ONLY : double
USE abundances, ONLY : yt
USE conditions, ONLY : rhot, t9t
USE controls, ONLY : iweak, idiag, iscrn
USE cross_sect_data, ONLY : n1i, n2i, n3i, irev1, irev2, irev3, iwk1, iwk2, iwk3, &
& iffn, ires1, ires2, ires3, nreac, csect1, csect2, csect3, rpf1, rpf2, n1i, rc1, &
& rc2, rc3, h2
USE ffn_data, ONLY : nffn, rf
USE nuclear_data, ONLY : nname, aa, zz
USE nuc_number, ONLY : ny
USE part_funct_data, ONLY : gg

IMPLICIT NONE

INTEGER            :: j, k, mu, ieos

REAL (KIND=double) :: t09(7)
REAL (KIND=double) :: ene, ye, yps, y2e, emas, eta
REAL (KIND=double) :: funct , rat, rhot2, zz1, zz2, zz3, zz12, aa12, h01, h02

emas               = 4.d0
ieos               = 0
eta                = 1.d+30

!-----------------------------------------------------------------------
!  Call EOS, necessary for screening
!-----------------------------------------------------------------------

CALL en( ny, yt, ye, yps, y2e )

!IF ( iscrn >= 1 ) THEN
!  V=1./rhot                                                       
!  CALL state( t9t, v, y, pa, ea, dlma, ye, pta, pva, eta, eva, ieos, ny, &
!&  emas,yps,beta)
!END IF

!-----------------------------------------------------------------------
!  Calculate necessary thermodynamic moments
!-----------------------------------------------------------------------

ene                = ye * rhot                                                               
t09(1)             = 1.d0                                                                
t09(2)             = t9t**(-1) 
t09(3)             = t9t**(-1./3.) 
t09(4)             = t9t**(+1./3.) 
t09(5)             = t9t          
t09(6)             = t9t**(5./3.)
t09(7)             = log(t9t)  
!     WRITE(50,"(a3,4es12.4)") 'THR',t9t,rhot,ene,ye

!-----------------------------------------------------------------------
!  Calculate partition functions for each nucleus at t9t
!-----------------------------------------------------------------------

CALL partf( t9t ) 

!-----------------------------------------------------------------------
!  If there are any FFN reactions, calculate their rates
!-----------------------------------------------------------------------

IF ( nffn > 0 ) THEN
  CALL ffnrate(t9t,ene)   
END IF

!-----------------------------------------------------------------------
!  Calculate csects for reactions with one reactant
!-----------------------------------------------------------------------

funct              = 1.0d0
rat                = 1.0d0
WHERE ( irev1 == 1 )                   ! If it's a reverse reaction
  rpf1(:)          = gg(n1i(2,:)) * gg(n1i(3,:))/gg(n1i(1,:))
ELSEWHERE
  rpf1             =1.0d0
END WHERE

IF ( iweak > 0 ) THEN             ! All Rates on

  WHERE ( iwk1 /= 2  .and.  iwk1 /= 3 )        ! If it's not an FFN reaction
    csect1         = rpf1 * DEXP(matmul(t09,rc1))
  ELSEWHERE
    iffn           = int(rc1(1,:))
    csect1         = rat * funct * rf(iffn)
  END WHERE
  WHERE ( iwk1 == 1 ) csect1 = ene * csect1 ! If it's a non-FFN EC reaction

ELSE IF ( iweak < 0 ) THEN         ! Only weak rates used

  WHERE ( iwk1 == 0 ) csect1 = 0.0d0
  WHERE ( iwk1 == 1 ) csect1 = ene * rpf1 * DEXP(MATMUL(t09,rc1))
  WHERE ( iwk1 == 2  .or.  iwk1 == 3 ) 
    iffn           = int(rc1(1,:))
    csect1         = rat * funct * rf(iffn)
  END WHERE
  WHERE ( iwk1 >= 4 ) csect1 = rpf1 * DEXP(MATMUL(t09,rc1))    

ELSE                         ! weak interactions are off (iweak==0)

  WHERE ( iwk1 == 0 ) 
    csect1         = rpf1 * DEXP(MATMUL(t09,rc1))  
  ELSEWHERE
    csect1         = 0.0d0
  END WHERE

END IF

IF ( idiag >= 5 ) THEN
  WRITE(50,"(a,i5)") 'CSect1',nreac(1)
  WRITE(50,"(i5,5a5,3i3,es17.9)") &
&    (k,nname(n1i(1,k)),'-->',(nname(n1i(j,k)),j=2,3),' ',&
&    iwk1(k),irev1(k),ires1(k),csect1(k), k=1,nreac(1))
END IF

!-----------------------------------------------------------------------
!  Calculate screening corrections for 2 reactant reactions
!-----------------------------------------------------------------------

IF ( iscrn >= 1 ) Then
!  CALL screen2
Else
  h2               = 0.0d0
Endif

!-----------------------------------------------------------------------
!  Calculate the csect for reactions with 2 reactants
!-----------------------------------------------------------------------

WHERE ( irev2 == 1 )   ! If it's a reverse reaction mult. by ratio of part. func.
  rpf2             = gg(n2i(3,:)) * gg(n2i(4,:))/( gg(n2i(1,:)) * gg(n2i(2,:)) )
ELSEWHERE
  rpf2             = 1.0d0
END WHERE

IF( iweak >0 ) THEN        ! All rates are on

  csect2           = rhot * rpf2 * DEXP(MATMUL(t09,rc2)+h2)
  WHERE ( iwk2 == 1 ) csect2 = ene * csect2  ! If it's a non FFN EC reaction

ELSE IF ( iweak < 0 ) THEN    ! Only weak rates

  WHERE ( iwk2 == 0 ) 
    csect2         = 0.0d0
  ELSEWHERE
    csect2         = rhot * rpf2 * DEXP(MATMUL(t09,rc2)+h2)
  END WHERE
  WHERE ( iwk2 == 1 ) csect2 = ene * csect2  

ELSE                    ! Weak interactions are off  (iweak=0)

  WHERE( iwk2 == 0 )
    csect2         = rhot * rpf2 * DEXP(MATMUL(t09,rc2)+h2)           
  ELSEWHERE
    csect2         = 0.0d0
  END WHERE

END IF

IF ( idiag >= 5 ) THEN

  WRITE (50,"(a,i5)") 'CSect2',nreac(2)
  WRITE (50,"(i5,5a5,3i3,es17.9)") &
&    (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,4), &
&    iwk2(k),irev2(k),ires2(k),csect2(k),k=1,nreac(2))
END IF

IF ( idiag >= 5 ) WRITE (50,"(a,i5)") 'CSect3',nreac(3)
rhot2              = rhot**2

DO mu=1,nreac(3)

  If(iscrn >= 1) Then
    zz1              = zz(n3i(1,mu))                                 
    zz2              = zz(n3i(2,mu))                                
    zz3              = zz(n3i(3,mu))                               
    zz12             = zz1 + zz2                              
    aa12             = aa(n3i(1,mu)) + aa(n3i(2,mu))                   

    IF( zz1*zz2 == 0.0d0 ) THEN
      h01            = 0.0d0
    ELSE
!    CALL screen(zz1,zz2,aa(n3i(1,mu)),aa(n3i(2,mu)),h01,fsr,fst) 
    END IF

    IF ( zz12 * zz3 == 0.0d0 ) THEN                          
      h02            = 0.0d0
    ELSE
!     CALL screen(zz12,zz3,aa12,aa(n3i(3,mu)),h02,fsr,fst)  
    END IF
  ELSE
    h01=0.0 ; h02=0.0
  ENDIF

  IF ( iweak > 0 ) THEN

    csect3(mu)     = t09(1) * rc3(1,mu) + t09(2) * rc3(2,mu) + &
&      t09(3) * rc3(3,mu) + t09(4) * rc3(4,mu) + t09(5) * rc3(5,mu) + &
&      t09(6) * rc3(6,mu) + t09(7) * rc3(7,mu)
    csect3(mu)     = rhot2 * DEXP(csect3(mu)+h01+h02)
    IF ( iwk3(mu) == 1 ) csect3(mu) = csect3(mu) * ene        
    IF (csect3(mu) < 1.d-20 ) csect3(mu) = 0.0d0                   

  ELSE IF ( iweak < 0 ) THEN  ! Only Weak Reactions on

    IF ( iwk3(mu) == 0 ) csect3(mu) = 0.d0 

  ELSE                  ! Weak Reactions off (iweak=0)

    IF ( iwk3(mu) == 0 ) THEN
      csect3(mu)   = t09(1) * rc3(1,mu) + t09(2) * rc3(2,mu) + &
&        t09(3) * rc3(3,mu) + t09(4) * rc3(4,mu) + t09(5) * rc3(5,mu) + &
&        t09(6) * rc3(6,mu) + t09(7) * rc3(7,mu)
      csect3(mu)   = rhot2 * DEXP(csect3(mu)+h01+h02)
      IF ( csect3(mu) < 1.d-20 ) csect3(mu) = 0.0d0
    ELSE
      csect3(mu)   = 0.d0
    END IF 

  END IF

  IF ( idiag >= 5 ) WRITE (50,"(i5,5a5,3i3,es17.9)") &
&    mu,(nname(n3i(j,mu)),j=1,3),'-->','+++', &
&    iwk3(mu),irev3(mu),ires3(mu),csect3(mu)

END DO

RETURN                                  
END SUBROUTINE cross_sect                                   
