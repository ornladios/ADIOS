SUBROUTINE poisson ( mode )
!
! -------------------------------------------------------------------------
! Solves the integral form of Poisson's equation in 2D spherical coordinates. 
!
!     Potential values are calculated at zone interfaces and
!     stored in array POT(1:imax+1,1:jmax+1).
!
!     POT(i,1):      0     1     2     3  ...
!     zone:          |--1--|--2--|--3--|  ...
!
!
!     MODE = 1  ==>  Calculate weights for azimuthal integration.
!                    Must be used in first call to this subroutine 
!                    or in the case of moving grid in an angular 
!                    direction. 
!                    Else  MODE = 0
!
!     INPUT            zph(ii,js)   density (g/cc)
!                      zxa(ii+1)    radius (in cm) of left zone interface
!                      zya(jj+1)    angle (in radians) of left zone interface
!
!     OUTPUT           zgx(ii,js)   acceleration in x direction
!                      zgy(ii,js)   acceleration in y direction
!                      zph(ii,js)   zone-centered potential
!
! ---------------------------------------------------------------------------

! GLOBALS

use zone
use global
      
include 'mpif.h'

! LOCALS

INTEGER, PARAMETER :: nleg = 10  ! highest Legendre polynomial to be used 

INTEGER :: i, j, l, ip1, ip2, imax1, imax2, mode, isym, mpierr
REAL :: fpg, gfac, fl1, fl2, gravk
REAL, DIMENSION(0:nleg+1,0:jj/pe) :: pleg
REAL, DIMENSION(0:nleg,jj/pe) :: pint
REAL, DIMENSION(0:ii) :: r2, phiin, phirout
REAL, DIMENSION(1:ii,0:nleg) :: atrm, btrm
REAL, DIMENSION(ii+1,js+1) :: pot

!---------------------------------------------------------------------

gravk = 6.6720e-8
fpg   = 4. * pi * gravk
gfac  = -0.5 * fpg


IF (mode .EQ. 1) THEN  ! For the first call calculate the weights for the azimuthal integration

  ! ... first calculate the Legendre-polynomials

  do j = 0, js-1
    pleg(0,j) = 1.0
    pleg(1,j) = cos(zya(mype*js+j+1))
  enddo
  pleg(0,js) = 1.0
  pleg(1,js) = cos(zya(mype*js+js)+zdy(mype*js+js))

  DO l = 2, nleg+1
    fl1 = 1. - 1./float(l)
    fl2 = 1. + fl1
    DO j = 0, js
      pleg(l,j) = fl2 * pleg(1,j) * pleg(l-1,j) - fl1 * pleg(l-2,j)
    ENDDO
  ENDDO

  ! ... now calculate the integrals of the Legendre-polynomials:
  !         \int_{th_{i-1}}^{th_i} sin(th)*pleg(cos(th)) dth

  DO j = 1, js
    pint(0,j) = pleg(1,j-1) - pleg(1,j)
  ENDDO

  DO l = 1 , nleg
    fl1 = 1.  / (2.*float(l) + 1.)
    DO j = 1, js
      pint(l,j) = ( pleg(l-1,j) - pleg(l-1,j-1) - pleg(l+1,j) + pleg(l+1,j-1) ) * fl1
    ENDDO
  ENDDO

ENDIF

! ... Now build inner and outer sum.

DO j = 1, js+1
 DO i = 1, imax+1
   pot(i,j) = 0.
 ENDDO
ENDDO

r2(0) = 0.
DO i = 1, imax
  r2(i) = (zxa(i)+zdx(i))**2 
ENDDO

phiin(0)= 0.
phirout(imax) = 0.
DO l = 0, nleg
 DO i = 1, imax
   btrm(i,l) = 0. 
   DO j = 1, js
     btrm(i,l) = btrm(i,l) + zph(i,j) * pint(l,j)
   ENDDO
 ENDDO
ENDDO

! add up all atrm's from all the processors

call MPI_ALLREDUCE(btrm, atrm, imax*(nleg+1), MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)

! --- loop for phiin:
DO l = 0, nleg

 DO i = 1, imax
   fl1 = ( zxa(i) / (zxa(i)+zdx(i)) )**float(l+1)
   phiin(i) = phiin(i-1) * fl1 + atrm(i,l) * ( r2(i) - r2(i-1)*fl1 ) / float(l+3)
 ENDDO 

 ! --- loop for phirout:
 IF (l .NE. 2) THEN
   DO i = imax-1, 1, -1
     ip1 = i + 1           
     fl1 = ( zxa(ip1) / (zxa(ip1)+zdx(ip1)) )**float(l)
     phirout(i) = phirout(ip1)* fl1 + atrm(ip1,l) * ( (zxa(ip1)+zdx(ip1))**2 * fl1 - zxa(ip1)*zxa(ip1) ) / float(2-l)
   ENDDO
 ELSE           
   DO i = imax-1, 1, -1             
     ip1 = i + 1
     phirout(i) = phirout(ip1)*(zxa(ip1)/(zxa(ip1)+zdx(ip1)))**float(l)+atrm(ip1,l)*zxa(ip1)**float(l)*log((zxa(ip1)+zdx(ip1))/zxa(ip1) )
   ENDDO
 ENDIF

 IF (l .EQ. 0) THEN
   phirout(0) = phirout(1) + 0.5*atrm(1,l)*zxa(2)**2
 ELSE
   phirout(0) = 0
 ENDIF

 ! ... now add all contributions together

 DO j = 0, js
   fl1 =  gfac * pleg(l,j)
   DO i = 0, imax
     pot(i+1,j+1) = pot(i+1,j+1) + fl1*(phiin(i)+phirout(i))
   ENDDO
 ENDDO

ENDDO

! ... compute differences in order to store accelerations

do j = 1, js
 do i = 1, imax
   zgx(i,j) = ((pot(i,j)+pot(i,j+1))-(pot(i+1,j)+pot(i+1,j+1)))/(2.0*zdx(i))
   zgy(i,j) = ((pot(i,j)+pot(i+1,j))-(pot(i,j+1)+pot(i+1,j+1)))/(2.0*zdy(j)*zxc(i))
   zph(i,j) = 0.25*(pot(i,j) + pot(i+1,j) + pot(i,j+1) + pot(i+1,j+1))
 enddo
enddo

RETURN
END


!________________________________________________________________________
!
! Ewald Mueller                        Internet: ewald@MPA-Garching.MPG.de  
! Max-Planck-Institut fuer Astrophysik                              
! Karl-Schwarzschild-Str. 1            phone :   +49 89 3299 3209 
! D-85748 Garching                     FAX:      +49 89 3299 3235
! Germany                              http://www.MPA-Garching.MPG.de
!_________________________________________________________________________

