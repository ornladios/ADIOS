SUBROUTINE poisson( mode, imin, imax, nx, ij_ray_dim, ik_ray_dim, ny, x_e, &
& x_c, dx_c, y_e, dy_c, rho_c, grav_x_c, grav_y_c, grav_pot_c, grav_x_e, &
& grav_y_e, grav_pot_e )
!
!-----------------------------------------------------------------------
!  Solves the integral form of Poisson's equation in 2D spherical coordinates. 
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
!     INPUT            rho_c(nx,ij_ray_dim,ik_ray_dim)    density (g cm^{-3})
!                      x_e(nx+1)              radius (in cm) of left zone interface
!                      y_e(ny+1)              angle (in radians) of left zone interface
!
!     OUTPUT           grav_x_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in x direction
!                      grav_y_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in y direction
!                      grav_pot_c(ii,ij_ray_dim,ik_ray_dim) zone-centered potential
!                      grav_x_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in x direction
!                      grav_y_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in y direction
!                      grav_pot_e(ii,ij_ray_dim,ik_ray_dim) zone-edged potential
!
!    Input arguments:
!  imin             : lower x-array index
!  imax             : upper x-array index
!  ny               : x-array extent
!  ij_ray_dim       : number of y-zones on a processor before swapping
!  ik_ray_dim       : number of z-zones on a processor before swapping
!  ny               : y-array extent
!  x_e              : x grid zone faces
!  x_c              ! x grid zone centers
!  dx_c             : x_e(i+1) - x_e(i)
!  y_e              : y grid zone left interfaces
!  dy_c             : y_e(i+1) - y_e(i)
!  rho_c            : density (g cm^{-3})
!
!    Output arguments:
!  grav_x_c         : zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_y_c         : zone-centered y-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_pot_c       : zone-centered zone-centered gravitational potential (erg g^{-1})
!  grav_x_e         : zone-edged x-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_y_e         : zone-edgedy-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_pot_e       : zone-edgedzone-centered gravitational potential (erg g^{-1})
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY : zero, half, frpi
USE physcnst_module, ONLY : g

USE parallel_module, ONLY : myid, ierr

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                 :: imin       ! minimum x-array index
INTEGER, INTENT(in)                                 :: imax       ! maximum x-array index
INTEGER, INTENT(in)                                 :: nx         ! x-array extent
INTEGER, INTENT(in)                                 :: ny         ! y-array extent
INTEGER, INTENT(in)                                 :: ij_ray_dim ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                                 :: ik_ray_dim ! number of z-zones on a processor before swapping

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)      :: x_e        ! x grid zone left interfaces (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: x_c        ! x grid zone centers (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: dx_c       ! x_e(i+1) - x_e(i) (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)      :: y_e        ! y grid zone left interfaces
REAL(KIND=double), INTENT(in), DIMENSION(ny)        :: dy_c       ! y_e(j+1) - y_e(j)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: rho_c      ! density (g cm^{-3})

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: grav_x_c   ! zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: grav_y_c   ! zone-centered y-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: grav_pot_c ! zone-centered zone-centered gravitational potential (erg g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: grav_x_e   ! zone-edged x-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim+1,ik_ray_dim) :: grav_y_e   ! zone-edged y-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: grav_pot_e ! zone-edged zone-centered gravitational potential (erg g^{-1})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER                                  :: nleg = 10         ! highest Legendre polynomial to be used 

INTEGER                                             :: mode              ! initialize flag
INTEGER                                             :: i                 ! radial index
INTEGER                                             :: ip1               ! i+1
INTEGER                                             :: j                 ! angular index
INTEGER                                             :: l                 ! order of Legendre polynomial

REAL(KIND=double)                                   :: fpg, gfac, fl1, fl2
REAL(KIND=double), DIMENSION(0:10+1,0:10), SAVE     :: pleg
REAL(KIND=double), DIMENSION(0:10,10), SAVE         :: pint
REAL(KIND=double), DIMENSION(0:imax)                :: r2, phiin, phirout
REAL(KIND=double), DIMENSION(1:imax,0:nleg)         :: atrm, btrm
REAL(KIND=double), DIMENSION(nx+1,ij_ray_dim+1)     :: pot

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' ik_ray_dim=',i4,' > 1 cannot expand grav in a Legendre series')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Stop if ik_ray_dim > 1
!-----------------------------------------------------------------------

IF ( ik_ray_dim > 1 ) THEN
  WRITE (nlog,1001) ik_ray_dim
  STOP
END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

fpg                    = frpi * g
gfac                   = -0.5d0 * fpg

!-----------------------------------------------------------------------
!  For the first call calculate the weights for the azimuthal integration
!-----------------------------------------------------------------------

IF ( mode == 1 ) THEN

!-----------------------------------------------------------------------
!  First calculate the Legendre-polynomials
!  Use the recurrence relation
!
!     (2l+1)yP (y) = (l+1)P   (y) + lP   (y)
!             l            l+1        l-1
!
! ==> P   (y) = (2l+1)yP (y)/(l+1) - lP   (y)/(l+1)
!      l+1              l              l-1
!
! ==> P (y) = (2l-1)yP   (y)/l - (l-1)P   (y)/l
!      l              l-1              l-2
!-----------------------------------------------------------------------

  DO j = 0, ij_ray_dim-1
    pleg(0,j)          = 1.0d0
    pleg(1,j)          = DCOS(y_e(myid*ij_ray_dim+j+1))
  END DO
  pleg(0,ij_ray_dim)   = 1.0d0
  pleg(1,ij_ray_dim)   = DCOS(y_e(myid*ij_ray_dim+ij_ray_dim+1))

  DO l = 2, nleg+1
    fl1                = 1.d0 - 1.d0/DBLE(l)
    fl2                = 1.d0 + fl1
    DO j = 0, ij_ray_dim
      pleg(l,j)        = fl2 * pleg(1,j) * pleg(l-1,j) - fl1 * pleg(l-2,j)
    END DO
  END DO

!-----------------------------------------------------------------------
!  Now calculate the integrals of the Legendre-polynomials:
!     \int_{th_{i-1}}^{th_i} sin(th)*pleg(cos(th)) dth
!
!     P'   (y) - P'   (y) = (2l+1)P (y)
!       l+1        l-1             l
!
!       j
!     /
! ==> |  P (y) dy = [ P   (j) - P   (j-1) ]/(2l+1) - [ P   (j) - P   (j-1) ]/(2l+1)
!     /   l            l-1       l-1                    l+1       l+1
!    j-1
!
!  = pint(l,j)
!-----------------------------------------------------------------------

 DO j = 1,ij_ray_dim
    pint(0,j)          = pleg(1,j-1) - pleg(1,j)
  ENDDO

  DO l = 1 , nleg
    fl1                = 1.d0/( 2.d0 * DBLE(l) + 1.d0 )
    DO j = 1,ij_ray_dim
      pint(l,j)        = ( pleg(l-1,j) - pleg(l-1,j-1) - pleg(l+1,j) + pleg(l+1,j-1) ) * fl1
    ENDDO
  ENDDO

RETURN
ENDIF ! mode == 1


!-----------------------------------------------------------------------
!                                        __
!              /                         \
!  atrm(i,l) = | dx pleg(l,x) rho(i,x) = /_   pint(l,j) rho(i,j)
!             /                           j
!
!  btrm(i,l) is the contribution to atrm(i,l) from eqch processor
!-----------------------------------------------------------------------

pot                    = zero

r2(0)                  = zero
r2(imin:imax)          = x_e(imin+1:imax+1)**2 

phiin(0)               = zero
phirout(imax)          = zero
DO l = 0,nleg
  DO i = 1,imax
    btrm(i,l)          = zero
    DO j = 1,ij_ray_dim
      btrm(i,l)        = btrm(i,l) + rho_c(i,j,k) * pint(l,j)
    ENDDO
  ENDDO
ENDDO

!-----------------------------------------------------------------------
!  Add up the btrm's from all the processors to get the atrm's
!-----------------------------------------------------------------------

atrm(imin:imax,0:nleg) = btrm(imin:imax,0:nleg)

!-----------------------------------------------------------------------
!  Compute phiin
!-----------------------------------------------------------------------

DO l = 0, nleg

  DO i = 1, imax
    fl1                = ( x_e(i) / x_e(i+1) )**float(l+1)
    phiin(i)           = phiin(i-1) * fl1 + atrm(i,l) * ( r2(i) - r2(i-1)*fl1 ) / float(l+3)
  END DO 

!-----------------------------------------------------------------------
!  Compute phirout
!-----------------------------------------------------------------------

  IF ( l /= 2 ) THEN
    DO i = imax-1,1,-1
      ip1              = i + 1           
      fl1              = ( x_e(ip1)/x_e(i+2) )**float(l)
      phirout(i)       = phirout(ip1) * fl1 + atrm(ip1,l) * ( x_e(i+2)**2 &
&                      * fl1 - x_e(ip1) * x_e(ip1) ) / ( 2.d0 - DBLE(l) )
    ENDDO
  ELSE           
    DO i = imax-1,1,-1             
      ip1              = i + 1
      phirout(i)       = phirout(ip1) * ( x_e(ip1)/x_e(i+2) )**float(l)   &
&                      + atrm(ip1,l) * x_e(ip1)**float(l) * DLOG( x_e(i+2)/x_e(ip1) )
    ENDDO
  ENDIF

  IF ( l == 0) THEN
    phirout(0)         = phirout(1) + half * atrm(1,l) * x_e(2)**2
  ELSE
    phirout(0)         = zero
  ENDIF

!-----------------------------------------------------------------------
!  Now add all contributions together
!-----------------------------------------------------------------------

  DO j = 0,ij_ray_dim
    fl1                =  gfac * pleg(l,j)
    DO i = 0, imax
      pot(i+1,j+1)     = pot(i+1,j+1) + fl1 * ( phiin(i) + phirout(i) )
    END DO ! i
  END DO ! j

END DO ! l

!-----------------------------------------------------------------------
!  Compute differences in order to store accelerations
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO i = 1,imax
    grav_x_c(i,j,1)    = ( ( pot(i,j) + pot(i,j+1) ) - ( pot(i+1,j) + pot(i+1,j+1) ) )/( 2.d0 * dx_c(i) )
    grav_y_c(i,j,1)    = ( ( pot(i,j) + pot(i+1,j) ) - ( pot(i,j+1) + pot(i+1,j+1) ) )/( 2.d0 * dy_c(j) * x_c(i) )
    grav_pot_c(i,j,1)  = 0.25d0 * ( pot(i,j) + pot(i+1,j) + pot(i,j+1) + pot(i+1,j+1) )
  END DO ! i
END DO ! j

grav_x_e  (1,:,1)      = zero
grav_y_e  (1,:,1)      = zero

DO j = 1,ij_ray_dim
  DO i = 2,imax
    grav_x_e(i,j,1)    = ( ( pot(i-1,j) + pot(i-1,j+1) ) - ( pot(i+1,j) + pot(i+1,j+1) ) )/( 2.d0 * dx_c(i) )
    grav_y_e(i,j,1)    = ( ( pot(i,j-1) + pot(i+1,j-1) ) - ( pot(i,j+1) + pot(i+1,j+1) ) )/( 2.d0 * dy_c(j) * x_c(i) )
    grav_pot_e(i,j,1)  = 0.5d0 * ( pot(i,j) + pot(i,j+1) )
  END DO ! i
END DO ! j

grav_x_e  (imax+1,:,1) = grav_x_e  (imax,:,1)
grav_y_e  (imax+1,:,1) = grav_y_e  (imax,:,1)
grav_pot_e(imax+1,:,1) = grav_pot_e(imax,:,1)

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

