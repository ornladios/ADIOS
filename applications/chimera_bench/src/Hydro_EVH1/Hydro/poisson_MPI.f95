SUBROUTINE poisson( mode, imin, imax, nx, ij_ray_dim, ik_ray_dim, ny, nz, x_e, 	&
	x_c, dx_c, y_e, y_c, dy_c, z_e, dz_c, rho_c, grav_x_c, grav_y_c, grav_z_c, &
	grav_pot_c, grav_x_e, grav_y_e, grav_z_e, grav_pot_e )

!-----------------------------------------------------------------------
!
!     File:         poisson_MPI_95
!     Type:         Subprogram
!     Author:       S. W. Bruenn, P. Marronetti
!		    Dept of Physics, FAU,
!                   Boca Raton, FL 33431-0991
!
!     Date:         1/22/08
!
!     Solves the integral form of Poisson's equation in 3D spherical coordinates. 
!
!     Potential values are calculated at zone interfaces and
!     stored in array POT(1:imax+1,1:jmax+1).
!
!     POT(i,1):      0     1     2     3  ...
!     zone:          |--1--|--2--|--3--|  ...
!
!
!     MODE = 1  ==>  Calculate weights for solid angle integration.
!                    Must be used in first call to this subroutine 
!                    or in the case of moving grid in an angular 
!                    direction. 
!                    Else  MODE = 0
!
!     INPUT            rho_c(nx,ij_ray_dim,ik_ray_dim)    density (g cm^{-3})
!                      x_e(nx+1)              radius (in cm) of left zone interface
!                      y_e(ny+1)              angle (in radians) of left zone interface
!                      z_e(nz+1)              angle (in radians) of left zone interface
!
!     OUTPUT           grav_x_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in x direction
!                      grav_y_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in y direction
!                      grav_z_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in z direction
!                      grav_pot_c(ii,ij_ray_dim,ik_ray_dim) zone-centered potential
!                      grav_x_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in x direction
!                      grav_y_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in y direction
!                      grav_z_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in z direction
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
!  y_c              ! y grid zone centers
!  dy_c             : y_e(i+1) - y_e(i)
!  z_e              : z grid zone left interfaces
!  dz_c             : z_e(i+1) - z_e(i)
!  rho_c            : density (g cm^{-3})
!
!    Output arguments:
!  grav_x_c         : zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_y_c         : zone-centered y-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_z_c         : zone-centered z-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_pot_c       : zone-centered zone-centered gravitational potential (erg g^{-1})
!  grav_x_e         : zone-edged x-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_y_e         : zone-edged y-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_z_e         : zone-edged z-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_pot_e       : zone-edgedzone-centered gravitational potential (erg g^{-1})
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY : zero, half, twpi, frpi
USE physcnst_module, ONLY : g, pi
USE edit_module, ONLY: nprint, nlog
USE parallel_module, ONLY : myid, myid_y, myid_z, ierr

USE mpi

IMPLICIT none

!include "mpif.h"
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                 :: imin       ! minimum x-array index
INTEGER, INTENT(in)                                 :: imax       ! maximum x-array index
INTEGER, INTENT(in)                                 :: nx         ! x-array extent
INTEGER, INTENT(in)                                 :: ny         ! y-array extent
INTEGER, INTENT(in)                                 :: nz         ! y-array extent
INTEGER, INTENT(in)                                 :: ij_ray_dim ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                                 :: ik_ray_dim ! number of z-zones on a processor before swapping

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)      :: x_e        ! x grid zone left interfaces (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: x_c        ! x grid zone centers (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: dx_c       ! x_e(i+1) - x_e(i) (cm)

REAL(KIND=double), INTENT(in), DIMENSION(ny+1)      :: y_e        ! y grid zone left interfaces
REAL(KIND=double), INTENT(in), DIMENSION(ny)        :: y_c        ! y grid zone centers
REAL(KIND=double), INTENT(in), DIMENSION(ny)        :: dy_c       ! y_e(j+1) - y_e(j)
REAL(KIND=double), INTENT(in), DIMENSION(nz+1)      :: z_e        ! z grid zone left interfaces
REAL(KIND=double), INTENT(in), DIMENSION(nz)        :: dz_c       ! z_e(k+1) - z_e(k)

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: rho_c      ! density (g cm^{-3})

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: grav_x_c   ! zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: grav_y_c   ! zone-centered y-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: grav_z_c   ! zone-centered z-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: grav_pot_c ! zone-centered zone-centered gravitational potential (erg g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: grav_x_e   ! zone-edged x-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: grav_y_e   ! zone-edged y-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: grav_z_e   ! zone-edged z-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: grav_pot_e ! zone-edged zone-centered gravitational potential (erg g^{-1})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER                          :: nleg_max = 12     ! highest Legendre polynomial degree to be used. 
                                                                 ! Keep this number under 2700 
INTEGER, PARAMETER                          :: ij_ray_max = 128, ik_ray_max = 128
INTEGER                                     :: mode              ! initialize flag
INTEGER                                     :: i                 ! radial index
INTEGER                                     :: ip1               ! i+1
INTEGER                                     :: j,k,jj,kk         ! angular indices
INTEGER                                     :: l,m          ! degree and order of Legendre polynomials

REAL(KIND=double)                                                    :: fl1, fl2, cos_m, sin_m
REAL(KIND=double), DIMENSION(0:nleg_max+1,0:nleg_max+1,0:ij_ray_max+1), SAVE :: pleg
REAL(KIND=double), DIMENSION(0:nleg_max+1,0:nleg_max+1), SAVE        :: fm    ! factor in potential series
REAL(KIND=double), DIMENSION(0:nleg_max,0:nleg_max,ij_ray_max), SAVE :: pint
REAL(KIND=double), DIMENSION(0:nleg_max,ik_ray_max), SAVE            :: cint, sint
REAL(KIND=double), DIMENSION(0:imax)                                 :: r2, phiin, phirout
COMPLEX(KIND=double), DIMENSION(1:imax,0:nleg_max,0:nleg_max)        :: atm, btm
REAL(KIND=double), DIMENSION(nx+1,0:ij_ray_dim+1,0:ik_ray_dim+1)     :: pot

!-----------------------------------------------------------------------
!        Local variables used by ucar subroutines
!-----------------------------------------------------------------------

INTEGER                                                              :: nrst, nrfn
INTEGER          , DIMENSION(ij_ray_dim+2)                           :: z0

REAL(KIND=double)                                                    :: sqrt2
REAL(KIND=double), DIMENSION(ij_ray_dim+1)                           :: thet1, thet2        
REAL(KIND=double), DIMENSION(nleg_max+1,ij_ray_dim+2)                :: pmm         
REAL(KIND=double), DIMENSION(nleg_max+1,ij_ray_dim+1)                :: pimm
REAL(KIND=double), DIMENSION(nleg_max+2,ij_ray_dim+1)                :: pnt,p1,p2

CHARACTER * 3 :: tail

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' ik_ray_dim=',i4,' > 1; cant expand grav in a Legendre series')

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
!  For the first call calculate the weights for the azimuthal integration
!-----------------------------------------------------------------------

IF ( mode == 1 ) THEN

!-----------------------------------------------------------------------
!
!  We calculate the fully normalized Associated Legendre Functions (ALFs)
!  and their integrals using subroutines developed by NGA (in 'ucar.f')
!  They are based on the algorithms from Paul (1978) and Gerstl (1980) 
!  IMPORTANT: These subroutines are valid for double precision and for 
!  Legendre pol. degrees of up to 2700.
!  More information can be found at (alf_sr_v121305):
!  http://earth-info.nga.mil/GandG/wgs84/gravitymod/new_egm/new_egm.html
!
!  pleg(l,m,j) := \tilde{P}^m_l(cos(theta_j))
!
!                  j
!                 /   
!  pint(l,m,j) := | \tilde{P}^m_l (y) dy 
!                 /     
!                 j-1
!
!-----------------------------------------------------------------------

  sqrt2 = DSQRT(2.d0)
  fl1   = 1.d0
  nrst  = 1
  nrfn  = ij_ray_dim+1
  
  DO j = 1, ij_ray_dim+1
    jj = myid_y*ij_ray_dim+j
    
    if (jj == 1) then
      thet2(j) = y_e(1)
    else
      thet2(j) = y_e(jj-1)
    endif
    
    thet1(j) = y_e(jj)
    
  END DO

! This subroutine calculates the sectorial (l=m) Assoc. Legendre
! functions which are used by the next sub.

  CALL alfisct(nrst,nrfn,thet1,thet2,nleg_max,fl1,pmm,pimm,  &
               z0,nleg_max)
  
  DO m = 0, nleg_max

!   This subr. calculates the ALFs. If you try to compare the 
!   output from this subroutine and published expressions for
!   P^m_l, keep in mind two things: 1) these are the fully-
!   normalized ALFs, 2) the pol. are defined module the "Condon
!   -Shortley" phase (-1)^m

    CALL alfiord(nrst,nrfn,thet1,thet2,m,nleg_max,pmm,pimm,  &
                 pnt,p1,p2,nleg_max)

    DO l = m, nleg_max

!   'fm' is the factor that enters in the final
!   assembly of the potential.

      if (m == 0) then
        fm(l,m) = - g / (2.d0 * DBLE(l) + 1.d0)
      else
        fm(l,m) = - g * sqrt2 / (2.d0 * DBLE(l) + 1.d0)
      endif
 
      DO j = 0, ij_ray_dim+1

        jj = myid_y*ij_ray_dim+j

        if (j == 0) then
          pleg(l,m,j) = p2(l+1,1)
        else
          pleg(l,m,j) = p1(l+1,j)
        endif
      
        if (jj == 0) pleg(l,m,j) = p1(l+1,j+2)
        if (jj == 1) pleg(l,m,j) = p2(l+1,1)

      END DO  ! j
      
      DO j = 1, ij_ray_dim
        pint(l,m,j) = pnt(l+1,j+1)
      END DO  ! j

    END DO  ! l
    
  END DO  ! m

!-----------------------------------------------------------------------
!  Now calculate the integrals of the azimuthal part:
!     A := \int_{phi_{k-1}}^{phi_k} exp(-I*m*phi) dphi
!
!  Dividing this into its real and imaginary parts, we get
!
!  Re(A) = \int_{phi_{k}}^{phi_{k+1}} cos(m*phi) dphi
!        = 1/m * \int_{k}^{k+1} du1
!        = sint(m,k)
!
!  where u1 := sin(m*phi). Similary
!
!  Im(A) = - \int_{phi_{k}}^{phi_{k+1}} sin(m*phi) dphi
!        = 1/m * \int_{k}^{k+1} du2
!        = cint(m,k)
!
!  where u2 := cos(m*phi). When m=0, we define
!	cint = 0 and sint = \delta \phi = \phi_{k+1} - \phi_{k}
!
!-----------------------------------------------------------------------

  cint(0,:)        = 0.d0
  
  DO k = 1, ik_ray_dim
    kk = myid_z*ik_ray_dim+k
    sint(0,k)        = dz_c(kk)
  ENDDO
  
  DO m = 1 , nleg_max
    fl1 = 1.d0 / DBLE(m) / sqrt2
    DO k = 1, ik_ray_dim
      kk = myid_z*ik_ray_dim+k
      cint(m,k)        = ( DCOS(m*z_e(kk+1)) - DCOS(m*z_e(kk)) ) * fl1
      sint(m,k)        = ( DSIN(m*z_e(kk+1)) - DSIN(m*z_e(kk)) ) * fl1
    ENDDO
  ENDDO
  
ENDIF ! mode == 1

!----------------------------------------------------------------------------------------
!                                                 __ __
!              /    /                             \  \
!  atm(i,l,m) =| dz | dy pleg(l,m,y) rho(i,y,z) = /_ /_  exp(-I m phi) pint(l,m,j) rho(i,j,k)
!              /    /                              k  j
!
!   The btm(i,l,m) are the contribution to atm(i,l,m) from each processor
!	
!----------------------------------------------------------------------------------------

DO l = 0,nleg_max
  DO m = 0,l
  
    DO i = 1, imax
    
      btm(i,l,m) = zero
      
      DO k = 1, ik_ray_dim
        DO j = 1, ij_ray_dim
          btm(i,l,m) = btm(i,l,m) +  &
                       CMPLX(sint(m,k) * pint(l,m,j) * rho_c(i,j,k), &
                             cint(m,k) * pint(l,m,j) * rho_c(i,j,k), double )
        ENDDO
      ENDDO
      
    ENDDO
    
  ENDDO
ENDDO

!-----------------------------------------------------------------------
!  Add up the btm's from all the processors to get the atm's
!-----------------------------------------------------------------------

CALL MPI_ALLREDUCE( btm, atm, imax*(nleg_max+1)**2, MPI_DOUBLE_COMPLEX, &
                    MPI_SUM, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Compute potential
!-----------------------------------------------------------------------

r2(0)                  = zero
r2(imin:imax)          = x_e(imin+1:imax+1)**2 

pot                    = zero
phiin(0)               = zero
phirout(imax)          = zero

DO l = 0, nleg_max
  DO m = 0, l

    DO k = 0, ik_ray_dim+1
    
      kk = myid_z*ik_ray_dim + k
      
      if (kk == 0) then
        cos_m  = DCOS(-m*z_e(2))
        sin_m  = DSIN(-m*z_e(2))
      else
        cos_m  = DCOS(m*z_e(kk))
        sin_m  = DSIN(m*z_e(kk))
      endif
      
!-----------------------------------------------------------------------
!  Compute phiin
!-----------------------------------------------------------------------

      DO i = 1, imax

        fl2       = cos_m * real(atm(i,l,m),double) - sin_m * aimag(atm(i,l,m))

        fl1       = ( x_e(i) / x_e(i+1) )**float(l+1)
        phiin(i)  = phiin(i-1) * fl1 + fl2 * ( r2(i) - r2(i-1)*fl1 ) / &
                    DBLE(l+3)

      END DO 

!-----------------------------------------------------------------------
!  Compute phirout
!-----------------------------------------------------------------------

      IF ( l /= 2 ) THEN
        DO i = imax-1,1,-1
          ip1              = i + 1           
          fl2              = cos_m * real(atm(ip1,l,m),double) - sin_m * aimag(atm(ip1,l,m))
 
          fl1              = ( x_e(ip1)/x_e(i+2) )**float(l)
          phirout(i)       = phirout(ip1) * fl1 + fl2            &
                       * ( x_e(i+2)**2 * fl1 - x_e(ip1) * x_e(ip1) )    &
                       / ( 2.d0 - DBLE(l) )
        ENDDO
      ELSE           
        DO i = imax-1,1,-1             
          ip1              = i + 1
          fl2              = cos_m * real(atm(ip1,l,m),double) - sin_m * aimag(atm(ip1,l,m))
  
          phirout(i)       = phirout(ip1) * ( x_e(ip1)/x_e(i+2) )**float(l) &
                       + fl2 * x_e(ip1)**float(l) * DLOG( x_e(i+2)/x_e(ip1) )
        ENDDO
      ENDIF

      IF ( l == 0) THEN
        fl2                = cos_m * real(atm(1,l,0)) - sin_m * aimag(atm(1,l,0))
        phirout(0)         = phirout(1) + half * fl2 * x_e(2)**2
      ELSE
        phirout(0)         = zero
      ENDIF

!-----------------------------------------------------------------------
!  Now add all contributions together
!-----------------------------------------------------------------------

      DO j = 0,ij_ray_dim+1
          fl1        = pleg(l,m,j) * fm(l,m) 
        DO i = 1, imax+1
          pot(i,j,k) = pot(i,j,k) + fl1 * ( phiin(i-1) + phirout(i-1) )
        END DO ! i
      END DO ! j

    END DO ! k
  
  END DO ! m
END DO ! l


!-----------------------------------------------------------------------
!  Compute differences in order to store accelerations
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
 
  fl1 = DSIN(y_c(myid_y*ij_ray_dim + j))

    DO i = 1,imax
  
      grav_x_c(i,j,k)  = ( ( pot(i,j,k)   + pot(i,j+1,k)   + pot(i,j,k+1)   + pot(i,j+1,k+1))       &
    			 - ( pot(i+1,j,k) + pot(i+1,j+1,k) + pot(i+1,j,k+1) + pot(i+1,j+1,k+1)) )	&
			 / ( 4.d0 * dx_c(i) )
			 
      grav_y_c(i,j,k)  = ( ( pot(i,j,k)   + pot(i+1,j,k)   + pot(i,j,k+1)   + pot(i+1,j,k+1) )  	&
    			 - ( pot(i,j+1,k) + pot(i+1,j+1,k) + pot(i,j+1,k+1) + pot(i+1,j+1,k+1)) ) 	&
			 /( 4.d0 * dy_c(j) * x_c(i) )
			 
      grav_z_c(i,j,k)  = ( ( pot(i,j,k)   + pot(i+1,j,k)   + pot(i,j+1,k)   + pot(i+1,j+1,k) )  	&
    			 - ( pot(i,j,k+1) + pot(i,j+1,k+1) + pot(i+1,j,k+1) + pot(i+1,j+1,k+1)) ) 	&
			 /( 4.d0 * dz_c(j) * x_c(i) * fl1 )
    
      grav_pot_c(i,j,k)= 0.125d0 * ( pot(i,j,k) + pot(i+1,j,k) + pot(i,j+1,k) + pot(i+1,j+1,k)  	&
    			+ pot(i,j,k+1) + pot(i+1,j,k+1) + pot(i,j+1,k+1) + pot(i+1,j+1,k+1) )
			
    END DO ! i = 1,imax
  END DO ! j = 1,ij_ray_dim
END DO ! k = 1,ik_ray_dim

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim

    DO i = 2,imax
      grav_x_e(i,j,k)  = ( ( pot(i-1,j,k) + pot(i-1,j+1,k) + pot(i-1,j,k+1) + pot(i-1,j+1,k+1)) 	&
    			-  ( pot(i+1,j,k) + pot(i+1,j+1,k) + pot(i+1,j,k+1) + pot(i+1,j+1,k+1)) ) 	&
			/( 4.d0 * (dx_c(i-1) + dx_c(i)) )

    END DO ! i = 2,imax
  END DO ! j = 1,ij_ray_dim
END DO ! k = 1,ik_ray_dim

grav_x_e  (1,:,:)      = zero
grav_x_e  (imax+1,:,:) = grav_x_e  (imax,:,:)
grav_y_e  = zero 
grav_z_e  = zero
grav_pot_e(imax+1,:,:) = grav_pot_e(imax,:,:)

RETURN
END
