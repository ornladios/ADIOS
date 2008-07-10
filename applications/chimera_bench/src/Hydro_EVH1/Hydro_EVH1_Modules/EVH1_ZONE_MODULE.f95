module evh1_zone 
!=======================================================================
! (formerly zone.h) global (3D) data arrays
!
! 09jun90 gbl
! 21jun90 jmb delete zeps
! jul 90 jfh 2D version
! jan 92 jmb 3D version
! oct 98 jmb mpi version
! oct 00 wrh composition added
! jan 01 wrh module-ized 
! feb 03 krd adding energy variable in account for gravity
!=======================================================================

USE kind_module
USE numerical_module, ONLY : zero
SAVE

INTEGER, PARAMETER                :: nmf =2   ! Number of constituents

!........Dimensions: kk and jj must be evenly divisible by pe...........

INTEGER, PARAMETER                :: pe=1                ! Designed number of PEs
INTEGER                           :: imax                ! Physical dimesions
INTEGER                           :: jmax                ! Physical dimesions
INTEGER                           :: kmax                ! Physical dimesions
INTEGER                           :: is,js               ! domain slice dimensions
INTEGER                           :: nmfx=nmf

!........State variables................................................

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zro     ! density: zone average (g cm^{-3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zei     ! internal energy: zone average (ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zegrav  ! grav. pot. energy: zone average (ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zeg     ! tot. energy: zone average (ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zux     ! zone average velocity x direction (cm s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zuy     ! zone average velocity y direction (cm s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zuz     ! zone average velocity z direction (cm s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zte     ! temperature (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: zye     ! electron fraction
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: znu_str ! neutrino stress (dynes/g)

!........Coordinates....................................................

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zxa     ! x grid zone edge location
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zya     ! y grid zone edge location
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zza     ! z grid zone edge location
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zdx     ! xa(i+1) - xa(i)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zdy     ! ya(i+1) - ya(i)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zdz     ! za(i+1) - za(i)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zxc     ! x grid zone midpoint locations
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zyc     ! y grid zone midpoint locations
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: zzc     ! z grid zone midpoint locations

!........PPM intepolation factors, for REAL(KIND=double) zones and......
!........ghosts.........................................................

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)       :: zparax
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)       :: zparay
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)       :: zparaz

END module evh1_zone
