MODULE evh1_global
!-----------------------------------------------------------------------
! (formerly global.h) 
! global variables accessible by (most) anything
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

INTEGER :: i_radial                 ! the unshifted radial zone corresponding to j_ray
INTEGER :: ncycle                   ! cycle number
INTEGER :: ndim                     ! # of Geometric dimensions
INTEGER :: ngeomx, ngeomy, ngeomz   ! XYZ Geometry flag
INTEGER :: nleftx, nlefty, nleftz   ! XYZ Lower Boundary Condition
INTEGER :: nrightx,nrighty,nrightz  ! XYZ Upper Boundary Condition
INTEGER :: restart                  ! =0 for a clean start
INTEGER :: nfile
INTEGER :: i_grav                   ! gravitation key

LOGICAL :: lagrangian               ! .true. is Lagrangian mode

!..Constants

REAL(KIND=double) :: time, dt, starttime
REAL(KIND=double) :: smallp, smallr, small  ! Minimum Values
REAL(KIND=double) :: svel                   ! maximum sound speed
REAL(KIND=double) :: courant                ! timestep fraction of courant limit
REAL(KIND=double) :: degen                  ! degeneracy - nondegeneracy criterion
REAL(KIND=double) :: v_diff                 ! grod-aligned shock smoothing parameter

!..Boundary values

REAL(KIND=double) :: uinflo,dinflo,vinflo,winflo,pinflo,einflo,yeinflo,gcinflo,geinflo 
REAL(KIND=double) :: uonflo,donflo,vonflo,wonflo,ponflo,eonflo,yeonflo,gconflo,geonflo 

!..Dump variables

!REAL :: gdump(42)

!..Common Block are still necessary for the dump variables to work
!      common /globe/ time,dt,starttime,courant,avgmass,boltzman,  
!    &       pi,G,bmass,smallp,smallr,small,svel,  
!    &       uinflo,dinflo,vinflo,winflo,pinflo,einflo, 
!    &       yeinflo,gcinflo,geinflo,uonflo,donflo,vonflo, 
!    &       wonflo,ponflo,eonflo,yeonflo,gconflo,geonflo, 
!    &       ndim,ncycle,ngeomx,ngeomy,ngeomz,nleftx,nlefty, 
!    &       nrightx,nrighty,nleftz,nrightz
!     equivalence (gdump(1),time)

END module evh1_global
