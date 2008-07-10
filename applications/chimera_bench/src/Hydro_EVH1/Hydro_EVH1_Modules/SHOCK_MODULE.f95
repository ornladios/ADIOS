!-----------------------------------------------------------------------
!    Module:       shock_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE shock_module

USE kind_module, ONLY : double

SAVE



!-----------------------------------------------------------------------
!  Pseudoviscous shock parameters
!-----------------------------------------------------------------------
!  ipq : pseudoviscosity switch
!
!     ipq = 0 : pseudoviscosity (pq_x(j)) computed normally.
!     ipq = 1 : pq_x(j)=0 unless u(jp) gt 0 for some jp, i.e., pq_x(j)=0 during infall.
!     ipq = 2 : pq_x(j)=0.
!
!  pq_x(j)  : pseudoviscous pressure in radial zone j.
!
!  pqr_x(j) : pseudoviscous pressure in radial zone j at the preceding 'hydro' timestep.
!
!  pqy_x(j) : pseudoviscous pressure in radial zone j due to angular shock.
!
!  pqz_x(j) : pseudoviscous pressure in radial zone j due to azimuthal shock.
!
!  q0_x(j)  : pseudoviscous pressure multiplyer for radial zone j.
!
!  pqcrit   : minimum value of q/p for presence of shock to be inferred.
!
!  pq_y(j)  : pseudoviscous pressure in angular zone j.
!
!  pqr_y(j) : pseudoviscous pressure in angular zone j at the preceding 'hydro' timestep.
!
!  q0_y(j)  : pseudoviscous pressure multiplier for angular zone j.
!
!  pq_z(j)  : pseudoviscous pressure in azimuthal zone k.
!
!  pqr_z(j) : pseudoviscous pressure in azimuthal zone k at the preceding 'hydro' timestep.
!
!  q0_z(j)  : pseudoviscous pressure multiplier for azimuthal zone k.
!
!  lshock   : shock is present in radial zone j if lshock = .true.
!-----------------------------------------------------------------------

LOGICAL, ALLOCATABLE, DIMENSION(:)               :: lshock

INTEGER                                          :: ipq

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pq_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pqr_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pqy_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pqz_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: q0_x
REAL(KIND=double)                                :: pqcrit

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pq_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pqr_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: q0_y

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pq_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: pqr_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: q0_z

!-----------------------------------------------------------------------
!  Radial location
!-----------------------------------------------------------------------
!  j_shk_radial_p     : outermost radial zone defining the shock
!   structure for the radial zones assigned to a processor.
!
!  j_shk_radial_all_p : outermost radial zone defining the shock
!   structure for all radial zones.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:)             :: j_shk_radial_p
INTEGER, ALLOCATABLE, DIMENSION(:,:)             :: j_shk_radial_all_p

!-----------------------------------------------------------------------
!  Shock boundaries
!-----------------------------------------------------------------------
!  jshockmn : Innermost radial zone defining the shock structure. Used in
!   the rezoning subroutines.
!
!  jshockmx : Outermost radial zone defining the shock structure. Used in
!   the rezoning subroutines.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                           :: jshockmn
INTEGER, DIMENSION(20)                           :: jshockmx

!-----------------------------------------------------------------------
!      nqpmax(j) :
!      qpmax(j)  :
!      qrpmax(j) :
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:)               :: nqpmax

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: qpmax
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: qrpmax


END module shock_module
