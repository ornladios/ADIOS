!-----------------------------------------------------------------------
!    Module:       incrmnt_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE incrmnt_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Temperature increments
!-----------------------------------------------------------------------
!  dtmpmn(j,i,ij_ray,ik_ray)  : Temperature change  of radial zone j
!   during the current cycle due to process i.
!
!  dtmpmn(j,1,ij_ray,ik_ray)  : temperature change of radial zone j due
!   to hydro x-sweep.
!
!  dtmpmn(j,2,ij_ray,ik_ray)  : temperature change of radial zone j due
!   to hydro y-sweep.
!
!  dtmpmn(j,3,ij_ray,ik_ray)  : temperature change of radial zone j due
!   to hydro z-sweep.
!
!  dtmpmn(j,4,ij_ray,ik_ray)  : temperature change of radial zone j due
!   to nuclear reactions.
!
!  dtmpmn(j,9,ij_ray,ik_ray)  : temperature change of radial zone j due
!   to ML convection.
!
!  dtmpmn(j,10,ij_ray,ik_ray) : temperature change of radial zone j due
!   to artificial launch of an explosion.
!
!  dtmpnn(j,i,ij_ray,ik_ray): Temperature change  of radial zone j during
!   the current cycle due neutrino transport process i.
!
!  dtmpnn(j,1,ij_ray,ik_ray): temperature change of radial zone j due to
!   source and transport.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: dtmpmn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: dtmpnn

!-----------------------------------------------------------------------
!  Internal energy increments
!-----------------------------------------------------------------------
!  denergy_int(j) : Internal energy cahnge of radial zone j due to hydro
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: denergy_int

!-----------------------------------------------------------------------
!  Electron fraction increments
!-----------------------------------------------------------------------
!  dye(j,n,ij_ray,ik_ray): Electron fraction change  of radial zone j
!   during the current cycle due to  n-neutrinos.
!
!  dyecnvt(j,ij_ray,ik_ray): change in ye in radial zone j due to ML
!    convection.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: dye
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: dyecnvt

END module incrmnt_module
