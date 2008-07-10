!-----------------------------------------------------------------------
!    Module:       scratch_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE scratch_module

USE kind_module

SAVE


!........Cube indices...................................................

!      idrse(j), itrse(j), and iyrse(j) are integers defining the location of log(rho), log(t), and ye
!       for radial zone j on the grid for neutrino-electron scattering functions, i.e.,

!           idrse(j)/dgrid < log(rho(j)) < ( idrse(j) + 1 )/dgrid 

!           itrse(j)/tgrid <  log(t(j))  < ( itrse(j) + 1 )/tgrid

!          0.5 - iyrse(j)/ygrid < ye < 0.5 - ( iyrse(j) + 1 )/ygrid

!      The eight grid points surrounding log(rho), log(t), and ye for radial zone j are referred to 
!       as the unit cube j. Neutrino-electron scattering rates for radial zone j are stored at 
!       the corners of unit cube j. Rates for zone j are interpolated from the rates stored at the
!       corners.

INTEGER, ALLOCATABLE, DIMENSION(:)                      :: idrse
INTEGER, ALLOCATABLE, DIMENSION(:)                      :: itrse
INTEGER, ALLOCATABLE, DIMENSION(:)                      :: iyrse


!........Neutrino-electron scattering functions at cube corners.........

!      scte0i(j,k,kp,id,it,iy), scte0ii(j,k,kp,id,it,iy): arrays of the zero moments of the 
!       neutrino-electron scattering functions as a function of the radial zone (j), incident
!       neutrino energy (k), final neutrino energy (kp) at the unit cube corners id, it, and iy 
!       (id, it, iy = 1,2). Interpolations of neutrino-electron scattering functions are computed
!       from this table.

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:)  :: scte0i
REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:)  :: scte0ii


!........Interpolated neutrino-electron scattering functions............

!      scef(j,k,n) : Neutrino-electron scattering function i

!                i = 1: a0w, i = 2: b0w, i = 3: c0w,
!                i = 4: a1w, i = 5: b1w, i = 6: c1w.

!      scefd(i,j,k,n)     = d(scef(i,j,k,n))/d(density)
!      sceft(i,j,k,n)     = d(scef(i,j,k,n))/d(temperature)
!      scefy(i,j,k,n)     = d(scef(i,j,k,n))/d(electron fraction)
!      scefp0(i,kp,j,k,n) = d(scef(i,j,k,n))/d(psi0(j,kp,n))
!      scefp1(i,kp,j,k,n) = d(scef(i,j,k,n))/d(psi1(j,kp,n))

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: scef
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: scefd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: sceft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: scefy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)    :: scefp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)    :: scefp1

END module scratch_module
