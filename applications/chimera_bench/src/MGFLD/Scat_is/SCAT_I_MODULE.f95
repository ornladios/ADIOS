!-----------------------------------------------------------------------
!    Module:       scat_i_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE scat_i_module

USE kind_module

SAVE


!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!
!  idrsi(j,ij_ray,ik_ray), itrsi(j,ij_ray,ik_ray), and iyrsi(j,ij_ray,ik_ray)
!   are integers defining the location of log(rho), log(t), and ye for
!   radial zone j on the grid  for isoenergetic scattering inverse mean
!   free paths, i.e.,
!
!     idrsi(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrsi(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itrsi(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrsi(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrsi(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrsi(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Isoenergetic scattering
!   inverse mean free paths for radial zone j are stored at the corners
!   of unit cube j. Rates for zone j are interpolated from the rates 
!   stored at the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: idrsi
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: itrsi
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: iyrsi


!-----------------------------------------------------------------------
!  Isoenergetic scattering functions at cube corners
!-----------------------------------------------------------------------
!
!  cohsct(j,k,id,it,iy,ij_ray,ik_ray) : array of isoenergetic scattering
!   functions (which are the same for neutrinos of any type) for neutrinos
!   of energy k in radial zone j at the unit cube corners id, it, and iy
!   (id, it, iy = 1,2). Interpolations of isoenergetic scattering inverse
!   mean free paths are computed from this table.
!-----------------------------------------------------------------------


REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: cohsct


!-----------------------------------------------------------------------
!  Isoenergetic antiscattering functions at cube corners
!-----------------------------------------------------------------------
!
!  cohbsct(j,k,id,it,iy,ij_ray,ik_ray) : array of isoenergetic scattering
!   functions (which are the same for antineutrinos of any type) for
!   antineutrinos of energy k in radial zone j at the unit cube corners
!   id, it, and iy (id, it, iy = 1,2). Interpolations of isoenergetic
!   scattering inverse mean free paths are computed from this table.
!-----------------------------------------------------------------------


REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: cohbsct


!-----------------------------------------------------------------------
!  Interpolated isoenergetic scattering inverse mean free paths
!-----------------------------------------------------------------------
!
!  scti(j,k,n) : Isoenergetic scattering inverse mean free path
!
!      sctid(j,k,n)     = d(scti(j,k))/d(density)
!      sctit(j,k,n)     = d(scti(j,k))/d(temperature)
!      sctiy(j,k,n)     = d(scti(j,k))/d(electron fraction)
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)         :: scti
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)         :: sctid
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)         :: sctit
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)         :: sctiy


!-----------------------------------------------------------------------
!  Zero and first moments of the isoenergetic scattering inverse mean
!   free paths
!-----------------------------------------------------------------------
!  rmdnpsi  : ith moment of the neutrino-proton scattering inverse mean
!   free path
!  rmdnnsi  : ith moment of the neutrino-neutron scattering inverse
!   mean free path
!  rmdnbpsi : ith moment of the antineutrino-proton scattering inverse
!   mean free path
!  rmdnbnsi : ith moment of the antineutrino-neutron scattering inverse
!   mean free path
!  rmdnnsi  : ith moment of the neutrino-helium scattering inverse mean
!   free path   
!  rmdnnsi  : ith moment of the neutrino-nucleus scattering inverse mean
!   free path   
!-----------------------------------------------------------------------

REAL(KIND=double)                                        :: rmdnps0
REAL(KIND=double)                                        :: rmdnns0
REAL(KIND=double)                                        :: rmdnbps0
REAL(KIND=double)                                        :: rmdnbns0
REAL(KIND=double)                                        :: rmdnhes0
REAL(KIND=double)                                        :: rmdnhs0
REAL(KIND=double)                                        :: rmdnps1
REAL(KIND=double)                                        :: rmdnns1
REAL(KIND=double)                                        :: rmdnbps1
REAL(KIND=double)                                        :: rmdnbns1
REAL(KIND=double)                                        :: rmdnhes1
REAL(KIND=double)                                        :: rmdnhs1

END module scat_i_module
