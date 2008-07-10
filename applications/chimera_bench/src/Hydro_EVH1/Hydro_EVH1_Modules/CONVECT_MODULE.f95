!-----------------------------------------------------------------------
!    Module:       convect_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE convect_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Convection instability parameter
!-----------------------------------------------------------------------
!  aledoux(j): Ledoux parameter:
!
!                   1     drho            drho
!     aledoux(j) = --- [( ---- )      - ( ---- )     ]
!                  rho     dr   star       dr   blob
!
!     aledoux(j) > 0: radial zones j and j+1 are unstable to mixing
!      by Ledoux convection.
!
!     aledoux(j) < 0: radial zones j and j+1 are stable to mixing by
!      Ledoux convection.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: aledoux

!-----------------------------------------------------------------------
!  Quantities needed to compute the convective flow
!-----------------------------------------------------------------------
!  psclht(j) : the pressure scale height at outer boundary of radial zone j.
!
!  usound(j) : the sound velocity at outer boundary of radial zone j.
!
!  lmix(j)   : the mixing length at outer boundary of radial zone j.
!
!  alphlmix  :  multiple of pressure scale height used to compute lmix.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: psclht
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: usound
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: lmix
REAL(KIND=double)                                           :: alphlmix

!-----------------------------------------------------------------------
!  Quantities computed from the convective flow
!-----------------------------------------------------------------------
!  ulmixl(j)  : the average convective velocity at outer boundary of radial
!   as computed by mixing length prescription.
!
!  ulcnvct(j) : the average convective velocity at outer boundary of radial
!   zone j.
!
!  dmcnvct(j) :  the mass mixed across outer boundary of radial zone j.
!
!  pcnvct(j)  : the convective, turbulent pressure in radial zone j.
!
!  scnvct(j)  : the rate at which thermal energy is transferred to the matter
!   due to viscous dissipation of turbulent eddies in radial zone j.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: ulmixl
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: ulcnvct
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: dmcnvct
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: pcnvct
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)                :: scnvct

!-----------------------------------------------------------------------
!  Neutrino advection parameter
!-----------------------------------------------------------------------
!  adrftcv : neutrino advection parameter -
!
!   if
!         vdriftk < adrftcv*ulcnvct
!   where vdriftk is the drift velocity of neutrinos of energy zone k,
!   ulcnvct is the convective velocity, and adrftcv is an arbitrary parameter
!   of order unity, neutrinos are advected with the matter.
!
!   if
!         vdriftk > adrftcv*ulcnvct
!   neutrinos are not advected with the matter. In this case psi0aa(k)
!   is set equal to psi0(jf,k,n) so the the neutrinos of energy k delivered
!   to the zone are equal to those removed, resulting in no net advection.
!-----------------------------------------------------------------------

REAL(KIND=double)                                           :: adrftcv

END module convect_module
