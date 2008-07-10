!-----------------------------------------------------------------------
!    Module:       nu_energy_grid_module
!    Author:       S. W. Bruenn
!    Date:         7/26/04
!-----------------------------------------------------------------------

MODULE nu_energy_grid_module

USE kind_module, ONLY : double

SAVE

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------
!  nnugp(n) : the number of energy zones for neutrinos of type n.
!
!  nnugpmx : the maximum number of energy zones for any neutrino flavor.
!
!  unui(k) : the value of the zone-centered neutrino energy at infinity
!   of energy zone k (MeV).
!
!  dunui(k) : the width of energy zone k at infinity (MeV),
!   dunui(k) = unubi(k+1) - unubi(k) 
!
!  unubi(k) : the value of the inner zone-edge neutrino energy at
!   infinity of energy zone k (MeV).
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:)                    :: nnugp
INTEGER                                               :: nnugpmx

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: unui
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: dunui
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: unubi


END module nu_energy_grid_module
