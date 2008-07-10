MODULE nuclear_data
!-----------------------------------------------------------------------------
!  This module contains the essential data for each included specie.
!  aa, zz, & nn are the mass, proton and neutron numbers, be is the binding
!  energy, nname is the 5 character nuclear name (nname(0) is placeholder for
!  non-nuclei).  Their array sizes are set in the routine format_nuclear_data.  
!-----------------------------------------------------------------------------

USE kind_module, ONLY : double

CHARACTER(len=5), ALLOCATABLE, DIMENSION(:)    :: nname

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: aa,zz,nn,be

END MODULE nuclear_data
