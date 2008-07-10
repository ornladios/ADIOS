MODULE cross_sect_data
!-----------------------------------------------------------------------------
!  This module contains the data needed to calculate the cross sections
!  nreac(i) are the # of reactions with i reactants.  
!  n?i list the nuclei affected by each reaction.  
!  The csect variables are the results.
!  The rc variables are data for the temperature dependant part of the csect.
!-----------------------------------------------------------------------------

USE kind_module, ONLY : double

INTEGER                                        :: nreac(3)

INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: n1i,n2i,n3i
INTEGER, ALLOCATABLE, DIMENSION(:)             :: iwk1,iwk2,iwk3    ! Weak Reaction
INTEGER, ALLOCATABLE, DIMENSION(:)             :: ires1,ires2,ires3 ! Resonant React
INTEGER, ALLOCATABLE, DIMENSION(:)             :: irev1,irev2,irev3 ! Reverse React
INTEGER, ALLOCATABLE, DIMENSION(:)             :: iffn ! Maps reaction to FFN list

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: csect1,csect2,csect3
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: rc1,rc2,rc3   ! dim(7,nreac)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rpf1,rpf2,h2,h3 ! ratios of PartF.
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: q1,q2,q3 ! reaction Q values

END MODULE cross_sect_data
