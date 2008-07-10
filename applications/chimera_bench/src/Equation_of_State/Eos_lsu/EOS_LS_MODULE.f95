!-----------------------------------------------------------------------
!    Module:       eos_ls_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE eos_ls_module

USE kind_module

SAVE


!.........Initial guess for Lattimer-Swesty EOS (cube corners)..........

!      inpvars(j,i,id,it,iy): initial guess of the ith input quantity for cube corner (id,it,iy) of 
!       radial zone j when calling  the Lattimer-Swesty equation of state. When possible, these 
!       quantities are saved from results of previous calls.

!          i = 1   : temperature
!          i = 2   : number density of nucleons in nuclei
!          i = 3   : (proton chemical potential)/kt
!          i = 4   : (neutron chemical potential)/kt

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: inpvars


!.........Initial guess for Lattimer-Swesty EOS (direct calls)..........

!      inpvarsa(j,i): initial guess of the ith input quantity of radial zone j
!       when calling  the Lattimer-Swesty equation of state. When possible, these 
!       quantities are saved from results of previous calls.

!          i = 1   : temperature
!          i = 2   : number density of nucleons in nuclei
!          i = 3   : (proton chemical potential)/kt
!          i = 4   : (neutron chemical potential)/kt

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)         :: inpvarsa

END module eos_ls_module
