SUBROUTINE nu_flux( jr_min, jr_max, nx, nez, nnu, psi1, f_nu )
!-----------------------------------------------------------------------
!
!    File:         nu_flux
!    Module:       nu_flux
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/25/01
!
!    Purpose:
!      To compute the total neutrino energy flux
!
!    Variables that must be passed through common:
!  ecoefa(j,k)    : 4.*pi/((h*c)**3)*w**3*dw
!  stwt(n)        : satistical weight of n-neutrinos 
!
!    Subprograms called:
!        pre_trans
!
!    Input arguments:
!  jr_min         : shifted lower radial zone index
!  jr_max         : shifted upper radial zone index
!  nx             : x-array extent
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!  psi1           : first angular moments of the neutrino occupation number
!
!    Output arguments:
!  f_nu           : neutrino energy flux (ergs cm^{-2} s^{-1})
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module
!  nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero
USE physcnst_module, ONLY: cvel

USE nu_dist_module, ONLY: ecoefa, stwt

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                  :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)                  :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)                  :: nx            ! x-array extent
INTEGER, INTENT(in)                  :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)                  :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu) :: psi1 ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nx)     :: f_nu          ! neutrino energy flux (ergs cm^{-2} s^{-1})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                              :: j             ! radial zone index
INTEGER                              :: n             ! neutrino flavor index

REAL(KIND=double), DIMENSION(nx,nnu) :: f_nu_n        ! neutrino energy flux/flavor (ergs cm^{-2} s^{-1})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE THE NEUTRINO ENERGY DENSITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize fnu; this will set f_nu(jr_min-1) to zero
!-----------------------------------------------------------------------

f_nu                          = zero

!-----------------------------------------------------------------------
!  Compute neutrino energy flux/flavor
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO j = jr_min,jr_max
    f_nu_n(j,n)               = SUM( ecoefa(j,:) * psi1(j,:,n) )
  END DO
END DO

!-----------------------------------------------------------------------
!  Compute total neutrino energy flux
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  f_nu(j)                     = SUM( f_nu_n(j,:) * stwt(:) ) * cvel
END DO

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE nu_flux
