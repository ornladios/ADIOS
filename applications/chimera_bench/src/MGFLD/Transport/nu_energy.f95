SUBROUTINE nu_energy( jr_min, jr_max, nx, nez, nnu, psi0, u_nu )
!-----------------------------------------------------------------------
!
!    File:         nu_energy
!    Module:       nu_energy
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/25/01
!
!    Purpose:
!      To compute the total neutrino energy density
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
!  psi0           : zero angular moments of the neutrino occupation number
!
!    Output arguments:
!  u_nu           : neutrino energy density (ergs cm^{-3})
!
!    Modules used:
!  kind_module, numerical_module
!  nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero

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

REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu) :: psi0 ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nx)     :: u_nu          ! neutrino energy density (ergs cm^{-3})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                              :: j             ! radial zone index
INTEGER                              :: n             ! neutrino flavor index

REAL(KIND=double), DIMENSION(nx,nnu) :: u_nu_n        ! neutrino energy density/flavor (ergs cm^{-3})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE THE NEUTRINO ENERGY DENSITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize enu
!-----------------------------------------------------------------------

u_nu                          = zero

!-----------------------------------------------------------------------
!  Compute neutrino energy density/flavor
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO j = jr_min,jr_max
    u_nu_n(j,n)               = SUM( ecoefa(j,:) * psi0(j,:,n) )
  END DO
END DO

!-----------------------------------------------------------------------
!  Compute total neutrino energy density
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  u_nu(j)                     = SUM( u_nu_n(j,:) * stwt(:) )
END DO

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE nu_energy
