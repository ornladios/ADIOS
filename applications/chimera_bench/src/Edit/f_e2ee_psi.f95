SUBROUTINE f_e2ee_psi( j, eta, f_e2ee )
!-----------------------------------------------------------------------
!
!    File:         f_e2ee_psi
!    Module:       f_e2ee_psi
!    Type:         Subrogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/22/02
!
!    Purpose:
!      To calculate F2*F4/F3**2 using the energy zoning of the mgfld code.
!
!    Subprograms called:
!  F_eta_psi
!
!    Input arguments:
!      none
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!      none
!
!    Modules:
!  kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: j            ! radial zone index

REAL(KIND=double), INTENT(in)  :: eta          ! chemical potential / kT

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT) :: f_e2ee       ! function value

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                        :: n            ! Fermi-function index

REAL(KIND=double)              :: f2           ! Fermi integral of order 2
REAL(KIND=double)              :: f3           ! Fermi integral of order 3
REAL(KIND=double)              :: f4           ! Fermi integral of order 4

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

n           = 2
CALL F_eta_psi( j, n, eta, f2 )
n           = 3
CALL F_eta_psi( j, n, eta, f3 )
n           = 4
CALL F_eta_psi( j, n, eta, f4 )
f_e2ee      = f2 * f4/f3**2

RETURN
END SUBROUTINE f_e2ee_psi
