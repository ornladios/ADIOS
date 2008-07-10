SUBROUTINE F_eta_psi( j, n_fermi, eta, f_n_eta )
!-----------------------------------------------------------------------
!
!    File:         F_eta_psi
!    Module:       F_eta_psi
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/21/02
!
!    Purpose:
!      To integrate the Fermi function using the energy zoning of the mgfld code
!
!                      *     n
!                     *     x dx
!          F (eta) =  *  ------------
!           n         *   x - eta
!                    *   e        + 1
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  eta     : Fermi function argument
!
!    Output arguments:
!  f_n_eta : computed value of the Fermi function of type n 
!             and argument eta
!
!    Variables that must be passed through common:
!      none
!
!    Modules:
!  kind_module, numerical_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one

USE nu_dist_module, ONLY : unu, dunu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)  :: n_fermi                 ! Fermi-function index
INTEGER, INTENT(IN)  :: j                       ! radial zone index

REAL(KIND=double), INTENT(IN)  :: eta           ! chemical potential / kT

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT) :: f_n_eta       ! value of Fermi integral

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER              :: k                       ! neutrino energy index

REAL(KIND=double), EXTERNAL    :: fexp          ! exp with upper and lower bounds

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

f_n_eta         = zero

!-----------------------------------------------------------------------
!  Return if nnugpmx = 0 for all n
!-----------------------------------------------------------------------

nnugpmx         = 0
IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!  Integrate the Fermi function
!-----------------------------------------------------------------------

Do k = 1,nnugpmx

f_n_eta         = f_n_eta + unu(j,k)**n_fermi * dunu(j,k) * ( one/( fexp( unu(j,k) - eta ) + one ) )

END DO

RETURN
END SUBROUTINE F_eta_psi
