SUBROUTINE F_eta( n_fermi, eta, f_n_eta )
!-----------------------------------------------------------------------
!
!    File:         F_eta
!    Module:       F_eta
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/21/02
!
!    Purpose:
!      To integrate the Fermi function
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
!  n       : Fermi function index
!  eta     : Fermi function argument
!
!    Output arguments:
!      f_n_eta : computed value of the Fermi function of type n 
!                 and argument eta
!
!    Variables that must be passed through common:
!      none
!
!    Modules:
!  kind_module, numerical_modlue
!  FD_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module

USE FD_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)  :: n_fermi                 ! Fermi-function index

REAL(KIND=double), INTENT(IN)  :: eta           ! chemical potential / kT

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT) :: f_n_eta       ! value of Fermi integral

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER   :: nleg = 32               ! order of Gauss-Legendre integration
INTEGER, PARAMETER   :: nlag = 24               ! order of Gauss-Laguerre integration

REAL(KIND=double)              :: f_n_eta1      ! value of Fermi integral for 0 < x < eta
REAL(KIND=double)              :: f_n_eta2      ! value of Fermi integral for eta < x < infinity
REAL(KIND=double), EXTERNAL    :: func_FD       ! x**n_fermi/( 1 + exp(x - eta) )
REAL(KIND=double)              :: xl            ! integration lower limit
REAL(KIND=double)              :: xu            ! integration upper limit

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

n_fD            = n_fermi
eta_FD          = eta

!-----------------------------------------------------------------------
!        eta > 0
!-----------------------------------------------------------------------

IF (eta > zero) THEN

!-----------------------------------------------------------------------
!  0 < x < eta
!-----------------------------------------------------------------------

xu              = eta
xl              = zero

CALL int_gauslegndr( func_FD, xl, xu, f_n_eta1, nleg )

!-----------------------------------------------------------------------
!  eta < x < nfinity
!-----------------------------------------------------------------------

xl              = eta

CALL int_gauslaguer( func_FD, xl, f_n_eta2, nlag )

f_n_eta         = f_n_eta1 + f_n_eta2

!-----------------------------------------------------------------------
!        eta <= 0
!-----------------------------------------------------------------------

ELSE

xl              = zero

CALL int_gauslaguer( func_FD, xl, f_n_eta, nlag )

END IF

RETURN
END SUBROUTINE F_eta
