SUBROUTINE dF_eta_deta( n_fermi, eta, df_n_deta )
!-----------------------------------------------------------------------
!
!    File:         dF_eta_deta
!    Module:       dF_eta_deta
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/05/02
!
!    Purpose:
!      To integrate the derivative of the Fermi function
!
!                             *   n    x - eta
!                            *   x dx e
!          dF (eta)/deta  =  *  ---------------
!            n               *    x - eta     2
!                           *   (e        + 1)
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  n         : Fermi function index
!  eta       : Fermi function argument
!
!    Output arguments:
!  df_n_deta : computed value of the Fermi function of type n 
!                 and argument eta
!
!    Variables that must be passed through common:
!      none
!
!    Modules:
!  numerical_modlue
!  FD_module
!
!-----------------------------------------------------------------------

USE numerical_module
USE FD_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)            :: n_fermi       ! Fermi-function index

REAL(KIND=double), INTENT(IN)  :: eta           ! chemical potential / kT

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT) :: df_n_deta     ! value of Fermi integral

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER             :: nleg = 32     ! order of Gauss-Legendre integration
INTEGER, PARAMETER             :: nlag = 24     ! order of Gauss-Laguerre integration

REAL(KIND=double)              :: df_n_deta1    ! value of Fermi integral for 0 < x < eta
REAL(KIND=double)              :: df_n_deta2    ! value of Fermi integral for eta < x < infinity
REAL(KIND=double)              :: xl            ! integration lower limit
REAL(KIND=double)              :: xu            ! integration upper limit

REAL(KIND=double), EXTERNAL    :: func_dFD_deta ! x**n_fermi*exp(x - eta)/( 1 + exp(x - eta) )**2

!-----------------------------------------------------------------------        
!-----------------------------------------------------------------------        

!-----------------------------------------------------------------------        
!  Initialize
!-----------------------------------------------------------------------        

n_fD            = n_fermi
eta_FD          = eta

!-----------------------------------------------------------------------        
!  eta > 0
!-----------------------------------------------------------------------        

IF (eta > zero) THEN

!........0 < x < eta

  xu            = eta
  xl            = zero
  CALL int_gauslegndr( func_dFD_deta, xl, xu, df_n_deta1, nleg )

!........eta < x < nfinity

  xl            = eta
  CALL int_gauslaguer( func_dFD_deta, xl, df_n_deta2, nlag ) 
  df_n_deta     = df_n_deta1 + df_n_deta2

!-----------------------------------------------------------------------        
!  eta <= 0
!-----------------------------------------------------------------------        

ELSE

  xl            = zero
  CALL int_gauslaguer( func_dFD_deta, xl, df_n_deta, nlag )

END IF

RETURN
END SUBROUTINE dF_eta_deta
