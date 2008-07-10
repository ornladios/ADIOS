SUBROUTINE enutcp( n, j, nu_t, nu_eta )
!-----------------------------------------------------------------------
!
!    File:         enutcp
!    Module:       enutcp
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/1/00
!
!    Purpose:
!      To compute the temperature and chemical potential of a Fermi-Dirac n-neutrino distribution
!       having the same zero, first, and second moments as the actual n-neutrino distribution.
!
!    Subprograms called:
!      f_e2ee_p, F_eta
!
!    Input arguments:
!  n           : neutrino type
!  j           : radial zone index
!
!    Output arguments:
!  nu_t        : temperature of the fit fermi-dirac n-neutrino distribution (MeV)
!  nu_eta      : (chemical potantial)/kT of the fit fermi-dirac n-neutrino distribution
!
!    Variables that must be passed through common:
!
!  nprint      : unit number of print file.
!  nnugp(n)    : number of energy zones for neutrinos of type n
!  unu(j,k)    : energy of energy zone k at radial zone j (MeV)
!  dunu(j,k)   : energy width of energy zone k at radial zone j (MeV)
!  psi0(j,k,n) : zeroth moment of of the neutrino occupation probability for neutrinos of type n, 
!                energy zone k, radial zone j
!
!    Modules used:
!  kind_module, numerical_module
!  edit_module, nu_dist_module, nu_energy_grid_module, 
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half, epsilon

USE edit_module, ONLY : nprint
USE nu_dist_module, ONLY : unu, psi0, ncoefa, enucpa, enuta
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)            :: j             ! radial zone index
INTEGER, INTENT(IN)            :: n             ! neutrino degree of freedom index

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT) :: nu_t          ! Fermi-Dirac temperature
REAL(KIND=double), INTENT(OUT) :: nu_eta        ! Fermi-Dirac eta

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER             :: itmax = 24    ! maximum number of iterations
INTEGER                        :: it            ! iteration index
INTEGER                        :: k             ! energy zone index
INTEGER                        :: nf            ! integer parameter for passing to F_eta

REAL(KIND=double), EXTERNAL    :: f_e2ee        ! F2*F4/F3**2
REAL(KIND=double)              :: f_e2ee_p      ! F2*F4/F3**2

REAL(KIND=double)              :: e0            ! zero energy momemt of ncoef * psi0
REAL(KIND=double)              :: e1            ! first energy momemt of ncoef * psi0
REAL(KIND=double)              :: e2            ! second energy momemt of ncoef * psi0
REAL(KIND=double)              :: e2ee          ! e2_mean/e_mean**2
REAL(KIND=double)              :: e2_mean       ! e2/e0
REAL(KIND=double)              :: e_mean        ! e1/e0
REAL(KIND=double)              :: f2            ! value returned by F_eta
REAL(KIND=double)              :: f3            ! value returned by F_eta
REAL(KIND=double)              :: dfeta         ! f_e2ee_p - e2ee
REAL(KIND=double)              :: enucmpmx      ! starting upper bound of nu_eta for bisection
REAL(KIND=double)              :: enucmpmn      ! starting lower bound of nu_eta for bisection

REAL(KIND=double), PARAMETER   :: e2ee_min = 1.0d+00    ! minimum e2ee for computing tolerance
REAL(KIND=double), PARAMETER   :: nu_eta_max = 40.d0   ! maximum value of nu_eta
REAL(KIND=double), PARAMETER   :: nu_eta_min = -40.d0  ! minimum value of nu_eta
REAL(KIND=double), PARAMETER   :: tol = 1.d-4           ! iteration tolerance

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 2001 format (' nu_eta will not converge, nu_eta=',1pe14.7,' enucmpmn=',1pe14.7,    &
&             ' enucmpmx=',1pe14.7,' f_e2ee_p',1pe14.7,' e2ee=',1pe14.7)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) THEN
   nu_t         = zero
   nu_eta       = zero
   RETURN
END IF ! nnugp(n) = 0

!-----------------------------------------------------------------------
!  Calculate e_mean and <e2>/<e><e>
!-----------------------------------------------------------------------

e0            = zero
e1            = zero
e2            = zero

DO k = 1,nnugp(n)
  e0          = e0 + ncoefa(j,k)               * psi0(j,k,n)
  e1          = e1 + ncoefa(j,k) * unu(j,k)    * psi0(j,k,n)
  e2          = e2 + ncoefa(j,k) * unu(j,k)**2 * psi0(j,k,n)
END DO

e_mean        = e1/( e0 + epsilon )
e2_mean       = e2/( e0 + epsilon )
e2ee          = e2_mean/( e_mean**2 + epsilon )

!-----------------------------------------------------------------------
!  Fit nu_t andf nu_eta
!-----------------------------------------------------------------------

SELECT CASE (n)

!-----------------------------------------------------------------------
!  n = 1 or 2, fit nu_eta
!-----------------------------------------------------------------------

CASE (1,2)

  enucmpmx     = nu_eta_max
  enucmpmn     = nu_eta_min

  DO it = 1,itmax

    nu_eta    = half * ( enucmpmn + enucmpmx )
    f_e2ee_p = f_e2ee(nu_eta)

    dfeta     = f_e2ee_p - e2ee

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

    IF ( dabs(dfeta) < tol * dmax1( dabs(e2ee), e2ee_min ) ) EXIT

!-----------------------------------------------------------------------
!  Change bounds if not converged
!-----------------------------------------------------------------------
    
    IF ( dfeta > zero ) enucmpmn = nu_eta
    IF ( dfeta < zero ) enucmpmx = nu_eta

!-----------------------------------------------------------------------
!  Test for iteration limit
!-----------------------------------------------------------------------

!    IF ( it == itmax) WRITE (nprint,2001) nu_eta,enucmpmn,enucmpmx,f_e2ee_p,e2ee

  END DO

  enucpa(j,n)  = nu_eta

!-----------------------------------------------------------------------
!  n = 1 or 2, fit nu_t
!-----------------------------------------------------------------------

  nf          = 2
  CALL F_eta( nf, nu_eta, f2 )
  nf          = 3
  CALL F_eta( nf, nu_eta, f3 )
  nu_t        = e_mean * f2/f3
  enuta(j,n)  = nu_t

!-----------------------------------------------------------------------
!  n = 3, fit nu_t and set nu_eta
!-----------------------------------------------------------------------

CASE (3)

  nu_eta      = zero
  enucpa(j,n) = zero
  nu_t        = e_mean/3.d0
  enuta(j,n)  = nu_t

END SELECT

RETURN
END SUBROUTINE enutcp
