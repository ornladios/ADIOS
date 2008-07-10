SUBROUTINE enutcp_2( n, j, nu_t, nu_eta )
!-----------------------------------------------------------------------
!
!    File:         enutcp_2
!    Module:       enutcp_2
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/4/02
!
!    Purpose:
!      To compute the temperature and chemical potential of a Fermi-Dirac n-neutrino distribution
!       having the same number density (zero moment) and the same energy density (first moment) as
!       the actual n-neutrino distribution.
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
!  kind_module, array_module, numerical_module
!  edit_module, nu_dist_module, nu_energy_grid_module, 
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY: nx, nez
USE numerical_module, ONLY: zero, half

USE edit_module, ONLY: nprint, nlog
USE nu_dist_module, ONLY: psi0, unu, dunu
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE
 
!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------
     
INTEGER, INTENT(in)               :: j             ! radial zone index
INTEGER, INTENT(in)               :: n             ! neutrino degree of freedom index

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)    :: nu_t          ! Fermi-Dirac temperature
REAL(KIND=double), INTENT(out)    :: nu_eta        ! Fermi-Dirac eta

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

INTEGER, PARAMETER                :: itmax = 30    ! maximum number of iterations
INTEGER                           :: i             ! Newton-Rhapson iteration index
INTEGER                           :: it            ! t iteration index
INTEGER                           :: ieta          ! eta iteration index
INTEGER                           :: k             ! energy zone index
INTEGER                           :: nf            ! integer parameter for passing to F_eta
INTEGER                           :: istat         ! allocation status
INTEGER, DIMENSION(2)             :: jc            ! integer array for gjr

REAL(KIND=double), DIMENSION(2,3) :: a             ! coefficient matrix in Newton-Rhapson iteration
REAL(KIND=double), DIMENSION(2)   :: v             ! vector of quantities for subroutine gjr

REAL(KIND=double)                 :: e0            ! zero energy momemt of ncoef * psi0
REAL(KIND=double)                 :: e1            ! first energy momemt of ncoef * psi0
REAL(KIND=double), PARAMETER      :: e0_min = 1.d-50 ! return if e0 < e0_min
REAL(KIND=double), PARAMETER      :: e1_min = 1.d-50 ! return if e1 < e1_min
REAL(KIND=double), PARAMETER      :: tol = 1.d-5   ! return if e1 < e1_min

REAL(KIND=double)                 :: e0_it         ! e0 calculated by iterated t and eta
REAL(KIND=double)                 :: e1_it         ! e1 calculated by iterated t and eta
REAL(KIND=double)                 :: de0_it_dt     ! d(e0)/dT
REAL(KIND=double)                 :: de0_it_deta   ! d(e0)/deta
REAL(KIND=double)                 :: de1_it_dt     ! d(e1)/dT
REAL(KIND=double)                 :: de1_it_deta   ! d(e1)/deta
REAL(KIND=double)                 :: f2            ! value returned by F_eta for index 2
REAL(KIND=double)                 :: f3            ! value returned by F_eta for index 3
REAL(KIND=double)                 :: df2_deta      ! value returned by dF_eta_deta for index 2
REAL(KIND=double)                 :: df3_deta      ! value returned by dF_eta_deta for index 3
REAL(KIND=double)                 :: dnu_t         ! correction to nu_t
REAL(KIND=double)                 :: dnu_eta       ! correction to nu_eta

REAL(KIND=double), PARAMETER      :: nu_tmax = 100.d0      ! Maximum value of nu_t
REAL(KIND=double), PARAMETER      :: nu_tmin = 1.d-1       ! Minimum value of nu_t
REAL(KIND=double), PARAMETER      :: nu_etamax = 200.d0    ! Maximum value of nu_eta
REAL(KIND=double), PARAMETER      :: nu_etamin = -30.d0    ! Minimum value of nu_eta

REAL(KIND=double)                 :: nu_tmaxb      ! upper bound of nu_t during bisection
REAL(KIND=double)                 :: nu_tminb      ! lower bound of nu_t during bisection
REAL(KIND=double)                 :: nu_etamaxb    ! upper bound of nu_eta during bisection
REAL(KIND=double)                 :: nu_etaminb    ! lower bound of nu_eta during bisection

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: nu_ta
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: nu_etaa

  201 FORMAT (' enut and enucmp will not converge in subroutine enutcp; n = ',i3, &
&              ' j=',i3,'; will try bisection method')
  301 FORMAT (' The  n-neutrino number will not converge; n=',i2,' j=',i3, &
&              ' nu_etaminb=',1pe11.4,' nu_etamaxb=',1pe11.4)
  303 FORMAT (' e0=',1pe11.4,' e0_it=',1pe11.4,' e1=',1pe11.4,' e1_it=',1pe11.4)
  305 FORMAT (' The n-neutrino energy will not converge;  n=',i2,' j=',i3, &
&              ' nu_tminb=',1pe11.4,' nu_tmaxb=',1pe11.4,' enucmp=',1pe11.4)
  307 FORMAT (' e0=',1pe11.4,' e0_it=',1pe11.4,' e1=',1pe11.4,' e1_it=',1pe11.4)
 1001 FORMAT (' Allocation problem for array ',a10,' in enutcp_2')
 2001 FORMAT (' Deallocation problem for array ',a10,' in enutcp_2')

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

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
!  Calculate the zero and first energy momemts of ncoef * psi0
!-----------------------------------------------------------------------

e0             = zero
e1             = zero

DO k = 1,nnugp(n)
  e0           = e0 + unu(j,k)**2 * dunu(j,k) * psi0(j,k,n)
  e1           = e1 + unu(j,k)**3 * dunu(j,k) * psi0(j,k,n)
END DO

!-----------------------------------------------------------------------
!  Return if values of enun or enue are too small
!-----------------------------------------------------------------------

IF ( e0 < e0_min  .or.  e1 < e1_min ) THEN
  nu_t         = zero
  nu_eta       = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (nu_ta(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_ta     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_etaa(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_etaa   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize arrays
!-----------------------------------------------------------------------

nu_ta            = zero
nu_etaa          = zero

!-----------------------------------------------------------------------
!
!                ||||| Newton-Rhapson iteration |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Try a Newton-Rhapson iteration if previous values of nu_t and nu_eta
!-----------------------------------------------------------------------

N_R_it: IF ( nu_ta(j,n)  /= zero  .and.  nu_etaa(j,n) /= zero ) THEN

!-----------------------------------------------------------------------
!  Initial guess of nu_t and nu_eta
!-----------------------------------------------------------------------

  nu_t           = nu_ta(j,n)
  nu_eta         = nu_etaa(j,n)

!-----------------------------------------------------------------------
!  Begin the Newton-Rhapson iteration
!-----------------------------------------------------------------------

  NR_iteration: DO i = 1,itmax

!-----------------------------------------------------------------------
!  Set up and solve simultaneous equations
!-----------------------------------------------------------------------

!         a(i,1)*denut + a(i,2)*denucp = a(i,3)
    nf          = 2
    CALL F_eta( nf, nu_eta, f2 )
    nf          = 3
    CALL F_eta( nf, nu_eta, f3 )
    nf          = 2
    CALL dF_eta_deta( nf, nu_eta, df2_deta )
    nf          = 3
    CALL dF_eta_deta( nf, nu_eta, df3_deta )

    e0_it        = nu_t**3 * f2
    e1_it        = nu_t**4 * f3
    de0_it_dt    = 3.d0 * nu_t**2 * f2
    de0_it_deta  = nu_t**3 * df2_deta
    de1_it_dt    = 4.d0 * nu_t**3 * f3
    de1_it_deta  = nu_t**4 * df3_deta

    a(1,1)       = de0_it_dt
    a(1,2)       = de0_it_deta
    a(1,3)       = e0 - e0_it
    a(2,1)       = de1_it_dt
    a(2,2)       = de1_it_deta
    a(2,3)       = e1 - e1_it

    v(1)         = 4.d0
    call gjr( a, 3, 2, 2, 3, *2050, jc, v )

    dnu_t        = a(1,3)
    dnu_eta      = a(2,3)

!-----------------------------------------------------------------------
!  Correct nu_t and nu_eta
!-----------------------------------------------------------------------

    nu_t         = nu_t   + dnu_t
    nu_eta       = nu_eta + dnu_eta
    nu_t         = min( nu_tmax  , max( nu_tmin , nu_t    ) )
    nu_eta       = min( nu_etamax, max( nu_etamin, nu_eta ) )

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

test: IF ( abs(dnu_t) <= tol * abs(nu_t)  .and.  abs(dnu_eta) <= tol * abs(nu_eta) ) THEN
      nu_ta(j,n)    = nu_t
      nu_etaa(j,n)  = nu_eta

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

      DEALLOCATE (nu_ta, STAT = istat)
        IF ( istat /= 0 ) THEN; var_name = 'nu_ta     '; WRITE (nlog,2001) var_name; END IF
      DEALLOCATE (nu_etaa, STAT = istat)
        IF ( istat /= 0 ) THEN; var_name = 'nu_etaa   '; WRITE (nlog,2001) var_name; END IF

      RETURN
    END IF test
      
  END DO NR_iteration
  
END IF N_R_it 

 2050 CONTINUE

!-----------------------------------------------------------------------
!
! ||||| Bisection method (Used if newton-rhapson fails to converge) |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Iterate to determine tnu
!-----------------------------------------------------------------------

nu_tmaxb     = nu_tmax
nu_tminb     = nu_tmin

nu_t_iteration: DO it = 1,itmax

  nu_t       = half * ( nu_tmaxb + nu_tminb )

!-----------------------------------------------------------------------
!  iterate to determine nu_eta ( = nu_cmp/kt )
!-----------------------------------------------------------------------

  nu_etamaxb = nu_etamax
  nu_etaminb = nu_etamin

nu_eta_iteration: DO ieta = 1,itmax

    nu_eta   = half * ( nu_etamaxb + nu_etaminb )

    nf          = 2
    CALL F_eta( nf, nu_eta, f2 )
    nf          = 3
    CALL F_eta( nf, nu_eta, f3 )

    e0_it        = nu_t**3 * f2
    e1_it        = nu_t**4 * f3

!-----------------------------------------------------------------------
!  Test for nu_eta convergence
!-----------------------------------------------------------------------

    IF ( abs(e0_it - e0) < tol * e0 ) EXIT

!-----------------------------------------------------------------------
!  Move eta boundary if not converged
!-----------------------------------------------------------------------

    IF ( e0_it < e0 ) then
      nu_etaminb   = nu_eta
    ELSE
      nu_etamaxb   = nu_eta
    END IF

!-----------------------------------------------------------------------
!  Print diagnostics if failure to converge
!-----------------------------------------------------------------------
    
    IF ( ieta == itmax ) THEN
!      WRITE (nprint,301) n,j,nu_etaminb,nu_etamaxb
!      WRITE (nprint,303) e0,e0_it,e1,e1_it
      nu_t    = zero
      nu_eta  = zero

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

      DEALLOCATE (nu_ta, STAT = istat)
        IF ( istat /= 0 ) WRITE (nlog,2001)
      DEALLOCATE (nu_etaa, STAT = istat)
        IF ( istat /= 0 ) WRITE (nlog,2001)

      RETURN
    END IF

  END DO nu_eta_iteration

!-----------------------------------------------------------------------
!  Test for nu_t convergence
!-----------------------------------------------------------------------

  IF ( abs(e1_it - e1) < tol * e1 ) EXIT

!-----------------------------------------------------------------------
!  Move t boundary if not converged
!-----------------------------------------------------------------------

  IF ( e1_it < e1 ) THEN
    nu_tminb       = nu_t
  ELSE
    nu_tmaxb       = nu_t
  END IF

!-----------------------------------------------------------------------
!  Print diagnostics if failure to converge
!-----------------------------------------------------------------------
  
  IF ( it == itmax ) THEN
    WRITE (nprint,305) n,j,nu_tminb,nu_tmaxb
    WRITE (nprint,307) e0,e0_it,e1,e1_it
    nu_t    = zero
    nu_eta  = zero

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

    DEALLOCATE (nu_ta, STAT = istat)
      IF ( istat /= 0 ) THEN; var_name = 'nu_ta     '; WRITE (nlog,2001) var_name; END IF
    DEALLOCATE (nu_etaa, STAT = istat)
      IF ( istat /= 0 ) THEN; var_name = 'nu_etaa   '; WRITE (nlog,2001) var_name; END IF

    RETURN
  END IF

END DO nu_t_iteration

nu_ta(j,n)   = nu_t 
nu_etaa(j,n) = nu_eta

!........Deallocate arrays..............................................

DEALLOCATE (nu_ta, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_ta     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (nu_etaa, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_etaa   '; WRITE (nlog,2001) var_name; END IF


RETURN
END SUBROUTINE enutcp_2
