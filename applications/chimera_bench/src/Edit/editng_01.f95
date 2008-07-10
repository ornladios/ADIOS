SUBROUTINE editng_01( n, jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editng_01
!    Module:       editng_01
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/05
!
!    Purpose:
!      To edit integral n-type neutrino transport data
!
!    Subprograms call:
!      date_and_time_print
!
!    Input arguments:
!  n          : neutrino flavor
!  jr_min     : inner radial zone of region for which configuration edit is to be made.
!  jr_max     : outer radial zone of region for which configuration edit is to be made.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!      none
!
!    Include files:
!      kind_module, array_module, numerical_module, physcnst_module
!      edit_module, mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nez, nnu
USE numerical_module, ONLY : zero, frpi, epsilon
USE physcnst_module, ONLY : ergfoe, rmu

USE edit_module, ONLY : prnttest, nprint, head, nlog
USE mdl_cnfg_module, ONLY : r, u, dmrst, ye
USE nu_dist_module, ONLY : unu, dunu, unue, dunue, psi0, psi1, unujcr, &
& nnujcr, unujrad, nnujrad, fluxnu, rjmh
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: n             ! neutrino flavor
INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status
INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: nn            ! neutrino flavor index

REAL(KIND=double)                :: psi0j2        ! second energy mioment of psi0
REAL(KIND=double)                :: psi0j3        ! third energy mioment of psi0
REAL(KIND=double)                :: psi0j4        ! forth energy mioment of psi0
REAL(KIND=double)                :: psi0j5        ! fifth energy mioment of psi0
REAL(KIND=double)                :: w2dw          ! unu(j)**2 * dunu(j)
REAL(KIND=double)                :: w3dw          ! unu(j)**3 * dunu(j)
REAL(KIND=double)                :: w4dw          ! unu(j)**4 * dunu(j)
REAL(KIND=double)                :: w5dw          ! unu(j)**5 * dunu(j)

REAL(KIND=double)                :: psi0j         ! psi0 interpolated to the zone edge
REAL(KIND=double)                :: e_density     ! the neutrino energy density
REAL(KIND=double)                :: e_flux        ! the neutrino flux
REAL(KIND=double)                :: e_mean        ! mean neutrino energy
REAL(KIND=double)                :: e_msq         ! mean square neutrino energy

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dunuj_trns   ! net energy in n-neutrinos transported into zone j (ergs)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dnnuj_trns   ! number in n-neutrinos transported into zone j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: denuj_trns   ! number in e-neutrinos transported into zone j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: danuj_trns   ! number in e-antineutrinos transported into zone j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: elecn        ! number in electrons residing in zone j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: totlpn       ! conserved number in leptons in zone j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: nnu_number   ! conserved number in leptons in zone j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: lum          ! luminosity of n-neutrinos (foes)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: e_rms        ! rms energy of n-neutrinos (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: inv_flxfct   ! inverse flux factor
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: pinch        ! pinch factor

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in editng_01')
 1101 FORMAT ('   j   n-nu e     n-nu n   n-nu eloss n-nu nloss  electrons   &
&totlpn    nunumber      lum       e-rms   /fluxfactor   pinch'/)
 1103 format (1x,i4,11(es11.3))
 2001 FORMAT (' Deallocation problem for array ',a10,' in editng_01')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ ALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (dunuj_trns(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunuj_trns'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dnnuj_trns(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dnnuj_trns'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (denuj_trns(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'denuj_trns'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (danuj_trns(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'denuj_trns'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (elecn(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'elecn     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (totlpn(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'totlpn    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnu_number(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnu_number'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (lum(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'lum       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_rms(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_rms     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (inv_flxfct(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'inv_flxfct'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pinch(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pinch     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

dunuj_trns          = zero
dnnuj_trns          = zero
denuj_trns          = zero
denuj_trns          = zero
elecn               = zero
totlpn              = zero
nnu_number          = zero
lum                 = zero
e_rms               = zero
inv_flxfct          = zero
pinch               = zero

!-----------------------------------------------------------------------
!
!                 \\\\\ COMPUTE EDIT QUANTITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Total number and energy of neutrinos in each zone
!-----------------------------------------------------------------------

DO nn = 1,nnu
  IF ( nnugp(nn) == 0 ) CYCLE
  CALL nu_number( jr_min, jr_max, nn, ij_ray, ik_ray, nx, nez, nnu, r, u, &
&  psi0, psi1 )
END DO

!-----------------------------------------------------------------------
!        Compute
!
!  dunuj_trns   : Net energy in n-neutrinos transported out of zone j
!  dnnuj_trns   : Net number in n-neutrinos transported out of  zone j
!  denuj_trns   : Net number in e-neutrinos transported out of  zone j
!  danuj_trns   : Net number in e-antineutrinos transported out of  zone j
!  elecn        : Net number in electrons residing in zone j
!  totlpn       : Net number in leptons residing in zone j - leptons
!                  transported in (Should remain constant)
!  nnu_number    : Number of n-neutrinos in zone j - net number transported in
!                  (should remain constant apart from production)
!-----------------------------------------------------------------------

DO j = jr_min,jr_max

  dunuj_trns(j)     = unujrad(j,n,ij_ray,ik_ray) - unujrad(j-1,n,ij_ray,ik_ray)
  dnnuj_trns(j)     = nnujrad(j,n,ij_ray,ik_ray) - nnujrad(j-1,n,ij_ray,ik_ray)
  denuj_trns(j)     = nnujrad(j,1,ij_ray,ik_ray) - nnujrad(j-1,1,ij_ray,ik_ray)
  danuj_trns(j)     = nnujrad(j,2,ij_ray,ik_ray) - nnujrad(j-1,2,ij_ray,ik_ray)
  elecn(j)          = dmrst(j) * ye(j)/rmu
  totlpn(j)         = nnujcr(j,1) + denuj_trns(j) + elecn(j) - nnujcr(j,2) - danuj_trns(j)
  nnu_number(j)     = nnujcr(j,n) + dnnuj_trns(j)

END DO

!-----------------------------------------------------------------------
!  n-neutrino luminosity
!-----------------------------------------------------------------------

CALL flux( jr_min, jr_max, n )
DO j = jr_min,jr_max
  lum(j)            = frpi * r(j) * r(j) * fluxnu(j,n) * ergfoe
END DO ! j

!-----------------------------------------------------------------------
!  n-neutrino rms energy
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  psi0j3            = zero
  psi0j5            = zero
  DO k = 1,nnugp(n)
    w3dw            = unu(j,k)**3 * dunu(j,k)
    w5dw            = unu(j,k)**5 * dunu(j,k)
    psi0j3          = psi0j3 + w3dw  * psi0(j,k,n)
    psi0j5          = psi0j5 + w5dw  * psi0(j,k,n)
  END DO ! k
  e_rms(j)          = DSQRT( DABS( psi0j5/( psi0j3 + epsilon ) ) )
END DO ! j

!-----------------------------------------------------------------------
!  inverse flux factors
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  e_density         = zero
  e_flux            = zero
  DO k = 1,nnugp(n)
    psi0j           = ( ( rjmh(j+1)**2 - r(j)**2 ) * psi0(j,k,n)         &
&                   + ( r(j)**2 - rjmh(j)**2 ) * psi0(j+1,k,n) )         &
&                   / ( rjmh(j+1)**2 - rjmh(j)**2 )
    e_density       = e_density + unue(j,k)**3 * dunue(j,k) * psi0j
    e_flux          = e_flux + unue(j,k)**3 * dunue(j,k) * psi1(j,k,n)
  END DO ! K
  inv_flxfct(j)     = e_density/( e_flux + epsilon )
END DO ! j

!-----------------------------------------------------------------------
!  pinch factors
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  psi0j2            = zero
  psi0j3            = zero
  psi0j4            = zero
  DO k = 1,nnugp(n)
    w2dw            = unu(j,k)**2 * dunu(j,k)
    w3dw            = unu(j,k)**3 * dunu(j,k)
    w4dw            = unu(j,k)**4 * dunu(j,k)
    psi0j2          = psi0j2 + w2dw  * psi0(j,k,n)
    psi0j3          = psi0j3 + w3dw  * psi0(j,k,n)
    psi0j4          = psi0j4 + w4dw  * psi0(j,k,n)
  END DO ! k
  e_mean            = psi0j3/( psi0j2 + epsilon )
  e_msq             = psi0j4/( psi0j2 + epsilon )
  pinch(j)          = 0.75d0 * e_msq/( e_mean * e_mean + epsilon )
END DO ! j

!-----------------------------------------------------------------------
!
!                       \\\\\ PRINT /////
!
!-----------------------------------------------------------------------

WRITE (nprint,1101)

DO j = jr_max,jr_min,-1
  WRITE (nprint,1103) j, unujcr(j,n), nnujcr(j,n), dunuj_trns(j), dnnuj_trns(j), &
&  elecn(j), totlpn(j), nnu_number(j), lum(j), e_rms(j), inv_flxfct(j), pinch(j)
END DO

!-----------------------------------------------------------------------
!
!                  \\\\\ DEALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

DEALLOCATE (dunuj_trns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunuj_trns'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dnnuj_trns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dnnuj_trns'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (denuj_trns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'denuj_trns'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (danuj_trns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'denuj_trns'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (elecn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'elecn     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (totlpn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'totlpn    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (nnu_number, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnu_number '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (lum, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'lum       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (e_rms, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_rms     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (inv_flxfct, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'inv_flxfct'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (pinch, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pinch     '; WRITE (nlog,2001) var_name; END IF



RETURN
END SUBROUTINE editng_01
