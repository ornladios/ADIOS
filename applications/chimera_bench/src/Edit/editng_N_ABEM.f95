SUBROUTINE editng_N_ABEM( n, jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editng_N_ABEM
!    Module:       editng_N_ABEM
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/05
!
!    Purpose:
!      To edit integral n-type neutrino transport data pertaining to
!       neujtrino absorption and emission on nucleons.
!
!    Subprograms call:
!  paired     : Computes quantities needed for editing PRA rates
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
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nez, nnu
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : cvel, rmu, kmev

USE abem_module, ONLY : emis, absor
USE edit_module, ONLY : prnttest, nprint, nlog, head, idxeng, dudt_ABEM
USE mdl_cnfg_module, ONLY : rho, t, ye
USE nu_dist_module, ONLY : ncoefa, ecoefa, psi0, psi1, unu
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor
INTEGER, INTENT(in)               :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

INTEGER                           :: istat         ! allocation status
INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy index

REAL(KIND=double)                 :: t_mev         ! temperature (MeV)

REAL(KIND=double)                 :: dpsirt        ! rate of change in psi0 due to PRA neutrino production
REAL(KIND=double)                 :: rentti        ! rate of neutrinos pruduced by PRA (g^{-1} s^{-1})
REAL(KIND=double)                 :: uentti        ! rate of energy chenge due to PRA (MeV g^{-1} s^{-1})

REAL(KIND=double)                 :: dpsira        ! rate of change in psi0 due to PRA neutrino annihilation
REAL(KIND=double)                 :: rentai        ! rate of neutrinos annihilation by PRA (g^{-1} s^{-1})
REAL(KIND=double)                 :: uentai        ! rate of energy chenge due to PRA neutrino annhihlation (MeV g^{-1} s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ravetp   ! neutrino production rate per barion due to PRA (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ravepa   ! neutrino annihilation rate per barion due to PRA (s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eavetp   ! Mean energy of neutrinos produced due to PRA (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eavepa   ! Mean energy of neutrinos annihilated due to PRA (MeV)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: urttp    ! energy transfer rate per baryon due to neutrino production by PRA (MeV s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: urtpa    ! energy transfer rate per baryon due to neutrino annihilation by PRA (MeV s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: urtpp    ! net energy transfer rate per baryon due to PRA (MeV s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dstpdt   ! d(entropy)/k-n)/dt by neutrino production due to PRA (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dspadt   ! d(entropy)/k-n)/dt by neutrino annihilation due to PRA (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dstadt   ! Net d(entropy)/k-n)/dt due to PRA (s^{-1})

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in editng_PR')
 1501 FORMAT (10x,'Enission and Absorption on Nucleons: Neutrino Number,&
 & Energy, and Matter Entropy Production Rates'/)
 1503 FORMAT ('   j  dndt_p/n   dndt_a/n   <e_nu>_p   <e_nu>_a   dudt_em/n &
 & dudt_ab/n dudt_net/n  dsdt_en/n  dsdt_ab/n dsdt_net/n'/)
 1505 FORMAT (1x,i4,10(es11.3))
 2001 FORMAT (' Deallocation problem for array ',a10,' in editng_PR')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) RETURN

!-----------------------------------------------------------------------
!
!                   \\\\\ ALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (ravetp(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravetp    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ravepa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravepa    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (eavetp(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavetp    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (eavepa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavepa    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (urttp(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'urttp     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (urtpa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'urtpa     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (urtpp(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'urtpp     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dstpdt(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dstpdt    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dspadt(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dspadt    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dstadt(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dstadt    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

ravetp                  = zero
ravepa                  = zero

eavetp                  = zero
eavepa                  = zero

urttp                   = zero
urtpa                   = zero
urtpp                   = zero

dstpdt                  = zero
dspadt                  = zero
dstadt                  = zero

!-----------------------------------------------------------------------
!
!                 \\\\\ COMPUTE EDIT QUANTITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute emission and absorption rates
!-----------------------------------------------------------------------

CALL abemrate( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )

DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

  rentti                = zero
  uentti                = zero
  rentai                = zero
  uentai                = zero

  t_mev                 = kmev * t(j)

  DO k = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Total paricle and energy change rate due to nucleon emission
!-----------------------------------------------------------------------

    dpsirt              = emis(j,k,n,ij_ray,ik_ray) * ( one - psi0(j,k,n) ) * cvel
    rentti              = rentti +            ncoefa(j,k) * dpsirt/rho(j)
    uentti              = uentti + unu(j,k) * ncoefa(j,k) * dpsirt/rho(j)

!-----------------------------------------------------------------------
!  Total paricle and energy change rate due to nucleon absorption
!-----------------------------------------------------------------------

    dpsira              = absor(j,k,n,ij_ray,ik_ray) * psi0(j,k,n) * cvel
    rentai              = rentai +            ncoefa(j,k) * dpsira/rho(j)
    uentai              = uentai + unu(j,k) * ncoefa(j,k) * dpsira/rho(j)

  END DO ! k

!-----------------------------------------------------------------------
!
!   ||||| Emission and absorption rates per baryon due to PRA |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino production rate per baryon by nucleon emission
!-----------------------------------------------------------------------

  ravetp(j)             = rentti * rmu

!-----------------------------------------------------------------------
!  Neutrino annihilation rate per baryon by nucleon absorption
!-----------------------------------------------------------------------

  ravepa(j)             = rentai * rmu

!-----------------------------------------------------------------------
!
! ||||| Mean neutrino energies in nucleon emission and absorption |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Mean energy of neutrinos produced due to nucleon emission
!-----------------------------------------------------------------------

  eavetp(j)             = uentti/SIGN( ( DABS( rentti ) + epsilon ), rentti )

!-----------------------------------------------------------------------
!  Mean energy of neutrinos annihilated due to nucleon absorption
!-----------------------------------------------------------------------

  eavepa(j)             = uentai/SIGN( ( DABS( rentai ) + epsilon ), rentai )

!-----------------------------------------------------------------------
!
!   ||||| Energy transfer rates per baryon by nucleon |||||
!   |||||           emission and absorption           |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon due to neutrino
!   production by PRA
!-----------------------------------------------------------------------

  urttp(j)              = - ravetp(j) * eavetp(j)

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon due to neutrino
!   annihilation by neustino absorption on nucleons
!-----------------------------------------------------------------------

  urtpa(j)              = ravepa(j) * eavepa(j)

!-----------------------------------------------------------------------
!  Net energy transfer rate to matter per baryon due to neustino emission
!   by nucleons
!-----------------------------------------------------------------------

  urtpp(j)              = urttp(j) + urtpa(j)
  dudt_ABEM(j,n)        = urtpp(j)

!-----------------------------------------------------------------------
!
!              ||||| d(entropy)/k-n)/dt by PRA |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  d(entropy)/k-n)/dt by neutrino production due to neutrino emission by
!   nucleons
!-----------------------------------------------------------------------

  dstpdt(j)             = urttp(j)/t_mev

!-----------------------------------------------------------------------
!  d(entropy)/k-n)/dt by neutrino annihilation due to neutrino
!   absorption on nucleons
!-----------------------------------------------------------------------

  dspadt(j)             = urtpa(j)/t_mev

!-----------------------------------------------------------------------
!  Net d(entropy)/k-n)/dt due to PRA
!-----------------------------------------------------------------------

  dstadt(j)             = dstpdt(j) + dspadt(j)

END DO ! j

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

WRITE (nprint,1501)
WRITE (nprint,1503)

DO j = jr_max,jr_min,-idxeng(8,n)
  WRITE (nprint,1505) j,ravetp(j),ravepa(j),eavetp(j),eavepa(j),urttp(j),urtpa(j), &
&  urtpp(j),dstpdt(j),dspadt(j),dstadt(j)
END DO ! j

!-----------------------------------------------------------------------
!
!                  \\\\\ DEALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------


DEALLOCATE (ravetp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravetp    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ravepa, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravepa    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (eavetp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavetp    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (eavepa, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavepa    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (urttp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'urttp     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (urtpa, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'urtpa     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (urtpp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'urtpp     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dstpdt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dstpdt    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dspadt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dspadt    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dstadt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dstadt    '; WRITE (nlog,2001) var_name; END IF


RETURN
END SUBROUTINE editng_N_ABEM
