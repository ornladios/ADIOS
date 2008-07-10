SUBROUTINE editng_NNS( n, jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editng_NNS
!    Module:       editng_NNS
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/05
!
!    Purpose:
!      To edit integral n-type neutrino transport data pertaining to
!       neutrino-nucleon inelastic scattering (NNS).
!
!    Subprograms call:
!  sctedns    : Computes quantities needed for editing NNS rates
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
USE numerical_module, ONLY : zero, half, epsilon
USE physcnst_module, ONLY : cvel, rmu, kmev

USE edit_module, ONLY : prnttest, nprint, nlog, head, idxeng, dudt_NNS
USE mdl_cnfg_module, ONLY : rho, t, ye
USE nu_dist_module, ONLY : ncoefa, ecoefa, psi0, psi1, unu
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

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

REAL(KIND=double)                 :: rmdns         ! neutrino inverse mean free path
REAL(KIND=double)                 :: rmdns0        ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                 :: rmdns1        ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                 :: a0w           ! coefficient for computing change in psi0
REAL(KIND=double)                 :: b0w           ! coefficient for computing change in psi0
REAL(KIND=double)                 :: c0w           ! coefficient for computing change in psi0
REAL(KIND=double)                 :: a1w           ! coefficient for computing change in psi1
REAL(KIND=double)                 :: b1w           ! coefficient for computing change in psi1
REAL(KIND=double)                 :: c1w           ! coefficient for computing change in psi1

REAL(KIND=double)                 :: arscte        ! coefficient for computing NNS scattering rate
REAL(KIND=double)                 :: brscte        ! coefficient for computing NNS scattering rate
REAL(KIND=double)                 :: crscte        ! coefficient for computing NNS scattering rate
REAL(KIND=double)                 :: arsctu        ! coefficient for computing NNS upscattering rate
REAL(KIND=double)                 :: brsctu        ! coefficient for computing NNS upscattering rate
REAL(KIND=double)                 :: crsctu        ! coefficient for computing NNS upscattering rate
REAL(KIND=double)                 :: arscti        ! coefficient for computing NNS isoscattering rate
REAL(KIND=double)                 :: brscti        ! coefficient for computing NNS isoscattering rate
REAL(KIND=double)                 :: crscti        ! coefficient for computing NNS isoscattering rate
REAL(KIND=double)                 :: arsctd        ! coefficient for computing NNS downscattering rate
REAL(KIND=double)                 :: brsctd        ! coefficient for computing NNS downscattering rate
REAL(KIND=double)                 :: crsctd        ! coefficient for computing NNS downscattering rate
REAL(KIND=double)                 :: aesctu        ! coefficient for computing NNS upscattering energy
REAL(KIND=double)                 :: besctu        ! coefficient for computing NNS upscattering energy
REAL(KIND=double)                 :: cesctu        ! coefficient for computing NNS upscattering energy
REAL(KIND=double)                 :: aescti        ! coefficient for computing NNS isoscattering energy
REAL(KIND=double)                 :: bescti        ! coefficient for computing NNS isoscattering energy
REAL(KIND=double)                 :: cescti        ! coefficient for computing NNS isoscattering energy
REAL(KIND=double)                 :: aesctd        ! coefficient for computing NNS downscattering energy
REAL(KIND=double)                 :: besctd        ! coefficient for computing NNS downscattering energy
REAL(KIND=double)                 :: cesctd        ! coefficient for computing NNS downscattering energy
REAL(KIND=double)                 :: arns          ! coefficient for computing NNS rate
REAL(KIND=double)                 :: brns          ! coefficient for computing NNS rate
REAL(KIND=double)                 :: crns          ! coefficient for computing NNS rate
REAL(KIND=double)                 :: aens          ! coefficient for computing NNS fractional eneergy transfer
REAL(KIND=double)                 :: bens          ! coefficient for computing NNS fractional eneergy transfer
REAL(KIND=double)                 :: cens          ! coefficient for computing NNS fractional eneergy transfer

REAL(KIND=double)                 :: t_mev         ! temperature (MeV)

REAL(KIND=double)                 :: dpsire        ! rate of change in psi0 due to NNS inscattering
REAL(KIND=double)                 :: dpsi0e        ! rate of change in psi0 due to NNS
REAL(KIND=double)                 :: rencti        ! rate of neutrinos scattered by NNS (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencti        ! rate of energy chenge due to NNS (MeV g^{-1} s^{-1})

REAL(KIND=double)                 :: dpsi0u        ! rate of change in psi0 due to NNS upscattering
REAL(KIND=double)                 :: dpsieu        ! rate of change in psi0 due to NNS upscattering
REAL(KIND=double)                 :: rencui        ! rate of neutrinos upscattered by NNS (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencui        ! rate of energy chenge due to NNS upscattere (MeV g^{-1} s^{-1})

REAL(KIND=double)                 :: dpsi0i        ! rate of change in psi0 due to NNS isoenergetic scattering
REAL(KIND=double)                 :: dpsiei        ! rate of change in psi0 due to NNS isoenergetic scattering
REAL(KIND=double)                 :: rencii        ! rate of neutrinos isoenergeticly scattered by NNS (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencii        ! rate of energy chenge due to NNS isoenergeticly scattered (MeV g^{-1} s^{-1})

REAL(KIND=double)                 :: dpsi0d        ! rate of change in psi0 due to NNS downscattering
REAL(KIND=double)                 :: dpsied        ! rate of change in psi0 due to NNS downscattering
REAL(KIND=double)                 :: rencdi        ! rate of neutrinos down scattered by NNS (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencdi        ! rate of energy chenge due to NNS down scattered (MeV g^{-1} s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: raveuw   ! upscattering rate of n-neutrinos per baryon by NNS (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: raveiw   ! isoenergetic scattering rate of n-neutrinos per baryon by NNS (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ravedw   ! downscattering rate of n-neutrinos per baryon by NNS (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ravetw   ! total cattering rate of n-neutrinos per baryon by NNS (s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eaveuw   ! mean energy change per event due to upscattering of n-neutrinos by NNS (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eaveiw   ! mean energy change per event due to isoenergetic scattering of n-neutrinos by NNS (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eavedw   ! mean energy change per event due to downscattering of n-neutrinos by NNS (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eavetw   ! mean energy change per event due to scattering of n-neutrinos by NNS (MeV)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: uraveuw  ! mean energy transfer rate by NNS upscattering of n-neutrinos (MeV s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: uravedw  ! mean energy transfer rate by NNS downscattering of n-neutrinos (MeV s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: uravew   ! mean energy transfer rate by NNS scattering of n-neutrinos (MeV s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dsnsutw  ! d(entropy)/k-n)/dt by NNS upscattering  of n-neutrinos (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dsnsdtw  ! d(entropy)/k-n)/dt by NNS downscattering  of n-neutrinos (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dsnstw   ! d(entropy)/k-n)/dt by NNS scattering  of n-neutrinos (s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)  :: psi1_jmh   ! Zone centered psi1

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in editng_NNS')
 1501 FORMAT (30x,'Neutrino-Nucleon Elastic Scattering: Neutrino Number,&
 & Energy, and Matter Entropy Production Rates'/)
 1503 FORMAT ('   j denu upsct   upsct/n   denu isct   isct/n   denu dnsct &
 &  dnsct/n  de-nu nncs  nncs rt/n  unsct u/n  dnsct u/n  nncs u/n &
 & ds NNSu/dt ds NNSd/dt ds NNSt/dt'/)
 1505 FORMAT (1x,i4,14(es11.3))
 2001 FORMAT (' Deallocation problem for array ',a10,' in editng_NNS')

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

ALLOCATE (raveuw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'raveuw    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (raveiw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'raveiw    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ravedw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravedw    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ravetw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravetw    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (eaveuw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eaveuw    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (eaveiw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eaveiw    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (eavedw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavedw    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (eavetw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavetw    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (uraveuw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uraveuw   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uravedw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uravedw   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uravew(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uravew    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dsnsutw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dsnsutw   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dsnsdtw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dsnsdtw   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dsnstw(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dsnstw    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi1_jmh(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1_jmh  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

raveuw                  = zero
raveiw                  = zero
ravedw                  = zero
ravetw                  = zero

eaveuw                  = zero
eaveiw                  = zero
eavedw                  = zero
eavetw                  = zero

uraveuw                 = zero
uravedw                 = zero
uravew                  = zero

dsnsutw                 = zero
dsnsdtw                 = zero
dsnstw                  = zero

!-----------------------------------------------------------------------
!
!                 \\\\\ COMPUTE EDIT QUANTITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Zone centered psi1
!-----------------------------------------------------------------------


DO j = jr_min,jr_max
  DO k = 1,nnugp(n)
    psi1_jmh(j,k,n)     = half * ( psi1(j-1,k,n) + psi1(j,k,n) )
  END DO ! k
END DO ! j

DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

  rencti                = zero
  uencti                = zero
  rencui                = zero
  uencui                = zero
  rencii                = zero
  uencui                = zero
  rencdi                = zero
  uencdi                = zero

  t_mev                 = kmev * t(j)

  DO k = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering info
!-----------------------------------------------------------------------

    CALL sctedns( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), a0w, b0w,    &
&    c0w, a1w, b1w, c1w, rmdns0, rmdns1, rmdns, arns, brns, crns, aens,      &
&    bens, cens, arscte, brscte, crscte, arsctu, brsctu, crsctu, arscti,     &
&    brscti, crscti, arsctd, brsctd, crsctd, aesctu, besctu, cesctu, aescti, &
&    bescti, cescti, aesctd, besctd, cesctd )

!-----------------------------------------------------------------------
!  Total scattering and energy change rate due to NNS
!-----------------------------------------------------------------------

    dpsire              = ( arscte * psi0(j,k,n) + brscte * psi1_jmh(j,k,n) + crscte ) * cvel
    dpsi0e              = ( a0w    * psi0(j,k,n) + b0w    * psi1_jmh(j,k,n) + c0w    ) * cvel
    rencti              = rencti +            ncoefa(j,k) * dpsire/rho(j)
    uencti              = uencti + unu(j,k) * ncoefa(j,k) * dpsi0e/rho(j)

!-----------------------------------------------------------------------
!  Upscattering of n-neutrinos by NNS
!-----------------------------------------------------------------------

    dpsi0u              = ( arsctu * psi0(j,k,n) + brsctu * psi1_jmh(j,k,n) + crsctu ) * cvel
    dpsieu              = ( aesctu * psi0(j,k,n) + besctu * psi1_jmh(j,k,n) + cesctu ) * cvel
    rencui              = rencui +            ncoefa(j,k) * dpsi0u/rho(j)
    uencui              = uencui + unu(j,k) * ncoefa(j,k) * dpsieu/rho(j)

!-----------------------------------------------------------------------
!  Isoenergetic scattering of n-neutrinos by NNS
!-----------------------------------------------------------------------

    dpsi0i              = ( arscti * psi0(j,k,n) + brscti * psi1_jmh(j,k,n) + crscti ) * cvel
    dpsiei              = ( aescti * psi0(j,k,n) + bescti * psi1_jmh(j,k,n) + cescti ) * cvel
    rencii              = rencii +            ncoefa(j,k) * dpsi0i/rho(j)
    uencii              = uencii + unu(j,k) * ncoefa(j,k) * dpsiei/rho(j)

!-----------------------------------------------------------------------
!  Down scattering of n-neutrinos by NNS
!-----------------------------------------------------------------------

    dpsi0d              = ( arsctd * psi0(j,k,n) + brsctd * psi1_jmh(j,k,n) + crsctd ) * cvel
    dpsied              = ( aesctd * psi0(j,k,n) + besctd * psi1_jmh(j,k,n) + cesctd ) * cvel
    rencdi              = rencdi +            ncoefa(j,k) * dpsi0d/rho(j)
    uencdi              = uencdi + unu(j,k) * ncoefa(j,k) * dpsied/rho(j)

  END DO ! k

!-----------------------------------------------------------------------
!
!          ||||| Scattering rates per baryon due to NNS |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Upscattering rate of n-neutrinos per baryon by NNS
!-----------------------------------------------------------------------

  raveuw(j)             = rencui * rmu

!-----------------------------------------------------------------------
!  Isoenergetic scattering rate of n-neutrinos per baryon by NNS
!-----------------------------------------------------------------------

  raveiw(j)             = rencii * rmu

!-----------------------------------------------------------------------
!  Downscattering rate of n-neutrinos per baryon by NNS
!-----------------------------------------------------------------------

  ravedw(j)             = rencdi * rmu

!-----------------------------------------------------------------------
!  Total cattering rate of n-neutrinos per baryon by NNS
!-----------------------------------------------------------------------

  ravetw(j)             = rencti * rmu

!-----------------------------------------------------------------------
!
!             ||||| Mean energy changes due to NNS |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Mean energy change per event due to upscattering of n-neutrinos
!   by NNS
!-----------------------------------------------------------------------

  eaveuw(j)             = uencui/SIGN( ( DABS( rencui ) + epsilon ), rencui )

!-----------------------------------------------------------------------
!  Mean energy change per event due to isoenergetic scattering of
!   n-neutrinos by NNS
!-----------------------------------------------------------------------

  eaveiw(j)             = uencii/SIGN( ( DABS( rencii ) + epsilon ), rencii )

!-----------------------------------------------------------------------
!  Mean energy change per event due to downscattering of n-neutrinos
!    by NNS
!-----------------------------------------------------------------------

  eavedw(j)             = uencdi/SIGN( ( DABS( rencdi ) + epsilon ), rencdi )

!-----------------------------------------------------------------------
!  Mean energy change per event due to scattering of n-neutrinos by NNS
!-----------------------------------------------------------------------

  eavetw(j)             = uencti/SIGN( ( DABS( rencti ) + epsilon ), rencti )

!-----------------------------------------------------------------------
!
!          ||||| Energy transfer rates per baryon by NNS |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon by NNS upscattering
!-----------------------------------------------------------------------

  uraveuw(j)            = - raveuw(j) * eaveuw(j)

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon by NNS downscattering
!-----------------------------------------------------------------------

  uravedw(j)            = - ravedw(j) * eavedw(j)

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon by NNS scattering
!-----------------------------------------------------------------------

  uravew(j)             = - ravetw(j) * eavetw(j)
  dudt_NNS(j,n)         = uravew(j)

!-----------------------------------------------------------------------
!        d(entropy)/k-n)/dt by NNS
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  d(entropy)/k-n)/dt by NNS downscattering
!-----------------------------------------------------------------------

  dsnsutw(j)            = ( -eaveuw(j)/t_mev ) * raveuw(j)

!-----------------------------------------------------------------------
!  d(entropy)/k-n)/dt by NAS downscattering
!-----------------------------------------------------------------------

  dsnsdtw(j)            = ( -eavedw(j)/t_mev ) * ravedw(j)

!-----------------------------------------------------------------------
!  d(entropy)/k-n)/dt by NAS scattering
!-----------------------------------------------------------------------

  dsnstw(j)             = ( -eavetw(j)/t_mev ) * ravetw(j)

END DO ! j

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

WRITE (nprint,1501)
WRITE (nprint,1503)

DO j = jr_max,jr_min,-idxeng(8,n)
  WRITE (nprint,1505) j,eaveuw(j),raveuw(j),eaveiw(j),raveiw(j),eavedw(j),ravedw(j), &
&  eavetw(j),ravetw(j),uraveuw(j),uravedw(j),uravew(j),dsnsutw(j),dsnsdtw(j),dsnstw(j)
END DO ! j

!-----------------------------------------------------------------------
!
!                  \\\\\ DEALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------


DEALLOCATE (raveuw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'raveuw    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (raveiw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'raveiw    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ravedw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravedw    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ravetw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ravetw    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (eaveuw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eaveuw    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (eaveiw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eaveiw    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (eavedw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavedw    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (eavetw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eavetw    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (uraveuw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uraveuw   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (uravedw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uravedw   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (uravew, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uravew    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dsnsutw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dsnsutw   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dsnsdtw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dsnsdtw   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dsnstw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dsnstw    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (psi1_jmh, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1_jmh  '; WRITE (nlog,2001) var_name; END IF


RETURN
END SUBROUTINE editng_NNS
