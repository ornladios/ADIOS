SUBROUTINE editng_NES( n, jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editng_NES
!    Module:       editng_NES
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/05
!
!    Purpose:
!      To edit integral n-type neutrino transport data pertaining to
!       neutrino-electron scattering (NES).
!
!    Subprograms call:
!  sctednes   : Computes quantities needed for editing NES rates
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

USE edit_module, ONLY : prnttest, nprint, nlog, head, idxeng, dudt_NES
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

REAL(KIND=double)                 :: rmdnes        ! neutrino inverse mean free path
REAL(KIND=double)                 :: rmdnes0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                 :: rmdnes1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                 :: a0w           ! coefficient for computing change in psi0
REAL(KIND=double)                 :: b0w           ! coefficient for computing change in psi0
REAL(KIND=double)                 :: c0w           ! coefficient for computing change in psi0
REAL(KIND=double)                 :: a1w           ! coefficient for computing change in psi1
REAL(KIND=double)                 :: b1w           ! coefficient for computing change in psi1
REAL(KIND=double)                 :: c1w           ! coefficient for computing change in psi1

REAL(KIND=double)                 :: arscte        ! coefficient for computing NES scattering rate
REAL(KIND=double)                 :: brscte        ! coefficient for computing NES scattering rate
REAL(KIND=double)                 :: crscte        ! coefficient for computing NES scattering rate
REAL(KIND=double)                 :: arsctu        ! coefficient for computing NES upscattering rate
REAL(KIND=double)                 :: brsctu        ! coefficient for computing NES upscattering rate
REAL(KIND=double)                 :: crsctu        ! coefficient for computing NES upscattering rate
REAL(KIND=double)                 :: arscti        ! coefficient for computing NES isoscattering rate
REAL(KIND=double)                 :: brscti        ! coefficient for computing NES isoscattering rate
REAL(KIND=double)                 :: crscti        ! coefficient for computing NES isoscattering rate
REAL(KIND=double)                 :: arsctd        ! coefficient for computing NES downscattering rate
REAL(KIND=double)                 :: brsctd        ! coefficient for computing NES downscattering rate
REAL(KIND=double)                 :: crsctd        ! coefficient for computing NES downscattering rate
REAL(KIND=double)                 :: aesctu        ! coefficient for computing NES upscattering energy
REAL(KIND=double)                 :: besctu        ! coefficient for computing NES upscattering energy
REAL(KIND=double)                 :: cesctu        ! coefficient for computing NES upscattering energy
REAL(KIND=double)                 :: aescti        ! coefficient for computing NES isoscattering energy
REAL(KIND=double)                 :: bescti        ! coefficient for computing NES isoscattering energy
REAL(KIND=double)                 :: cescti        ! coefficient for computing NES isoscattering energy
REAL(KIND=double)                 :: aesctd        ! coefficient for computing NES downscattering energy
REAL(KIND=double)                 :: besctd        ! coefficient for computing NES downscattering energy
REAL(KIND=double)                 :: cesctd        ! coefficient for computing NES downscattering energy
REAL(KIND=double)                 :: arnes         ! coefficient for computing NES rate
REAL(KIND=double)                 :: brnes         ! coefficient for computing NES rate
REAL(KIND=double)                 :: crnes         ! coefficient for computing NES rate
REAL(KIND=double)                 :: aenes         ! coefficient for computing NES fractional eneergy transfer
REAL(KIND=double)                 :: benes         ! coefficient for computing NES fractional eneergy transfer
REAL(KIND=double)                 :: cenes         ! coefficient for computing NES fractional eneergy transfer

REAL(KIND=double)                 :: t_mev         ! temperature (MeV)

REAL(KIND=double)                 :: dpsire        ! rate of change in psi0 due to NES inscattering
REAL(KIND=double)                 :: dpsi0e        ! rate of change in psi0 due to NES
REAL(KIND=double)                 :: rencti        ! rate of neutrinos scattered by NES (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencti        ! rate of energy chenge due to NES (MeV g^{-1} s^{-1})

REAL(KIND=double)                 :: dpsi0u        ! rate of change in psi0 due to NES upscattering
REAL(KIND=double)                 :: dpsieu        ! rate of change in psi0 due to NES upscattering
REAL(KIND=double)                 :: rencui        ! rate of neutrinos upscattered by NES (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencui        ! rate of energy chenge due to NES upscattere (MeV g^{-1} s^{-1})

REAL(KIND=double)                 :: dpsi0i        ! rate of change in psi0 due to NES isoenergetic scattering
REAL(KIND=double)                 :: dpsiei        ! rate of change in psi0 due to NES isoenergetic scattering
REAL(KIND=double)                 :: rencii        ! rate of neutrinos isoenergeticly scattered by NES (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencii        ! rate of energy chenge due to NES isoenergeticly scattered (MeV g^{-1} s^{-1})

REAL(KIND=double)                 :: dpsi0d        ! rate of change in psi0 due to NES downscattering
REAL(KIND=double)                 :: dpsied        ! rate of change in psi0 due to NES downscattering
REAL(KIND=double)                 :: rencdi        ! rate of neutrinos down scattered by NES (g^{-1} s^{-1})
REAL(KIND=double)                 :: uencdi        ! rate of energy chenge due to NES down scattered (MeV g^{-1} s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: raveuw   ! upscattering rate of n-neutrinos per baryon by NES (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: raveiw   ! isoenergetic scattering rate of n-neutrinos per baryon by NES (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ravedw   ! downscattering rate of n-neutrinos per baryon by NES (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ravetw   ! total cattering rate of n-neutrinos per baryon by NES (s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eaveuw   ! mean energy change per event due to upscattering of n-neutrinos by NES (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eaveiw   ! mean energy change per event due to isoenergetic scattering of n-neutrinos by NES (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eavedw   ! mean energy change per event due to downscattering of n-neutrinos by NES (MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eavetw   ! mean energy change per event due to scattering of n-neutrinos by NES (MeV)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: uraveuw  ! mean energy transfer rate by NES upscattering of n-neutrinos (MeV s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: uravedw  ! mean energy transfer rate by NES downscattering of n-neutrinos (MeV s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: uravew   ! mean energy transfer rate by NES scattering of n-neutrinos (MeV s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dsnsutw  ! d(entropy)/k-n)/dt by NES upscattering  of n-neutrinos (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dsnsdtw  ! d(entropy)/k-n)/dt by NES downscattering  of n-neutrinos (s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dsnstw   ! d(entropy)/k-n)/dt by NES scattering  of n-neutrinos (s^{-1})

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)  :: psi1_jmh   ! Zone centered psi1

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in editng_NES')
 1501 FORMAT (30x,'Neutrino-Electron Scattering: Neutrino Number, Energy,&
 & and Matter Entropy Production Rates'/)
 1503 FORMAT ('   j denu upsct   upsct/n   denu isct   isct/n   denu dnsct &
 &  dnsct/n   de-nu nes  nes rt/n   unsct u/n  dnsct u/n   nes un &
 &  ds NESu/dt ds NESd/dt ds NESt/dt'/)
 1505 FORMAT (1x,i4,14(es11.3))
 2001 FORMAT (' Deallocation problem for array ',a10,' in editng_NES')

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

psi1_jmh                = zero

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
!  Neutrino-electron elastic scattering info
!-----------------------------------------------------------------------

    CALL sctednes( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), a0w, b0w,    &
&    c0w, a1w, b1w, c1w, rmdnes0, rmdnes1, rmdnes, arnes, brnes, crnes,       &
&    aenes, benes, cenes, arscte, brscte, crscte, arsctu, brsctu, crsctu,     &
&    arscti, brscti, crscti,  arsctd, brsctd, crsctd, aesctu, besctu, cesctu, &
&    aescti, bescti, cescti, aesctd, besctd, cesctd )

!-----------------------------------------------------------------------
!  Total scattering and energy change rate due to NES
!-----------------------------------------------------------------------

    dpsire              = ( arscte * psi0(j,k,n) + brscte * psi1_jmh(j,k,n) + crscte ) * cvel
    dpsi0e              = ( a0w    * psi0(j,k,n) + b0w    * psi1_jmh(j,k,n) + c0w    ) * cvel
    rencti              = rencti +            ncoefa(j,k) * dpsire/rho(j)
    uencti              = uencti + unu(j,k) * ncoefa(j,k) * dpsi0e/rho(j)

!-----------------------------------------------------------------------
!  Upscattering of n-neutrinos by NES
!-----------------------------------------------------------------------

    dpsi0u              = ( arsctu * psi0(j,k,n) + brsctu * psi1_jmh(j,k,n) + crsctu ) * cvel
    dpsieu              = ( aesctu * psi0(j,k,n) + besctu * psi1_jmh(j,k,n) + cesctu ) * cvel
    rencui              = rencui +            ncoefa(j,k) * dpsi0u/rho(j)
    uencui              = uencui + unu(j,k) * ncoefa(j,k) * dpsieu/rho(j)

!-----------------------------------------------------------------------
!  Isoenergetic scattering of n-neutrinos by NES
!-----------------------------------------------------------------------

    dpsi0i              = ( arscti * psi0(j,k,n) + brscti * psi1_jmh(j,k,n) + crscti ) * cvel
    dpsiei              = ( aescti * psi0(j,k,n) + bescti * psi1_jmh(j,k,n) + cescti ) * cvel
    rencii              = rencii +            ncoefa(j,k) * dpsi0i/rho(j)
    uencii              = uencii + unu(j,k) * ncoefa(j,k) * dpsiei/rho(j)

!-----------------------------------------------------------------------
!  Down scattering of n-neutrinos by NES
!-----------------------------------------------------------------------

    dpsi0d              = ( arsctd * psi0(j,k,n) + brsctd * psi1_jmh(j,k,n) + crsctd ) * cvel
    dpsied              = ( aesctd * psi0(j,k,n) + besctd * psi1_jmh(j,k,n) + cesctd ) * cvel
    rencdi              = rencdi +            ncoefa(j,k) * dpsi0d/rho(j)
    uencdi              = uencdi + unu(j,k) * ncoefa(j,k) * dpsied/rho(j)

  END DO ! k

!-----------------------------------------------------------------------
!
!          ||||| Scattering rates per baryon due to NES |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Upscattering rate of n-neutrinos per baryon by NES
!-----------------------------------------------------------------------

  raveuw(j)             = rencui * rmu

!-----------------------------------------------------------------------
!  Isoenergetic scattering rate of n-neutrinos per baryon by NES
!-----------------------------------------------------------------------

  raveiw(j)             = rencii * rmu

!-----------------------------------------------------------------------
!  Downscattering rate of n-neutrinos per baryon by NES
!-----------------------------------------------------------------------

  ravedw(j)             = rencdi * rmu

!-----------------------------------------------------------------------
!  Total cattering rate of n-neutrinos per baryon by NES
!-----------------------------------------------------------------------

  ravetw(j)             = rencti * rmu

!-----------------------------------------------------------------------
!
!             ||||| Mean energy changes due to NES |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Mean energy change per event due to upscattering of n-neutrinos
!   by NES
!-----------------------------------------------------------------------

  eaveuw(j)             = uencui/SIGN( ( DABS( rencui ) + epsilon ), rencui )

!-----------------------------------------------------------------------
!  Mean energy change per event due to isoenergetic scattering of
!   n-neutrinos by NES
!-----------------------------------------------------------------------

  eaveiw(j)             = uencii/SIGN( ( DABS( rencii ) + epsilon ), rencii )

!-----------------------------------------------------------------------
!  Mean energy change per event due to downscattering of n-neutrinos
!    by NES
!-----------------------------------------------------------------------

  eavedw(j)             = uencdi/SIGN( ( DABS( rencdi ) + epsilon ), rencdi )

!-----------------------------------------------------------------------
!  Mean energy change per event due to scattering of n-neutrinos by NES
!-----------------------------------------------------------------------

  eavetw(j)             = uencti/SIGN( ( DABS( rencti ) + epsilon ), rencti )

!-----------------------------------------------------------------------
!
!          ||||| Energy transfer rates per baryon by NES |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon by NES upscattering
!-----------------------------------------------------------------------

  uraveuw(j)            = - raveuw(j) * eaveuw(j)

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon by NES downscattering
!-----------------------------------------------------------------------

  uravedw(j)            = - ravedw(j) * eavedw(j)

!-----------------------------------------------------------------------
!  Mean energy transfer rate to matter per baryon by NES scattering
!-----------------------------------------------------------------------

  uravew(j)             = - ravetw(j) * eavetw(j)
  dudt_NES(j,n)         = uravew(j)

!-----------------------------------------------------------------------
!
!              ||||| d(entropy)/k-n)/dt by NES |||||
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  d(entropy)/k-n)/dt by NES downscattering
!-----------------------------------------------------------------------

  dsnsutw(j)            = ( -eaveuw(j)/t_mev ) * raveuw(j)

!-----------------------------------------------------------------------
!  d(entropy)/k-n)/dt by NES downscattering
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
END SUBROUTINE editng_NES
