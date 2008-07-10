SUBROUTINE sctnArate( n, jmin, jmax, ij_ray, ik_ray, rho, t, ye, nx )
!-----------------------------------------------------------------------
!
!    File:         sctnArate
!    Module:       sctnArate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/09/02
!
!    Purpose:
!      To compute the zero and first legendre coefs for the n-type
!       neutrino-nucleus inelastic scattering functions and their
!       derivatives with respect to temperature and electron fraction.
!
!        The scattering functions are included in the multi-group
!         diffusion equations, which have the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2 = xw + yw*psi0 + zw*dpsi0/dr
!
!         where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!        inelastic neutrino scattering contributes to the terms a0w,
!          a1w, b0w, b1w, c0w, c1w as follows:
!
!            a0w  = -K   Int{ w2'dw'[ phi0in(w,w')psi0(w') + phi0out(w,w')( 1 - psi0(w') ) ] }
!            b0w  = -K/3 Int{ w2'dw'[ phi1in(w,w') - phi1out(w,w') ]psi1(w') }
!            c0w  =  K   Int{ w2'dw'[ phi0in(w,w')psi0(w') ] }
!            a1w  = -K   Int{ w2'dw'[ phi1in(w,w') - phi1out(w,w') ]psi1(w') }
!            b1w  = -K   Int{ w2'dw'[ phi0in(w,w')psi0(w') + phi0out(w,w')( 1 - psi0(w') ) ] }
!            c1w  =  K   Int{ w2'dw'[ phi1in(w,w')psi1(w') ] }
!
!        Isoenergetic scattering contributes to the term b1w as follows:
!
!            b1w  = K*w2*[ phi1(w,w) - phi0(w,w) ]
!
!    Input arguments:
!
!  n                   : neutrino type
!  jmin                : inner radial zone number
!  jmax                : outer radial zone number
!  ij_ray              : j-index of a radial ray
!  ik_ray              : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!
!  rho                 : matter density (g/cm**3)
!  t                   : matter temperature (K)
!  ye                  : matter electron fraction
!  iscat               : 0, all neutrino scattering processes omitted 1,
!                         neutrino scattering processes not necessarily omitted
!   isctnA             : 0, neutrino-nucleus inelastic scattering omitted
!                        1, neutrino-nucleus inelastic scattering included
!  rhosctnAemn         : density below which e-neutrino-nucleus inelastic
!                         scattering is turned off.
!  rhosctnAemx         : density above which e-neutrino-nucleus inelastic
!                         scattering is turned off.
!  rhosctnAtmn         : density below which t-neutrino-nucleus inelastic
!                         scattering is turned off.
!  rhosctnAtmx         : density above which t-neutrino-nucleus inelastic
!                         scattering is turned off.
!  nnugp(n)            : number of energy zones for neutrinos of type n
!  unu(k)              : energy of energy zone k (MeV)
!  dunu(k)             : energy width of energy zone k (MeV)
!  psi0(j,k,n)         : zeroth moment of of the neutrino occupation probability
!                         for neutrinos of type n, energy zone k, radial zone j
!  psi1(j,k,n)         : first moment of of the neutrino occupation probability
!                         for neutrinos of type n, energy zone k, radial zone j
!
!    Output arguments (common):
!
!  scnAf(i,j,k,n)      : neutrino-electron scattering function i
!                       i = 1: a0w, i = 2: b0w, i = 3: c0w,
!                       i = 4: a1w, i = 5: b1w, i = 6: c1w,
!
!  scnAfd(i,j,k,n)     : d(scnf(i,j,k,n))/d(density)
!  scnAft(i,j,k,n)     : d(scnf(i,j,k,n))/d(temperature)
!  scnAfy(i,j,k,n)     : d(scnf(i,j,k,n))/d(electron fraction)
!  scnAfp0(i,kp,j,k,n) : d(scnf(i,j,k,n))/d(psi0(j,kp,n))
!  scnAfp1(i,kp,j,k,n) : d(scnf(i,j,k,n))/d(psi1(j,kp,n))
!
!    Subprograms called:
!      sctnAkrnl
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module,
!  scat_n_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nez, nnu
USE numerical_module, ONLY : zero, one
USE physcnst_module

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : iscat, isctnA, rhosctnAemn, rhosctnAemx, rhosctnAtmn, rhosctnAtmx
USE scat_nA_module, ONLY : scnAf, scnAfd, scnAft, scnAfy, scnAfp0, scnAfp1

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: jmin          ! minimum radial zone
INTEGER, INTENT(in)               :: jmax          ! maximum radial zone
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)               :: nx            ! radial array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rho    ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: t      ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: ye     ! electron fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

LOGICAL                           :: sctnA_off

INTEGER                           :: i             ! NNNS scattering function index
INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! incoming neutrino energy zone index
INTEGER                           :: kp            ! outcoming neutrino energy zone index
INTEGER                           :: ke            ! integrated scattering function index

INTEGER                           :: istat         ! allocation status

REAL(KIND=double)                 :: c1            ! integrated scattering function
REAL(KIND=double)                 :: zt            ! integrated scattering function
REAL(KIND=double)                 :: asp           ! integrated scattering function
REAL(KIND=double)                 :: gp            ! integrated scattering function
REAL(KIND=double)                 :: c1d           ! d(c1)/d(rho)
REAL(KIND=double)                 :: ztd           ! d(zt)/d(rho)
REAL(KIND=double)                 :: aspd          ! d(asp)/d(rho)
REAL(KIND=double)                 :: gpd           ! d(gp)/d(rho)
REAL(KIND=double)                 :: c1t           ! d(c1)/d(t)
REAL(KIND=double)                 :: ztt           ! d(zt)/d(t)
REAL(KIND=double)                 :: aspt          ! d(asp)/d(t)
REAL(KIND=double)                 :: gpt           ! d(gp)/d(t)
REAL(KIND=double)                 :: c1y           ! d(c1)/d(ye)
REAL(KIND=double)                 :: zty           ! d(zt)/d(ye)
REAL(KIND=double)                 :: aspy          ! d(asp)/d(ye)
REAL(KIND=double)                 :: gpy           ! d(gp)/d(ye)

REAL(KIND=double)                 :: w2dw          ! w2 * neutrino energy bin width
REAL(KIND=double)                 :: x0            ! psi0
REAL(KIND=double)                 :: x1            ! psi1

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: c1p0          ! d(c1)/d(psi0(k))
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ztp0          ! d(zt)/d(psi0(k))
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: aspp1         ! d(asp)/d(psi1(k))
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: gpp1          ! d(gp)/d(psi1(k))

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0_ot     ! zero moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1_ot     ! first moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0d_ot    ! d(scant0)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1d_ot    ! d(scant1)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0t_ot    ! d(scant0)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1t_ot    ! d(scant1)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0y_ot    ! d(scant0)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1y_ot    ! d(scant1)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0_in     ! zero moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1_in     ! first moment of the NNS scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0d_in    ! d(scant0)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1d_in    ! d(scant1)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0t_in    ! d(scant0)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1t_in    ! d(scant1)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn0y_in    ! d(scant0)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scatn1y_in    ! d(scant1)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0in          ! zero moment of the NNS in scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ind         ! d(f0in)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0int         ! d(f0in)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0iny         ! d(f0in)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ot          ! zero moment of the NNS out scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0otd         ! d(fotn)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ott         ! d(fotn)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0oty         ! d(fotn)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1in          ! first moment of the NNS in scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ind         ! d(f0in)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1int         ! d(f0in)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1iny         ! d(f0in)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ot          ! first moment of the NNS out scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1otd         ! d(fotn)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ott         ! d(fotn)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1oty         ! d(fotn)/d(ye)

EXTERNAL fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in sctnArate')
 2001 FORMAT (' Deallocation problem for array ',a10,' in sctnArate')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set neutrino-nucleus inelastic scattering functions to zero if
!     nnugp(n) =  0
!   or
!     isctnA = 0
!   or
!     iscat = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0  .or.  iscat  == 0  .or.  isctnA == 0 ) RETURN

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (c1p0(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c1p0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ztp0(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ztp0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aspp1(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aspp1     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gpp1(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gpp1      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scatn0_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_ot '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_ot '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0d_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1d_ot(nez), STAT = istat) 
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0t_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1t_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0y_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1y_ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_ot'; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scatn0_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_in '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_in '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0d_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1d_in(nez), STAT = istat) 
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0t_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1t_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn0y_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scatn1y_in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_in'; WRITE (nlog,1001) var_name; END IF

ALLOCATE (f0in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0in      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0ind(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ind     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0int(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0int     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0iny(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0iny     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (f0ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ot      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0otd(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0otd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0ott(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ott     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0oty(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0oty     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (f1in(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1in      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1ind(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ind     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1int(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1int     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1iny(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1iny     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (f1ot(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ot      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1otd(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1otd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1ott(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ott     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1oty(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1oty     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

c1p0                         = zero
ztp0                         = zero
aspp1                        = zero
gpp1                         = zero

scatn0_ot                    = zero
scatn1_ot                    = zero
scatn0d_ot                   = zero
scatn1d_ot                   = zero
scatn0t_ot                   = zero
scatn1t_ot                   = zero
scatn0y_ot                   = zero
scatn1y_ot                   = zero

scatn0_in                    = zero
scatn1_in                    = zero
scatn0d_in                   = zero
scatn1d_in                   = zero
scatn0t_in                   = zero
scatn1t_in                   = zero
scatn0y_in                   = zero
scatn1y_in                   = zero

f0in                         = zero
f0ind                        = zero
f0int                        = zero
f0iny                        = zero

f0ot                         = zero
f0otd                        = zero
f0ott                        = zero
f0oty                        = zero

f1in                         = zero
f1ind                        = zero
f1int                        = zero
f1iny                        = zero

f1ot                         = zero
f1otd                        = zero
f1ott                        = zero
f1oty                        = zero

!-----------------------------------------------------------------------
!                      ||||| Loop over j |||||
!-----------------------------------------------------------------------

DO j = jmin,jmax

!-----------------------------------------------------------------------        
!  Set neutrino-nucleus inelastic scattering functions to zero if rho
!   outside specified boundaries.
!-----------------------------------------------------------------------

  sctnA_off                       = .false.
  IF ( n < 3   .and.  rho(j) < rhosctnAemn ) sctnA_off = .true.
  IF ( n < 3   .and.  rho(j) > rhosctnAemx ) sctnA_off = .true.
  IF ( n == 3  .and.  rho(j) < rhosctnAtmn ) sctnA_off = .true.
  IF ( n == 3  .and.  rho(j) > rhosctnAtmx ) sctnA_off = .true.

  IF ( sctnA_off ) THEN
    DO k = 1,nez
      DO i = 1,6
        scnAf(i,j,k,n)            = zero
        scnAfd(i,j,k,n)           = zero
        scnAft(i,j,k,n)           = zero
        scnAfy(i,j,k,n)           = zero
        DO kp = 1,nez
          scnAfp0(i,kp,j,k,n)     = zero
          scnAfp1(i,kp,j,k,n)     = zero
        END DO
      END DO
    END DO
    CYCLE
  END IF ! sctnA_off

!-----------------------------------------------------------------------
!                      ||||| Loop over k |||||
!-----------------------------------------------------------------------

  DO k = 1,nnugp(n)

!-----------------------------------------------------------------------        
!  NAS neutrino scattering functions
!-----------------------------------------------------------------------

    c1             = zero                                                        
    zt             = zero
    asp            = zero
    gp             = zero
    c1d            = zero                                                        
    ztd            = zero
    aspd           = zero
    gpd            = zero
    c1t            = zero                                                        
    ztt            = zero
    aspt           = zero
    gpt            = zero
    c1y            = zero                                                        
    zty            = zero
    aspy           = zero
    gpy            = zero

    DO ke = 1,nnugp(n)
      c1p0(ke)     = zero
      ztp0(ke)     = zero
      aspp1(ke)    = zero
      gpp1(ke)     = zero
    END DO

    CALL sctnAkrnl( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), scatn0_ot, &
    & scatn1_ot, scatn0d_ot, scatn1d_ot, scatn0t_ot, scatn1t_ot, scatn0y_ot, &
    & scatn1y_ot, scatn0_in, scatn1_in, scatn0d_in, scatn1d_in, scatn0t_in, &
    & scatn1t_in, scatn0y_in, scatn1y_in, nez )

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      f0ot(kp)     = scatn0_ot (kp)
      f0in(kp)     = scatn0_in (kp)
      f0otd(kp)    = scatn0d_ot(kp)
      f0ind(kp)    = scatn0d_in(kp)
      f0ott(kp)    = scatn0t_ot(kp)
      f0int(kp)    = scatn0t_in(kp)
      f0oty(kp)    = scatn0y_ot(kp)
      f0iny(kp)    = scatn0y_in(kp)
      f1ot(kp)     = scatn1_ot (kp)
      f1in(kp)     = scatn1_in (kp)
      f1otd(kp)    = scatn1d_ot(kp)
      f1ind(kp)    = scatn1d_in(kp)
      f1ott(kp)    = scatn1t_ot(kp)
      f1int(kp)    = scatn1t_in(kp)
      f1oty(kp)    = scatn1y_ot(kp)
      f1iny(kp)    = scatn1y_in(kp)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      w2dw         = 2.d0 * pi * unu(j,kp) * unu(j,kp) * dunu(j,kp)
      x0           = psi0(j,kp,n)
      x1           = psi1(j,kp,n)

      zt           = zt   + ( f0in(kp) * x0 + f0ot(kp)  * ( one - x0 ) )  * w2dw
      asp          = asp  + ( f1in(kp)  - f1ot(kp)  ) * x1                * w2dw
      c1           = c1   +   f0in(kp)                * x0                * w2dw
      gp           = gp   +   f1in(kp)                * x1                * w2dw

      ztd          = ztd  + ( f0ind(kp) * x0 + f0otd(kp) * ( one - x0 ) ) * w2dw
      aspd         = aspd + ( f1ind(kp) - f1otd(kp) ) * x1                * w2dw
      c1d          = c1d  +   f0ind(kp)               * x0                * w2dw
      gpd          = gpd  +   f1ind(kp)               * x1                * w2dw

      ztt          = ztt  + ( f0int(kp) * x0 + f0ott(kp) * ( one - x0 ) ) * w2dw
      aspt         = aspt + ( f1int(kp) - f1ott(kp) ) * x1                * w2dw
      c1t          = c1t  +   f0int(kp)               * x0                * w2dw
      gpt          = gpt  +   f1int(kp)               * x1                * w2dw

      zty          = zty  + ( f0iny(kp) * x0 + f0oty(kp) * ( one - x0 ) ) * w2dw
      aspy         = aspy + ( f1iny(kp) - f1oty(kp) ) * x1                * w2dw
      c1y          = c1y  +   f0iny(kp)               * x0                * w2dw
      gpy          = gpy  +   f1iny(kp)               * x1                * w2dw

      ztp0(kp)     = ( f0in(kp) - f0ot(kp) )                              * w2dw
      c1p0(kp)     =   f0in(kp)                                           * w2dw
      aspp1(kp)    = ( f1in(kp) - f1ot(kp) )                              * w2dw
      gpp1(kp)     =   f1in(kp)                                           * w2dw

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO

!-----------------------------------------------------------------------
!  Compute a0w, b0w, c0w, a1w, b1w, c1w, and derivatives.
!
!  scnAf(1) = a0w   scnAfi(1) = da0w/di   scnAfpi(1,k) = da0w/dpsii(k)
!  scnAf(2) = b0w   scnAfi(2) = db0w/di   scnAfpi(2,k) = db0w/dpsii(k)
!  scnAf(3) = c0w   scnAfi(3) = dc0w/di   scnAfpi(3,k) = dc0w/dpsii(k)
!  scnAf(4) = a1w   scnAfi(4) = da1w/di   scnAfpi(4,k) = da1w/dpsii(k)
!  scnAf(5) = b1w   scnAfi(5) = db1w/di   scnAfpi(5,k) = db1w/dpsii(k)
!  scnAf(6) = c1w   scnAfi(6) = dc1w/di   scnAfpi(6,k) = dc1w/dpsii(k)
!-----------------------------------------------------------------------

    scnAf(1,j,k,n)  = -zt
    scnAf(2,j,k,n)  = -asp/3.d0
    scnAf(3,j,k,n)  = c1
    scnAf(4,j,k,n)  = -asp
    scnAf(5,j,k,n)  = -zt
    scnAf(6,j,k,n)  = gp

    scnAfd(1,j,k,n) = -ztd
    scnAfd(2,j,k,n) = -aspd/3.d0
    scnAfd(3,j,k,n) = c1d
    scnAfd(4,j,k,n) = -aspd
    scnAfd(5,j,k,n) = -ztd
    scnAfd(6,j,k,n) = gpd

    scnAft(1,j,k,n) = -ztt
    scnAft(2,j,k,n) = -aspt/3.d0
    scnAft(3,j,k,n) = c1t
    scnAft(4,j,k,n) = -aspt
    scnAft(5,j,k,n) = -ztt
    scnAft(6,j,k,n) = gpt

    scnAfy(1,j,k,n) = -zty
    scnAfy(2,j,k,n) = -aspy/3.d0
    scnAfy(3,j,k,n) = c1y
    scnAfy(4,j,k,n) = -aspy
    scnAfy(5,j,k,n) = -zty
    scnAfy(6,j,k,n) = gpy

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      scnAfp0(1,kp,j,k,n) = -ztp0(kp)
      scnAfp0(2,kp,j,k,n) = zero
      scnAfp0(3,kp,j,k,n) = c1p0(kp)
      scnAfp0(4,kp,j,k,n) = zero
      scnAfp0(5,kp,j,k,n) = -ztp0(kp)
      scnAfp0(6,kp,j,k,n) = zero

      scnAfp1(1,kp,j,k,n) = zero
      scnAfp1(2,kp,j,k,n) = -aspp1(kp)/3.d0
      scnAfp1(3,kp,j,k,n) = zero
      scnAfp1(4,kp,j,k,n) = -aspp1(kp)
      scnAfp1(5,kp,j,k,n) = zero
      scnAfp1(6,kp,j,k,n) = gpp1(kp)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO

!-----------------------------------------------------------------------
!                 ||||| End loop over j and k |||||
!-----------------------------------------------------------------------

  END DO
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (c1p0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c1p0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ztp0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ztp0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (aspp1, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aspp1     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (gpp1, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gpp1      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (scatn0_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_ot '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_ot '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0d_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1d_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0t_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1t_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0y_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1y_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_ot'; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (scatn0_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0_in '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1_in '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0d_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0d_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1d_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1d_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0t_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0t_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1t_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1t_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn0y_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn0y_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scatn1y_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scatn1y_in'; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (f0in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0in      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0ind, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ind     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0int, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0int     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0iny, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0iny     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (f0ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ot      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0otd, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0otd     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0ott, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ott     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0oty, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0oty     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (f1in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1in      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1ind, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ind     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1int, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1int     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1iny, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1iny     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (f1ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ot      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1otd, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1otd     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1ott, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ott     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1oty, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1oty     '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE sctnArate
