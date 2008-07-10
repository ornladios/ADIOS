SUBROUTINE scterate( n , jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )
!-----------------------------------------------------------------------
!
!    File:         scterate
!    Module:       scterate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/01/03
!
!    Purpose:
!      To compute the zero and first legendre coefs for the n-type neutrino-electron 
!       scattering functions and their derivatives with respect to temperature
!       and electron fraction are computed here.
!
!        The scattering functions are included in the multi-group
!          diffusion equations, which have the form
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
!        Inelastic neutrino scattering contributes to the terms a0w,
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
!  n           : neutrino type
!  jr_min      : inner radial zone number
!  jr_max      : outer radial zone number
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!
!  rho         : matter density (g/cm**3)
!  t           : matter temperature (K)
!  ye          : matter electron fraction
!  iscat       : 0, all neutrino scattering processes omitted
!                1, neutrino scattering processes not necessarily omitted
!  nes         : 0, neutrino-electron scattering omitted
!                1, neutrino-electron scattering included
!  rhonesmn    : rho < rhonesmn: neutrino-electron scattering omitted
!  rhonesmx    : rho > rhonesmx: neutrino-electron scattering omitted
!  nnugp(n)    : number of energy zones for neutrinos of type n
!  unu(k)      : energy of energy zone k (MeV)
!  dunu(k)     : energy width of energy zone k (MeV)
!  psi0(j,k,n) : zeroth moment of of the neutrino occupation probability for neutrinos 
!                of type n, energy zone k, radial zone j
!  psi1(j,k,n) : first moment of of the neutrino occupation probability for neutrinos
!                 of type n, energy zone k, radial zone j
!
!    Output arguments (common):
!
!  scef(i,j,k,n)
!              : neutrino-electron scattering function i
!                i = 1: a0w, i = 2: b0w, i = 3: c0w,
!                i = 4: a1w, i = 5: b1w, i = 6: c1w,
!  scefd(i,j,k,n)
!              : d(scef(i,j,k,n))/d(density)
!  sceft(i,j,k,n)
!              : d(scef(i,j,k,n))/d(temperature)
!  scefy(i,j,k,n)
!              : d(scef(i,j,k,n))/d(electron fraction)
!  scefp0(i,kp,j,k,n)
!              : d(scef(i,j,k,n))/d(psi0(j,kp,n))
!  scefp1(i,kp,j,k,n)
!              : d(scef(i,j,k,n))/d(psi1(j,kp,n))
!
!    Subprograms called:
!      sctekrnl
!
!    Include files:
!  array_module, kind_module, numerical_module, physcnst_module
!  edit_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module,
!  scat_e_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nez, nnu
USE numerical_module, ONLY : zero, one
USE physcnst_module

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : nes, iscat, rhonesmn, rhonesmx
USE scat_e_module, ONLY : scef, scefd, sceft, scefy, scefp0, scefp1

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: jr_min        ! minimum radial zone
INTEGER, INTENT(in)               :: jr_max        ! maximum radial zone
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

LOGICAL                           :: scte_off

INTEGER                           :: i             ! NES scattering function index
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

REAL(KIND=double)                 :: enuin         ! incoming neutrino energy/kt
REAL(KIND=double)                 :: enuout        ! outgoing neutrino energy/kt
REAL(KIND=double)                 :: tmev          ! temperature (MeV)

REAL(KIND=double)                 :: arg           ! argument of an exponential
REAL(KIND=double)                 :: rinv          ! factor relating in and out scattering

REAL(KIND=double)                 :: w2            ! square of the incoming neutrino energy
REAL(KIND=double)                 :: w2dw          ! w2 * neutrino energy bin width
REAL(KIND=double)                 :: x0            ! psi0
REAL(KIND=double)                 :: x1            ! psi1

REAL(KIND=double)                 :: fexp          ! exponential function

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: c1p0          ! d(c1)/d(psi0(k))
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ztp0          ! d(zt)/d(psi0(k))
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: aspp1         ! d(asp)/d(psi1(k))
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: gpp1          ! d(gp)/d(psi1(k))

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate0        ! zero moment of the NES scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate1        ! first moment of the NES scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate0d       ! d(scat0)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate1d       ! d(scat1)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate0t       ! d(scat0)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate1t       ! d(scat1)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate0y       ! d(scat0)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: scate1y       ! d(scat1)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0in          ! zero moment of the NES in scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ind         ! d(f0in)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0int         ! d(f0in)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0iny         ! d(f0in)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ot          ! zero moment of the NES out scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0otd         ! d(fotn)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ott         ! d(fotn)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0oty         ! d(fotn)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1in          ! first moment of the NES in scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ind         ! d(f0in)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1int         ! d(f0in)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1iny         ! d(f0in)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ot          ! first moment of the NES out scattering function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1otd         ! d(fotn)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ott         ! d(fotn)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1oty         ! d(fotn)/d(ye)

EXTERNAL fexp

 1001 FORMAT (' Allocation problem for array ',a10,' in scterate')
 2001 FORMAT (' Deallocation problem for array ',a10,' in scterate')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------        
!  Set neutrino-electron scattering functions to zero if
!     nnugp(n) =  0
!   or
!     nes = 0
!   or
!     iscat = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0  .or.  nes == 0  .or.  iscat == 0 ) RETURN

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

ALLOCATE (scate0(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scate1(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scate0d(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0d   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scate1d(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1d   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scate0t(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0t   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scate1t(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1t   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scate0y(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0y   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scate1y(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1y   '; WRITE (nlog,1001) var_name; END IF

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

scate0                       = zero
scate1                       = zero
scate0d                      = zero
scate1d                      = zero
scate0t                      = zero
scate1t                      = zero
scate0y                      = zero
scate1y                      = zero

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

outer: DO j = jr_min,jr_max

!-----------------------------------------------------------------------        
!  Set neutrino-electron functions to zero if rho outside specified
!   boundaries.
!-----------------------------------------------------------------------

  scte_off           = .false.
  IF ( rho(j) < rhonesmn ) scte_off = .true.
  IF ( rho(j) > rhonesmx ) scte_off = .true.
  IF ( scte_off ) THEN
    DO k = 1,nez
      DO i = 1,6
        scef(i,j,k,n)        = zero
        scefd(i,j,k,n)       = zero
        sceft(i,j,k,n)       = zero
        scefy(i,j,k,n)       = zero
        DO kp = 1,nez
          scefp0(i,kp,j,k,n) = zero
          scefp1(i,kp,j,k,n) = zero
        END DO ! kp
      END DO ! i
    END DO ! k
    CYCLE ! outer
  END IF ! scte_off

!-----------------------------------------------------------------------
!                      ||||| Loop over k |||||
!-----------------------------------------------------------------------

  DO k = 1,nnugp(n)

    tmev           = t(j) * kmev
    enuin          = unu(j,k)/tmev

!-----------------------------------------------------------------------        
!  NES neutrino scattering functions
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
    END DO ! ke

    CALL sctekrnl( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j),  scate0, &
&    scate1, scate0d, scate1d, scate0t, scate1t, scate0y, scate1y, nez )

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      enuout       = unu(j,kp)/tmev

      IF ( kp < k ) THEN

        arg        = enuout - enuin
        rinv       = fexp( arg )

        f0ot(kp)   = scate0(kp)
        f0in(kp)   = scate0(kp)                                  * rinv
        f0otd(kp)  = scate0d(kp)
        f0ind(kp)  = scate0d(kp)                                 * rinv
        f0ott(kp)  = scate0t(kp)
        f0int(kp)  = scate0t(kp)                                 * rinv
        f0oty(kp)  = scate0y(kp)
        f0iny(kp)  = scate0y(kp)                                 * rinv
        f1ot(kp)   = scate1(kp)
        f1in(kp)   = scate1(kp)                                  * rinv
        f1otd(kp)  = scate1d(kp)
        f1ind(kp)  = scate1d(kp)                                 * rinv
        f1ott(kp)  = scate1t(kp)
        f1int(kp)  = scate1t(kp)                                 * rinv
        f1oty(kp)  = scate1y(kp)
        f1iny(kp)  = scate1y(kp)                                 * rinv

      ELSE IF ( kp .eq. k ) THEN

        f0ot(kp)   = scate0 (kp)
        f0in(kp)   = scate0 (kp)
        f0otd(kp)  = scate0d(kp)
        f0ind(kp)  = scate0d(kp)
        f0ott(kp)  = scate0t(kp)
        f0int(kp)  = scate0t(kp)
        f0oty(kp)  = scate0y(kp)
        f0iny(kp)  = scate0y(kp)
        f1ot(kp)   = scate1 (kp)
        f1in(kp)   = scate1 (kp)
        f1otd(kp)  = scate1d(kp)
        f1ind(kp)  = scate1d(kp)
        f1ott(kp)  = scate1t(kp)
        f1int(kp)  = scate1t(kp)
        f1oty(kp)  = scate1y(kp)
        f1iny(kp)  = scate1y(kp)

      ELSE IF ( kp .gt. k ) THEN

        arg        = enuin - enuout
        rinv       = fexp( arg )

        f0in(kp)   = scate0(kp)
        f0ot(kp)   = scate0(kp)                                  * rinv
        f0ind(kp)  = scate0d(kp)
        f0otd(kp)  = scate0d(kp)                                 * rinv
        f0int(kp)  = scate0t(kp)
        f0ott(kp)  = scate0t(kp)                                 * rinv
        f0iny(kp)  = scate0y(kp)
        f0oty(kp)  = scate0y(kp)                                 * rinv
        f1in(kp)   = scate1(kp)
        f1ot(kp)   = scate1(kp)                                  * rinv
        f1ind(kp)  = scate1d(kp)
        f1otd(kp)  = scate1d(kp)                                 * rinv
        f1int(kp)  = scate1t(kp)
        f1ott(kp)  = scate1t(kp)                                 * rinv
        f1iny(kp)  = scate1y(kp)
        f1oty(kp)  = scate1y(kp)                                 * rinv

      END IF

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO ! kp

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      w2           = unu(j,kp) * unu(j,kp)
      w2dw         = w2 * dunu(j,kp)
      x0           = psi0(j,kp,n)
      x1           = psi1(j,kp,n)

      zt           = zt   + ( f0in(kp) * x0 + f0ot(kp)   * ( one - x0 ) ) * w2dw
      asp          = asp  + ( f1in(kp)  - f1ot(kp)  )    * x1             * w2dw
      c1           = c1   +   f0in(kp)                   * x0             * w2dw
      gp           = gp   +   f1in(kp)                   * x1             * w2dw

      ztd          = ztd  + ( f0ind(kp) * x0 + f0otd(kp) * ( one - x0 ) ) * w2dw
      aspd         = aspd + ( f1ind(kp) - f1otd(kp) )    * x1             * w2dw
      c1d          = c1d  +   f0ind(kp)                  * x0             * w2dw
      gpd          = gpd  +   f1ind(kp)                  * x1             * w2dw

      ztt          = ztt  + ( f0int(kp) * x0 + f0ott(kp) * ( one - x0 ) ) * w2dw
      aspt         = aspt + ( f1int(kp) - f1ott(kp) )    * x1             * w2dw
      c1t          = c1t  +   f0int(kp)                  * x0             * w2dw
      gpt          = gpt  +   f1int(kp)                  * x1             * w2dw

      zty          = zty  + ( f0iny(kp) * x0 + f0oty(kp) * ( one - x0 ) ) * w2dw
      aspy         = aspy + ( f1iny(kp) - f1oty(kp) )    * x1             * w2dw
      c1y          = c1y  +   f0iny(kp)                  * x0             * w2dw
      gpy          = gpy  +   f1iny(kp)                  * x1             * w2dw

      ztp0(kp)     = ( f0in(kp) - f0ot(kp) )                              * w2dw
      c1p0(kp)     =   f0in(kp)                                           * w2dw
      aspp1(kp)    = ( f1in(kp) - f1ot(kp) )                              * w2dw
      gpp1(kp)     =   f1in(kp)                                           * w2dw

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO ! kp

!-----------------------------------------------------------------------
!  Compute a0w, b0w, c0w, a1w, b1w, c1w, and derivatives
!
!   scef(1) = a0w    scefi(1) = da0w/di    scefpi(1,k) = da0w/dpsii(k)
!   scef(2) = b0w    scefi(2) = db0w/di    scefpi(2,k) = db0w/dpsii(k)
!   scef(3) = c0w    scefi(3) = dc0w/di    scefpi(3,k) = dc0w/dpsii(k)
!   scef(4) = a1w    scefi(4) = da1w/di    scefpi(4,k) = da1w/dpsii(k)
!   scef(5) = b1w    scefi(5) = db1w/di    scefpi(5,k) = db1w/dpsii(k)
!   scef(6) = c1w    scefi(6) = dc1w/di    scefpi(6,k) = dc1w/dpsii(k)
!-----------------------------------------------------------------------

    scef(1,j,k,n)  = -zt
    scef(2,j,k,n)  = -asp/3.d0
    scef(3,j,k,n)  = c1
    scef(4,j,k,n)  = -asp
    scef(5,j,k,n)  = -zt
    scef(6,j,k,n)  = gp

    scefd(1,j,k,n) = -ztd
    scefd(2,j,k,n) = -aspd/3.d0
    scefd(3,j,k,n) = c1d
    scefd(4,j,k,n) = -aspd
    scefd(5,j,k,n) = -ztd
    scefd(6,j,k,n) = gpd

    sceft(1,j,k,n) = -ztt
    sceft(2,j,k,n) = -aspt/3.d0
    sceft(3,j,k,n) = c1t
    sceft(4,j,k,n) = -aspt
    sceft(5,j,k,n) = -ztt
    sceft(6,j,k,n) = gpt

    scefy(1,j,k,n) = -zty
    scefy(2,j,k,n) = -aspy/3.d0
    scefy(3,j,k,n) = c1y
    scefy(4,j,k,n) = -aspy
    scefy(5,j,k,n) = -zty
    scefy(6,j,k,n) = gpy

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      scefp0(1,kp,j,k,n) = -ztp0(kp)
      scefp0(2,kp,j,k,n) = zero
      scefp0(3,kp,j,k,n) = c1p0(kp)
      scefp0(4,kp,j,k,n) = zero
      scefp0(5,kp,j,k,n) = -ztp0(kp)
      scefp0(6,kp,j,k,n) = zero

      scefp1(1,kp,j,k,n) = zero
      scefp1(2,kp,j,k,n) = -aspp1(kp)/3.d0
      scefp1(3,kp,j,k,n) = zero
      scefp1(4,kp,j,k,n) = -aspp1(kp)
      scefp1(5,kp,j,k,n) = zero
      scefp1(6,kp,j,k,n) = gpp1(kp)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO ! kp

!-----------------------------------------------------------------------
!                 ||||| End loop over j and k |||||
!-----------------------------------------------------------------------

  END DO ! k
END DO outer

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

DEALLOCATE (scate0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scate1, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scate0d, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0d   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scate1d, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1d   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scate0t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0t   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scate1t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1t   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scate0y, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate0y   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (scate1y, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scate1y   '; WRITE (nlog,2001) var_name; END IF

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
END SUBROUTINE scterate
