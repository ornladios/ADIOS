SUBROUTINE sweepbc( nleft, nright, nmin, nmax, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         sweepbc
!    Module:       sweepbc
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/02/04
!
!    Purpose:
!      To impose boundary conditions at the ends of the 1-D swweep arrays.
!
!    Subprograms called:
!  coord_bc         : computes the ghost coordinates
!  tgvndsye_x       : determines the temperature given rho, s, and ye
!  eqstt_x          : interpolates quantities in the local EOS table
!  esrgnz_x         : reganerates local EOS table if necessary
!
!    Input arguments:
!  nleft            : left-hand boundary flags
!  nright           : right-hand boundary flags
!  nmin             : minimum paddded array index
!  nmin             : maximum paddded array index
!  ij_ray           : j-index of a radial ray
!  ik_ray           : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  eos_snc_x_module, evh1_global, evh1_sweep, evh1_bound, evh1_zone,
!  mdl_cnfg_module, mgfld_remap_module, nucbrn_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY: nnc
USE numerical_module, ONLY : zero, half, third
USE physcnst_module, ONLY: pi, G => g

USE eos_snc_x_module, ONLY: nse, aesv
USE evh1_global, ONLY : smallp, ngeomx
USE evh1_sweep, ONLY : r, u, v, w, p, ei, e, egrav, temp, ye, ge, gc, dvol, &
& xa, dx, xa0, dx0, s=>entrop, lapse_c, lapse_e
USE evh1_bound, ONLY : r_bcl, u_bcl, v_bcl, w_bcl, p_bcl, ei_bcl, temp_bcl, &
& ye_bcl, ge_bcl, gc_bcl, r_bcr, u_bcr, v_bcr, w_bcr, ei_bcr, temp_bcr, ye_bcr
USE evh1_zone, ONLY: imax
USE mdl_cnfg_module, ONLY: rho_m=>rho,t_m=>t,ye_m=>ye
USE mgfld_remap_module, ONLY : eb
USE nucbrn_module, ONLY: xn, a_nuc_rep, z_nuc_rep, be_nuc_rep

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: nleft    ! left boundary condition key
INTEGER, INTENT(in)                    :: nright   ! right boundary condition key
INTEGER, INTENT(in)                    :: nmin     ! minimum paddded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum paddded array index
INTEGER, INTENT(in)                    :: ij_ray   ! j-index of a radial ray
INTEGER, INTENT(in)                    :: ik_ray   ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                                :: first = .true.

INTEGER                                :: n        ! zone index
INTEGER                                :: i        ! composition index
INTEGER                                :: j        ! mgfld radial zone
INTEGER                                :: nmax1n   ! nmax+1-n
INTEGER                                :: nminn1   ! nmin+n-1

REAL(KIND=double)                      :: massi    ! enclosed mass
REAL(KIND=double)                      :: mass_co  ! central mass
REAL(KIND=double)                      :: gcr      ! 1/gc(nmax
REAL(KIND=double)                      :: dpdr     ! pressure gradient
REAL(KIND=double)                      :: dp       ! change in pressure
REAL(KIND=double)                      :: prat     ! pressure ratio
REAL(KIND=double)                      :: rrat     ! density ratio
REAL(KIND=double)                      :: xscale   ! velocity scaling in outer zones to enforce homology
REAL(KIND=double)                      :: t_s      ! temperature obtained assuming s_nmax(nmax+n) = s_nmax(nmax)
REAL(KIND=double)                      :: s_nmax   ! entropy of outer zone
REAL(KIND=double)                      :: ei_j     ! internal energy
REAL(KIND=double)                      :: dedd     ! derivative of the energy wrt density
REAL(KIND=double)                      :: dedt     ! derivative of the temperature wrt density
REAL(KIND=double)                      :: dedy     ! derivative of the electron fraction wrt density

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : fixed (eg, u_bcl,p_bcr,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find ghost coordinates
!-----------------------------------------------------------------------

! print *,' Calling coord_bc from sweepbc'
CALL coord_bc( nleft, nright, nmin, nmax, xa, dx, xa0, dx0, nmax+6 )

!-----------------------------------------------------------------------
!  Load left (inner) ghosts
!-----------------------------------------------------------------------

IF( nleft == 0 ) THEN          ! symmetric accross left (inner) edge

  DO n = 1, 6
    nminn1               = MIN( nmin + n - 1, nmax )
    r      (nmin-n)      = r      (nminn1)
    u      (nmin-n)      = -u     (nminn1)
    v      (nmin-n)      = v      (nminn1)
    w      (nmin-n)      = w      (nminn1)
    p      (nmin-n)      = p      (nminn1)
    s      (nmin-n)      = s      (nminn1)
    ei     (nmin-n)      = ei     (nminn1)
    eb     (nmin-n)      = eb     (nminn1)
    e      (nmin-n)      = e      (nminn1)
    egrav  (nmin-n)      = egrav  (nminn1)
    temp   (nmin-n)      = temp   (nminn1)
    ye     (nmin-n)      = ye     (nminn1)
    ge     (nmin-n)      = ge     (nminn1)
    gc     (nmin-n)      = gc     (nminn1)
    lapse_c(nmin-n)      = lapse_c(nminn1)
    lapse_e(nmin-n)      = lapse_e(nminn1)
  END DO

ELSE IF ( nleft == 1 ) THEN    ! Zero Gradient 

  DO n = 1, 6
    r      (nmin-n)      = r      (nmin)
    u      (nmin-n)      = u      (nmin)
    v      (nmin-n)      = v      (nmin)
    w      (nmin-n)      = w      (nmin)
    p      (nmin-n)      = p      (nmin)
    s      (nmin-n)      = s      (nmin)
    ei     (nmin-n)      = ei     (nmin)
    eb     (nmin-n)      = eb     (nmin)
    e      (nmin-n)      = e      (nmin)
!         egrav (nmin-n) = egrav (nmin)
    ye     (nmin-n)      = ye     (nmin)
    temp   (nmin-n)      = temp   (nmin)
    ge     (nmin-n)      = ge     (nmin)
    gc     (nmin-n)      = gc     (nmin)
    lapse_c(nmin-n)      = lapse_c(nmin)
    lapse_e(nmin-n)      = lapse_e(nmin)
  END DO

ELSE IF ( nleft == 2 ) THEN   ! Externally Fixed

  DO n = 1, 6
    r      (nmin-n)      = r_bcl
    u      (nmin-n)      = u_bcl
    v      (nmin-n)      = v_bcl
    w      (nmin-n)      = w_bcl
    p      (nmin-n)      = p_bcl
    ei     (nmin-n)      = ei_bcl
    eb     (nmin-n)      = eb(nmin)
    e      (nmin-n)      = ei_bcl + 0.5d0 * ( u_bcl**2 + v_bcl**2 + w_bcl**2 )
    temp   (nmin-n)      = temp_bcl
    ye     (nmin-n)      = ye_bcl
    ge     (nmin-n)      = ge_bcl
    gc     (nmin-n)      = gc_bcl
    lapse_c(nmin-n)      = lapse_c(nmin)
    lapse_e(nmin-n)      = lapse_e(nmin)
  END DO

ELSE IF ( nleft == 3 ) THEN   ! Periodic

  DO n = 1, 6
    r      (nmin-n)      = r      (nmax+1-n)
    u      (nmin-n)      = u      (nmax+1-n)
    v      (nmin-n)      = v      (nmax+1-n)
    w      (nmin-n)      = w      (nmax+1-n)
    p      (nmin-n)      = p      (nmax+1-n)
    ei     (nmin-n)      = ei     (nmax+1-n)
    eb     (nmin-n)      = eb     (nmax+1-n)
    e      (nmin-n)      = e      (nmax+1-n)
    temp   (nmin-n)      = temp   (nmax+1-n)
    ye     (nmin-n)      = ye     (nmax+1-n)
!         mf(:,nmin-n) = mf(:,nmax+1-n)
    ge     (nmin-n)      = ge     (nmax+1-n)
    gc     (nmin-n)      = gc     (nmax+1-n)
    lapse_c(nmin-n)      = lapse_c(nmax+1-n)
    lapse_e(nmin-n)      = lapse_e(nmax+1-n)
  END DO

END IF

!-----------------------------------------------------------------------
!  Load right (outer) ghosts
!-----------------------------------------------------------------------

IF ( nright == 0 ) THEN          ! symmetric accross right (outer) edge

  DO n = 1, 6
    nmax1n               = MAX( nmax + 1 - n, nmin )
    r      (nmax+n)      = r      (nmax1n)
    u      (nmax+n)      = -u     (nmax1n)
    v      (nmax+n)      = v      (nmax1n)
    w      (nmax+n)      = w      (nmax1n)
    p      (nmax+n)      = p      (nmax1n)
    s      (nmax+n)      = s      (nmax1n)
    ei     (nmax+n)      = ei     (nmax1n)
    eb     (nmax+n)      = eb     (nmax1n)
    e      (nmax+n)      = e      (nmax1n)
    temp   (nmax+n)      = temp   (nmax1n)
    ye     (nmax+n)      = ye     (nmax1n)
    ge     (nmax+n)      = ge     (nmax1n)
    gc     (nmax+n)      = gc     (nmax1n)
    lapse_c(nmax+n)      = lapse_c(nmax1n)
    lapse_e(nmax+n)      = lapse_e(nmax1n)
  END DO

ELSE IF ( nright == 1 ) THEN    ! Zero Gradient 

  DO n = 1, 6
    r      (nmax+n)      = r (nmax)
    u      (nmax+n)      = u (nmax)
    v      (nmax+n)      = v (nmax)
    w      (nmax+n)      = w (nmax)
    p      (nmax+n)      = p (nmax)
    s      (nmax+n)      = s (nmax)
    ei     (nmax+n)      = ei(nmax)
    eb     (nmax+n)      = eb(nmax)
    egrav  (nmax+n)      = egrav(nmax)
    e      (nmax+n)      = e (nmax)
    temp   (nmax+n)      = temp(nmax)
    ye     (nmax+n)      = ye(nmax)
!          mf(:,nmax+n) = mf(:,nmax)
    ge     (nmax+n)      = ge(nmax)
    gc     (nmax+n)      = gc(nmax)
    lapse_c(nmax+n)      = lapse_c(nmax)
    lapse_e(nmax+n)      = lapse_e(nmax)
  END DO

ELSE IF ( nright == 2 ) THEN   ! Externally Fixed

  DO n = 1, 6
    r      (nmax+n)      = r_bcr
    rho_m  (nmax+n-5)    = r_bcr
    u      (nmax+n)      = u_bcr
    v      (nmax+n)      = v_bcr
    w      (nmax+n)      = w_bcr
    egrav  (nmax+n)      = egrav(nmax)
    ei     (nmax+n)      = ei_bcr
    eb     (nmax+n)      = eb(nmax)
    e      (nmax+n)      = ei_bcr + 0.5d0 * ( u_bcr**2 + v_bcr**2 + w_bcr**2 )
    temp   (nmax+n)      = temp_bcr
    lapse_c(nmax+n)      = lapse_c(nmax)
    lapse_e(nmax+n)      = lapse_e(nmax)
    t_m    (nmax+n-5)    = temp_bcr
    ye     (nmax+n)      = ye_bcr
    ye_m   (nmax+n-5)    = ye_bcr
    nse    (nmax+n-5,ij_ray,ik_ray) = nse(nmax-5,ij_ray,ik_ray)
    DO i = 1,nnc
      xn   (nmax+n-5,i)  = xn(nmax-5,i)
    END DO
    a_nuc_rep(nmax+n-5)  = a_nuc_rep(nmax-5)
    z_nuc_rep(nmax+n-5)  = z_nuc_rep(nmax-5)
    be_nuc_rep(nmax+n-5) = be_nuc_rep(nmax-5)
!         mf(:,nmax+n) = mf_bcr
!         p (nmax+n) = p_bcr
!         ge(nmax+n) = ge_bcr
!         gc(nmax+n) = gc_bcr
  END DO

!  print *,' Calling esrgnz_x'
  CALL esrgnz_x( nmax+1-5, nmax+6-5, rho_m, t_m, ye_m, ij_ray, ik_ray )

!  print *,' Calling tgvndeye_sweep_x'
  CALL tgvndeye_sweep_x( nmax+1, nmax+6, ij_ray, ik_ray, r, r )
!  print *,' Returning from tgvndeye_sweep_x'

ELSE IF ( nright == 3 ) THEN   ! Periodic

  DO n = 1, 6
    r      (nmax+n)      = r      (nmin+n-1)
    u      (nmax+n)      = u      (nmin+n-1)
    v      (nmax+n)      = v      (nmin+n-1)
    w      (nmax+n)      = w      (nmin+n-1)
    ei     (nmax+n)      = ei     (nmin+n-1)
    eb     (nmax+n)      = eb     (nmin+n-1)
    e      (nmax+n)      = e      (nmin+n-1)
    temp   (nmax+n)      = temp   (nmin+n-1)
    ye     (nmax+n)      = ye     (nmin+n-1)
!          mf(:,nmax+n) = mf(:,nmin+n-1)
    p      (nmax+n)      = p      (nmin+n-1)
    ge     (nmax+n)      = ge     (nmin+n-1)
    gc     (nmax+n)      = gc     (nmin+n-1)
    lapse_c(nmax+n)      = lapse_c(nmin+n-1)
    lapse_e(nmax+n)      = lapse_e(nmin+n-1)
  END DO

ELSE IF ( nright == 4 ) THEN   ! Balance pressure with gravity  

!........Compute volume elements

  CALL volume ( ngeomx ) 

!........Calculate mass interior

  massi                  = mass_co
  DO n = nmin,nmax
    massi                = massi + 4.0d0 * pi * r(n) * dvol(n)
  END DO
  gcr                    = 1.0d0/gc(nmax)
  dpdr                   = ( p(nmax) - p(nmax-1) )/( xa(nmax) - xa(nmax-1 ) &
&                        + 0.5d0 * ( dx(nmax) - dx(nmax-1) ) )

  DO n = 1, 6

!........Gamma constant

    ge(nmax+n)           = ge(nmax)
    gc(nmax+n)           = gc(nmax)

!........Approximate the hydrostatic dp, using the previous density

    dp                   = -G * massi * dx(nmax+n) * r(nmax+n-1)/xa(nmax+n)**2
    p(nmax+n)            = p(nmax+n-1) + dp
    p(nmax+n)            = DMAX1( smallp, p(nmax+n) )

!  Calculate P assuming constant dP/dr
!         p (nmax+n) = p (nmax)+ dpdr*
!    &      (xa(nmax+n)-xa(nmax)+.5*(dx(nmax+n)-dx(nmax)))
!	  p (nmax+n) = max(smallp,p(nmax+n))

!........Set rho to keep P/rho**gamma constant

    r(nmax+n)            = r(nmax+n-1) * (p(nmax+n)/p(nmax+n-1))**gcr
!    ei(nmax+n)        = p(nmax+n)/(r(nmax+n) * ( ge(nmax+n)-1.0d0) )

!-----------------------------------------------------------------------
!  Enforce homology in the velocity
!         xscale=xa(nmax+n)/xa(nmax) 
!         u (nmax+n) = u (nmax)*xscale
!  Constant velocity
!-----------------------------------------------------------------------

    xscale               = xa(nmax+n)/xa(nmax)
    u      (nmax+n)      = u      (nmax) * xscale
    v      (nmax+n)      = v      (nmax)
    w      (nmax+n)      = w      (nmax)
    egrav  (nmax+n)      = egrav  (nmax)
    ei     (nmax+n)      = ei     (nmax)
    eb     (nmax+n)      = eb     (nmax)
    e      (nmax+n)      = ei(nmax+n) + 0.5d0 * ( u(nmax+n)**2 + v(nmax+n)**2 + w(nmax+n)**2 )
    temp   (nmax+n)      = temp   (nmax)
    ye     (nmax+n)      = ye     (nmax)
    lapse_c(nmax+n)      = lapse_c(nmax)
    lapse_e(nmax+n)      = lapse_e(nmax)
!          mf(:,nmax+n) = mf(:,nmax)
  END DO

ELSE IF ( nright == 5 ) THEN  ! Scaled Down P and rho

  prat                    = p(nmax)/p(nmax-1)
  rrat                    = r(nmax)/r(nmax-1)
  gcr                     = 1.0d0/gc(nmax)                  

  DO n=1,6

!........Constant gamma

    ge(nmax+n)            = ge(nmax)
    gc(nmax+n)            = gc(nmax)

!........Scale down the pressure by the ratio of the last real zones

    p(nmax+n)             = p(nmax+n-1) * prat
    p(nmax+n)             = DMAX1( smallp, p(nmax+n) )

!........Set rho to keep P/rho**gamma constant

    r(nmax+n)             = r(nmax+n-1) * (p(nmax+n)/p(nmax+n-1))**gcr
    ei(nmax+n)            = p(nmax+n)/( r(nmax+n) * ( ge(nmax+n) - 1.0d0 ) )

!-----------------------------------------------------------------------
!  Enforce homology in the velocity
!         xscale=xa(nmax+n)/xa(nmax) 
!         u (nmax+n) = u (nmax)*xscale
!  Constant velocity
!-----------------------------------------------------------------------

    u      (nmax+n)      = u      (nmax)
    v      (nmax+n)      = v      (nmax)
    w      (nmax+n)      = w      (nmax)
    eb     (nmax+n)      = eb     (nmax)
    e      (nmax+n)      = ei     (nmax+n) + 0.5d0 * ( u(nmax+n)**2 + v(nmax+n)**2 + w(nmax+n)**2 )
    temp   (nmax+n)      = temp   (nmax)
    ye     (nmax+n)      = ye     (nmax)
    lapse_c(nmax+n)      = lapse_c(nmax)
    lapse_e(nmax+n)      = lapse_e(nmax)

  END DO

ELSE IF ( nright == 6 ) THEN    ! Zero Gradient (in s rather than t)

!........Compute volume elements

! print *,' Calling volume from sweepbc'
  CALL volume ( ngeomx )

!........Calculate mass interior

  massi                   = mass_co
  DO n = nmin,nmax
    massi                 = massi + 4.0d0 * pi * r(n) * dvol(n)
  END DO

!........Calculate P assuming constant dP/dr

    dpdr                  = ( p(nmax) - p(nmax-1) )/( half * ( xa(nmax+1) - xa(nmax-1 ) ) )
    DO n = 1, 5
      p(nmax+n)           = p(nmax) + dpdr * ( half * ( xa(nmax+n+1) + xa(nmax+n) ) &
&                         - half * ( xa(nmax+1) + xa(nmax) ) )
      p(nmax+n)           = DMAX1( 1.d-1 * p(nmax), p(nmax+n) )
    END DO

!........Calculate egrav

    DO n = 1, 5
      egrav(nmax+n)       = - G * massi/( half * ( xa(nmax+n)**3 + xa(nmax+n+1)**3 ) )**third
    END DO

!-----------------------------------------------------------------------
!        Calculate u
!
!  Approximate ghost u as following free-fall profile if u(nmax) < 0
!  Set u(nmax+n) = u(nmax) if u(nmax) > 0
!-----------------------------------------------------------------------

    DO n = 1, 5
      IF ( u(nmax) < zero ) THEN
        u(nmax+n)        = u(nmax-10) * ( half * ( xa(nmax-10+1) + xa(nmax-10) ) &
&                        /( half * ( xa(nmax+n+1) + xa(nmax+n) ) ) )**half
      ELSE
        u(nmax+n)        = u(nmax)
      END IF
    END DO

!-----------------------------------------------------------------------
!        Entropy boundary condition
!
!  s_nmax is initially set to the initial entropy of the outer zone
!  If u(nmax) < 0, set s(nmax) to s_nmax
!  If u(nmax) > 0, reset s_nmax to s(nmax)
!
!  If u(nmax) < zero and r(nmax) < r_bcr, or if u(nmax) >= zero,
!   set r_bcr to r(nmax)
!-----------------------------------------------------------------------

  IF ( first ) THEN
    first                = .false.
    s_nmax               = aesv(nmax-5,3,ij_ray,ik_ray)
  END IF


  IF ( u(nmax) < zero  .and.  r(nmax) < r_bcr ) THEN
    r_bcr                = r(nmax)
  ELSE IF ( u(nmax) >= zero ) THEN
    r_bcr                = r(nmax)
  END IF

  IF ( u(nmax) < zero ) THEN
    CALL tgvndsye_x(nmax-5, ij_ray, ik_ray, r_bcr, s_nmax, ye(nmax),    &
&    temp(nmax), t_s)
    CALL eqstt_x(2, nmax-5, ij_ray, ik_ray, r_bcr, t_s, ye(nmax), ei_j, &
&    dedd, dedt, dedy)
  ELSE
    s_nmax               = aesv(nmax-5,3,ij_ray,ik_ray)
    ei_j                 = ei(nmax)
    t_s                  = temp(nmax)
  END IF
    
  DO n = 1, 6
    DO i = 1,nnc
      xn(nmax+n-5,i)     = xn(nmax-5,i)
    END DO
    a_nuc_rep(nmax+n-5)  = a_nuc_rep(nmax-5)
    z_nuc_rep(nmax+n-5)  = z_nuc_rep(nmax-5)
    be_nuc_rep(nmax+n-5) = be_nuc_rep(nmax-5)
    j                    = nmax + n - 5

    r      (nmax+n)      = r_bcr
    v      (nmax+n)      = v_bcr
    w      (nmax+n)      = w_bcr
    ei     (nmax+n)      = ei_j
    eb     (nmax+n)      = eb(nmax)
    e      (nmax+n)      = ei(nmax+n) + 0.5d0 * ( u(nmax+n)**2 + v(nmax+n)**2 + w(nmax+n)**2 )
    temp   (nmax+n)      = t_s
    s      (nmax+n)      = s_nmax
    ye     (nmax+n)      = ye     (nmax)
    ge     (nmax+n)      = ge     (nmax)
    gc     (nmax+n)      = gc     (nmax)
    lapse_c(nmax+n)      = lapse_c(nmax)
    lapse_e(nmax+n)      = lapse_e(nmax)

  END DO !  n = 1, 6

END IF

RETURN
END SUBROUTINE sweepbc

