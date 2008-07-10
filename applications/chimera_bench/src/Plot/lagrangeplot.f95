SUBROUTINE lagrangeplot( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         lagrangeplot
!    Module:       lagrangeplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To creates a file for specified variables as a function
!       of time for a specified Lagrange point.
!
!    Subprograms called:
!  date_and_time_print : prints the cycle number, date, and time
!  eqstta_x            : computes EOS quantities directly
!  abem_cal            : computes absorption and emission inverse mean free paths
!
!    Input arguments:
!  ij_ray              : index denoting the j-index of a specific radial ray
!  ik_ray              : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module mdl_cnfg_module,  nucbrn_module,
!  nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nx
USE numerical_module, ONLY : zero, half, one, epsilon, frpi, ncoef
USE physcnst_module, ONLY : kmev, msolar, cvel, rmu

USE edit_module, ONLY : ilagplt, dtimeplot, nlagplt, nlagdump, nprint,  &
& msslag, head, data_path
USE eos_snc_x_module, ONLY : nse, aesv, nuc_number
USE mdl_cnfg_module, ONLY : jr_max, u, r, rho, t, ye, rstmss
USE nucbrn_module, ONLY : a_name, xn
USE nu_dist_module, ONLY : unu, dunu, stwt, psi0
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY : time, dtnph

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)               :: ik_ray        ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)               :: lagfile

LOGICAL                           :: lprint
LOGICAL                           :: first_c = .true.

INTEGER                           :: i             ! composition index
INTEGER                           :: j             ! radial zone index
INTEGER                           :: jd            ! particular radial zone index
INTEGER                           :: jt            ! radial index for a particular temperature
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index
INTEGER                           :: istat         ! open and CLOSE file flag
INTEGER                           :: n_time        ! unit number to print the cycle number, date, and time

INTEGER                           :: itime         ! parameter for criterion for writing to nuradplot files
INTEGER                           :: itimeprev     ! parameter for criterion for writing to nuradplot files

INTEGER                           :: it            ! iteration index
INTEGER, PARAMETER                :: itmax = 30    ! maximum number of iterations

INTEGER                           :: i_n           ! neutron abundance index
INTEGER                           :: i_p           ! proton abundance index
INTEGER                           :: i_4He         ! 4He abundance index

REAL(KIND=double)                 :: T_up = 8.d-1
REAL(KIND=double)                 :: T_down = 2.d-1
REAL(KIND=double)                 :: tol = 1.d-3

REAL(KIND=double), DIMENSION(nx)  :: x_n           ! neutron mass fraction
REAL(KIND=double), DIMENSION(nx)  :: x_p           ! proton mass fraction
REAL(KIND=double), DIMENSION(nx)  :: x_He          ! helium mass fraction
REAL(KIND=double), DIMENSION(nx)  :: x_heavy       ! heavy nuclei mass fraction

REAL(KIND=double)                 :: tmult         ! multiplier for criterion for writing to nuradplot files

REAL(KIND=double)                 :: t_msslag      ! temperature of the Lagrangian point
REAL(KIND=double)                 :: rstmssjd      ! rest mass enclosed by radial index jd
REAL(KIND=double)                 :: rstmssjdm1    ! rest mass enclosed by radial index jd-1
REAL(KIND=double)                 :: u_msslag      ! radial velocity of the Lagrangian point
REAL(KIND=double)                 :: r_msslag      ! radius of the Lagrangian point
REAL(KIND=double)                 :: rho_msslag    ! density of the Lagrangian point
REAL(KIND=double)                 :: s_msslag      ! entropy of the Lagrangian point
REAL(KIND=double)                 :: ye_msslag     ! electron fraction of the Lagrangian point
REAL(KIND=double)                 :: xneut_msslag  ! neutron mass fraction of the Lagrangian point
REAL(KIND=double)                 :: xprot_msslag  ! proton mass fraction of the Lagrangian point
REAL(KIND=double)                 :: xhe_msslag    ! helium mass fraction of the Lagrangian point
REAL(KIND=double)                 :: xhv_msslag    ! heacy nuclei mass fraction of the Lagrangian point
REAL(KIND=double)                 :: tau_jd        ! expansion time scale of material at radial index jd
REAL(KIND=double)                 :: tau_jdm1      ! expansion time scale of material at radial index jd-1
REAL(KIND=double)                 :: tau_msslag    ! expansion time scale of material at Lagrangian point
REAL(KIND=double)                 :: mdot_jd       ! mass luminosity of material at radial index jd
REAL(KIND=double)                 :: mdot_jdm1     ! mass luminosity of material at radial index jd-1
REAL(KIND=double)                 :: mdot_msslag   ! mass luminosity of Lagrangian point
REAL(KIND=double)                 :: ye_max        ! upper value of ye for bisection iteration
REAL(KIND=double)                 :: ye_min        ! lower value of ye for bisection iteration
REAL(KIND=double)                 :: ye_test       ! test value of ye
REAL(KIND=double)                 :: dye_dt        ! d(ye)/d(time) for the lagrangian mass point
REAL(KIND=double)                 :: dye_dtmax     ! maximum d(ye)/d(time) for convergence tenting
REAL(KIND=double)                 :: xneut         ! neutron mass fraction
REAL(KIND=double)                 :: xprot         ! proton mass fraction
REAL(KIND=double)                 :: cmpn          ! neutron chemical potential
REAL(KIND=double)                 :: cmpp          ! proton chemical potential
REAL(KIND=double)                 :: cmpe          ! electron chemical potential
REAL(KIND=double)                 :: xh            ! heavy nucleus mass fraction
REAL(KIND=double)                 :: ah            ! heavy nucleus mass number
REAL(KIND=double)                 :: zh            ! heavy nucleus charge number
REAL(KIND=double)                 :: dum1          ! dummy index
REAL(KIND=double)                 :: dum2          ! dummy index
REAL(KIND=double)                 :: dum3          ! dummy index
REAL(KIND=double)                 :: absornp       ! absorption inverse mean free path
REAL(KIND=double)                 :: emitnp        ! emission inverse mean free path
REAL(KIND=double)                 :: ncoefap       ! coefficient for computing neutrino number
REAL(KIND=double)                 :: dpsi_dt       ! d(psi0)/d(time)
REAL(KIND=double)                 :: yeeq_jd       ! equilibrium value of ye at radial zone jd
REAL(KIND=double)                 :: yeeq_jdm1     ! equilibrium value of ye at radial zone jd-1
REAL(KIND=double)                 :: yeeq_msslag   ! equilibrium value of ye at Lagrangian point
REAL(KIND=double)                 :: yeeqab_jd     ! kinetic equilibrium value of ye at radial zone jd
REAL(KIND=double)                 :: yeeqab_jdm1   ! kinetic equilibrium value of ye at radial zone jd-1
REAL(KIND=double)                 :: yeeqab_msslag ! kinetic equilibrium value of ye at Lagrangian point

!      double precision ah,absornp,cmpe,cmpn,cmpp,
!     *                 dpsi_dt,dum1,dum2,dum3,dye_dt,dye_dtmax,
!     *                 emitnp,mdot_jd,mdot_jdm1,mdot_msslag,ncoefap,
!     *                 rho_msslag,r_msslag,rstmssjd,rstmssjdm1,
!     *                 s_msslag,tau_jd,tau_jdm1,tau_msslag,T_down,tol,
!     *                 t_msslag,tmult,T_up,u_msslag,xh,
!     *                 x_He,xhe_msslag,x_heavy,xhv_msslag,
!     *                 x_n,xneut,xneut_msslag,xprot,xprot_msslag,x_p,
!     *                 yeeq_jd,yeeq_jdm1, yeeqab_jd,yeeqab_jdm1,
!     *                 yeeq_msslag, yeeqab_msslag,
!     *                 ye_max,ye_min,ye_msslag,ye_test,zh
!                                                                      c
!                                                                      c

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    3 format (1x,a128)
    5 format (1x,128('*'))
    7 format (1x,'time =',1pe15.8,' mass_enclosed =',1pe15.8)
  101 format ('     time         u          r         rho         T     &
&     s         ye         x_n        x_p       x_He      x_heavy')
  103 format ('     sec         cm/s       cm        g/cm3       MeV    /baryon'/)
  105 format (1x,11(1pe11.3))
  201 format ('     time        tau       m_dot       ye        ye_eq   &
&    ye_eq')
  203 format ('     sec         sec     sol mass/s             em & ab  &
&   ab only'/)
  205 format (1x,6(1pe11.3))
 1001 format (' jt cannot be found in subroutine lagrangeplot, T_up=',  &
& 1pe10.3,' t(2)=',1pe10.3,' t(jr_max)=',1pe10.3)
 1201 format (' jd cannot be found in subroutine lagrangeplot, msslag=',&
& 1pe10.3,' rstmss(jr_max)=',1pe10.3)
 8001 format (' File lagplot1.d cannot be opened in subroutime lagrangeplot')
 8501 format (' File lagplot2.d cannot be opened in subroutime lagrangeplot')
 9001 format (' File lagplot1.d cannot be closed in subroutime lagrangeplot')
 9501 format (' File lagplot2.d cannot be closed in subroutime lagrangeplot')

!-----------------------------------------------------------------------
!  Return IF ivarplt = 0
!-----------------------------------------------------------------------

IF ( ilagplt == 0 ) RETURN

!-----------------------------------------------------------------------
!  Time criterion
!-----------------------------------------------------------------------

lprint             = .false.
tmult              = 1.d+3/dtimeplot
itime              = INT( time * tmult )
itimeprev          = INT( ( time - dtnph ) * tmult )
IF ( itime > itimeprev ) lprint = .true.

IF ( .not. lprint ) RETURN

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

n_time             = nlagplt

!-----------------------------------------------------------------------
!  If nlagdump = 0, find msslag such that t = T_max at the radial shell
!   enclosing a rest mass = msslag; set newfile = true.
!-----------------------------------------------------------------------

IF ( nlagdump == 0 ) then
  DO j = 2,jr_max
    IF ( t(j) * kmev < T_up ) THEN
      jt           = j
      GO TO 1100
    END IF ! t(j) < T_up
  END DO ! j = 2,jr_max
  WRITE (nprint,1001) T_up, t(2), t(jr_max)
  RETURN

 1100 CONTINUE
  msslag            = rinterp( rstmss(jt)/msolar, rstmss(jt-1)/msolar,  &
&                     t(jt) * kmev, T_up, t(jt-1) * kmev )
  nlagdump          = nlagdump + 1
END IF ! nlagdump = 0

!-----------------------------------------------------------------------
!  Find jd such that rstmss(jd) > msslag > rstmss(jd-1)
!-----------------------------------------------------------------------

DO j = 2,jr_max
  IF ( rstmss(j)/msolar >= msslag ) THEN
    jd              = j
    GO TO 1300
  END IF ! rstmss(j)/msolar > msslag
  WRITE (nprint,1201) msslag,rstmss(jr_max)/msolar
  RETURN
END DO ! j = 2,jr_max

 1300 CONTINUE

!-----------------------------------------------------------------------
!  If t_msslag < T_down, find msslag such that t = T_max at the radial
!   shell enclosing a rest mass = msslag; set newfile = true
!-----------------------------------------------------------------------

t_msslag            = rinterp( t(jd) * kmev, t(jd-1) * kmev,            &
&                     rstmss(jd)/msolar, msslag, rstmss(jd-1)/msolar )
IF ( t_msslag < T_down ) THEN

  DO j = 2,jr_max
    IF ( t(j) * kmev < T_up ) THEN
      jt            = j
      GO TO 1500
    END IF ! t(j) < T_up
  END DO ! j = 2,jr_max
  WRITE (nprint,1001) T_up,t(2)*kmev, t(jr_max)*kmev
  RETURN
 1500  CONTINUE
  msslag            = rinterp( rstmss(jt)/msolar, rstmss(jt-1)/msolar,   &
&                     t(jt)*kmev, T_up, t(jt-1)*kmev )
  nlagdump          = nlagdump + 1

!-----------------------------------------------------------------------
!  Find jd such that rstmss(jd) > msslag > rstmss(jd-1)
!-----------------------------------------------------------------------

  DO j = 2,jr_max
    IF ( rstmss(j)/msolar >= msslag ) THEN
      jd            = j
      GO TO 1700
    END IF ! rstmss(j)/msolar > msslag
  END DO
  WRITE (nprint,1201) msslag,rstmss(jr_max)/msolar
  RETURN
 1700  CONTINUE

END IF ! t_msslag < T_down

!-----------------------------------------------------------------------
!  Composition indexing
!-----------------------------------------------------------------------

IF ( first_c ) THEN
  first_c           = .false.
  i_n               = nuc_number + 1
  i_p               = nuc_number + 1
  i_4He             = nuc_number + 1
  DO i = 1,nuc_number
    IF (      a_name(i) == '  n  ' ) THEN
      i_n           = i
    END IF 
    IF (      a_name(i) == '  p  ' ) THEN
      i_p           = i
    END IF 
    IF (      a_name(i) == '  4He' ) THEN
      i_4He         = i
    END IF 
  END DO ! i
END IF ! first_c

!-----------------------------------------------------------------------
!  Determine abundance of free neutrons, free protons, helium nuclei
!   and heavies
!-----------------------------------------------------------------------

IF ( nse(jd-1,ij_ray,ik_ray) == 0 ) THEN
  x_n(jd-1)         = xn(jd-1,i_n)
  x_p(jd-1)         = xn(jd-1,i_p)
  x_He(jd-1)        = xn(jd-1,i_4He)
  x_heavy(jd-1)     = DMAX1( 1.d0 - x_n(jd-1) - x_p(jd-1) - x_He(jd-1), zero )
ELSE
  x_n(jd-1)         = aesv(jd-1,7,ij_ray,ik_ray)
  x_p(jd-1)         = aesv(jd-1,8,ij_ray,ik_ray)
  x_heavy(jd-1)     = aesv(jd-1,9,ij_ray,ik_ray)
  x_He(jd-1)        = DMAX1( one - x_n(jd-1) - x_p(jd-1) - x_heavy(jd-1), zero )
END IF ! nse(jd-1,ij_ray,ik_ray) = 0

IF ( nse(jd,ij_ray,ik_ray) == 0 ) THEN
  x_n(jd)           = xn(jd,i_n)
  x_p(jd)           = xn(jd,i_p)
  x_He(jd)          = xn(jd,i_4He)
  x_heavy(jd)       = DMAX1( 1.d0 - x_n(jd) - x_p(jd) - x_He(jd), zero )
ELSE
  x_n(jd)           = aesv(jd,7,ij_ray,ik_ray)
  x_p(jd)           = aesv(jd,8,ij_ray,ik_ray)
  x_heavy(jd)       = aesv(jd,9,ij_ray,ik_ray)
  x_He(jd)          = DMAX1( one - x_n(jd) - x_p(jd) - x_heavy(jd), zero )
END IF

!-----------------------------------------------------------------------
!  Interpolate variables in mass to msslag
!-----------------------------------------------------------------------

rstmssjd            = rstmss(jd  )/msolar
rstmssjdm1          = rstmss(jd-1)/msolar

u_msslag            = rinterp( u(jd), u(jd-1),                        &
&                     rstmssjd, msslag, rstmssjdm1 )
r_msslag            = rinterp( r(jd), r(jd-1),                        &
&                     rstmssjd, msslag, rstmssjdm1 )
rho_msslag          = rinterp( rho(jd), rho(jd-1),                    &
&                     rstmssjd, msslag, rstmssjdm1 )
t_msslag            = rinterp( t(jd)*kmev, t(jd-1) * kmev,            &
&                     rstmssjd, msslag, rstmssjdm1 )
s_msslag            = rinterp( aesv(jd,3,ij_ray, ik_ray),             &
&                     aesv(jd-1,3,ij_ray, ik_ray),                    &
&                     rstmssjd, msslag, rstmssjdm1 )
ye_msslag           = rinterp( ye(jd), ye(jd-1),                      &
&                     rstmssjd, msslag, rstmssjdm1 )
xneut_msslag        = rinterp( x_n(jd), x_n(jd-1),                    &
&                     rstmssjd, msslag, rstmssjdm1 )
xprot_msslag        = rinterp( x_p(jd), x_p(jd-1),                    &
&                     rstmssjd, msslag, rstmssjdm1 )
xhe_msslag          = rinterp( x_He(jd), x_He(jd-1),                  &
&                     rstmssjd, msslag, rstmssjdm1 )
xhv_msslag          = rinterp( x_heavy(jd), x_heavy(jd-1),            &
&                     rstmssjd, msslag, rstmssjdm1 )

!-----------------------------------------------------------------------
!  Write to nlagplt ("variable lagplot files")
!-----------------------------------------------------------------------

WRITE (lagfile,'(a20,i3.3,a2)') '/Plot_Files/lagplot1_',nlagdump,'.d'
lagfile             = TRIM(data_path)//TRIM(lagfile)

!-----------------------------------------------------------------------
!  Open lagplot1.d file.                                         c
!-----------------------------------------------------------------------

OPEN (UNIT=nlagplt, FILE=TRIM(lagfile), STATUS='new', IOSTAT=istat)
IF ( istat == 0 ) THEN
  WRITE (nlagplt,3) head
  WRITE (nlagplt,5)
  WRITE (nlagplt,7) time,msslag
  CALL date_and_time_print( n_time )
  WRITE (nlagplt,101)
  WRITE (nlagplt,103)
END IF ! istat == 0
      
IF ( istat /= 0 ) OPEN (UNIT=nlagplt,FILE=lagfile, STATUS='old',        &
&  POSITION='append', IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,8001)
  STOP
END IF ! istat /= 0

WRITE (nlagplt,105) time, u_msslag, r_msslag, rho_msslag, t_msslag,     &
& s_msslag, ye_msslag, xneut_msslag, xprot_msslag, xhe_msslag, xhv_msslag

!-----------------------------------------------------------------------
!  Close lagplot1.d file
!-----------------------------------------------------------------------

CLOSE (UNIT=nlagplt, STATUS='keep', ERR=9000)

!-----------------------------------------------------------------------
!  Compute the expansion timescale
!-----------------------------------------------------------------------

tau_jd              = r(jd  )/( 2.d0 * u(jd  ) + epsilon )
tau_jdm1            = r(jd-1)/( 2.d0 * u(jd-1) + epsilon )

!-----------------------------------------------------------------------
!  Compute the mass loss rate
!-----------------------------------------------------------------------

mdot_jd             = frpi * r(jd  )**2 * u(jd  ) * rho(jd  )/msolar
mdot_jdm1           = frpi * r(jd-1)**2 * u(jd-1) * rho(jd-1)/msolar

!-----------------------------------------------------------------------
!  Compute the kinetic equilibrium value of ye at radial zone jd with
!   both emission and absorption.
!-----------------------------------------------------------------------

ye_max              = ye(jd) + 0.02d0
ye_min              = ye(jd) - 0.02d0
DO it = 1,itmax

  ye_test           = half * ( ye_min + ye_max )
  dye_dt            = zero
  dye_dtmax         = zero
  CALL eqstta_x( 7 , jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, xneut, dum1, dum2, dum3 )
  CALL eqstta_x( 8 , jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, xprot, dum1, dum2, dum3 )
  CALL eqstta_x( 4 , jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, cmpn , dum1, dum2, dum3 )
  CALL eqstta_x( 5 , jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, cmpp , dum1, dum2, dum3 )
  CALL eqstta_x( 6 , jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, cmpe , dum1, dum2, dum3 )
  CALL eqstta_x( 9 , jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, xh   , dum1, dum2, dum3 )
  CALL eqstta_x( 10, jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, ah   , dum1, dum2, dum3 )
  CALL eqstta_x( 11, jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, zh   , dum1, dum2, dum3 )

  n                 = 1 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd,k), rho(jd), t(jd), xneut, xprot, xh, ah,    &
&    zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd,k)**2 * dunu(jd,k)
    dpsi_dt         = cvel * ( emitnp * ( one - psi0(jd,k,n) )            &
&                   - absornp * psi0(jd,k,n) )
    dye_dt          = dye_dt - ( ncoefap * rmu/rho(jd) ) * dpsi_dt
    dye_dtmax       = DMAX1( DABS( ( ncoefap * rmu/rho(jd) ) * dpsi_dt ), & 
&                     dye_dtmax )
  END DO ! k = 1,nnugp(n)

  n                 = 2 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd,k), rho(jd), t(jd), xneut, xprot, xh, ah,    &
&    zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd,k)**2 * dunu(jd,k)
    dpsi_dt         = cvel * ( emitnp * ( one - psi0(jd,k,n) )            &
&                   - absornp * psi0(jd,k,n) )
    dye_dt          = dye_dt + ( ncoefap * rmu/rho(jd) ) * dpsi_dt
    dye_dtmax       = DMAX1( DABS( ( ncoefap * rmu/rho(jd) ) * dpsi_dt ), &
&                     dye_dtmax )
  END DO ! k = 1,nnugp(n)

  IF ( DABS(dye_dt) <= tol * dye_dtmax ) EXIT
  IF ( dye_dt <= zero ) THEN
    ye_max          = ye_test
  ELSE
    ye_min          = ye_test
  END IF ! dye_dt < 0

END DO ! it = 1,itmax
yeeq_jd             = ye_test

!-----------------------------------------------------------------------
!  Compute the kinetic equilibrium value of ye at radial zone jd - 1
!   with both emission and absorption
!-----------------------------------------------------------------------

ye_max              = ye(jd-1) + 0.02d0
ye_min              = ye(jd-1) - 0.02d0
DO it = 1,itmax

  ye_test           = half * ( ye_min + ye_max )
  dye_dt            = zero
  dye_dtmax         = zero
  CALL eqstta_x( 7, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, xneut, dum1, dum2, dum3 )
  CALL eqstta_x( 8, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, xprot, dum1, dum2, dum3 )
  CALL eqstta_x( 4, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, cmpn , dum1, dum2, dum3 )
  CALL eqstta_x( 5, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, cmpp , dum1, dum2, dum3 )
  CALL eqstta_x( 6, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, cmpe , dum1, dum2, dum3 )

  n                 = 1 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd-1,k), rho(jd-1), t(jd-1), xneut, xprot,      &
&    xh, ah, zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd-1,k)**2 * dunu(jd-1,k)
    dpsi_dt         = cvel * ( emitnp * ( one - psi0(jd-1,k,n) )          &
&                   - absornp * psi0(jd-1,k,n) )
    dye_dt          = dye_dt - ( ncoefap * rmu/rho(jd-1) ) * dpsi_dt
    dye_dtmax       = DMAX1( DABS( ( ncoefap*rmu/rho(jd) )*dpsi_dt ),     &
&                     dye_dtmax )
  END DO ! k = 1,nnugp(n)

  n                 = 2 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd-1,k), rho(jd-1), t(jd-1), xneut, xprot,      &
&    xh, ah, zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd-1,k)**2 * dunu(jd-1,k)
    dpsi_dt         = cvel * ( emitnp * ( one - psi0(jd-1,k,n) )          &
&                   - absornp * psi0(jd-1,k,n) )
    dye_dt          = dye_dt + ( ncoefap * rmu/rho(jd-1) ) * dpsi_dt
    dye_dtmax       = DMAX1( DABS( ( ncoefap * rmu/rho(jd) ) * dpsi_dt ), &
&                     dye_dtmax )
  END DO ! k = 1,nnugp(n)

  IF ( DABS(dye_dt) <= tol*dye_dtmax ) EXIT
  IF ( dye_dt <= zero ) THEN
    ye_max          = ye_test
  ELSE
    ye_min          = ye_test
  END IF ! dye_dt < 0

END DO ! it = 1,itmax
yeeq_jdm1           = ye_test

!-----------------------------------------------------------------------
!  Compute the kinetic equilibrium value of ye at radial zone jd due
!   to neutrino absorption only.
!-----------------------------------------------------------------------

ye_max              = ye(jd) + 0.02d0
ye_min              = ye(jd) - 0.02d0
DO it = 1,itmax

  ye_test           = half * ( ye_min + ye_max )
  dye_dt            = zero
  CALL eqstta_x( 7, jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, xneut, dum1, dum2, dum3)
  CALL eqstta_x( 8, jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, xprot, dum1, dum2, dum3)
  CALL eqstta_x( 4, jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, cmpn , dum1, dum2, dum3)
  CALL eqstta_x( 5, jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, cmpp , dum1, dum2, dum3)
  CALL eqstta_x( 6, jd, ij_ray, ik_ray, rho(jd), t(jd), ye_test, cmpe , dum1, dum2, dum3)

  n                 = 1 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd,k), rho(jd), t(jd), xneut, xprot, xh,      &
&    ah, zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd,k)**2 * dunu(jd,k)
    dpsi_dt         = cvel * ( - absornp * psi0(jd,k,n) )
    dye_dt          = dye_dt - ( ncoefap * rmu/rho(jd) ) * dpsi_dt
  END DO !k = 1,nnugp(n)

  n                 = 2 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd,k), rho(jd), t(jd), xneut, xprot, xh,      &
&    ah, zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd,k)**2 * dunu(jd,k)
    dpsi_dt         = cvel * ( - absornp * psi0(jd,k,n) )
    dye_dt          = dye_dt + ( ncoefap * rmu/rho(jd) ) * dpsi_dt
  END DO ! k = 1,nnugp(n)

  IF ( DABS(dye_dt) <= tol * ye(jd) ) EXIT
  IF ( dye_dt <= zero ) THEN
    ye_max          = ye_test
  ELSE
    ye_min          = ye_test
  END IF ! dye_dt

END DO ! it = 1,itmax
yeeqab_jd           = ye_test

!-----------------------------------------------------------------------
!  Compute the kinetic equilibrium value of ye at radial zone jd - 1
!   due to neutrino absorption only
!-----------------------------------------------------------------------

ye_max              = ye(jd-1) + 0.02d0
ye_min              = ye(jd-1) - 0.02d0
DO it = 1,itmax

  ye_test           = half * ( ye_min + ye_max )
  dye_dt            = zero
  CALL eqstta_x( 7, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, xneut ,dum1 ,dum2 ,dum3 )
  CALL eqstta_x( 8, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, xprot ,dum1 ,dum2 ,dum3 )
  CALL eqstta_x( 4, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, cmpn  ,dum1 ,dum2 ,dum3 )
  CALL eqstta_x( 5, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, cmpp  ,dum1 ,dum2 ,dum3 )
  CALL eqstta_x( 6, jd-1, ij_ray, ik_ray, rho(jd-1), t(jd-1), ye_test, cmpe  ,dum1, dum2 ,dum3 )

  n                 = 1 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd-1,k), rho(jd-1), t(jd-1), xneut, xprot, xh, &
&    ah, zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd-1,k)**2 * dunu(jd-1,k)
    dpsi_dt         = cvel * ( - absornp * psi0(jd-1,k,n) )
    dye_dt          = dye_dt - ( ncoefap*rmu/rho(jd-1) ) * dpsi_dt
  END DO ! k = 1,nnugp(n)

  n                 = 2 
  DO k = 1,nnugp(n)
    CALL abem_cal( n, unu(jd-1,k), rho(jd-1), t(jd-1), xneut, xprot, xh, &
&    ah, zh, cmpn, cmpp, cmpe, absornp, emitnp, ye_test )
    ncoefap         = ncoef * stwt(n) * unu(jd-1,k)**2 * dunu(jd-1,k)
    dpsi_dt         = cvel * ( - absornp * psi0(jd-1,k,n) )
    dye_dt          = dye_dt + ( ncoefap * rmu/rho(jd-1) ) * dpsi_dt
  END DO ! k = 1,nnugp(n)

  IF ( DABS(dye_dt) <= tol * ye(jd-1) ) EXIT
  IF ( dye_dt <= zero ) THEN
    ye_max          = ye_test
  ELSE
    ye_min          = ye_test
  END IF ! dye_dt

END DO ! it = 1,itmax
yeeqab_jdm1         = ye_test

!-----------------------------------------------------------------------
!  Interpolate variables in mass to msslag
!-----------------------------------------------------------------------

rstmssjd            = rstmss(jd  )/msolar
rstmssjdm1          = rstmss(jd-1)/msolar

tau_msslag          = rinterp( tau_jd   , tau_jdm1   , rstmssjd, msslag, rstmssjdm1 )
mdot_msslag         = rinterp( mdot_jd  , mdot_jdm1  , rstmssjd, msslag, rstmssjdm1 )
yeeq_msslag         = rinterp( yeeq_jd  , yeeq_jdm1  , rstmssjd, msslag, rstmssjdm1 )
yeeqab_msslag       = rinterp( yeeqab_jd, yeeqab_jdm1, rstmssjd, msslag, rstmssjdm1 )

!-----------------------------------------------------------------------
!  Write to nlagplt ("variable lagplot files")
!-----------------------------------------------------------------------

WRITE (lagfile,'(a25,i3.3,a2)') '/Plot_Files/lagplot2_',nlagdump,'.d'

!-----------------------------------------------------------------------
!  Open lagplot2.d file
!-----------------------------------------------------------------------

OPEN (UNIT=nlagplt,FILE=TRIM(lagfile),STATUS='new',IOSTAT=istat)
IF ( istat == 0 ) THEN
  WRITE (nlagplt,3) head
  WRITE (nlagplt,5)
  WRITE (nlagplt,7) time,msslag
  CALL date_and_time_print( n_time )
  WRITE (nlagplt,201)
  WRITE (nlagplt,203)
END IF ! istat == 0
      
IF ( istat .ne. 0 ) OPEN ( UNIT=nlagplt, FILE=lagfile, STATUS='old',    &
& POSITION='append', IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,8501)
  STOP
END IF ! istat /= 0

WRITE (nlagplt,205) time, tau_msslag, mdot_msslag, ye_msslag,           &
& yeeq_msslag, yeeqab_msslag

!-----------------------------------------------------------------------
!  Close lagplot2.d file
!-----------------------------------------------------------------------

CLOSE (unit=nlagplt,status='keep',err=9500)

RETURN

 9000 CONTINUE
WRITE (nprint,9001)
STOP

 9500 CONTINUE
WRITE (nprint,9501)
STOP

CONTAINS
REAL (KIND=double) FUNCTION rinterp(a,b,x,y,z)

REAL (KIND=double) :: a
REAL (KIND=double) :: b
REAL (KIND=double) :: x
REAL (KIND=double) :: y
REAL (KIND=double) :: z

rinterp             = b + ( a - b ) * ( y - z )/( x - z )
END FUNCTION rinterp

END SUBROUTINE lagrangeplot
