SUBROUTINE bnuplot( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         bnuplot
!    Module:       bnuplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/11/00
!
!    Purpose:
!      To dump luminosity, mean energy, shock, convection, and
!       boundary to files at selected intervals.
!
!    Subprograms called:
!      eqstta_x
!
!    Input arguments:
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!      kind_module, array_module, numerical_module, physcnst_module
!      abem_module, convect_module, cycle_module, edit_module,
!      eos_snc_x_module, mdl_cnfg_module, nu_dist_module,
!      nu_energy_grid_module, prb_cntl_module, scat_i_module,
!      shock_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nnu
USE numerical_module, ONLY : zero, half, one, epsilon, frpi, third
USE physcnst_module, ONLY : kmev, msolar

USE abem_module, ONLY : absor
USE convect_module, ONLY : ulcnvct
USE edit_module, ONLY : dtimeplot, iplotinnerb, rinnerb, nprint, nplotinnerb, &
& iplotouterb, routerb, nplotouterb, iplotlum, r_lumerms, nplotlum, iplotshk, &
& nplotshk, iplotcnv, nplotcnv, iplotmss, nplotmss, data_path
USE eos_snc_x_module, ONLY : aesv, nse
USE mdl_cnfg_module, ONLY : r, jr_min, jr_max, rstmss, u, rho, t, ye, dr, grvmss
USE nu_dist_module, ONLY : fluxnu, stwt, ecoefa, psi0, unu, psi1, dunujeadt, &
& ecoefae
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE scat_i_module, ONLY : scti
USE shock_module, ONLY : pq_x
USE t_cntrl_module, ONLY : time, dtnph

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                   :: ij_ray         ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                   :: ik_ray         ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                               :: itaut

INTEGER                               :: i              ! do index
INTEGER                               :: j              ! radial zone index
INTEGER                               :: ju             ! upper radial zone index for interpolation
INTEGER                               :: jd             ! lower radial zone index for interpolation
INTEGER                               :: k              ! neutrino energy zone index
INTEGER                               :: n              ! neutrino flavor zone index
INTEGER                               :: ic             ! convective boundary insex
INTEGER                               :: itime          ! integer of ( time * tmult )
INTEGER                               :: itimeprev      ! previous value of int( time * tmult )
INTEGER                               :: istat          ! open - close flag

INTEGER                               :: jsph           ! radial zone just outside Planck averaged neutrinospehre
INTEGER                               :: j0             ! radial index just outside shock
INTEGER                               :: jc_max         ! used for determining the MLT convective boundaries
INTEGER                               :: jc_max1        ! used for determining the MLT convective boundaries

REAL(KIND=double)                     :: tmult          ! used to determine when to write a file
REAL(KIND=double)                     :: rsave          ! saved value of r(1)

REAL(KIND=double)                     :: rjmhlu         ! log of rjmh at ju
REAL(KIND=double)                     :: rjmhld         ! log of rjmh at jd
REAL(KIND=double)                     :: rbdaryl        ! log of rinnerb

REAL(KIND=double)                     :: rholu          ! log of rho at ju
REAL(KIND=double)                     :: rhold          ! log of rho at jd
REAL(KIND=double)                     :: robdaryl       ! log of rho interpolated to rinnerb
REAL(KIND=double)                     :: robdary        ! rho interpolated to rinnerb

REAL(KIND=double)                     :: tmev           ! temperature (MeV)
REAL(KIND=double)                     :: tlu            ! log of temperature at ju
REAL(KIND=double)                     :: tld            ! log of temperature at jd
REAL(KIND=double)                     :: tbdaryl        ! log of temperature interpolated to rinnerb
REAL(KIND=double)                     :: tbdary         ! temperature interpolated to rinnerb

REAL(KIND=double)                     :: yelu           ! log of electron fraction at ju
REAL(KIND=double)                     :: yeld           ! log of electron fraction at jd
REAL(KIND=double)                     :: yebdaryl       ! log of electron fraction interpolated to rinnerb
REAL(KIND=double)                     :: yebdary        ! electron fraction interpolated to rinnerb

REAL(KIND=double)                     :: fluxcjd        ! neutrino luminosity coefficient at jd
REAL(KIND=double)                     :: fluxcjdm       ! neutrino luminosity coefficient at jd - 1
REAL(KIND=double)                     :: r_lumjd        ! luminosity at jd
REAL(KIND=double)                     :: r_lumjdm       ! luminosity at jd - 1
REAL(KIND=double)                     :: taujk          ! neutrino mean free path as a function of j,k
REAL(KIND=double)                     :: taut           ! Planck averaged neutrino mean free path at j
REAL(KIND=double)                     :: tautp          ! Planck averaged neutrino mean free path at j+1
REAL(KIND=double)                     :: ymdavt         ! used for  computing Planck averaged neutrino mean free path
REAL(KIND=double)                     :: ymdave         ! used for  computing Planck averaged neutrino mean free path
REAL(KIND=double)                     :: flux           ! neutrino flux
REAL(KIND=double)                     :: remsf          ! Planck averaged neutrinosphere radius
REAL(KIND=double), DIMENSION(nnu)     :: rnusph         ! Planck averaged neutrinosphere radius

REAL(KIND=double)                     :: rstmssj        ! interpolated enclosed rest mass (g)
REAL(KIND=double)                     :: uj             ! interpolated velocity (cm s^{-1})

REAL(KIND=double)                     :: pp             ! pressure at robdary, tbdary, yebdary (dynes cm^{-2})
REAL(KIND=double)                     :: ee             ! pressure at robdary, tbdary, yebdary (ergs g^{-1})
REAL(KIND=double)                     :: ss             ! entropy at robdary, tbdary, yebdary
REAL(KIND=double)                     :: u1             ! dummy variable
REAL(KIND=double)                     :: u2             ! dummy variable
REAL(KIND=double)                     :: u3             ! dummy variable

REAL(KIND=double)                     :: coefe2         ! e4 de at jd
REAL(KIND=double)                     :: coefe2m        ! e4 de at jd - 1
REAL(KIND=double)                     :: dennjd         ! integral of psi0 e2 de at jd
REAL(KIND=double)                     :: dennjdm        ! integral of psi0 e2 de at jd - 1
REAL(KIND=double)                     :: dene2jd        ! integral of psi0 e4 de at jd
REAL(KIND=double)                     :: dene2jdm       ! integral of psi0 e4 de at jd - 1
REAL(KIND=double)                     :: enuv2jd        ! n neutrinmo rms energy at jd
REAL(KIND=double)                     :: enuv2jdm       ! n neutrinmo rms energy at jd - 1
REAL(KIND=double)                     :: enuv2          ! n neutrinmo rms energy at r_lumerms

REAL(KIND=double), DIMENSION(nx)      :: rjmh           ! zone centered radii
REAL(KIND=double), DIMENSION(nnu)     :: lum_r          ! n neutrino luminosity at radius r_lumerms
REAL(KIND=double), DIMENSION(nnu)     :: enuvrms        ! n neutrino rms energy at radius r_lumerms

REAL(KIND=double), PARAMETER          :: pqcrit = 0.2d0 ! minimum pq_x/p for presence of shock

REAL(KIND=double), DIMENSION(4)       :: rconvmn        ! inner boundary of a convective region
REAL(KIND=double), DIMENSION(4)       :: rconvmx        ! outer boundary of a convective region

REAL(KIND=double)                     :: dmdt           ! mass accretion rate into shock (g s^{-1}}
REAL(KIND=double)                     :: shkmassgr      ! gravitational mass accretion rate into shock (M_solar s^{-1}}
REAL(KIND=double)                     :: shkmassnewt    ! rest mass accretion rate into shock (M_solar s^{-1}}

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (a1)
  201 FORMAT (es15.8,9(es12.4))
  301 FORMAT (es15.8,12(es12.4))
  401 FORMAT (es15.8,10(es12.4))
  501 FORMAT (es15.8,8(es12.4))
  601 FORMAT (es15.8,3(es12.4))
 1001 FORMAT (' ju cannot be found in subroutine bnuplot for rinnerb=',es10.3,' for r(ju)')
 2001 FORMAT (' ju cannot be found in subroutine bnuplot for rinnerb=',es10.3,' for rjmh(ju)')
 3001 FORMAT (' ju cannot be found in subroutine bnuplot for routerb=',es10.3,' for r(ju)')
 4001 FORMAT (' ju cannot be found in subroutine bnuplot for routerb=',es10.3,' for rjmh(ju)')
 5001 FORMAT (' jd cannot be found in subroutine bnuplot for r_lumerms=', es10.3,' for r(jd)')
 5501 FORMAT (' jd cannot be found in subroutine bnuplot for r_lumerms=', es10.3,' for rjmh(jd)')
 8001 FORMAT (' File nplotinnerb cannot be opened in subroutime bnuplot')
 8101 FORMAT (' File nplotouterb cannot be opened in subroutime bnuplot')
 8201 FORMAT (' File nplotlum cannot be opened in subroutime bnuplot')
 8301 FORMAT (' File nplotshk cannot be opened in subroutime bnuplot')
 8401 FORMAT (' File nplotcnv cannot be opened in subroutime bnuplot')
 8501 FORMAT (' File nplotmss cannot be opened in subroutime bnuplot')

!-----------------------------------------------------------------------
!        Evaluate criterion for writing to plot files.
!-----------------------------------------------------------------------

tmult              = 1.d+3/dtimeplot
itime              = INT( time * tmult )
itimeprev          = INT( ( time - dtnph ) * tmult )
if ( itime == itimeprev ) RETURN

!-----------------------------------------------------------------------
!        Compute the zone-centered radii.
!-----------------------------------------------------------------------

rsave              = r(1)
r(1)               = epsilon
DO j = 2,jr_max
  rjmh(j)          = ( half * ( r(j)**3 + r(j-1)**3 ) )**third
END DO
r(1)               = rsave

!-----------------------------------------------------------------------
!
!                     \\\\\ INNER BOUNDARY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Evaluate criterion for writing entries to the inner
!         boundary file.
!-----------------------------------------------------------------------

IF ( iplotinnerb > 0 ) THEN

!-----------------------------------------------------------------------
!        Interpolate zone-edged quantities to the inner boundary.
!-----------------------------------------------------------------------

  ju               = 0
  jd               = 0
  DO j = 2,jr_max
    IF ( r(j) > rinnerb ) THEN
      ju           = j
      jd           = j - 1
      EXIT
    END IF ! r(j) > rinnerb
  END DO

  IF ( ju == 0 ) THEN
    WRITE (nprint,1001) rinnerb
    GO TO 2900
  END IF ! ju == 0

  rstmssj          = rinterp( rstmss(ju), rstmss(jd), r(ju), rinnerb, r(jd) )
  uj               = rinterp(u(ju)      , u(jd)     , r(ju), rinnerb, r(jd) )

!-----------------------------------------------------------------------
!        Interpolate zone-centered quantities to the inner boundary.
!-----------------------------------------------------------------------

  ju               = 0
  jd               = 0
  DO j = 3,jr_max
    IF ( rjmh(j) > rinnerb ) THEN
      ju           = j
      jd           = j - 1
      EXIT
    END IF ! rjmh(j) > rinnerb
  END DO
  
  IF ( ju == 0 ) THEN
     WRITE (nprint,2001) rinnerb
     GO TO 2900
  END IF ! ju == 0

  rjmhlu           = DLOG10( rjmh(ju) )
  rjmhld           = DLOG10( rjmh(jd) )
  rbdaryl          = DLOG10( rinnerb )

  rholu            = DLOG10( rho(ju) )
  rhold            = DLOG10( rho(jd) )
  robdaryl         = rinterp( rholu, rhold, rjmhlu, rbdaryl, rjmhld )
  robdary          = 10.d+00**robdaryl

  tlu              = DLOG10(   t(ju) )
  tld              = DLOG10(   t(jd) )
  tbdaryl          = rinterp(  tlu,  tld,rjmhlu,rbdaryl,rjmhld)
  tbdary           = 10.d+00** tbdaryl

  yelu             = DLOG10(  ye(ju) )
  yeld             = DLOG10(  ye(jd) )
  yebdaryl         = rinterp( yelu, yeld,rjmhlu,rbdaryl,rjmhld)
  yebdary          = 10.d+00**yebdaryl

!-----------------------------------------------------------------------
!        Compute the boundary pressure, energy, and entropy.
!-----------------------------------------------------------------------

  CALL eqstta_x( 1, j, ij_ray, ik_ray, robdary, tbdary, yebdary, pp, u1, &
&  u2, u3 )
  CALL eqstta_x( 2, j, ij_ray, ik_ray, robdary, tbdary, yebdary, ee, u1, &
&  u2, u3 )
  CALL eqstta_x( 3, j, ij_ray, ik_ray, robdary, tbdary, yebdary, ss, u1, &
&  u2, u3 )
  
!-----------------------------------------------------------------------
!        Open inner boundary file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nplotinnerb,FILE=TRIM(data_path)//'/Plot_Files/ibound.d', STATUS='new', &
&  IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=nplotinnerb, FILE=TRIM(data_path)//'/Plot_Files/ibound.d', &
& STATUS='old',POSITION='append', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8001)
    STOP
  END IF

!-----------------------------------------------------------------------
!        Write to inner boundary file.
!-----------------------------------------------------------------------

  tmev             = kmev * tbdary
  WRITE (nplotinnerb,201) time, rinnerb, rstmssj, uj, robdary, tmev, &
& yebdary, pp, ee, ss

!-----------------------------------------------------------------------
!        Close inner boundary file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nplotinnerb,STATUS='keep')

END IF ! iplotinnerb > 0

!-----------------------------------------------------------------------
!
!                     \\\\\ OUTER BOUNDARY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Evaluate criterion for writing entries to the outer
!         boundary file.
!-----------------------------------------------------------------------

 2900 CONTINUE
IF ( iplotouterb > 0 ) THEN

!-----------------------------------------------------------------------
!        Interpolate zone-edged quantities to the outer boundary.
!-----------------------------------------------------------------------

  ju               = 0
  jd               = 0
  DO j = 2,jr_max
    IF ( r(j) > routerb ) THEN
      ju           = j
      jd           = j - 1
      EXIT
    END IF ! r(j) > rinnerb
  END DO ! j = 2,jr_max

  IF ( ju == 0 ) THEN
    WRITE (nprint,3001) routerb
    GO TO 4900
  END IF ! ju == 0

  rstmssj          = rinterp( rstmss(ju), rstmss(jd), r(ju), routerb, r(jd) )
  uj               = rinterp( u(ju)     , u(jd)     , r(ju), routerb, r(jd) )

!----------------------------------------------------------------------!
!        Interpolate zone-centered quantities to the outer boundary.   !
!----------------------------------------------------------------------!

  ju               = 0
  jd               = 0
  DO j = 3,jr_max
    IF ( rjmh(j) > routerb ) THEN
      ju           = j
      jd           = j - 1
      EXIT
    END IF ! rjmh(j) > routerb
  END DO ! 3,jr_max
  
  IF ( ju == 0 ) THEN
     WRITE (nprint,4001) routerb
     GO TO 4900
  END IF ! ju == 0

  rjmhlu           = DLOG10( rjmh(ju) )
  rjmhld           = DLOG10( rjmh(jd) )
  rbdaryl          = DLOG10( routerb )

  rholu            = DLOG10( rho(ju) )
  rhold            = DLOG10( rho(jd) )
  robdaryl         = rinterp( rholu, rhold, rjmhlu, rbdaryl, rjmhld )
  robdary          = 10.d+00**robdaryl

  tlu              = DLOG10(   t(ju) )
  tld              = DLOG10(   t(jd) )
  tbdaryl          = rinterp(  tlu,  tld,rjmhlu,rbdaryl,rjmhld)
  tbdary           = 10.d+00** tbdaryl

  yelu             = DLOG10(  ye(ju) )
  yeld             = DLOG10(  ye(jd) )
  yebdaryl         = rinterp( yelu, yeld,rjmhlu,rbdaryl,rjmhld)
  yebdary          = 10.d+00**yebdaryl

!-----------------------------------------------------------------------
!        Compute the boundary pressure, energy, and entropy.
!-----------------------------------------------------------------------

  CALL eqstta_x( 1, j, ij_ray, ik_ray, robdary, tbdary, yebdary, pp, u1, &
&  u2, u3 )
  CALL eqstta_x( 2, j, ij_ray, ik_ray, robdary, tbdary, yebdary, ee, u1, &
&  u2, u3 )
  CALL eqstta_x( 3, j, ij_ray, ik_ray, robdary, tbdary, yebdary, ss, u1, &
&  u2, u3 )

!----------------------------------------------------------------------!
!        Open outer boundary file.                                     !
!----------------------------------------------------------------------!

  OPEN (UNIT=nplotouterb,FILE=TRIM(data_path)//'/Plot_Files/obound.d', STATUS='new', &
&  IOSTAT=istat)
  IF ( istat .ne. 0 ) OPEN (UNIT=nplotouterb, FILE=TRIM(data_path)//'/Plot_Files/obound.d', &
&  STATUS='old', POSITION='append', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8101)
    STOP
  END IF

!-----------------------------------------------------------------------
!        Write to outer boundary file.
!-----------------------------------------------------------------------

  tmev             = kmev * tbdary
  WRITE (nplotouterb,201) time, routerb, rstmssj, uj, robdary, tmev, &
& yebdary, pp, ee, ss

!-----------------------------------------------------------------------
!        Close outer boundary file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nplotouterb, STATUS='keep')

END IF ! iplotouterb > 0

!-----------------------------------------------------------------------
!
!                   \\\\\ LUMINOSITY BOUNDARY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Evaluate criterion for writing entries to luminosity
!         boundary file.
!-----------------------------------------------------------------------

 4900 CONTINUE
 
IF ( nnugpmx /= 0  .and. iplotlum < 0 ) THEN

!-----------------------------------------------------------------------
!        Neutrino luminosities, rms energies, and emission surfaces.
!
!        Find jd such that r(jd  ) > r_lumerms
!                    and   r(jd-1) < r_lumerms.
!----------------------------------------------------------------------!

  jd               = 0
  DO j = 2,jr_max
    IF ( r(j) >= r_lumerms ) THEN
      jd           = j
      EXIT
    END IF ! r(j) > r_lumerms
  END DO

  IF ( jd == 0 ) THEN
    WRITE (nprint,5001) r_lumerms
    GO TO 5500
  END IF

!-----------------------------------------------------------------------
!        Compute the n-neutrino luminosities at radius r_lumerms.
!-----------------------------------------------------------------------

  fluxcjd     = frpi * r(jd  ) * r(jd  ) * 1.d-30
  fluxcjdm    = frpi * r(jd-1) * r(jd-1) * 1.d-30
  DO  n = 1,nnu
    r_lumjd     = fluxnu(jd  ,n) * fluxcjd  * 1.d-2/stwt(n)
    r_lumjdm    = fluxnu(jd-1,n) * fluxcjdm * 1.d-2/stwt(n)
    lum_r(n)    = rinterp( r_lumjd, r_lumjdm, r(jd), r_lumerms, r(jd-1) )
  END DO

!----------------------------------------------------------------------!
!        Find jd such that rjmh(jd  ) > r_lumerms                       !
!                    and   rjmh(jd-1) < r_lumerms.                      !
!----------------------------------------------------------------------!

 5500  CONTINUE

  jd               = 0
  DO j = 2,jr_max
    IF ( rjmh(j) > r_lumerms ) THEN
      jd           = j
      EXIT
    END IF ! rjmh(j) > r_lumerms
  END DO

  IF ( jd == 0 ) THEN
    WRITE (nprint,5501) r_lumerms
    GO TO 6900
  END IF

!-----------------------------------------------------------------------
!
!                   \\\\\ RMS ENERGY BOUNDARY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Compute the n-neutrino rms energies at radius r_lumerms.
!
!          <erms> = int( psi0*unu**5 )/int( psi0*unu**3 )
!-----------------------------------------------------------------------
 
   DO n = 1,nnu

    dennjd         = zero
    dene2jd        = zero
    dennjdm        = zero
    dene2jdm       = zero

    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)
      coefe2       = ecoefa(jd  ,k) * unu(jd  ,k) * unu(jd  ,k)
      coefe2m      = ecoefa(jd-1,k) * unu(jd-1,k) * unu(jd-1,k)
      dennjd       = dennjd    + ecoefa(jd  ,k) * psi0(jd  ,k,n)
      dennjdm      = dennjdm   + ecoefa(jd-1,k) * psi0(jd-1,k,n)
      dene2jd      = dene2jd   + coefe2  * psi0(jd  ,k,n)
      dene2jdm     = dene2jdm  + coefe2m * psi0(jd-1,k,n)
    END DO

    enuv2jd        = dene2jd /( dennjd   + epsilon )
    enuv2jdm       = dene2jdm/( dennjdm  + epsilon )
    enuv2          = rinterp( enuv2jd, enuv2jdm, rjmh(jd), r_lumerms, rjmh(jd-1) )
    enuvrms(n)     = dsqrt( dabs(enuv2) + epsilon )

  END DO

!-----------------------------------------------------------------------
!        Plank averaged emission surfaces.
!-----------------------------------------------------------------------

 6900  CONTINUE

  DO n = 1,nnu

    taut           = zero
    itaut          = .false.

    DO j = jr_max,2,-1

      tautp        = taut
      ymdavt       = zero
      flux         = zero

      DO k = 1,nnugpmx
        taujk      = ( absor(j,k,n,ij_ray,ik_ray) + scti(j,k,n) ) * dr(j)
        ymdavt     = ymdavt + taujk * psi1(j,k,n) * ecoefae(j,k)
        flux       = flux   +         psi1(j,k,n) * ecoefae(j,k)
      END DO !  k = 1,nnugpmx

      ymdave       = ymdavt/( flux + epsilon )
      taut         = taut + ymdave

      IF ( taut >= one ) THEN
        itaut      = .true.
        jsph       = j
        EXIT
      END IF ! taut > 1

    END DO ! j = jr_max,2,-1

    IF ( itaut ) THEN
      IF ( jsph == 2 ) THEN
        remsf      = r(2)
      ELSE ! jsph /= 2
        remsf      = r(jsph) + ( r(jsph-1) - r(jsph) ) &
&                  * ( one - tautp )/( taut - tautp + 1.d-20 )
      END IF ! jsph = 2
    ELSE ! itaut = false
      remsf        = zero
    END IF ! itaut

    rnusph(n)      = remsf

  END DO ! n = 1,nnu

!-----------------------------------------------------------------------
!        Open luminosity boundary file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nplotlum, FILE=TRIM(data_path)//'/Plot_Files/lbound.d', STATUS='new', IOSTAT=istat )
  IF ( istat /= 0 ) OPEN (UNIT=nplotlum, FILE=TRIM(data_path)//'/Plot_Files/lbound.d', &
& STATUS='old', POSITION='append', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8201)
    STOP
  END IF ! istat /= 0
  
!-----------------------------------------------------------------------
!        Write to luminosity boundary file.
!-----------------------------------------------------------------------

  WRITE (nplotlum,301) time, (lum_r(n), n = 1,nnu), (enuvrms(n), n = 1,nnu), &
& (rnusph(n), n = 1,nnu)

!-----------------------------------------------------------------------
!        Close luminosity boundary file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nplotlum, STATUS='keep')

END IF ! iplotlum > 0

!-----------------------------------------------------------------------
!
!                 \\\\\ MLT CONVECTIVE BOUNDARIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Evaluate criterion for writing entries to convection file.
!-----------------------------------------------------------------------

IF ( iplotcnv > 0 ) THEN

!-----------------------------------------------------------------------
!        Boundaries of active convection zones.
!-----------------------------------------------------------------------

  rconvmn          = zero
  rconvmx          = zero

!........Determine outer boundary

  jc_max           = jr_max

  jc_max1          = jr_max
  c_outer: DO ic = 1,4
    DO j = jc_max-1,2,-1
      IF ( ulcnvct(j) > zero  .and.  ulcnvct(j+1) == zero ) THEN
         rconvmx(ic) = r(j)
         jc_max1   = j
         EXIT
      END IF ! ulcnvct(j) > 0, ulcnvct(j+1) = 0
    END DO ! j = jr_max-1,2,-1

    IF ( jc_max1 == jc_max ) EXIT c_outer
    jc_max         = jc_max1

!........Determine inner boundary

    DO j = jc_max-1,1,-1
      IF ( ulcnvct(j) == zero  .and.  ulcnvct(j+1) > zero ) THEN
        rconvmn(ic) = r(j)
        jc_max1    = j
        EXIT
      END IF ! ulcnvct(j) = 0, ulcnvct(jj+1) > 0
    END DO ! j = jr_max-1,1,-1

    IF ( jc_max1 == jc_max ) EXIT c_outer
    jc_max         = jc_max1

  END DO c_outer ! ic = 1,4

!-----------------------------------------------------------------------
!        Open convection file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nplotcnv, FILE=TRIM(data_path)//'/Plot_Files/conv.d', STATUS='new', IOSTAT=istat)
  IF ( istat .ne. 0 ) OPEN (UNIT=nplotcnv, FILE=TRIM(data_path)//'/Plot_Files/conv.d', &
&  STATUS='old', POSITION='append', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8401)
    STOP
  END IF ! istat /= 0
!-----------------------------------------------------------------------
!        Write to convection file.
!-----------------------------------------------------------------------

  WRITE (nplotcnv,501) time,(rconvmn(i),rconvmx(i),i = 1,4)

!-----------------------------------------------------------------------
!        Close convection file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nplotcnv, STATUS='keep')

END IF ! iplotcnv > 0

!-----------------------------------------------------------------------
!
!                       \\\\\ MASS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Evaluate criterion for writing entries to mass file.
!-----------------------------------------------------------------------

IF ( iplotmss > 0 ) THEN

!-----------------------------------------------------------------------
!        Mass accretion rate and mass below shock.
!-----------------------------------------------------------------------

  shkmassnewt      = zero
  shkmassgr        = zero
  dmdt             = zero

  DO j = jr_max-1,2,-1
    j0             = j + 1
    IF ( pq_x(j,ij_ray,ik_ray)/aesv(j,1,ij_ray,ik_ray) >= pqcrit ) EXIT
  END DO ! j = jr_max-1,2,-1

  IF ( j0 /= 3 ) THEN
    dmdt           = frpi * r(j0)**2 * u(j0) * half * ( rho(j0) + rho(j0+1) )/msolar
    shkmassnewt    = rstmss(j0)/msolar
    shkmassgr      = grvmss(j0)/msolar
  END IF ! j0 /= 3

!-----------------------------------------------------------------------
!        Open mass file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nplotmss, FILE=TRIM(data_path)//'/Plot_Files/mass.d', STATUS='new', IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=nplotmss, FILE=TRIM(data_path)//'/Plot_Files/mass.d', &
&  STATUS='old', POSITION='append', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8501)
    STOP
  END IF

!-----------------------------------------------------------------------
!        Write to mass file.
!-----------------------------------------------------------------------

  WRITE (nplotmss,601) time, dmdt, shkmassnewt, shkmassgr

!-----------------------------------------------------------------------
!        Close mass file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nplotmss, STATUS='keep')

END IF ! iplotmss > 0

!-----------------------------------------------------------------------
!        Return.
!-----------------------------------------------------------------------

RETURN


CONTAINS
REAL (KIND=double) FUNCTION rinterp(a,b,x,y,z)

REAL (KIND=double) :: a
REAL (KIND=double) :: b
REAL (KIND=double) :: x
REAL (KIND=double) :: y
REAL (KIND=double) :: z

rinterp            = b + ( a - b ) * ( y - z )/( x - z )
END FUNCTION rinterp

END
