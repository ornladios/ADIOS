SUBROUTINE eqstta_x( i, j, ij_ray, ik_ray, rhoj, tj, yej, v, vd, vt, vy )
!-----------------------------------------------------------------------
!
!    File:         eqstlsa_x
!    Module:       eqstta_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/15/03
!
!    Purpose:
!      To call the Lattimer-Swesty equation of state directly. The
!       equation of state quantities returned are obtained directly
!       in contradistinction to being interpolaed from a grid set
!       up around the point in question. The Cooperstein-BCK
!       equation of state is called if the independent variables
!       lie outside the range permitted by the Lattimer-Swesty
!       equation of state.
!
!    Variables that must be passed through common:
!  nse(j,i_ray) : nuclear statistical equilibrium flag for radial zone j.
!  gamhv        :  BCK gamma.
!  wnm          :  BCK nucl matter energy; - 16.0 MeV.
!  ws           :  BCK bulk surface coefficient.
!  xk0          :  BCK K_0(x=.5) symmettric matter.
!  xkzafac      :  BCK drop in K_0 with x :
!              k(x) = k(1/2)(1-xkzafac*(1-2x)**2).
!
!    Subprograms called:
!  loadmx       : loads eos data into the LS eos
!  inveos       : invokes the LS eos
!  eos          : invokes the BCK eos
!
!    Input arguments:
!  i            : thermodynamic function index.
!   i = 1       : pressure.
!   i = 2       : energy.
!   i = 3       : entropy.
!   i = 4       : neutron chemical potential.
!   i = 5       : proton chemical potential.
!   i = 6       : electron chemical potential.
!   i = 7       : free neutron mass fraction.
!   i = 8       : free proton mass fraction.
!   i = 9       : heavy nucleus mass fraction.
!   i = 10      : heavy nucleus mass number.
!   i = 11      : heavy nucleus charge number.
!   i = 12      : gamma1.
!   i = 13      : equation of state quantities are obtained from the
!                  Cooperstein-BCK equation of state.                  !
!  j            : radial zone index.
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  rhoj         : matter density (g/cm**3).
!  tj           : matter temperature (K).
!  yej          : matter electron fraction.
!
!    Output arguments:
!  v            : value of thermodynamic variable i.
!  vd           : derivative of v wrt rhoj.
!  vt           : derivative of v wrt tj.
!  vy           : derivative of v wrt yej.
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, el_eos_module, eos_ls_module, eos_m4c_module,
!  eos_bck_module, eos_snc_x_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx
USE numerical_module, ONLY: zero, half, one
USE physcnst_module, ONLY: dmnp, kmev, rmu, cm3fm3, ergmev

USE edit_module, ONLY: nprint, nlog
USE el_eos_module, ONLY: MUSUBE, ES, PS, EPRESS, EU
USE eos_ls_module, ONLY : inpvarsa
USE eos_m4c_module, ONLY: inpvar, iflag, forflg, brydns, eosflg, xprev, PTOT, &
& UTOT, STOT, XNUT, XPROT, XH, MUN, MUPROT, A, X, GAM_S, BS, BSNUC, DPDN, DPDT, &
& DPDY, DUDN, DUDT, DUDY, DSDN, DSDT, DSDY
USE eos_bck_module, ONLY: jshel, dbck, tbck, yebck, dtran, ue, xnbck, xpbck, xhbck, ptotbck, uhat, etot,   &
& stotbck, un, zabck, ahbck, uea, una, uhata, thetaa, zaa, xaa, dtrana, theta, xabck, se, sd, sh, sneu, &
& ueaa, unaa, uhataa, thetaaa, zaaa, xaaa, dtranaa, ee_m=>ee, pe_m=>pe
USE eos_snc_x_module, ONLY: eos_i, nse, eosrho

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: i             ! eos dependent variable index
INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)    :: rhoj          ! density (g cm^{-3})
REAL(KIND=double), INTENT(in)    :: tj            ! temperature (K)
REAL(KIND=double), INTENT(in)    :: yej           ! electron fraction

!-----------------------------------------------------------------------
!        Out variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: v             ! i'th EOS variable
REAL(KIND=double), INTENT(out)   :: vd            ! derivative wrt demsity of the i'th EOS variable
REAL(KIND=double), INTENT(out)   :: vt            ! derivative wrt temperature of the i'th EOS variable
REAL(KIND=double), INTENT(out)   :: vy            ! derivative wrt electron fraction of the i'th EOS variable

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER               :: nxd = 5000

CHARACTER (LEN=1)                :: eos_ii

LOGICAL                          :: first = .true.
LOGICAL                          :: firstbck = .true.
LOGICAL, DIMENSION(nxd)          :: firstls

INTEGER                          :: ii            ! do index
INTEGER                          :: jj            ! radial zone do index
INTEGER                          :: sf            ! eos flag

REAL(KIND=double), PARAMETER     :: eost = 0.05d0 ! lower boundary of LS EOS table
REAL(KIND=double), PARAMETER     :: UTOT0 = 8.9d0 ! change in the zero of energy (MeV)
REAL(KIND=double), PARAMETER     :: delta = 0.01d+00 ! increment used for computing derivatives

REAL(KIND=double), DIMENSION(6)  :: rhod          ! density array for computing derivatives
REAL(KIND=double), DIMENSION(6)  :: td            ! temperature array for computing derivatives
REAL(KIND=double), DIMENSION(6)  :: yed           ! electron fraction array for computing derivatives
REAL(KIND=double), DIMENSION(6)  :: pderiv        ! pressure array for computing derivatives
REAL(KIND=double), DIMENSION(6)  :: uderiv        ! energy array for computing derivatives

REAL(KIND=double)                :: kfm           ! ( # nucleons/gram )( cm3/fm3 )
REAL(KIND=double)                :: kp            ! ( erg/cm3 ) / ( mev/fm3 )
REAL(KIND=double)                :: ku            ! ( # nucleons/gram )( erg/MeV )
REAL(KIND=double)                :: tmev          ! temperature (MeV)

REAL(KIND=double)                :: pprev         ! previous value of proton density
REAL(KIND=double)                :: t_old         ! old value of temperature

REAL(KIND=double)                :: pe            ! electron pressure (dynes cm^{-2})
REAL(KIND=double)                :: ee            ! electron energy (ergs cm^{-3})
REAL(KIND=double)                :: cmpe          ! electron chemical potential (MeV)
REAL(KIND=double)                :: yeplus        ! positron fraction
REAL(KIND=double)                :: rel           ! relativity parameter

REAL(KIND=double)                :: ye_inveos     ! min( ye, 0.5 )

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT ('Equation of state identifier does not match either L or B')
  201 FORMAT (' nx=',i3,' exceeds nxd=',i3,' in subroutine eqstlsa')
 1001 FORMAT (' sf = 0 in eqstta_x; j=',i3,' ij_ray,ik_ray=',i3,        &
& ' rhod=',es10.3,' td=',es10.3,' yed=',es10.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Test that nxd > nx
!-----------------------------------------------------------------------

IF ( nx > nxd ) THEN
  WRITE (nprint,201) nx,nxd
  WRITE (nlog,201) nx,nxd
  STOP
END IF

!-----------------------------------------------------------------------
!
!               \\\\\ DETERMINE WHICH EOS TO USE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Examine equation of state identifier
!-----------------------------------------------------------------------

IF ( eos_i /= 'L'  .and.  eos_i /= 'B'  .and.  eos_i /= 'S' ) THEN
  WRITE (nprint,101)
  STOP
END IF

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN

  kfm              = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
  kp               = ergmev/cm3fm3 ! ( erg/cm3 ) / ( mev/fm3 )
  ku               = ergmev/rmu    ! ( # nucleons/gram )( erg/mev )
  iflag            = 1
  forflg           = 0
  CALL loadmx()

  DO ii = 1,nx
    firstls(ii)    = .true.
  END DO

  first            = .false.

END IF

!-----------------------------------------------------------------------
!  Convert to independent variables of the Lattimer-Swesty and
!   Cooperstein BCK equation of state
!-----------------------------------------------------------------------

brydns             = rhoj * kfm
tmev               = tj * kmev

!-----------------------------------------------------------------------
!  If     brydns       > eosrho         and
!         tmev         > eost           and
!         nse(j,ij_ray,ik_ray) = 1      and
!         i            ne 12
!  then            use Lattimer-Swesty eos
!  otherwise,      use Cooperstein bck eos
!-----------------------------------------------------------------------

IF ( brydns >= eosrho           .and.  &
&    tmev >= eost               .and.  &
&    nse(j,ij_ray,ik_ray) == 1  .and.  &
&    i /= 12                    .and.  &
&    eos_i /= 'S' )                     THEN
  eos_ii           = 'L'
ELSE
  eos_ii           = 'B'
END IF ! brydns >= eosrho, tmev >= eost, nse(j,ij_ray,ik_ray) = 1, i /= 12, eos_i /= 'S'

!-----------------------------------------------------------------------
!
!                  \\\\\ LATTIMER-SWEATY EOS /////
!
!-----------------------------------------------------------------------

IF ( eos_ii == 'L' ) THEN

!-----------------------------------------------------------------------
!  Load Lattimer-Swesty equation of state parameters
!-----------------------------------------------------------------------

  IF ( firstls(j) ) THEN
    inpvarsa(j,1)  = yej * brydns
    inpvarsa(j,2)  = 0.155d+00
    inpvarsa(j,3)  = -15.d+00
    inpvarsa(j,4)  = -10.d+00
    firstls(j)     = .false.
  END IF ! firstls(j)

!-----------------------------------------------------------------------
!  Call the Lattimer-Swesty equation of state
!-----------------------------------------------------------------------

  pprev            = DMIN1( inpvarsa(j,1), brydns * yej )
  inpvar(1)        = tmev
  inpvar(2)        = inpvarsa(j,2)
  inpvar(3)        = inpvarsa(j,3)
  inpvar(4)        = inpvarsa(j,4)

  ye_inveos        = DMIN1( yej, 0.5d0 )

  CALL inveos( inpvar, t_old, ye_inveos, brydns, iflag, eosflg, forflg, &
&  sf, xprev, pprev )
  IF (sf == 0 ) WRITE (nprint,1001) j, ij_ray, ik_ray, rhoj, tj, ye_inveos

!-----------------------------------------------------------------------
!  Store equation of state parameters for zone j
!-----------------------------------------------------------------------

  inpvarsa(j,1)    = pprev
  inpvarsa(j,2)    = inpvar(2)
  inpvarsa(j,3)    = inpvar(3)
  inpvarsa(j,4)    = inpvar(4)

!-----------------------------------------------------------------------
!  Load individual contributions to the entropy
!-----------------------------------------------------------------------

  stotbck          = STOT
  se               = ES
  sneu             = PS
  sd               = BS - BSNUC
  sh               = BSNUC

!-----------------------------------------------------------------------
!  Compute electron equation of state
!-----------------------------------------------------------------------

  CALL e_p_eos( brydns, tmev, yej, pe, ee, se, cmpe, yeplus, rel )

  pe_m             = pe
  ee_m             = ee

!-----------------------------------------------------------------------
!  Load v, vd, vt, and vr, changing from MeV and fm to cgs units, as
!   appropriate.
!-----------------------------------------------------------------------

  i_LS: SELECT CASE (i)

    CASE(1)
    v              = ( PTOT - EPRESS + pe ) * kp
    vd             = DPDN * kp * kfm
    vt             = DPDT * kp * kmev
    vy             = DPDY * kp
    RETURN

    CASE(2)
    v              = ( UTOT + UTOT0 - EU + ee ) * ku
    vd             = DUDN * ku * kfm
    vt             = DUDT * ku * kmev
    vy             = DUDY * ku
    RETURN

    CASE(3)
    v              = STOT - ES + se
    vd             = DSDN * kfm
    vt             = DSDT * kmev
    vy             = DSDY
    RETURN

    CASE(4)
    v              = MUN - dmnp
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

    CASE(5)
    v              = MUPROT
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

    CASE(6)
    v              = cmpe
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

    CASE(7)
    v              = XNUT
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

    CASE(8)
    v              = XPROT
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

    CASE(9)
    v              = XH
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

    CASE(10)
    v              = A
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

    CASE(11)
    v              = X * A
    vd             = zero
    vt             = zero
    vy             = zero
    RETURN

  END SELECT i_LS

!-----------------------------------------------------------------------
!  End Lattimer-Swesty eos
!-----------------------------------------------------------------------

  RETURN
END IF

!-----------------------------------------------------------------------
!
!                   \\\\\ COOPERSTEIN-BCK EOS /////
!
!-----------------------------------------------------------------------

IF ( eos_ii == 'B'  .and.  i /= 12 ) then

!-----------------------------------------------------------------------
!  Initialize BCK inputs
!-----------------------------------------------------------------------

  IF ( firstbck ) then
    DO jj = 1,nx
      ueaa   (jj)  = zero
      unaa   (jj)  = zero
      uhataa (jj)  = zero
      thetaaa(jj)  = zero
      zaaa   (jj)  = zero
      xaaa   (jj)  = zero
      dtranaa(jj)  = zero
    END DO
  firstbck         = .false.
  END IF

!-----------------------------------------------------------------------
!  Initialize BCK inputs
!-----------------------------------------------------------------------

  ue               = ueaa(j)
  un               = unaa(j)
  uhat             = uhataa(j)
  theta            = thetaaa(j)
  zabck            = zaaa(j)
  xabck            = xaaa(j)
  dtran            = dtranaa(j)
  jshel            = j
  dbck             = brydns
  tbck             = tmev
  yebck            = yej

!-----------------------------------------------------------------------
!  Call the BCK equation of state
!-----------------------------------------------------------------------

  CALL eos( ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Store equation of state parameters for zone j
!-----------------------------------------------------------------------

  ueaa   (j)       = ue
  unaa   (j)       = un
  uhataa (j)       = uhat
  thetaaa(j)       = theta
  zaaa   (j)       = zabck
  xaaa   (j)       = xabck
  dtranaa(j)       = dtran

!-----------------------------------------------------------------------
!  Convert thermodynamic quantities from MeV and fm to cgs units, where
!   appropriate.
!-----------------------------------------------------------------------

  i_BCK: SELECT CASE (i)

    CASE(1)
    v              = ptotbck * kp

    CASE(2)
    v              = etot * ku

    CASE(3)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = stotbck
    RETURN

    CASE(4)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = un
    RETURN

    CASE(5)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = un - uhat
    RETURN

    CASE(6)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = ue
    RETURN

    CASE(7)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = xnbck
    RETURN

    CASE(8)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = xpbck
    RETURN

    CASE(9)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = xhbck
    RETURN

    CASE(10)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = ahbck
    RETURN

    CASE(11)
    vd             = zero
    vt             = zero
    vy             = zero
    v              = zabck * ahbck
    RETURN

  END SELECT i_BCK

!-----------------------------------------------------------------------
!  Set up variables for computing derivatives by finite differences.
!-----------------------------------------------------------------------

  DO ii = 1,6
    rhod(ii)       = rhoj
    td(ii)         = tj
    yed(ii)        = yej
  END DO

  rhod(1)          = ( one - delta ) * rhoj
  rhod(2)          = ( one + delta ) * rhoj
  td(3)            = ( one - delta ) * tj
  td(4)            = ( one + delta ) * tj
  yed(5)           = ( one - delta ) * yej
  yed(6)           = ( one + delta ) * yej

!-----------------------------------------------------------------------
!  Compute derivatives of the pressure.
!-----------------------------------------------------------------------
            
  i_deriv: SELECT CASE (i)

    CASE(1)

    DO ii = 1,6
      dbck         = rhod(ii) * kfm
      tbck         = td(ii) * kmev
      yebck        = yed(ii)
      CALL eos( ij_ray, ik_ray )
      pderiv(ii)   = ptotbck * kp
    END DO

    vd             = ( pderiv(2) - pderiv(1) )/( 2.d+00 * delta * rhoj )
    vt             = ( pderiv(4) - pderiv(3) )/( 2.d+00 * delta * tj   )
    vy             = ( pderiv(6) - pderiv(5) )/( 2.d+00 * delta * yej  )
    RETURN

!-----------------------------------------------------------------------
!  Compute derivatives of the energy.
!-----------------------------------------------------------------------

    CASE(2)

    DO ii = 1,6
      dbck         = rhod(ii) * kfm
      tbck         = td(ii) * kmev
      yebck        = yed(ii)
      CALL eos( ij_ray, ik_ray )
      uderiv(ii)   = etot * ku
    END DO

    vd             = ( uderiv(2) - uderiv(1) )/( 2.d+00 * delta * rhoj )
    vt             = ( uderiv(4) - uderiv(3) )/( 2.d+00 * delta * tj   )
    vy             = ( uderiv(6) - uderiv(5) )/( 2.d+00 * delta * yej  )
    RETURN

  END SELECT i_deriv

!-----------------------------------------------------------------------
!
!              \\\\\ END LS OR COOPERSTEIN-BCK EOS /////
!
!-----------------------------------------------------------------------

END IF ! brydns > eosrho, tmev > eost, nse(j,ij_ray,ik_ray) = 1, eos_i = 'L'

IF ( i == 12 ) THEN

!-----------------------------------------------------------------------
!  Initialize BCK inputs.
!-----------------------------------------------------------------------

  IF ( firstbck ) then
    DO jj = 1,nx
      ueaa   (jj)  = zero
      unaa   (jj)  = zero
      uhataa (jj)  = zero
      thetaaa(jj)  = zero
      zaaa   (jj)  = zero
      xaaa   (jj)  = zero
      dtranaa(jj)  = zero
    END DO
  firstbck         = .false.
  END IF

!-----------------------------------------------------------------------
!  Initialize BCK inputs.
!-----------------------------------------------------------------------

  ue               = ueaa(j)
  un               = unaa(j)
  uhat             = uhataa(j)
  theta            = thetaaa(j)
  zabck            = zaaa(j)
  xabck            = xaaa(j)
  dtran            = dtranaa(j)
  jshel            = j
  dbck             = brydns
  tbck             = tmev
  yebck            = yej

!-----------------------------------------------------------------------
!  Call the BCK equation of state.
!-----------------------------------------------------------------------

  CALL eos( ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Store equation of state parameters for zone j.
!-----------------------------------------------------------------------

  ueaa   (j)       = ue
  unaa   (j)       = un
  uhataa (j)       = uhat
  thetaaa(j)       = theta
  zaaa   (j)       = zabck
  xaaa   (j)       = xabck
  dtranaa(j)       = dtran

  RETURN
END IF ! i = 12

END SUBROUTINE eqstta_x