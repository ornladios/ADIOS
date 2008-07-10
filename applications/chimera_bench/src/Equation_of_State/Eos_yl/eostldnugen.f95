SUBROUTINE eostldnugen( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya, sf, rho )
!-----------------------------------------------------------------------
!
!    File:         eostldnugen
!    Module:       eostldnugen
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/25/02
!
!    Purpose:
!      To call the Lattimer-Swesty or the Cooperstein-BCK equation
!       of state with independent variables rho, t, and ye-yl, and to
!       store the thermodynamic quantities in arrays for
!       interpolation. E-neutrino-antineutrino lepton number
!       (leptons/fm**3) and neutrino energy (MeV/fm**3) are stored
!       in i = 10 and i = 11, respectively. If rho > rho_couple, the eos
!       is iterated on ye until the yl computed equals the yl inputed.
!
!    Variables that must be passed through common:
!  nse(j)         : nuclear statistical equilibrium flag for radial zone j.
!  dgrid(idty(j)) : number of table entries per decade in rho for zone j.
!  tgrid(idty(j)) : number of table entries per decade in t for zone j.
!  ygrid(idty(j)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j.
!  idty(j)        : index for dgrid, tgrid, and ygrid for zone j.
!  rhoes(i)       : density boundary between idty=i and idty=i+1
!  tblt(i,ida,ita,iya)
!                 : temporary equation of state table array.
!  gamhv          : BCK gamma.
!  wnm            : BCK nucl matter energy; - 16.0 MeV.
!  ws             : BCK bulk surface coefficient.
!  xk0            : BCK K_0(x=.5) symmettric matter.
!  xkzafac        : BCK drop in K_0 with x :
!                  k(x) = k(1/2)(1-xkzafac*(1-2x)**2).
!
!    Subprograms called:
!        loadmx, inveos, eos
!
!    Input arguments:
!  j              : radial zone index.
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!  idd            : density grid index
!  itt            : temperature grid index
!  iyy            : electron fraction grid index
!  ida            : equation of state table density index
!  ita            : equation of state table temperature index
!  iya            : equation of state table electron fraction index
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  el_eos_module, eos_m4c_module, eos_bck_module, eos_drv_module, eos_ls_module,
!  eos_snc_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module
USE numerical_module, ONLY: zero, half, one
USE physcnst_module, ONLY: pi, kmev, hbarc, cm3fm3, rmu, dmnp, ergmev

USE edit_module, ONLY: nprint, nlog
USE el_eos_module, ONLY:  EPRESS, EU, ES, MUSUBE
USE eos_m4c_module, ONLY: inpvar, iflag, forflg, brydns, eosflg, xprev, &
& PTOT, UTOT, STOT, XNUT, XPROT, XH, MUN, MUPROT 
USE eos_bck_module, ONLY: jshel, dbck, tbck, yebck, dtran, ue, xnbck, xpbck, &
& xhbck, ptotbck, uhat, etot, stotbck
USE eos_drv_module, ONLY: tblt, rho_couple
USE eos_ls_module, ONLY: inpvars
USE eos_snc_x_module, ONLY: eos_i, idty, dgrid, tgrid, ygrid, eosrho, nse

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: idd           ! density index
INTEGER, INTENT(in)              :: itt           ! temperature index
INTEGER, INTENT(in)              :: iyy           ! lepton fraction index
INTEGER, INTENT(in)              :: ida           ! density index in EOS cube array
INTEGER, INTENT(in)              :: ita           ! temperature index in EOS cube array
INTEGER, INTENT(in)              :: iya           ! lepton fraction index in EOS cube array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i             ! do index

CHARACTER (LEN=1)                :: eos_ls

LOGICAL, DIMENSION(5000)         :: firstls = (/ (.true., i = 1,5000) /)
LOGICAL                          :: first = .true.
LOGICAL                          :: yl_nonconv

INTEGER                          :: it            ! iteration index
INTEGER, PARAMETER               :: itmax = 40    ! maximum number of bisection iterations
INTEGER                          :: sf            ! equation of state flag

REAL(KIND=double), PARAMETER     :: eost = 0.05d0 ! lower boundary of LS EOS table

REAL(KIND=double)                :: pi2           ! pi**2
REAL(KIND=double)                :: pi4           ! pi**4
REAL(KIND=double)                :: kfm           ! ( # nucleons/gram )( cm3/fm3 )
REAL(KIND=double)                :: kp            ! ( erg/cm3 ) / ( mev/fm3 )
REAL(KIND=double)                :: ku            ! ( # nucleons/gram )( erg/mev )
REAL(KIND=double)                :: rhod          ! density at the cube corners
REAL(KIND=double)                :: td            ! temperature at the cube corners
REAL(KIND=double)                :: yld           ! lepton fraction at the cube corners
REAL(KIND=double)                :: ye_yl         ! electron or lepton at the cube corners
REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double), DIMENSION(1000,4) :: guess       ! initial values for the EOS
REAL(KIND=double)                :: pprev         ! previous value of the proton number
REAL(KIND=double)                :: rho           ! density [gm cm^{-3}]

REAL(KIND=double)                :: ye            ! electron fraction
REAL(KIND=double)                :: yemin         ! minimum electron fraction (guess)
REAL(KIND=double)                :: yemax         ! maximum electron fraction (guess)

REAL(KIND=double)                :: pe            ! electron pressure (dynes cm^{-2})
REAL(KIND=double)                :: ee            ! electron energy (ergs cm^{-3})
REAL(KIND=double)                :: se            ! electron entropy
REAL(KIND=double)                :: cmpe          ! electron chemical potential (MeV)
REAL(KIND=double)                :: yeplus        ! positron fraction
REAL(KIND=double)                :: rel           ! relativity parameter

REAL(KIND=double)                :: etae          ! electron chemical potential / kT
REAL(KIND=double)                :: etahat        ! neutron - proton chemical potential / kT
REAL(KIND=double)                :: etanu         ! electron neutrino chemical potential / kT
REAL(KIND=double)                :: nnuex         ! electron neutrino - electron antineutrino number
REAL(KIND=double)                :: yneu          ! nnuex per baryon
REAL(KIND=double)                :: yltest        ! iterated value of yl
REAL(KIND=double), PARAMETER     :: tol = 1.d-4   ! iteration tolerence

REAL(KIND=double), PARAMETER     :: UTOT0 = 8.9d0 ! change in the zero of energy (MeV)
REAL(KIND=double)                :: eneue         ! electron neutrino + antineutrino energy (MeV/fm3)
REAL(KIND=double)                :: eneumt        ! mu and tau neutrino pair energy (MeV/fm3)
REAL(KIND=double)                :: eneu          ! total neutrino energy (MeV/fm3)
REAL(KIND=double)                :: eneub         ! total neutrino energy per baryon (MeV)
REAL(KIND=double)                :: pneu          ! total neutrino pressure (MeV/fm3)
REAL(KIND=double)                :: uneu          ! electron neutrino chemical potential (MeV)
REAL(KIND=double)                :: sneu          ! total neutrino entropy (/fm3)
REAL(KIND=double)                :: sneub         ! dimensionless total neutrino entropy per baryon
REAL(KIND=double)                :: told          ! 

REAL(KIND=double)                :: ps            ! pressure
REAL(KIND=double)                :: us            ! energy
REAL(KIND=double)                :: ss            ! entropy

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

   51 FORMAT (' yl will not converge in subroutine eostldnugen')
   52 FORMAT (' ye=',1pe12.5,' yemin=',1pe12.5,' yemax=',1pe12.5,' yltest=',1pe12.5,' ylj=',1pe12.5)
  101 FORMAT (' Equation of state identifier does not match either L or B')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Examine equation of state identifier
!-----------------------------------------------------------------------

IF ( eos_i /= 'L'  .and.  eos_i /= 'B'  .and.  eos_i /= 'S' ) THEN
  WRITE (nprint,101)
  WRITE (nlog,101)
  STOP
END IF

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN

  pi2              = pi * pi       ! pi**2
  pi4              = pi2 * pi2     ! pi**4
  kfm              = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
  kp               = ergmev/cm3fm3 ! ( erg/cm3 ) / ( mev/fm3 )
  ku               = ergmev/rmu    ! ( # nucleons/gram )( erg/mev )
  iflag            = 1
  forflg           = 0
  CALL loadmx()
  first            = .false.
END IF

!-----------------------------------------------------------------------
!  Compute independent variables at grid point
!-----------------------------------------------------------------------

rhod               = 10.d0**( DBLE(idd)/dgrid(idty(j,ij_ray,ik_ray)) )
td                 = 10.d0**( DBLE(itt)/tgrid(idty(j,ij_ray,ik_ray)) )
yld                = one  - ( DBLE(iyy)/ygrid(idty(j,ij_ray,ik_ray)) )

!-----------------------------------------------------------------------
!  Convert to independent variables of the Lattimer-Swesty and
!   Cooperstein BCK equation of stat
!-----------------------------------------------------------------------

brydns             = rhod * kfm
tmev               = td * kmev
ye_yl              = yld

!-----------------------------------------------------------------------
!  Test the low density, low temperature, high ye corner of
!   the unit cell to determine which equation of state to use.
!   The other corners must use the same equation of state.
!
!     If     brydns       > eosrho         and
!            tmev         > eost           and
!            nse(j,ij_ray,ik_ray) = 1      and
!            eos_i        = 'L'
!     then            use Lattimer-Swesty eos
!     otherwise,      use Cooperstein bck eos
!-----------------------------------------------------------------------

IF ( ida == 1  .and.  ita == 1  .and.  iya == 1 ) THEN
  yl_nonconv       = .false.
  IF ( brydns >= eosrho           .and.  &
&      tmev >= eost               .and.  &
&      nse(j,ij_ray,ik_ray) == 1  .and.  &
&      eos_i /= 'S' )                             THEN
    eos_ls         = 'L'
  ELSE
    eos_ls         = 'B'
  END IF ! brydns > eosrho, tmev > eost, nse(j,ij_ray,ik_ray) = 1, eos_i /= 'S'
END IF ! ida = 1, ita = 1, iya = 1

IF ( eos_ls == 'L' ) THEN

!-----------------------------------------------------------------------
!
!       \\\\\ CALL THE LATTIMER-SWESTY EQUATION OF STATE /////
!
!-----------------------------------------------------------------------

  IF ( rho > rho_couple ) THEN

!-----------------------------------------------------------------------
!  Initialize for bisection to determine ye.
!-----------------------------------------------------------------------

    yemin          = 0.02d0    ! minimum temperature (guess)
    yemax          = 0.5d0     ! maximum temperature (guess)

!-----------------------------------------------------------------------
!  Bisection iteration
!-----------------------------------------------------------------------

  biesction_1: DO it = 1,itmax

      ye           = half * ( yemin + yemax )

!-----------------------------------------------------------------------
!  Load Lattimer-Swesty equation of state parameters
!-----------------------------------------------------------------------

      IF ( firstls(j) ) THEN

        guess(j,1) = ye * brydns
        guess(j,2) = 0.155d0
        guess(j,3) = -15.d0
        guess(j,4) = -10.d0
        firstls(j) = .false.

      END IF ! firstls(j)

      pprev        = dmin1( guess(j,1), brydns * ye )
      inpvar(1)    = tmev
      inpvar(2)    = guess(j,2)
      inpvar(3)    = guess(j,3)
      inpvar(4)    = guess(j,4)

!-----------------------------------------------------------------------
!  Call the Lattimer-Swesty equation of state
!-----------------------------------------------------------------------

      CALL inveos( inpvar, told, ye, brydns, iflag, eosflg, forflg, sf, &
&      xprev, pprev)

!-----------------------------------------------------------------------
!  Store equation of state parameters for zone j
!-----------------------------------------------------------------------

      guess(j,1)   = pprev
      guess(j,2)   = inpvar(2)
      guess(j,3)   = inpvar(3)
      guess(j,4)   = inpvar(4)

!-----------------------------------------------------------------------
!  Compute electron equation of state
!-----------------------------------------------------------------------

      CALL e_p_eos( brydns, tmev, ye, pe, ee, se, cmpe, yeplus, rel )

!-----------------------------------------------------------------------
!  Compute yneu and yltest
!-----------------------------------------------------------------------

      etae         = cmpe/tmev
      etahat       = ( MUN - MUPROT )/tmev
      etanu        = etae - etahat
      nnuex        = 0.5d0 * tmev**3/( pi2 * hbarc**3 ) * ( etanu**3 + pi2 * etanu )/3.0d0
      yneu         = nnuex/brydns
      yltest       = ye + yneu

!-----------------------------------------------------------------------
!  Test for convergence of yltest to yl
!-----------------------------------------------------------------------

      IF ( DABS(yltest - ye_yl) <= tol * DMAX1( DABS(ye_yl), 0.05d+00 ) ) EXIT

!-----------------------------------------------------------------------
!  Bisect the interval if not converged
!-----------------------------------------------------------------------

      IF ( yltest <= ye_yl ) THEN
        yemin      = ye
      ELSE
        yemax      = ye
      END IF ! yltest <= yl

!-----------------------------------------------------------------------
!  Not converged? WRITE diagnostics and go on.
!-----------------------------------------------------------------------

      IF ( it == itmax ) THEN
        yl_nonconv = .true.
      END IF

    END DO biesction_1

    IF ( ida == 4  .and.  ita == 4  .and.  iya == 4  .and.  yl_nonconv ) THEN
      WRITE (nlog,51)
      WRITE (nlog,52) ye, yemin, yemax, yltest, ye_yl
      yl_nonconv   = .false.
    END IF

!-----------------------------------------------------------------------
!  Compute thermodynamic functions for neutrinos.
!
!     eneue  : electron neutrino + antineutrino energy (MeV/fm3).
!     eneumt : mu and tau neutrino pair energy (MeV/fm3).
!     eneu   : total neutrino energy (MeV/fm3).
!     eneub  : total neutrino energy per baryon (MeV).
!     pneu   : total neutrino pressure (MeV/fm3).
!     uneu   : electron neutrino chemical potential (MeV).
!     sneu   : total neutrino entropy (/fm3).
!     sneub  : dimensionless total neutrino entropy per baryon.
!-----------------------------------------------------------------------

    eneue          = half * tmev**4/( pi2 * hbarc**3 ) * 0.25d0           &
&                  * ( etanu**4 + 2.d0 * pi2 * etanu**2 + 7.d0 * pi4/15.0d0 )
    eneumt         = 0.5d0 * tmev**4/( pi2 * hbarc**3 ) * 0.25d0 * ( 2.d0 * 7.d0 * pi4/15.0d0 )
    eneu           = eneue + eneumt
    eneub          = eneu/brydns
    pneu           = eneu/3.0d0
    uneu           = etanu*tmev
    sneu           = ( eneu + pneu - uneu * nnuex )/tmev
    sneub          = sneu/brydns

  ELSE ! rho <= rho_couple

!-----------------------------------------------------------------------
!  Load Lattimer-Swesty equation of state parameters
!-----------------------------------------------------------------------

    IF ( firstls(j) ) THEN

      guess(j,1)   = ye * brydns
      guess(j,2)   = 0.155d0
      guess(j,3)   = -15.d0
      guess(j,4)   = -10.d0
      firstls(j)   = .false.

    END IF ! firstls(j)

    pprev          = dmin1( guess(j,1), brydns * ye )
    inpvar(1)      = tmev
    inpvar(2)      = guess(j,2)
    inpvar(3)      = guess(j,3)
    inpvar(4)      = guess(j,4)
    ye             = ye_yl

!-----------------------------------------------------------------------
!  Call the Lattimer-Swesty equation of state
!-----------------------------------------------------------------------

    CALL inveos( inpvar, told, ye, brydns, iflag, eosflg, forflg, sf, &
&    xprev, pprev)

!-----------------------------------------------------------------------
!  Store equation of state parameters for zone j
!-----------------------------------------------------------------------

    guess(j,1)     = pprev
    guess(j,2)     = inpvar(2)
    guess(j,3)     = inpvar(3)
    guess(j,4)     = inpvar(4)

!-----------------------------------------------------------------------
!  Compute electron equation of state
!-----------------------------------------------------------------------

    CALL e_p_eos( brydns, tmev, ye, pe, ee, se, cmpe, yeplus, rel )

!-----------------------------------------------------------------------
!  Set neutrino thermodynamic quantities to zero
!-----------------------------------------------------------------------

    pneu           = zero
    eneub          = zero
    sneub          = zero

  END IF ! rho > rho_couple

!-----------------------------------------------------------------------
!  Convert thermodynamic quantities from units of mev and fm to cgs
!   units and store table entries.
!-----------------------------------------------------------------------

  tblt(1 ,ida,ita,iya) = ( PTOT - EPRESS + pe + pneu ) * kp
  tblt(2 ,ida,ita,iya) = ( UTOT + UTOT0 - EU + ee + eneub ) * ku
  tblt(3 ,ida,ita,iya) = STOT - ES + se + sneub

!-----------------------------------------------------------------------
!  End Lattimer-Swesty eos
!-----------------------------------------------------------------------

ELSE ! eos_ls /= 'L'


!-----------------------------------------------------------------------
!
!        \\\\\ USE THE COOPERSTEIN BCK EQUATION OF STATE /////
!
!-----------------------------------------------------------------------

  IF ( rho > rho_couple ) THEN

  yemin            = 0.02d0    ! minimum electron fraction (guess)
  yemax            = 0.5d0     ! maximum electron fraction (guess)

!-----------------------------------------------------------------------
!  Initialize for bisection to determine ye.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Bisection iteration
!-----------------------------------------------------------------------

    biesction_2: DO it = 1,itmax

      ye           = half * ( yemin + yemax )

!-----------------------------------------------------------------------
!  Load BCK equation of state parameters
!-----------------------------------------------------------------------

      jshel        = j
      dbck         = brydns
      tbck         = tmev
      yebck        = ye

!-----------------------------------------------------------------------
!  Call the bck equation of state
!-----------------------------------------------------------------------

      CALL eos( ij_ray, ik_ray )

      IF ( xnbck<  zero  .or.  xpbck < zero  .or.  stot < zero ) THEN
        dtran        = dbck
        CALL eos( ij_ray, ik_ray )
      END IF

!-----------------------------------------------------------------------
!  Compute yneu and yltest
!-----------------------------------------------------------------------

      etae         = ue/tmev
      etahat       = ( uhat + dmnp )/tmev
      etanu        = etae - etahat
      nnuex        = half * tmev**3/( pi2 * hbarc**3 ) * ( etanu**3 + pi2 * etanu )/3.0d0
      yneu         = nnuex/dbck
      yltest       = ye + yneu

!-----------------------------------------------------------------------
!  Test for convergence of yltest to ye_yl
!-----------------------------------------------------------------------

      IF ( DABS(yltest - ye_yl) <= tol * DMAX1( DABS(ye_yl), 0.05d+00 ) ) EXIT

!-----------------------------------------------------------------------
!  Bisect the interval if not converged
!-----------------------------------------------------------------------

      IF ( yltest <= ye_yl ) THEN
        yemin      = ye
      ELSE
        yemax      = ye
      END IF

      IF ( it == itmax ) THEN
        yl_nonconv = .true.
      END IF

    END DO biesction_2

    IF ( ida == 4  .and.  ita == 4  .and.  iya == 4  .and.  yl_nonconv ) THEN
      WRITE (nprint,51)
      WRITE (nlog,51)
      WRITE (nprint,52) ye, yemin, yemax, yltest, ye_yl
      WRITE (nlog,52) ye, yemin, yemax, yltest, ye_yl
      yl_nonconv   = .false.
    END IF

!-----------------------------------------------------------------------
!  Compute thermodynamic functions for neutrinos.
!
!     eneue  : electron neutrino + antineutrino energy (MeV/fm3).
!     eneumt : mu and tau neutrino pair energy (MeV/fm3).
!     eneu   : total neutrino energy (MeV/fm3).
!     eneub  : total neutrino energy per baryon (MeV).
!     pneu   : total neutrino pressure (MeV/fm3).
!     uneu   : electron neutrino chemical potential (MeV).
!     sneu   : total neutrino entropy (/fm3).
!     sneub  : dimensionless total neutrino entropy per baryon.
!-----------------------------------------------------------------------

    eneue          = half * tbck**4/( pi2 * hbarc**3 ) * 0.25d0            &
&                  * ( etanu**4 + 2.d0 * pi2 * etanu**2 + 7.d0*pi4/15.0d0 )
    eneumt         = half * tbck**4/( pi2 * hbarc**3 ) * 0.25d0 * ( 2.d0 * 7.d0 * pi4/15.0d0 )
    eneu           = eneue + eneumt
    eneub          = eneu/dbck
    pneu           = eneu/3.0d0
    uneu           = etanu * tbck
    sneu           = ( eneu + pneu - uneu * nnuex )/tbck
    sneub          = sneu/dbck

  ELSE ! rho <= rho_couple

!-----------------------------------------------------------------------
!  Load BCK equation of state parameters
!-----------------------------------------------------------------------

    jshel          = j
    dbck           = brydns
    tbck           = tmev
    yebck          = ye_yl

!-----------------------------------------------------------------------
!  Call the bck equation of state
!-----------------------------------------------------------------------

    CALL eos( ij_ray, ik_ray )

    IF ( xnbck  <  zero  .or.  xpbck < zero  .or.  stot < zero ) THEN
      dtran        = dbck
      CALL eos( ij_ray, ik_ray )
    END IF

!-----------------------------------------------------------------------
!  Set neutrno thermodynamic quantities to zero
!-----------------------------------------------------------------------

    pneu           = zero
    eneub          = zero
    sneub          = zero

  END IF ! rho > rho_couple

!-----------------------------------------------------------------------
!  Convert thermodynamic quantities from units of mev and fm  to cgs
!   units
!-----------------------------------------------------------------------

  ps               = ( ptotbck + pneu ) * kp
  us               = ( etot + eneub ) * ku
  ss               = stotbck + sneub

!-----------------------------------------------------------------------
!  Store table entries
!-----------------------------------------------------------------------

  tblt(1 ,ida,ita,iya) = nnuex
  tblt(2 ,ida,ita,iya) = eneue
  tblt(3 ,ida,ita,iya) = ss

!-----------------------------------------------------------------------
!  End Cooperstein-BCK equation of state
!-----------------------------------------------------------------------

END IF ! eosls == 'L'


RETURN
END SUBROUTINE eostldnugen

