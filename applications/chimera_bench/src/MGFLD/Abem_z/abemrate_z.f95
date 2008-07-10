SUBROUTINE abemrate_z( kmin, kmax, ki_ray, kj_ray, rho_in, t_in, ye_in, &
& nz )
!-----------------------------------------------------------------------
!
!    File:         abemrate_z
!    Module:       abemrate_z
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/11/07
!
!    Purpose:
!      To interpolate, for all energy zones ie, all azimuthal zones k,
!       and all neutrino types ie, the logs of the neutrino and
!       antineutrino absorption and emission inverse mean free paths
!       in a local table of nearest entries created for each zone.
!
!      The table consists of the eight nearest-neighbor entries of
!       the log(rho), log(t), and ye grid.
!      Derivatives of the absorption and emission inverse mean free
!       paths are obtained from the interpolation expressions for
!       these rates by direct differentiation of these expressions
!
!      The absorption and emission rates are included in the
!       multi-group diffusion equations, which have the form
!
!
!                      1                  d psi0(ie)
!          psi1(ie)  = - - lamt(ie) F(ie) ----------
!                      3                      dr
!                                 2
!          1  d psi0(ie)    1  d r psi1(ie)
!          -  ---------- +  -  ------------  = RHS0(ie)
!          !      dr        2       dr
!                          r
!        Absorption and emission contribute to the terms RHS0(ie)
!         and lamt(ie) as
!
!
!          RHS0(ie) = emis(ie) [1 - psi0(ie)] - absor(ie) psi0(ie)
!
!                     psi0(ie)
!          lamt(ie) = --------
!                     emis(ie)
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  kmin           : inner azimuthal zone number
!  kmax           : outer azimuthal zone number
!  ki_ray         : x (radial) index of a specific z (azimuthal) ray
!  kj_ray         : y (angular) index of a specific z (azimuthal) ray
!  rho_in         : density (g cm^{-3})
!  t_in           : temperature (K)
!  ye_in          : electron fraction
!  nz             : y-array extent
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  nnugp(n)       : number of energy zones for neutrinos of type n
!  dgrid(idty(k)) : number of table entries per decade in rho for zone k
!  tgrid(idty(k)) : number of table entries per decade in t for zone k
!  ygrid(idty(k)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone k
!  idty(k)        : index for dgrid, tgrid, and ygrid for zone k
!  em             : zero moment of the neutrino emission rate
!  ab             : zero moment of the neutrino absorption rate
!  rho(k)         : density of zone k
!  t(k)           : temperature of zone k
!  ye(k)          : electron fraction of zone k
!  idrae(k,kj_ray,ki_ray) : rho grid index for zone k
!  itrae(k,kj_ray,ki_ray) : t grid index for zone k
!  iyrae(k,kj_ray,ki_ray) : ye grid index for zone k
!
!    Output arguments (common):
!
!  emis(k,ie,n)    : emission inverse mean free path (1/cm) for radial zone k, energy zone ie,
!                    neutrino type n
!  emisd(k,ie,n)   : d(emis(k,ie,n))/d(rho(k))
!  emist(k,ie,n)   : d(emis(k,ie,n))/d(t  (k))
!  emisy(k,ie,n)   : d(emis(k,ie,n))/d(ye (k))
!  absor(k,ie,n)   : absorption inverse mean free path (1/cm) for radial zone k, energy zone ie,
!                    neutrino type n
!  absord(k,ie,n)  : d(absor(k,ie,n))/d(rho(k))
!  absort(k,ie,n)  : d(absor(k,ie,n))/d(t  (k))
!  absory(k,ie,n)  : d(absor(k,ie,n))/d(ye (k))
!
!    Modules used:
!  kind_module, array_module, numerical_module, physcnst_module,
!  abem_module, eos_snc_x_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, ik_ray_dim
USE numerical_module, ONLY : zero, one
USE physcnst_module

USE abem_z_module, ONLY : em, ab, emis, emisd, emist, emisy, absor, &
& absord, absort, absory, idrae, itrae, iyrae
USE eos_snc_z_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : iaefnp, iaencnu, iaence

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: kmin            ! minimum azimuthal zone index
INTEGER, INTENT(in)              :: kmax            ! maximum azimuthal zone index
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: nz              ! azimuthal array extent

REAL(KIND=double), INTENT(in), DIMENSION(nz)    :: rho_in ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz)    :: t_in   ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nz)    :: ye_in  ! electron fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: k             ! azimuthal zone index
INTEGER                          :: ie            ! incoming neutrino energy zone index
INTEGER                          :: n             ! neutrino flavor index

INTEGER                          :: id            ! EOS density index
INTEGER                          :: it            ! EOS temperature index
INTEGER                          :: iy            ! EOS electron fraction index

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

INTEGER, PARAMETER               :: ida = 1       ! lower cube table density index
INTEGER, PARAMETER               :: idap1 = 2     ! upper cube table density index
INTEGER, PARAMETER               :: ita = 1       ! lower cube table temperature index
INTEGER, PARAMETER               :: itap1 = 2     ! upper cube table temperature index
INTEGER, PARAMETER               :: iya = 1       ! lower cube table electron fraction index
INTEGER, PARAMETER               :: iyap1 = 2     ! upper cube table electron fraction index

REAL(KIND=double)                :: log_e         ! log10(e)
REAL(KIND=double)                :: ln_10         ! ln(10)

REAL(KIND=double)                :: rhomin        ! minimum rho of EOS box
REAL(KIND=double)                :: rhomax        ! maximum rho of EOS box
REAL(KIND=double)                :: tmin          ! minimum t of EOS box
REAL(KIND=double)                :: tmax          ! maximum t of EOS box
REAL(KIND=double)                :: yemin         ! minimum ye of EOS box
REAL(KIND=double)                :: yemax         ! maximum ye of EOS box

REAL(KIND=double)                :: rho           ! density (g/cm**3)
REAL(KIND=double)                :: t             ! temperature (K)
REAL(KIND=double)                :: ye            ! electron fraction

REAL(KIND=double)                :: fd            ! position of rho in grid
REAL(KIND=double)                :: fdp           ! position of rho in grid wrt lower cube index
REAL(KIND=double)                :: fdm           ! position of rho in grid wrt upper cube index
REAL(KIND=double)                :: fdd           ! d(fd)/d(rho)
REAL(KIND=double)                :: ft            ! position of t in grid
REAL(KIND=double)                :: ftp           ! position of t in grid wrt lower cube index
REAL(KIND=double)                :: ftm           ! position of t in grid wrt upper cube index
REAL(KIND=double)                :: ftt           ! d(ft)/d(t)
REAL(KIND=double)                :: fy            ! position of ye in grid
REAL(KIND=double)                :: fyp           ! position of ye in grid wrt lower cube index
REAL(KIND=double)                :: fym           ! position of ye in grid wrt upper cube index
REAL(KIND=double)                :: fyy           ! d(fy)/d(ye)

REAL(KIND=double)                :: em111         ! scalar table entry for interpolation
REAL(KIND=double)                :: em211         ! scalar table entry for interpolation
REAL(KIND=double)                :: em121         ! scalar table entry for interpolation
REAL(KIND=double)                :: em112         ! scalar table entry for interpolation
REAL(KIND=double)                :: em221         ! scalar table entry for interpolation
REAL(KIND=double)                :: em212         ! scalar table entry for interpolation
REAL(KIND=double)                :: em122         ! scalar table entry for interpolation
REAL(KIND=double)                :: em222         ! scalar table entry for interpolation

REAL(KIND=double)                :: ab111         ! scalar table entry for interpolation
REAL(KIND=double)                :: ab211         ! scalar table entry for interpolation
REAL(KIND=double)                :: ab121         ! scalar table entry for interpolation
REAL(KIND=double)                :: ab112         ! scalar table entry for interpolation
REAL(KIND=double)                :: ab221         ! scalar table entry for interpolation
REAL(KIND=double)                :: ab212         ! scalar table entry for interpolation
REAL(KIND=double)                :: ab122         ! scalar table entry for interpolation
REAL(KIND=double)                :: ab222         ! scalar table entry for interpolation

REAL(KIND=double)                :: emisl         ! log(emis)
REAL(KIND=double)                :: emisdl        ! log(emisd)
REAL(KIND=double)                :: emistl        ! log(emist)
REAL(KIND=double)                :: emisyl        ! log(emisy)

REAL(KIND=double)                :: absorl        ! log(absor)
REAL(KIND=double)                :: absordl       ! log(absord)
REAL(KIND=double)                :: absortl       ! log(absort)
REAL(KIND=double)                :: absoryl       ! log(absory)

REAL(KIND=double)                :: f10           ! 10** function
EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Set rates to zero and return if
!
! all nnugp(n) = 0, or ( iaefnp = 0 and iaencnu = 0 and iaence = 0 )
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  ln_10            = DLOG(10.d0)
  log_e            = DLOG10( dexp(1.d0) )
END IF ! first


IF ( nnugpmx == 0  .or.  ( iaefnp == 0  .and.  iaencnu == 0  .and.  iaence == 0 ) ) THEN

  emis  (kmin:kmax,:,:) = zero
  emisd (kmin:kmax,:,:) = zero
  emist (kmin:kmax,:,:) = zero
  emisy (kmin:kmax,:,:) = zero
  absor (kmin:kmax,:,:) = zero
  absord(kmin:kmax,:,:) = zero
  absort(kmin:kmax,:,:) = zero
  absory(kmin:kmax,:,:) = zero

  RETURN
END IF

i_ray              = ik_ray_dim * ( ki_ray - 1 ) + kj_ray

!-----------------------------------------------------------------------
!  Loop over k
!-----------------------------------------------------------------------

DO k = kmin,kmax

!-----------------------------------------------------------------------
!  Ensure that rho, t, and ye lie within "unit cube"
!-----------------------------------------------------------------------

  id               = idrae(k,kj_ray,ki_ray)
  it               = itrae(k,kj_ray,ki_ray)
  iy               = iyrae(k,kj_ray,ki_ray)

  rhomin           = 10.d+00**( DBLE(id  )/dgrid(idty(k,kj_ray,ki_ray)) )
  rhomax           = 10.d+00**( DBLE(id+1)/dgrid(idty(k,kj_ray,ki_ray)) )
  tmin             = 10.d+00**( DBLE(it  )/tgrid(idty(k,kj_ray,ki_ray)) )
  tmax             = 10.d+00**( DBLE(it+1)/tgrid(idty(k,kj_ray,ki_ray)) )
  yemax            = one     -( DBLE(iy  )/ygrid(idty(k,kj_ray,ki_ray)) )
  yemin            = one     -( DBLE(iy+1)/ygrid(idty(k,kj_ray,ki_ray)) )

  rho              = DMAX1( dmin1( rho_in(k), rhomax ), rhomin )
  t                = DMAX1( dmin1( t_in(k)  , tmax   ), tmin   )
  ye               = DMAX1( dmin1( ye_in(k) , yemax  ), yemin  )

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

  fd               = dgrid(idty(k,kj_ray,ki_ray)) * LOG10( rho  )
  id               = idrae(k,kj_ray,ki_ray)
  fdp              = fd - DBLE( id )
  fdm              = one - fdp
  fdd              = log_e * dgrid(idty(k,kj_ray,ki_ray))/( rho  )

  ft               = tgrid(idty(k,kj_ray,ki_ray)) * LOG10( t    )
  it               = itrae(k,kj_ray,ki_ray)
  ftp              = ft - DBLE( it )
  ftm              = one - ftp
  ftt              = log_e * tgrid(idty(k,kj_ray,ki_ray))/( t    )

  fy               = ygrid(idty(k,kj_ray,ki_ray)) * ( one - ye )
  iy               = iyrae(k,kj_ray,ki_ray)
  fyp              = fy - DBLE( iy )
  fym              = one - fyp
  fyy              = - ygrid(idty(k,kj_ray,ki_ray))

!-----------------------------------------------------------------------
!  Loop over n and ie
!-----------------------------------------------------------------------

  DO n = 1,2
    DO ie = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store emission and absorption functionsin temporary scalar variables
!   for interpolation
!-----------------------------------------------------------------------

      em111        = em(k,ie,n,ida  ,ita  ,iya  ,i_ray)
      em211        = em(k,ie,n,idap1,ita  ,iya  ,i_ray)
      em121        = em(k,ie,n,ida  ,itap1,iya  ,i_ray)
      em112        = em(k,ie,n,ida  ,ita  ,iyap1,i_ray)
      em221        = em(k,ie,n,idap1,itap1,iya  ,i_ray)
      em212        = em(k,ie,n,idap1,ita  ,iyap1,i_ray)
      em122        = em(k,ie,n,ida  ,itap1,iyap1,i_ray)
      em222        = em(k,ie,n,idap1,itap1,iyap1,i_ray)

      ab111        = ab(k,ie,n,ida  ,ita  ,iya  ,i_ray)
      ab211        = ab(k,ie,n,idap1,ita  ,iya  ,i_ray)
      ab121        = ab(k,ie,n,ida  ,itap1,iya  ,i_ray)
      ab112        = ab(k,ie,n,ida  ,ita  ,iyap1,i_ray)
      ab221        = ab(k,ie,n,idap1,itap1,iya  ,i_ray)
      ab212        = ab(k,ie,n,idap1,ita  ,iyap1,i_ray)
      ab122        = ab(k,ie,n,ida  ,itap1,iyap1,i_ray)
      ab222        = ab(k,ie,n,idap1,itap1,iyap1,i_ray)


!-----------------------------------------------------------------------
!  Interpolate emission inv mfp's
!-----------------------------------------------------------------------

      emisl        =  fym * ( fdm * ( ftm * em111 + ftp * em121 )   &
&                  +          fdp * ( ftm * em211 + ftp * em221 ) ) &
&                  +  fyp * ( fdm * ( ftm * em112 + ftp * em122 )   &
&                  +          fdp * ( ftm * em212 + ftp * em222 ) )

      emis(k,ie,n)  = f10(emisl)

!-----------------------------------------------------------------------
!  Interpolate d(emission rate)/d(density)
!-----------------------------------------------------------------------

      emisdl       = fdd * ( fym * ( ftm * ( -em111 + em211 )   &
&                  +                 ftp * ( -em121 + em221 ) ) &
&                  +         fyp * ( ftm * ( -em112 + em212 )   &
&                  +                 ftp * ( -em122 + em222 ) ) )

      emisd(k,ie,n) = ln_10 * emisdl * emis(k,ie,n)

!-----------------------------------------------------------------------
!  Interpolate d(emission rate)/d(temperature)
!-----------------------------------------------------------------------

      emistl       = ftt * ( fym * ( fdm * ( -em111 + em121 )   &
&                  +                 fdp * ( -em211 + em221 ) ) &
&                  +         fyp * ( fdm * ( -em112 + em122 )   &
&                  +                 fdp * ( -em212 + em222 ) ) )

      emist(k,ie,n) = ln_10 * emistl * emis(k,ie,n)

!-----------------------------------------------------------------------
!  Interpolate d(emission rate)/d(ye)
!-----------------------------------------------------------------------

      emisyl       = fyy * ( fdm * ( ftm * ( -em111 + em112 )   &
&                  +                 ftp * ( -em121 + em122 ) ) &
&                  +         fdp * ( ftm * ( -em211 + em212 )   &
&                  +                 ftp * ( -em221 + em222 ) ) )

      emisy(k,ie,n) = ln_10 * emisyl * emis(k,ie,n)

!-----------------------------------------------------------------------
!  Interpolate absortion inv mfp's
!-----------------------------------------------------------------------

      absorl       =  fym * ( fdm * ( ftm * ab111 + ftp * ab121 )   &
&                  +          fdp * ( ftm * ab211 + ftp * ab221 ) ) &
&                  +  fyp * ( fdm * ( ftm * ab112 + ftp * ab122 )   &
&                  +          fdp * ( ftm * ab212 + ftp * ab222 ) )

      absor(k,ie,n) = f10(absorl)

!-----------------------------------------------------------------------
!  Interpolate d(absortion rate)/d(density)
!-----------------------------------------------------------------------

      absordl      = fdd * ( fym * ( ftm * ( -ab111 + ab211 )   &
&                  +                 ftp * ( -ab121 + ab221 ) ) &
&                  +         fyp * ( ftm * ( -ab112 + ab212 )   &
&                  +                 ftp * ( -ab122 + ab222 ) ) )

      absord(k,ie,n)= ln_10 * absordl * absor(k,ie,n)

!-----------------------------------------------------------------------
!  Interpolate d(absortion rate)/d(temperature)
!-----------------------------------------------------------------------

      absortl      = ftt * ( fym * ( fdm * ( -ab111 + ab121 )   &
&                  +                 fdp * ( -ab211 + ab221 ) ) &
&                  +         fyp * ( fdm * ( -ab112 + ab122 )   &
&                  +                 fdp * ( -ab212 + ab222 ) ) )

      absort(k,ie,n)= ln_10 * absortl * absor(k,ie,n)

!-----------------------------------------------------------------------
!  Interpolate d(absortion rate)/d(ye)
!-----------------------------------------------------------------------

      absoryl      = fyy * ( fdm * ( ftm * ( -ab111 + ab112 )   &
&                  +                 ftp * ( -ab121 + ab122 ) ) &
&                  +         fdp * ( ftm * ( -ab211 + ab212 )   &
&                  +                 ftp * ( -ab221 + ab222 ) ) )

      absory(k,ie,n)= ln_10 * absoryl * absor(k,ie,n)

!-----------------------------------------------------------------------
!  End k,n,ie loop
!-----------------------------------------------------------------------

    END DO
  END DO
END DO

RETURN
END SUBROUTINE abemrate_z
