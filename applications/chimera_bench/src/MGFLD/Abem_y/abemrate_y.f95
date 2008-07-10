SUBROUTINE abemrate_y( jmin, jmax, ji_ray, jk_ray, rho_in, t_in, ye_in, &
& ny )
!-----------------------------------------------------------------------
!
!    File:         abemrate_y
!    Module:       abemrate_y
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/11/00
!
!    Purpose:
!      To interpolate, for all energy zones k, all angular zones j,
!       and all neutrino types k, the logs of the neutrino and
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
!                      1               d psi0(k)
!          psi1(k)  = - - lamt(k) F(k) --------
!                      3                  dr
!                              2
!          1  d psi0(k)    1  d r psi1(k)
!          -  --------- +  -  -----------  = RHS0(k)
!          !     dr         2     dr
!                          r
!        Absorption and emission contribute to the terms RHS0(k)
!         and lamt(k) as
!
!
!          RHS0(k) = emis(k) [1 - psi0(k)] - absor(k) psi0(k)
!
!                    psi0(k)
!          lamt(k) = -------
!                    emis(k)
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  jmin           : inner angular zone number
!  jmax           : outer angular zone number
!  ji_ray         : x (radial) index of a specific angular ray
!  jk_ray         : z (azimuthal) index of a specific angular ray
!  rho_in         : density (g cm^{-3})
!  t_in           : temperature (K)
!  ye_in          : electron fraction
!  ny             : y-array extent
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  nnugp(n)       : number of energy zones for neutrinos of type n
!  dgrid(idty(j)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j)) : number of table entries per decade in t for zone j
!  ygrid(idty(j)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j)        : index for dgrid, tgrid, and ygrid for zone j
!  em             : zero moment of the neutrino emission rate
!  ab             : zero moment of the neutrino absorption rate
!  rho(j)         : density of zone j
!  t(j)           : temperature of zone j
!  ye(j)          : electron fraction of zone j
!  idrae(j,ji_ray,jk_ray) : rho grid index for zone j
!  itrae(j,ji_ray,jk_ray) : t grid index for zone j
!  iyrae(j,ji_ray,jk_ray) : ye grid index for zone j
!
!    Output arguments (common):
!
!  emis(j,k,n)    : emission inverse mean free path (1/cm) for angular zone j, energy zone k,
!                    neutrino type n
!  emisd(j,k,n)   : d(emis(j,k,n))/d(rho(j))
!  emist(j,k,n)   : d(emis(j,k,n))/d(t  (j))
!  emisy(j,k,n)   : d(emis(j,k,n))/d(ye (j))
!  absor(j,k,n)   : absorption inverse mean free path (1/cm) for angular zone j, energy zone k,
!                    neutrino type n
!  absord(j,k,n)  : d(absor(j,k,n))/d(rho(j))
!  absort(j,k,n)  : d(absor(j,k,n))/d(t  (j))
!  absory(j,k,n)  : d(absor(j,k,n))/d(ye (j))
!
!    Modules used:
!  kind_module, array_module, numerical_module, physcnst_module,
!  abem_module, eos_snc_x_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : zero, one
USE physcnst_module

USE abem_y_module, ONLY : em, ab, emis, emisd, emist, emisy, absor, &
& absord, absort, absory, idrae, itrae, iyrae
USE eos_snc_y_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : iaefnp, iaencnu, iaence

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jmin            ! minimum angular zone index
INTEGER, INTENT(in)              :: jmax            ! maximum angular zone index
INTEGER, INTENT(in)              :: ji_ray          ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray          ! z (azimuthal) index of a specific angular ray
INTEGER, INTENT(in)              :: ny              ! angular array extent

REAL(KIND=double), INTENT(in), DIMENSION(ny)    :: rho_in ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny)    :: t_in   ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny)    :: ye_in  ! electron fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: j             ! angular zone index
INTEGER                          :: k             ! incoming neutrino energy zone index
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

  emis  (jmin:jmax,:,:) = zero
  emisd (jmin:jmax,:,:) = zero
  emist (jmin:jmax,:,:) = zero
  emisy (jmin:jmax,:,:) = zero
  absor (jmin:jmax,:,:) = zero
  absord(jmin:jmax,:,:) = zero
  absort(jmin:jmax,:,:) = zero
  absory(jmin:jmax,:,:) = zero

  RETURN
END IF

i_ray              = ij_ray_dim * ( jk_ray - 1 ) + ji_ray

!-----------------------------------------------------------------------
!  Loop over j
!-----------------------------------------------------------------------

DO j = jmin,jmax

!-----------------------------------------------------------------------
!  Ensure that rho, t, and ye lie within "unit cube"
!-----------------------------------------------------------------------

  id               = idrae(j,ji_ray,jk_ray)
  it               = itrae(j,ji_ray,jk_ray)
  iy               = iyrae(j,ji_ray,jk_ray)

  rhomin           = 10.d+00**( DBLE(id  )/dgrid(idty(j,ji_ray,jk_ray)) )
  rhomax           = 10.d+00**( DBLE(id+1)/dgrid(idty(j,ji_ray,jk_ray)) )
  tmin             = 10.d+00**( DBLE(it  )/tgrid(idty(j,ji_ray,jk_ray)) )
  tmax             = 10.d+00**( DBLE(it+1)/tgrid(idty(j,ji_ray,jk_ray)) )
  yemax            = one     -( DBLE(iy  )/ygrid(idty(j,ji_ray,jk_ray)) )
  yemin            = one     -( DBLE(iy+1)/ygrid(idty(j,ji_ray,jk_ray)) )

  rho              = DMAX1( dmin1( rho_in(j), rhomax ), rhomin )
  t                = DMAX1( dmin1( t_in(j)  , tmax   ), tmin   )
  ye               = DMAX1( dmin1( ye_in(j) , yemax  ), yemin  )

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

  fd               = dgrid(idty(j,ji_ray,jk_ray)) * LOG10( rho  )
  id               = idrae(j,ji_ray,jk_ray)
  fdp              = fd - DBLE( id )
  fdm              = one - fdp
  fdd              = log_e * dgrid(idty(j,ji_ray,jk_ray))/( rho  )

  ft               = tgrid(idty(j,ji_ray,jk_ray)) * LOG10( t    )
  it               = itrae(j,ji_ray,jk_ray)
  ftp              = ft - DBLE( it )
  ftm              = one - ftp
  ftt              = log_e * tgrid(idty(j,ji_ray,jk_ray))/( t    )

  fy               = ygrid(idty(j,ji_ray,jk_ray)) * ( one - ye )
  iy               = iyrae(j,ji_ray,jk_ray)
  fyp              = fy - DBLE( iy )
  fym              = one - fyp
  fyy              = - ygrid(idty(j,ji_ray,jk_ray))

!-----------------------------------------------------------------------
!  Loop over n and k
!-----------------------------------------------------------------------

  DO n = 1,2
    DO k = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store emission and absorption functionsin temporary scalar variables
!   for interpolation
!-----------------------------------------------------------------------

      em111        = em(j,k,n,ida  ,ita  ,iya  ,i_ray)
      em211        = em(j,k,n,idap1,ita  ,iya  ,i_ray)
      em121        = em(j,k,n,ida  ,itap1,iya  ,i_ray)
      em112        = em(j,k,n,ida  ,ita  ,iyap1,i_ray)
      em221        = em(j,k,n,idap1,itap1,iya  ,i_ray)
      em212        = em(j,k,n,idap1,ita  ,iyap1,i_ray)
      em122        = em(j,k,n,ida  ,itap1,iyap1,i_ray)
      em222        = em(j,k,n,idap1,itap1,iyap1,i_ray)

      ab111        = ab(j,k,n,ida  ,ita  ,iya  ,i_ray)
      ab211        = ab(j,k,n,idap1,ita  ,iya  ,i_ray)
      ab121        = ab(j,k,n,ida  ,itap1,iya  ,i_ray)
      ab112        = ab(j,k,n,ida  ,ita  ,iyap1,i_ray)
      ab221        = ab(j,k,n,idap1,itap1,iya  ,i_ray)
      ab212        = ab(j,k,n,idap1,ita  ,iyap1,i_ray)
      ab122        = ab(j,k,n,ida  ,itap1,iyap1,i_ray)
      ab222        = ab(j,k,n,idap1,itap1,iyap1,i_ray)


!-----------------------------------------------------------------------
!  Interpolate emission inv mfp's
!-----------------------------------------------------------------------

      emisl        =  fym * ( fdm * ( ftm * em111 + ftp * em121 )   &
&                  +          fdp * ( ftm * em211 + ftp * em221 ) ) &
&                  +  fyp * ( fdm * ( ftm * em112 + ftp * em122 )   &
&                  +          fdp * ( ftm * em212 + ftp * em222 ) )

      emis(j,k,n)  = f10(emisl)

!-----------------------------------------------------------------------
!  Interpolate d(emission rate)/d(density)
!-----------------------------------------------------------------------

      emisdl       = fdd * ( fym * ( ftm * ( -em111 + em211 )   &
&                  +                 ftp * ( -em121 + em221 ) ) &
&                  +         fyp * ( ftm * ( -em112 + em212 )   &
&                  +                 ftp * ( -em122 + em222 ) ) )

      emisd(j,k,n) = ln_10 * emisdl * emis(j,k,n)

!-----------------------------------------------------------------------
!  Interpolate d(emission rate)/d(temperature)
!-----------------------------------------------------------------------

      emistl       = ftt * ( fym * ( fdm * ( -em111 + em121 )   &
&                  +                 fdp * ( -em211 + em221 ) ) &
&                  +         fyp * ( fdm * ( -em112 + em122 )   &
&                  +                 fdp * ( -em212 + em222 ) ) )

      emist(j,k,n) = ln_10 * emistl * emis(j,k,n)

!-----------------------------------------------------------------------
!  Interpolate d(emission rate)/d(ye)
!-----------------------------------------------------------------------

      emisyl       = fyy * ( fdm * ( ftm * ( -em111 + em112 )   &
&                  +                 ftp * ( -em121 + em122 ) ) &
&                  +         fdp * ( ftm * ( -em211 + em212 )   &
&                  +                 ftp * ( -em221 + em222 ) ) )

      emisy(j,k,n) = ln_10 * emisyl * emis(j,k,n)

!-----------------------------------------------------------------------
!  Interpolate absortion inv mfp's
!-----------------------------------------------------------------------

      absorl       =  fym * ( fdm * ( ftm * ab111 + ftp * ab121 )   &
&                  +          fdp * ( ftm * ab211 + ftp * ab221 ) ) &
&                  +  fyp * ( fdm * ( ftm * ab112 + ftp * ab122 )   &
&                  +          fdp * ( ftm * ab212 + ftp * ab222 ) )

      absor(j,k,n) = f10(absorl)

!-----------------------------------------------------------------------
!  Interpolate d(absortion rate)/d(density)
!-----------------------------------------------------------------------

      absordl      = fdd * ( fym * ( ftm * ( -ab111 + ab211 )   &
&                  +                 ftp * ( -ab121 + ab221 ) ) &
&                  +         fyp * ( ftm * ( -ab112 + ab212 )   &
&                  +                 ftp * ( -ab122 + ab222 ) ) )

      absord(j,k,n)= ln_10 * absordl * absor(j,k,n)

!-----------------------------------------------------------------------
!  Interpolate d(absortion rate)/d(temperature)
!-----------------------------------------------------------------------

      absortl      = ftt * ( fym * ( fdm * ( -ab111 + ab121 )   &
&                  +                 fdp * ( -ab211 + ab221 ) ) &
&                  +         fyp * ( fdm * ( -ab112 + ab122 )   &
&                  +                 fdp * ( -ab212 + ab222 ) ) )

      absort(j,k,n)= ln_10 * absortl * absor(j,k,n)

!-----------------------------------------------------------------------
!  Interpolate d(absortion rate)/d(ye)
!-----------------------------------------------------------------------

      absoryl      = fyy * ( fdm * ( ftm * ( -ab111 + ab112 )   &
&                  +                 ftp * ( -ab121 + ab122 ) ) &
&                  +         fdp * ( ftm * ( -ab211 + ab212 )   &
&                  +                 ftp * ( -ab221 + ab222 ) ) )

      absory(j,k,n)= ln_10 * absoryl * absor(j,k,n)

!-----------------------------------------------------------------------
!  End j,n,k loop
!-----------------------------------------------------------------------

    END DO
  END DO
END DO

RETURN
END SUBROUTINE abemrate_y
