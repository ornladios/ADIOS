SUBROUTINE sctikrnl( j, ij_ray, ik_ray, k, n, rho, t, ye, coh, cohd, coht, &
& cohy )
!-----------------------------------------------------------------------
!
!    File:         sctikrnl
!    Module:       sctikrnl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To interpolate isoenergetic neutrino-nucleus (nucleon) scattering
!       kernals in a local table of nearest entries created for each zone.
!      The table consists of the eight nearest-neighbor entries in rho,
!       t, and ye of the state point.
!      Derivatives of the scattering rates are obtained from the expressions
!       for these rates by direct differentiation of these expressions.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  j                      : radial zone number
!  ij_ray                 : j-index of a radial ray
!  ik_ray                 : k-index of a radial ray
!  n                      : neutrino flavor index
!  rho                    : matter density (g/cm**3)
!  t                      : matter temperature (K)
!  ye                     : matter electron fraction
!  in                     : 0, neutrino-neutron isoenergetic scattering omitted
!                           1, neutrino-neutron isoenergetic scattering included
!  ip                     : 0, neutrino-proton isoenergetic scattering omitted
!                           1, neutrino-proton isoenergetic scattering included
!  ihe                    : 0, neutrino-helium isoenergetic scattering omitted
!                           1, neutrino-helium isoenergetic scattering included
!  iheavy                 : 0, neutrino-heavy nucleus isoenergetic scattering omitted
!                           1, neutrino-heavy nucleus isoenergetic scattering included
!  iscat                  : 0, all neutrino scattering processes omitted
!                           1, neutrino scattering processes not necessarily omitted
!
!    Output arguments:
!
!  coh                    : isoenergetic scattering function
!  cohd                   : d(coh)/d(rho)
!  coht                   : d(coh)/d(t)
!  cohy                   : d(coh)/d(ye)
!
!    Input arguments (common):
!
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5
!                           and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)  : index for dgrid, tgrid, and ygrid for zone j
!  idrsi(j,ij_ray,ik_ray) : rho grid index for zone j
!  itrsi(j,ij_ray,ik_ray) : t grid index for zone j
!  iyrsi(j,ij_ray,ik_ray) : ye grid index for zone j
!   cohsct                : neutrino-nucleus (nucleon) isoenergetic scattering kernal array
!
!    Include files:
!  numerical_module, numerical_module
!  eos_snc_x_module, scat_i.cmn!
!
!----------------------------------------------------------------------- 

USE kind_module
USE numerical_module, ONLY : zero, half, one

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty
USE scat_i_module, ONLY : idrsi, itrsi, iyrsi, cohsct, cohbsct

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: k             ! incoming neutrino energy zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: n             ! neutrino flavor index

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)                 :: coh     ! isoenergetic scattering function
REAL(KIND=double), INTENT(out)                 :: cohd    ! d(coh)/d(rho)
REAL(KIND=double), INTENT(out)                 :: coht    ! d(coh)/d(t)
REAL(KIND=double), INTENT(out)                 :: cohy    ! d(coh)/d(ye)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: id            ! EOS density index
INTEGER                          :: it            ! EOS temperature index
INTEGER                          :: iy            ! EOS electron fraction index

INTEGER, PARAMETER               :: ida = 1       ! lower cube table density index
INTEGER, PARAMETER               :: idap1 = 2     ! upper cube table density index
INTEGER, PARAMETER               :: ita = 1       ! lower cube table temperature index
INTEGER, PARAMETER               :: itap1 = 2     ! upper cube table temperature index
INTEGER, PARAMETER               :: iya = 1       ! lower cube table electron fraction index
INTEGER, PARAMETER               :: iyap1 = 2     ! upper cube table electron fraction index

REAL(KIND=double)                :: log_e         ! log10(e)
REAL(KIND=double)                :: ln_10         ! ln(10)

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

REAL(KIND=double)                :: ch111         ! scalar table entry for interpolation
REAL(KIND=double)                :: ch211         ! scalar table entry for interpolation
REAL(KIND=double)                :: ch121         ! scalar table entry for interpolation
REAL(KIND=double)                :: ch112         ! scalar table entry for interpolation
REAL(KIND=double)                :: ch221         ! scalar table entry for interpolation
REAL(KIND=double)                :: ch212         ! scalar table entry for interpolation
REAL(KIND=double)                :: ch122         ! scalar table entry for interpolation
REAL(KIND=double)                :: ch222         ! scalar table entry for interpolation

REAL(KIND=double)                :: cohl          ! interpolated log of coh
REAL(KIND=double)                :: cohdl         ! interpolated log of d(coh)/d(rho)
REAL(KIND=double)                :: cohtl         ! interpolated log of d(coh)/d(t)
REAL(KIND=double)                :: cohyl         ! interpolated log of d(coh)/d(ye)

REAL(KIND=double)                :: f10           ! 10** function

EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF (first) THEN
  first            = .false.

  log_e            = DLOG10( DEXP( one ) )
  ln_10            = DLOG( 1.d+01 )

END IF

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho )
id                 = idrsi(j,ij_ray,ik_ray)
fdp                = DMAX1( DMIN1( fd - DBLE( id ), one ), zero )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho )

ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   )
it                 = itrsi(j,ij_ray,ik_ray)
ftp                = DMAX1( DMIN1( ft - DBLE( it ), one ), zero )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   ) ! log_e = log(e)

fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )
iy                 = iyrsi(j,ij_ray,ik_ray)
fyp                = DMAX1( DMIN1( fy - DBLE( iy ), one ), zero )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

!-----------------------------------------------------------------------
!  Store coherent scattering functions in temporary scalar variables
!   for interpolation.
!-----------------------------------------------------------------------

IF ( n == 1  .or.  n == 3 ) THEN
  ch111            = cohsct(j,k,ida  ,ita  ,iya  ,ij_ray,ik_ray)
  ch211            = cohsct(j,k,idap1,ita  ,iya  ,ij_ray,ik_ray)
  ch121            = cohsct(j,k,ida  ,itap1,iya  ,ij_ray,ik_ray)
  ch112            = cohsct(j,k,ida  ,ita  ,iyap1,ij_ray,ik_ray)
  ch221            = cohsct(j,k,idap1,itap1,iya  ,ij_ray,ik_ray)
  ch212            = cohsct(j,k,idap1,ita  ,iyap1,ij_ray,ik_ray)
  ch122            = cohsct(j,k,ida  ,itap1,iyap1,ij_ray,ik_ray)
  ch222            = cohsct(j,k,idap1,itap1,iyap1,ij_ray,ik_ray)
ELSE
  ch111            = cohbsct(j,k,ida  ,ita  ,iya  ,ij_ray,ik_ray)
  ch211            = cohbsct(j,k,idap1,ita  ,iya  ,ij_ray,ik_ray)
  ch121            = cohbsct(j,k,ida  ,itap1,iya  ,ij_ray,ik_ray)
  ch112            = cohbsct(j,k,ida  ,ita  ,iyap1,ij_ray,ik_ray)
  ch221            = cohbsct(j,k,idap1,itap1,iya  ,ij_ray,ik_ray)
  ch212            = cohbsct(j,k,idap1,ita  ,iyap1,ij_ray,ik_ray)
  ch122            = cohbsct(j,k,ida  ,itap1,iyap1,ij_ray,ik_ray)
  ch222            = cohbsct(j,k,idap1,itap1,iyap1,ij_ray,ik_ray)
END IF ! n == 1  .or.  n == 3

!-----------------------------------------------------------------------
!  Interpolate coherent scattering kernal
!-----------------------------------------------------------------------

cohl               = fym * ( fdm * ( ftm * ch111 + ftp * ch121 )   &
&                  +         fdp * ( ftm * ch211 + ftp * ch221 ) ) &
&                  + fyp * ( fdm * ( ftm * ch112 + ftp * ch122 )   &
&                  +         fdp * ( ftm * ch212 + ftp * ch222 ) )

coh                = f10(cohl)

!-----------------------------------------------------------------------
!  Interpolate d(coherent scattering kernal)/d(density)
!-----------------------------------------------------------------------

cohdl              = fdd * ( fym * ( ftm * ( -ch111 + ch211 )   &
&                  +                 ftp * ( -ch121 + ch221 ) ) &
&                  +         fyp * ( ftm * ( -ch112 + ch212 )   &
&                  +                 ftp * ( -ch122 + ch222 ) ) )

cohd               = ln_10 * cohdl * coh

!-----------------------------------------------------------------------
!  Interpolate d(coherent scattering kernal)/d(temperature)
!-----------------------------------------------------------------------

cohtl              = ftt * ( fym * ( fdm * ( -ch111 + ch121 )   &
&                  +                 fdp * ( -ch211 + ch221 ) ) &
&                  +         fyp * ( fdm * ( -ch112 + ch122 )   &
&                  +                 fdp * ( -ch212 + ch222 ) ) )

coht               = ln_10 * cohtl * coh

!-----------------------------------------------------------------------
!  Interpolate d(coherent scattering kernal/d(ye)
!-----------------------------------------------------------------------

cohyl              = fyy * ( fdm * ( ftm * ( -ch111 + ch112 )   &
&                  +                 ftp * ( -ch121 + ch122 ) ) &
&                  +         fdp * ( ftm * ( -ch211 + ch212 )   &
&                  +                 ftp * ( -ch221 + ch222 ) ) )

cohy               = ln_10 * cohyl * coh

RETURN
END SUBROUTINE sctikrnl
