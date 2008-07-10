SUBROUTINE eosidx_x( j, ij_ray, ik_ray, rho, t, ye, fdm, fdp, fdd, ftm, &
& ftp, ftt, fym, fyp, fyy )
!-----------------------------------------------------------------------
!
!    File:         eosidx_x
!    Module:       eosidx_x
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!        To compute the independent variable grid indices and
!         interpolation coefficients for subroutine eqstt_x.
!
!    Input arguments:
!
!  j                    : radial zone index.
!  ij_ray               : j-index of a radial ray
!  ik_ray               : k-index of a radial ray
!  rho                  : matter density (g/cm**3).
!  t                    : matter temperature (K).
!  ye                   : electron fraction.
!
!    Output arguments:
!
!  fdp                  : distance in log10(rho) from the nearest
!                          lower value of log10(rho) on the rho grid;
!                          used for interpolating equation of state quantities.
!  fdm                  : 1 - fdp; distance in log10(rho) from the nearest
!                          upper value of log10(rho) on the rho grid; used
!                          for interpolating equation of state quantities.
!  fdd                  : derivative of fdp wrt rho; used for computing
!                          derivatives of equation of state quantities wrt rho.
!  ftp                  : distance in log10(t) from the nearest lower
!                          value of log10(t) on the t grid; used for
!                          interpolating equation of state quantities.
!  ftm                  : 1 - ftp; distance in log10(t) from the nearest
!                          upper value of log10(t) on the t grid; used for
!                          interpolating equation of state quantities.
!  ftt                  : derivative of ftp wrt t; used for computing
!                          derivatives of equation of state quantities wrt t.
!  fyp                  : distance in ye from the nearest upper value of
!                          ye on the ye grid; used for interpolating equation
!                          of state quantities.
!  fym                  : 1 - fyp; distance in ye from the nearest lower
!                          value of ye on the ye grid; used for interpolating
!                          equation of state quantities.
!  fyy                  : derivative of ftp wrt ye; used for computing
!                          derivatives of equation of state quantities wrt ye.
!
!    Input arguments (common):
!
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j.
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j.
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5
!                          and ye = 0 for zone j.
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j.
!  idr(j,ij_ray,ik_ray)         : rho grid index for zone j.
!  itr(j,ij_ray,ik_ray)         : t grid index for zone j.
!  iyr(j,ij_ray,ik_ray)         : ye grid index for zone j.
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, numerical_module  
!  eos_snc_x_module!
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : one

USE eos_snc_x_module, ONLY : idr, itr, iyr, idty, dgrid, tgrid, ygrid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: fdp           ! position of rho in grid wrt lower cube index
REAL(KIND=double), INTENT(out)   :: fdm           ! position of rho in grid wrt upper cube index
REAL(KIND=double), INTENT(out)   :: fdd           ! d(fd)/d(rho)
REAL(KIND=double), INTENT(out)   :: ftp           ! position of t in grid wrt lower cube index
REAL(KIND=double), INTENT(out)   :: ftm           ! position of t in grid wrt upper cube index
REAL(KIND=double), INTENT(out)   :: ftt           ! d(ft)/d(t)
REAL(KIND=double), INTENT(out)   :: fyp           ! position of ye in grid wrt lower cube index
REAL(KIND=double), INTENT(out)   :: fym           ! position of ye in grid wrt upper cube index
REAL(KIND=double), INTENT(out)   :: fyy           ! d(fy)/d(ye)

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: id            ! EOS density index
INTEGER                          :: it            ! EOS temperature index
INTEGER                          :: iy            ! EOS electron fraction index

REAL(KIND=double)                :: log_e         ! log10(e)
REAL(KIND=double)                :: fd            ! position of rho in grid
REAL(KIND=double)                :: ft            ! position of t in grid
REAL(KIND=double)                :: fy            ! position of ye in grid

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  log_e            = DLOG10( DEXP(1.d+00) )
END IF

!-----------------------------------------------------------------------
!        Compute independent variable grid indices (id, it, iy)
!         and interpolation coefficients (fdp, fdm, fdd, etc.)
!-----------------------------------------------------------------------

fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10(rho )
ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10(t   )
fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )

id                 = idr(j,ij_ray,ik_ray)
fdp                = fd - DBLE( id )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho )

it                 = itr(j,ij_ray,ik_ray)
ftp                = ft - DBLE( it )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   )

iy                 = iyr(j,ij_ray,ik_ray)
fyp                = fy - DBLE( iy )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

RETURN
END SUBROUTINE eosidx_x
