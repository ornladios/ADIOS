SUBROUTINE sctnkrnl( n, j, ij_ray, ik_ray, k, rho, t, ye, scatn0, scatn1, &
& scatn0d, scatn1d, scatn0t, scatn1t, scatn0y, scatn1y, nez )
!-----------------------------------------------------------------------
!
!    File:         sctnkrnl
!    Module:       sctnkrnl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/09/03
!
!    Purpose:
!      To interpolate neutrino-nucleon inelastic down- and isoenergetic
!       scattering kernals from a local table of nearest entries created
!       for each zone.
!      The table consists of the eight nearest-neighbor entries in rho, t,
!       and ye of the state point.
!      Derivatives of the scattering rates are obtained from the
!       interpolation formula for these rates by direct differentiation
!       of the interpolation formula.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!
!  n             : neutrino type
!  j             : radial zone number
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!  k             : neutrino energy grouip
!  rho           : density (g/cm**3)
!  t             : temperature (K)
!  ye            : electron fraction
!
!    Output arguments:
!
!  scatn0        : zero moment of the neutrino-nucleon inelastic down or
!                   isoenergetic scattering function
!  scatn1        : first moment of the neutrino-nucleon inelastic down or
!                   isoenergetic scattering function
!  scatn0d       : derivative with respect to density of the zero moment
!                   of the neutrino-nucleon inelastic down or isoenergetic
!                   scattering function
!  scatn1d       : derivative with respect to density of the first moment
!                   of the neutrino-nucleon inelastic down or isoenergetic
!                   scattering function
!  scatn0t       : derivative with respect to temperature of the zero
!                   moment of the neutrino-nucleon inelastic down or
!                   isoenergetic scattering function
!  scatn1t       : derivative with respect to temperature of the first
!                   moment of the neutrino-nucleon inelastic down or
!                  isoenergetic scattering function
!  scatn0y       : derivative with respect to electron fraction of the
!                   zero moment of the neutrino-nucleon inelastic down or
!                  isoenergetic scattering function
!  scatn1y       : derivative with respect to electron fraction of the
!                   first moment of the neutrino-nucleon inelastic down
!                   or isoenergetic cattering function
!
!    Variables that must be passed through common:
!
!  sctn0         : zero moment of the neutrino-nucleon inelastic out
!                   scattering function array
!  dgrid(idty(j,ij_ray,ik_ray))
!                : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray))
!                : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray))
!                : number of table entries in ye between ye = 0.5 and
!                   ye = 0 for zone j
!  idty(j,ij_ray,ik_ray) : index for dgrid, tgrid, and ygrid for zone j
!  cv            : weak interaction coefficient
!  ca            : weak interaction coefficient
!  idrn(j,ij_ray,ik_ray) : rho grid index for zone j
!  itrn(j,ij_ray,ik_ray) : t grid index for zone j
!  iyrn(j,ij_ray,ik_ray) : ye grid index for zone j
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  eos_snc_x_module, nu_energy_grid_module, scat_n_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : zero, half, one
USE physcnst_module, ONLY : cv, ca

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp
USE scat_n_module, ONLY : sctn0, sctnb0, idrn, itrn, iyrn

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: n             ! neutrino flavor index
INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: k             ! incoming neutrino energy zone index
INTEGER, INTENT(in)              :: nez           ! dimension of neutrino energy arrays

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn0  ! zero moment of the NNS function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn1  ! first moment of the NNS function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn0d ! d(scat0)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn1d ! d(scat1)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn0t ! d(scat0)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn1t ! d(scat1)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn0y ! d(scat0)/d(ye)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scatn1y ! d(scat1)/d(ye)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: kp            ! outcoming neutrino energy zone index

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

REAL(KIND=double)                :: sn111         ! scalar table entry for interpolation
REAL(KIND=double)                :: sn211         ! scalar table entry for interpolation
REAL(KIND=double)                :: sn121         ! scalar table entry for interpolation
REAL(KIND=double)                :: sn112         ! scalar table entry for interpolation
REAL(KIND=double)                :: sn221         ! scalar table entry for interpolation
REAL(KIND=double)                :: sn212         ! scalar table entry for interpolation
REAL(KIND=double)                :: sn122         ! scalar table entry for interpolation
REAL(KIND=double)                :: sn222         ! scalar table entry for interpolation

REAL(KIND=double)                :: scatn0l       ! interpolated log of zero-m NNS funciton i
REAL(KIND=double)                :: scatn0dl      ! interpolated log of d(scatn0)/d(rho)
REAL(KIND=double)                :: scatn0tl      ! interpolated log of d(scatn0)/d(t)
REAL(KIND=double)                :: scatn0yl      ! interpolated log of d(scatn0)/d(ye)

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

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho )
id                 = idrn(j,ij_ray,ik_ray)
fdp                = DMAX1( DMIN1( fd - DBLE( id ), one ), zero )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho )

ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   )
it                 = itrn(j,ij_ray,ik_ray)
ftp                = DMAX1( DMIN1( ft - DBLE( it ), one ), zero )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   ) ! log_e = log(e)

fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )
iy                 = iyrn(j,ij_ray,ik_ray)
fyp                = DMAX1( DMIN1( fy - DBLE( iy ), one ), zero )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

DO kp = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store neutrino-nucleon down- and iso-e scattering functions in
!   temporary scalar variables for interpolation
!-----------------------------------------------------------------------

  IF ( n == 1  .or.  n == 3 ) THEN
    sn111          = sctn0 (j,k,kp,ida  ,ita  ,iya  ,i_ray)
    sn211          = sctn0 (j,k,kp,idap1,ita  ,iya  ,i_ray)
    sn121          = sctn0 (j,k,kp,ida  ,itap1,iya  ,i_ray)
    sn112          = sctn0 (j,k,kp,ida  ,ita  ,iyap1,i_ray)
    sn221          = sctn0 (j,k,kp,idap1,itap1,iya  ,i_ray)
    sn212          = sctn0 (j,k,kp,idap1,ita  ,iyap1,i_ray)
    sn122          = sctn0 (j,k,kp,ida  ,itap1,iyap1,i_ray)
    sn222          = sctn0 (j,k,kp,idap1,itap1,iyap1,i_ray)  
  ELSE
    sn111          = sctnb0 (j,k,kp,ida  ,ita  ,iya  ,i_ray)
    sn211          = sctnb0 (j,k,kp,idap1,ita  ,iya  ,i_ray)
    sn121          = sctnb0 (j,k,kp,ida  ,itap1,iya  ,i_ray)
    sn112          = sctnb0 (j,k,kp,ida  ,ita  ,iyap1,i_ray)
    sn221          = sctnb0 (j,k,kp,idap1,itap1,iya  ,i_ray)
    sn212          = sctnb0 (j,k,kp,idap1,ita  ,iyap1,i_ray)
    sn122          = sctnb0 (j,k,kp,ida  ,itap1,iyap1,i_ray)
    sn222          = sctnb0 (j,k,kp,idap1,itap1,iyap1,i_ray)  
  END IF ! n == 1  .or.  n == 3

!-----------------------------------------------------------------------
!  Interpolate NNS down- and iso-e  scattering kernals
!-----------------------------------------------------------------------

  scatn0l          = fym * ( fdm * ( ftm * sn111  + ftp * sn121  )   &
&                  +         fdp * ( ftm * sn211  + ftp * sn221  ) ) &
&                  + fyp * ( fdm * ( ftm * sn112  + ftp * sn122  )   &
&                  +         fdp * ( ftm * sn212  + ftp * sn222  ) )

  scatn0(kp)       = f10(scatn0l)
  scatn1(kp)       = zero

!-----------------------------------------------------------------------
!  Interpolate d(nns down- and iso-e scattering kernal)/d(rho)
!-----------------------------------------------------------------------

  scatn0dl          = fdd * ( fym * ( ftm * ( -sn111  + sn211  )   &
&                   +                 ftp * ( -sn121  + sn221  ) ) &
&                   +         fyp * ( ftm * ( -sn112  + sn212  )   &
&                   +                 ftp * ( -sn122  + sn222  ) ) )

  scatn0d(kp)       = ln_10 * scatn0dl * scatn0(kp)
  scatn1d(kp)       = zero

!-----------------------------------------------------------------------
!  Interpolate d(nns down- and iso-e scattering kernal)/d(t)
!-----------------------------------------------------------------------

  scatn0tl          = ftt * ( fym * ( fdm * ( -sn111  + sn121  )   &
&                   +                 fdp * ( -sn211  + sn221  ) ) &
&                   +         fyp * ( fdm * ( -sn112  + sn122  )   &
&                   +                 fdp * ( -sn212  + sn222  ) ) )

  scatn0t(kp)       = ln_10 * scatn0tl * scatn0(kp)
  scatn1t(kp)       = zero

!-----------------------------------------------------------------------
!  Interpolate d(nns down- and iso-e scattering kernal)/d(ye)
!-----------------------------------------------------------------------

  scatn0yl          = fyy * ( fdm * ( ftm * ( -sn111  + sn112  )   &
&                   +                 ftp * ( -sn121  + sn122  ) ) &
&                   +         fdp * ( ftm * ( -sn211  + sn212  )   &
&                   +                 ftp * ( -sn221  + sn222  ) ) )

  scatn0y(kp)       = ln_10 * scatn0yl * scatn0(kp)
  scatn1y(kp)       = zero

END DO

RETURN
END SUBROUTINE sctnkrnl
