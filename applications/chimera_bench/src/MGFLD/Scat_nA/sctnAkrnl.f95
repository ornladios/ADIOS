SUBROUTINE sctnAkrnl( n, j, ij_ray, ik_ray, k, rho, t, ye, scata0_out, &
& scata1_out, scata0d_out, scata1d_out, scata0t_out, scata1t_out, &
& scata0y_out, scata1y_out, scata0_in, scata1_in, scata0d_in, scata1d_in, &
& scata0t_in, scata1t_in, scata0y_in, scata1y_in, nez )
!-----------------------------------------------------------------------
!
!    File:         sctnAkrnl
!    Module:       sctnAkrnl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/28/05
!
!    Purpose:
!      To interpolate neutrino-nucleus inelastic down- and isoenergetic
!       scattering kernals from a local table of nearest entries created
!       for each zone.
!      The table consists of the eight nearest-neighbor entries in rho,
!       t, and ye of the state point.
!      Derivatives of the scattering rates are obtained from the interpolation
!       formula for these rates by direct differentiation of the interpolation
!       formula.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!
!  n              : neutrino type
!  j              : radial zone number
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!  k              : neutrino energy grouip
!  rho            : density (g/cm**3)
!  t              : temperature (K)
!  ye             : electron fraction
!
!    Output arguments:
!
!  scata0_out     : zero moment of the neutrino-nucleus inelastic out
!                    scattering function
!  scata1_out     : first moment of the neutrino-nucleus inelastic out
!                    scattering function
!  scata0d_out    : derivative with respect to density of the zero moment of the
!                    neutrino-nucleus inelastic out scattering function
!  scata1d_out    : derivative with respect to density of the first moment of the
!                    neutrino-nucleus inelastic out scattering function
!  scata0t_out    : derivative with respect to temperature of the zero moment of the
!                    neutrino-nucleus inelastic out scattering function
!  scata1t_out    : derivative with respect to temperature of the first moment of the
!                    neutrino-nucleus inelastic out scattering function
!  scata0y_out    : derivative with respect to electron fraction of the zero moment of the
!                    neutrino-nucleus inelastic out scattering function
!  scata1y_out    : derivative with respect to electron fraction of the first moment of
!                    the neutrino-nucleus inelastic out cattering function
!  scata0_in      : zero moment of the neutrino-nucleus inelastic in
!                    scattering function
!  scata1_in      : first moment of the neutrino-nucleus inelastic in
!                    scattering function
!  scata0d_in     : derivative with respect to density of the zero moment of the
!                    neutrino-nucleus inelastic in scattering function
!  scata1d_in     : derivative with respect to density of the first moment of the
!                    neutrino-nucleus inelastic in scattering function
!  scata0t_in     : derivative with respect to temperature of the zero moment of the
!                    neutrino-nucleus inelastic in scattering function
!  scata1t_in     : derivative with respect to temperature of the first moment of the
!                    neutrino-nucleus inelastic in scattering function
!  scata0y_in     : derivative with respect to electron fraction of the zero moment of the
!                    neutrino-nucleus inelastic in scattering function
!  scata1y_in     : derivative with respect to electron fraction of the first moment of
!                    the neutrino-nucleus inelastic in cattering function
!
!    Variables that must be passed through common:
!
!  sctnA0          : zero moment of the neutrino-nucleus inelastic out scattering function array
!  dgrid(idty(j,ij_ray,ik_ray))
!                 : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray))
!                 : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray))
!                 : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)
!                 : index for dgrid, tgrid, and ygrid for zone j
!  cv             : weak interaction coefficient
!  ca             : weak interaction coefficient
!  idrnA(j,ij_ray,ik_ray)
!                 : rho grid index for zone j
!  itrnA(j,ij_ray,ik_ray)
!                 : t grid index for zone j
!  iyrnA(j,ij_ray,ik_ray)
!                 : ye grid index for zone j
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
USE scat_nA_module, ONLY : sctnA0, idrnA, itrnA, iyrnA

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

REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0_out  ! zero moment of the NNS function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1_out  ! first moment of the NNS function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0d_out ! d(scat0)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1d_out ! d(scat1)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0t_out ! d(scat0)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1t_out ! d(scat1)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0y_out ! d(scat0)/d(ye)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1y_out ! d(scat1)/d(ye)

REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0_in   ! zero moment of the NNS function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1_in   ! first moment of the NNS function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0d_in  ! d(scat0)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1d_in  ! d(scat1)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0t_in  ! d(scat0)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1t_in  ! d(scat1)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata0y_in  ! d(scat0)/d(ye)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scata1y_in  ! d(scat1)/d(ye)

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

REAL(KIND=double)                :: sa111         ! scalar table entry for interpolation
REAL(KIND=double)                :: sa211         ! scalar table entry for interpolation
REAL(KIND=double)                :: sa121         ! scalar table entry for interpolation
REAL(KIND=double)                :: sa112         ! scalar table entry for interpolation
REAL(KIND=double)                :: sa221         ! scalar table entry for interpolation
REAL(KIND=double)                :: sa212         ! scalar table entry for interpolation
REAL(KIND=double)                :: sa122         ! scalar table entry for interpolation
REAL(KIND=double)                :: sa222         ! scalar table entry for interpolation

REAL(KIND=double)                :: scata0l       ! interpolated log of zero-m NNS funciton i
REAL(KIND=double)                :: scata0dl      ! interpolated log of d(scata0)/d(rho)
REAL(KIND=double)                :: scata0tl      ! interpolated log of d(scata0)/d(t)
REAL(KIND=double)                :: scata0yl      ! interpolated log of d(scata0)/d(ye)

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
id                 = idrnA(j,ij_ray,ik_ray)
fdp                = DMAX1( DMIN1( fd - DBLE( id ), one ), zero )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho )

ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   )
it                 = itrnA(j,ij_ray,ik_ray)
ftp                = DMAX1( DMIN1( ft - DBLE( it ), one ), zero )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   ) ! log_e = log(e)

fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )
iy                 = iyrnA(j,ij_ray,ik_ray)
fyp                = DMAX1( DMIN1( fy - DBLE( iy ), one ), zero )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

!-----------------------------------------------------------------------
!
!         \\\\\ NEUTRINO-NUCLEUS OUT-SCATTERING KERNELS /////
!
!-----------------------------------------------------------------------

DO kp = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store neutrino-nucleus out scattering functions in temporary scalar
!   variables for interpolation
!-----------------------------------------------------------------------

  sa111            = sctnA0 (j,k,kp,ida  ,ita  ,iya  ,i_ray)
  sa211            = sctnA0 (j,k,kp,idap1,ita  ,iya  ,i_ray)
  sa121            = sctnA0 (j,k,kp,ida  ,itap1,iya  ,i_ray)
  sa112            = sctnA0 (j,k,kp,ida  ,ita  ,iyap1,i_ray)
  sa221            = sctnA0 (j,k,kp,idap1,itap1,iya  ,i_ray)
  sa212            = sctnA0 (j,k,kp,idap1,ita  ,iyap1,i_ray)
  sa122            = sctnA0 (j,k,kp,ida  ,itap1,iyap1,i_ray)
  sa222            = sctnA0 (j,k,kp,idap1,itap1,iyap1,i_ray)  

!-----------------------------------------------------------------------
!  Interpolate NAS out scattering kernals
!-----------------------------------------------------------------------

  scata0l          = fym * ( fdm * ( ftm * sa111  + ftp * sa121  )   &
&                  +         fdp * ( ftm * sa211  + ftp * sa221  ) ) &
&                  + fyp * ( fdm * ( ftm * sa112  + ftp * sa122  )   &
&                  +         fdp * ( ftm * sa212  + ftp * sa222  ) )

  scata0_out(kp)   = f10(scata0l)
  scata1_out(kp)   = zero

!-----------------------------------------------------------------------
!  Interpolate d(NAS out scattering kernal)/d(rho)
!-----------------------------------------------------------------------

  scata0dl          = fdd * ( fym * ( ftm * ( -sa111  + sa211  )   &
&                   +                 ftp * ( -sa121  + sa221  ) ) &
&                   +         fyp * ( ftm * ( -sa112  + sa212  )   &
&                   +                 ftp * ( -sa122  + sa222  ) ) )

  scata0d_out(kp)   = ln_10 * scata0dl * scata0_out(kp)
  scata1d_out(kp)   = zero

!-----------------------------------------------------------------------
!  Interpolate d(NAS out scattering kernal)/d(t)
!-----------------------------------------------------------------------

  scata0tl          = ftt * ( fym * ( fdm * ( -sa111  + sa121  )   &
&                   +                 fdp * ( -sa211  + sa221  ) ) &
&                   +         fyp * ( fdm * ( -sa112  + sa122  )   &
&                   +                 fdp * ( -sa212  + sa222  ) ) )

  scata0t_out(kp)   = ln_10 * scata0tl * scata0_out(kp)
  scata1t_out(kp)   = zero

!-----------------------------------------------------------------------
!  Interpolate d(NAS out scattering kernal)/d(ye)
!-----------------------------------------------------------------------

  scata0yl          = fyy * ( fdm * ( ftm * ( -sa111  + sa112  )   &
&                   +                 ftp * ( -sa121  + sa122  ) ) &
&                   +         fdp * ( ftm * ( -sa211  + sa212  )   &
&                   +                 ftp * ( -sa221  + sa222  ) ) )

  scata0y_out(kp)   = ln_10 * scata0yl * scata0_out(kp)
  scata1y_out(kp)   = zero

END DO

!-----------------------------------------------------------------------
!
!         \\\\\ NEUTRINO-NUCLEUS IN-SCATTERING KERNELS /////
!
!-----------------------------------------------------------------------

DO kp = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store neutrino-nucleus in scattering functions in temporary scalar
!   variables for interpolation
!-----------------------------------------------------------------------

  sa111            = sctnA0 (j,kp,k,ida  ,ita  ,iya  ,i_ray)
  sa211            = sctnA0 (j,kp,k,idap1,ita  ,iya  ,i_ray)
  sa121            = sctnA0 (j,kp,k,ida  ,itap1,iya  ,i_ray)
  sa112            = sctnA0 (j,kp,k,ida  ,ita  ,iyap1,i_ray)
  sa221            = sctnA0 (j,kp,k,idap1,itap1,iya  ,i_ray)
  sa212            = sctnA0 (j,kp,k,idap1,ita  ,iyap1,i_ray)
  sa122            = sctnA0 (j,kp,k,ida  ,itap1,iyap1,i_ray)
  sa222            = sctnA0 (j,kp,k,idap1,itap1,iyap1,i_ray)  

!-----------------------------------------------------------------------
!  Interpolate NAS out scattering kernals
!-----------------------------------------------------------------------

  scata0l          = fym * ( fdm * ( ftm * sa111  + ftp * sa121  )   &
&                  +         fdp * ( ftm * sa211  + ftp * sa221  ) ) &
&                  + fyp * ( fdm * ( ftm * sa112  + ftp * sa122  )   &
&                  +         fdp * ( ftm * sa212  + ftp * sa222  ) )

  scata0_in(kp)    = f10(scata0l)
  scata1_in(kp)    = zero

!-----------------------------------------------------------------------
!  Interpolate d(NAS out scattering kernal)/d(rho)
!-----------------------------------------------------------------------

  scata0dl          = fdd * ( fym * ( ftm * ( -sa111  + sa211  )   &
&                   +                 ftp * ( -sa121  + sa221  ) ) &
&                   +         fyp * ( ftm * ( -sa112  + sa212  )   &
&                   +                 ftp * ( -sa122  + sa222  ) ) )

  scata0d_in(kp)    = ln_10 * scata0dl * scata0_in(kp)
  scata1d_in(kp)    = zero

!-----------------------------------------------------------------------
!  Interpolate d(NAS out scattering kernal)/d(t)
!-----------------------------------------------------------------------

  scata0tl          = ftt * ( fym * ( fdm * ( -sa111  + sa121  )   &
&                   +                 fdp * ( -sa211  + sa221  ) ) &
&                   +         fyp * ( fdm * ( -sa112  + sa122  )   &
&                   +                 fdp * ( -sa212  + sa222  ) ) )

  scata0t_in(kp)    = ln_10 * scata0tl * scata0_in(kp)
  scata1t_in(kp)    = zero

!-----------------------------------------------------------------------
!  Interpolate d(NAS out scattering kernal)/d(ye)
!-----------------------------------------------------------------------

  scata0yl          = fyy * ( fdm * ( ftm * ( -sa111  + sa112  )   &
&                   +                 ftp * ( -sa121  + sa122  ) ) &
&                   +         fdp * ( ftm * ( -sa211  + sa212  )   &
&                   +                 ftp * ( -sa221  + sa222  ) ) )

  scata0y_in(kp)    = ln_10 * scata0yl * scata0_in(kp)
  scata1y_in(kp)    = zero

END DO

RETURN
END SUBROUTINE sctnAkrnl
