SUBROUTINE sctekrnl( n, j, ij_ray, ik_ray, k, rho, t, ye, scat0, scat1, &
& scat0d, scat1d, scat0t, scat1t, scat0y, scat1y, nez )
!-----------------------------------------------------------------------
!
!    File:         sctekrnl
!    Module:       sctekrnl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/01/03
!
!    Purpose:
!      To interpolate neutrino-electron down- and isoenergetic scattering kernals
!       from a local table of nearest entries created for each zone.
!      The table consists of the eight nearest-neighbor entries in
!       rho, t, and ye of the state point.
!      Derivatives of the scattering rates are obtained from the interpolation formula
!       for these rates by direct differentiation of the interpolation formula.
!
!    Subprograms called:
!      none
!
!    Input arguments:
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
!  scat0          : zero moment of the NES down or isoenergetic scattering function
!  scat1          : first moment of the NES down or isoenergetic scattering function
!  scat0d         : derivative with respect to density of the zero moment
!                    of the NES down or isoenergetic scattering function
!  scat1d         : derivative with respect to density of the first moment
!                    of the NES down or isoenergetic scattering function
!  scat0t         : derivative with respect to temperature of the zero moment
!                    of the NES down or isoenergetic scattering function
!  scat1t         : derivative with respect to temperature of the first moment
!                    of the NES down or isoenergetic scattering function
!  scat0y         : derivative with respect to electron fraction of the zero moment
!                    of the NES down or isoenergetic scattering function
!  scat0y         : derivative with respect to electron fraction of the first moment
!                    of the NES down or isoenergetic scattering function
!
!    Variables that must be passed through common:
!  sctin0         : zero moment of the nonconservative neutrino in-scattering kernal array
!  sctot0         : zero moment of the nonconservative neutrino out-scattering kernal array
!  sctin1         : first moment of the nonconservative neutrino in-scattering kernal array
!  sctot1         : first moment of the nonconservative neutrino out-scattering kernal array
!  dgrid(idty(j)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j)) : number of table entries per decade in t for zone j
!  ygrid(idty(j)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j)        : index for dgrid, tgrid, and ygrid for zone j
!  cv             : weak interaction coefficient
!  ca             : weak interaction coefficient
!  idrse(j,ij_ray,ik_ray) : rho grid index for zone j
!  itrse(j,ij_ray,ik_ray) : t grid index for zone j
!  iyrse(j,ij_ray,ik_ray) : ye grid index for zone j
!
!        The variables that should be saved are
!
!  cnes1(nnu)     : first nes coefficiant for neutrinos of type nnu
!  cnes2(nnu)     : second nes coefficiant for neutrinos of type nnu
!
!    Include files:
!  array_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, nu_energy_grid_module, scat_e_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : zero, half, one
USE physcnst_module, ONLY : cv, ca

USE edit_module, ONLY : nprint, nlog
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp
USE scat_e_module, ONLY : scte0i, scte0ii, idrse, itrse, iyrse

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
INTEGER, INTENT(in)              :: nez           ! neutrino energy array dimension

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat0  ! zero moment of the NES function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat1  ! first moment of the NES function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat0d ! d(scat0)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat1d ! d(scat1)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat0t ! d(scat0)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat1t ! d(scat1)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat0y ! d(scat0)/d(ye)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: scat1y ! d(scat1)/d(ye)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                            :: first = .true.

INTEGER                            :: kp            ! outcoming neutrino energy zone index

INTEGER                            :: id            ! EOS density index
INTEGER                            :: it            ! EOS temperature index
INTEGER                            :: iy            ! EOS electron fraction index

INTEGER                            :: i_ray         ! f(ij_ray,ik_ray)

INTEGER, PARAMETER                 :: ida = 1       ! lower cube table density index
INTEGER, PARAMETER                 :: idap1 = 2     ! upper cube table density index
INTEGER, PARAMETER                 :: ita = 1       ! lower cube table temperature index
INTEGER, PARAMETER                 :: itap1 = 2     ! upper cube table temperature index
INTEGER, PARAMETER                 :: iya = 1       ! lower cube table electron fraction index
INTEGER, PARAMETER                 :: iyap1 = 2     ! upper cube table electron fraction index

INTEGER, PARAMETER                 :: nnud=4        ! upper cube table electron fraction index

REAL(KIND=double), DIMENSION(nnud) :: cnes1         ! neutrino flavor coefficient
REAL(KIND=double), DIMENSION(nnud) :: cnes2         ! neutrino flavor coefficient

REAL(KIND=double)                  :: log_e         ! log10(e)
REAL(KIND=double)                  :: ln_10         ! ln(10)

REAL(KIND=double)                  :: fd            ! position of rho in grid
REAL(KIND=double)                  :: fdp           ! position of rho in grid wrt lower cube index
REAL(KIND=double)                  :: fdm           ! position of rho in grid wrt upper cube index
REAL(KIND=double)                  :: fdd           ! d(fd)/d(rho)
REAL(KIND=double)                  :: ft            ! position of t in grid
REAL(KIND=double)                  :: ftp           ! position of t in grid wrt lower cube index
REAL(KIND=double)                  :: ftm           ! position of t in grid wrt upper cube index
REAL(KIND=double)                  :: ftt           ! d(ft)/d(t)
REAL(KIND=double)                  :: fy            ! position of ye in grid
REAL(KIND=double)                  :: fyp           ! position of ye in grid wrt lower cube index
REAL(KIND=double)                  :: fym           ! position of ye in grid wrt upper cube index
REAL(KIND=double)                  :: fyy           ! d(fy)/d(ye)

REAL(KIND=double)                  :: si111         ! scalar table entry for interpolation
REAL(KIND=double)                  :: si211         ! scalar table entry for interpolation
REAL(KIND=double)                  :: si121         ! scalar table entry for interpolation
REAL(KIND=double)                  :: si112         ! scalar table entry for interpolation
REAL(KIND=double)                  :: si221         ! scalar table entry for interpolation
REAL(KIND=double)                  :: si212         ! scalar table entry for interpolation
REAL(KIND=double)                  :: si122         ! scalar table entry for interpolation
REAL(KIND=double)                  :: si222         ! scalar table entry for interpolation

REAL(KIND=double)                  :: sii111        ! scalar table entry for interpolation
REAL(KIND=double)                  :: sii211        ! scalar table entry for interpolation
REAL(KIND=double)                  :: sii121        ! scalar table entry for interpolation
REAL(KIND=double)                  :: sii112        ! scalar table entry for interpolation
REAL(KIND=double)                  :: sii221        ! scalar table entry for interpolation
REAL(KIND=double)                  :: sii212        ! scalar table entry for interpolation
REAL(KIND=double)                  :: sii122        ! scalar table entry for interpolation
REAL(KIND=double)                  :: sii222        ! scalar table entry for interpolation

REAL(KIND=double)                  :: scate0il      ! interpolated log of zero-m NES funciton i
REAL(KIND=double)                  :: scate0iil     ! interpolated log of zero-m NES funciton ii
REAL(KIND=double)                  :: scate0i       ! interpolated zero-m NES funciton i
REAL(KIND=double)                  :: scate0ii      ! interpolated zero-m NES funciton ii
REAL(KIND=double)                  :: scate1i       ! interpolated first-m NES funciton i
REAL(KIND=double)                  :: scate1ii      ! interpolated first-m NES funciton ii

REAL(KIND=double)                  :: scate0idl     ! interpolated log of d(scte0i)/d(rho)
REAL(KIND=double)                  :: scate0iidl    ! interpolated log of d(scte0ii)/d(rho)
REAL(KIND=double)                  :: scate0id      ! interpolated d(scte0i)/d(rho)
REAL(KIND=double)                  :: scate0iid     ! interpolated d(scte0ii)/d(rho)
REAL(KIND=double)                  :: scate1id      ! interpolated d(scte1i)/d(rho)
REAL(KIND=double)                  :: scate1iid     ! interpolated d(scte1ii)/d(rho)

REAL(KIND=double)                  :: scate0itl     ! interpolated log of d(scte0i)/d(t)
REAL(KIND=double)                  :: scate0iitl    ! interpolated log of d(scte0ii)/d(t)
REAL(KIND=double)                  :: scate0it      ! interpolated d(scte0i)/d(t)
REAL(KIND=double)                  :: scate0iit     ! interpolated d(scte0ii)/d(t)
REAL(KIND=double)                  :: scate1it      ! interpolated d(scte1i)/d(t)
REAL(KIND=double)                  :: scate1iit     ! interpolated d(scte1ii)/d(t)

REAL(KIND=double)                  :: scate0iyl     ! interpolated log of d(scte0i)/d(ye)
REAL(KIND=double)                  :: scate0iiyl    ! interpolated log of d(scte0ii)/d(ye)
REAL(KIND=double)                  :: scate0iy      ! interpolated d(scte0i)/d(ye)
REAL(KIND=double)                  :: scate0iiy     ! interpolated d(scte0ii)/d(ye)
REAL(KIND=double)                  :: scate1iy      ! interpolated d(scte1i)/d(ye)
REAL(KIND=double)                  :: scate1iiy     ! interpolated d(scte1ii)/d(ye)

REAL(KIND=double)                 :: f10           ! 10** function
EXTERNAL f10

!-----------------------------------------------------------------------
!        Format
!-----------------------------------------------------------------------

 1001 FORMAT (' nnu=',i3,' exceeds nnud=',i3,' in subroutine sctekrnl')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute standard model coefficients of the neutrino-electron
!   scattering functions.
!-----------------------------------------------------------------------

IF (first) THEN
  first            = .false.

  IF ( nnu > nnud ) THEN
    WRITE (nprint,1001) nnu,nnud
    WRITE (nlog,1001) nnu,nnud
    STOP
  END IF

  cnes1(1)         = ( cv + ca )**2
  cnes2(1)         = ( cv - ca )**2
  IF ( nnu >= 2 ) THEN
    cnes1(2)       = cnes2(1)
    cnes2(2)       = cnes1(1)
  END IF
  IF ( nnu >= 3 ) THEN
    cnes1(3)       = ( cv + ca - 2.d+0 )**2
    cnes2(3)       = ( cv - ca )**2
  END IF
  IF ( nnu >= 4 ) THEN
    cnes1(4)       = cnes2(3)
    cnes2(4)       = cnes1(3)
  END IF

  log_e            = DLOG10( DEXP( one ) )
  ln_10            = DLOG( 1.d+01 )

END IF

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho )
id                 = idrse(j,ij_ray,ik_ray)
fdp                = DMAX1( DMIN1( fd - DBLE( id ), one ), zero )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho )

ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   )
it                 = itrse(j,ij_ray,ik_ray)
ftp                = DMAX1( DMIN1( ft - DBLE( it ), one ), zero )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   ) ! log_e = log(e)

fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )
iy                 = iyrse(j,ij_ray,ik_ray)
fyp                = DMAX1( DMIN1( fy - DBLE( iy ), one ), zero )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

DO kp = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store neutrino-electron down- and iso-e scattering functions in
!   temporary scalar variables for interpolation
!-----------------------------------------------------------------------

  si111            = scte0i (j,k,kp,ida  ,ita  ,iya  ,i_ray)
  si211            = scte0i (j,k,kp,idap1,ita  ,iya  ,i_ray)
  si121            = scte0i (j,k,kp,ida  ,itap1,iya  ,i_ray)
  si112            = scte0i (j,k,kp,ida  ,ita  ,iyap1,i_ray)
  si221            = scte0i (j,k,kp,idap1,itap1,iya  ,i_ray)
  si212            = scte0i (j,k,kp,idap1,ita  ,iyap1,i_ray)
  si122            = scte0i (j,k,kp,ida  ,itap1,iyap1,i_ray)
  si222            = scte0i (j,k,kp,idap1,itap1,iyap1,i_ray)

  sii111           = scte0ii(j,k,kp,ida  ,ita  ,iya  ,i_ray)
  sii211           = scte0ii(j,k,kp,idap1,ita  ,iya  ,i_ray)
  sii121           = scte0ii(j,k,kp,ida  ,itap1,iya  ,i_ray)
  sii112           = scte0ii(j,k,kp,ida  ,ita  ,iyap1,i_ray)
  sii221           = scte0ii(j,k,kp,idap1,itap1,iya  ,i_ray)
  sii212           = scte0ii(j,k,kp,idap1,ita  ,iyap1,i_ray)
  sii122           = scte0ii(j,k,kp,ida  ,itap1,iyap1,i_ray)
  sii222           = scte0ii(j,k,kp,idap1,itap1,iyap1,i_ray)

!-----------------------------------------------------------------------
!  Interpolate NES down- and iso-e  scattering kernals
!-----------------------------------------------------------------------

  scate0il         = fym * ( fdm * ( ftm * si111  + ftp * si121  )   &
&                  +         fdp * ( ftm * si211  + ftp * si221  ) ) &
&                  + fyp * ( fdm * ( ftm * si112  + ftp * si122  )   &
&                  +         fdp * ( ftm * si212  + ftp * si222  ) )

  scate0i          = f10(scate0il)
  scate1i          = zero

  scate0iil        = fym * ( fdm * ( ftm * sii111 + ftp * sii121 )   &
&                  +         fdp * ( ftm * sii211 + ftp * sii221 ) ) &
&                  + fyp * ( fdm * ( ftm * sii112 + ftp * sii122 )   &
&                  +         fdp * ( ftm * sii212 + ftp * sii222 ) )

  scate0ii         = f10(scate0iil)
  scate1ii         = zero

!-----------------------------------------------------------------------
!  Interpolate d(nes down- and iso-e scattering kernal)/d(rho)
!-----------------------------------------------------------------------

  scate0idl        = fdd * ( fym * ( ftm * ( -si111  + si211  )   &
&                  +                 ftp * ( -si121  + si221  ) ) &
&                  +         fyp * ( ftm * ( -si112  + si212  )   &
&                  +                 ftp * ( -si122  + si222  ) ) )

  scate0id         = ln_10 * scate0idl * scate0i
  scate1id         = zero

  scate0iidl       = fdd * ( fym * ( ftm * ( -sii111 + sii211 )   &
&                  +                 ftp * ( -sii121 + sii221 ) ) &
&                  +         fyp * ( ftm * ( -sii112 + sii212 )   &
&                  +                 ftp * ( -sii122 + sii222 ) ) )

  scate0iid        = ln_10 * scate0iidl * scate0ii
  scate1iid        = zero

!-----------------------------------------------------------------------
!  Interpolate d(nes down- and iso-e scattering kernal)/d(t)
!-----------------------------------------------------------------------

  scate0itl        = ftt * ( fym * ( fdm * ( -si111  + si121  )   &
&                  +                 fdp * ( -si211  + si221  ) ) &
&                  +         fyp * ( fdm * ( -si112  + si122  )   &
&                  +                 fdp * ( -si212  + si222  ) ) )

  scate0it         = ln_10 * scate0itl * scate0i
  scate1it         = zero

  scate0iitl       = ftt * ( fym * ( fdm * ( -sii111 + sii121 )   &
&                  +                 fdp * ( -sii211 + sii221 ) ) &
&                  +         fyp * ( fdm * ( -sii112 + sii122 )   &
&                  +                 fdp * ( -sii212 + sii222 ) ) )

  scate0iit        = ln_10 * scate0iitl * scate0ii
  scate1iit        = zero

!-----------------------------------------------------------------------
!  Interpolate d(nes down- and iso-e scattering kernal)/d(ye)
!-----------------------------------------------------------------------

  scate0iyl        = fyy * ( fdm * ( ftm * ( -si111  + si112  )   &
&                  +                 ftp * ( -si121  + si122  ) ) &
&                  +         fdp * ( ftm * ( -si211  + si212  )   &
&                  +                 ftp * ( -si221  + si222  ) ) )

  scate0iy         = ln_10 * scate0iyl * scate0i
  scate1iy         = zero

  scate0iiyl       = fyy * ( fdm * ( ftm * ( -sii111 + sii112 )   &
&                  +                 ftp * ( -sii121 + sii122 ) ) &
&                  +         fdp * ( ftm * ( -sii211 + sii212 )   &
&                  +                 ftp * ( -sii221 + sii222 ) ) )

  scate0iiy        = ln_10 * scate0iiyl * scate0ii
  scate1iiy        = zero

!-----------------------------------------------------------------------
!  Construct net scattering functions
!-----------------------------------------------------------------------

  scat0(kp)        = cnes1(n) * scate0i  + cnes2(n) * scate0ii
  scat0d(kp)       = cnes1(n) * scate0id + cnes2(n) * scate0iid
  scat0t(kp)       = cnes1(n) * scate0it + cnes2(n) * scate0iit
  scat0y(kp)       = cnes1(n) * scate0iy + cnes2(n) * scate0iiy
  scat1(kp)        = cnes1(n) * scate1i  + cnes2(n) * scate1ii
  scat1d(kp)       = cnes1(n) * scate1id + cnes2(n) * scate1iid
  scat1t(kp)       = cnes1(n) * scate1it + cnes2(n) * scate1iit
  scat1y(kp)       = cnes1(n) * scate1iy + cnes2(n) * scate1iiy
  
END DO !  kp = 1,nnugp(n)

RETURN
END SUBROUTINE sctekrnl
