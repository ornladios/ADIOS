SUBROUTINE edit_Global_mean( imin, imax, nx, time, t_tb, x_in, mass_a_ray, &
&  rho_min, rho_max, rho_mean, sigma_rho, t_min, t_max, t_mean, sigma_t, &
&  u_min, u_max, u_mean, sigma_u, v_min, v_max, v_mean, sigma_v, w_min, &
&  w_max, w_mean, sigma_w, s_min, s_max, s_mean, sigma_s, ye_min, ye_max, &
&  ye_mean, sigma_ye, mach_min, mach_max, mach_mean, sigma_mach, c_shock, &
& nprint )
!-----------------------------------------------------------------------
!
!    File:         edit_Global_mean
!    Module:       edit_Global_mean
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/2/00
!
!    Purpose:
!      To edit the minimum, maximum, and mean values of quantities as a
!      function of the radius.
!
!    Subprograms called:
!  date_and_time_print : fetches and prints the date and time
!
!    Input arguments:
!  imin         : minimum x-array index for the edit
!  imax         : maximum x-array index for the edit
!  nx           : x-array extent
!  time         : elapsed time
!  t_tb         : time from core bounce
!  x_in         : radial midpoint of zone (cm)
!  mass_a_ray   : mass/angular ray (g)
!  rho_min      : minimum density/(radial zone) (g cm^{-3})
!  rho_max      : maximum density/(radial zone) (g cm^{-3})
!  rho_mean     : mass averaged density/radius (g cm^{-3})
!  sigma_rho    : RMS density/(radial zone) (g cm^{-3})
!  t_min        : minimum temperature/(radial zone) (MeV)
!  t_max        : maximum temperature/(radial zone) (MeV)
!  t_mean       : mass averaged temperature/radius (MeV)
!  sigma_t      : RMS temperature/(radial zone) (MeV)
!  u_min        : minimum radial velocity/(radial zone) (cm s^{-1})
!  u_max        : maximum radial velocity/(radial zone) (cm s^{-1})
!  u_mean       : mass averaged radial velocity/radius (cm s^{-1})
!  sigma_u      : RMS radial velocity/(radial zone) (MeV)
!  v_min        : minimum y-velocity/(radial zone) (cm s^{-1})
!  v_max        : maximum y-velocity/(radial zone) (cm s^{-1})
!  v_mean       : mass averaged y-velocity/radius (cm s^{-1})
!  sigma_v      : RMS y-velocity velocity/(radial zone) (MeV)
!  w_min        : minimum z-velocity/(radial zone) (cm s^{-1})
!  w_max        : maximum z-velocity/(radial zone) (cm s^{-1})
!  w_mean       : mass averaged z-velocity/radius (cm s^{-1})
!  sigma_w      : RMS z-velocity velocity/(radial zone) (MeV)
!  s_min        : minimum entropy(radial zone)
!  s_max        : maximum entropy/(radial zone)
!  s_mean       : mass averaged entropy/radius
!  sigma_s      : RMS entropy velocity/(radial zone) (MeV)
!  ye_min       : minimum electron fraction(radial zone)
!  ye_max       : maximum electron fraction/(radial zone)
!  ye_mean      : mass averaged electron fraction/radius
!  sigma_ye     : RMS electron fraction velocity/(radial zone) (MeV)
!  mach_min     : minimum mach number/radius
!  mach_max     : maximum mach number/radius
!  mach_mean    : mass averaged mach number/radius
!  sigma_mach   : mass averaged mach number/radius
!  c_shock      : shock location
!  nprint       : unit numberto print data
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_head, shock_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half
USE physcnst_module, ONLY : msolar

USE edit_module, ONLY : head
USE shock_module, ONLY: j_shk_radial_all_p


IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                              :: imin           ! minimum unshifted radial zone index
INTEGER, INTENT(in)                              :: imax           ! maximum unshifted radial zone index
INTEGER, INTENT(in)                              :: nx             ! x-array extent
INTEGER, INTENT(in)                              :: nprint         ! unit numberto print data

CHARACTER (len=2), INTENT(in), DIMENSION(nx)     :: c_shock        ! shock location

REAL(KIND=double), INTENT(in)                    :: time           ! elapsed time
REAL(KIND=double), INTENT(in)                    :: t_tb           ! time from core bounce
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: x_in           ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mass_a_ray     ! mass/(angular ray) (g)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: rho_mean       ! mass averaged density/radius (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: rho_min        ! minimum density/(radial zone) (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: rho_max        ! maximum density/(radial zone) (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_rho      ! RMS density/(radial zone) (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: t_mean         ! mass averaged temperature/radius (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: t_min          ! minimum temperature/(radial zone) (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: t_max          ! maximum temperature/(radial zone) (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_t        ! RMS temperature/(radial zone) (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: u_mean         ! mass averaged radial velocity/radius (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: u_min          ! minimum radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: u_max          ! maximum radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_u        ! RMS radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: v_mean         ! mass averaged y-velocity/radius (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: v_min          ! minimum y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: v_max          ! maximum y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_v        ! RMS y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: w_mean         ! mass averaged z-velocity/radius (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: w_min          ! minimum z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: w_max          ! maximum z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_w        ! RMS z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: s_mean         ! mass averaged entropy/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: s_min          ! minimum entropy/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: s_max          ! maximum entropy/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_s        ! RMS entropy/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ye_mean        ! mass averaged electron fraction/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ye_min         ! minimum electron fraction/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ye_max         ! maximum electron fraction/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_ye       ! RMS electron fraction/(radial zone) (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mach_min       ! minimum mach number/radius
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mach_max       ! maximum mach number/radius
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mach_mean      ! mass averaged mach number/(radial zone)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: sigma_mach     ! RMS mach number/radius

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i             ! do index

REAL(KIND=double), DIMENSION(nx) :: mass          ! enclosed mass to radial zone center

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1a)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (9x,'Elapsed time=',es14.7,10x,' Time from bounce=',es14.7/)
    9 FORMAT (41x,'Minimum, maximum, and mean data')
   11 FORMAT (40x,33('-')/)
  101 FORMAT ('   i      r            mass          rho      rho_min    rho_max   sigma_rho &
&    t        t_min      t_max      sigma_t       u        u_min      u_max     sigma_t   &
&    v        v_min      v_max     sigma_v       w        w_min      w_max     sigma_w   &
&  mach     mach_min   mach_max  sigma_mach      s        s_min      s_max     sigma_s   &
&   ye       ye_min     ye_max    sigma_ye'/)
  103 FORMAT (1x,i4,es11.3,a2,es15.7,32(es11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                          \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,7) time,t_tb
WRITE (nprint,9)
WRITE (nprint,11)
WRITE (nprint,101)

mass              = zero
mass(imin)        = half * mass_a_ray(imin)
DO i = imin+1,imax
  mass(i)         = mass(i-1) + half * ( mass_a_ray(i-1) + mass_a_ray(i) )
END DO

DO i = imax, imin, -1
  WRITE (nprint,103) i, x_in(i), c_shock(i), mass(i)/msolar, rho_mean(i), &
&  rho_min(i), rho_max(i), sigma_rho(i), t_mean(i), t_min(i),t_max(i), &
&  sigma_t(i), u_mean(i), u_min(i), u_max(i), sigma_u(i), v_mean(i), &
&  v_min(i), v_max(i), sigma_v(i), w_mean(i), w_min(i), w_max(i), sigma_w(i), &
&  mach_mean(i), mach_min(i), mach_max(i), sigma_mach(i), s_mean(i), &
&  s_min(i), s_max(i), sigma_s(i), ye_mean(i), ye_min(i), ye_max(i), &
&  sigma_ye(i)
END DO ! i

RETURN
END SUBROUTINE edit_Global_mean
