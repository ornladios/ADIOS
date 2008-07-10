SUBROUTINE eqstld( jr_min, jr_max, rho, t, yl, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         eqstld
!    Module:       eqstld
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/14/02
!
!    Purpose:
!      To interpolate the logs of all equation of state quantities
!       from a local table of nearest entries created for each
!       radial shell midpoint, and to interpolate their derivatives
!       with respect to rho, t, and yl from nearest neighbor entries
!       of the derivatives. This subroutine has been optimized for a
!       vector machine. This subroutine shold be called at the end
!       of each cycle of the code to update all the thermodynamic
!       variables.
!
!    Variables that must be passed through common:
!
!  dgrid(idty(j,ij_ray,ik_ray))
!                : number of table entries per decade in rho
!                   for zone j.
!  tgrid(idty(j,ij_ray,ik_ray))
!                : number of table entries per decade in t
!                   for zone j.
!  ygrid(idty(j,ij_ray,ik_ray))
!                : number of table entries in ye between
!                   ye = 0.5 and ye = 0 for zone j.
!  idty(j,ij_ray,ik_ray) : index for dgrid, tgrid, and ygrid for zone j.
!  rhoes(i)      : density boundary between idty=i and idty=i+1
!  idr(j,ij_ray,ik_ray)  : rho grid index for zone j.
!  itr(j,ij_ray,ik_ray)  : t grid index for zone j.
!  iyr(j,ij_ray,ik_ray)  : ye grid index for zone j.
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray):
!               : equation of state table array
!  estbled(i,j,ida,ita,iya,ij_ray,ik_ray):
!               : equation of state dervatives with respect to rho
!                  table array.
!  estblet(i,j,ida,ita,iya,ij_ray,ik_ray):
!               : equation of state dervatives with respect to t
!                  table array.
!  estbley(i,j,ida,ita,iya,ij_ray,ik_ray):
!               : equation of state dervatives with respect to yl
!                  table array.
!  escnst(i,j,ij_ray,ik_ray) 
!               : additive constant for estble
!  escnstd(i,j,ij_ray,ik_ray)
!               : additive constant for estbled
!  escnstt(i,j,ij_ray,ik_ray)
!               : additive constant for estblet
!  escnsty(i,j,ij_ray,ik_ray)
!               : additive constant for estbley
!  i            : thermodynamic function index
!  i = 1        : pressure
!  i = 2        : energy
!  i = 3        : entropy
!
!    Subprograms called:
!        none
!
!
!    Input arguments:
!  jr_min        : minimum radial zone for which thermodynamic
!                   variables are to be evaluated
!  jr_max        : maximum radial zone for which thermodynamic
!                   variables are to be evaluated
!  rho           : shifted matter density array (g/cm**3).
!  t             : shifted matter matter temperature array (K).
!  yl            : shifted matter matter electron fraction array.
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  eos_drv_module, edit_module, eos_snc_x_module, lep_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx
USE numerical_module, ONLY : zero, one, epsilon

USE eos_drv_module, ONLY: idr, itr, iyr, escnst, estble, estbled, estblet, &
& estbley, escnstd, escnstt, escnsty, aesv, aesvd, aesvt, aesvy
USE edit_module, ONLY : nprint, nlog
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min          ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max          ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t   ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: yl  ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: j             ! radial zone index
INTEGER                          :: i             ! equation of state variable index
INTEGER                          :: istat         ! allocate file flag

INTEGER                          :: id            ! density index
INTEGER                          :: it            ! temperature index
INTEGER                          :: iy            ! lepton fraction index
INTEGER, PARAMETER               :: ida = 1       ! density lower table index
INTEGER, PARAMETER               :: idap1 = 2     ! density upper table index
INTEGER, PARAMETER               :: ita = 1       ! temperature lower table index index
INTEGER, PARAMETER               :: itap1 = 2     ! temperature upper table index index
INTEGER, PARAMETER               :: iya = 1       ! lepton fraction lower table index index
INTEGER, PARAMETER               :: iyap1 = 2     ! lepton fraction upper table index index

REAL(KIND=double)                :: av111         ! table entry for interpolation
REAL(KIND=double)                :: av211         ! table entry for interpolation
REAL(KIND=double)                :: av121         ! table entry for interpolation
REAL(KIND=double)                :: av112         ! table entry for interpolation
REAL(KIND=double)                :: av221         ! table entry for interpolation
REAL(KIND=double)                :: av212         ! table entry for interpolation
REAL(KIND=double)                :: av122         ! table entry for interpolation
REAL(KIND=double)                :: av222         ! table entry for interpolation
REAL(KIND=double)                :: besv          ! interpolated quantity in the log
REAL(KIND=double)                :: besvd         ! interpolated d(quantity)/d(rho) in the log
REAL(KIND=double)                :: besvt         ! interpolated d(quantity)/d(t) in the log
REAL(KIND=double)                :: besvy         ! interpolated d(quantity)/d(yl) in the log
REAL(KIND=double)                :: esvtmp        ! interpolated quantity
REAL(KIND=double)                :: esvtmpd       ! interpolated d(quantity)/d(rho)
REAL(KIND=double)                :: esvtmpt       ! interpolated d(quantity)/d(t)
REAL(KIND=double)                :: esvtmpy       ! interpolated d(quantity)/d(yl)

REAL(KIND=double)                :: fd            ! location of log(rho) in the table
REAL(KIND=double)                :: ft            ! location of log(t) in the table
REAL(KIND=double)                :: fy            ! location of yl in the table
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fdm  ! distance of log(rho) from top of cube
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fdp  ! distance of log(rho) from bottom of cube
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ftm  ! distance of log(t) from top of cube
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ftp  ! distance of log(t) from bottom of cube
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fym  ! distance of ye from top of cube
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fyp  ! distance of ye from bottom of cube

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in eqstld')
 2001 FORMAT (' Deallocation problem for array ',a10,' in eqstld')

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (fdp(nx), STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fdp       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fdm(nx), STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fdm       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ftp(nx), STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'ftp       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ftm(nx), STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'ftm       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fyp(nx), STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fyp       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fym(nx), STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fym       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

fdp                  = zero
fdm                  = zero
ftp                  = zero
ftm                  = zero
fyp                  = zero
fym                  = zero

DO j = jr_min, jr_max

!-----------------------------------------------------------------------
!  Compute and store interpolation coefficients
!-----------------------------------------------------------------------

  fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10(rho(j) )
  ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10(t(j)   )
  fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - yl(j) )

  id                 = idr(j,ij_ray,ik_ray)
  fdp(j)             = fd - DBLE( id )
  fdm(j)             = one - fdp(j)

  it                 = itr(j,ij_ray,ik_ray)
  ftp(j)             = ft - DBLE( it )
  ftm(j)             = one - ftp(j)

  iy                 = iyr(j,ij_ray,ik_ray)
  fyp(j)             = fy - DBLE( iy )
  fym(j)             = one - fyp(j)

END DO !  j = jr_min, jr_max

DO i = 1, 3
  DO j = jr_min, jr_max

!-----------------------------------------------------------------------
!  Store equation of state data for quantity i in temporary scalar
!   variables for interpolation.
!-----------------------------------------------------------------------

    av111            = estble(i,j,ida  ,ita  ,iya  ,ij_ray,ik_ray)
    av211            = estble(i,j,idap1,ita  ,iya  ,ij_ray,ik_ray)
    av121            = estble(i,j,ida  ,itap1,iya  ,ij_ray,ik_ray)
    av112            = estble(i,j,ida  ,ita  ,iyap1,ij_ray,ik_ray)
    av221            = estble(i,j,idap1,itap1,iya  ,ij_ray,ik_ray)
    av212            = estble(i,j,idap1,ita  ,iyap1,ij_ray,ik_ray)
    av122            = estble(i,j,ida  ,itap1,iyap1,ij_ray,ik_ray)
    av222            = estble(i,j,idap1,itap1,iyap1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Interpolate equation of state variable i (esv)
!-----------------------------------------------------------------------


    besv             =  fym(j) * ( fdm(j) * ( ftm(j) * av111 + ftp(j) * av121 )   &
&                     +            fdp(j) * ( ftm(j) * av211 + ftp(j) * av221 ) ) &
&                     + fyp(j) * ( fdm(j) * ( ftm(j) * av112 + ftp(j) * av122 )   &
&                     +            fdp(j) * ( ftm(j) * av212 + ftp(j) * av222 ) )
    esvtmp           = 10.d+00**besv
    aesv(j,i,ij_ray,ik_ray)  = esvtmp - escnst(i,j,ij_ray,ik_ray) - epsilon
    IF ( i >= 7 ) aesv(j,i,ij_ray,ik_ray) = DMAX1( aesv(j,i,ij_ray,ik_ray), zero )

!-----------------------------------------------------------------------
!  Store density derivative of equation of state data for quantity i
!   in temporary scalar variables for interpolation.
!-----------------------------------------------------------------------

    av111            = estbled(i,j,ida  ,ita  ,iya  ,ij_ray,ik_ray)
    av211            = estbled(i,j,idap1,ita  ,iya  ,ij_ray,ik_ray)
    av121            = estbled(i,j,ida  ,itap1,iya  ,ij_ray,ik_ray)
    av112            = estbled(i,j,ida  ,ita  ,iyap1,ij_ray,ik_ray)
    av221            = estbled(i,j,idap1,itap1,iya  ,ij_ray,ik_ray)
    av212            = estbled(i,j,idap1,ita  ,iyap1,ij_ray,ik_ray)
    av122            = estbled(i,j,ida  ,itap1,iyap1,ij_ray,ik_ray)
    av222            = estbled(i,j,idap1,itap1,iyap1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Interpolate equation of state variable d(esv(i))/d(density)
!-----------------------------------------------------------------------

    besvd            = fym(j) * ( fdm(j) * ( ftm(j) * av111 + ftp(j) * av121 )   &
&                    +            fdp(j) * ( ftm(j) * av211 + ftp(j) * av221 ) ) &
&                    + fyp(j) * ( fdm(j) * ( ftm(j) * av112 + ftp(j) * av122 )   &
&                    +            fdp(j) * ( ftm(j) * av212 + ftp(j) * av222 ) )
    esvtmpd          = 10.d0**besvd
    aesvd(j,i,ij_ray,ik_ray) = esvtmpd - escnstd(i,j,ij_ray,ik_ray) - epsilon

!-----------------------------------------------------------------------
!  Store temperature derivative of equation of state data for quantity i
!   in temporary scalar variables for interpolation.
!-----------------------------------------------------------------------

    av111            = estblet(i,j,ida  ,ita  ,iya  ,ij_ray,ik_ray)
    av211            = estblet(i,j,idap1,ita  ,iya  ,ij_ray,ik_ray)
    av121            = estblet(i,j,ida  ,itap1,iya  ,ij_ray,ik_ray)
    av112            = estblet(i,j,ida  ,ita  ,iyap1,ij_ray,ik_ray)
    av221            = estblet(i,j,idap1,itap1,iya  ,ij_ray,ik_ray)
    av212            = estblet(i,j,idap1,ita  ,iyap1,ij_ray,ik_ray)
    av122            = estblet(i,j,ida  ,itap1,iyap1,ij_ray,ik_ray)
    av222            = estblet(i,j,idap1,itap1,iyap1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Interpolate equation of state variable d(esv(i))/d(temperature)
!-----------------------------------------------------------------------

    besvt            = fym(j) * ( fdm(j) * ( ftm(j) * av111 + ftp(j) * av121 )   &
&                    +            fdp(j) * ( ftm(j) * av211 + ftp(j) * av221 ) ) &
&                    + fyp(j) * ( fdm(j) * ( ftm(j) * av112 + ftp(j) * av122 )   &
&                    +            fdp(j) * ( ftm(j) * av212 + ftp(j) * av222 ) )
    esvtmpt          = 10.d0**besvt
    aesvt(j,i,ij_ray,ik_ray) = esvtmpt - escnstt(i,j,ij_ray,ik_ray) - epsilon

!-----------------------------------------------------------------------
!  Store lepton fraction derivative of equation of state data for quantity i
!   in temporary scalar variables for interpolation.
!-----------------------------------------------------------------------

    av111            = estbley(i,j,ida  ,ita  ,iya  ,ij_ray,ik_ray)
    av211            = estbley(i,j,idap1,ita  ,iya  ,ij_ray,ik_ray)
    av121            = estbley(i,j,ida  ,itap1,iya  ,ij_ray,ik_ray)
    av112            = estbley(i,j,ida  ,ita  ,iyap1,ij_ray,ik_ray)
    av221            = estbley(i,j,idap1,itap1,iya  ,ij_ray,ik_ray)
    av212            = estbley(i,j,idap1,ita  ,iyap1,ij_ray,ik_ray)
    av122            = estbley(i,j,ida  ,itap1,iyap1,ij_ray,ik_ray)
    av222            = estbley(i,j,idap1,itap1,iyap1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Interpolate equation of state variable d(esv(i))/d(lepton fraction)
!-----------------------------------------------------------------------

    besvy            = fym(j) * ( fdm(j) * ( ftm(j) * av111 + ftp(j) * av121 )   &
&                    +            fdp(j) * ( ftm(j) * av211 + ftp(j) * av221 ) ) &
&                    + fyp(j) * ( fdm(j) * ( ftm(j) * av112 + ftp(j) * av122 )   &
&                    +            fdp(j) * ( ftm(j) * av212 + ftp(j) * av222 ) )
    esvtmpy          = 10.d0**besvy
    aesvy(j,i,ij_ray,ik_ray) = esvtmpy - escnsty(i,j,ij_ray,ik_ray) - epsilon

  END DO ! j = jr_min, jr_max
END DO ! i = 1, 12

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (fdp, STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fdp       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fdm, STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fdm       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ftp, STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'ftp       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ftm, STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'ftm       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fyp, STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fyp       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fym, STAT = istat)
 IF ( istat /= 0 ) THEN; var_name = 'fym       '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE eqstld