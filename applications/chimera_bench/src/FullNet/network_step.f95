SUBROUTINE network_step( dtime_start, dtime_stop, rho, temp, ye, &
& xn, nuc_min, nuc_max, nx, nncext, x_active, dtime )
!----------------------------------------------------------------------
!  This routine evolves the abundances from dtime_start to dtime_stop
!----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnc
USE numerical_module, ONLY : zero, one

USE abundances, ONLY : y
USE nuclear_data, ONLY : aa
USE nuc_number, ONLY : ny
USE thermo_data, ONLY : nh, tstart, tstop, t9h, rhoh, th

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                     :: nx          ! radial array extent
INTEGER, INTENT(in)                     :: nncext      ! abundance array extent
INTEGER, INTENT(in)                     :: nuc_min     ! shifted index of first zone to compute nuclear burn
INTEGER, INTENT(in)                     :: nuc_max     ! shifted index of last zone to compute nuclear burn

REAL(KIND=double), INTENT(in)           :: dtime_start ! time at beginning of time step
REAL(KIND=double), INTENT(in)           :: dtime_stop  ! time at end of time step
REAL(KIND=double), INTENT(in)           :: dtime       ! time step of cycle or subcycle

REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: rho         ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: temp        ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: ye          ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx)        :: x_active    ! mass fraction of nuclei participating in reactions

!-----------------------------------------------------------------------
!        Input-output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,nncext) :: xn       ! abundance mass fractions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                 :: i, k, n     ! loop variables

REAL(KIND=double)                       :: mfnorm      ! mass fraction normalization parameter
REAL(KIND=double)                       :: mfsum       ! sum of mass fractions

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

nh                = 2

!-----------------------------------------------------------------------
!  Set starting and stopping time for network
!-----------------------------------------------------------------------

tstart            = dtime_start
tstop             = dtime_stop

!-----------------------------------------------------------------------
!  Loop over all zones
!-----------------------------------------------------------------------

n                 = nuc_min - 1

DO i = nuc_min,nuc_max

!-----------------------------------------------------------------------
!  Normalize mass fractions, convert to abundances
!-----------------------------------------------------------------------

  mfnorm          = one/x_active(i)
  y(1:ny)         = mfnorm * xn(i,1:ny)/aa(1:ny)

!-----------------------------------------------------------------------
!  Specify thermo conditions for each zone
!-----------------------------------------------------------------------

  t9h(1:2)        = temp(i)/1.0d9
  rhoh(1:2)       = rho(i)
  th(1:2)         = (/tstart,tstop/)

!-----------------------------------------------------------------------
!  Call network
!-----------------------------------------------------------------------

  n               = n + 1
  CALL full_net( n, dtime )

!-----------------------------------------------------------------------
!  Renormalize mass fractions
!-----------------------------------------------------------------------

  mfsum           = SUM(y*aa)
  mfnorm          = x_active(i)/mfsum

  DO k = 1,ny
    y(k)          = mfnorm * y(k)
  END DO

!-----------------------------------------------------------------------
!  Update internal energy and mass fractions
!-----------------------------------------------------------------------

  DO k = 1,ny
    xn(i,k)       = aa(k) * y(k)
  END DO ! k
               
END DO ! i

RETURN
END
