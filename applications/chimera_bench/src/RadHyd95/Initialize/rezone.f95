SUBROUTINE rezone( c_radhyd_data, i_radhyd_data, d_radhyd_data, &
& l_rezone_data, d_rezone_data )
!-----------------------------------------------------------------------
!
!    File:         rezone
!    Module:       rezone
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/13/04
!
!    Purpose:
!      To generate a y (angular and a z (azimuthal) grid and to laod
!       the rezone_array buffers.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  grid           : generates a y (angular) and a z (azimuthal) grid
!
!    Input-output arguments:
!  c_radhyd_data  : character array of radhyd keys
!  i_radhyd_data  : integer array of radhyd keys
!  d_radhyd_data  : 64 bit real array of radhyd keys
!  l_rezone_data  : logical array of zone and rezone data
!  d_rezone_data  : 64 bit real array of zone and rezone data
!
!    Output arguments:
!        none
!
!    Output arguments (common):
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: nx, ny, nz, nez, nnu, nnc, ij_ray_dim, ik_ray_dim
USE numerical_module, ONLY: zero, half
USE physcnst_module, ONLY: pi

USE edit_module, ONLY : nlog


IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(inout), DIMENSION(20) :: c_radhyd_data ! character array of radhyd keys

LOGICAL, INTENT(inout), DIMENSION(1)            :: l_rezone_data ! lagrangian flag

INTEGER, INTENT(inout), DIMENSION(50)           :: i_radhyd_data ! integer array of radhyd keys

REAL(KIND=double), INTENT(inout), DIMENSION(50)                  :: d_radhyd_data ! 64 bit real array of radhyd keys
REAL(KIND=double), INTENT(inout), DIMENSION(20+3*ny+3*nz+2)      :: d_rezone_data ! 64 bit real array of rezone data

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=2)                              :: lagr       ! Lagrangian - Eulerian toggle
CHARACTER (len=2)                              :: rezn       !  Regridding toggle.

LOGICAL                                        :: lagrangian ! lagrangian flag

INTEGER                                        :: i          ! do index
INTEGER                                        :: ngeomy     ! y-geometry flag
INTEGER                                        :: ngeomz     ! z-geometry flag

INTEGER                                        :: ndim       ! problem dimension
INTEGER                                        :: jmin       ! minimum y-array index
INTEGER                                        :: jmax       ! maximum y-array index
INTEGER                                        :: kmin       ! minimum z-array index
INTEGER                                        :: kmax       ! maximum z-array index


REAL(KIND=double), DIMENSION(ny+1)             :: y_ei       ! y edge coordinate
REAL(KIND=double), DIMENSION(ny)               :: y_ci       ! y center coordinate
REAL(KIND=double), DIMENSION(ny)               :: dy_ci      ! y_ei(j+1) - y_ei(j)
REAL(KIND=double), DIMENSION(nz+1)             :: z_ei       ! z edge coordinate
REAL(KIND=double), DIMENSION(nz)               :: z_ci       ! z center coordinate
REAL(KIND=double), DIMENSION(nz)               :: dz_ci      ! z_ei(j+1) - z_ei(j)
REAL(KIND=double)                              :: xmin       ! minimum x
REAL(KIND=double)                              :: xmax       ! maximum x
REAL(KIND=double)                              :: ymin       ! minimum y
REAL(KIND=double)                              :: ymax       ! maximum y
REAL(KIND=double)                              :: zmin       ! minimum z
REAL(KIND=double)                              :: zmax       ! maximum z
REAL(KIND=double)                              :: xmin_i     ! minimum x before modified in rezone
REAL(KIND=double)                              :: xmax_i     ! maximum x before modified in rezone
REAL(KIND=double)                              :: ymin_i     ! minimum y before modified in rezone
REAL(KIND=double)                              :: ymax_i     ! maximum y before modified in rezone
REAL(KIND=double)                              :: zmin_i     ! minimum z before modified in rezone
REAL(KIND=double)                              :: zmax_i     ! maximum z before modified in rezone
REAL(KIND=double)                              :: zoom       ! grid zoom factor
REAL(KIND=double)                              :: courant    ! maximum courant condition

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!         \\\\\ UNPACK DATA NEEDED FOR SETTING UP GRID /////
!
!-----------------------------------------------------------------------

lagr                 = c_radhyd_data(1)
rezn                 = c_radhyd_data(2)

ndim                 = i_radhyd_data(3)
ngeomy               = i_radhyd_data(5)
ngeomz               = i_radhyd_data(6)
jmin                 = i_radhyd_data(15)
jmax                 = i_radhyd_data(16)
kmin                 = i_radhyd_data(17)
kmax                 = i_radhyd_data(18)

!-----------------------------------------------------------------------
!  Set EVH1 globals
!-----------------------------------------------------------------------

courant              = 0.60d0
IF ( lagr == 'ye' ) lagrangian = .true.
IF ( lagr == 'no' ) lagrangian = .false.

!-----------------------------------------------------------------------
!  Get radial limits from MGFLD
!-----------------------------------------------------------------------

xmin                 = d_radhyd_data(24)
xmax                 = d_radhyd_data(25)
ymin                 = d_radhyd_data(26)
ymax                 = d_radhyd_data(27)
zmin                 = d_radhyd_data(28)
zmax                 = d_radhyd_data(29)

xmin_i               = d_radhyd_data(24)
xmax_i               = d_radhyd_data(25)
ymin_i               = d_radhyd_data(26)
ymax_i               = d_radhyd_data(27)
zmin_i               = d_radhyd_data(28)
zmax_i               = d_radhyd_data(29)
   
!-----------------------------------------------------------------------
!  Build j grid
!-----------------------------------------------------------------------

y_ei                 = zero
dy_ci                = zero
y_ci                 = zero

IF ( ndim == 1 ) THEN
  IF ( ndim == 1  .and.  ngeomy >= 3 ) THEN         ! For angular dimensions,
    ymin             = ymin * pi  ! convert to radians
    ymax             = ymax * pi
  END IF ! ngeomy >= 3
  y_ei(1)            = ymin
  y_ei(2)            = ymax
  dy_ci(1)           = ymax - ymin
  y_ci(1)            = half * ( ymin + ymax )
END IF ! ndim == 1

zoom                 = 1.0d0
IF ( ndim >= 2 ) THEN
  IF ( ngeomy >= 3 ) THEN         ! For angular dimensions,
    ymin             = ymin * pi  ! convert to radians
    ymax             = ymax * pi
  END IF ! geomy >= 3
  CALL grid( jmin, jmax, ymin, ymax, y_ei, y_ci, dy_ci, zoom, ny )
END IF ! ndim >= 2
 
!-----------------------------------------------------------------------
!  Build k grid
!-----------------------------------------------------------------------

z_ei                 = zero
dz_ci                = zero
z_ci                 = zero

zoom                 = 1.0d0
IF ( ndim == 3 ) THEN
  IF ( ngeomz >= 3 ) THEN         ! For angular dimensions,
    zmin             = zmin * pi  ! convert to radians
    zmax             = zmax * pi
  END IF ! ngeomz == 3
  CALL grid( kmin, kmax, zmin, zmax, z_ei, z_ci, dz_ci, zoom, nz )
END IF ! ndim == 3
      
!-----------------------------------------------------------------------
!
!                       \\\\\ PACK DATA /////
!
!-----------------------------------------------------------------------

l_rezone_data(1)     = lagrangian

d_rezone_data(1)     = xmin
d_rezone_data(2)     = xmax
d_rezone_data(3)     = ymin
d_rezone_data(4)     = ymax
d_rezone_data(5)     = zmin
d_rezone_data(6)     = zmax
d_rezone_data(7)     = xmin_i
d_rezone_data(8)     = xmax_i
d_rezone_data(9)     = ymin_i
d_rezone_data(10)    = ymax_i
d_rezone_data(11)    = zmin_i
d_rezone_data(12)    = zmax_i
d_rezone_data(13)    = courant

DO i = 1,ny
  d_rezone_data(20+i)             = y_ei(i)
  d_rezone_data(20+ny+1+i)        = y_ci(i)
  d_rezone_data(20+2*ny+1+i)      = dy_ci(i)
END DO
d_rezone_data(20+ny+1)            = y_ei(ny+1)

DO i = 1,nz
  d_rezone_data(20+3*ny+1+i)      = z_ei(i)
  d_rezone_data(20+3*ny+nz+2+i)   = z_ci(i)
  d_rezone_data(20+3*ny+2*nz+2+i) = dz_ci(i)
END DO
d_rezone_data(20+3*ny+nz+2)       = z_ei(nz+1)

RETURN
END SUBROUTINE rezone 
