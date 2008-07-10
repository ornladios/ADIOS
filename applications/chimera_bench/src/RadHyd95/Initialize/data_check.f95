SUBROUTINE data_check( c_radhyd_data, i_radhyd_data, d_radhyd_data, &
& d_hydro_data, nx, ny, nz, nnu, n_proc )
!-----------------------------------------------------------------------
!
!    File:         data_check
!    Module:       data_check
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/04
!
!    Purpose:
!      Check data for consistency.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  c_radhyd_data   : Lagrangian - Eulerian toggle
!  i_radhyd_data   : integer array of radhyd keys
!  d_radhyd_data   : 64 bit real array of radhyd keys
!  d_hydro_data    : 64 bit real array of transport keys
!  nx              : x_array extent
!  ny              : y_array extent
!  nz              : z_array extent
!  nnu             : neutrino flavor array extent
!  n_proc          : number of processors assigned to the run
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE edit_module, ONLY : nprint, nlog

IMPLICIT none
SAVE


!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(in), DIMENSION(20) :: c_radhyd_data ! Lagrangian - Eulerian toggle

INTEGER, INTENT(in), DIMENSION(50)           :: i_radhyd_data ! integer array of edit keys

INTEGER, INTENT(in)                          :: nx            ! x-array extent
INTEGER, INTENT(in)                          :: ny            ! y-array extent
INTEGER, INTENT(in)                          :: nz            ! z-array extent
INTEGER, INTENT(in)                          :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)                          :: n_proc        ! number of processors assigned to the run

REAL(KIND=double), DIMENSION(50)             :: d_radhyd_data ! 64 bit real array of radhyd keys
REAL(KIND=double), DIMENSION((30+nx))        :: d_hydro_data  ! 64 bit real array of transport keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                                  :: lagrangian ! T: lagr = 'ye'; F: lagr = 'no'

CHARACTER (len=2)                        :: lagr       ! Lagrangian - Eulerian toggle
CHARACTER (len=2)                        :: regrid     ! re-grid switch
CHARACTER (len=2)                        :: G_trns     ! Galilean transformation switch

INTEGER                                  :: ndim       ! number of geometric dimensions
INTEGER                                  :: ngeomx     ! x-geometry flag
INTEGER                                  :: ngeomy     ! y-geometry flag
INTEGER                                  :: ngeomz     ! z-geometry flag
INTEGER                                  :: nleftx     ! lower x-boundary condition flag
INTEGER                                  :: nlefty     ! lower y-boundary condition flag
INTEGER                                  :: nleftz     ! lower z-boundary condition flag
INTEGER                                  :: nrightx    ! upper x-boundary condition flag
INTEGER                                  :: nrighty    ! upper y-boundary condition flag
INTEGER                                  :: nrightz    ! upper z-boundary condition flag
INTEGER                                  :: imin       ! minimum x-array index
INTEGER                                  :: imax       ! maximum x-array index
INTEGER                                  :: jmin       ! minimum y-array index
INTEGER                                  :: jmax       ! maximum y-array index
INTEGER                                  :: kmin       ! minimum z-array index
INTEGER                                  :: kmax       ! maximum z-array index

REAL(KIND=double)                        :: grid_frac  ! fraction of a grid width grid is allowed to move per time step
REAL(KIND=double)                        :: courant    ! fraction of a grid width sound is allowed to propagate per time step

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1007 FORMAT (' imax + 2=',i4,' > nx=',i4,' not enough radial zones to fit model')
 1009 FORMAT (' jmax=',i4,' /= ny=',i4,' jmax must equal ny in this version of the code')
 1011 FORMAT (' kmax=',i4,' /= nz=',i4,' kmax must equal nz in this version of the code')
 1013 FORMAT (' imax-imin+1=',i4,' is incomensurate with n_proc=',i4)
 1015 FORMAT (' jmax-jmin+1=',i4,' kmax-kmin+1=',i4, &
 & ' (jmax-jmin+1)(kmax-kmin+1) is incomensurate with n_proc=',i4)
 1017 FORMAT (' imax-imin+1=',i4,' is incomensurate with jmax-jmin+1=',i4)
 1019 FORMAT (' imax-imin+1=',i4,' is incomensurate with kmax-kmin+1=',i4)
 1031 FORMAT (' ndim not set')
 1033 FORMAT (' ngeomx not set')
 1035 FORMAT (' ngeomy not set')
 1037 FORMAT (' ngeomz not set')
 1039 FORMAT (' nleftx not set')
 1041 FORMAT (' nrightx not set')
 1043 FORMAT (' nlefty not set')
 1045 FORMAT (' nrighty not set')
 1047 FORMAT (' nleftz not set')
 1049 FORMAT (' nrightz not set')
 1061 FORMAT (' ndim =',i3,' and jmax=',i3,'> 1, 1D run with jmax > 1 not allowed')
 1063 FORMAT (' ndim =',i3,' and kmax=',i3,'> 1, 1D run with kmax > 1 not allowed')
 1065 FORMAT (' ndim =',i3,', jmax=', i3,', and kmax=',i3,'; MD run with jmax = kmax = 1 not allowed')
 1067 FORMAT (', jmax=', i3,', kmax=',i3,' and ndim =',i3,'; MD run with jmax > 1, kmax > 1, and ndim < 3 inconsistent')
 1069 FORMAT (' G_trns =',a2,', ndim=', i3,' need MD run to do Galilean transformations')
 1071 FORMAT (' grid_frac=',es11.3,' + courant=',es11.3,' > 0.9d0. Not allowed with re-grid option turned on.')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

lagr                = c_radhyd_data(1)
regrid              = c_radhyd_data(4)
G_trns              = c_radhyd_data(10)

ndim                = i_radhyd_data(3)
ngeomx              = i_radhyd_data(4)
ngeomy              = i_radhyd_data(5)
ngeomz              = i_radhyd_data(6)
nleftx              = i_radhyd_data(7)
nrightx             = i_radhyd_data(8)
nlefty              = i_radhyd_data(9)
nrighty             = i_radhyd_data(10)
nleftz              = i_radhyd_data(11)
nrightz             = i_radhyd_data(12)
imin                = i_radhyd_data(13)
imax                = i_radhyd_data(14)
jmin                = i_radhyd_data(15)
jmax                = i_radhyd_data(16)
kmin                = i_radhyd_data(17)
kmax                = i_radhyd_data(18)

grid_frac           = d_radhyd_data(19)
courant             = d_hydro_data(7)

!-----------------------------------------------------------------------
!
!               \\\\\ CHECK CONSISTENCY OF DATA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Test the compatability of the model extents and the array extents
!-----------------------------------------------------------------------

IF ( imax + 2 > nx ) THEN
  WRITE (nlog,1007) imax + 2, nx
  WRITE (nprint,1007) imax + 2, nx
  STOP
END IF ! imax + 2 > nx

IF ( jmax /= ny ) THEN
  WRITE (nlog,1009) jmax, ny
  WRITE (nprint,1009) jmax, ny
  STOP
END IF ! jmax /= ny

IF ( kmax /= nz ) THEN
  WRITE (nlog,1011) kmax, nz
  WRITE (nprint,1011) kmax, nz
  STOP
END IF ! kmax /= nz

!-----------------------------------------------------------------------
!  Test for commensorability of array extents with tne number of
!   processors
!-----------------------------------------------------------------------

IF ( imax-imin+1 >= n_proc ) THEN
  IF ( MOD( imax-imin+1, n_proc ) /= 0 ) THEN
    WRITE (nlog,1013) imax-imin+1, n_proc
    WRITE (nprint,1013) imax-imin+1, n_proc
    STOP
  END IF ! MOD( imax-imin+1, n_proc ) /= 0
ELSE ! imax-imin+1 < n_proc
  IF ( MOD( n_proc, imax-imin+1 ) /= 0 ) THEN
    WRITE (nlog,1013) imax-imin+1, n_proc
    WRITE (nprint,1013) imax-imin+1, n_proc
    STOP
  END IF ! MOD( imax-imin+1, n_proc ) /= 0
END IF ! imax-imin+1 >= n_proc

IF ( MOD( ( jmax-jmin+1 ) * ( kmax-kmin+1 ), n_proc ) /= 0 ) THEN
  WRITE (nlog,1015) jmax-jmin+1, n_proc
  WRITE (nprint,1015) jmax-jmin+1, kmax-kmin+1, n_proc
  STOP
END IF ! MOD( imax-imin+1, n_proc ) /= 0

!-----------------------------------------------------------------------
!  Test for commensorability of array extents with each other
!-----------------------------------------------------------------------

IF ( MOD( imax-imin+1, jmax-jmin+1 ) /= 0 ) THEN
  WRITE (nlog,1017) imax-imin+1, jmax-jmin+1
  WRITE (nprint,1017) imax-imin+1, jmax-jmin+1
  STOP
END IF ! MOD( imax-imin+1, jmax-jmin+1 ) /= 0

IF ( MOD( imax-imin+1, kmax-kmin+1 ) /= 0 ) THEN
  WRITE (nlog,1019) imax-imin+1, kmax-kmin+1
  WRITE (nprint,1019) imax-imin+1, kmax-kmin+1
  STOP
END IF ! MOD( imax-imin+1, jmax-jmin+1 ) /= 0

!-----------------------------------------------------------------------
!  Test if this is a Lagrangian run, which requires 1D
!-----------------------------------------------------------------------

IF (lagr/='no') THEN
  IF ( jmax - jmin > 1 .or. kmax - kmin > 1  .or.  ndim /= 1 ) THEN
    WRITE (nlog,*) 'Multi-D Lagrangian Runs impossible',jmin,jmax,kmin,kmax,ndim
    STOP
  ELSE
    lagrangian       =.true.
  END IF
ELSE
  lagrangian         =.false.
END IF
WRITE (nlog,*) 'lagrangian=',lagrangian

!-----------------------------------------------------------------------
!  Test for the setting of geometry and boundary conditions
!-----------------------------------------------------------------------

IF ( ndim == -1 ) THEN
  WRITE (nlog,1031) 
  STOP
END IF

IF ( ngeomx == -1 ) THEN
  WRITE (nlog,1033) 
  STOP
END IF

IF ( ndim <= 1  .and.  ngeomy == -1 ) THEN
  WRITE (nlog,1035) 
  STOP
END IF

IF ( ndim <= 2  .and.  ngeomz == -1 ) THEN
  WRITE (nlog,1037) 
  STOP
END IF

IF ( nleftx == -1 ) THEN
  WRITE (nlog,1039) 
  STOP
END IF

IF ( nrightx == -1 ) THEN
  WRITE (nlog,1041) 
  STOP
END IF

IF ( ndim <= 1  .and.  nlefty == -1 ) THEN
  WRITE (nlog,1043) 
  STOP
END IF

IF ( ndim <= 1  .and.  nrighty == -1 ) THEN
  WRITE (nlog,1045) 
  STOP
END IF

IF ( ndim <= 2  .and.  nleftz == -1 ) THEN
  WRITE (nlog,1047) 
  STOP
END IF

IF ( ndim <= 2  .and.  nrightz == -1 ) THEN
  WRITE (nlog,1049) 
  STOP
END IF

!-----------------------------------------------------------------------
!  Test for consistency of dimensionality and array extents
!-----------------------------------------------------------------------

IF ( ndim == 1  .and.  jmax > 1 ) THEN
  WRITE (nlog,1061) ndim, jmax
  STOP
END IF

IF ( ndim == 1  .and.  kmax > 1 ) THEN
  WRITE (nlog,1063) ndim, kmax
  STOP
END IF

IF ( ndim > 1  .and.  jmax == 1  .and.  kmax == 1 ) THEN
  WRITE (nlog,1065) ndim, jmax, kmax
  STOP
END IF

IF ( jmax > 1  .and.  kmax > 1  .and.  ndim < 3 ) THEN
  WRITE (nlog,1067) jmax, kmax, ndim
  STOP
END IF

!-----------------------------------------------------------------------
!  Test for consistency of Galilean transformations and dimensionality
!-----------------------------------------------------------------------

IF ( G_trns == 'ye'  .and.  ndim == 1 ) THEN
  WRITE (nlog,1069) G_trns, ndim
  STOP
END IF

!-----------------------------------------------------------------------
!  Test for adequate grid movement restrictions
!-----------------------------------------------------------------------

IF ( regrid == 'ye' ) THEN
  IF ( grid_frac + courant > 0.9d0 ) THEN
    WRITE (nlog,1071) grid_frac, courant
    STOP
  END IF ! grid_frac + courant > 0.9d0
END IF ! regrid == 'ye'

RETURN
END SUBROUTINE data_check
