SUBROUTINE read_pack_initial_model( nreadp, nprint, iskip, nx, nez, nnu, &
& ij_ray_dim, ik_ray_dim, i_model_data, d_model_data1, d_model_data2,    &
& d_model_data3, d_psi_data2, nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_initial_model
!    Module:       read_pack_initial_model
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the initial model configuration and to pack them into an
!       integer and a real*8 array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp        : unit number from which to read
!  nprint        : unit number from which to print
!  iskip         : echo data read flag
!  nx            : radial array dimension
!  nez           : neutrino energy array extent
!  ij_ray_dim    : number of y-zones on a processor before swapping
!  ik_ray_dim    : number of z-zones on a processor before swapping
!  nnu           : neutrino flavor array extent
!  nrst          : cycle number at start or restart
!
!    Output arguments:
!  i_model_data  : integer array of initial model data
!  d_model_data1 : 64 bit real array of initial model data
!  d_model_data2 : 64 bit real array of initial model data
!  d_model_data3 : 64 bit real array of initial model data
!  d_psi_data2   : 64 bit real array of initial model data
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : dnurad, unukrad, unujrad, nnukrad, nnujrad, &
& e_rad, unurad, elec_rad, nnurad

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(2)                         :: i_model_data  ! integer array of initial model data


REAL(KIND=double), INTENT(out), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)         :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx)                               :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx+1)                             :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: d_psi_data2   ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line
CHARACTER (len=10)               :: var_name

INTEGER                          :: i             ! do index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

INTEGER                          :: k1            ! integer data variable to read in an interger datum
INTEGER                          :: k2            ! integer data variable to read in an interger datum
INTEGER                          :: k3            ! integer data variable to read in an interger datum
INTEGER                          :: k4            ! integer data variable to read in an interger datum
INTEGER                          :: k5            ! integer data variable to read in an interger datum

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl1           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl2           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl3           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl4           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl5           ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!  Radial mass zoning
!-----------------------------------------------------------------------
!  jr_max  : the innermost shifted radial zone for thvar data
!
!  jr_max  : the outermost shifted radial zone for thvar data
!
!  jrx_max : the innermost shifted radial zone for stvar data
!
!  jrx_max : the outermost shifted radial zone for stvar data
!
!  imin    : the innermost unshifted radial zone
!
!  imax    : the outermost unshifted radial zone
!-----------------------------------------------------------------------

INTEGER                                          :: jr_min
INTEGER                                          :: jr_max
INTEGER                                          :: jrs_min
INTEGER                                          :: jrs_max
INTEGER                                          :: imin
INTEGER                                          :: imax

!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------
!  rho(j) : the density of radial zone j at timestep m (current time)
!   (g/cm**3).
!
!  t(j)   : the temperature of radial zone j at timestep m (K).
!
!  ye(j)  : the electron fraction of radial zone j at timestep m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye

!-----------------------------------------------------------------------
!  Mechanical state
!-----------------------------------------------------------------------
!  u(j)  : the zone averaged x-velocity of radial zone j at timestep m
!   (cm/s).
!
!  v(j)  : the zone averaged y-velocity of radial zone j at timestep m
!   (cm/s).
!
!  e(j)  : the zone averaged z-velocity of radial zone j at timestep m
!   (cm/s).
!
!  r(j)  : the radius of outer boundary of radial zone j at timestep m
!   (cm).
!
!  dr(j) : r(j) - r(j-1).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: u
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: v
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: w
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dr

!-----------------------------------------------------------------------
!  Specific energy density arrays for rezonin
!-----------------------------------------------------------------------
!  uad(j) : the specific energy added to radial zone j to set the GR
!   gravity equal to the Newtonian gravity at problem initiation (ergs/g).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: uad

!-----------------------------------------------------------------------
!  Radiation Variables
!-----------------------------------------------------------------------
!  psi0(j,k,n) : the zero moment of the occupation distribution for
!   neutrinos at the midpoint of radial zone j, of energy zone k, and of
!   type n.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi0

    1 FORMAT (a128)
    2 FORMAT (1x,a128)
  101 FORMAT (a128)
  103 FORMAT (a6)
  105 FORMAT (1x,a128)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)
  115 FORMAT (28x,2a)
  131 FORMAT (20x,i10,5x,e15.8)
  133 FORMAT (1x,a6,14x,i10,5x,es15.8,22x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  161 FORMAT (10x,2i10,5x,e15.8)
  163 FORMAT (1x,a6,4x,2i10,5x,es15.8,22x,a16)
  211 FORMAT (10x,i10,5x,e15.8)
  213 FORMAT (1x,a6,4x,i10,5x,es15.8)
  231 FORMAT (10x,i10,3(5x,e15.8))
  233 FORMAT (1x,a6,4x,i10,3(5x,es15.8))
  241 FORMAT (10x,i10,5(5x,e15.8))
  243 FORMAT (1x,a6,4x,i10,5(5x,es15.8))
  311 FORMAT (10x,2i5,5(i3,1x,e15.8))
  313 FORMAT (1x,a6,4x,2i5,5(i3,1x,es15.8))
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)
 1001 FORMAT (' Allocation problem for array ',a10,' in read_pack_initial_model')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_pack_initial_model')
 3001 FORMAT (' jr_min (thvar data) =',i4,' /= jrs_min (stvar data) =',i4)
 3003 FORMAT (' jr_max (thvar data) =',i4,' /= jrs_max (stvar data) =',i4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rho(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uad(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uad       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

jr_min                     = 1000000
jr_max                     = 0
jrs_min                    = 1000000
jrs_max                    = 0
imin                       = 0
imax                       = 0

rho                        = zero
t                          = zero
ye                         = zero
u                          = zero
r                          = zero
dr                         = zero
uad                        = zero
psi0                       = zero

!-----------------------------------------------------------------------
!
!                \\\\\ READ INITIAL MODEL DATA /////
!
!-----------------------------------------------------------------------

REWIND (nreadp)

READ: DO

!-----------------------------------------------------------------------
!  Read line
!-----------------------------------------------------------------------

  READ (nreadp,101,END=5000) line

  type                     = line(1:6)
  name                     = line(73:88)

!-----------------------------------------------------------------------
!  Ignore comment cards in the data stream
!-----------------------------------------------------------------------

  IF ( type == 'cccccc'  .or.  type(1:1) == '!'  .or.  type == '      ') THEN
    IF ( iskip == 0 ) WRITE (nprint,105) line
    CYCLE
  END IF ! type = 'cccccc' or type = '!'

!-----------------------------------------------------------------------
!  Model configuration arrays
!-----------------------------------------------------------------------

!........thvar

  IF ( type == 'thvar ' ) THEN
    READ (line ,231) j,rl1,rl2,rl3
    rho(j)                 = rl1
    t(j)                   = rl2
    ye(j)                  = rl3
    IF ( j > jr_max ) jr_max = j
    IF ( j < jr_min ) jr_min = j
    IF ( iskip == 0 ) WRITE (nprint,233) type,j,rho(j),t(j),ye(j)
    CYCLE
  END IF ! type = 'thvar '

!........stvar

  IF ( type == 'stvar ' ) THEN
    READ (line ,241) j,rl1,rl2,rl3,rl4,rl5
    u(j)                   = rl1
    v(j)                   = rl2
    w(j)                   = rl3
    dr(j)                  = rl4
    r(j)                   = rl5
    IF ( j > jrs_max ) jrs_max = j
    IF ( j < jrs_min ) jrs_min = j
    IF ( iskip == 0 ) WRITE (nprint,243) type,j,u(j),v(j),w(j),dr(j),r(j)
    CYCLE
  END IF ! type = 'stvar '

!........uad

  IF ( type == 'uad   ' ) THEN
    READ (line ,211) j,rl
    uad(j)                 = rl
    IF ( iskip == 0 ) WRITE (nprint,213) type,j,uad(j)
    CYCLE
  END IF ! type = 'uad   '

!........psi0

  IF ( type == 'psi0  ' ) THEN
    READ (line ,311) n,j,k1,rl1,k2,rl2,k3,rl3,k4,rl4,k5,rl5
    psi0(j,k1,n)           = rl1
    psi0(j,k2,n)           = rl2
    psi0(j,k3,n)           = rl3
    psi0(j,k4,n)           = rl4
    psi0(j,k5,n)           = rl5
    IF ( iskip == 0 ) WRITE (nprint,313) type,n,j,k1,psi0(j,k1,n),k2,psi0(j,k2,n),k3, &
&    psi0(j,k3,n),k4,psi0(j,k4,n),k5,psi0(j,k5,n)
    CYCLE
  END IF ! type = 'psi0  '

!........unucr

  IF ( type == 'unucr ' ) THEN

    IF ( name == 'e_rad   ') THEN
      READ (line ,141) rl
      e_rad(:,:)           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,e_rad(1,1),name
      CYCLE
    END IF ! name = 'e_rad   '

    IF ( name == 'unurad  ') THEN
      READ (line ,131) n,rl
      unurad(n,:,:)        = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,unurad(n,1,1),name
      CYCLE
    END IF ! name = 'unurad  '

    IF ( name == 'unukrad ') THEN
      READ (line ,161) n,k,rl
      unukrad(k,n,:,:)     = rl
      IF ( iskip == 0 ) WRITE (nprint,163) type,n,k,unukrad(k,n,1,1),name
      CYCLE
    END IF ! name = 'unukrad '

    IF ( name == 'unujrad ') THEN
      READ (line ,161) n,j,rl
      unujrad(j,n,:,:)     = rl
      IF ( iskip == 0 ) WRITE (nprint,163) type,n,j,unujrad(j,n,1,1),name
      CYCLE
    END IF ! name = 'unujrad '

  END IF ! type = 'unucr '

!........nnucr

  IF ( type == 'nnucr ' ) THEN

    IF ( name == 'elec_rad') THEN
      READ (line ,141) rl
      elec_rad(:,:)        = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,elec_rad(1,1),name
      CYCLE
    END IF ! name = 'elec_rad'

    IF ( name == 'nnukrad ') THEN
      READ (line ,161) n,k,rl
      nnukrad(k,n,:,:)     = rl
      IF ( iskip == 0 ) WRITE (nprint,163) type,n,k,nnukrad(k,n,1,1),name
      CYCLE
    END IF ! name = 'nnukrad '

    IF ( name == 'nnujrad ') THEN
      READ (line ,161) n,j,rl
      nnujrad(j,n,:,:)     = rl
      IF ( iskip == 0 ) WRITE (nprint,163) type,n,j,nnujrad(j,n,1,1),name
      CYCLE
    END IF ! name = 'nnujrad '

    IF ( name == 'nnurad  ') THEN
      READ (line ,131) n,rl
      nnurad(n,:,:)        = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,nnurad(n,1,1),name
      CYCLE
    END IF ! name = 'nnurad  '

  END IF ! type = 'nnucr '

!-----------------------------------------------------------------------
!  Unrecognized card
!-----------------------------------------------------------------------

  IF ( nrst == 0 ) THEN
    WRITE (nprint,401)
    WRITE (nprint,403) line
  END IF

END DO READ

 5000 CONTINUE

!-----------------------------------------------------------------------
!
!             \\\\\ CHECK DATA EXTENTS FOR CONSISTENCY /////
!
!  jr_min and jrs_min must be the same; jr_max and jrs_max must be the same
!-----------------------------------------------------------------------

IF ( jr_min /= jrs_min ) THEN
  WRITE (nlog,3001) jr_min, jrs_min
  WRITE (nprint,3001) jr_min, jrs_min
  STOP
END IF ! jr_min /= jrs_min

IF ( jr_max /= jrs_max ) THEN
  WRITE (nlog,3003) jr_max, jrs_max
  WRITE (nprint,3003) jr_max, jrs_max
  STOP
END IF ! jr_min /= jrs_min

!-----------------------------------------------------------------------
!
!                \\\\\ PACK INITIAL MODEL DATA /////
!
!-----------------------------------------------------------------------

DO i = 1,nx-1
  r(i+1)                   = r(i) + dr(i+1)
END DO ! i = 1,nx-1

imin                       = jr_min - 1
imax                       = jr_max - 1

i_model_data(1)            = imin
i_model_data(2)            = imax

DO i = 1,nx-1
  d_model_data1(1,i,:,:)   = rho(i+1)
  d_model_data1(2,i,:,:)   = t(i+1)
  d_model_data1(3,i,:,:)   = ye(i+1)
  d_model_data1(4,i,:,:)   = u(i+1)
  d_model_data1(5,i,:,:)   = v(i+1)
  d_model_data1(6,i,:,:)   = w(i+1)
  d_model_data1(7,i,:,:)   = zero
  d_model_data2(1,i)       = dr(i+1)
  d_model_data3(1,i)       = r(i)
END DO ! i

DO i = 1,nx-1
  DO k = 1,nez
    DO n = 1,nnu
      d_psi_data2(1,i,k,n,:,:) = psi0(i+1,k,n)
    END DO ! n = 1,nnu
  END DO ! k = 1,nez
END DO ! i = 1,nx-1

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (rho, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (v, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (w, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (r, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (uad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uad       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE read_pack_initial_model
