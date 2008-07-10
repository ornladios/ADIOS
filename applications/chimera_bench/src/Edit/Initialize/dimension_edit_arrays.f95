SUBROUTINE dimension_edit_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         dimension_edit_arrays
!    Module:       dimension_edit_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate the dimensions to certain of the evh1 arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array extent
!  ny         : y (angular) array extent
!  nz         : z (azimuthal) array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x (radial) array extent
INTEGER, INTENT(in)              :: ny            ! y (angular) array extent
INTEGER, INTENT(in)              :: nz            ! z (azimuthal) array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' edit_arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_edit_arrays')

!-----------------------------------------------------------------------
!
!              \\\\\ ALLOCATE EDIT ARRAY DIMENSIONS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  General edit parameters
!-----------------------------------------------------------------------

ALLOCATE (nprt(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nprt      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nmodel(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nmodel    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Parameters for diferential neutrino edits
!-----------------------------------------------------------------------

ALLOCATE (intedn(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'intedn    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nedn(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedn      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (idxedn(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idxedn    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Parameters for integral neutrino edits
!-----------------------------------------------------------------------

ALLOCATE (intdng(40,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'intdng    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nedng(40,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedng     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (idxeng(40,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idxeng    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Parameters for nurad plot edits
!-----------------------------------------------------------------------

ALLOCATE (nu_r(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_r      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_rt(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rt     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_rho(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rho    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_rhot(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rhot   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unu_r(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unu_r     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unu_rt(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unu_rt    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unu_rho(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unu_rho   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unu_rhot(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unu_rhot  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Parameters for neutrino distribution plot edits
!-----------------------------------------------------------------------

ALLOCATE (psi0dat(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0dat   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1dat(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1dat   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Work
!-----------------------------------------------------------------------

ALLOCATE (pdv(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pdv       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (twrk(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'twrk      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Stresses
!-----------------------------------------------------------------------

ALLOCATE (nustrss(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nustrss   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pstrss(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pstrss    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gstrss(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gstrss    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gstrss_cx(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gstrss_cx '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gstrss_cy(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gstrss_cy '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rstrss(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rstrss    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Gravotational potential
!-----------------------------------------------------------------------

ALLOCATE (g_pot(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'g_pot     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Solid angles subtended by rays
!-----------------------------------------------------------------------

ALLOCATE (d_omega(ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_omega   '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Angularly averaged quantities
!-----------------------------------------------------------------------

ALLOCATE (rhobar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhobar    '; WRITE  (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-matter energy transfer rates
!-----------------------------------------------------------------------

ALLOCATE (dudt_ABEM(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_ABEM '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dudt_NES(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_NES  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dudt_NNS(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_NNS  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dudt_NAS(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_NAS  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dudt_PR(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_PR   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dudt_Brem(nx,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_Brem '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dudt_NET(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dudt_NET  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE EDIT PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  General edit parameters
!-----------------------------------------------------------------------

nread                     = 0
nrrst                     = 0
iprint                    = 0
nprint                    = 0
intprt                    = 0
nprt                      = 0
iprtbgn                   = 0
iflprt                    = 0
nmodel                    = 0
iplot                     = 0
nplot                     = 0
intplf                    = 0
npltf                     = 0
nplotc                    = 0
nplote                    = 0
nplota                    = 0
nplott                    = 0
i_edit                    = 0

rhoprint                  = -1.d+00

!-----------------------------------------------------------------------
!  Parameters for configuration edits
!-----------------------------------------------------------------------

intedc                    = 9000000
nedc                      = 0
idxedc                    = 1

!-----------------------------------------------------------------------
!  Parameters for kinetic, internal, and gravitational energy edits
!-----------------------------------------------------------------------

intede                    = 9000000
nede                      = 0
idxede                    = 1

!-----------------------------------------------------------------------
!  Parameters for editmi edits
!-----------------------------------------------------------------------

intdmi                    = 9000000
nedmi                     = 0
idxemi                    = 1

!-----------------------------------------------------------------------
!  Parameters for mass average edits
!-----------------------------------------------------------------------

intdma                    = 9000000
nedma                     = 0
idxema                    = 1

!-----------------------------------------------------------------------
!  Parameters for hydro edits
!-----------------------------------------------------------------------

intedh                    = 9000000
nedh                      = 0
idxedh                    = 1

!-----------------------------------------------------------------------
!  Parameters for pressure - stress edits
!-----------------------------------------------------------------------

intdps                    = 9000000
nedps                     = 0
idxeps                    = 1

!-----------------------------------------------------------------------
!  Parameters for energy edits
!-----------------------------------------------------------------------

intedu                    = 9000000
nedu                      = 0
idxedu                    = 1

!-----------------------------------------------------------------------
!  Parameters for composition edits
!-----------------------------------------------------------------------

intedy                    = 9000000
nedy                      = 0
idxedy                    = 1

!-----------------------------------------------------------------------
!  Parameters for entropy and chemical potential edits
!-----------------------------------------------------------------------

intdsc                    = 9000000
nedsc                     = 0
idxesc                    = 1

!-----------------------------------------------------------------------
!  Parameters for diferential neutrino edits
!-----------------------------------------------------------------------

intedn                    = 9000000
nedn                      = 0
idxedn                    = 1
niedn                     = 9000000
neden                     = 0

!-----------------------------------------------------------------------
!  Parameters for integral neutrino edits
!-----------------------------------------------------------------------

intdng                    = 9000000
nedng                     = 0
idxeng                    = 1

!-----------------------------------------------------------------------
!  Parameters for integral neutrino edits
!-----------------------------------------------------------------------

ilumplt                   = 0
nlumplt1                  = 0
nlumplt2                  = 0
nlum                      = 0
intlum                    = 0
ncylum                    = 0

r_lum                      = 1.d+08
d_lum                    = 1.d+12

!-----------------------------------------------------------------------
!  Parameters for rms neutrino energy plot edits
!-----------------------------------------------------------------------

ienuplt                   = 0
nenuplt1                  = 0
nenuplt2                  = 0
nenu                      = 0
intenu                    = 0
ncyenu                    = 0

r_e_rms                      = 1.d+08
d_e_rms                    = 1.d+12

!-----------------------------------------------------------------------
!  Parameters for selected radii plot edits
!-----------------------------------------------------------------------

irnuplt                   = 0
nrnuplt                   = 0
nrnu                      = 0
intrnu                    = 0
ncyrnu                    = 0

!-----------------------------------------------------------------------
!  Parameters for configuration plot edits
!-----------------------------------------------------------------------

ivarplt                   = 0
nvar                      = 0
nvarint                   = 0
nvarplt                   = 0
nvardump                  = 0

dtvarplot                 = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for comparison plot edits
!-----------------------------------------------------------------------

icomplt                   = 0
ncomplt                   = 0
ncomdump                  = 0

dtcomplot                 = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for bounday plot edits
!-----------------------------------------------------------------------

nplotinnerb               = 0
iplotinnerb               = 0
nplotouterb               = 0
iplotouterb               = 0
nplotlum                  = 0
iplotlum                  = 0
nplotshk                  = 0
iplotshk                  = 0
nplotcnv                  = 0
iplotcnv                  = 0
nplotmss                  = 0
iplotmss                  = 0

dtimeplot                 = 1.d+100
rinnerb                   = 1.d+100
routerb                   = 1.d+100
r_lumerms                 = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for nurad plot edits
!-----------------------------------------------------------------------

nplotnurad                = 0
iplotnurad                = 0

dtnuradplot               = 1.d+100
r_nurad                   = 1.d+100
rho_nurad                 = 1.d+100
nu_r                      = zero
nu_rt                     = zero
nu_rho                    = zero
nu_rhot                   = zero
unu_r                     = zero
unu_rt                    = zero
unu_rho                   = zero
unu_rhot                  = zero

!-----------------------------------------------------------------------
!  Parameters for neutrino distribution plot edits
!-----------------------------------------------------------------------

nnudata                   = 0
inudata                   = 0

dtnudata                  = 1.d+100
r_nudata                  = 1.d+100
t_nudata                  = 1.d+100
psi0dat                   = zero
psi1dat                   = zero

!-----------------------------------------------------------------------
!  Parameters for lagrangian plot edits
!-----------------------------------------------------------------------

nlagplt                   = 0
ilagplt                   = 0
nlagdump                  = 0

msslag                    = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for r-lagrangian plot edits
!-----------------------------------------------------------------------

nrlagplt                  = 0
irlagplt                  = 0

dmlag                     = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for energy-check plot edits
!-----------------------------------------------------------------------

n_eplt                    = 0
i_eplt                    = 0

dt_eplot                  = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional u (radial) edits
!-----------------------------------------------------------------------

iedMDu                    = 0
nedMDu                    = 0
intedMDu                  = 200
n_editMDu                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional v (angular) edits
!-----------------------------------------------------------------------

iedMDv                    = 0
nedMDv                    = 0
intedMDv                  = 200
n_editMDv                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional v (azimuthal) edits
!-----------------------------------------------------------------------

iedMDw                    = 0
nedMDw                    = 0
intedMDw                  = 200
n_editMDw                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional s (entropy) edits
!-----------------------------------------------------------------------

iedMDs                    = 0
nedMDs                    = 0
intedMDs                  = 200
n_editMDs                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional d (density) edits
!-----------------------------------------------------------------------

iedMDd                    = 0
nedMDd                    = 0
intedMDd                  = 200
n_editMDd                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional e (internal energy) edits
!-----------------------------------------------------------------------

iedMDe                    = 0
nedMDe                    = 0
intedMDe                  = 200
n_editMDe                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional p (pressure) edits
!-----------------------------------------------------------------------

iedMDp                    = 0
nedMDp                    = 0
intedMDp                  = 200
n_editMDp                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional enu (neutrino energy density) edits
!-----------------------------------------------------------------------

iedMDenu                  = 0
nedMDenu                  = 0
intedMDenu                = 200
n_editMDenu               = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional enu (neutrino energy flux) edits
!-----------------------------------------------------------------------

iedMDfnu                  = 0
nedMDfnu                  = 0
intedMDfnu                = 200
n_editMDfnu               = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional a (mean mass number) edits
!-----------------------------------------------------------------------

iedMDa                    = 0
nedMDa                    = 0
intedMDa                  = 200
n_editMDa                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional x (parameter) edits
!-----------------------------------------------------------------------

iedMDx                    = 0
nedMDx                    = 0
intedMDx                  = 200
n_editMDx                 = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional yw (electron fraction) edits
!-----------------------------------------------------------------------

iedMDye                   = 0
nedMDye                   = 0
intedMDye                 = 200
n_editMDye                = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional cm (composition) edits
!-----------------------------------------------------------------------

iedMDcm                   = 0
nedMDcm                   = 0
intedMDcm                 = 200
n_editMDcm                = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional nu (neutrino energy deposition) edits
!-----------------------------------------------------------------------

iedMDnu                   = 0
nedMDnu                   = 0
intedMDnu                 = 200
n_editMDnu                = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional nu (nuclear energy deposition) edits
!-----------------------------------------------------------------------

iedMDnc                   = 0
nedMDnc                   = 0
intedMDnc                 = 200
n_editMDnc                = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional nl (neutrino luminosity) edits
!-----------------------------------------------------------------------

iedMDnl                   = 0
nedMDnl                   = 0
intedMDnl                 = 200
n_editMDnl                = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional ne (neutrino rms energy) edits
!-----------------------------------------------------------------------

iedMDne                   = 0
nedMDne                   = 0
intedMDne                 = 200
n_editMDne                = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional gx (x-gravitational acceleration) edits
!-----------------------------------------------------------------------

iedMDgx                   = 0
nedMDgx                   = 0
intedMDgx                 = 200
n_editMDgx                = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional gy (y-gravitational acceleration) edits
!-----------------------------------------------------------------------

iedMDgy                   = 0
nedMDgy                   = 0
intedMDgy                 = 200
n_editMDgy                = 0

!-----------------------------------------------------------------------
!  Parameters for 2-D BVw-edits (Brunt-Vaisala frequency)
!-----------------------------------------------------------------------

iedMDBVw                  = 0
nedMDBVw                  = 0
intedMDBVw                = 200
n_editMDBVw               = 0

!-----------------------------------------------------------------------
!  Parameters for MD-dimensional gy (y-gravitational acceleration) edits
!-----------------------------------------------------------------------

iedMDyl                   = 0
nedMDyl                   = 0
intedMDyl                 = 200
n_editMDyl                = 0

!-----------------------------------------------------------------------
!  Time parameters for MD edits
!-----------------------------------------------------------------------

i_editMD                  = 0
n_editMD                  = 0
dt_MDedit1                = 1.d+20
dt_MDedit2                = 1.d+20

!-----------------------------------------------------------------------
!  Parameters for multidimensional HDF edits
!-----------------------------------------------------------------------

i_HDFedit                 = 0
nd_HDFedit                = 31
n_HDFedit                 = 0
dt_HDFedit1               = 1.d+20
dt_HDFedit2               = 1.d+20

!-----------------------------------------------------------------------
!  Parameters for temperary restart dumps
!-----------------------------------------------------------------------

intrst                    = 9000000
nnrst                     = 0
nrstd1                    = 0
nrstd2                    = 0
nrstfl                    = 0
nouttmp                   = 0

!-----------------------------------------------------------------------
!  Parameters for permanent restart dumps
!-----------------------------------------------------------------------

noutpmt                   = 0
intprm                    = 9000000
nprm                      = 0
ncyrst                    = -1
ncychg                    = 0
irstbgn                   = 0
noutfl                    = 0

!-----------------------------------------------------------------------
!  Work
!-----------------------------------------------------------------------

pdv                       = zero
twrk                      = zero

!-----------------------------------------------------------------------
!  Stresses
!-----------------------------------------------------------------------

nustrss                   = zero
pstrss                    = zero
gstrss                    = zero
gstrss_cx                 = zero
gstrss_cy                 = zero
rstrss                    = zero

!-----------------------------------------------------------------------
!  Gravitational potential
!-----------------------------------------------------------------------

g_pot                     = zero

!-----------------------------------------------------------------------
!  Neutron star parameters
!-----------------------------------------------------------------------

mass_ns                   = zero
vel_ns                    = zero

!-----------------------------------------------------------------------
!  Solid angles subtended by rays
!-----------------------------------------------------------------------

d_omega                   = zero

!-----------------------------------------------------------------------
!  Angularly averaged quantities
!-----------------------------------------------------------------------

rhobar                    = zero

!-----------------------------------------------------------------------
!  Neutrino-matter energy transfer rates
!-----------------------------------------------------------------------

dudt_ABEM                 = zero
dudt_NES                  = zero
dudt_NNS                  = zero
dudt_NAS                  = zero
dudt_PR                   = zero
dudt_Brem                 = zero
dudt_NET                  = zero

WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_edit_arrays
