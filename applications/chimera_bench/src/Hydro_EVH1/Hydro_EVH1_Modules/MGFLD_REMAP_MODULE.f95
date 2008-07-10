MODULE mgfld_remap_module
!=======================================================================
! Stores variables needed for MGFLD to remap to Eulerian grid
!=======================================================================

USE kind_module

SAVE

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: r0i         ! initial density
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dvoli       ! initial zone volume
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dmi         ! initial zone masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: t0i         ! initial zone temperatures
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ei0i        ! initial zone internal energies
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dei         ! internal energy increment (in ergs/g)  

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: r0l         ! density before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dvoll       ! zone volume before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dxl         ! zone widths before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: xal         ! zone positions before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dml         ! zone masses before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: u0l         ! zone velocities before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: t0l         ! zone temperatures before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ei0l        ! zone internal energies before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ye0l        ! zone electron fractions before remap

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: r0e         ! density after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dvole       ! zone volume after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dxe         ! zone widths after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: xae         ! zone positions after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dme         ! zone masses after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: u0e         ! zone velocities after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: t0e         ! zone temperatures after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ei0e        ! zone internal energies after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ye0e        ! zone electron fractions after remap

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)  :: psi0_re     ! remapped psi0's
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)    :: comp        ! remapped xn's

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: r           ! density (g cm^{-3}) (working array)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: temp        ! temperature (K) (working array)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: ye          ! electron fraction (working array)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: xa          ! radial coordinates after Lagr update (working array)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dx          ! zone thicknesses after Lagr update (working array)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: xa0         ! radial coordinates after Eul remap (working array)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: dx0         ! zone thicknesses after Eul remap (working array)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: eb          ! nuclear binding energy (ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: fluxbe      ! total nuclear binding energy transferred during remap (ergs)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: e_bind_zn0  ! total nuclear binding energy initially in each mass shell (ergs)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)      :: fluxye_comp ! total electron fraction transferred during remap (ergs)


END MODULE mgfld_remap_module
