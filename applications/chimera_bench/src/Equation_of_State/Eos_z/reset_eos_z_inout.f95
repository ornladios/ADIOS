SUBROUTINE reset_eos_z_inout( kmin, kmax, nz, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, ls, le, ldim, nprintp, idtyp, nsep, rhop, tp, yep, eip, xnp, &
& be_nuc_repp, a_nuc_repp, z_nuc_repp, reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         reset_eos_z_inout
!    Module:       reset_eos_z_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To load the variables passed from radial_ray_module and to call
!       the subroutines that reload the eos rate tables.
!
!    Input arguments:
!
!  kmin              : inner z-array index
!  kmax              : outer z-array index
!  nz                : z-array extent
!  ki_ray            : x (radial) index of a specific z (azimuthal) ray
!  kj_ray            : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim        : the number of y-zones on a processor before swapping with y
!  k_ray_dim         : the number of radial zones on a processor after swapping with z
!  ls                : inner abundance array index
!  le                : outer abundance array index
!  ldim              : abundance array extent
!  nprintp           : unit number to print diagnostics
!  nsep              : nuclear statistical equilibrium flag
!  rhop              : density after hydro advance (cm^{-3})
!  tp                : initial temperatures (K)
!  yep               : initial electron fractions
!  eip               : initial internal energies (ergs g^{-1})
!  xnp               : mass fractions
!  be_nuc_repp       : binding energy of mean auxiliary nucleus
!  a_nuc_repp        : mass number of mean auxiliary nucleus
!  z_nuc_repp        : charge number of mean auxiliary nucleus
!  reset_comp_eos    : composition EOS reset flag
!
!    Output arguments:
!
!  idtyp             : grid selectgor index
!  eip               : updated internal energy (ergs g^{-1})
!
!    Subprograms called:
!  eos_z_reset       : resets eos tables
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  eos_snc_z_module, mdl_cnfg_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : msolar

USE eos_snc_z_module, ONLY: xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, idty, &
& aesv, nse
USE mdl_cnfg_z_module, ONLY : rho, ye, t
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos   ! composition EOS reset flag

INTEGER, INTENT(in)              :: kmin             ! minimum y-zone index
INTEGER, INTENT(in)              :: kmax             ! maximum y-zone index
INTEGER, INTENT(in)              :: nz               ! z-array extent

INTEGER, INTENT(in)              :: ki_ray           ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray           ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! the number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim        ! the number of radial zones on a processor after swapping with z

INTEGER                          :: ls               ! minimum abundance array index
INTEGER                          :: le               ! maximum abundance array index
INTEGER                          :: ldim             ! abundance array extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)          :: rhop        ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)          :: yep         ! electron fraction

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)                   :: idtyp       ! grid selector index

!-----------------------------------------------------------------------
!        Input - Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)                 :: nsep        ! nse-nonnse flag

REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: tp          ! temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: eip         ! internal energy (ergs g^{-1})
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ldim,ij_ray_dim,k_ray_dim)  :: xnp         ! mass fractions
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: be_nuc_repp ! binding energy of the representative heavy nucleus (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: a_nuc_repp  ! mass number of the representative heavy nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: z_nuc_repp  ! charge number of the representative heavy nucleus

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!          \\\\\ SET UP FOR RESETING EOS AND RATE TABLES /////
!
!        Load variables received from radial_ray_module into
!         mdl_cnfg_module and eos_snc_x_module
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to mgfld arrays
!-----------------------------------------------------------------------

rho(kmin:kmax)               = rhop(kmin:kmax,kj_ray,ki_ray)
t  (kmin:kmax)               = tp  (kmin:kmax,kj_ray,ki_ray)
ye (kmin:kmax)               = yep (kmin:kmax,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Transfer abundances
!-----------------------------------------------------------------------

be_nuc_rep(kmin:kmax)        = be_nuc_repp(kmin:kmax,kj_ray,ki_ray)
a_nuc_rep (kmin:kmax)        = a_nuc_repp (kmin:kmax,kj_ray,ki_ray)
z_nuc_rep (kmin:kmax)        = z_nuc_repp (kmin:kmax,kj_ray,ki_ray)
nse(kmin:kmax,kj_ray,ki_ray) = nsep(kmin:kmax,kj_ray,ki_ray)

xn(kmin:kmax,ls:le)          = xnp(kmin:kmax,ls:le,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!
!           \\\\\ RESET EOS AND NEUTRINO RATE TABLES /////
!
!-----------------------------------------------------------------------

CALL eos_z_reset( kmin, kmax, rho, t, ye, ki_ray, kj_ray, reset_comp_eos )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return updated variables
!-----------------------------------------------------------------------

idtyp(kmin:kmax,kj_ray,ki_ray)       = idty(kmin:kmax,kj_ray,ki_ray)
nsep (kmin:kmax,kj_ray,ki_ray)       = nse (kmin:kmax,kj_ray,ki_ray)
tp   (kmin:kmax,kj_ray,ki_ray)       = t   (kmin:kmax)
eip  (kmin:kmax,kj_ray,ki_ray)       = aesv(kmin:kmax,2,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Return abundances
!-----------------------------------------------------------------------

be_nuc_repp(kmin:kmax,kj_ray,ki_ray) = be_nuc_rep(kmin:kmax)
a_nuc_repp (kmin:kmax,kj_ray,ki_ray) = a_nuc_rep (kmin:kmax)
z_nuc_repp (kmin:kmax,kj_ray,ki_ray) = z_nuc_rep (kmin:kmax)

xnp(kmin:kmax,ls:le,kj_ray,ki_ray)   = xn(kmin:kmax,ls:le)

RETURN
END SUBROUTINE reset_eos_z_inout
