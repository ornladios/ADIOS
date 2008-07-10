SUBROUTINE reset_eos_y_inout( jmin, jmax, ny, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, ls, le, ldim, nprintp, idtyp, nsep, rhop, tp, yep, eip, xnp, &
& be_nuc_repp, a_nuc_repp, z_nuc_repp, reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         reset_eos_y_inout
!    Module:       reset_eos_y_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/05
!
!    Purpose:
!      To load the variables passed from radial_ray_module and to call
!       the subroutines that reload the eos rate tables.
!
!    Input arguments:
!
!  jmin              : inner y-array index
!  jmax              : outer y-array index
!  ny                : y-array extent
!  ji_ray            : x (radial) index of a specific y (angular) ray
!  jk_ray            : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim         : the number of radial zones on a processor after swapping with y
!  ik_ray_dim        : the number of z-zones on a processor before swapping with z
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
!  eos_y_reset       : resets eos tables
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  eos_snc_x_module, mdl_cnfg_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : msolar

USE eos_snc_y_module, ONLY: xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, idty, &
& aesv, nse
USE mdl_cnfg_y_module, ONLY : rho, ye, t
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos   ! composition EOS reset flag

INTEGER, INTENT(in)              :: jmin             ! minimum y-zone index
INTEGER, INTENT(in)              :: jmax             ! maximum y-zone index
INTEGER, INTENT(in)              :: ny               ! y-array extent

INTEGER, INTENT(in)              :: ji_ray           ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray           ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: j_ray_dim        ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of radial zones on a processor before swapping with z

INTEGER                          :: ls               ! minimum abundance array index
INTEGER                          :: le               ! maximum abundance array index
INTEGER                          :: ldim             ! abundance array extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: rhop        ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)          :: yep         ! electron fraction

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)                   :: idtyp       ! grid selector index

!-----------------------------------------------------------------------
!        Input - Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)                 :: nsep        ! nse-nonnse flag

REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: tp          ! temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: eip         ! internal energy (ergs g^{-1})
REAL(KIND=double), INTENT(inout), DIMENSION(ny,ldim,j_ray_dim,ik_ray_dim)  :: xnp         ! mass fractions
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: be_nuc_repp ! binding energy of the representative heavy nucleus (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: a_nuc_repp  ! mass number of the representative heavy nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: z_nuc_repp  ! charge number of the representative heavy nucleus

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

rho(jmin:jmax)               = rhop(jmin:jmax,ji_ray,jk_ray)
t  (jmin:jmax)               = tp  (jmin:jmax,ji_ray,jk_ray)
ye (jmin:jmax)               = yep (jmin:jmax,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Transfer abundances
!-----------------------------------------------------------------------

be_nuc_rep(jmin:jmax)        = be_nuc_repp(jmin:jmax,ji_ray,jk_ray)
a_nuc_rep (jmin:jmax)        = a_nuc_repp (jmin:jmax,ji_ray,jk_ray)
z_nuc_rep (jmin:jmax)        = z_nuc_repp (jmin:jmax,ji_ray,jk_ray)
nse(jmin:jmax,ji_ray,jk_ray) = nsep(jmin:jmax,ji_ray,jk_ray)

xn(jmin:jmax,ls:le)          = xnp(jmin:jmax,ls:le,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!
!           \\\\\ RESET EOS AND NEUTRINO RATE TABLES /////
!
!-----------------------------------------------------------------------

CALL eos_y_reset( jmin, jmax, rho, t, ye, ji_ray, jk_ray, reset_comp_eos )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return updated variables
!-----------------------------------------------------------------------

idtyp(jmin:jmax,ji_ray,jk_ray)       = idty(jmin:jmax,ji_ray,jk_ray)
nsep (jmin:jmax,ji_ray,jk_ray)       = nse (jmin:jmax,ji_ray,jk_ray)
tp   (jmin:jmax,ji_ray,jk_ray)       = t   (jmin:jmax)
eip  (jmin:jmax,ji_ray,jk_ray)       = aesv(jmin:jmax,2,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Return abundances
!-----------------------------------------------------------------------

be_nuc_repp(jmin:jmax,ji_ray,jk_ray) = be_nuc_rep(jmin:jmax)
a_nuc_repp (jmin:jmax,ji_ray,jk_ray) = a_nuc_rep (jmin:jmax)
z_nuc_repp (jmin:jmax,ji_ray,jk_ray) = z_nuc_rep (jmin:jmax)

xnp(jmin:jmax,ls:le,ji_ray,jk_ray)   = xn(jmin:jmax,ls:le)

RETURN
END SUBROUTINE reset_eos_y_inout
