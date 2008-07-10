SUBROUTINE reset_eos_x_inout( imin, imax, idim, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, ls, le, ldim, nprintp, idtyp, rhop, tp, yep, eip, xnp, &
& be_nuc_repp, a_nuc_repp, z_nuc_repp, nsep, reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         reset_eos_x_inout
!    Module:       reset_eos_x_inout
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
!  imin              : inner x-array index
!  imax              : outer x-array index
!  idim              : x-array extent
!  ij_ray            : j-index of a radial ray
!  ik_ray            : k-index of a radial ray
!  ij_ray_dim        : number of y-zones on a processor before swapping
!  ik_ray_dim        : number of z-zones on a processor before swapping
!  ls                : inner abundance array index
!  le                : outer abundance array index
!  ldim              : abundance array extent
!  nprintp           : unit number to print diagnostics
!  rhop              : density after hydro advance (cm^{-3})
!  tp                : initial temperatures (K)
!  yep               : initial electron fractions
!  eip               : initial internal energies (ergs g^{-1})
!  xnp               : mass fractions
!  be_nuc_repp       : binding energy of mean auxiliary nucleus
!  a_nuc_repp        : mass number of mean auxiliary nucleus
!  z_nuc_repp        : charge number of mean auxiliary nucleus
!  nsep              : nuclear statistical equilibrium flag
!  reset_comp_eos    : composition EOS reset flag
!
!    Output arguments:
!
!  idtyp             : grid selector index
!  eip               : updated internal energy
!  aesvp             : unshifted array of equation of state quantities
!
!    Subprograms called:
!  nse_test          : test for flashing from non_NSE tp NSE and vice versa
!  eos_reset         : resets eos tables
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  eos_snc_x_module, mdl_cnfg_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : msolar

USE eos_snc_x_module, ONLY: xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, idty, aesv, nse
USE mdl_cnfg_module, ONLY : rho, ye, t, jr_min, jr_max
USE nucbrn_module, ONLY: xn_n=>xn, be_nuc_rep_n=>be_nuc_rep, &
& a_nuc_rep_n=>a_nuc_rep, z_nuc_rep_n=>z_nuc_rep, nse_n=>nse

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos   ! composition EOS reset flag

INTEGER, INTENT(in)              :: imin             ! minimum x-zone index
INTEGER, INTENT(in)              :: imax             ! maximum x-zone index
INTEGER, INTENT(in)              :: idim             ! x-array extent

INTEGER, INTENT(in)              :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray           ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)              :: ls               ! minimum abundance array index
INTEGER, INTENT(in)              :: le               ! maximum abundance array index
INTEGER, INTENT(in)              :: ldim             ! abundance array extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)          :: rhop        ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)          :: yep         ! electron fraction

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(idim,ij_ray_dim,ik_ray_dim)                   :: idtyp       ! grid selector index

REAL(KIND=double), INTENT(out), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: eip         ! specific internal energy (ergs g^{-1})

!-----------------------------------------------------------------------
!        Input - Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(inout), DIMENSION(idim,ij_ray_dim,ik_ray_dim)                 :: nsep        ! nse-nonnse flag

REAL(KIND=double), INTENT(inout), DIMENSION(idim,ij_ray_dim,ik_ray_dim)       :: tp          ! temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(idim,ldim,ij_ray_dim,ik_ray_dim)  :: xnp         ! mass fractions
REAL(KIND=double), INTENT(inout), DIMENSION(idim,ij_ray_dim,ik_ray_dim)       :: be_nuc_repp ! binding energy of the representative heavy nucleus (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(idim,ij_ray_dim,ik_ray_dim)       :: a_nuc_repp  ! mass number of the representative heavy nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(idim,ij_ray_dim,ik_ray_dim)       :: z_nuc_repp  ! charge number of the representative heavy nucleus

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

jr_min                               = imin + 1
jr_max                               = imax + 1

!-----------------------------------------------------------------------
!
!          \\\\\ SET UP FOR RESETING EOS AND RATE TABLES /////
!
!        Load variables received from radial_ray_module into
!         mdl_cnfg_module and eos_snc_x_module
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to zone-shifted mgfld
!   arrays
!-----------------------------------------------------------------------

rho(jr_min:jr_max)                   = rhop(imin:imax,ij_ray,ik_ray)
t  (jr_min:jr_max)                   = tp  (imin:imax,ij_ray,ik_ray)
ye (jr_min:jr_max)                   = yep (imin:imax,ij_ray,ik_ray)

rho(jr_max+1)                        = rho(jr_max)
t  (jr_max+1)                        = t  (jr_max)
ye (jr_max+1)                        = ye (jr_max)

!-----------------------------------------------------------------------
!  Transfer abundances to zone-shifted mgfld arrays
!-----------------------------------------------------------------------

be_nuc_rep  (jr_min:jr_max)          = be_nuc_repp(imin:imax,ij_ray,ik_ray)
be_nuc_rep_n(jr_min:jr_max)          = be_nuc_repp(imin:imax,ij_ray,ik_ray)
a_nuc_rep   (jr_min:jr_max)          = a_nuc_repp (imin:imax,ij_ray,ik_ray)
a_nuc_rep_n (jr_min:jr_max)          = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep   (jr_min:jr_max)          = z_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep_n (jr_min:jr_max)          = z_nuc_repp (imin:imax,ij_ray,ik_ray)

be_nuc_rep  (jr_max+1)               = be_nuc_rep(jr_max)
be_nuc_rep_n(jr_max+1)               = be_nuc_rep(jr_max)
a_nuc_rep   (jr_max+1)               = a_nuc_rep (jr_max)
a_nuc_rep_n (jr_max+1)               = a_nuc_rep (jr_max)
z_nuc_rep   (jr_max+1)               = z_nuc_rep (jr_max)
z_nuc_rep_n (jr_max+1)               = z_nuc_rep (jr_max)

xn  (jr_min:jr_max,ls:le)            = xnp(imin:imax,ls:le,ij_ray,ik_ray)
xn_n(jr_min:jr_max,ls:le)            = xnp(imin:imax,ls:le,ij_ray,ik_ray)

xn  (jr_max+1,ls:le)                 = xn(jr_max,ls:le)
xn_n(jr_max+1,ls:le)                 = xn(jr_max,ls:le)

nse  (jr_min:jr_max,ij_ray,ik_ray)   = nsep(imin:imax,ij_ray,ik_ray)
nse_n(jr_min:jr_max)                 = nsep(imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!
!             \\\\\ TEST FOR FLASHING OR DEFLASHING /////
!
!-----------------------------------------------------------------------

CALL nse_test( jr_min, jr_max, rho, t, ye, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim )

!-----------------------------------------------------------------------
!
!           \\\\\ RESET EOS AND NEUTRINO RATE TABLES /////
!
!-----------------------------------------------------------------------

CALL eos_x_reset( idim, jr_min, jr_max, rho, t, ye, ij_ray, ik_ray, &
& reset_comp_eos )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return updated variables
!-----------------------------------------------------------------------

idtyp(imin:imax,ij_ray,ik_ray)       = idty(jr_min:jr_max,ij_ray,ik_ray)
nsep (imin:imax,ij_ray,ik_ray)       = nse (jr_min:jr_max,ij_ray,ik_ray)
tp   (imin:imax,ij_ray,ik_ray)       = t   (jr_min:jr_max)
eip  (imin:imax,ij_ray,ik_ray)       = aesv(jr_min:jr_max,2,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Return abundances
!-----------------------------------------------------------------------

be_nuc_repp(imin:imax,ij_ray,ik_ray) = be_nuc_rep(jr_min:jr_max)
a_nuc_repp (imin:imax,ij_ray,ik_ray) = a_nuc_rep (jr_min:jr_max)
z_nuc_repp (imin:imax,ij_ray,ik_ray) = z_nuc_rep (jr_min:jr_max)

xnp(imin:imax,ls:le,ij_ray,ik_ray)   = xn(jr_min:jr_max,ls:le)

RETURN
END SUBROUTINE reset_eos_x_inout
