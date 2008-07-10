SUBROUTINE pre_remap_comp_x( nleft, nright, imin, imax, rho, t, ye, &
& ij_ray, ik_ray, i_nnse, nnc )
!-----------------------------------------------------------------------
!
!    File:         pre_remap_comp_x
!    Module:       pre_remap_comp_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/03/05
!
!    Purpose:
!      To load the abundance mass fractions, xn, into a padded array, comp, with
!       ghost zones at either end, and to fill these ghost zone with values
!       depending on the boundary condition flags.
!      xn arrays are stored in nucbrn_module.
!
!    Subprograms called:
!  coordbc   : laods x-grid ghost zones
!  sweepbc_r : loads density, r, ghost zones
!  deflash_x : deflashes material in NSE to be advected into a zone in non-NSE
!
!    Input arguments:
!  nleft     : left-hand boundary condition flag
!  nright    : right-hand boundary condition flag
!  imin      : minimum x-array index
!  imax      : maximim x-array index
!  rho       : shifted matter density array (g/cm**3).
!  t         : shifted matter matter temperature array (K).
!  ye        : shifted matter matter electron fraction array.
!  ij_ray    : index denoting the j-index of a specific radial ray
!  ik_ray    : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!  i_nnse    : minimum padded zone index for composition
!
!    Subprograms called:
!  coordbc   : laods x-grid ghost zones
!  sweepbc_r : loads density, r, ghost zones
!
!    Include files:
!  kind_module, array_module, numerical_module
!  evh1_bound, evh1_sweep, eos_snc_x_module, mgfld_remap_module,
!  nucbrn_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY : nx
USE numerical_module, ONLY : zero

USE evh1_bound, ONLY : comp_bcl, comp_bcr, r_bcl, r_bcr
USE evh1_sweep, ONLY: u
USE eos_snc_x_module, ONLY : nse
USE mgfld_remap_module, ONLY : comp, r, xa, dx, xa0, dx0
USE nucbrn_module, ONLY : xn, nuc_number

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------
     
INTEGER, INTENT(in)                    :: nleft     ! left boundary geometry index
INTEGER, INTENT(in)                    :: nright    ! right boundary geometry index
INTEGER, INTENT(in)                    :: imin      ! minimum x-array index
INTEGER, INTENT(in)                    :: imax      ! maximim x-array index
INTEGER, INTENT(in)                    :: ij_ray    ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                    :: ik_ray    ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                    :: nnc       ! composition array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t   ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ye  ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out)                   :: i_nnse    ! ndex of first non-nse zone (1,imax)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1)                      :: eos_reset ! EOS reset flag

INTEGER                                :: n_nucp1   ! nuc_number + 1
INTEGER                                :: i         ! zone index
INTEGER                                :: nmin      ! minimum padded index
INTEGER                                :: nmax      ! maximum padded index
INTEGER                                :: ntot      ! total number of padded indices
INTEGER                                :: n         ! padded zone index
INTEGER                                :: nminc     ! minimum padded zone index for non_nse composition
INTEGER                                :: nmincm1   ! nminc - 1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nmin                          = imin + 6
nmax                          = imax + 6
ntot                          = imax + 12
n_nucp1                       = nuc_number + 1

comp(nmin:nmax,1:n_nucp1)     = 1.d-30

!-----------------------------------------------------------------------
!  Find nse - nonnse boundary - nse has been loaded during a prior call
!   to eos_reset. It is the radial zone-shifted nse flag.
!  i_nnse is the first unshifted radial for which nse = 0
!  If i_nnse = 0, all zones are in NSE, so return
!  If i_nnse = 1, all zones are in non-NSE. Set i_nnse = 2 to avoid 
!   an index of 0 at the lower boundary
!-----------------------------------------------------------------------

i_nnse                        = 0
DO i = 1,imax
  IF ( nse(i+1,ij_ray,ik_ray) == 0 ) THEN
    i_nnse                    = i
    EXIT
  END IF
END DO

IF ( i_nnse == 0 ) RETURN
IF ( i_nnse == 1 ) i_nnse = 2

!-----------------------------------------------------------------------
!  nminc is the first padded radial index for which nse = 0
!-----------------------------------------------------------------------

nminc                         = i_nnse + 6
nmincm1                       = nminc - 1

!-----------------------------------------------------------------------
!  Load padded comp array, xn --> comp
!
!  If nse(j) = 0, transfer xn(j,nc) to comp.
!  If nse(j) = 1, deflash zone and transfer xn(j,nc) to comp.
!  Deflash i_nnse-1 zone for left boundary values of xn(j,nc) in order
!   "feed" rightward advection
!-----------------------------------------------------------------------

DO i = i_nnse-1,imax
  n                           = i + 6
  IF ( nse(i+1,ij_ray,ik_ray) == 1 ) THEN
    eos_reset                 = 'n'
    CALL deflash_x( i+1, rho, t, ye, ij_ray, ik_ray, eos_reset )
  END IF
  comp(n,1:n_nucp1)           = xn(i+1,1:n_nucp1)
END DO

!-----------------------------------------------------------------------
!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : fixed (eg, u_bcl,p_bcr,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))
!
!  Find ghost coordinates
!-----------------------------------------------------------------------

CALL coord_bc( nleft, nright, nmin, nmax, xa, dx, xa0, dx0, imax+12 )
CALL sweepbc_r( nleft ,nright, nmin, nmax )

!-----------------------------------------------------------------------
!  Add six ghost zones to xn --> comp array
!-----------------------------------------------------------------------

!........Load left (inner) ghosts

IF ( nleft == 0 ) THEN       ! Zero Gradient 

  DO n = 1, 6
    r (nmin-n)                = r (nmincm1)
    comp(nmincm1-n,1:n_nucp1) = comp(nmincm1,1:n_nucp1)
  END DO

ELSE IF ( nleft == 1 ) THEN  ! Zero Gradient 

  DO n = 1, 6
    r (nmin-n)                = r (nmincm1)
    comp(nmincm1-n,1:n_nucp1) = comp(nmincm1,1:n_nucp1)
  END DO

ELSE IF ( nleft == 2 ) THEN  ! Externally Fixed

  DO n = 1, 6
    r (nmin-n)                = r_bcl
    comp(nmincm1-n,1:n_nucp1) = comp_bcl(1:n_nucp1)
  END DO

ELSE IF ( nleft == 3 ) THEN  ! Periodic

  DO n = 1, 6
    r (nmin-n)                = r (nmax+1-n)
    comp(nmincm1-n,1:n_nucp1) = comp(nmax+1-n,1:n_nucp1)
  END DO

END IF

!........ Load right (outer) ghosts

if ( nright == 0 ) THEN      ! symmetric accross right (outer) edge

  DO n = 1, 6
    comp(nmax+n,1:n_nucp1)    = comp(nmax+1-n,1:n_nucp1)
  END DO

ELSE IF ( nright == 1 ) THEN ! Zero Gradient 

  DO n = 1, 6
    r (nmax+n)                = r (nmax)
    comp(nmax+n,1:n_nucp1)    = comp(nmax,1:n_nucp1)
  END DO

ELSE IF ( nright == 2 ) THEN ! Externally Fixed

  DO n = 1, 6
    r (nmax+n)                = r_bcr
    comp(nmax+n,1:n_nucp1)    = comp_bcr(1:n_nucp1)
  END DO

!  CALL tgvndeye_sweep(nmax+1,nmax+6)

ELSE IF ( nright == 3 ) THEN ! Periodic

  DO n = 1, 6
    r (nmax+n)                = r (nmin+n-1)
    comp(nmax+n,1:n_nucp1)    = comp(nminc+n-1,1:n_nucp1)
  END DO

ELSE IF ( nright == 4 ) THEN ! Externally Fixed

  DO n = 1, 6
    r (nmax+n)                = r_bcr
    comp(nmax+n,1:n_nucp1)    = comp(nmax,1:n_nucp1)
  END DO

ELSE IF ( nright == 5 ) THEN ! Externally Fixed

  DO n = 1, 6
    r (nmax+n)                = r_bcr
    comp(nmax+n,1:n_nucp1)    = comp(nmax,1:n_nucp1)
  END DO

ELSE IF ( nright == 6 ) THEN ! Externally Fixed


  IF ( u(nmax) < zero  .and.  r(nmax) < r_bcr ) THEN
    r_bcr                     = r(nmax)
  ELSE IF ( u(nmax) >= zero ) THEN
    r_bcr                     = r(nmax)
  END IF

  DO n = 1, 6
    r (nmax+n)                = r_bcr
    comp(nmax+n,1:n_nucp1)    = comp(nmax,1:n_nucp1)
  END DO

END IF

RETURN
END SUBROUTINE pre_remap_comp_x
