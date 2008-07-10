SUBROUTINE pre_remap_comp_y( nleft, nright, jmin, jmax, rho, t, ye, &
& ji_ray, jk_ray, nnc )
!-----------------------------------------------------------------------
!
!    File:         pre_remap_comp_y
!    Module:       pre_remap_comp_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/03/05
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
!  deflash_y : switches material from NES to nonNSE
!
!    Input arguments:
!  nleft     : left-hand boundary condition flag
!  nright    : right-hand boundary condition flag
!  jmin      : minimum y-array index
!  jmax      : maximim y-array index
!  rho       : angular array of matter density (g/cm**3)
!  t         : angular array of matter temperature (K)
!  ye        : angular array of matter electron fraction
!  ji_ray    : x (radial) index of a specific y (angular) ray
!  jk_ray    : z (azimuthal) index of a specific y (angular) ray
!  nnc       : composition array extent
!
!    Output arguments:
!      none
!
!    Subprograms called:
!  coordbc   : laods x-grid ghost zones
!  sweepbc_r : loads density, r, ghost zones
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_y_module, evh1_bound, mgfld_remap_module,
!  nucbrn_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY : ny
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE eos_snc_y_module, ONLY : xn, nse
USE evh1_bound, ONLY : comp_bcl, comp_bcr, r_bcl, r_bcr
USE mgfld_remap_module, ONLY : comp, r, xa, dx, xa0, dx0, e_bind_zn0, eb, fluxbe
USE nucbrn_module, ONLY : nuc_number

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------
     
INTEGER, INTENT(in)                    :: nleft     ! left boundary geometry index
INTEGER, INTENT(in)                    :: nright    ! right boundary geometry index
INTEGER, INTENT(in)                    :: jmin      ! minimum y-array index
INTEGER, INTENT(in)                    :: jmax      ! maximim y-array index
INTEGER, INTENT(in)                    :: ji_ray    ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)                    :: jk_ray    ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)                    :: nnc       ! composition array extent

REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rho ! angular array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: t   ! angular array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: ye  ! angular array of matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: n_nucp1   ! nuc_number + 1
INTEGER                                :: nmin      ! minimum padded index
INTEGER                                :: nmax      ! maximum padded index
INTEGER                                :: ntot      ! total number of padded indices
INTEGER                                :: j         ! y-zone index
INTEGER                                :: n         ! padded zone index
INTEGER                                :: nminc     ! minimum padded zone index for composition
INTEGER                                :: nminn1    ! MIN( nmin + n - 1, nmax )
INTEGER                                :: nmax1n    ! MAX( nmax + 1 - n, nmin )
INTEGER                                :: nse_zn    ! 1 if all zones or in nse; otherwise 0

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' nrighty=',i3,' in subroutine pre_remap_comp_y; not an allowed value.')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nmin                       = jmin + 6
nmax                       = jmax + 6
ntot                       = jmax + 12
n_nucp1                    = nuc_number + 1

!-----------------------------------------------------------------------
!  Determine nse_zn (1 if all zones or in nse; otherwise 0)
!-----------------------------------------------------------------------

nse_zn                     = 1
DO j = jmin,jmax
  IF ( nse(j,ji_ray,jk_ray) == 0 ) THEN
    nse_zn                 = 0
    EXIT
  END IF
END DO

!-----------------------------------------------------------------------
!  If nse_zn = 1, all zones are in NSE
!  Set e_bind_zn0(n) = eb = fluxbe = 0 and return
!-----------------------------------------------------------------------

IF ( nse_zn == 1 ) THEN
  eb                       = zero
  e_bind_zn0               = zero
  fluxbe                   = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Initialize comp
!-----------------------------------------------------------------------

comp(nmin:nmax,1:n_nucp1) = 1.d-30

!-----------------------------------------------------------------------
! Store xn --> comp array padded with six ghost zones
!
!  If nse(j) = 0, transfer xn(j,nc) to comp.
!  If nse(j) = 1, deflash zone and transfer xn(j,nc) to comp.
!-----------------------------------------------------------------------

DO j = jmin,jmax
  n                        = j + 6
  IF ( nse(j,ji_ray,jk_ray) == 1 ) CALL deflash_y( j, rho, t, ye, ji_ray, &
&  jk_ray)
  comp(n,1:n_nucp1)        = xn(j,1:n_nucp1)
END DO

!-----------------------------------------------------------------------
!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : fixed (eg, u_bcl,p_bcr,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find ghost coordinates
!-----------------------------------------------------------------------

CALL coord_bc( nleft, nright, nmin, nmax, xa, dx, xa0, dx0, ntot )
CALL sweepbc_r( nleft, nright, nmin, nmax )

!........Load left (inner) ghosts

IF ( nleft == 0 ) THEN       ! symmetric accross right (outer) edge

  DO n = 1, 6
    nminn1                 = MIN( nmin + n - 1, nmax )
    r (nmin-n)             = r (nminn1)
    comp(nmin-n,1:n_nucp1) = comp(nminn1,1:n_nucp1)
  END DO

ELSE IF ( nleft == 1 ) THEN  ! Zero Gradient 

  DO n = 1, 6
    r (nmin-n)             = r (nmin)
    comp(nmin-n,1:n_nucp1) = comp(nmin,1:n_nucp1)
  END DO

ELSE IF ( nleft == 2 ) THEN  ! Externally Fixed

  DO n = 1, 6
    r (nmin-n)             = r_bcl
    comp(nmin-n,1:n_nucp1) = comp_bcl(1:n_nucp1)
  END DO

ELSE IF ( nleft == 3 ) THEN  ! Periodic

  DO n = 1, 6
    nmax1n                 = MAX( nmax+1-n, nmin )
    r (nmin-n)             = r (nmax1n)
    comp(nmin-n,1:n_nucp1) = comp(nmax1n,1:n_nucp1)
  END DO

END IF

!........ Load right (outer) ghosts

IF ( nright == 0 ) THEN      ! symmetric accross right (outer) edge

  DO n = 1, 6
    nmax1n                 = MAX( nmax+1-n, nmin )
    r (nmax+n)             = r (nmax1n)
    comp(nmax+n,1:n_nucp1) = comp(nmax1n,1:n_nucp1)
  END DO

ELSE IF ( nright == 1 ) THEN ! Zero Gradient 

  DO n = 1, 6
    r (nmax+n)             = r (nmax)
    comp(nmax+n,1:n_nucp1) = comp(nmax,1:n_nucp1)
  END DO

ELSE IF ( nright == 2 ) THEN ! Externally Fixed

  DO n = 1, 6
    r (nmax+n)             = r_bcr
    comp(nmax+n,1:n_nucp1) = comp_bcr(1:n_nucp1)
  END DO

ELSE IF ( nright == 3 ) THEN ! Periodic

  DO n = 1, 6
    r (nmax+n)             = r (nmin+n-1)
    comp(nmax+n,1:n_nucp1) = comp(nminc+n-1,1:n_nucp1)
  END DO

ELSE

  WRITE (nlog,1001) nright
  STOP

END IF

RETURN
END SUBROUTINE pre_remap_comp_y
