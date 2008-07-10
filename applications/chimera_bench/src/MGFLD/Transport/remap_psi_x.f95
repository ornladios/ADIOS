SUBROUTINE remap_psi_x( ngeom, imin, imax, ij_ray, ik_ray, nx, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         remap_psi_x
!    Module:       remap_psi_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/03/05
!
!    Purpose:
!      To remap the zero moment of the neutrino distribution function
!       to the updated grid.
!
!    Input arguments:
!  ngeom       : problem geometry flag
!  imin        : minimum x-array index
!  imax        : maximim x-array index
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!  nx          : x-array extent
!  nez         : neutrino energy array extent
!  nnu         : neutrino flavor array extent
!
!    Output arguments:
!      none
!
!    Subprograms called:
!  parabola    : Computes piecewise parabolic fits to psi0
!  volume_zone : Computes volumes in the Lagrangian and Eulerian grid
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, evh1_sweep, evh1_zone, mgfld_remap_module,
!  nu_dist_module, nu_energy_grid_module, prb_cntl_module
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY : zero, frpi
USE physcnst_module, ONLY : cvel

USE edit_module, ONLY : nlog
USE evh1_sweep, ONLY: radius, lapse_c, u
USE evh1_zone, ONLY : zparax
USE mgfld_remap_module, ONLY : psi0_re, xa, dx, xa0, dx0
USE nu_dist_module, ONLY : psi0, unukrad, nnukrad, nnurad, unurad, ncoefa, &
& ecoefa, stwt, psi1
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx, unubi, dunui
USE prb_cntl_module, ONLY: ireltrns

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                              :: ngeom    ! geometry index
INTEGER, INTENT(in)                              :: imin     ! minimum x-array index
INTEGER, INTENT(in)                              :: imax     ! maximim x-array index
INTEGER, INTENT(in)                              :: ij_ray   ! j-index of a radial ray
INTEGER, INTENT(in)                              :: ik_ray   ! k-index of a radial ray
INTEGER, INTENT(in)                              :: nx       ! x-array extent
INTEGER, INTENT(in)                              :: nez      ! neutrino energy array extent
INTEGER, INTENT(in)                              :: nnu      ! neutrino energy flavor extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                               :: var_name

LOGICAL                                          :: first = .true.

INTEGER                                          :: istat      ! allocation status flag
INTEGER                                          :: nn         ! neutrino flavor index
INTEGER                                          :: k          ! neutrino energy index
INTEGER                                          :: nmin       ! minimum padded index
INTEGER                                          :: nmax       ! maximum padded index
INTEGER                                          :: ntot       ! total number of padded indices
INTEGER                                          :: jr_min     ! minimum shifted radial index
INTEGER                                          :: jr_max     ! maximum shifted radial index
INTEGER                                          :: n          ! padded aray index
INTEGER                                          :: nm1        ! n-1 padded aray index
INTEGER                                          :: i          ! unshifted aray index
INTEGER, PARAMETER                               :: n_ei = 4   ! number of inner zones to remap internal energy

INTEGER                                          :: nk         ! neutrino energy padded aray index
INTEGER                                          :: nmink      ! minimum padded neutrino energy zone index
INTEGER                                          :: nmaxk      ! maximum padded neutrino energy zone index
INTEGER                                          :: ntotk      ! total number of neutrino energy padded indices
INTEGER                                          :: nminkn1    ! nmink+nk-1
INTEGER                                          :: nmaxk1n    ! nmaxk+1-nk

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi0l      ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi06      ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dps0i      ! izero neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi0l_re   ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi06_re   ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dps0i_re   ! izero neutrino moment slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol       ! volume after Lagr step
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol0      ! volume after remap to final grid
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol0_inv  ! inverse volume after remap to final grid

REAL(KIND=double), DIMENSION(300)                :: evol       ! volume of energy bins
REAL(KIND=double), DIMENSION(300)                :: evol_inv   ! 1/( volume of energy bins )
REAL(KIND=double), DIMENSION(300)                :: e3         ! (neutrino energy bin edge)**3

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dm         ! dummy array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: delta      ! volume of overlapping subshells

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi        ! 1-D storage array for zero neutrino moments
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: fluxpsi_ot ! neutrinos contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: fluxpsi_in ! neutrinos contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: fluxpsi_it ! neutrinos contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi0_tmp   ! zero neutrino moment for advecting in energy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: agrjmh     ! zone-centered lapse of transfered netrinos before remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: agrajmh    ! zone-centered lapse of transfered netrinos after remap

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: para       ! parabolic coefficients

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: deltx      ! grid displacement during Lagrangian step
REAL(KIND=double)                                :: fractn     ! half the zone Lagrangian and Eulerian frid overlap
REAL(KIND=double)                                :: fractn2    ! 1 - 4/3*fractn

REAL(KIND=double)                                :: third      ! 1/3
REAL(KIND=double)                                :: twothd     ! 2/3
REAL(KIND=double)                                :: fourthd    ! 4/3

REAL(KIND=double)                                :: beta_rel   ! (relative velocity of acceptor zone to donar zone)/c
REAL(KIND=double)                                :: cvel_inv   ! 1/c_vel

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi1_1D    ! 1-D storage array for the first neutrino moments
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: xk         ! padded neutrino energy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dk         ! padded neutrino energy width
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi1l      ! first neutrino moment at left edge of energy zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi16      ! first neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dps1i      ! first neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi1l_re   ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi16_re   ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dps1i_re   ! izero neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: para1      ! array of parabolic interpolation coefficients

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in remap_psi_x')
 2001 FORMAT (' Deallocation problem for array ',a10,' in remap_psi_x')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (psi0l(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0l     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi06(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi06     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dps0i(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dps0i     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0l_re(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0l_re  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi06_re(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi06_re  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dps0i_re(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dps0i_re  '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dvol(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol0_inv(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0_inv '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dm(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (delta(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (deltx(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'deltx     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxpsi_ot(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi_ot'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxpsi_in(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi_in'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxpsi_it(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi_it'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0_tmp(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_tmp  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrjmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrjmh    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agrajmh(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrajmh   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (para(10,nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi1_1D(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1_1D   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xk(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xk        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dk(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dk        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1l(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1l     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi16(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi16     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dps1i(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dps1i     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1l_re(nx+1,nez+12,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1l_re  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi16_re(nx+1,nez+12,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi16_re  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dps1i_re(nx+1,nez+12,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dps1i_re  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (para1(10,nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para1     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

psi0l                   = zero
psi06                   = zero
dps0i                   = zero
psi0l_re                = zero
psi06_re                = zero
dps0i_re                = zero

dvol                    = zero
dvol0                   = zero
dvol0_inv               = zero

dm                      = zero
delta                   = zero

psi                     = zero
fluxpsi_ot              = zero
fluxpsi_in              = zero

para                    = zero

IF ( first ) THEN
  first                 = .false.
  third                 = 1.d0/3.d0
  twothd                = 2.d0 * third
  fourthd               = 4.d0 * third
  cvel_inv              = 1.d0/cvel
  e3(1:nnugpmx+1)       = unubi(1:nnugpmx+1) * unubi(1:nnugpmx+1) * unubi(1:nnugpmx+1)
  evol(1:nnugpmx)       = third * ( e3(2:nnugpmx+1) - e3(1:nnugpmx) )
  evol_inv(1:nnugpmx)   = 1.d0/evol(1:nnugpmx)
END IF

nmin                    = imin + 6
nmax                    = imax + 6
ntot                    = imax + 12
jr_min                  = imin + 1
jr_max                  = imax + 1

nmink                   = 7
nmaxk                   = nnugpmx + 6
ntotk                   = nnugpmx + 12

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
!
!                \\\\\ ENERGY EDGE VALUES OF PSI1 /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Put neutrino energies in padded arrays
!-----------------------------------------------------------------------

xk                      = zero
dk                      = zero
xk(7:nnugpmx+6)         = unubi(1:nnugpmx)
dk(7:nnugpmx+6)         = dunui(1:nnugpmx)

!-----------------------------------------------------------------------
!  Impose lower boundary conditions
!-----------------------------------------------------------------------

DO nk = 1,6
  nminkn1               = MIN( nmink + nk - 1, nmaxk )
  dk  (nmink-nk)        = dk  (nminkn1)
  xk  (nmink-nk)        = xk  (nminkn1) - dk  (nmink-nk)
END DO

!-----------------------------------------------------------------------
!  Impose upper boundary conditions
!-----------------------------------------------------------------------

DO nk = 1,6
  nmaxk1n               = MAX( nmaxk + 1 - nk, nmink )
  dk  (nmaxk+nk)        = dk  (nmaxk1n)
  xk  (nmaxk+nk)        = xk  (nmaxk1n) + dk (nmaxk+nk-1)
END DO

!-----------------------------------------------------------------------
!  Generate interpolation functions.
!-----------------------------------------------------------------------

CALL paraset_nu( ntotk, para1, dk, xk, nmink-4, nmaxk+4, ngeom )

!-----------------------------------------------------------------------
!  Generate energy edge values of psi1
!-----------------------------------------------------------------------

DO nn = 1,nnu
  IF ( nnugp(nn)  ==  0 ) CYCLE
  DO i = imin, imax+1

    psi1_1D(nmink:nmaxk) = psi1(i,1:nnugpmx,nn)
    DO nk = 1,6
      nminkn1            = MIN( nmink + nk - 1, nmaxk )
      nmaxk1n            = MAX( nmaxk + 1 - nk, nmink )
      psi1_1D (nmink-nk) = psi1_1D(nminkn1)
      psi1_1D (nmaxk+nk) = psi1_1D(nmaxk1n)
    END DO ! nk = 1,6

    CALL parabola_nu( nmink-4, nmaxk+4, ntotk, para1, psi1_1D, dps1i,  &
& psi16, psi1l, dm, 0, 0, ngeom )

    dps1i_re(i,1:nnugpmx,nn) = dps1i(nmink:nmaxk)
    psi16_re(i,1:nnugpmx,nn) = psi16(nmink:nmaxk)
    psi1l_re(i,1:nnugpmx,nn) = psi1l(nmink:nmaxk)

  END DO ! i = imin, imax+1
END DO ! nn = 1,nnu

!-----------------------------------------------------------------------
!
!  Generate interpolation functions for psi0.
!
!
!     a(x) = a    + x[ da  = a   (1 - x) ]
!             L,i        i    6,i
!
!     x = ( xi - xi     )/dxi       xi     < xi < xi
!                  i-1/2     i        i-1/2         i+1/2
!
!-----------------------------------------------------------------------

DO nn = 1,nnu
  IF ( nnugp(nn)  ==  0 ) CYCLE
  DO k = 1,nnugp(nn)

    psi(1:imax+12)   = psi0_re(1:imax+12,k,nn)

    CALL parabola( nmin-3, nmax+3, ntot, zparax, psi, dps0i, psi06, psi0l, &
&    dm, 0, 0, ngeom )

    dps0i_re(1:imax+12,k,nn)  = dps0i(1:imax+12)
    psi06_re(1:imax+12,k,nn)  = psi06(1:imax+12)
    psi0l_re(1:imax+12,k,nn)  = psi0l(1:imax+12)

  END DO ! k = 1,nnugp(nn)
END DO ! nn = 1,nnu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
!  Calculate the volume of the overlapping subshells (delta)
!-----------------------------------------------------------------------

IF ( ngeom == 0 ) THEN

  DO n = nmin, nmax+1
    delta(n)            = xa(n) - xa0(n)
  END DO

ELSE IF ( ngeom == 1 ) THEN

  DO n = nmin, nmax+1
    delta(n)            = xa(n) - xa0(n)
    delta(n)            = delta(n) * ( xa0(n) + .5d0 * delta(n) )
  END DO

ELSE IF ( ngeom == 2) THEN ! delta(n) = third * ( xa**3 - xa0**3 )

  DO n = nmin, nmax+1
    delta(n)            = xa(n) - xa0(n)
    delta(n)            = delta(n) * ( xa0(n) * ( xa0(n) + delta(n) )  + third * delta(n)**2 )
  END DO

ELSE IF( ngeom == 3 ) THEN

  DO n = nmin, nmax+1
    delta(n)            = ( xa(n) - xa0(n) )  * radius
  END DO

ELSE IF( ngeom == 4 ) THEN

  DO n = nmin, nmax+1
    delta(n)            = ( cos(xa0(n)) - cos(xa(n)) ) * radius
  END DO

ELSE IF( ngeom == 5 ) THEN

DO n = nmin, nmax+1
  delta(n)              = ( xa(n) - xa0(n) ) * radius 
END DO

END IF

!-----------------------------------------------------------------------
!
!  Calculate. fluxpsi, proportional to the the total number in the 
!   subshell created by the overlap of the Lagrangian and Eulerian
!   grids.
!
!  If the zone face has moved to the right (deltx(n) > 0), use the
!   integral from the left side of Lagrangian zone n (i.e., nm1 = n-1).
!
!  If the zone face has moved to the left (deltx(n) < 0), use the 
!   integral from the right side of zone n.
!
!-----------------------------------------------------------------------

DO n = nmin, nmax + 1
  deltx(n)              = xa(n) - xa0(n)
  i                     = n - 6
  i                     = n - 5
  nm1                   = n - 1

!-----------------------------------------------------------------------
!
!        ...deltx(n) > 0.0...
!
!
!        a   (x) = al    + z(da    + a6 (1-z)),    z = (x - x      )/dx
!         n-1        n-1       n-1     n-1                    n-1      n-1
!                                                                Lagr     Lagr
!
!             x
!              n
!           /   Lagr
!     1    |
!   _____  |          a   (x) dx = al    + da    - del z [ da    - a6    ( 1 -  4/3 del z ) ]
!   del x  |           n-1           n-1     n-1             n-1     n-1
!          /
!            x        - del x
!             n
!              Lagr
!
!........Material is remaped outwards
!
!                  |   n-1 --> n   |
!                  |               |
!                  |  delta(n) > 0 |
!         n-1      |       --------|-->         n
!                  |               |          i=n-5
!
!                  i               i
!                   final           Lagr
!
!                  <---- del x ---->
!
!       n-1 material in overlapping shell added to n
!
!-----------------------------------------------------------------------

  IF ( deltx(n) >= 0.0d0 ) THEN

    fractn              = 0.5d0 * deltx(n)/dx(nm1)
    fractn2             = 1.d0 - fourthd * fractn
    beta_rel            = ( u(nm1) - u(n) ) * cvel_inv
!    beta_rel            = zero

!-----------------------------------------------------------------------
!  psi0_tmp is the interpolated psi0 in the region to be remapped
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        psi0_tmp(i,k,nn) = psi0l_re(nm1,k,nn) + dps0i_re(nm1,k,nn) &
&                        - fractn * ( dps0i_re(nm1,k,nn) - fractn2 * psi06_re(nm1,k,nn) )
      END DO ! k = 1,nnugp(nn)
    END DO !  nn = 1,nnu

!-----------------------------------------------------------------------
!  fluxpsi_ot(n-1,k,nn) is the number of (k,nn) neutrinos removed from
!   radial zone n-1
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        fluxpsi_ot(n-1,k,nn) = fluxpsi_ot(n-1,k,nn) + psi0_tmp(i,k,nn) * delta(n)
      END DO ! k = 1,nnugp(nn)
    END DO ! nn = 1,nnu

!-----------------------------------------------------------------------
!  Correct psi0_tmp for the change in lapse in being transferred from
!   radial zone n-1 to radial zone n if ireltrns > 0
!-----------------------------------------------------------------------

    IF ( ireltrns > 0 ) THEN

      agrjmh(i)            = lapse_c(nm1)
      agrajmh(i)           = lapse_c(n)
      CALL nu_energy_agr_advct_inout_x( i, i, nx, nez, nnu, agrjmh, &
&      agrajmh, psi0_tmp )

    END IF ! ireltrns > 0

!-----------------------------------------------------------------------
!  fluxpsi_in(n,k,nn) is the number of (k,nn) neutrinos added to
!   radial zone n
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        fluxpsi_it(n,k,nn) = psi0_tmp(i,k,nn) * delta(n)
      END DO ! nn = 1,nnu
      IF ( beta_rel >= zero ) THEN
!-----------------------------------------------------------------------
!             beta_rel = ( u(nm1) - u(n) ) * cvel_inv > 0
!                !   u(n-1)      !  !      u(n)     !
!                !    <--        ! -!->  <------    !
!                !               !  !               !
!                    neutrino spectrum shifted up
!             Use psi1 at right-hand side of donor zone
!-----------------------------------------------------------------------
        k                      = 1
        fluxpsi_in(n,k,nn)     = fluxpsi_in(n,k,nn)                   &
&                              + fluxpsi_it(n,k,nn)                   &
&                              - beta_rel * ( e3(k+1) * ( psi1l_re(i,k  ,nn) + dps1i_re(i,k  ,nn) ) ) &
&                              * delta(n) * evol_inv(k)
        DO k = 2,nnugp(nn) - 1
          fluxpsi_in(n,k,nn)   = fluxpsi_in(n,k,nn)                   &
&                              + fluxpsi_it(n,k,nn)                   &
&                              - beta_rel * ( e3(k+1) * ( psi1l_re(i,k  ,nn) + dps1i_re(i,k  ,nn) )   &
&                                           - e3(k  ) * ( psi1l_re(i,k-1,nn) + dps1i_re(i,k-1,nn) ) ) &
&                              * delta(n) * evol_inv(k)
        END DO ! k = 1,nnugp(nn) - 1
        k                      = nnugp(nn)
        fluxpsi_in(n,k,nn)     = fluxpsi_in(n,k,nn)                   &
&                              + fluxpsi_it(n,k,nn)                   &
&                              - beta_rel * ( - e3(k) * ( psi1l_re(i,k-1,nn) + dps1i_re(i,k-1,nn) ) ) &
&                              * delta(n) * evol_inv(k)
      ELSE ! beta_rel < zero
!-----------------------------------------------------------------------
!             beta_rel = ( u(nm1) - u(n) ) * cvel_inv < 0
!                !   u(n-1)      !  !      u(n)     !
!                ! <--------     ! -!->    <--      !
!                !               !  !               !
!                   neutrino spectrum shifted down
!              Use psi1 at laft-hand side of donor zone
!-----------------------------------------------------------------------
        k                      = 1
        fluxpsi_in(n,k,nn)     = fluxpsi_in(n,k,nn)                   &
&                              + fluxpsi_it(n,k,nn)                   &
&                              + beta_rel * ( e3(k+1) * psi1l_re(i,k+1,nn) )                          &
&                              * delta(n) * evol_inv(k)
        DO k = 2,nnugp(nn) - 1
          fluxpsi_in(n,k,nn)   = fluxpsi_in(n,k  ,nn)                 &
&                              + fluxpsi_it(n,k,nn)                   &
&                              + beta_rel * ( e3(k+1) * psi1l_re(i,k+1,nn)                            &
&                                           - e3(k  ) * psi1l_re(i,k  ,nn) )                          &
&                              * delta(n) * evol_inv(k)
        END DO ! k = 1,nnugp(nn)
        k                      = nnugp(nn)
        fluxpsi_in(n,k,nn)     = fluxpsi_in(n,k,nn)                   &
&                              + fluxpsi_it(n,k,nn)                   &
&                              + beta_rel * ( - e3(k) * psi1l_re(i,k  ,nn) )                          &
&                              * delta(n) * evol_inv(k)
      END IF ! beta_rel >= zero
    END DO

  ELSE ! deltx(n) < 0.0d0

!-----------------------------------------------------------------------
!
!        ...deltx(n) < 0.0...
!
!        a (x) = al  + z(da  + a6 (1-z))    z = (x - x      )/dx
!         n        n       n     n                    n         n
!                                                      Lagr      Lagr
!
!             x      + del x
!              n
!           /   Lagr
!     1    |
!   _____  |          a (x) dx = al  + del z [ da  + a6  ( 1 -  4/3 del z ) ]
!   del x  |           n           n             n     n
!          |  
!          /
!            x
!             n
!              Lagr
!
!
!........Material is remaped inwards
!
!                  |    n-1 <-- n    |
!                  |                 |
!                  |  delta(n) <  0  |
!        n-1    <--|---------        |            n
!       i=n-6      |                 |          i=n-5
!
!                  i                 i
!                   Lagr              final
!
!
!                  <----- del x ----->
!
!         n material in overlapping shell added to n-1
!
!-----------------------------------------------------------------------

    fractn              = 0.5d0 * ( - deltx(n) )/dx(n)
    fractn2             = 1.d0 - fourthd * fractn
    beta_rel            = ( u(n) - u(nm1) ) * cvel_inv
!    beta_rel            = zero

!-----------------------------------------------------------------------
!  For deltx(n) < 0.0, fluxpsi_ot(n) are the neutrinos remapped into n-1
!   from n, fluxpsi_in(n) are the neutrinos remapped out of n.
!   Redistribute the neutrinos in fluxpsi_ot because of the different
!   lapses in n-1 and n.
!
!  fluxpsi_ot(n,k,nn) : number of (k,nn) neutrinos removed from zone n
!  fluxpsi_in(n,k,nn) : number of (k,nn) neutrinos added to zone n
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  psi0_tmp is the interpolated psi0 in the region to be remapped
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        psi0_tmp(i,k,nn) = psi0l_re(n,k,nn) &
&                        + fractn * ( dps0i_re(n,k,nn) + fractn2 * psi06_re(n,k,nn) )
      END DO ! k = 1,nnugp(nn)
    END DO ! nn = 1,nnu

!-----------------------------------------------------------------------
!  fluxpsi_ot(n,k,nn) is the number of (k,nn) neutrinos removed from
!   radial zone n
!  delta(n) < 0
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        fluxpsi_ot(n,k,nn) = fluxpsi_ot(n,k,nn) + psi0_tmp(i,k,nn) * ( - delta(n) )
      END DO ! k = 1,nnugp(nn)
    END DO ! nn = 1,nnu

!-----------------------------------------------------------------------
!  Correct psi0_tmp for the change in lapse in being transferred from
!   radial zone n to radial zone n-1 if ireltrns > 0
!-----------------------------------------------------------------------

    IF ( ireltrns > 0 ) THEN

      agrjmh(i)            = lapse_c(n)
      agrajmh(i)           = lapse_c(n-1)
      CALL nu_energy_agr_advct_inout_x( i, i, nx, nez, nnu, agrjmh, &
&      agrajmh, psi0_tmp )

    END IF ! ireltrns > 0

!-----------------------------------------------------------------------
!  fluxpsi_in(n-1,k,nn) is the number of (k,nn) neutrinos added to
!   radial zone n-1
!  delta(n) < 0
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        fluxpsi_it(n,k,nn) = psi0_tmp(i,k,nn) * ( - delta(n) )
      END DO ! nn = 1,nnu
      IF ( beta_rel >= zero ) THEN
!-----------------------------------------------------------------------
!             beta_rel = ( u(n) - u(nm1) ) * cvel_inv > 0
!                !   u(n-1)      !  !      u(n)     !
!                !   <-----    <-!- !      <--      !
!                !               !  !               !
!                    neutrino spectrum shifted up
!             Use psi1 at right-hand side of donor zone
!-----------------------------------------------------------------------
        k                      = 1
        fluxpsi_in(n-1,k,nn)   = fluxpsi_in(n-1,k,nn)                 &
&                              + fluxpsi_it(n,k,nn)                   &
&                              - beta_rel * ( e3(k+1) * ( psi1l_re(i,k  ,nn) + dps1i_re(i,k  ,nn) ) ) &
&                              * ( - delta(n) ) * evol_inv(k)
        DO k = 2,nnugp(nn) - 1
          fluxpsi_in(n-1,k,nn) = fluxpsi_in(n-1,k,nn)                 &
&                              + fluxpsi_it(n,k,nn)                   &
&                              - beta_rel * ( e3(k+1) * ( psi1l_re(i,k  ,nn) + dps1i_re(i,k  ,nn) )   &
&                                           - e3(k  ) * ( psi1l_re(i,k-1,nn) + dps1i_re(i,k-1,nn) ) ) &
&                              * ( - delta(n) ) * evol_inv(k)
        END DO ! k = 1,nnugp(nn) - 1
        k                      = nnugp(nn)
        fluxpsi_in(n-1,k,nn)   = fluxpsi_in(n-1,k,nn)                 &
&                              + fluxpsi_it(n,k,nn)                   &
&                              - beta_rel * ( - e3(k) * ( psi1l_re(i,k-1,nn) + dps1i_re(i,k-1,nn) ) ) &
&                              * ( - delta(n) ) * evol_inv(k)
      ELSE ! beta_rel < zero
!-----------------------------------------------------------------------
!             beta_rel = ( u(n) - u(nm1) ) * cvel_inv < 0
!                !   u(n-1)      !  !      u(n)     !
!                !   ----->    <-!- !      <--      !
!                !               !  !               !
!                   neutrino spectrum shifted down
!              Use psi1 at laft-hand side of donor zone
!-----------------------------------------------------------------------
        k                      = 1
        fluxpsi_in(n-1,k,nn)   = fluxpsi_in(n-1,k,nn)                 &
&                              + fluxpsi_it(n,k,nn)                   &
&                              + beta_rel * ( e3(k+1) * psi1l_re(i,k+1,nn) )                          &
&                              * ( - delta(n) )  * evol_inv(k)
        DO k = 2,nnugp(nn) - 1
          fluxpsi_in(n-1,k,nn) = fluxpsi_in(n-1,k,nn)                 &
&                              + fluxpsi_it(n,k,nn)                   &
&                              + beta_rel * ( e3(k+1) * psi1l_re(i,k+1,nn)                            &
&                                           - e3(k  ) * psi1l_re(i,k  ,nn) )                          &
&                              * ( - delta(n) ) * evol_inv(k)
        END DO ! k = 1,nnugp(nn)
        k                      = nnugp(nn)
        fluxpsi_in(n-1,k,nn)   = fluxpsi_in(n-1,k,nn)                 &
&                              + fluxpsi_it(n,k,nn)                   &
&                              + beta_rel * ( - e3(k) * psi1l_re(i,k  ,nn) )                          &
&                              * delta(n) * evol_inv(k)
      END IF ! beta_rel >= zero
    END DO ! nn = 1,nnu

  END IF ! deltx(n) >= 0.0d0
END DO ! n = nmin, nmax + 1

!-----------------------------------------------------------------------
!  Calculate volumes before and after Eul remap
!-----------------------------------------------------------------------

CALL volume_zone( ngeom, imin, imax, xa, dx, xa0, dx0, dvol, dvol0 )

!-----------------------------------------------------------------------
!  Advect psi0 by moving the subshell quantities into the appropriate
!   Eulerian zone. 
!-----------------------------------------------------------------------

dvol0_inv(nmin:nmax) = 1.0d0/dvol0(nmin:nmax)

DO n = nmin, nmax
  DO k = 1, nnugpmx
    DO nn = 1, nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      psi0_re(n,k,nn)   = ( psi0_re(n,k,nn) * dvol(n) &
&             + fluxpsi_in(n,k,nn) - fluxpsi_ot(n,k,nn) ) * dvol0_inv(n)
    END DO ! nn = 1,nnu
  END DO ! k = 1,nnugp(nn)
END DO ! n = nmin,nmax

!-----------------------------------------------------------------------
!  Restore psi0
!-----------------------------------------------------------------------

DO i = 1,imax
  n = i + 6
  psi0(i+1,:,:)            = psi0_re(n,:,:)
END DO ! i = 1,imax

!-----------------------------------------------------------------------
!  Prevent overfilling of neutrino states
!-----------------------------------------------------------------------

agrjmh(imin+1:imax+1)  = lapse_c(nmin:nmax)
CALL rebal_agr( imin+1, imax+1, nx, nez, nnu, psi0, agrjmh )

!-----------------------------------------------------------------------
!  Keep track of neutrinos entering or leaving the grid by updating
!   unukrad(k,n,ij_ray,ik_ray), nnukrad(k,n,ij_ray,ik_ray),
!   nnurad(n,ij_ray,ik_ray), unulum(n), and unurad(n,ij_ray,ik_ray).
!-----------------------------------------------------------------------

DO nn = 1,nnu
  IF ( nnugp(nn)  ==  0 ) CYCLE
  DO k = 1,nnugp(nn)
    IF ( deltx(nmax+1) >= zero ) THEN
      nnukrad(k,nn,ij_ray,ik_ray) = nnukrad(k,nn,ij_ray,ik_ray)         &
&                                 + stwt(nn) * frpi * ncoefa(imax+1,k) * fluxpsi_ot(nmax,k,nn)
      unukrad(k,nn,ij_ray,ik_ray) = unukrad(k,nn,ij_ray,ik_ray)         &
&                                 + stwt(nn) * frpi * ecoefa(imax+1,k) * fluxpsi_ot(nmax,k,nn)
      nnurad(nn,ij_ray,ik_ray)    = nnurad(nn,ij_ray,ik_ray)            &
&                                 + stwt(nn) * frpi * ncoefa(imax+1,k) * fluxpsi_ot(nmax,k,nn)
      unurad(nn,ij_ray,ik_ray)    = unurad(nn,ij_ray,ik_ray)            &
&                                 + stwt(nn) * frpi * ecoefa(imax+1,k) * fluxpsi_ot(nmax,k,nn)
    ELSE
      nnukrad(k,nn,ij_ray,ik_ray) = nnukrad(k,nn,ij_ray,ik_ray)         &
&                                 - stwt(nn) * frpi * ncoefa(imax+1,k) * fluxpsi_in(nmax,k,nn)
      unukrad(k,nn,ij_ray,ik_ray) = unukrad(k,nn,ij_ray,ik_ray)         &
&                                 - stwt(nn) * frpi * ecoefa(imax+1,k) * fluxpsi_in(nmax,k,nn)
      nnurad(nn,ij_ray,ik_ray)    = nnurad(nn,ij_ray,ik_ray)            &
&                                 - stwt(nn) * frpi * ncoefa(imax+1,k) * fluxpsi_in(nmax,k,nn)
      unurad(nn,ij_ray,ik_ray)    = unurad(nn,ij_ray,ik_ray)            &
&                                 - stwt(nn) * frpi * ecoefa(imax+1,k) * fluxpsi_in(nmax,k,nn)
    END IF ! deltx(nmax+1) >= zero
  END DO !  k = 1,nnugp(nn)
END DO ! nn = 1,nnu

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (psi0l, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0l     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi06, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi06     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dps0i, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0l     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi0l_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0l_re  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi06_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi06_re  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dps0i_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dps0i_re  '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dvol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol0_inv, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0_inv '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dm, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (delta, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (deltx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'deltx     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (psi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxpsi_ot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi_ot'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxpsi_in, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi_in'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxpsi_it, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi_it'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi0_tmp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_tmp  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (agrjmh, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrjmh    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (agrajmh, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agrajmh   '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (para, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (psi1_1D, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1_1D   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (xk, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xk        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dk, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dk        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi1l, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1l     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi16, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi16     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dps1i, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dps1i     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi1l_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1l_re  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi16_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi16_re  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dps1i_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dps1i_re  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (para1, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para1     '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE remap_psi_x
