SUBROUTINE remap_psi_z( ngeom, kmin, kmax, ki_ray, kj_ray, nz, nez, nnu, &
& k_shock, time, t_bounce, tb_dy_shift )
!-----------------------------------------------------------------------
!
!    File:         remap_psi_z
!    Module:       remap_psi_z
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
!  ngeom       : problem geometry glag
!  kmin        : minimum z-array index
!  kmax        : maximim z-array index
!  ki_ray      : x (radial) index of a specific z (azimuthal) ray
!  kj_ray      : y (angular) index of a specific z (azimuthal) ray
!  nz          : z-array extent
!  nez         : neutrino energy array extent
!  nnu         : neutrino flavor array extent
!  k_shock     : zones marked for additional diffusion
!  time        : elapsed time
!  t_bounce    : time from bounce
!  tb_dy_shift : time from bounce to trun off grid wiggle and/or psi0 diffusion
!
!    Output arguments:
!      none
!
!    Subprograms called:
!  parabola    : Computes piecewise parabolic fits to psi0
!  volume_zone : Computes volumes in the Lagrangian and Eulerian grid
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, evh1_global, evh1_sweep, evh1_zone,
!  mgfld_remap_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY : zero, epsilon

USE edit_module, ONLY : nlog
USE evh1_global, ONLY: v_diff
USE evh1_sweep, ONLY: radius
USE evh1_zone, ONLY : zparay
USE mgfld_remap_module, ONLY : psi0_re, xa, dx, xa0, dx0
USE nu_dist_module, ONLY : psi0
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                              :: ngeom       ! geometry index
INTEGER, INTENT(in)                              :: kmin        ! minimum z-array index
INTEGER, INTENT(in)                              :: kmax        ! maximim z-array index
INTEGER, INTENT(in)                              :: ki_ray      ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                              :: kj_ray      ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                              :: nz          ! z-array extent
INTEGER, INTENT(in)                              :: nez         ! neutrino energy array extent
INTEGER, INTENT(in)                              :: nnu         ! neutrino energy flavor extent
INTEGER, INTENT(in), DIMENSION(nz)               :: k_shock     ! zones marked for added y-diffusion

REAL(KIND=double), INTENT(in)                    :: time        ! elapsed time
REAL(KIND=double), INTENT(in)                    :: t_bounce    ! time from bounce
REAL(KIND=double), INTENT(in)                    :: tb_dy_shift ! time from bounce to trun off grid wiggle and/or psi0 diffusion

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                               :: var_name

LOGICAL                                          :: first = .true.
LOGICAL                                          :: l_diff_off

INTEGER                                          :: istat       ! allocation status flag
INTEGER                                          :: nn          ! neutrino flavor index
INTEGER                                          :: k           ! neutrino energy index
INTEGER                                          :: iz          ! angular zone index
INTEGER                                          :: nmin        ! minimum padded index
INTEGER                                          :: nmax        ! maximum padded index
INTEGER                                          :: ntot        ! total number of padded indices
INTEGER                                          :: n           ! padded aray index
INTEGER                                          :: nm1         ! n-1 padded aray index
INTEGER                                          :: i           ! unpadded aray index
INTEGER, PARAMETER                               :: n_ei = 4    ! number of inner zones to remap internal energy

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psil        ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi6        ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dpsi        ! izero neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psil_re     ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi6_re     ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dpsi_re     ! izero neutrino moment slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol        ! volume after Lagr step
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol0       ! volume after Eul remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol0_inv   ! inverse volume after Eul remap

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dm          ! dummy array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: delta       ! volume of overlapping subshells

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi         ! 1-D storage array for zero neutrino moments
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: fluxpsi     ! neutrinos contained in overlapping subshells

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: para        ! parabolic coefficients

REAL(KIND=double)                                :: deltx       ! grid displacement during Lagrangian step
REAL(KIND=double)                                :: fractn      ! half the zone Lagrangian and Eulerian frid overlap
REAL(KIND=double)                                :: fractn2     ! 1 - 4/3*fractn

REAL(KIND=double)                                :: third       ! 1/3
REAL(KIND=double)                                :: twothd      ! 2/3
REAL(KIND=double)                                :: fourthd     ! 4/3

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in remap_psi_z')
 2001 FORMAT (' Deallocation problem for array ',a10,' in remap_psi_z')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (psil(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psil_re(nz+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil_re   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi6_re(nz+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6_re   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi_re(nz+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi_re   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dvol(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol0(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol0_inv(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0_inv '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dm(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (delta(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxpsi(nz+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (para(10,nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

psil                    = zero
psi6                    = zero
dpsi                    = zero
psil_re                 = zero
psi6_re                 = zero
dpsi_re                 = zero

dvol                    = zero
dvol0                   = zero
dvol0_inv               = zero

dm                      = zero
delta                   = zero

psi                     = zero
fluxpsi                 = zero

para                    = zero

IF ( first ) THEN
  first              = .false.
  third              = 1.d0/3.d0
  twothd             = 2.d0 * third
  fourthd            = 4.d0 * third
END IF

fluxpsi              = zero

nmin                 = kmin + 6
nmax                 = kmax + 6
ntot                 = kmax + 12

!-----------------------------------------------------------------------
!
!  Generate interpolation functions for psi0.
!
!
!     a(x) = a     + x[ da   = a    (1 - x) ]
!             L,iz        iz    6,iz
!
!     x = ( xi - xi      )/dxi        xi      < xi < xi
!                  iz-1/2     iz        iz-1/2         iz+1/2
!
!-----------------------------------------------------------------------

DO nn = 1,nnu
  IF ( nnugp(nn)  ==  0 ) CYCLE
  DO k = 1,nnugp(nn)

    psi(1:kmax+12)   = psi0_re(1:kmax+12,k,nn)

    CALL parabola( nmin-3, nmax+3, ntot, zparay, psi, dpsi, psi6, psil, &
&    dm, 0, 0, ngeom )

    dpsi_re(1:kmax+12,k,nn) = dpsi(1:kmax+12)
    psi6_re(1:kmax+12,k,nn) = psi6(1:kmax+12)
    psil_re(1:kmax+12,k,nn) = psil(1:kmax+12)

  END DO ! k = 1,nnugp(nn)
END DO ! nn = 1,nnu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
!  Calculate the volume of the overlapping subshells (delta)
!-----------------------------------------------------------------------

IF ( ngeom == 0 ) THEN

  DO n = nmin, nmax+1
    delta(n)         = xa(n) - xa0(n)
  END DO

ELSE IF ( ngeom == 1 ) THEN

  DO n = nmin, nmax+1
    delta(n)         = xa(n) - xa0(n)
    delta(n)         = delta(n) * ( xa0(n) + .5d0 * delta(n) )
  END DO

ELSE IF ( ngeom == 2) THEN ! delta(n) = third * ( xa**3 - xa0**3 )

  DO n = nmin, nmax+1
    delta(n)         = xa(n) - xa0(n)
    delta(n)         = delta(n) * ( xa0(n) * ( xa0(n) + delta(n))  + third * delta(n)**2 )
  END DO

ELSE IF( ngeom == 3 ) THEN

  DO n = nmin, nmax+1
    delta(n)         = ( xa(n) - xa0(n) )  * radius
  END DO

ELSE IF( ngeom == 4 ) THEN

  DO n = nmin, nmax+1
    delta(n)         = ( cos(xa0(n)) - cos(xa(n)) ) * radius
  END DO

ELSE IF( ngeom == 5 ) THEN

DO n = nmin, nmax+1
  delta(n)           = ( xa(n) - xa0(n) ) * radius 
END DO

END IF

!-----------------------------------------------------------------------
!
!  Calculate. fluxpsi, proportional to the the total number in the 
!   subshell created by the overlap of the Lagrangian and Eulerian
!   grids.
!
!  If the zone face has moved to the right (deltx > 0), use the
!   integral from the left side of Lagrangian zone n (i.e., nm1 = n-1).
!
!  If the zone face has moved to the left (deltx < 0), use the 
!   integral from the right side of zone n.
!
!-----------------------------------------------------------------------

DO n = nmin+1, nmax
  deltx              = xa(n) - xa0(n)

!-----------------------------------------------------------------------
!
!        ...deltx > 0.0...
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
!        (iz-6)    |               |          (iz-5)
!                  n               n
!                   Eulr            Lagr
!
!                  <---- del x ---->
!
!       n-1 material in overlapping shell added to n
!
!-----------------------------------------------------------------------

  IF ( deltx >= 0.0d0 ) THEN

    nm1              = n - 1
    fractn           = 0.5d0 * deltx/dx(nm1)
    fractn2          = 1.d0 - fourthd * fractn

    DO k = 1,nnugpmx
      DO nn = 1,nnu
        IF ( nnugp(nn)  ==  0 ) CYCLE
        fluxpsi(n,k,nn) = ( psil_re(nm1,k,nn) + dpsi_re(nm1,k,nn) &
&                       - fractn * ( dpsi_re(nm1,k,nn) - fractn2 * psi6_re(nm1,k,nn) ) ) * delta(n)

      END DO ! nn = 1,nnu
    END DO ! k = 1,nnugpmx

  ELSE

!-----------------------------------------------------------------------
!
!        ...deltx < 0.0...
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
!       (iz-6)     |                 |          (iz-5)
!                  n                 n
!                   Lagr              Eulr
!
!
!                  <----- del x ----->
!
!         n material in overlapping shell added to n-1
!
!-----------------------------------------------------------------------

    fractn           = 0.5d0 * deltx/dx(n)
    fractn2          = 1.d0 + fourthd * fractn

    DO k = 1,nnugpmx
      DO nn = 1,nnu
        IF ( nnugp(nn)  ==  0 ) CYCLE
        fluxpsi(n,k,nn) = ( psil_re(n,k,nn) &
&                       - fractn * ( dpsi_re(n,k,nn) + fractn2 * psi6_re(n,k,nn) ) ) * delta(n)

      END DO ! nn = 1,nnu
    END DO ! k = 1,nnugpmx

  END IF
END DO

!-----------------------------------------------------------------------
!  Calculate volumes before and after Eul remap
!-----------------------------------------------------------------------

CALL volume_zone( ngeom, kmin, kmax, xa, dx, xa0, dx0, dvol, dvol0 )

!-----------------------------------------------------------------------
!  Add diffusion to a shock aligned parallel to the y-axis
!-----------------------------------------------------------------------

l_diff_off           = ( t_bounce > epsilon )  .and.  ( time - t_bounce > tb_dy_shift )

DO n = nmin,nmax
  iz                 = n - 6
  DO k = 1,nnugpmx
    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      IF ( k_shock(iz) == 1  .or.  ( .not. l_diff_off ) ) THEN
        fluxpsi(n+1,k,nn) = fluxpsi(n+1,k,nn) &
&                      - v_diff * ( psi0_re(n+1,k,nn) - psi0_re(n,k,nn) ) * DMIN1( dvol(n), dvol(n+1) )
      END IF ! k_shock(iz) == 1
    END DO ! nn = 1,nnu
  END DO ! k = 1,nnugpmx
END DO ! n = nmin,nmax

!-----------------------------------------------------------------------
!  Advect psi0 by moving the subshell quantities into the appropriate
!   Eulerian zone. 
!-----------------------------------------------------------------------

dvol0_inv(nmin:nmax) = 1.0d0/dvol0(nmin:nmax)

DO n = nmin, nmax
  DO k = 1,nnugpmx
    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      psi0_re(n,k,nn) = ( psi0_re(n,k,nn) * dvol(n) &
&             + fluxpsi(n,k,nn) - fluxpsi(n+1,k,nn) ) * dvol0_inv(n)
    END DO !  nn = 1,nnu
  END DO ! k = 1,nnugpmx
END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
!  Restore psi0
!-----------------------------------------------------------------------

DO i = 1,kmax
  DO k = 1,nnugpmx
    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      n = i + 6
      psi0(i+1,k,nn) = psi0_re(n,k,nn)
    END DO ! nn = 1,nnu
  END DO ! k = 1,nnugpmx
END DO ! i = 1,kmax

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (psil, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dpsi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psil_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil_re   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi6_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6_re   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dpsi_re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi_re   '; WRITE (nlog,2001) var_name; END IF

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

DEALLOCATE (psi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxpsi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi   '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (para, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE remap_psi_z















