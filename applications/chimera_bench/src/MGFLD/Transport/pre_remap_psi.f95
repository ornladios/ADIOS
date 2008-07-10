SUBROUTINE pre_remap_psi( nleft, nright, imin, imax, nnu )
!-----------------------------------------------------------------------
!
!    File:         pre_remap_psi
!    Module:       pre_remap_psi
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/03/05
!
!    Purpose:
!      To load psi0 ghost zones in preparation of the remap step.
!
!    Input arguments:
!  nleft     : left-hand boundary condition flag
!  nright    : right-hand boundary condition flag
!  imin      : minimum y-array index
!  imax      : maximim y-array index
!  nnu       : neutrino energy flavor extent
!
!    Output arguments:
!      none
!
!    Subprograms called:
!  coordbc   : laods x-grid ghost zones
!
!    Include files:
!  kind_module
!  evh1_bound, nu_energy_grid_module, mgfld_remap_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double

USE evh1_bound, ONLY : psi0_bcl
USE nu_energy_grid_module, ONLY : nnugp
USE mgfld_remap_module, ONLY : psi0_re, xa, dx, xa0, dx0

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------
     
INTEGER, INTENT(in)                    :: nleft    ! left boundary geometry flag
INTEGER, INTENT(in)                    :: nright   ! right boundary geometry flag
INTEGER, INTENT(in)                    :: imin     ! minimum x-array index
INTEGER, INTENT(in)                    :: imax     ! maximim x-array index
INTEGER, INTENT(in)                    :: nnu      ! neutrino energy flavor extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: nn       ! neutrino flavor index
INTEGER                                :: kk       ! neutrino energy index
INTEGER                                :: n        ! padded zone index
INTEGER                                :: nmin     ! minimum padded index
INTEGER                                :: nmax     ! maximum padded index
INTEGER                                :: ntot     ! total number of padded indices
INTEGER                                :: nmax1n   ! nmax+1-n
INTEGER                                :: nminn1   ! nmin+n-1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : fixed (eg, u_bcl,p_bcr,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nmin                          = imin + 6
nmax                          = imax + 6
ntot                          = imax + 12

!-----------------------------------------------------------------------
!  Find ghost coordinates
!-----------------------------------------------------------------------

CALL coord_bc( nleft, nright, nmin, nmax, xa, dx, xa0, dx0, imax+12 )

!-----------------------------------------------------------------------
!  Load left (inner) ghosts
!-----------------------------------------------------------------------

IF ( nleft == 0 ) THEN       ! symmetric accross left (inner) edge

  DO n = 1, 6
    nminn1                    = MIN( nmin + n - 1, nmax )
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmin-n,kk,nn) = psi0_re(nminn1,kk,nn)
      END DO
    END DO
  END DO

ELSE IF ( nleft == 1 ) THEN  ! Zero Gradient 

  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmin-n,kk,nn) = psi0_re(nmin,kk,nn)
      END DO
    END DO
  END DO

ELSE IF ( nleft == 2 ) THEN  ! Externally Fixed

  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmin-n,kk,nn) = psi0_bcl
      END DO
    END DO
  END DO

ELSE IF ( nleft == 3 ) THEN  ! Periodic

  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmin-n,kk,nn) = psi0_re(nmax+1-n,kk,nn)
      END DO
    END DO
  END DO

END IF

!-----------------------------------------------------------------------
!  Load right (outer) ghosts
!-----------------------------------------------------------------------

if ( nright == 0 ) THEN      ! symmetric accross right (outer) edge

  DO n = 1, 6
    nmax1n                    = MAX( nmax + 1 - n, nmin )
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmax+n,kk,nn) = psi0_re(nmax1n,kk,nn)
      END DO
    END DO
  END DO

ELSE IF ( nright == 1 ) THEN  ! Radially diluted

  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmax+n,kk,nn) = psi0_re(nmax,kk,nn) * xa(nmax)**2/xa(nmax+n)**2
      END DO
    END DO
  END DO

ELSE IF ( nright == 2 ) THEN ! Radially diluted


  DO n = 1, 6
     DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmax+n,kk,nn) = psi0_re(nmax,kk,nn) * xa(nmax)**2/xa(nmax+n)**2
      END DO
    END DO
  END DO

ELSE IF ( nright == 3 ) THEN ! Periodic

  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmax+n,kk,nn) = psi0_re(nmin+n-1,kk,nn)
      END DO
    END DO
  END DO

ELSE IF ( nright == 4 ) THEN ! Radially diluted 

  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmax+n,kk,nn) = psi0_re(nmax,kk,nn) * xa(nmax)**2/xa(nmax+n)**2
      END DO
    END DO
  END DO

ELSE IF ( nright == 5 ) THEN ! Radially diluted
  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmax+n,kk,nn) = psi0_re(nmax,kk,nn) * xa(nmax)**2/xa(nmax+n)**2
      END DO
    END DO
  END DO

ELSE IF ( nright == 6 ) THEN ! Radially diluted
  DO n = 1, 6
    DO nn = 1,nnu
      IF ( nnugp(nn) == 0 ) CYCLE
      DO kk = 1,nnugp(nn)
        psi0_re(nmax+n,kk,nn) = psi0_re(nmax,kk,nn) * xa(nmax)**2/xa(nmax+n)**2
      END DO
    END DO
  END DO

END IF

RETURN
END SUBROUTINE pre_remap_psi
