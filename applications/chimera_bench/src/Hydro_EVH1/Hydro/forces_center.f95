SUBROUTINE forces_center( ij_ray, ik_ray, xa, xf, uf, vf, wf, r, umid, &
& grav, fict, nuf, pseudo_fict )
!-----------------------------------------------------------------------
!
!    File:         forces_center
!    Module:       forces_center
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To calculate real (ie, gravity and radiation) and fictitious (i.e., 
!       coriolis and centrifugal) forces at the zone centers.
!
!    Input arguments:
!
!  ij_ray           : y-index of an x ray, x-index of an y ray, x index of a z array
!  ik_ray           : z-index of an x ray, z-index of an y ray, y index of a z array
!  xa               : padded zone edged coordinate
!  xf               : padded zone centered coordinate
!  uf               : x-velocity
!  vf               : y-velocity
!  wf               : z-velocity
!  umid             : padded zone-edged velocity
!  r                : density
!
!    Output arguments:
!
!  grav             : padded gravitational force
!  fict             : padded fictitious force
!  nuf              : padded neutrino stress
!  pseudo_fict      : padded pseudo fictitious force
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, evh1_global, evh1_sweep
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: max_12
USE numerical_module, ONLY : zero, third, half, one, epsilon
USE physcnst_module, ONLY : pi, G => g

USE edit_module, ONLY : nlog
USE evh1_global, ONLY : ngeomx, ngeomy, ngeomz, nleftx, small, i_radial
USE evh1_sweep, ONLY : sweep, theta, ctheta, stheta, radius, nmin, nmax, &
& dvol0, rhobar, nu_strc, g_force_c, g_pot_c, e_nu_c, f_nu_e

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                  :: ij_ray      ! y-index of an x ray, x-index of an y ray, x index of a z array
INTEGER, INTENT(in)                  :: ik_ray      ! z-index of an x ray, z-index of an y ray, y index of a z array

REAL(KIND=double), DIMENSION(max_12) :: xa          ! zone edge position
REAL(KIND=double), DIMENSION(max_12) :: xf          ! mass averaged mid zone position
REAL(KIND=double), DIMENSION(max_12) :: uf          ! x-velocity
REAL(KIND=double), DIMENSION(max_12) :: vf          ! y-velocity
REAL(KIND=double), DIMENSION(max_12) :: wf          ! z-velocity
REAL(KIND=double), DIMENSION(max_12) :: r           ! density
REAL(KIND=double), DIMENSION(max_12) :: umid        ! x-edge velocity

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(max_12) :: grav        ! gravitational force
REAL(KIND=double), DIMENSION(max_12) :: fict        ! fictitious force
REAL(KIND=double), DIMENSION(max_12) :: nuf         ! neutrino stress
REAL(KIND=double), DIMENSION(max_12) :: pseudo_fict ! pseudo fictitious force

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                              :: n           ! neutrino flavor index
INTEGER                              :: nsym        ! ghost zone index

REAL(KIND=double), DIMENSION(max_12) :: area        ! area of zone interfaces
REAL(KIND=double), DIMENSION(max_12) :: poteng      ! ptential energy

REAL(KIND=double)                    :: massi       ! mass interior to zone edge
REAL(KIND=double)                    :: massih      ! mass interior to zone center
REAL(KIND=double)                    :: massz       ! mass of shell
REAL(KIND=double)                    :: mass_co     ! central mass
REAL(KIND=double)                    :: xf0         ! minimum radius
REAL(KIND=double)                    :: sinxf0      ! DSIN(xf0)

REAL(KIND=double), PARAMETER         :: fopi = 4.d0 * pi

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

grav                            = zero
fict                            = zero
nuf                             = zero
pseudo_fict                     = zero

!-----------------------------------------------------------------------
!
!                  \\\\\ X-DIRECTION FORCES /////
!
!-----------------------------------------------------------------------

IF ( sweep == 'x') THEN

  IF ( ngeomx == 0  .and.  ngeomy == 0 ) then       ! Cartesian

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = zero
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomx == 0  .and.  ngeomy == 1 ) THEN  ! Cylindrical radial

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = wf(n) * wf(n)/xf(n)
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomx == 1  .and.  ngeomy == 3 ) THEN  ! Cylindrical polar

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = vf(n) * vf(n)/xf(n)
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomx == 2 ) THEN                      ! Spherical

!-----------------------------------------------------------------------
!  This is where the Poisson solution should get used
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!        ||||| Problems with a reflecting inner boundary |||||
!        |||||          (assumed to be at r=0)           |||||
!
!-----------------------------------------------------------------------

    IF ( nleftx == 0 ) THEN

!-----------------------------------------------------------------------
!  Initialize, fill inner ghosts and compute area (modulo 4pi)
!-----------------------------------------------------------------------

      massi                     = zero
      grav(nmin-4:nmin-1)       = zero
      fict(nmin-4:nmin-1)       = zero
      poteng(nmin-4:nmin-1)     = zero
      area(nmin-4:nmin-1)       = xf(nmin-4:nmin-1) * xf(nmin-4:nmin-1)
      area(nmin:nmax+5)         = xf(nmin:nmax+5) * xf(nmin:nmax+5)
      nuf                       = zero
      fict                      = zero

!-----------------------------------------------------------------------
!  Compute gravity
!
!  For zone center coordinates, interior mass includes portion of
!   current zone
!
!  For volume centered coordinates this is half of zone mass
!  For linear center, the fraction is smaller
!     hdx     = dx0(n) * 0.5
!     hdvol   = hdx*(xa0(n)*(xa0(n)+hdx)+hdx*hdx/3.)*4.0*pi
!     massih = massi + rhobar(n) * hdvol 
!
!  Gravitational forces
!-----------------------------------------------------------------------

      DO n = nmin, nmax
        grav(n)                 = g_force_c(n)
        poteng(n)               = g_pot_c(n)
      END DO
      grav  (nmax+1:nmax+5)     = grav  (nmax)
      poteng(nmax+1:nmax+5)     = poteng(nmax)

!-----------------------------------------------------------------------
!  Compute nongravitational forces
!
!  nuf(n)    : zone-centered neutrino stress
!  fict(n)   : fictitious forces
!-----------------------------------------------------------------------

      DO n = nmin, nmax+5
        nuf(n)                  = nu_strc(n)
        IF ( xf(n) /= zero ) THEN
          fict(n)               = ( wf(n) * wf(n) + vf(n) * vf(n) )/xf(n)
        END IF
      END DO ! n = nmin, nmax+5

!-----------------------------------------------------------------------
!  Set force at inner zone edge =0
!-----------------------------------------------------------------------

      fict(nmin)    = zero
             
!-----------------------------------------------------------------------
!  Impose symetric forces on interior ghosts
!-----------------------------------------------------------------------

      DO n = nmin-4,nmin-1
        nsym                    = nmin - n + nmin - 1
        grav(n)                 = -grav(nsym)
        fict(n)                 = -fict(nsym)
      END DO ! n = nmin-4, nmax-1

!-----------------------------------------------------------------------
!
!       ||||| Spherical Problems with the possibility of a |||||
!       |||||                central object                |||||
!
!-----------------------------------------------------------------------

    ELSE ! nleftx /= 0
             
!-----------------------------------------------------------------------
!  Start with the mass of the core or central object
!-----------------------------------------------------------------------

      massi                     = mass_co
             
!-----------------------------------------------------------------------
!  Subtract mass of inner ghosts
!-----------------------------------------------------------------------

      DO n = nmin-4,nmin-1
        massi                   = massi - rhobar(n) * dvol0(n) * fopi
      END DO ! n = nmin-4, nmax-1
             
!-----------------------------------------------------------------------
!  Gravity, including inner half of zone
!-----------------------------------------------------------------------

      DO n = nmin-4, nmax+5

        massz                   = rhobar(n) * dvol0(n) * fopi
             
!-----------------------------------------------------------------------
!  For zone center coordinates, interior mass includes half of current
!   zone
!-----------------------------------------------------------------------

        massih                  = massi + 0.5d0 * massz
             
!-----------------------------------------------------------------------
!  Update mass interior
!-----------------------------------------------------------------------

        massi                   = massi + massz  !  Update mass interior
             
!........Gravity

!        grav(n)     = -G * massih/xf(n)**2
        grav(n)                 = zero
             
!........Fictitios Forces

        fict(n)                 = zero
        fict(n)                 = ( wf(n) * wf(n) + vf(n) * vf(n) )/xf(n)

      END DO ! n = nmin-4, nmax+5

!-----------------------------------------------------------------------
!  for a star in which one may just USE the mass of the proto-neutron
!  star, the following may be used instead of the above
!
!macro       DO n = nmin-4, nmax+5
!macro         grav(n) = -GM_co / xf(n)**2
!macro         fict(n) = (wf(n)*wf(n)+vf(n)*vf(n))/xf(n)
!macro       END DO
!-----------------------------------------------------------------------

    END IF ! nleftx == 0
  END IF ! ngeomx == 2
END IF ! sweep == 'x'

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-----------------------------------------------------------------------
!
!                  \\\\\ Y-DIRECTION FORCES /////
!
!-----------------------------------------------------------------------

IF ( sweep == 'y' ) THEN

  IF( ngeomy == 0 ) THEN                            ! Cartesian

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      grav(n)                   = - grav(nsym)
      fict(n)                   = zero
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 1 ) THEN                      ! Cylindrical radial

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = wf(n) * wf(n)/xf(n)
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 3  .and.  ngeomx == 1 ) THEN  ! Cylindrical polar (2D)

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = - uf(n) * wf(n)/radius
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 3  .and.  ngeomx == 2 ) THEN  ! Spherical equator (2D)

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = - uf(n) * wf(n)/radius
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 4  .and.  ngeomx == 2 ) THEN  ! Spherical

    nuf                         = zero
    DO n = nmin-4, nmax+5
      IF ( xf(n) == zero ) THEN
        xf0                     = small
      ELSE
        xf0                     = xf(n)
      END IF ! xf(n) == zero 

!-----------------------------------------------------------------------
!  Set Grav to g_force_c
!-----------------------------------------------------------------------

      grav(n)                   = g_force_c(n)

!-----------------------------------------------------------------------
!  Zone-centered fictitious forces
!-----------------------------------------------------------------------

      fict(n)                   = -uf(n) * wf(n)/radius
      sinxf0                    = DSIN(xf0)
      IF ( DABS(sinxf0) > 1.d-3 ) THEN
        fict(n)                 = fict(n) + vf(n) * vf(n)/radius * DCOS(xf(n))/sinxf0
      END IF ! DABS(sinxf0) > 1.d-3
      IF ( xf(n) == zero ) fict(n) = zero

!-----------------------------------------------------------------------
!  Zone-centered neutrino stress
!-----------------------------------------------------------------------

      nuf(n)                    = nu_strc(n)

    END DO ! n = nmin-4, nmax+5

  END IF ! ngeomy == 0
END IF ! sweep == 'y'

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-----------------------------------------------------------------------
!
!                  \\\\\ Z-DIRECTION FORCES /////
!
!-----------------------------------------------------------------------

IF ( sweep == 'z' ) THEN

  pseudo_fict                   = zero

  IF ( ngeomz == 0 ) THEN                            ! Cartesian

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = zero
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomz == 3 ) THEN                       ! Cylindrical

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = - uf(n) * wf(n)/radius
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomz == 5 ) THEN                       ! Spherical

    DO n = nmin-4, nmax+5
      grav(n)                   = zero
      fict(n)                   = - uf(n) * wf(n) * stheta/radius
      IF ( stheta > 1.d-3 ) THEN
        fict(n)                 = fict(n) - uf(n) * vf(n) * ctheta/radius
      END IF ! stheta > 1.d-3
      IF ( xf(n) == zero ) fict(n) = zero

!-----------------------------------------------------------------------
!  Zone-centered neutrino stress
!-----------------------------------------------------------------------

      nuf(n)                    = nu_strc(n)

    END DO ! n = nmin-4, nmax+5

  END IF ! ngeomz == 0
END IF ! sweep == 'z'

RETURN
END SUBROUTINE forces_center
