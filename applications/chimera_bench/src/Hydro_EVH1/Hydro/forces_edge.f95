SUBROUTINE forces_edge( ij_ray, ik_ray, xf, uf, vf, wf, r, grav, fict, &
& nuf )
!-----------------------------------------------------------------------
!
!    File:         forces_edge
!    Module:       forces_edge
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To calculate real (ie, gravity and radiation) and fictitious (i.e., 
!       coriolis and centrifugal) forces at the zone interfaces.
!
!    Input arguments:
!
!  ij_ray           : y-index of an x ray, x-index of an y ray, x index of a z array
!  ik_ray           : z-index of an x ray, z-index of an y ray, y index of a z array
!  xf               : padded zone edged coordinate
!  uf               : x-velocity
!  vf               : y-velocity
!  wf               : z-velocity
!  r                : density
!
!    Output arguments:
!
!  grav             : padded gravitational force
!  fict             : padded fictitious force
!  nuf              : padded neutrino stress
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  e_advct_module, evh1_global, evh1_sweep, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: max_12
USE numerical_module, ONLY : zero, epsilon, third
USE physcnst_module, ONLY : pi, G => g

USE evh1_global, ONLY : ngeomx, ngeomy, ngeomz, nleftx, small, i_radial
USE evh1_sweep, ONLY : sweep, ctheta, stheta, radius, nmin, nmax, dvol0, &
& rhobar, nu_stre, g_force_e, g_pot_e, e_nu_c, f_nu_e
USE parallel_module, ONLY : myid_y

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                  :: ij_ray   ! y-index of an x ray, x-index of an y ray, x index of a z array
INTEGER, INTENT(in)                  :: ik_ray   ! z-index of an x ray, z-index of an y ray, y index of a z array

REAL(KIND=double), DIMENSION(max_12) :: xf       ! zone position
REAL(KIND=double), DIMENSION(max_12) :: uf       ! x-velocity
REAL(KIND=double), DIMENSION(max_12) :: vf       ! y-velocity
REAL(KIND=double), DIMENSION(max_12) :: wf       ! z-velocity
REAL(KIND=double), DIMENSION(max_12) :: r        ! density

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(max_12) :: grav     ! gravitational force
REAL(KIND=double), DIMENSION(max_12) :: fict     ! fictitious force
REAL(KIND=double), DIMENSION(max_12) :: nuf      ! neutrino stress

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                              :: n        ! neutrino flavor index
INTEGER                              :: nsym     ! ghost zone index

REAL(KIND=double), DIMENSION(max_12) :: area     ! area of zone interfaces
REAL(KIND=double), DIMENSION(max_12) :: poteng   ! ptential energy

REAL(KIND=double)                    :: massi    ! mass interior to zone edge
REAL(KIND=double)                    :: massih   ! mass interior to zone center
REAL(KIND=double)                    :: massz    ! mass of shell
REAL(KIND=double)                    :: mass_co  ! central mass
REAL(KIND=double)                    :: xf0      ! minimum radius
REAL(KIND=double)                    :: sinxf0   ! DSIN(xf0)

REAL(KIND=double), PARAMETER         :: fopi = 4.d0 * pi

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

grav                          = zero
fict                          = zero
nuf                           = zero

!-----------------------------------------------------------------------
!
!                  \\\\\ X-DIRECTION FORCES /////
!
!-----------------------------------------------------------------------

IF ( sweep == 'x') THEN

  IF ( ngeomx == 0  .and.  ngeomy == 0 ) then       ! Cartesian

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = zero
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomx == 0  .and.  ngeomy == 1 ) THEN  ! Cylindrical radial

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = wf(n) * wf(n)/xf(n)
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomx == 1  .and.  ngeomy == 3 ) THEN  ! Cylindrical polar

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = vf(n) * vf(n)/xf(n)
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

      massi                   = zero
      grav(nmin-4:nmin-1)     = zero
      fict(nmin-4:nmin-1)     = zero
      poteng(nmin-4:nmin-1)   = zero
      area(nmin-4:nmin-1)     = xf(nmin-4:nmin-1) * xf(nmin-4:nmin-1)
      area(nmin:nmax+5)       = xf(nmin:nmax+5) * xf(nmin:nmax+5)
      nuf                     = zero
      fict                    = zero

!-----------------------------------------------------------------------
!  Compute gravity
!-----------------------------------------------------------------------

      DO n = nmin, nmax
        grav(n)               = g_force_e(n)
        poteng(n)             = g_pot_e(n)
      END DO
      grav  (nmax+1:nmax+5)   = grav  (nmax)
      poteng(nmax+1:nmax+5)   = poteng(nmax)

!-----------------------------------------------------------------------
!  Compute nongravitational forces
!
!  nuf(n)    : zone-edged neutrino stress
!  fict(n)   : fictitious forces
!-----------------------------------------------------------------------

      DO n = nmin, nmax+5
        nuf(n)                = nu_stre(n)
        IF ( xf(n) /= zero ) THEN
          fict(n)             = ( wf(n) * wf(n) + vf(n) * vf(n) )/xf(n)
        END IF
      END DO ! n
             
!-----------------------------------------------------------------------
!  Set force at inner zone edge =0
!-----------------------------------------------------------------------

      fict(nmin)              = zero
             
!-----------------------------------------------------------------------
!  Impose symetric forces on interior ghosts
!-----------------------------------------------------------------------

      DO n = nmin-4,nmin-1
        nsym                  = nmin - n + nmin
        grav(n)               = -grav(nsym)
        fict(n)               = -fict(nsym)
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

      massi                   = mass_co
             
!-----------------------------------------------------------------------
!  Subtract mass of inner ghosts
!-----------------------------------------------------------------------

      DO n = nmin-4,nmin-1
        massi                 = massi - rhobar(n) * dvol0(n) * fopi
      END DO ! n = nmin-4, nmax-1
             
!-----------------------------------------------------------------------
!  Gravity, including inner half of zone
!-----------------------------------------------------------------------

      DO n = nmin-4, nmax+5

        massz                 = rhobar(n) * dvol0(n) * fopi

!-----------------------------------------------------------------------
!  Zone-edge masses
!-----------------------------------------------------------------------

        massih                = massi
             
!-----------------------------------------------------------------------
!  Update mass interior
!-----------------------------------------------------------------------

        massi                 = massi + massz  !  Update mass interior
             
!........Gravity

        grav(n)               = -G * massih/xf(n)**2
             
!........Fictitios Forces

        fict(n)               = zero
        fict(n)               = ( wf(n) * wf(n) + vf(n) * vf(n) )/xf(n)

      END DO

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
      grav(n)                 = zero
      grav(n)                 = - grav(nsym)
      fict(n)                 = zero
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 1 ) THEN                      ! Cylindrical radial

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = wf(n) * wf(n)/xf(n)
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 3  .and.  ngeomx == 1 ) THEN  ! Cylindrical polar (2D)

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = - uf(n) * wf(n)/radius
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 3  .and.  ngeomx == 2 ) THEN  ! Spherical equator (2D)

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = - uf(n) * wf(n)/radius
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomy == 4  .and.  ngeomx == 2 ) THEN  ! Spherical

    nuf                       = zero
    DO n = nmin-4, nmax+5
      IF ( xf(n) == zero ) THEN
        xf0                   = small
      ELSE
        xf0                   = xf(n)
      END IF ! xf(n) = zero

!-----------------------------------------------------------------------
!  Set Grav to zero
!-----------------------------------------------------------------------

      grav(n)                 = g_force_e(n)

!-----------------------------------------------------------------------
!  Zone-centered fictitious forces
!-----------------------------------------------------------------------

      fict(n)                 = -uf(n) * wf(n)/radius
      sinxf0                  = DSIN(xf0)
      IF ( DABS(sinxf0) > 1.d-3 ) THEN
        fict(n)               = fict(n) + vf(n) * vf(n)/radius * DCOS(xf(n))/sinxf0
      END IF ! DABS(sinxf0) > 1.d-3
      IF ( xf(n) == zero ) fict(n) = zero

!-----------------------------------------------------------------------
!  Zone-edge neutrino stress
!-----------------------------------------------------------------------

      nuf(n)                  = nu_stre(n)

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

  IF ( ngeomz == 0 ) THEN                            ! Cartesian

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = zero
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomz == 3 ) THEN                       ! Cylindrical

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = - uf(n) * wf(n) * stheta/radius
    END DO ! n = nmin-4, nmax+5

  ELSE IF ( ngeomz == 5 ) THEN                       ! Spherical

    DO n = nmin-4, nmax+5
      grav(n)                 = zero
      fict(n)                 = - uf(n) * wf(n)/radius
      IF ( stheta > 1.d-3 ) THEN
        fict(n)               = fict(n) - uf(n) * vf(n) * ctheta/radius
      END IF ! stheta > 1.d-3
      IF ( xf(n) == zero ) fict(n) = zero

!-----------------------------------------------------------------------
!  Zone-edge neutrino stress
!-----------------------------------------------------------------------

      nuf(n)                    = nu_stre(n)

    END DO ! n = nmin-4, nmax+5

  END IF ! sweep == 'z'
END IF ! END IF

RETURN
END SUBROUTINE forces_edge
