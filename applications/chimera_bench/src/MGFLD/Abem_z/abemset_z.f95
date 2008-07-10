SUBROUTINE abemset_z( k, ki_ray, kj_ray, rho, t, ye, agr )
!-----------------------------------------------------------------------
!
!    File:         abemset_z
!    Module:       abemset_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/19/06
!
!    Purpose:
!      To direct the recalculation of the emission and absorption inverse
!       mean free paths at the corners of the unit cell for azimuthal zone k.
!      The state rho(k), t(k), and ye(k) for azimuthal zone k is first tested
!       to determine if it still lies within the umit cell for zone k. If
!       so, the program returns. If not, the unit cell surrounding the state
!       is determined, and the emission and absorption inverse mean free
!       paths at the corners of this unit cell are determined and stored.
!
!    Variables that must be passed through common:
!      none
!
!    Subprograms called:
!  abemrgn     : fills the nearrest neighbor table entries of emission and absorption rates
!
!    Input arguments:
!
!  k           : azimuthal zone index
!  ki_ray      : x (radial) index of a specific z (azimuthal) ray
!  kj_ray      : y (angular) index of a specific z (azimuthal) ray
!  rho         : matter density (g/cm**3)
!  t           : matter temperature (K)
!  ye          : electron fraction
!  agr         : lapse function
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!
!  iaefnp      : 0 (1); absorption and emission on free neutrons and
!                 protons set to zero (included)
!  rhoaefnp    : density above which e-neutrino and e-antineutrino
!                 absorption on free neutrons and protons is turned off
!  iaencnu     : 0 (1); absorption and emission on nuclei with cross
!                 sections given by Haxton set to zero (included)
!  roaencnu    : density above which e-neutrino and e-antineutrino
!                 absorption on nuclei (Haxton's rates) is turned off
!  iaence      : 0 (1); absorption and emission of e-neutrinos on nuclei
!                 with cross sections given by Fuller, et al.
!                 set to zero (included)
!  edmpe       : ( en - cmpn ) - ( ep - cmpp) ; mean difference between
!                 the exitation energy of daughter and parent
!                 nucleus for e-neutrino absorption on nuclei
!  iaenca      : 0 (1); absorption and emission of e-antineutrinos
!                 on nuclei set to zero (included)
!  edmpa       : ( en - cmpn ) - ( ep - cmpp) ; mean difference
!                 between the exitation energy of daughter and parent
!                 nucleus for e-antineutrino absorption on nuclei
!
!    Modules used:
!  kind_module, numerical_module
!  abem_module, eos_snc_x_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : one

USE abem_z_module, ONLY : idrae, itrae, iyrae
USE eos_snc_z_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : iaefnp, iaencnu, rhoaefnp, roaencnu


IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: k             ! azimuthal zone index
INTEGER, INTENT(in)              :: kj_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: ki_ray        ! y (angular) index of a specific z (azimuthal) ray

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction
REAL(KIND=double), INTENT(in)    :: agr           ! lapse function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: id            ! density grid index
INTEGER                          :: it            ! temperature grid index
INTEGER                          :: iy            ! electron fraction grid index

INTEGER                          :: idd           ! density do index
INTEGER                          :: itt           ! temperature do index
INTEGER                          :: iyy           ! electron fraction do index

INTEGER                          :: idp           ! density cube index
INTEGER                          :: itp           ! temperature cube index
INTEGER                          :: iyp           ! electron fraction cube index

INTEGER                          :: iaefnpp       ! local value of iaefnp
INTEGER                          :: iaencnup      ! local value of iaencnu

REAL(KIND=double)                :: rhotest       ! value of rho at upper cube corner

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if no neutrinos
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!  Compute independent variable grid indices, (id,it,iy). These define
!   the corner of a "unit cube" in thermodynamic state space surrounding
!   the system point.
!
!   id: In a grid of log(density) of dgrid points per decade,
!
!          id < log(rho)*dgrid < id + 1
!
!   it: In a grid of log(temperature) of dgrid points per decade,
!
!          it < log(t)*tgrid < it + 1
!
!   iy: In a grid of electron fraction of ygrid points between 0.5 and 0.0,
!
!          iy < ( 1.0 - ye*ygrid ) < iy + 1
!-----------------------------------------------------------------------

id                 = INT( dgrid(idty(k,kj_ray,ki_ray)) * DLOG10( rho ) )
it                 = INT( tgrid(idty(k,kj_ray,ki_ray)) * DLOG10( t   ) )
iy                 = INT( ygrid(idty(k,kj_ray,ki_ray)) * ( one - ye ) )

!-----------------------------------------------------------------------
!  Test whether zone is still within local table
!-----------------------------------------------------------------------

IF ( id == idrae(k,kj_ray,ki_ray)  .and.                                &
&    it == itrae(k,kj_ray,ki_ray)  .and.                                &
&    iy == iyrae(k,kj_ray,ki_ray) )                                      RETURN

!-----------------------------------------------------------------------
!  Recompute local table if zone state is no longer inside prior table.
!-----------------------------------------------------------------------

iaefnpp            = iaefnp
iaencnup           = iaencnu
rhotest            = 10.d0**( DBLE(id+1)/dgrid(idty(k,kj_ray,ki_ray)) )

IF ( rhotest > rhoaefnp ) iaefnpp  = 0
iF ( rhotest > roaencnu ) iaencnup = 0

DO idd = id,id+1
  DO itt = it,it+1
    DO iyy = iy,iy+1

      idp          = ( idd - id + 1 )
      itp          = ( itt - it + 1 )
      iyp          = ( iyy - iy + 1 )
      CALL abemrgn_z( k, ki_ray, kj_ray, idd, itt, iyy, idp, itp, iyp, agr )

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Save the indices
!-----------------------------------------------------------------------

idrae(k,kj_ray,ki_ray) = id
itrae(k,kj_ray,ki_ray) = it
iyrae(k,kj_ray,ki_ray) = iy

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE abemset_z
