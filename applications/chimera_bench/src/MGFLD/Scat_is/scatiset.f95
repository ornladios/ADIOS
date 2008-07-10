SUBROUTINE scatiset( j, ij_ray, ik_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         scatiset
!    Module:       scatiset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To direct the calculation or recalculation of the isoenergetic
!       scattering inverse mean free paths at the corners of the unit
!       cell for radial shell j.
!
!    Subprograms called:
!        scatirgn
!
!    Input arguments:
!
!  j                    : radial zone number
!  ij_ray               : j-index of a radial ray
!  ik_ray               : k-index of a radial ray
!  in                   : 0, neutrino-neutron isoenergetic scattering omitted
!                         1, neutrino-neutron isoenergetic scattering included
!  ip                   : 0, neutrino-proton isoenergetic scattering omitted
!                         1, neutrino-proton isoenergetic scattering included
!  ihe                  : 0, neutrino-helium isoenergetic scattering omitted
!                         1, neutrino-helium isoenergetic scattering included
!  iheavy               : 0, neutrino-heavy nucleus isoenergetic scattering omitted
!                         1, neutrino-heavy nucleus isoenergetic scattering included
!  iscat                : 0, all neutrino scattering processes omitted
!                         1, neutrino scattering processes not necessarily omitted
!  idd                  : density grid index
!  itt                  : temperature grid index
!  iyy                  : electron fraction grid index
!  ida                  : absorption-emission function table density index
!  ita                  : absorption-emission function table temperature index
!  iya                  : absorption-emission function table electron fraction index
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  nnugp(n)             : number of energy zones for neutrinos of type n
!
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and
!                          ye = 0 for zone j
!
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  idrsi(j,ij_ray,ik_ray)       : rho grid index for zone j
!  itrsi(j,ij_ray,ik_ray)       : t grid index for zone j
!  iyrsi(j,ij_ray,ik_ray)       : ye grid index for zone j
!
!    Include files:
!  kind_module, numerical_module
!  eos_snc_x_module, nu_energy_grid_module, prb_cntl_module, scat_i_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : one

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid,idty
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : iscat, in, ip, ihe, iheavy
USE scat_i_module, ONLY : idrsi, itrsi, iyrsi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction

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

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  No isoenergetic scattering
!     if iscat eq 0
!   or if
!     in = 0  and  ip = 0  and  ihe = 0  and  iheavy = 0
!   or if
!     nnugp(n) = 0 for all n
!-----------------------------------------------------------------------

IF ( iscat == 0  .or.  nnugpmx == 0 ) RETURN
IF ( in == 0  .and.  ip == 0  .and.  ihe == 0  .and.  iheavy == 0 ) RETURN


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

id                 = INT( dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho ) )
it                 = INT( tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   ) )
iy                 = INT( ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye ) )

!----------------------------------------------------------------------!
!  Test whether zone is still within local table                 !
!----------------------------------------------------------------------!

IF ( id == idrsi(j,ij_ray,ik_ray)  .and.  &
&    it == itrsi(j,ij_ray,ik_ray)  .and.  &
&    iy == iyrsi(j,ij_ray,ik_ray) ) RETURN

!-----------------------------------------------------------------------
!  Recompute local table if zone state is no longer inside prior table.
!-----------------------------------------------------------------------

DO idd = id,id+1
  DO itt = it,it+1
    DO iyy = iy,iy+1

      idp          = ( idd - id + 1 )
      itp          = ( itt - it + 1 )
      iyp          = ( iyy - iy + 1 )
      CALL scatirgn( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, iyp )

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Save the indices
!-----------------------------------------------------------------------

idrsi(j,ij_ray,ik_ray) = id
itrsi(j,ij_ray,ik_ray) = it
iyrsi(j,ij_ray,ik_ray) = iy

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE scatiset
