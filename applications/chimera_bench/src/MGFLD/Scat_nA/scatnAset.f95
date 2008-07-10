SUBROUTINE scatnAset( j, ij_ray, ik_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         scatnAset
!    Module:       scatnAset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/10/03
!
!    Purpose:
!      To test whether the rho-T-ye state point has moved out of
!       the unit rho-T-ye cell, or if the unit cube has not been set
!       up. If so, subroutine scatargn is called to recompute (or
!       compute for the first time) the neutrino-nucleus inelastic
!       scattering  functions at corners of the unit cube.
!
!    Subprograms called:
!      scatnArgn
!
!    Input arguments:
!
!  j                      : radial zone number
!  ij_ray                 : j-index of a radial ray
!  ik_ray                 : k-index of a radial ray
!  rho                    : density (g/cm**3)
!  t                      : temperature (K)
!  ye                     : electron fraction
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!
!  isctnA                 : 0 - neutrino-nucleus inelastic scattering turned off;
!                            sctnA subroutines are bypassed; sctnA scattering
!                            function arrays, if used, must be zeroed elsewhere
!                         : 1 - neutrino-nucleus inelastic scattering turned on
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and
!                          ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)  : index for dgrid, tgrid, and ygrid for zone j
!  idrnA(j,ij_ray,ik_ray) : rho grid index for zone j
!  itrnA(j,ij_ray,ik_ray) : t grid index for zone j
!  iyrnA(j,ij_ray,ik_ray) : ye grid index for zone j
!
!    Include files:
!  kind_module, array_module, numerical_module
!  eos_snc_modul, nu_energy_grid_module, prb_cntl_module, scat_a_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu
USE numerical_module, ONLY : one

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : isctnA
USE scat_nA_module, ONLY : idrnA, itrnA, iyrnA

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
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: nnugpmx       ! maximum neutrino group number

INTEGER                          :: idd           ! density do index
INTEGER                          :: itt           ! temperature do index
INTEGER                          :: iyy           ! electron fraction do index

INTEGER                          :: idp           ! density cube index
INTEGER                          :: itp           ! temperature cube index
INTEGER                          :: iyp           ! electron fraction cube index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------        
!  No neutrino-nucleus inelastic scattering if
!     iscta = 0
!   or if
!     nnugp(n) = 0 for all n
!-----------------------------------------------------------------------

nnugpmx            = 0
DO n = 1,nnu
  nnugpmx          = MAX( nnugpmx, nnugp(n) )
END DO

IF ( isctnA == 0  .or.  nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!  Compute independent variable grid indices, (id,it,iy). These define
!   the corner of a "unit cube" in thermodynamic state space surrounding
!   the system point.
!
!  id : In a grid of log(density) of dgrid points per decade,
!
!          id < log(rho)*dgrid < id + 1
!
!  it : In a grid of log(temperature) of dgrid points per decade,
!
!          it < log(t)*tgrid < it + 1
!
!  iy : In a grid of electron fraction of ygrid points between
!        0.5 and 0.0,
!
!          iy < ( 1.0 - ye*ygrid ) < iy + 1
!-----------------------------------------------------------------------

id                 = INT( dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho ) )
it                 = INT( tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   ) )
iy                 = INT( ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye ) )

!-----------------------------------------------------------------------
!  Test whether zone is still within local table
!-----------------------------------------------------------------------

IF ( id == idrnA(j,ij_ray,ik_ray)  .and.  &
&    it == itrnA(j,ij_ray,ik_ray)  .and.  &
&    iy == iyrnA(j,ij_ray,ik_ray) ) RETURN

!-----------------------------------------------------------------------
!  Recompute local table if zone state is no longer inside prior table.
!-----------------------------------------------------------------------

DO idd = id,id+1
  DO itt = it,it+1
    DO iyy = iy,iy+1

      idp          = ( idd - id + 1 )
      itp          = ( itt - it + 1 )
      iyp          = ( iyy - iy + 1 )
      CALL scatnArgn( j, ij_ray,ik_ray, idd, itt, iyy, idp, itp, iyp )

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Save the indices
!-----------------------------------------------------------------------

idrnA(j,ij_ray,ik_ray) = id
itrnA(j,ij_ray,ik_ray) = it
iyrnA(j,ij_ray,ik_ray) = iy

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE scatnAset
