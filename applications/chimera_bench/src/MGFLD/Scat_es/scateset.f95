SUBROUTINE scateset( j, ij_ray, ik_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         scateset
!    Module:       scateset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/29/03
!
!    Purpose:
!      To test whether the rho-T-ye state point has moved out of
!       the unit rho-T-ye cell, or if the unit cube has not been set
!       up. If so, subroutine scatergn is called to recompute (or
!       compute for the first time) the neutrino-electron scattering 
!       functions at corners of the unit cube.
!
!    Subprograms called:
!      scatergn
!
!    Input arguments:
!
!  j              : radial zone number
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!  rho            : density (g/cm**3)
!  t              : temperature (K)
!  ye             : electron fraction
! 
!    Output arguments:
!      none
!
!    Input arguments (common):
!
!  nes            : 0 - neutrino-electron scattering turned off; nes
!                    subroutines are bypassed; nes scattering function
!                    arrays, if used, must be zeroed elsewhere
!                 : 1 - neutrino-electron scattering turned on
!  rhonesmn       : density below which n-neutrino-electron scattering is
!                    turned off (nes scatering function arrays are zeroed
!  rhonesmx       : density above which n-neutrino-electron scattering is
!                    turned  off (nes scatering function arrays are zeroed
!  dgrid(idty(j,i_ray))
!                 : number of table entries per decade in rho for zone j
!  tgrid(idty(j,i_ray))
!                 : number of table entries per decade in t for zone j
!  ygrid(idty(j,i_ray))
!                 : number of table entries in ye between ye = 0.5 and
!                    ye = 0 for zone j
!  idty(j,i_ray)  : index for dgrid, tgrid, and ygrid for zone j
!  idrse(j,i_ray) : rho grid index for zone j
!  itrse(j,i_ray) : t grid index for zone j
!  iyrse(j,i_ray) : ye grid index for zone j
!
!    Include files:
!  kind_modul, array_module, numerical_module
!  eos_snc_modul, nu_energy_grid_module, prb_cntl.cmn, scat_e.cmn
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : ij_ray_dim
USE numerical_module, ONLY : one

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid,idty
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : nes
USE scat_e_module, ONLY : idrse, itrse, iyrse, scte0i, scte0ii

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

INTEGER                          :: k             ! incomiong neutrino energy zone index
INTEGER                          :: kp            ! outgoing neutrino energy zone index

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

INTEGER                          :: id            ! density grid index
INTEGER                          :: it            ! temperature grid index
INTEGER                          :: iy            ! electron fraction grid index
INTEGER                          :: idtr_test     ! cumulative difference between current and proir table indices

INTEGER                          :: idd           ! density do index
INTEGER                          :: itt           ! temperature do index
INTEGER                          :: iyy           ! electron fraction do index

INTEGER                          :: idp           ! density cube index
INTEGER                          :: itp           ! temperature cube index
INTEGER                          :: iyp           ! electron fraction cube index

INTEGER                          :: ida           ! transposed density cube index
INTEGER                          :: ita           ! transposed temperature cube index
INTEGER                          :: iya           ! transposed electron fraction cube index

INTEGER, DIMENSION(2,2,2)        :: id_trans      ! old density index to be transposes
INTEGER, DIMENSION(2,2,2)        :: it_trans      ! old temperature index to be transposes
INTEGER, DIMENSION(2,2,2)        :: iy_trans      ! old electron fraction index to be transposes

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------        
!  No neutrino-electron scattering if
!     nes = 0
!   or if
!     nnugp(n) = 0 for all n
!-----------------------------------------------------------------------

IF ( nes == 0  .or.  nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!  Compute independent variable grid indices, (id,it,iy). These define
!   the corner of a "unit cube" in thermodynamic state space surrounding
!   the system point.
!
!  id : In a grid of log(density) of dgrid points per decade,
!
!     id < log(rho)*dgrid < id + 1
!
!  it : In a grid of log(temperature) of dgrid points per decade,
!
!     it < log(t)*tgrid < it + 1
!
!  iy : In a grid of electron fraction of ygrid points between
!        0.5 and 0.0,
!
!     iy < ( 1.0 - ye*ygrid ) < iy + 1
!-----------------------------------------------------------------------

id                 = INT( dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho ) )
it                 = INT( tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   ) )
iy                 = INT( ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye ) )

!-----------------------------------------------------------------------
!  Test whether zone is still within local table
!-----------------------------------------------------------------------
idtr_test          = ABS( id - idrse(j,ij_ray,ik_ray) )                  &
&                  + ABS( it - itrse(j,ij_ray,ik_ray) )                  &
&                  + ABS( iy - iyrse(j,ij_ray,ik_ray) )

IF ( idtr_test == 0 ) THEN

  RETURN

!-----------------------------------------------------------------------
!  Recompute entire local table.
!-----------------------------------------------------------------------

ELSE IF ( idtr_test > 1 ) THEN

  DO idd = id,id+1
    DO itt = it,it+1
      DO iyy = iy,iy+1

        idp        = ( idd - id + 1 )
        itp        = ( itt - it + 1 )
        iyp        = ( iyy - iy + 1 )
        CALL scatergn( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, iyp )

      END DO ! iyy = iy,iy+1
    END DO ! itt = it,it+1
  END DO ! idd = id,id+1

!-----------------------------------------------------------------------
!  Transpose one side of the table to the other, recompute the first
!   side.
!-----------------------------------------------------------------------

ELSE

  CALL index_transpose( id, it, iy, idrse(j,ij_ray,ik_ray), itrse(j,ij_ray,ik_ray), &
&  iyrse(j,ij_ray,ik_ray), id_trans, it_trans, iy_trans )

!-----------------------------------------------------------------------
!  Transpose one side of the table to the other.
!  Combine indices to uniquely index the radial ray
!-----------------------------------------------------------------------

  i_ray            = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

  DO idd = id,id+1
    DO itt = it,it+1
      DO iyy = iy,iy+1

        idp        = ( idd - id + 1 )
        itp        = ( itt - it + 1 )
        iyp        = ( iyy - iy + 1 )
        
        ida        = id_trans(idp,itp,iyp)
        IF ( ida == 0 ) CYCLE
        ita        = it_trans(idp,itp,iyp)
        iya        = iy_trans(idp,itp,iyp)
        DO k = 1,nnugpmx
          DO kp = 1,k
            scte0i (j,k,kp,idp,itp,iyp,i_ray) = scte0i (j,k,kp,ida,ita,iya,i_ray)
            scte0i (j,kp,k,idp,itp,iyp,i_ray) = scte0i (j,kp,k,ida,ita,iya,i_ray)
            scte0ii(j,k,kp,idp,itp,iyp,i_ray) = scte0ii(j,k,kp,ida,ita,iya,i_ray)
            scte0ii(j,kp,k,idp,itp,iyp,i_ray) = scte0ii(j,kp,k,ida,ita,iya,i_ray)
          END DO ! kp = 1,k
        END DO ! k = 1,nnugpmx

      END DO ! iyy = iy,iy+1
    END DO ! itt = it,it+1
  END DO ! idd = id,id+1

  DO idd = id,id+1
    DO itt = it,it+1
      DO iyy = iy,iy+1

        idp        = ( idd - id + 1 )
        itp        = ( itt - it + 1 )
        iyp        = ( iyy - iy + 1 )
        
        ida        = id_trans(idp,itp,iyp)
        IF ( ida == 0 ) THEN
          CALL scatergn( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, iyp )
        END IF ! ida == 0

      END DO ! iyy = iy,iy+1
    END DO ! itt = it,it+1
  END DO ! idd = id,id+1

END IF ! idtr_test == 0

IF ( id == idrse(j,ij_ray,ik_ray)  .and.  &
&    it == itrse(j,ij_ray,ik_ray)  .and.  &
&    iy == iyrse(j,ij_ray,ik_ray) ) RETURN

!-----------------------------------------------------------------------
!  Save the indices
!-----------------------------------------------------------------------

idrse(j,ij_ray,ik_ray) = id
itrse(j,ij_ray,ik_ray) = it
iyrse(j,ij_ray,ik_ray) = iy

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE scateset
