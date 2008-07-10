SUBROUTINE pairset( j, ij_ray, ik_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         pairset
!    Module:       pairset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/5/03
!
!    Purpose:
!      To test whether the rho-T-ye state point has moved out of the unit
!       rho-T-ye cell, or if the unit cube has not been set up. If so,
!       subroutine pairrgn is called to recompute (or compute for the first
!       time) the pair annihilation functions at corners of the unit cube.
!
!    Subprograms called:
!  pairrgn              : sets up the calculation of the pair annihilation
!                          kernels at the cube corners
!
!    Input arguments:
!
!  j                    : radial zone number
!  ij_ray               : j-index of a radial ray
!  ik_ray               : k-index of a radial ray
!  rho                  : density (g/cm**3)
!  t                    : temperature (K)
!  ye                   : electron fraction
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!
!  ipair                : 0 - pair annihilation turned off; pair annihilation
!                          subroutines are bypassed;  pair annihilation function
!                          arrays, if used, must be zeroed elsewhere
!                       : 1 - pair annihilaiton turned on
!  rhopairemn           : density below which n_e - n_ebar pair annihilation is
!                          not computed (pair annihilation functions are zeroed).
!  rhopairemn           : density above which n_e - n_ebar pair annihilation is
!                          not computed (pair annihilation functions are zeroed).
!  rhopairxmn           : density below which n_x - n_xbar pair annihilation is
!                          not computed (pair annihilation functions are zeroed).
!  rhopairxmn           : density above which n_x - n_xbar pair annihilation is
!                          not computed (pair annihilation functions are zeroed).
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  idrpp(j,ij_ray,ik_ray)       : rho grid index for zone j
!  itrpp(j,ij_ray,ik_ray)       : t grid index for zone j
!  iyrpp(j,ij_ray,ik_ray)       : ye grid index for zone j
!
!    Include files:
!  kind_modul, array_module, numerical_module
!  eos_snc_x_module, nu_energy_grid_module, pair_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY : ij_ray_dim
USE numerical_module, ONLY : one

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid,idty
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE pair_module, ONLY : idrpp, itrpp, iyrpp, paira0i, paira0ii
USE prb_cntl_module, ONLY : ipair

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
!  No pair-annihilation if
!     ipair = 0
!   or if
!     nnugp(n) = 0 for all n
!-----------------------------------------------------------------------

IF ( ipair == 0  .or.  nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!  Compute independent variable grid indices, (id,it,iy). These define
!   the corner of a "unit cube" in thermodynamic state space
!   surrounding the system point.
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
idtr_test          = ABS( id - idrpp(j,ij_ray,ik_ray) )                 &
&                  + ABS( it - itrpp(j,ij_ray,ik_ray) )                 &
&                  + ABS( iy - iyrpp(j,ij_ray,ik_ray) )

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
        CALL pairrgn( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, iyp )

      END DO ! iyy = iy,iy+1
    END DO ! itt = it,it+1
  END DO ! idd = id,id+1

!-----------------------------------------------------------------------
!  Transpose one side of the table to the other, recompute the first
!   side.
!-----------------------------------------------------------------------

ELSE

  CALL index_transpose( id, it, iy, idrpp(j,ij_ray,ik_ray), itrpp(j,ij_ray,ik_ray), &
&  iyrpp(j,ij_ray,ik_ray), id_trans, it_trans, iy_trans )

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
          DO kp = 1,nnugpmx
            paira0i (j,k,kp,idp,itp,iyp,i_ray) = paira0i (j,k,kp,ida,ita,iya,i_ray)
            paira0ii(j,k,kp,idp,itp,iyp,i_ray) = paira0ii(j,k,kp,ida,ita,iya,i_ray)
          END DO ! kp = 1,nnugpmx
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
          CALL pairrgn( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, iyp )
        END IF ! ida == 0

      END DO ! iyy = iy,iy+1
    END DO ! itt = it,it+1
  END DO ! idd = id,id+1

END IF ! idtr_test == 0

IF ( id == idrpp(j,ij_ray,ik_ray)  .and.  &
&    it == itrpp(j,ij_ray,ik_ray)  .and.  &
&    iy == iyrpp(j,ij_ray,ik_ray) ) RETURN

!-----------------------------------------------------------------------
!  Save the indices
!-----------------------------------------------------------------------

idrpp(j,ij_ray,ik_ray) = id
itrpp(j,ij_ray,ik_ray) = it
iyrpp(j,ij_ray,ik_ray) = iy

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE pairset
