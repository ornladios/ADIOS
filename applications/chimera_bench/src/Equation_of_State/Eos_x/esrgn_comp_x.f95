SUBROUTINE esrgn_comp_x( j, ij_ray, ik_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         esrgn_comp_x
!    Module:       esrgn_comp_x
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!       To reevaluate the eos quantities at the corners of the same eos cube,
!        i.e., without testing whether a new local cube is needed.
!       Subroutine eosdtgen_x is called to fill this table with eos quantities.
!
!    Input arguments:
!
!  j                   : radial zone index.
!  ij_ray              : j-index of a radial ray
!  ik_ray              : k-index of a radial ray
!  rho                 : matter density (g/cm**3).
!  t                   : matter temperature (K).
!  ye                  : matter electron fraction.
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  rhoes(i)                     : density boundary between idty=i and idty=i+1
!  idr(j,ij_ray,ik_ray)         : rho grid index for zone j
!  itr(j,ij_ray,ik_ray)         : t grid index for zone j
!  iyr(j,ij_ray,ik_ray)         : ye grid index for zone j
!
!    Output arguments (common):
!
!   estble(i,j,idp,itp,iyp,ij_ray,ik_ray) : equation of state table array
!
!    Subprograms called:
!  eosdtgen_x          : computes the EOS quantities
!
!    Include files:
!  kind_module, numerical_module,
!  cycle_module, edit_module, eos_snc_x_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, epsilon

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint
USE eos_snc_x_module, ONLY : idr, itr, iyr, estble, escnst

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

INTEGER                          :: i             ! thermodynamic variable index

INTEGER                          :: id            ! density grid index
INTEGER                          :: it            ! temperature grid index
INTEGER                          :: iy            ! electron fraction grid index

INTEGER                          :: idd           ! density do index
INTEGER                          :: itt           ! temperature do index
INTEGER                          :: iyy           ! electron fraction do index

INTEGER                          :: idp           ! density cube index
INTEGER                          :: itp           ! temperature cube index
INTEGER                          :: iyp           ! electron fraction cube index

INTEGER                          :: sf            ! eos flag

REAL(KIND=double), DIMENSION(2,2) :: estmn
REAL(KIND=double), DIMENSION(2,2) :: estmx

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Table generation unsuccessful in esrgn_comp_x, ncycle=',i7, &
& ' j=',i4,' ij_ray=',i4,' ik_ray=',i3,' rho, t, ye=',3(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Get independent variable grid indices
!-----------------------------------------------------------------------

id                          = idr(j,ij_ray,ik_ray)
it                          = itr(j,ij_ray,ik_ray)
iy                          = iyr(j,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Recompute local table.
!-----------------------------------------------------------------------

DO idd = id,id+1
  DO itt = it,it+1
    DO iyy = iy,iy+1

      idp                   = ( idd - id + 1 )
      itp                   = ( itt - it + 1 )
      iyp                   = ( iyy - iy + 1 )
      CALL eosdtgen_x( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, iyp, sf )

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  If sf = 0. table generation was unsuccessful, use old table entries
!   and print message.
!-----------------------------------------------------------------------

IF ( sf == 0 ) THEN
  WRITE (nprint,101) ncycle, j, ij_ray, ik_ray, rho, t, ye
  RETURN
END IF

!-----------------------------------------------------------------------
!  Save the thermodymaic state indices for radial zone j
!-----------------------------------------------------------------------

idr(j,ij_ray,ik_ray)        = id
itr(j,ij_ray,ik_ray)        = it
iyr(j,ij_ray,ik_ray)        = iy

!-----------------------------------------------------------------------
!  Compute escnst(i,j), which is added to the thermodynamic table
!   entries for the ith quantity of radial zone j to ensure that all
!   entries are positive so that logarithms can be taken. This quantity
!   will be subtracted from the interpolated value of the ith
!   thermodynmaic quantity to restore its original value.
!-----------------------------------------------------------------------

DO i = 1,12
  escnst(i,j,ij_ray,ik_ray) = zero
  DO id = 1,2
    DO it = 1,2
      DO iy = 1,2
        escnst(i,j,ij_ray,ik_ray)                                        &
&                           = DMIN1( escnst(i,j,ij_ray,ik_ray),          &
&                             DBLE( estble(i,j,id,it,iy,ij_ray,ik_ray) ) )
      END DO
    END DO
  END DO
  escnst(i,j,ij_ray,ik_ray) = -2.d+00 * escnst(i,j,ij_ray,ik_ray)
END DO

!-----------------------------------------------------------------------
!  Store logarithms of table entries.
!-----------------------------------------------------------------------

DO id = 1,2
  DO it = 1,2
    DO iy = 1,2
      DO i = 1,12
        estble(i,j,id,it,iy,ij_ray,ik_ray)                               &
&                           = DLOG10( estble(i,j,id,it,iy,ij_ray,ik_ray) &
&                           + escnst(i,j,ij_ray,ik_ray) + epsilon )
      END DO
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Rearrange for monotonicity of energy.
!-----------------------------------------------------------------------

DO id = 1,2
  DO iy = 1,2
    estmn(id,iy)            = AMIN1( estble(2,j,id,1,iy,ij_ray,ik_ray),  &
&                                    estble(2,j,id,2,iy,ij_ray,ik_ray) )
    estmx(id,iy)            = AMAX1( estble(2,j,id,1,iy,ij_ray,ik_ray),  &
&                                    estble(2,j,id,2,iy,ij_ray,ik_ray) )
    estble(2,j,id,1,iy,ij_ray, ik_ray) = estmn(id,iy)
    estble(2,j,id,2,iy,ij_ray, ik_ray) = estmx(id,iy)
  END DO
END DO

!-----------------------------------------------------------------------
!  Table generation complete.
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE esrgn_comp_x
