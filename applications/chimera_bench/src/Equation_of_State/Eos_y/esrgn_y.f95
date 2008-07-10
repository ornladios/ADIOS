SUBROUTINE esrgn_y( j, ji_ray, jk_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         esrgn_y
!    Module:       esrgn_y
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/18/05
!
!    Purpose:
!      To test radial zone j to determine if its thermodynamic
!       state point is still within the local table of nearest
!       neighbor entries created for that zone in a prior cycle.
!       If not, a new local table of nearest neighbor entries is
!       created. Subroutine eosdtgen_y is called to fill this table
!       with eos quantities.
!
!    Input arguments:
!
!  j                    : y (angular) zone index
!  ji_ray               : x (radial) index of a specific angular ray
!  jk_ray               : z (azimuthal) index of a specific angular ray
!  rho                  : matter density (g/cm**3).
!  t                    : matter temperature (K).
!  ye                   : matter electron fraction.
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  dgrid(idty(j,ji_ray,jk_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ji_ray,jk_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ji_ray,jk_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ji_ray,jk_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  rhoes(i)                     : density boundary between idty=i and idty=i+1
!  idr(j,ji_ray,jk_ray)         : rho grid index for zone j
!  itr(j,ji_ray,jk_ray)         : t grid index for zone j
!  iyr(j,ji_ray,jk_ray)         : ye grid index for zone j
! 
!    Output arguments (common):
!
!   estble(i,j,idp,itp,iyp,ji_ray,jk_ray) : equation of state table array
!
!    Subprograms called:
!  eosdtgen_y           : computes the EOS quantities
!
!    Include files:
!  kind_module, numerical_module,
!  cycle_module, edit_module, eos_snc_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, epsilon

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nlog
USE eos_snc_y_module, ONLY : idty, dgrid, tgrid, ygrid, rhoes, idr, itr, &
& iyr, estble, escnst

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! y (angular) zone index
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific angular ray

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

REAL(KIND=double)                :: ye_min        ! minimim upper cube corner value of ye
REAL(KIND=double)                :: yed           ! max( ye, ye_min )

REAL(KIND=double)                :: fd            ! location of rho in rho-t-ye grid
REAL(KIND=double)                :: ft            ! location of t in rho-t-ye grid
REAL(KIND=double)                :: fy            ! location of ye in rho-t-ye grid

REAL(KIND=double), DIMENSION(2,2) :: estmn
REAL(KIND=double), DIMENSION(2,2) :: estmx

  101 FORMAT (' Table generation unsuccessful in esrgn_y, ncycle=',i7, &
& ' j=',i4,' ji_ray,jk_ray=',i4,' rho, t, ye=',3(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine which grid to use, then store in the grid index array
!   idty(j,ji_ray,jk_ray).
!-----------------------------------------------------------------------

IF ( rho < rhoes(1) ) THEN
  idty(j,ji_ray,jk_ray)           = 1
ELSE IF ( rho < rhoes(2) ) THEN
  idty(j,ji_ray,jk_ray)           = 2
ELSE
  idty(j,ji_ray,jk_ray)           = 3
END IF

!-----------------------------------------------------------------------
!  The value of ye_min is set so that the lower cube corner value of
!   ye is at least one grid point away from zero.
!-----------------------------------------------------------------------

ye_min                            = 2.d0/ygrid(idty(j,ji_ray,jk_ray))
yed                               = DMAX1( ye, ye_min )

!-----------------------------------------------------------------------
!  Compute independent variable grid indices.
!-----------------------------------------------------------------------

fd                                = dgrid(idty(j,ji_ray,jk_ray)) * dlog10(rho)
ft                                = tgrid(idty(j,ji_ray,jk_ray)) * dlog10(t  )
fy                                = ygrid(idty(j,ji_ray,jk_ray)) * (one - yed)

id                                = INT(fd)
it                                = INT(ft)
iy                                = INT(fy)

!-----------------------------------------------------------------------
!  Test whether radial the thermodynamic state of radial zone j is
!   still within its local table of nearest neighbor entries.
!-----------------------------------------------------------------------

IF ( id == idr(j,ji_ray,jk_ray)  .and.      &
&    it == itr(j,ji_ray,jk_ray)  .and.      &
&    iy == iyr(j,ji_ray,jk_ray) )             RETURN

!-----------------------------------------------------------------------
!  Recompute local table if the thermodynamic state of radial zone j
!   is no longer inside its local table.
!-----------------------------------------------------------------------

DO idd = id,id+1
  DO itt = it,it+1
    DO iyy = iy,iy+1

      idp                         = ( idd - id + 1 )
      itp                         = ( itt - it + 1 )
      iyp                         = ( iyy - iy + 1 )
      CALL eosdtgen_y( j, ji_ray,jk_ray, idd, itt, iyy, idp, itp, iyp, sf )

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  If sf = 0. table generation was unsuccessful, use old table  entries
!   and print message.
!-----------------------------------------------------------------------

IF ( sf == 0 ) THEN
  WRITE (nprint,101) ncycle, j, ji_ray, jk_ray, rho, t, ye
  WRITE (nlog,101) ncycle, j, ji_ray, jk_ray, rho, t, ye
  RETURN
END IF

!-----------------------------------------------------------------------
!  Save the thermodymaic state indices for radial zone j
!-----------------------------------------------------------------------

idr(j,ji_ray,jk_ray)              = id
itr(j,ji_ray,jk_ray)              = it
iyr(j,ji_ray,jk_ray)              = iy

!-----------------------------------------------------------------------
!  Compute escnst(i,j), which is added to the thermodynamic table
!   entries for the ith quantity of radial zone j to ensure that all
!   entries are positive so that logarithms can be taken. This quantity
!   will be subtracted from the interpolated value of the ith
!   thermodynmaic quantity to restore its original value.
!-----------------------------------------------------------------------

DO i = 1,12
  escnst(i,j,ji_ray,jk_ray)       = zero
  DO id = 1,2
    DO it = 1,2
      DO iy = 1,2
        escnst(i,j,ji_ray,jk_ray) = DMIN1( escnst(i,j,ji_ray,jk_ray),   &
&                                   DBLE( estble(i,j,id,it,iy,ji_ray,jk_ray) ) )
      END DO
    END DO
  END DO
  escnst(i,j,ji_ray,jk_ray)       = -2.d+00 * escnst(i,j,ji_ray,jk_ray)
END DO

!-----------------------------------------------------------------------
!  Store logarithms of table entries
!-----------------------------------------------------------------------

DO id = 1,2
  DO it = 1,2
    DO iy = 1,2
      DO i = 1,12
        estble(i,j,id,it,iy,ji_ray,jk_ray)                              &
&                                 = DLOG10( estble(i,j,id,it,iy,ji_ray,jk_ray) &
&                                 + escnst(i,j,ji_ray,jk_ray) + epsilon )
      END DO
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Rearrange for monotonicity of energy
!-----------------------------------------------------------------------

DO id = 1,2
  DO iy = 1,2
    estmn(id,iy)                  = AMIN1( estble(2,j,id,1,iy,ji_ray,jk_ray), &
&                                   estble(2,j,id,2,iy,ji_ray,jk_ray) )
    estmx(id,iy)                  = AMAX1( estble(2,j,id,1,iy,ji_ray,jk_ray), &
&                                   estble(2,j,id,2,iy,ji_ray,jk_ray) )
    estble(2,j,id,1,iy,ji_ray,jk_ray) = estmn(id,iy)
    estble(2,j,id,2,iy,ji_ray,jk_ray) = estmx(id,iy)
  END DO
END DO

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE esrgn_y
