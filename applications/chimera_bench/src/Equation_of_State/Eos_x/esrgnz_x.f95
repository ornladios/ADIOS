SUBROUTINE esrgnz_x( jr_min, jr_max, rho, t, ye, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         esrgnz_x
!    Module:       esrgnz_x
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To test each radial zone j between jr_min and jr_max to
!       determine if its thermodynamic state point is still within
!       the local table of nearest neighbor entries created for
!       that zone in a prior cycle. If not, a new local table of
!       nearest neighbor entries is created. Subroutine eosdtgen_x
!       is called to fill this table with eos quantities.
!
!    Input arguments:
!  jr_min                       : minimum radial zone for which thermodynamic 
!                                  variables are to be evaluated
!  jr_max                       : minimum radial zone for which thermodynamic
!                                  variables are to be evaluated
!  rho                          : shifted matter density array (g/cm**3).
!  t                            : shifted matter matter temperature array (K).
!  ye                           : shifted matter matter electron fraction array.
!  ij_ray                       : j-index of a radial ray
!  ik_ray                       : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  rho(j)                       : matter density (g/cm**3) of radial zone j
!  t(j)                         : matter temperature (K) of radial zone j
!  ye(j)                        : matter electron fraction of radial zone j
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in
!   rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in
!   t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between
!   ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  rhoes(i)                     : density boundary between idty=i and idty=i+1
!  idr(j,ij_ray,ik_ray)         : rho grid index for zone j
!  itr(j,ij_ray,ik_ray)         : t grid index for zone j
!  iyr(j,ij_ray,ik_ray)         : ye grid index for zone j
!
!    Output arguments (common):
!
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray) : equation of state table array
!  idty(j,ij_ray,ik_ray)                 : grid selector
!
!    Subprograms called:
!  eosdtgen_x           : computes EOS quantities
!
!    Include files:
!  kind_module, array_module, numerical_module,
!  cycle_module, edit_module, eos_snc_x_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx
USE numerical_module, ONLY : zero, one, epsilon

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nlog
USE eos_snc_x_module, ONLY : idty, dgrid, tgrid, ygrid, rhoes, idr, itr, &
& iyr, estble, escnst
USE parallel_module, ONLY : myid

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t   ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ye  ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i             ! thermodynamic variable index
INTEGER                          :: j             ! radial zone index

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

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Table generation unsuccessful in esrgnz_x, ncycle=',i7, &
& ' myid=',i4,' j=',i4,' ij_ray=',i4,' ij_ray=',i4,' rho, t, ye=',3(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Begin radial loop
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Determine which grid to use, then store in the grid index array
!   idty(j,ij_ray,ik_ray).
!-----------------------------------------------------------------------

  IF ( rho(j) < rhoes(1) ) THEN
    idty(j,ij_ray,ik_ray) = 1
  ELSE IF ( rho(j) < rhoes(2) ) THEN
    idty(j,ij_ray,ik_ray) = 2
  ELSE
    idty(j,ij_ray,ik_ray) = 3
  END IF

!-----------------------------------------------------------------------
!  The value of ye_min is set so that the lower cube corner value
!   of ye is at least one grid point away from zero.
!-----------------------------------------------------------------------

  ye_min                  = 2.d0/ygrid(idty(j,ij_ray,ik_ray))
  yed                     = DMAX1( ye(j), ye_min )

!-----------------------------------------------------------------------
!  Compute independent variable grid indices
!-----------------------------------------------------------------------

  fd                      = dgrid(idty(j,ij_ray,ik_ray)) * dlog10(rho(j))
  ft                      = tgrid(idty(j,ij_ray,ik_ray)) * dlog10(t(j)  )
  fy                      = ygrid(idty(j,ij_ray,ik_ray)) * (one - yed   )

  id                      = IDINT(fd)
  it                      = IDINT(ft)
  iy                      = IDINT(fy)

!-----------------------------------------------------------------------
!  Test whether the thermodynamic state of radial zone j is still
!   within its local table of nearest neighbor entries.
!-----------------------------------------------------------------------

  IF ( id == idr(j,ij_ray,ik_ray)  .and.   &
&      it == itr(j,ij_ray,ik_ray)  .and.   &
&      iy == iyr(j,ij_ray,ik_ray)              ) CYCLE

!-----------------------------------------------------------------------
!  Recompute local table if the thermodynamic state of radial zone j is
!   no longer inside its local table.
!-----------------------------------------------------------------------

  DO idd = id,id+1
    DO itt = it,it+1
      DO iyy = iy,iy+1

        idp               = ( idd - id + 1 )
        itp               = ( itt - it + 1 )
        iyp               = ( iyy - iy + 1 )
        CALL eosdtgen_x( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, &
&        iyp, sf )

      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!  If sf = 0. table generation was unsuccessful, use old table
!   entries and print message.
!-----------------------------------------------------------------------

  IF ( sf == 0 ) THEN
    WRITE (nprint,101) ncycle, myid, j, ij_ray, ik_ray, rho(j), t(j), ye(j)
    WRITE (nlog,101) ncycle, myid, j, ij_ray, ik_ray, rho(j), t(j), ye(j)
    CYCLE
  END IF

!-----------------------------------------------------------------------
!  Save the thermodymaic state indices for radial zone j
!-----------------------------------------------------------------------

  idr(j,ij_ray,ik_ray)    = id
  itr(j,ij_ray,ik_ray)    = it
  iyr(j,ij_ray,ik_ray)    = iy

!-----------------------------------------------------------------------
!  Compute escnst(i,j,ij_ray,ik_ray), which is added to the thermodynamic
!   table entries for the ith quantity of radial zone j to ensure that
!   all entries are positive so that logarithms can be taken. This
!   quantity will be subtracted from the interpolated value of the ith
!   thermodynmaic quantity to restore its original value.
!-----------------------------------------------------------------------

  DO i = 1,12
    escnst(i,j,ij_ray,ik_ray) = zero
    DO id = 1,2
      DO it = 1,2
        DO iy = 1,2
          escnst(i,j,ij_ray,ik_ray)                                        &
&                             = DMIN1( escnst(i,j,ij_ray,ik_ray),          &
&                               DBLE( estble(i,j,id,it,iy,ij_ray,ik_ray) ) )
        END DO
      END DO
    END DO
    escnst(i,j,ij_ray,ik_ray) = -2.d+00 * escnst(i,j,ij_ray,ik_ray)
  END DO

!-----------------------------------------------------------------------
!  Store logarithms of table entries
!-----------------------------------------------------------------------

  DO id = 1,2
    DO it = 1,2
      DO iy = 1,2
        DO i = 1,12
          estble(i,j,id,it,iy,ij_ray,ik_ray)                               &
&                             = DLOG10( estble(i,j,id,it,iy,ij_ray,ik_ray) &
&                             + escnst(i,j,ij_ray,ik_ray) + epsilon )
        END DO
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!  Rearrange for monotonicity of energy
!-----------------------------------------------------------------------

  DO id = 1,2
    DO iy = 1,2
      estmn(id,iy)            = AMIN1( estble(2,j,id,1,iy,ij_ray,ik_ray), &
&                                      estble(2,j,id,2,iy,ij_ray,ik_ray) )
      estmx(id,iy)            = AMAX1( estble(2,j,id,1,iy,ij_ray,ik_ray), &
&                                      estble(2,j,id,2,iy,ij_ray,ik_ray) )
      estble(2,j,id,1,iy,ij_ray,ik_ray) = estmn(id,iy)
      estble(2,j,id,2,iy,ij_ray,ik_ray) = estmx(id,iy)
    END DO
  END DO

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  End radial loop
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END DO !  j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE esrgnz_x
