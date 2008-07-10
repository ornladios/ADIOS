SUBROUTINE esrgnz_comp_z( kmin, kmax, rho, t, ye, ki_ray, kj_ray )
!-----------------------------------------------------------------------
!
!    File:         esrgnz_comp_z
!    Module:       esrgnz_comp_z
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To test each radial zone j between kmin and kmax to determine
!       if its thermodynamic state point is still within the local table
!       of nearest neighbor entries created for that zone in a prior
!       cycle. If not, or if nse = 0, a new local table of nearest
!       neighbor entries is created. Subroutine eosdtgen_z is called to
!       fill this table with eos quantities.
!
!    Input arguments:
!  kmin                 : minimum z (azimuthal) zone for which thermodynamic variables are to be evaluated
!  kmax                 : maximum z (azimuthal) zone for which thermodynamic variables are to be evaluated
!  rho                  : azimuthal array of matter density (g/cm**3)
!  t                    : azimuthal array of matter temperature (K)
!  ye                   : azimuthal array of matter electron fraction
!  ki_ray               : x (radial) index of a specific z (azimuthal) ray
!  kj_ray               : y (angular) index of a specific z (azimuthal) ray
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  dgrid(idty(j,kj_ray,ki_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,kj_ray,ki_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,kj_ray,ki_ray)) : number of table entries in ye between ye = 0.5 and
!                                  ye = 0 for zone j
!  idty(j,kj_ray,ki_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  rhoes(i)                     : density boundary between idty=i and idty=i+1
!  idr(j,kj_ray,ki_ray)         : rho grid index for zone j
!  itr(j,kj_ray,ki_ray)         : t grid index for zone j
!  iyr(j,kj_ray,ki_ray)         : ye grid index for zone j
!
!    Output arguments (common):
!
!  estble(i,j,ida,ita,iya,kj_ray,ki_ray) : equation of state table array
!  idty(j,kj_ray,ki_ray)                 : grid selector
!
!    Subprograms called:
!  eosdtgen_z           : computes EOS quantities
!
!    Include files:
!  kind_module, numerical_module,
!  cycle_module, edit_module, eos_snc_z_module, mdl_cnfg_z_module,
!  parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nz
USE numerical_module, ONLY : zero, one, epsilon

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nlog
USE eos_snc_z_module, ONLY : idty, dgrid, tgrid, ygrid, rhoes, idr, itr, &
& iyr, estble, escnst, nse
USE parallel_module, ONLY : myid

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: kmin            ! minimum z (azimuthal) zone index
INTEGER, INTENT(in)              :: kmax            ! maximum z (azimuthal) zone index
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray

REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rho ! azimuthal array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: t   ! azimuthal array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: ye  ! azimuthal array of matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! thermodynamic variable index
INTEGER                          :: j               ! radial zone index

INTEGER                          :: id              ! density grid index
INTEGER                          :: it              ! temperature grid index
INTEGER                          :: iy              ! electron fraction grid index

INTEGER                          :: idd             ! density do index
INTEGER                          :: itt             ! temperature do index
INTEGER                          :: iyy             ! electron fraction do index

INTEGER                          :: idp             ! density cube index
INTEGER                          :: itp             ! temperature cube index
INTEGER                          :: iyp             ! electron fraction cube index

INTEGER                          :: sf              ! eos flag

REAL(KIND=double)                :: ye_min          ! minimim upper cube corner value of ye
REAL(KIND=double)                :: yed             ! max( ye, ye_min )

REAL(KIND=double)                :: fd              ! location of rho in rho-t-ye grid
REAL(KIND=double)                :: ft              ! location of t in rho-t-ye grid
REAL(KIND=double)                :: fy              ! location of ye in rho-t-ye grid

REAL(KIND=double), DIMENSION(2,2) :: estmn
REAL(KIND=double), DIMENSION(2,2) :: estmx

  101 FORMAT (' Table generation unsuccessful in esrgnz_comp_z, ncycle=',i7, &
& ' myid=',i4,' ki_ray=',i4,' kj_ray=',i4,' j=',i4,' rho, t, ye=',3(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Begin angular loop
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO j = kmin,kmax

!-----------------------------------------------------------------------
!  Determine which grid to use, then store in the grid index array
!   idty(j,i_ray).
!-----------------------------------------------------------------------

  IF ( rho(j) < rhoes(1) )               THEN
    idty(j,kj_ray,ki_ray)             = 1
  ELSE IF ( rho(j) < rhoes(2) )          THEN
    idty(j,kj_ray,ki_ray)             = 2
  ELSE
    idty(j,kj_ray,ki_ray)             = 3
  END IF

!-----------------------------------------------------------------------
!  The value of ye_min is set so that the lower cube corner value
!   of ye is at least one grid point away from zero.
!-----------------------------------------------------------------------

  ye_min                              = 2.d0/ygrid(idty(j,kj_ray,ki_ray))
  yed                                 = DMAX1( ye(j), ye_min )

!-----------------------------------------------------------------------
!  Compute independent variable grid indices
!-----------------------------------------------------------------------

  fd                                  = dgrid(idty(j,kj_ray,ki_ray)) * dlog10(rho(j))
  ft                                  = tgrid(idty(j,kj_ray,ki_ray)) * dlog10(t(j)  )
  fy                                  = ygrid(idty(j,kj_ray,ki_ray)) * (one - yed   )

  id                                  = INT(fd)
  it                                  = INT(ft)
  iy                                  = INT(fy)

!-----------------------------------------------------------------------
!  Test whether the thermodynamic state of angular zone j is still
!   within its local table of nearest neighbor entries.
!-----------------------------------------------------------------------

  IF ( id == idr(j,kj_ray,ki_ray)  .and.   &
&      it == itr(j,kj_ray,ki_ray)  .and.   &
&      iy == iyr(j,kj_ray,ki_ray)  .and.   &
&      nse(j,kj_ray,ki_ray) == 1               ) CYCLE

!-----------------------------------------------------------------------
!  Recompute local table if the thermodynamic state of angular zone j is
!   no longer inside its local table or if nse = 0.
!-----------------------------------------------------------------------

  DO idd = id,id+1
    DO itt = it,it+1
      DO iyy = iy,iy+1

        idp                           = ( idd - id + 1 )
        itp                           = ( itt - it + 1 )
        iyp                           = ( iyy - iy + 1 )
        CALL eosdtgen_z( j, ki_ray, kj_ray, idd, itt, iyy, idp, itp, iyp, sf )

      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!  If sf = 0. table generation was unsuccessful, use old table entries
!   and print message.
!-----------------------------------------------------------------------

  IF ( sf == 0 ) THEN
    WRITE (nprint,101) ncycle, myid, j, ki_ray, kj_ray, rho(j), t(j), ye(j)
    WRITE (nlog,101) ncycle, myid, j, ki_ray, kj_ray, rho(j), t(j), ye(j)
    CYCLE
  END IF

!-----------------------------------------------------------------------
!  Save the thermodymaic state indices for angular zone j
!-----------------------------------------------------------------------

  idr(j,kj_ray,ki_ray)              = id
  itr(j,kj_ray,ki_ray)              = it
  iyr(j,kj_ray,ki_ray)              = iy

!-----------------------------------------------------------------------
!  Compute escnst(i,j,kj_ray,ki_ray), which is added to the thermodynamic
!   table entries for the ith quantity of radial zone j to ensure that
!   all entries are positive so that logarithms can be taken. This 
!   quantity will be subtracted from the interpolated value of the ith
!   thermodynmaic quantity to  restore its original value.
!-----------------------------------------------------------------------

  DO i = 1,12
    escnst(i,j,kj_ray,ki_ray)       = zero
    DO id = 1,2
      DO it = 1,2
        DO iy = 1,2
          escnst(i,j,kj_ray,ki_ray) = DMIN1( escnst(i,j,kj_ray,ki_ray), DBLE( estble(i,j,id,it,iy,kj_ray,ki_ray) ) )
        END DO
      END DO
    END DO
    escnst(i,j,kj_ray,ki_ray)       = -2.d+00 * escnst(i,j,kj_ray,ki_ray)
  END DO

!-----------------------------------------------------------------------
!  Store logarithms of table entries
!-----------------------------------------------------------------------

  DO id = 1,2
    DO it = 1,2
      DO iy = 1,2
        DO i = 1,12
          estble(i,j,id,it,iy,kj_ray,ki_ray)                            &
&                                   = dlog10( estble(i,j,id,it,iy,kj_ray,ki_ray) &
&                                   + escnst(i,j,kj_ray,ki_ray) + epsilon )
        END DO
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!  Rearrange for monotonicity of energy
!-----------------------------------------------------------------------

  DO id = 1,2
    DO iy = 1,2
      estmn(id,iy)                  = AMIN1( estble(2,j,id,1,iy,kj_ray,ki_ray), &
&                                     estble(2,j,id,2,iy,kj_ray,ki_ray) )
      estmx(id,iy)                  = AMAX1( estble(2,j,id,1,iy,kj_ray,ki_ray), &
&                                     estble(2,j,id,2,iy,kj_ray,ki_ray) )
      estble(2,j,id,1,iy,kj_ray,ki_ray) = estmn(id,iy)
      estble(2,j,id,2,iy,kj_ray,ki_ray) = estmx(id,iy)
    END DO
  END DO

!-----------------------------------------------------------------------
!  End radial loop
!-----------------------------------------------------------------------

END DO

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE esrgnz_comp_z
