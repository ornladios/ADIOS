SUBROUTINE esrgn_comp_z( k, ki_ray, kj_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         esrgn_comp_z
!    Module:       esrgn_comp_z
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To reevaluate the eos quantities at the corners of the same eos cube.
!       Subroutine eosdtgen_z is called to fill this table with eos quantities.
!
!    Input arguments:
!
!  k                    : z (azimuthal) zone index
!  ki_ray               : x (radial) index of a specific z (azimuthaal) ray
!  kj_ray               : y (angular) index of a specific z (azimuthaal) ray
!  rho                  : matter density (g/cm**3).
!  t                    : matter temperature (K).
!  ye                   : matter electron fraction.
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  dgrid(idty(k,kj_ray,ki_ray)) : number of table entries per decade in rho for zone k
!  tgrid(idty(k,kj_ray,ki_ray)) : number of table entries per decade in t for zone k
!  ygrid(idty(k,kj_ray,ki_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone k
!  idty(k,kj_ray,ki_ray)        : index for dgrid, tgrid, and ygrid for zone k
!  rhoes(i)                     : density boundary between idty=i and idty=i+1
!  idr(k,kj_ray,ki_ray)         : rho grid index for zone k
!  itr(k,kj_ray,ki_ray)         : t grid index for zone k
!  iyr(k,kj_ray,ki_ray)         : ye grid index for zone k
!
!    Output arguments (common):
!
!  estble(i,k,idp,itp,iyp,kj_ray,ki_ray) : equation of state table array
!
!    Subprograms called:
!  eosdtgen_z           : computes the EOS quantities
!
!    Include files:
!  kind_module, numerical_module,
!  cycle_module, edit_module, eos_snc_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, epsilon

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nlog
USE eos_snc_z_module, ONLY : idr, itr, iyr, estble, escnst

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: k             ! z (azimuthal) zone index
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthaal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthaal) ray

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

  101 FORMAT (' Table generation unsuccessful in esrgn_comp_y, ncycle=',i7, &
& ' k=',i4,' ki_ray=',i3,' kj_ray=',i4,' rho, t, ye=',3(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Get independent variable grid indices
!-----------------------------------------------------------------------

id                          = idr(k,kj_ray,ki_ray)
it                          = itr(k,kj_ray,ki_ray)
iy                          = iyr(k,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Recompute local table.
!-----------------------------------------------------------------------

DO idd = id,id+1
  DO itt = it,it+1
    DO iyy = iy,iy+1

      idp                   = ( idd - id + 1 )
      itp                   = ( itt - it + 1 )
      iyp                   = ( iyy - iy + 1 )
      CALL eosdtgen_z( k, ki_ray, kj_ray, idd, itt, iyy, idp, itp, iyp, sf )

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  If sf = 0. table generation was unsuccessful, use old table  entries
!  and print message.
!-----------------------------------------------------------------------

IF ( sf == 0 ) THEN
  WRITE (nprint,101) ncycle, k, ki_ray, kj_ray, rho, t, ye
  WRITE (nlog,101) ncycle, k, ki_ray, kj_ray, rho, t, ye
  RETURN
END IF

!-----------------------------------------------------------------------
!  Save the thermodymaic state indices for radial zone k
!-----------------------------------------------------------------------
 
idr(k,kj_ray,ki_ray)        = id
itr(k,kj_ray,ki_ray)        = it
iyr(k,kj_ray,ki_ray)        = iy

!-----------------------------------------------------------------------
!  Compute escnst(i,k), which is added to the thermodynamic table
!   entries for the ith quantity of radial zone k to ensure that all
!   entries are positive so that logarithms can be taken. This quantity
!   will be subtracted from the interpolated value of the ith
!   thermodynmaic quantity to restore its original value.
!-----------------------------------------------------------------------

DO i = 1,12
  escnst(i,k,kj_ray,ki_ray) = zero
  DO id = 1,2
    DO it = 1,2
      DO iy = 1,2
        escnst(i,k,kj_ray,ki_ray) = DMIN1( escnst(i,k,kj_ray,ki_ray),   &
&                                   DBLE( estble(i,k,id,it,iy,kj_ray,ki_ray) ) )
      END DO
    END DO
  END DO
  escnst(i,k,kj_ray,ki_ray) = -2.d+00 * escnst(i,k,kj_ray,ki_ray)
END DO

!-----------------------------------------------------------------------
!  Store logarithms of table entries
!-----------------------------------------------------------------------

DO id = 1,2
  DO it = 1,2
    DO iy = 1,2
      DO i = 1,12
        estble(i,k,id,it,iy,kj_ray,ki_ray) = DLOG10( estble(i,k,id,it,iy,kj_ray,ki_ray) &
&                           + escnst(i,k,kj_ray,ki_ray) + epsilon )
      END DO
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!  Rearrange for monotonicity of energy
!-----------------------------------------------------------------------

DO id = 1,2
  DO iy = 1,2
    estmn(id,iy)            = AMIN1( estble(2,k,id,1,iy,kj_ray,ki_ray), &
&                             estble(2,k,id,2,iy,kj_ray,ki_ray) )
    estmx(id,iy)            = AMAX1( estble(2,k,id,1,iy,kj_ray,ki_ray), &
&                             estble(2,k,id,2,iy,kj_ray,ki_ray) )
    estble(2,k,id,1,iy,kj_ray,ki_ray) = estmn(id,iy)
    estble(2,k,id,2,iy,kj_ray,ki_ray) = estmx(id,iy)
  END DO
END DO

!-----------------------------------------------------------------------
!  Table generation complete
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE esrgn_comp_z
