SUBROUTINE edit_Global_composition( imin, imax, nx, nnc, time, t_tb, &
& x_c, xn_a_ray, xn_mean, a_mean, nse_min, nse_max, c_shock, nprint )
!-----------------------------------------------------------------------
!
!    File:         edit_Global_composition
!    Module:       edit_Global_composition
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/2/00
!
!    Purpose:
!      To edit the minimum, maximum, and mean values of quantities as a
!      function of the radius.
!
!    Subprograms called:
!  date_and_time_print : fetches and prints the date and time
!
!    Input arguments:
!  imin         : minimum x-array index for the edit
!  imax         : maximum x-array index for the edit
!  nx           : x-array extent
!  nnc          : composition array extent
!  time         : elapsed time
!  t_tb         : time from core bounce
!  x_c          : radial midpoint of zone (cm)
!  xn_a_ray     : composition nasses (solar masses)
!  xn_mean      : ngular average of the composition mass fractions
!  a_mean       : angular average of the mean nuclear mass number along an angular ray
!  nse_min      : minimum value of NSE-nonNSE boundary
!  nse_max      : maximum value of NSE-nonNSE boundary
!  c_shock      : shock location
!  nprint       : unit numberto print data
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_head, nucbrn_module, shock_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half
USE physcnst_module, ONLY : msolar

USE edit_module, ONLY : head
USE nucbrn_module, ONLY : a_name

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                              :: imin           ! minimum unshifted radial zone index
INTEGER, INTENT(in)                              :: imax           ! maximum unshifted radial zone index
INTEGER, INTENT(in)                              :: nx             ! x-array extent
INTEGER, INTENT(in)                              :: nnc            ! composition array extent
INTEGER, INTENT(in)                              :: nse_min        ! minimum value of NSE-nonNSE boundary
INTEGER, INTENT(in)                              :: nse_max        ! maximum value of NSE-nonNSE boundary
INTEGER, INTENT(in)                              :: nprint         ! unit numberto print data

CHARACTER (len=2), INTENT(in), DIMENSION(nx)     :: c_shock        ! shock location

REAL(KIND=double), INTENT(in)                    :: time           ! elapsed time
REAL(KIND=double), INTENT(in)                    :: t_tb           ! time from core bounce
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: x_c            ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc) :: xn_a_ray       ! composition nasses (solar masses)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc) :: xn_mean        ! angular average of the composition mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: a_mean         ! angular average of the mean nuclear mass number along an angular ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=2), DIMENSION(nx)                 :: c_nse

INTEGER                                          :: i             ! do index

REAL(KIND=double), DIMENSION(nx,nnc)             :: xn_a_ray_sum  ! composition nasses sumed from outer edge

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1a)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (9x,'Elapsed time=',es14.7,10x,' Time from bounce=',es14.7/)
    9 FORMAT (41x,'Angular Average of the Composition Mass Fractions')
   11 FORMAT (40x,51('-')/)
   19 FORMAT (41x,'Composition Masses per Zone')
   21 FORMAT (40x,29('-')/)
   29 FORMAT (41x,'Composition Masses Summed from the Surface')
   31 FORMAT (40x,44('-')/)
  101 FORMAT ('   j           r     ',17(3x,a5,3x),'   a_mean'/)
  103 FORMAT (1x,i4,1x,a2,es11.3,a2,18(es11.3))
  111 FORMAT ('   j           r     ',17(3x,a5,3x)/)
  113 FORMAT (1x,i4,1x,a2,es11.3,a2,17(es11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                          \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

DO i = imin,imax
  IF ( i < nse_min  .and.  i < nse_max ) THEN
    c_nse(i)            = '  '
  ELSE IF ( i < nse_max ) THEN
    c_nse(i)            = '* '
  ELSE
    c_nse(i)            = '**'
  END IF
END DO ! i = imax, imin, -1

!-----------------------------------------------------------------------
!  Angular average of the composition mass fractions 
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,7) time,t_tb
WRITE (nprint,9)
WRITE (nprint,11)
WRITE (nprint,101) a_name

DO i = imax, imin, -1
  WRITE (nprint,103) i, c_nse(i), x_c(i), c_shock(i), xn_mean(i,:), a_mean(i)
END DO ! i = imax, imin, -1

!-----------------------------------------------------------------------
!  Composition masses per radial zone summed over angular zones
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,7) time,t_tb
WRITE (nprint,19)
WRITE (nprint,21)
WRITE (nprint,111) a_name

DO i = imax, imin, -1
  WRITE (nprint,113) i, c_nse(i), x_c(i), c_shock(i), xn_a_ray(i,:)
END DO ! i = imax, imin, -1

!-----------------------------------------------------------------------
!  Composition masses summed from the surface
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,7) time,t_tb
WRITE (nprint,29)
WRITE (nprint,31)
WRITE (nprint,111) a_name

xn_a_ray_sum            = zero
DO i = imax,imin,-1
  xn_a_ray_sum(i,1:nnc) = xn_a_ray_sum(i+1,1:nnc) + xn_a_ray(i,1:nnc)
END DO ! i = imax, imin, -1

DO i = imax, imin, -1
  WRITE (nprint,113) i, c_nse(i), x_c(i), c_shock(i), xn_a_ray_sum(i,:)
END DO ! i

RETURN
END SUBROUTINE edit_Global_composition
