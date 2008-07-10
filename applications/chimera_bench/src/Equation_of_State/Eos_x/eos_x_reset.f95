SUBROUTINE eos_x_reset( nx, jr_min, jr_max, rho, t, ye, ij_ray, ik_ray, &
& reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         eos_x_reset
!    Module:       eos_x_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/28/00
!
!    Purpose:
!      To reset the thermodymamic and neutrino rate tables.
!
!    Subprograms called:
!  eqstt_x        : interpolates quantities in local EOS table
!  esrgnz_comp_x  : regenerates (qlwaqys) the local EOS table
!  esrgnz_x       : reganerates (if necessary) local EOS table
!  eqstz_x        : interpolates quantities in local EOS table
!  gammaz_x       : computes the EOS gammas
!
!    Input arguments:
!  nx             : x-array extent
!  jr_min         : minimum radial zone for which thermodynamic variables are to be evaluated
!  jr_max         : maximum radial zone for which thermodynamic variables are to be evaluated
!  rho            : shifted matter density array (g/cm**3).
!  t              : shifted matter matter temperature array (K).
!  ye             : shifted matter matter electron fraction array.
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!  reset_comp_eos : composition EOS reset flag
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_x_module, nucbrn_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one

USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY : duesrc, aesv, xn, nuc_number, nse, a_nuc_rep, &
& z_nuc_rep
USE nucbrn_module, ONLY: a_name

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos   ! composition EOS reset flag

LOGICAL                          :: first_c = .true.

INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: jr_min           ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max           ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray           ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho  ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t    ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ye   ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: i                ! composition index
INTEGER                          :: j                ! x-array zone index
INTEGER                          :: jr_maxp          ! jr_max+1

INTEGER                          :: i_n              ! neutron abundance index
INTEGER                          :: i_p              ! proton abundance index
INTEGER                          :: i_4He            ! 4He abundance index

REAL(KIND=double)                :: duddp            ! dudd dummy variable
REAL(KIND=double)                :: dudtp            ! dudt dummy variable
REAL(KIND=double)                :: dudyp            ! dudy dummy variable

REAL(KIND=double), DIMENSION(nx) :: umat             ! internal energy; used for comparing regrids

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

umat                                 = zero

!-----------------------------------------------------------------------
!
!             \\\\\ RESET EOS TABLES /////
!
!-----------------------------------------------------------------------

jr_maxp                              = jr_max + 1

!-----------------------------------------------------------------------
!  Recompute thermodynamic quantities
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), umat(j), &
&  duddp, dudtp, dudyp )
END DO

IF ( reset_comp_eos == 'ye' ) THEN
  CALL esrgnz_comp_x( jr_min, jr_maxp, rho, t, ye, ij_ray, ik_ray )
ELSE
  CALL esrgnz_x( jr_min, jr_maxp, rho, t, ye, ij_ray, ik_ray )
END IF
CALL eqstz_x( jr_min, jr_maxp, rho, t, ye, ij_ray, ik_ray )
CALL gammaz_x( jr_min, jr_maxp, rho, t, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Compute difference in internal energy due to possible EOS regridding
!-----------------------------------------------------------------------

duesrc(jr_min:jr_max,ij_ray, ik_ray) = duesrc(jr_min:jr_max,ij_ray,ik_ray) &
&                                    + aesv(jr_min:jr_max,2,ij_ray,ik_ray) - umat(jr_min:jr_max)

!-----------------------------------------------------------------------
!  Put nse composition in composition arrays
!-----------------------------------------------------------------------

IF ( first_c ) THEN
  first_c                 = .false.
  i_n                     = nuc_number + 1
  i_p                     = nuc_number + 1
  i_4He                   = nuc_number + 1
  DO i = 1,nuc_number
    IF ( a_name(i) == '  n  ' ) THEN
      i_n                 = i
    END IF ! a_name(i) == '  n  '
    IF ( a_name(i) == '  p  ' ) THEN
      i_p                 = i
    END IF ! a_name(i) == '  p  '
    IF ( a_name(i) == '  4He' ) THEN
      i_4He               = i
    END IF ! a_name(i) == '  4He'
  END DO ! i = 1,nuc_number
END IF ! first_c

DO j = jr_min,jr_max
  IF ( nse(j,ij_ray,ik_ray) == 0 ) CYCLE
  xn(j,:)                 = zero
  xn(j,i_n)               = aesv(j,7,ij_ray,ik_ray)
  xn(j,i_p)               = aesv(j,8,ij_ray,ik_ray)
  xn(j,nuc_number+1)      = aesv(j,9,ij_ray,ik_ray)
  xn(j,i_4He)             = DMAX1( one - xn(j,i_n) - xn(j,i_p) - xn(j,nuc_number+1), zero )
  a_nuc_rep(j)            = aesv(j,10,ij_ray,ik_ray)
  z_nuc_rep(j)            = aesv(j,11,ij_ray,ik_ray)
END DO ! j = jr_min,jr_max

RETURN
END SUBROUTINE eos_x_reset
