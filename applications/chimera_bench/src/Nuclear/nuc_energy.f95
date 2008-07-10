SUBROUTINE nuc_energy( j, enb_mev, enm )
!===============================================================================
!  This routine finds moments of the abundance distribution useful for
!  hydrodynamics, including the total abundance, electron fraction, binding 
!  energy, and rest mass excess energy and outputs ergs/g. 
!===============================================================================

USE kind_module, ONLY : double
USE array_module, ONLY : nnc
USE numerical_module, ONLY : zero, one, epsilon
USE physcnst_module, ONLY : ergmev, rmu

USE eos_snc_x_module, ONLY : xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, nuc_number, &
& a_nuc, z_nuc, be_nuc

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: j        ! radial zone index

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)    :: enb_mev  ! mean binding energy per particle
REAL(KIND=double), INTENT(out)    :: enm      ! rest mass energy

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                           :: l        ! nuclear abundance index
INTEGER                           :: n_nucp1  ! nuc_number + 1

REAL(KIND=double), PARAMETER      :: mex_p=7.28899d0, mex_n=8.07144d0
REAL(KIND=double), DIMENSION(nnc) :: y
REAL(KIND=double)                 :: ztot     ! electron fraction

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

n_nucp1            = nuc_number + 1

be_nuc(n_nucp1)    = be_nuc_rep(j)
a_nuc(n_nucp1)     = a_nuc_rep(j)
z_nuc(n_nucp1)     = z_nuc_rep(j)
enb_mev            = zero
ztot               = zero
      
DO l = 1,n_nucp1
  y(l)             = xn(j,l)/( a_nuc(l) + epsilon )
  ztot             = ztot + z_nuc(l) * y(l)
  enb_mev          = enb_mev  + be_nuc(l) * y(l)
END DO

enm                = mex_p * ztot + mex_n * ( one - ztot ) - enb_mev

!-----------------------------------------------------------------------
!  Change units from MeV/nucleon to erg/g
!-----------------------------------------------------------------------

enm                = ergmev * enm/rmu

RETURN
END SUBROUTINE nuc_energy
