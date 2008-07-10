SUBROUTINE unpack_eos_keys( c_eos_data, d_eos_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_eos_keys
!    Module:       unpack_eos_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/23/04
!
!    Purpose:
!      To unpack the edit key arrays and restore the values
!       to the appropriate variables in edit_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  c_eos_data : integer array of edit keys
!  d_eos_data : real*8 array of edit keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  eos_snc_x_module, eos_snc_y_module, prb_cntl_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, rhoes, eos_i, eosrho
USE eos_snc_y_module, ONLY : dgrid_y=>dgrid, tgrid_y=>tgrid, ygrid_y=>ygrid, &
& rhoes_y=>rhoes, eos_i_y=>eos_i, eosrho_y=>eosrho
USE eos_snc_z_module, ONLY : dgrid_z=>dgrid, tgrid_z=>tgrid, ygrid_z=>ygrid, &
& rhoes_z=>rhoes, eos_i_z=>eos_i, eosrho_z=>eosrho
USE prb_cntl_module, ONLY : tnse, tdnse
USE radial_ray_module, ONLY : dgrid_r=>dgrid, tgrid_r=>tgrid, ygrid_r=>ygrid, &
& rhoes_r=>rhoes

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER (len=1), INTENT(in), DIMENSION(1)  :: c_eos_data  ! character array of edit keys

REAL(KIND=double), INTENT(in), DIMENSION(14) :: d_eos_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i             ! do i ndex

!-----------------------------------------------------------------------
!
!                   \\\\\ UNPACK EOS KEYS /////
!
!-----------------------------------------------------------------------

eos_i                     = c_eos_data(1)
eos_i_y                   = c_eos_data(1)
eos_i_z                   = c_eos_data(1)

DO i = 1,3
  dgrid(i)                = d_eos_data(0*3+i)
  tgrid(i)                = d_eos_data(1*3+i)
  ygrid(i)                = d_eos_data(2*3+i)
  dgrid_r(i)              = d_eos_data(0*3+i)
  tgrid_r(i)              = d_eos_data(1*3+i)
  ygrid_r(i)              = d_eos_data(2*3+i)
  dgrid_y(i)              = d_eos_data(0*3+i)
  tgrid_y(i)              = d_eos_data(1*3+i)
  ygrid_y(i)              = d_eos_data(2*3+i)
  dgrid_z(i)              = d_eos_data(0*3+i)
  tgrid_z(i)              = d_eos_data(1*3+i)
  ygrid_z(i)              = d_eos_data(2*3+i)
END DO

rhoes(1)                  = d_eos_data(10)
rhoes(2)                  = d_eos_data(11)
rhoes_r(1)                = d_eos_data(10)
rhoes_r(2)                = d_eos_data(11)
rhoes_y(1)                = d_eos_data(10)
rhoes_y(2)                = d_eos_data(11)
rhoes_z(1)                = d_eos_data(10)
rhoes_z(2)                = d_eos_data(11)
eosrho                    = d_eos_data(12)
eosrho_y                  = d_eos_data(12)
eosrho_z                  = d_eos_data(12)
tnse                      = d_eos_data(13)
tdnse                     = d_eos_data(14)

RETURN
END SUBROUTINE unpack_eos_keys
