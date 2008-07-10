SUBROUTINE unpack_hydro_keys( nx, i_hydro_data, d_hydro_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_hydro_keys
!    Module:       unpack_hydro_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the hydro keys defining themode of hydrodynamics,
!       tolerances, and initial time step and hydro time step controls.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp       : unit number from which to read
!  nprint       : unit number from which to print
!  nx           : radial array dimension
!
!    Output arguments:
!  i_hydro_data : integer array of hydro keys
!  d_hydro_data : real*8 array of hydro keys
!
!    Include files:
!  kind_module
!  bomb_module, boundary_module, convect_module, evh1_global,
!  prb_cntl_module, shock_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE bomb_module, ONLY : e_bomb, bomb_time, t_start_bomb, jexpl_min, jexpl_max
USE boundary_module, ONLY : iubcjmn, ubcjmn, iubcjmx, ubcjmx, iuset, uset, r_u, &
& ipbnd, pbound
USE convect_module, ONLY : adrftcv, alphlmix
USE evh1_global, ONLY : degen
USE prb_cntl_module, ONLY : ihydro, irelhy, ilapsehy, ilcnvct, iemix, iyemix, &
& ipsimix, ipcnvct, iscnvct, i_bomb
USE shock_module, ONLY : ipq, q0_x, q0_y
USE t_cntrl_module, ONLY : idtj, jtv, jtvshk, ncyshk, ncyrho, rhojtv, tcntrl, &
& rdtmax, dtst1, dtst2, dtst3, ttst1

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension

INTEGER, INTENT(in), DIMENSION(30)                :: i_hydro_data  ! integer array of transport keys

REAL(KIND=double), INTENT(in), DIMENSION((30+nx)) :: d_hydro_data  ! 64 bit real array of transport keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i             ! do i ndex

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!
!                  \\\\\ UNPACK HYDRO KEYS /////
!
!-----------------------------------------------------------------------

ihydro                   = i_hydro_data(1)
irelhy                   = i_hydro_data(2)
ilapsehy                 = i_hydro_data(3)
iubcjmn                  = i_hydro_data(4)
iubcjmx                  = i_hydro_data(5)
iuset                    = i_hydro_data(6)
ipbnd                    = i_hydro_data(7)
idtj                     = i_hydro_data(8)
jtv                      = i_hydro_data(9)
jtvshk                   = i_hydro_data(10)
ncyshk                   = i_hydro_data(11)
ncyrho                   = i_hydro_data(12)
ipq                      = i_hydro_data(13)
ilcnvct                  = i_hydro_data(14)
iemix                    = i_hydro_data(15)
iyemix                   = i_hydro_data(16)
ipsimix                  = i_hydro_data(17)
ipcnvct                  = i_hydro_data(18)
iscnvct                  = i_hydro_data(19)
i_bomb                   = i_hydro_data(20)
jexpl_min                = i_hydro_data(21)
jexpl_max                = i_hydro_data(22)

ubcjmn                   = d_hydro_data(1)
ubcjmx                   = d_hydro_data(2)
r_u                      = d_hydro_data(3)
pbound                   = d_hydro_data(4)
rhojtv                   = d_hydro_data(6)
tcntrl(1)                = d_hydro_data(7)
tcntrl(2)                = d_hydro_data(8)
tcntrl(3)                = d_hydro_data(9)
tcntrl(4)                = d_hydro_data(10)
tcntrl(5)                = d_hydro_data(11)
tcntrl(6)                = d_hydro_data(12)
tcntrl(7)                = d_hydro_data(13)
tcntrl(8)                = d_hydro_data(14)
tcntrl(9)                = d_hydro_data(15)
tcntrl(9)                = d_hydro_data(15)
tcntrl(10)               = d_hydro_data(16)
dtst1                    = d_hydro_data(17)
dtst2                    = d_hydro_data(18)
dtst3                    = d_hydro_data(19)
ttst1                    = d_hydro_data(20)
q0_x(1)                  = d_hydro_data(21)
q0_y(1)                  = d_hydro_data(21)
adrftcv                  = d_hydro_data(22)
alphlmix                 = d_hydro_data(23)
e_bomb                   = d_hydro_data(24)
bomb_time                = d_hydro_data(25)
t_start_bomb             = d_hydro_data(26)
degen                    = d_hydro_data(27)
tcntrl(31)               = d_hydro_data(28)

DO i = 1,nx
  uset(i)                = d_hydro_data(30+i)
END DO


END SUBROUTINE unpack_hydro_keys
