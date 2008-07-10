SUBROUTINE unpack_radial_ray_keys( c_radhyd_data, i_radhyd_data, d_radhyd_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_radial_ray_keys
!    Module:       unpack_radial_ray_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/02/04
!
!    Purpose:
!      To unpack the radhyd keys defining the dimensions, geometry,
!       and principal parameters of the probkem.
!
!    Subprograms called:
!
!    Input arguments:
!  c_radhyd_data : character array of radhyd_ray keys
!  i_radhyd_data : integer array of radhyd_ray keys
!  d_radhyd_data : 64 bit real array of radhyd_ray keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  radial_ray_module, evh1_zone, prb_cntl_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : ncycle, ncymax, time, t_start, t_bounce, &
& t_stop, tb_stop, dtnph, dtnmh, ndim, ngeomx, ngeomy, ngeomz, nleftx, &
& nlefty, nleftz, nrightx, nrighty, nrightz, lagr, rezn, v_trans_0, imin, &
& imax, jmin, jmax, kmin, kmax, n_lgrgrid, n_eulgrid, n1zoom, n2zoom, &
& n3zoom, r_1, r_2, r_3, m_1, m_2, m_3, zoome1, zoome2, zoome3, m_grid, &
& regrid, t_bounce_lagr_chg, t_bounce_mgrd_chg, grid_frac, int_pre_b, &
& int_post_b, rho_regrid, dtnph_trans, y_shft, dy_shift, ncy_shift, &
& tb_dy_shift, sub_cy_yz, t_step_xyz, xmin, xmax, ymin, ymax, zmin, zmax, &
& nu_equil, rho_equilibrate, t_equilibrate, i_grav, v_diff, G_trns, vel_ns, &
& rot, A, beta

USE evh1_global, ONLY : i_grav_g=>i_grav, v_diff_g=>v_diff
USE evh1_zone, ONLY : imax_e=>imax, jmax_e=>jmax, kmax_e=>kmax
USE prb_cntl_module, ONLY : lagrp=>lagr, reznp=>rezn, v_trans_0p=>v_trans_0, &
& i_gravp=>i_grav
USE t_cntrl_module, ONLY : time_tcntl=>time, t_start_tcntl=>t_start, &
& t_bounce_tcntl=>t_bounce, t_stop_tcntl=>t_stop, tb_stop_tcntl=>tb_stop, &
& dtnph_tcntl=>dtnph, dtnmh_tcntl=>dtnmh, dtnph_trans_t=>dtnph_trans, &
& t_step_xyz_t=>t_step_xyz

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(in), DIMENSION(20) :: c_radhyd_data ! Lagrangian - Eulerian toggle

INTEGER, INTENT(in), DIMENSION(50)           :: i_radhyd_data ! integer array of edit keys

REAL(KIND=double), INTENT(in), DIMENSION(50) :: d_radhyd_data ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!
!                  \\\\\ UNPACK RADHYD KEYS /////
!
!-----------------------------------------------------------------------

lagr                = c_radhyd_data(1)
rezn                = c_radhyd_data(2)
m_grid              = c_radhyd_data(3)
regrid              = c_radhyd_data(4)
y_shft              = c_radhyd_data(5)
v_trans_0           = c_radhyd_data(6)
sub_cy_yz           = c_radhyd_data(7)
t_step_xyz          = c_radhyd_data(8)
nu_equil            = c_radhyd_data(9)
G_trns              = c_radhyd_data(10)
rot                 = c_radhyd_data(11)

ncycle              = i_radhyd_data(1)
ncymax              = i_radhyd_data(2)
ndim                = i_radhyd_data(3)
ngeomx              = i_radhyd_data(4)
ngeomy              = i_radhyd_data(5)
ngeomz              = i_radhyd_data(6)
nleftx              = i_radhyd_data(7)
nrightx             = i_radhyd_data(8)
nlefty              = i_radhyd_data(9)
nrighty             = i_radhyd_data(10)
nleftz              = i_radhyd_data(11)
nrightz             = i_radhyd_data(12)
imin                = i_radhyd_data(13)
imax                = i_radhyd_data(14)
jmin                = i_radhyd_data(15)
jmax                = i_radhyd_data(16)
kmin                = i_radhyd_data(17)
kmax                = i_radhyd_data(18)
n_lgrgrid           = i_radhyd_data(19)
n_eulgrid           = i_radhyd_data(20)
n1zoom              = i_radhyd_data(21)
n2zoom              = i_radhyd_data(22)
n3zoom              = i_radhyd_data(23)
ncy_shift           = i_radhyd_data(24)
i_grav              = i_radhyd_data(25)
int_pre_b           = i_radhyd_data(26)
int_post_b          = i_radhyd_data(27)

time                = d_radhyd_data(1)
t_start             = d_radhyd_data(2)
t_bounce            = d_radhyd_data(3)
t_stop              = d_radhyd_data(4)
tb_stop             = d_radhyd_data(5)
dtnph               = d_radhyd_data(6)
dtnmh               = d_radhyd_data(7)
r_1                 = d_radhyd_data(8)
r_2                 = d_radhyd_data(9)
r_3                 = d_radhyd_data(10)
m_1                 = d_radhyd_data(11)
m_2                 = d_radhyd_data(12)
m_3                 = d_radhyd_data(13)
zoome1              = d_radhyd_data(14)
zoome2              = d_radhyd_data(15)
zoome3              = d_radhyd_data(16)
t_bounce_lagr_chg   = d_radhyd_data(17)
t_bounce_mgrd_chg   = d_radhyd_data(18)
grid_frac           = d_radhyd_data(19)
rho_regrid          = d_radhyd_data(20)
dy_shift            = d_radhyd_data(21)
tb_dy_shift         = d_radhyd_data(22)
v_diff              = d_radhyd_data(23)
xmin                = d_radhyd_data(24)
xmax                = d_radhyd_data(25)
ymin                = d_radhyd_data(26)
ymax                = d_radhyd_data(27)
zmin                = d_radhyd_data(28)
zmax                = d_radhyd_data(29)
rho_equilibrate     = d_radhyd_data(30)
t_equilibrate       = d_radhyd_data(31)
vel_ns              = d_radhyd_data(32)
A                   = d_radhyd_data(33)
beta                = d_radhyd_data(34)

time_tcntl          = d_radhyd_data(1)
t_start_tcntl       = d_radhyd_data(2)
t_bounce_tcntl      = d_radhyd_data(3)
t_stop_tcntl        = d_radhyd_data(4)
tb_stop_tcntl       = d_radhyd_data(5)
dtnph_tcntl         = d_radhyd_data(6)
dtnmh_tcntl         = d_radhyd_data(7)

!-----------------------------------------------------------------------
!  Transfer lagr and rezn to prb_cntl_module
!-----------------------------------------------------------------------
 
lagrp               = lagr
reznp               = rezn
v_trans_0p          = v_trans_0
i_gravp             = i_grav

!-----------------------------------------------------------------------
!  Transfer index limits to evh1_zone_module
!-----------------------------------------------------------------------

imax_e              = imax
jmax_e              = jmax
kmax_e              = kmax

!-----------------------------------------------------------------------
!  Transfer t_step_xyz to t_cntrl_module
!-----------------------------------------------------------------------

t_step_xyz_t        = t_step_xyz

!-----------------------------------------------------------------------
!  Transfer i_grav to evh1_global_module
!-----------------------------------------------------------------------

i_grav_g            = i_grav
v_diff_g            = v_diff

!-----------------------------------------------------------------------
!  Set dtnph_trans
!-----------------------------------------------------------------------

dtnph_trans         = dtnph
dtnph_trans_t       = dtnph


RETURN
END SUBROUTINE unpack_radial_ray_keys
