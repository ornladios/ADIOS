
!-----------------------------------------------------------------------
!         Computational cycles
!
!      ncycle: the cycle number of the calculation.
!
!      ncymax: problem termination criterion. Calculation terminated when ncycle > ncymax.
!-----------------------------------------------------------------------

cycle                        0                                          ncycle
cycle                  9000000                                          ncymax

!-----------------------------------------------------------------------
!         Times
!
!      time: the elapsed time since the initiation of the calculation (i.e., since ncycle = 0)
!
!      tstart: the time at the initiation of the calculation (typically tstart = 0.0.)
!
!      t_stop: the elapsed time at which the simulation is stopped.
!-----------------------------------------------------------------------

sptime                              0.00000000e+00                      time
sptime                              0.00000000e+00                      t_start
sptime                              0.00000000e+00                      t_bounce
sptime                              1.00000000E+03                      t_stop
sptime                              5.00000000E-01                      tb_stop

!-----------------------------------------------------------------------
!         Time steps
!
!      dtnph: the coordinate 'hydro' time step, set by hydro and nuclear reactions (if the material
!       is not in nse), between time cycle m and time cycle m + 1.
!
!      dtnmh:- the coordinate 'hydro' time step between time cycle m - 1 and time cycle m.
!-----------------------------------------------------------------------

sptime                              5.00000000e-06                      dtnph
sptime                              5.00000000e-06                      dtnmh

!-----------------------------------------------------------------------
!      lagr = ye: Lagrangian hydrodynamics.
!      lagr = no: Eulerian hydrodynamics.
!-----------------------------------------------------------------------

lagr                        no                                          lagr
lagr                                1.00000000E-03                      tb_lagr

!-----------------------------------------------------------------------
!      m_grid = ye: Moving radial grid if lagr = no.
!      m_grid = no: Eulerian if lagr = no.
!-----------------------------------------------------------------------

m_grid                      ye                                          m_grid
m_grid                              1.00000000E-03                      tb_mgrd

!-----------------------------------------------------------------------
!  regrid     : regrid grid switch
!
!     re_grid = ye : regrid every dt_regrid seconds beginning grid_frac
!                    seconds after bounce
!     re_grid = no : regrid option off
!
!  int_pre_b  : number of cycles between successive regrids before the
!                central density reaches 10^{14} g cm^{-3}
!  int_post_b : number of cycles between successive regrids after the
!                central density reaches 10^{14} g cm^{-3}
!  grid_frac  : fraction of a grid width the grid is allowed to move per
!                time step
!  rho_regrid : regrid up to density rho_regrid or 2 zones behind shock,
!                whichever is larger [g cm^{-3}]
!-----------------------------------------------------------------------

regrid                      ye                                          regrid

regrid                      25                                          int_pre_b
regrid                      25                                          int_post_b
regrid                              1.00000000E-01                      grid_frac
regrid                              1.00000000E+06                      rho_regrid

!-----------------------------------------------------------------------
!      rot  : imposed rotation switch
!
!      rot = ye : impart to the initial model a rotation with constant
!       angular velocity on cylinders according to the rotation law
!
!                                _      _ -1
!                               |      2 |
!                               |     r  |
!          Omega(r) = Omega_{0} | 1 + -- |
!                               |      2 |
!                               |_    A _|
!
!       where Omega(r) is the angular velocity, r the distance from the
!       rotation axis. A and beta (the ratio of the rotational energy to
!       the gravitational binding energy are free parameters that determine
!       the  rotational energy of the model and the distribution of angular
!       momentum. Omega_{0} is iterated until beta achieves the value
!       selected.
!
!     rot = no : No rotation imposed on the initial model.
!     m_grid = no : Eulerian if lagr = no.
!
!  A    : differential rotation parameter [km]
!  beta : ratio of the rotational energy to the gravitational binding
!   energy
!-----------------------------------------------------------------------
!      rot = ye: Rotation imposed on the initial model
!      rot = no: No rotation imposed
!-----------------------------------------------------------------------

rot                         no                                          rot
rot                                 1.00000000E+03                      A
rot                                 2.50000000E-03                      beta

!-----------------------------------------------------------------------
!      y_shft      : zone wiggle switch
!      dy_shift    : fraction of a zone width to wiggle grid
!      ncy_shift   : cycle number at which to commence zone wiggling
!      tb_dy_shift : time after bounce at which to stop zone wiggling
!-----------------------------------------------------------------------

y_shft                      no                                          y_shft
y_shft                              1.00000000E-01                      dy_shift
y_shft                     100                                          ncy_shift
y_shft                              5.00000000E-03                      tb_dy_shift

!-----------------------------------------------------------------------
!      v_diff : grid-aligned shock smoothing parameter
!-----------------------------------------------------------------------

v_diff                              0.07500000E+00                      v_diff

!-----------------------------------------------------------------------
!      v_trans_0 = ye: transverse velocities 0 above shock.
!      v_trans_0 = no: transverse velocities 0 above shock computed.
!-----------------------------------------------------------------------

v_tran                      ye                                          v_trans_0

!-----------------------------------------------------------------------
!      sub_cy_yz = ye : subcycle yz hydrodynamics relative to x hydrodynamics.
!      sub_cy_yz = no : yz subvycle option off.
!-----------------------------------------------------------------------

sub_cy                      no                                          sub_cy_yz

!-----------------------------------------------------------------------
!      t_step_xyz = ye : hydro time step set to minimum of xyz hydro.
!      t_step_xyz = no : hydro time step set to minimum of x hydro.
!-----------------------------------------------------------------------

g_tstp                      no                                          t_step_xyz

!-----------------------------------------------------------------------
!      i_grav = 1 : spherical symmetric Newtonian gravity
!      i_grav = 2 : Newtonian Poisson solver
!-----------------------------------------------------------------------

i_grav                       1                                          i_grav

!-----------------------------------------------------------------------
!      nu_equil = ye   : neutrino equilibration on.
!      nu_equil = no   : neutrino equilibration off.
!      rho_equilibrate : density above which to equilibrate neutrinos
!      t_equilibrate   : time after bounce for nu_equil = ye
!-----------------------------------------------------------------------

nequil                      no                                          nu_equil
nequil                              1.00000000E+14                      rho_equilibrate 
nequil                              1.00000000E+01                      t_equilibrate

!-----------------------------------------------------------------------
!      G_trns = ye : Perform Galilean transformation at each time step
!      G_trns = no : Omit Galilean transformation
!-----------------------------------------------------------------------

G_trns                      no                                          G_trns

!-----------------------------------------------------------------------
!      rezn = ye: Rezone radial grid upon start or restart.
!      rezn = no: No rezoning of radial grid.
!-----------------------------------------------------------------------

rezn                        no                                          rezn

rezn                         2                                          n_eulgrid
rezn                        70                                          n1zoom
rezn                        90                                          n2zoom
rezn                       115                                          n3zoom
rezn                                2.00000000E-03                      r_1   
rezn                                8.00000000E-01                      r_2   
rezn                                1.25000000E+00                      r_3   
rezn                                0.97143000E+00                      zoome1
rezn                                1.00001000E+00                      zoome2
rezn                                0.89000000E+00                      zoome3

!-----------------------------------------------------------------------
!         Geometry
!
!      ndim    : number of dimensions.
!
!      ngeomx  : geometry of first dimension
!      ngeomy  : geometry of second dimension
!      ngeomz  : geometry of third dimension
!
!      nleftx  : left boundary condition of first dimension
!      nlefty  : left boundary condition of second dimension
!      nleftz  : left boundary condition of third dimension
!
!      nrightx : right boundary condition of first dimension
!      nrighty : right boundary condition of second dimension
!      nrightz : right boundary condition of third dimension
!
!-----------------------------------------------------------------------

geom                         3                                          ndim
geom                         2                                          ngeomx
geom                         4                                          ngeomy
geom                         5                                          ngeomz
geom                         0                                          nleftx
geom                         6                                          nrightx
geom                         0                                          nlefty
geom                         0                                          nrighty
geom                         0                                          nleftz
geom                         0                                          nrightz

!-----------------------------------------------------------------------
!      imin: innerrmost radial zone of evh1 modules; (EVH1 only)
!      imax: outermost radial zone of evh1 modules; (EVH1 only)
!
!-----------------------------------------------------------------------

evh1zn                       1                                          imin
evh1zn                     256                                          imax

!-----------------------------------------------------------------------
!      jmin: innerrmost spatial y-zone (theta in cylindrical spherical
!       coordinates). (EVH1 only)
!      jmax: outermost spatial y-zone (theta in cylindrical spherical
!       coordinates). (EVH1 only)
!-----------------------------------------------------------------------

evh1zn                       1                                          jmin
evh1zn                      16                                          jmax

!-----------------------------------------------------------------------
!      kmin: outermost spatial z-zone (phi in spherical coordinates)
!        (EVH1 only)
!      kmax: outermost spatial z-zone (phi in spherical coordinates)
!        (EVH1 only)
!-----------------------------------------------------------------------

evh1zn                       1                                          kmin
evh1zn                       4                                          kmax

!-----------------------------------------------------------------------
!      vel_ns : Neutron star velocity (cm s^{-2})
!-----------------------------------------------------------------------

vel_ns                              0.00000000E+00                      vel_ns

!-----------------------------------------------------------------------
!     xmin : minimum value of x-coordinate (used only if not given in
!      the initial model data)
!     xmax : maximum value of x-coordinate (used only if not given in
!      the initial model data)
!-----------------------------------------------------------------------

pb_dim                              0.00000000E+00                      xmin
pb_dim                              2.00000000E+09                      xmax

!-----------------------------------------------------------------------
!      ymin: minimum value of y-coordinate (1/pi radians for ngeomy >= 3)
!      ymax: maximum value of y-coordinate (1/pi radians for ngeomy >= 3)
!-----------------------------------------------------------------------

pb_dim                              0.00000000E+00                      ymin
pb_dim                              1.00000000E+00                      ymax

!-----------------------------------------------------------------------
!      zmin: minimum value of y-coordinate (1/pi radians for ngeomy >= 3)
!      zmax: maximum value of y-coordinate (1/pi radians for ngeomy >= 3)
!-----------------------------------------------------------------------

pb_dim                              0.00000000E+00                      zmin
pb_dim                              2.00000000E+00                      zmax

