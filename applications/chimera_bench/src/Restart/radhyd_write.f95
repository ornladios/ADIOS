SUBROUTINE radhyd_write( ndump )
!-----------------------------------------------------------------------
!
!    File:         radhyd_write
!    Module:       radhyd_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/00
!
!    Purpose:
!      To dump the radhyd keys.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!   ndump       : unit number to dump the data
!
!    Output arguments:
!        none
!
!    Include files:
!      cycle_module, nucbrn_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : nse_c, nprint, ndim, ngeomx, ngeomy, ngeomz, &
& nleftx, nlefty, nleftz, nrightx, nrighty, nrightz, imin, imax, jmin, jmax, &
& kmin, kmax, ncycle, ncymax, time, t_start, t_bounce, t_stop, tb_stop, &
& dtnph, dtnmh, lagr, rezn, m_grid, t_bounce_lagr_chg, t_bounce_mgrd_chg, &
& regrid, grid_frac, int_pre_b, int_post_b, rho_regrid, rot, A, beta, y_shft, &
& dy_shift, ncy_shift, tb_dy_shift, v_trans_0, xmin_i, xmax_i, ymin_i, ymax_i, &
& ymax_i, zmin_i, zmax_i, sub_cy_yz, t_step_xyz, nu_equil, rho_equilibrate, &
& t_equilibrate, i_grav, v_diff, G_trns, vel_ns
USE rezone_module, ONLY : n_lgrgrid, n_eulgrid, n1zoom, n2zoom, n3zoom, m_1, &
& m_2, m_3, r_1, r_2, r_3, zoome1, zoome2, zoome3

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ndump           ! unit number to write restart file

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!  Header
!-----------------------------------------------------------------------

    1 FORMAT ('!                  \\\\\                ///// ')
    2 FORMAT ('!                         RADHYD KEYS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!-----------------------------------------------------------------------
!  Cycle number
!-----------------------------------------------------------------------

   21 FORMAT ('!      ncycle : the cycle number of the calculation.')
   22 FORMAT ('!      ncymax : probelm termination crierion. Calculation terminated when ncycle > ncymax.')
   23 FORMAT ('cycle ',14x,i10,42x,'ncycle')
   24 FORMAT ('cycle ',14x,i10,42x,'ncymax')

!-----------------------------------------------------------------------
!  Times
!-----------------------------------------------------------------------

   31 FORMAT ('!      time : the elapsed time since the initiation of the calculation.')
   32 FORMAT ('!      t_start : the time at the initiation of the calculation.')
   33 FORMAT ('!      t_bounce : the time of core bounce.')
   34 FORMAT ('!      t_stop : the elapsed time after which calculation is terminated.')
   35 FORMAT ('!      tb_stop : the elapsed time from bounce after which calculation is terminated.')
   36 FORMAT ('sptime',29x,1pe15.8,22x,'time')
   37 FORMAT ('sptime',29x,1pe15.8,22x,'t_start')
   38 FORMAT ('sptime',29x,1pe15.8,22x,'t_bounce')
   39 FORMAT ('sptime',29x,1pe15.8,22x,'t_stop')
   40 FORMAT ('sptime',29x,1pe15.8,22x,'tb_stop')

!-----------------------------------------------------------------------
!  Time steps
!-----------------------------------------------------------------------

   41 FORMAT ('!      dtnph : the current time step.')
   42 FORMAT ('!      dtmph : the preceding time step.')
   43 FORMAT ('sptime',29x,1pe15.8,22x,'dtnph')
   44 FORMAT ('sptime',29x,1pe15.8,22x,'dtnmh')

!-----------------------------------------------------------------------
!  Lagrangian-Eulerian switch
!-----------------------------------------------------------------------

   51 FORMAT ('!      lagr              : Lagrangian-Eulerian switch.')
   52 FORMAT ('!      t_bounce_lagr_chg : time after bounce for lagr = no')
   53 FORMAT ('lagr  ',22x,a2,42x,'lagr')
   54 FORMAT ('lagr  ',29x,1pe15.8,22x,'tb_lagr')

!-----------------------------------------------------------------------
!  Moving radial grid switch
!-----------------------------------------------------------------------

   61 FORMAT ('!      m_grid            : Moving radial grid switch.')
   62 FORMAT ('!      t_bounce_mgrd_chg : time after bounce for m_grid = no.')
   63 FORMAT ('m_grid',22x,a2,42x,'m_grid')
   64 FORMAT ('m_grid',29x,1pe15.8,22x,'tb_mgrd')

!-----------------------------------------------------------------------
!  Regridding switch
!-----------------------------------------------------------------------

   71 FORMAT ('!      regrid     : Regridding switch.')
   72 FORMAT ('!      grid_frac  : time after bounce to resume regridding [ms].')
   73 FORMAT ('!      int_pre_b  : number of cycles between successive regrids before bounce.')
   74 FORMAT ('!      int_post_b : number of cycles between successive regrids after bounce.')
   75 FORMAT ('!      rho_regrid : minimum density to regrid [g cm^{-3}].')
   76 FORMAT ('regrid',22x,a2,42x,'regrid')
   77 FORMAT ('regrid',14x,i10,42x,'int_pre_b')
   78 FORMAT ('regrid',14x,i10,42x,'int_post_b')
   79 FORMAT ('regrid',29x,1pe15.8,22x,'grid_frac')
   80 FORMAT ('regrid',29x,1pe15.8,22x,'rho_regrid')

!-----------------------------------------------------------------------
!  Regridding switch
!-----------------------------------------------------------------------

   81 FORMAT ('!      rot  : Impose initial rotation switch.')
   82 FORMAT ('!      A    : differential rotation parameter.')
   83 FORMAT ('!      beta : rotation magnitude parameter.')
   84 FORMAT ('rot   ',22x,a2,42x,'rot   ')
   85 FORMAT ('rot   ',29x,1pe15.8,22x,'A     ')
   86 FORMAT ('rot   ',29x,1pe15.8,22x,'beta  ')

!-----------------------------------------------------------------------
!  Zone wiggling parameters
!-----------------------------------------------------------------------

   91 FORMAT ('!      y_shft = ye : zone wiggle on.')
   92 FORMAT ('!      y_shft = no : zone wiggle off.')
   93 FORMAT ('!      dy_shift    : fraction of a zone width to wiggle grid.')
   94 FORMAT ('!      ncy_shift   : cycle number at which to commence zone wiggling.')
   95 FORMAT ('!      tb_dy_shift : time after bounce at which to stop zone wiggling.')
   96 FORMAT ('y_shft',22x,a2,42x,'y_shft')
   97 FORMAT ('y_shft',29x,1pe15.8,22x,'dy_shift')
   98 FORMAT ('y_shft',14x,i10,42x,'ncy_shift')
   99 FORMAT ('y_shft',29x,1pe15.8,22x,'tb_dy_shift')

!-----------------------------------------------------------------------
!  Grid-aligned ahock smoothing parameter
!-----------------------------------------------------------------------

  101 FORMAT ('!      v_diff: grod-aligned shock smoothing parameter.')
  102 FORMAT ('v_diff',29x,1pe15.8,22x,'v_diff')

!-----------------------------------------------------------------------
!  Transverse velocity zero above shock switch
!-----------------------------------------------------------------------

  106 FORMAT ('!      v_trans_0 = ye : transverse velocities 0 above shock.')
  107 FORMAT ('!      v_trans_0 = no : transverse velocities 0 above shock computed.')
  108 FORMAT ('v_tran',22x,a2,42x,'v_trans_0')

!-----------------------------------------------------------------------
!  yz-subcycling switch
!-----------------------------------------------------------------------

  111 FORMAT ('!      sub_cy_yz = ye : subcycle yz hydrodynamics relative to x hydrodynamics.')
  112 FORMAT ('!      sub_cy_yz = no : yz subvycle option off')
  113 FORMAT ('sub_cy',22x,a2,42x,'sub_cy_yz')

!-----------------------------------------------------------------------
!  Global hydro time-step switch switch
!-----------------------------------------------------------------------

  116 FORMAT ('!      t_step_xyz = ye : hydro time step set to minimum of xyz hydro.')
  117 FORMAT ('!      t_step_xyz = no : hydro time step set to minimum of x hydro.')
  118 FORMAT ('g_tstp',22x,a2,42x,'t_step_xyz')

!-----------------------------------------------------------------------
!  Gravitation switch
!-----------------------------------------------------------------------

  121 FORMAT ('!      i_grav = 1 : spherical symmetric gravity.')
  122 FORMAT ('!      i_grav = 2 : nonspherical components computed from by a Newtonian Poisson solver.')
  123 FORMAT ('i_grav',14x,i10,42x,'i_grav')

!-----------------------------------------------------------------------
!  Neutrino equilibration switch
!-----------------------------------------------------------------------

  131 FORMAT ('!      nu_equil = ye   : neutrino equilibration on.')
  132 FORMAT ('!      nu_equil = no   : neutrino equilibration off.')
  133 FORMAT ('!      rho_equilibrate : density above which to equilibrate neutrinos with matter.')
  134 FORMAT ('!      t_equilibrate   : time after bounce for nu_equil = ye.')
  135 FORMAT ('nequil',22x,a2,42x,'nu_equil')
  136 FORMAT ('nequil',29x,1pe15.8,22x,'rho_equilibrate')
  137 FORMAT ('nequil',29x,1pe15.8,22x,'t_equilibrate')

!-----------------------------------------------------------------------
!  Rezoning switch
!-----------------------------------------------------------------------

  141 FORMAT ('!      G_trns : Galilean transformation.')
  142 FORMAT ('G_trns',22x,a2,42x,'G_trns')

!-----------------------------------------------------------------------
!  Rezoning switch
!-----------------------------------------------------------------------

  148 FORMAT ('!      rezn : Regridding switch.')
  149 FORMAT ('rezn  ',22x,a2,42x,'rezn')

!-----------------------------------------------------------------------
!  Rezoning parameters
!-----------------------------------------------------------------------

  151 FORMAT ('!      n_lgrgrid : Lagrangian regridding toggle.')
  152 FORMAT ('!      n_eulgrid : Eulerian regridding toggle.')
  153 FORMAT ('!      m_1 : mass of first zone in a Lagrangian regrid.')
  154 FORMAT ('!      m_2 : mass separating grid 1 from grid 2 in a Lagrangian regrid.')
  155 FORMAT ('!      m_3 : mass separating grid 2 from grid 3 in a Lagrangian regrid.')
  156 FORMAT ('!      r_1 : radius of first zone in an Eularian regrid.')
  157 FORMAT ('!      r_2 : radius separating grid 1 from grid 2 in an Eularian regrid.')
  158 FORMAT ('!      r_3 : radius separating grid 2 from grid 3 in an Eularian regrid.')
  159 FORMAT ('!      n1zoom : number of zones between m_1 and m_2, Lagrangian (r_1 and r_2, Eularian).')
  160 FORMAT ('!      n2zoom : number of zones between m_2 and m_3, Lagrangian (r_2 and r_3, Eularian).')
  161 FORMAT ('!      n3zoom : number of zones between m_3 and edge, Lagrangian (r_2 and edge, Eularian).')
  162 FORMAT ('!      zoome1 : mass ratio, Lagrangian (radius ratio, Eulerian) for first n1 zones.')
  163 FORMAT ('!      n2zoom : mass ratio, Lagrangian (radius ratio, Eulerian) for next n2zoom zones.')
  164 FORMAT ('!      n3zoom : mass ratio, Lagrangian (radius ratio, Eulerian) for remaining n3zoom zones.')
  165 FORMAT ('rezn  ',14x,i10,42x,'n_lgrgrid')
  166 FORMAT ('rezn  ',14x,i10,42x,'n_eulgrid')
  167 FORMAT ('rezn  ',29x,1pe15.8,22x,'m_1   ')
  168 FORMAT ('rezn  ',29x,1pe15.8,22x,'m_2   ')
  169 FORMAT ('rezn  ',29x,1pe15.8,22x,'m_3   ')
  170 FORMAT ('rezn  ',29x,1pe15.8,22x,'r_1   ')
  171 FORMAT ('rezn  ',29x,1pe15.8,22x,'r_2   ')
  172 FORMAT ('rezn  ',29x,1pe15.8,22x,'r_3   ')
  173 FORMAT ('rezn  ',14x,i10,42x,'n1zoom')
  174 FORMAT ('rezn  ',14x,i10,42x,'n2zoom')
  175 FORMAT ('rezn  ',14x,i10,42x,'n3zoom')
  176 FORMAT ('rezn  ',29x,1pe15.8,22x,'zoome1')
  177 FORMAT ('rezn  ',29x,1pe15.8,22x,'zoome2')
  178 FORMAT ('rezn  ',29x,1pe15.8,22x,'zoome3')

!-----------------------------------------------------------------------
!  Geometry
!-----------------------------------------------------------------------

  181 FORMAT ('!      ndim : number of geometric dimensions.')
  182 FORMAT ('!      ngeomx : x-geometry flag.')
  183 FORMAT ('!      ngeomy : y-geometry flag.')
  184 FORMAT ('!      ngeomz : z-geometry flag.')
  185 FORMAT ('!      nleftx : lower x-boundary condition flag.')
  186 FORMAT ('!      nlefty : lower y-boundary condition flag.')
  187 FORMAT ('!      nleftz : lower z-boundary condition flag.')
  188 FORMAT ('!      nrightx : upper x-boundary condition flag.')
  189 FORMAT ('!      nrighty : upper y-boundary condition flag.')
  190 FORMAT ('!      nrightz : upper z-boundary condition flag.')
  191 FORMAT ('geom  ',14x,i10,42x,'ndim')
  192 FORMAT ('geom  ',14x,i10,42x,'ngeomx')
  193 FORMAT ('geom  ',14x,i10,42x,'ngeomy')
  194 FORMAT ('geom  ',14x,i10,42x,'ngeomz')
  195 FORMAT ('geom  ',14x,i10,42x,'nleftx')
  196 FORMAT ('geom  ',14x,i10,42x,'nrightx')
  197 FORMAT ('geom  ',14x,i10,42x,'nlefty')
  198 FORMAT ('geom  ',14x,i10,42x,'nrighty')
  199 FORMAT ('geom  ',14x,i10,42x,'nleftz')
  200 FORMAT ('geom  ',14x,i10,42x,'nrightz')

!-----------------------------------------------------------------------
!  Zoning
!-----------------------------------------------------------------------

  201 FORMAT ('!      imin : inner physical x (radial) index.')
  202 FORMAT ('!      imax : outer physical x (radial) index.')
  203 FORMAT ('!      jmin : inner physical y (angular) index.')
  204 FORMAT ('!      jmax : outer physical y (angular) index.')
  205 FORMAT ('!      kmin : inner physical z index.')
  206 FORMAT ('!      kmax : outer physical z index.')
  207 FORMAT ('evh1zn',14x,i10,42x,'imin')
  208 FORMAT ('evh1zn',14x,i10,42x,'imax')
  209 FORMAT ('evh1zn',14x,i10,42x,'jmin')
  210 FORMAT ('evh1zn',14x,i10,42x,'jmax')
  211 FORMAT ('evh1zn',14x,i10,42x,'kmin')
  212 FORMAT ('evh1zn',14x,i10,42x,'kmax')

!-----------------------------------------------------------------------
!  Neutron star velocity
!-----------------------------------------------------------------------

  216 FORMAT ('!      vel_ns : Neutron star velocity (cm s^{-2}).')
  217 FORMAT ('vel_ns',29x,1pe15.8,22x,'vel_ns')

!-----------------------------------------------------------------------
!  Problem size
!-----------------------------------------------------------------------

  231 FORMAT ('!      xmin : minimum value of x-coordinate.')
  232 FORMAT ('!      xmax : maximum value of x-coordinate.')
  233 FORMAT ('!      ymin : minimum value of y-coordinate.')
  234 FORMAT ('!      ymax : maximum value of y-coordinate.')
  235 FORMAT ('!      zmin : minimum value of z-coordinate.')
  236 FORMAT ('!      zmax : maximum value of z-coordinate.')
  237 FORMAT ('pb_dim',29x,1pe15.8,22x,'xmin')
  238 FORMAT ('pb_dim',29x,1pe15.8,22x,'xmax')
  239 FORMAT ('pb_dim',29x,1pe15.8,22x,'ymin')
  240 FORMAT ('pb_dim',29x,1pe15.8,22x,'ymax')
  241 FORMAT ('pb_dim',29x,1pe15.8,22x,'zmin')
  242 FORMAT ('pb_dim',29x,1pe15.8,22x,'zmax')

!-----------------------------------------------------------------------
!  Document the dump
!-----------------------------------------------------------------------

 1001 FORMAT (' ***Radhyd keys dump written at cycle     ',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!
!                      \\\\\ WRITE KEYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Header
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!-----------------------------------------------------------------------
!  Cycle number
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,22)
WRITE (ndump,15)
WRITE (ndump,23) ncycle
WRITE (ndump,24) ncymax

!-----------------------------------------------------------------------
!  Times
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,31)
WRITE (ndump,32)
WRITE (ndump,33)
WRITE (ndump,34)
WRITE (ndump,35)
WRITE (ndump,15)
WRITE (ndump,36) time
WRITE (ndump,37) t_start
WRITE (ndump,38) t_bounce
WRITE (ndump,39) t_stop
WRITE (ndump,40) tb_stop

!-----------------------------------------------------------------------
!  Time steps
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,41)
WRITE (ndump,42)
WRITE (ndump,15)
WRITE (ndump,43) dtnph
WRITE (ndump,44) dtnmh

!-----------------------------------------------------------------------
!  Lagrangian-Eulerian switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,51)
WRITE (ndump,52)
WRITE (ndump,15)
WRITE (ndump,53) lagr
WRITE (ndump,54) t_bounce_lagr_chg

!-----------------------------------------------------------------------
!  Moving radial grid switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,61)
WRITE (ndump,62)
WRITE (ndump,15)
WRITE (ndump,63) m_grid
WRITE (ndump,64) t_bounce_mgrd_chg

!-----------------------------------------------------------------------
!  Moving radial grid switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,71)
WRITE (ndump,72)
WRITE (ndump,73)
WRITE (ndump,74)
WRITE (ndump,75)
WRITE (ndump,15)
WRITE (ndump,76) regrid
WRITE (ndump,77) int_pre_b
WRITE (ndump,78) int_post_b
WRITE (ndump,79) grid_frac
WRITE (ndump,80) rho_regrid

!-----------------------------------------------------------------------
!  Impose initial rotation switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,81)
WRITE (ndump,82)
WRITE (ndump,83)
WRITE (ndump,15)
WRITE (ndump,84) rot
WRITE (ndump,85) A
WRITE (ndump,86) beta

!-----------------------------------------------------------------------
!  Zone wiggling parameters
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,91)
WRITE (ndump,92)
WRITE (ndump,93)
WRITE (ndump,94)
WRITE (ndump,95)
WRITE (ndump,15)
WRITE (ndump,96) y_shft
WRITE (ndump,97) dy_shift
WRITE (ndump,98) ncy_shift
WRITE (ndump,99) tb_dy_shift

!-----------------------------------------------------------------------
!  Grid-aligned ahock smoothing parameter
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,101)
WRITE (ndump,15)
WRITE (ndump,102) v_diff

!-----------------------------------------------------------------------
!  Transverse velocity zero above shock switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,106)
WRITE (ndump,107)
WRITE (ndump,15)
WRITE (ndump,108) v_trans_0

!-----------------------------------------------------------------------
!  yz-subcycling switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,111)
WRITE (ndump,112)
WRITE (ndump,15)
WRITE (ndump,113) sub_cy_yz

!-----------------------------------------------------------------------
!  Global hydro time-step switch switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,116)
WRITE (ndump,117)
WRITE (ndump,15)
WRITE (ndump,118) t_step_xyz

!-----------------------------------------------------------------------
!  Gravitation switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,121)
WRITE (ndump,122)
WRITE (ndump,15)
WRITE (ndump,123) i_grav

!-----------------------------------------------------------------------
!  Neutrino equilibration switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,131)
WRITE (ndump,132)
WRITE (ndump,133)
WRITE (ndump,134)
WRITE (ndump,15)
WRITE (ndump,135) nu_equil
WRITE (ndump,136) rho_equilibrate
WRITE (ndump,137) t_equilibrate

!-----------------------------------------------------------------------
!  Galilean transformation switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,141)
WRITE (ndump,15)
WRITE (ndump,142) G_trns

!-----------------------------------------------------------------------
!  Rezoning switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,148)
WRITE (ndump,15)
WRITE (ndump,149) rezn

!-----------------------------------------------------------------------
!  Rezoning parameters
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,151)
WRITE (ndump,152)
WRITE (ndump,153)
WRITE (ndump,154)
WRITE (ndump,155)
WRITE (ndump,156)
WRITE (ndump,157)
WRITE (ndump,158)
WRITE (ndump,159)
WRITE (ndump,160)
WRITE (ndump,161)
WRITE (ndump,162)
WRITE (ndump,163)
WRITE (ndump,164)
WRITE (ndump,15)
WRITE (ndump,165) n_lgrgrid
WRITE (ndump,166) n_eulgrid
WRITE (ndump,167) m_1
WRITE (ndump,168) m_2
WRITE (ndump,169) m_3
WRITE (ndump,170) r_1
WRITE (ndump,171) r_2
WRITE (ndump,172) r_3
WRITE (ndump,173) n1zoom
WRITE (ndump,174) n2zoom
WRITE (ndump,175) n3zoom
WRITE (ndump,176) zoome1
WRITE (ndump,177) zoome2
WRITE (ndump,178) zoome3

!-----------------------------------------------------------------------
!  Geometry
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,181)
WRITE (ndump,182)
WRITE (ndump,183)
WRITE (ndump,184)
WRITE (ndump,185)
WRITE (ndump,186)
WRITE (ndump,187)
WRITE (ndump,188)
WRITE (ndump,189)
WRITE (ndump,190)
WRITE (ndump,15)
WRITE (ndump,191) ndim
WRITE (ndump,192) ngeomx
WRITE (ndump,193) ngeomy
WRITE (ndump,194) ngeomz
WRITE (ndump,195) nleftx
WRITE (ndump,196) nrightx
WRITE (ndump,197) nlefty
WRITE (ndump,198) nrighty
WRITE (ndump,199) nleftz
WRITE (ndump,200) nrightz

!-----------------------------------------------------------------------
!  Zoning
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,201)
WRITE (ndump,202)
WRITE (ndump,203)
WRITE (ndump,204)
WRITE (ndump,205)
WRITE (ndump,206)
WRITE (ndump,15)
WRITE (ndump,207) imin
WRITE (ndump,208) imax
WRITE (ndump,209) jmin
WRITE (ndump,210) jmax
WRITE (ndump,211) kmin
WRITE (ndump,212) kmax

!-----------------------------------------------------------------------
!  Neutron star velocity
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,216)
WRITE (ndump,15)
WRITE (ndump,217) vel_ns

!-----------------------------------------------------------------------
!  Problem size
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,231)
WRITE (ndump,232)
WRITE (ndump,233)
WRITE (ndump,234)
WRITE (ndump,235)
WRITE (ndump,236)
WRITE (ndump,15)
WRITE (ndump,237) xmin_i
WRITE (ndump,238) xmax_i
WRITE (ndump,239) ymin_i
WRITE (ndump,240) ymax_i
WRITE (ndump,241) zmin_i
WRITE (ndump,242) zmax_i

!-----------------------------------------------------------------------
!  Record the dump
!-----------------------------------------------------------------------

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE radhyd_write
