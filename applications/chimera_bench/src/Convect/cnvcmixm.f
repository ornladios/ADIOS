c----------------------------------------------------------------------c
c                                                                      c
c    File:         cnvcmixm                                            c
c    Module:       cnvcmixm                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/4/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the mass, dmcnvct(j), that is exhanged between       c
c       radial zones j and j+1 when convective velocity between        c
c       these zones is nonzero.                                        c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c  jmin       : inner radial zone of region for which convective       c
c                stability is to be tested.                            c
c  jmax       : outer radial zone of region for which convective       c
c                stability is to be tested.                            c
c                                                                      c
c    Output arguments:                                                 c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      numerical_module, physcnst_module                               c
c      convect_module, mdl_cnfg_module, t_cntrl_module                 c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvcmixm(jmin,jmax)

      USE numerical_module
      USE physcnst_module

      USE convect_module
      USE mdl_cnfg_module
      USE t_cntrl_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,jmax,jmin
      double precision rhoj
c----------------------------------------------------------------------c
c        If ulcnvct(j) > 0, compute the mass, dmcnvct(j),  that is     c
c         exhanged between radial zones j and j+1. A factor of 2*pi    c
c         rather than 4*pi is used to compute the mass flux through    c
c         a zone boundary to take into account the presence of both    c
c         downstreams and upstreams.                                   c
c----------------------------------------------------------------------c
      do 1000 j = jmin,jmax-1
       rhoj        = half * ( rho(j) + rho(j+1) )
       dmcnvct(j)  = 2.d+00 * pi * r(j) * r(j) * rhoj
     *             * ulcnvct(j) * dtnph
 1000 continue
c                                                                      c
c                                                                      c
      return
      end