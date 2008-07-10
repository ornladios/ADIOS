c----------------------------------------------------------------------c
c                                                                      c
c    File:         cvmixye                                             c
c    Module:       cvmixye                                             c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/4/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the change in the electron fraction in each          c
c       radial zone due to convective mixing.                          c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c  jmin       : inner radial zone of region for which convective       c
c                stability is to be tested.                            c
c  jmax       : outer radial zone of region for which convective       c
c                stability is to be tested.                            c
c  i_ray      : index denoting a specific radial ray                   c
c                                                                      c
c    Output arguments:                                                 c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      numerical_module                                                c
c      convect_module, eos_snc_x_module, incrmnt_module,                 c
c      mdl_cnfg_module, prb_cntl_module                                c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cvmixye(jmin,jmax,i_ray)

      USE numerical_module

      USE convect_module
      USE eos_snc_x_module
      USE incrmnt_module
      USE mdl_cnfg_module
      USE prb_cntl_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,jmax,jmin
c----------------------------------------------------------------------c
c        Initialize                                                    c
c----------------------------------------------------------------------c
      do 100 j = jmin,jmax
       dyecnvt(j,i_ray)  = zero
  100 continue 
c----------------------------------------------------------------------c
c        Return if iyemix = 0.                                         c
c----------------------------------------------------------------------c
      if ( iyemix .eq. 0 ) return
c----------------------------------------------------------------------c
c        Compute dyecnvt(jmin,i_ray) for the innermost radial zone.    c
c----------------------------------------------------------------------c
      dyecnvt(jmin,i_ray)= ( ye(jmin+1) - ye(jmin) ) * dmcnvct(jmin)
     *             / dmrst(jmin)
c----------------------------------------------------------------------c
c        Compute dyecnvt(jmax,i_ray) for the outermost radial zone.    c
c----------------------------------------------------------------------c
      dyecnvt(jmax,i_ray)= ( ye(jmax-1) - ye(jmax) ) * dmcnvct(jmax-1)
     *             / dmrst(jmax)
c----------------------------------------------------------------------c
c        Compute dyecnvt(j,i_ray) for the rest of the radial zones.    c
c----------------------------------------------------------------------c
      do 1000 j = jmin+1,jmax-1 
c                                                                      c
       dyecnvt(j,i_ray) = ( ( ye(j+1) - ye(j  ) ) * dmcnvct(j  )
     *            -  ( ye(j  ) - ye(j-1) ) * dmcnvct(j-1) )
     *            / dmrst(j)
 1000 continue
c                                                                      c
c                                                                      c
      return
      end