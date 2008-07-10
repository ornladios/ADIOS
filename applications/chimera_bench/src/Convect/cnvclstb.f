c----------------------------------------------------------------------c
c                                                                      c
c    File:         cnvclstb                                            c
c    Module:       cnvclstb                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/1/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To test instability for Ledoux convection by computing          c
c       the quantity aledoux(j) for each radial zone j, where          c
c                                                                      c
c                     1     drho            drho                       c
c       aledoux(j) = --- [( ---- )      - ( ---- )     ]               c
c                    rho     dr   star       dr   blob                 c
c                                                                      c
c       aledoux(j) > 0: radial zones j and j+1 are unstable to         c
c                        mixing bu Ledoux convection.                  c
c       aledoux(j) < 0: radial zones j and j+1 are stable to           c
c                        mixing bu Ledoux convection.                  c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c  jmin       : inner radial zone of region for which convective       c
c                stability is to be tested.                            c
c  jmax       : outer radial zone of region for which convective       c
c                stability is to be tested.                            c
!  i_ray      : index denoting a specific radial ray
c                                                                      c
c    Output arguments:                                                 c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      numerical_module                                                c
c      convect_module, eos_snc_x_module, mdl_cnfg_module,                c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvclstb(jmin,jmax,i_ray)

      USE numerical_module

      USE convect_module
      USE eos_snc_x_module
      USE mdl_cnfg_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,i_ray,jmax,jmin
      double precision dpj,drhoj,drj,gamj,pj,rhoj,rhopgj
c----------------------------------------------------------------------c
c        Compute aledoux(j).                                           c
c----------------------------------------------------------------------c
      aledoux(jmin-1)
     *             = zero
      aledoux(jmax)= zero
      do 1000 j = jmin,jmax-1
       rhoj        = half * ( rho(j) + rho(j+1) )
c      rhopgj      = half * ( rho(j  )/( aesv(j  ,1,i_ray) * aesv(j  ,12,i_ray) )
c    *             +      rho(j+1)/( aesv(j+1,1,i_ray)*aesv(j+1,12,i_ray) ) )
       pj          = half * ( aesv(j,1,i_ray) + aesv(j+1,1),i_ray )
c      gamj        = half * ( aesv(j,12,i_ray) + aesv(j+1,12,i_ray) )
       gamj        = half * ( gam1(j,i_ray) + gam1(j+1,i_ray) )
       rhopgj      = rhoj/( pj * gamj )
       drj         = half * ( r(j+1) - r(j-1) )
       drhoj       = rho(j+1) - rho(j)
       dpj         = aesv(j+1,1,i_ray) - aesv(j,1,i_ray)
       aledoux(j)  = ( drhoj - rhopgj * dpj )/( drj * rhoj )
 1000 continue
c                                                                      c
c                                                                      c
      return
      end