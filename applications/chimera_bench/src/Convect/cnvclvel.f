c----------------------------------------------------------------------c
c                                                                      c
c    File:         cnvclvel                                            c
c    Module:       cnvclvel                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/3/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the mixing length and convective velocities,         c
c       ulmixl(j) and ulcnvct(j), for those radial zone boundaries     c
c       for which aledoux(j) > 0.                                      c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c      jmin : inner radial zone of region for which convective         c
c              stability is to be tested.                              c
c      jmax : outer radial zone of region for which convective         c
c              stability is to be tested.                              c
c      i_ray: index denoting a specific radial ray                   c
c                                                                      c
c    Output arguments:                                                 c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      numerical_module, physcnst_module                               c
c      convect_module, eos_snc_x_module, mdl_cnfg_module,                c
c      shock_module                                                    c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvclvel(jmin,jmax,i_ray)

      USE numerical_module
      USE physcnst_module

      USE convect_module
      USE eos_snc_x_module
      USE mdl_cnfg_module
      USE shock_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,jmax,jmin,i_ray,jshockmin
      double precision gam1j,pj,pqcrit,rhoj,usound2
c                                                                      c
c                                                                      c
      data pqcrit  /0.3d+00/
c----------------------------------------------------------------------c
c        Initialize.                                                   c
c----------------------------------------------------------------------c
      do 1000 j = 1,nz
       ulmixl(j)   = zero
       ulcnvct(j)  = zero
 1000 continue
c----------------------------------------------------------------------c
c        Compute the sound velocities usound(j) at the radial zone     c
c         boundaries.                                                  c
c----------------------------------------------------------------------c
      do 2000 j = jmin,jmax-1
       pj          = half * ( aesv(j,1,i_ray) + aesv(j+1,1,i_ray) )
       rhoj        = half * ( rho(j) + rho(j+1) )
c      gam1j       = half * ( aesv(j,12,i_ray) + aesv(j+1,12,i_ray) )
       gam1j       = half * ( gam1(j,i_ray) + gam1(j+1,i_ray) )
       usound2     = gam1j * pj/rhoj
       usound(j)   = dsqrt( dabs(usound2) + epsilon )
 2000 continue
c----------------------------------------------------------------------c
c        Compute the mixing length velocities ulmixl(j) at the         c
c         radial zone boundaries.                                      c
c----------------------------------------------------------------------c
      do 3000 j = jmin,jmax-1
       ulmixl(j)   = zero
       if ( aledoux(j) .gt. zero ) then
        ulmixl(j)  = half * lmix(j)
     *              * dsqrt( g * rstmss(j) * aledoux(j)/r(j)**2 )
        ulmixl(j)  = dmin1( ulmixl(j), usound(j) )
       end if ! aledoux(j) > 0
 3000 continue
c----------------------------------------------------------------------c
c        Determine the shock location in order to avoid computing      c
c         convection in the negative entropy profile of the shock.     c
c                                                                      c
c jshock   :  pq_x(jshock)/p(shock) > pqcrit                           c
c jshockmin:  min(jshock)                                              c
c----------------------------------------------------------------------c
      do 4000 j = jmin,jmax-1
       jshockmin   = j 
       if ( pq_x(j,i_ray)/aesv(j,1,i_ray) .gt. pqcrit ) go to 4500
 4000 continue
 4500 continue
c----------------------------------------------------------------------c
c        Compute the convective velocities ulcnvct(j) at the           c
c         radial zone boundaries below the shock.                      c
c----------------------------------------------------------------------c
      do 5000 j = jmin+10,jmax-1
c      uldt        = lmix(j)/dtnph
c      ulcnvctr    = ulcnvct(j)
c      ulcnvct(j)  = dsqrt( uldt**2 + ulmixl(j)**2
c    *             + 2.d+00 * ulcnvct(j) * uldt ) - uldt
       if ( j .le. jshockmin - 3  .and.  rho(j) .ge. 1.d+10 ) then
        ulcnvct(j)   = ulmixl(j)
       else
        ulcnvct(j)   = zero
       end if ! j le jshockmn - 3  and  rho(j) ge 1.d+10
 5000 continue
c                                                                      c
c                                                                      c
      return
      end