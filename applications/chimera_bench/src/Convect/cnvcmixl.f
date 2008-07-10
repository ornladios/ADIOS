c----------------------------------------------------------------------c
c                                                                      c
c    File:         cnvcmixl                                            c
c    Module:       cnvcmixl                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/3/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the array psclht(j), the pressure scale height       c
c       evaluated at the boundary between radial zone j and j+1,       c
c       and from that the mixing length array lmix(j).                 c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c  jmin       : inner radial zone of region for which convective       c
c              stability is to be tested.                              c
c  jmax       : outer radial zone of region for which convective       c
c              stability is to be tested.                              c
c  i_ray      : index denoting a specific radial ray                   c
c                                                                      c
c    Output arguments:                                                 c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      array_module, numerical_module, physcnst_module                 c
c      convect_module, eos_snc_x_module, mdl_cnfg_module,                c
c      shock_module                                                    c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvcmixl(jmin,jmax,i_ray)

      USE array_module
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
      integer j,i_ray,jj,jlower,jmax,jmin,jshockmin,jupper
      double precision gravinv,lmixf,lmixb,pj,pqcrit,rhoj
c                                                                      c
c                                                                      c
      dimension lmixf(nx),lmixb(nx)
c                                                                      c
c                                                                      c
      data pqcrit /0.3d+00/
c----------------------------------------------------------------------c
c        Initialize.                                                   c
c----------------------------------------------------------------------c
      do 1000 j = 1,nx
       lmixf(j)    = zero
       lmixb(j)    = zero
       lmix(j)     = zero
 1000 continue
c----------------------------------------------------------------------c
c        Compute psclht(j).                                            c
c----------------------------------------------------------------------c
      do 2000 j = jmin,jmax-1
       rhoj        = half * ( rho(j+1) + rho(j) )
       pj          = half * ( aesv(j+1,1,i_ray) + aesv(j,1,i_ray) )
       gravinv     = r(j)**2/( g * rstmss(j) )
       psclht(j)   = gravinv * pj/rhoj
 2000       continue
c----------------------------------------------------------------------c
c        Determine the shock location in order to avoid computing      c
c         convection in the negative entropy profile of the shock.     c
c                                                                      c
c jshock   :  pq_x(jshock)/p(shock) > pqcrit                           c
c jshockmin:  min(jshock)                                              c
c----------------------------------------------------------------------c
      do 3000 j = jmin,jmax-1
       jshockmin   = j 
       if ( pq_x(j,i_ray)/aesv(j,1,i_ray) .gt. pqcrit ) go to 3500
 3000 continue
 3500 continue
c----------------------------------------------------------------------c
c        Compute lmixf(j), the formal mixing length.                   c
c----------------------------------------------------------------------c
      do 4000 j = jmin,jmax-1
       lmixf(j)    = alphlmix * psclht(j)
 4000 continue
c----------------------------------------------------------------------c
c        Compute lmixb(j), the width of the convectively unstable      c
c         region.                                                      c
c----------------------------------------------------------------------c
      do 5000 j = jmin,jshockmin-1
       if ( aledoux(j) .gt. zero  .or.  aledoux(j-1) .gt. zero ) then
c----------------------------------------------------------------------c
c            Find outer boundary of convectively unstable region.      c
c----------------------------------------------------------------------c
        do 500 jj = j,jshockmin-1
         jupper    = jj
         if ( aledoux(jj-1) .le. zero ) go to 510
  500   continue
  510   continue
c----------------------------------------------------------------------c
c            Find inner boundary of convectively unstable region.      c
c----------------------------------------------------------------------c
        do 520 jj = j,jmin,-1
         jlower    = jj
         if ( aledoux(jj) .le. zero ) go to 530
  520   continue
  530   continue
c----------------------------------------------------------------------c
c            Compute lmixb(j).                                         c
c----------------------------------------------------------------------c
        lmixb(j)   = r(jupper) - r(jlower-1)
       end if ! aledoux(j) > 0  or  aledoux(j-1) > zero
 5000 continue
c----------------------------------------------------------------------c
c        Compute lmix(j), the minimum of lmixf(j) and lmixb(j)         c
c----------------------------------------------------------------------c
      do 6000 j = jmin,jmax-1
       lmix(j)     = dmin1( lmixf(j), lmixb(j) )
 6000 continue
c                                                                      c
c                                                                      c
      return
      end