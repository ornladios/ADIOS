c----------------------------------------------------------------------c
c                                                                      c
c    File:         cvmixe                                              c
c    Module:       cvmixe                                              c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/7/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the change in the energy and, in turn,               c
c       temperature in each radial zone due to convective              c
c       mixing.                                                        c
c                                                                      c
c    Subprograms called:                                               c
c      cnvccmpe, cnvccmpt                                              c
c                                                                      c
c    Input arguments:                                                  c
c  jmin       : inner radial zone of region for which convective       c
c                stability is to be tested.                            c
c  jmax       : outer radial zone of region for which convective       c
c                stability is to be tested.                            c
c  i_ray      : index denoting a specific radial ray                   c
c    Output arguments:                                                 c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      numerical_module                                                c
c      convect_module, eos_snc_x_module, incrmnt_module,                 c
c      mdl_cnfg_module, prb_cntl_module, shock_module                  c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cvmixe(jmin,jmax,i_ray)

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
      integer j,i_ray,jj,jm1,jmax,jmaxm1,jmin,jminp1,jp1
      double precision de,deinjjm1,deinjjp1,deoutjjm1,deoutjjp1,
     *                 du,efnl,tfnl,tp,yep
c                                                                      c
c                                                                      c
      dimension du(nz),efnl(nz),de(nz)
      dimension deoutjjp1(nz),deinjjp1(nz),deoutjjm1(nz),deinjjm1(nz)
c----------------------------------------------------------------------c
c        Return if iemix = 0.                                          c
c----------------------------------------------------------------------c
      if ( iemix .eq. 0 ) return
      jj           = 1
c----------------------------------------------------------------------c
c  Compute the specific internal energy of a blob after it moves       c
c   outward or inward between radial zones j and j+1 undergoing a      c
c   density change equal to the density difference between the two     c
c   zones.                                                             c
c                                                                      c
c    deoutjjp1(j) : specific internal energy of blob as it leaves      c
c                    radial zone j on its way to radial zone j+1.      c
c    deinjjp1(j)  : specific internal energy of blob as it arrives     c
c                    at radial zone j from radial zone j+1.            c
c    deoutjjm1(j) : specific internal energy of blob as it leaves      c
c                    radial zone j on its way to radial zone j-1.      c
c    deinjjm1(j)  : specific internal energy of blob as it arrives     c
c                    at radial zone j from radial zone j-1.            c
c                                                                      c
c  Note that the innermost zone (jmin) exchanges matter only with the  c
c   zone adjacent to it on the outside, and the outermost zone (jmax)  c
c   exchanges matter only with the zone adjacent to it on the inside.  c
c                                                                      c
c  Subroutine cnvccmpe(j,rho,s,ye,t,e) computes the specific internal  c
c   energy e of matter having density rho, entropy s, and electron     c
c   fraction ye, using t as an initial guess and the equation of state c
c   grid indices of radial zone j.                                     c
c----------------------------------------------------------------------c
      do 500 j = jmin,jmax
       deoutjjp1(j)   = zero
       deinjjp1(j)    = zero
       deoutjjm1(j)   = zero
       deinjjm1(j)    = zero
  500 continue
c                                                                      c
      if ( ulcnvct(jmin) .ne. zero ) then
       deoutjjp1(jmin)= aesv(jmin  ,2,i_ray)
       jminp1         = jmin + 1
       tp             = t(jmin)
       call cnvccmpe(jmin,rho(jmin),aesv(jminp1,3,i_ray),ye(jminp1),tp,
     *               deinjjp1(jmin))
      end if ! ulcnvct(jmin) ne 0
c                                                                      c
      if ( ulcnvct(jmax-1) .ne. zero ) then
       deoutjjm1(jmax)= aesv(jmax  ,2,i_ray)
       jmaxm1         = jmax - 1
       tp             = t(jmax)
       call cnvccmpe(jmax,rho(jmax),aesv(jmaxm1,3,i_ray),ye(jmaxm1),tp,
     *               deinjjm1(jmax))
      end if ! ulcnvct(jmax-1) ne 0
c                                                                      c
      do 1000 j = jmin+1,jmax-1
       if ( ulcnvct(j) .ne. zero ) then
        deoutjjp1(j)  = aesv(j  ,2,i_ray)
        jp1           = j + 1
        tp            = t(j)
        call cnvccmpe(j,rho(j),aesv(jp1,3,i_ray),ye(jp1),tp,deinjjp1(j))
       end if ! ulcnvct(j) ne 0
       if ( ulcnvct(j-1) .ne. zero ) then
        deoutjjm1(j)  = aesv(j  ,2,i_ray)
        jm1           = j - 1
        tp            = t(j)
        call cnvccmpe(j,rho(j),aesv(jm1,3,i_ray),ye(jm1),tp,deinjjm1(j))
       end if ! ulcnvct(j-1) ne 0
 1000 continue
c----------------------------------------------------------------------c
c  The change in the specific internal energy in a zone will be given  c
c   by the NET specific internal energy brought into the zone by       c
c   convection between the zone of interest and its two adjacent       c
c   zones. For the innermost zone, there is no convective source or    c
c   sink of specific internal energy below it. Consequently, the NET   c
c   change in the specific internal energy of this zone is given by    c
c   what specific internal energy it exchanges with the zone just      c
c   above it. For the outermost zone, its NET change in specific       c
c   internal energy is given by what it exchanges with the zone below  c
c   it. Below, we first compute the net change in the internal energy  c
c   (du) of a zone j, then its new specific internal energy (efnl).    c
c----------------------------------------------------------------------c
      du(jmin)        = ( deinjjp1(jmin) - deoutjjp1(jmin) )
     *                * dmcnvct(jmin  )
      du(jmax)        = ( deinjjm1(jmax) - deoutjjm1(jmax) )
     *                * dmcnvct(jmax-1)
c                                                                      c
      do 2000 j=jmin+1,jmax-1
       du(j)          = ( deinjjm1(j) - deoutjjm1(j) )*dmcnvct(j-1)
     *                +( deinjjp1(j) - deoutjjp1(j) )*dmcnvct(j  )
 2000 continue
c                                                                      c
      do 3000 j=jmin,jmax
       de(j)          = du(j)/dmrst(j)
       efnl(j)        = aesv(j  ,2,i_ray) + de(j)
 3000 continue 
c----------------------------------------------------------------------c
c  Compute the change in the temperature corresponding to the          c
c   convective change, dyecnvt(j,i_ray), in the electron fraction and  c
c   the convective change, de(j), in the specific internal energy.     c
c----------------------------------------------------------------------c
      do 4000 j=jmin,jmax
       if ( dyecnvt(j,i_ray) .ne. zero  .or.  de(j) .ne. zero ) then
        yep           = ye(j) + dyecnvt(j,i_ray)
        tfnl          = t(j)
        call cnvccmpt(j,rho(j),efnl(j),yep,tfnl)
        dtmpmn(j,9,i_ray) = tfnl - t(j)
       end if ! dyecnvt(j,i_ray) ne 0  or  de(j) ne 0
 4000 continue
c                                                                      c
c                                                                      c
      return
      end