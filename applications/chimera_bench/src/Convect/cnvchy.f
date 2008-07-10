c----------------------------------------------------------------------c
c                                                                      c
c    File:         cnvchy                                              c
c    Module:       cnvchy                                              c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         2/16/96                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the turbulent pressure, pcnvct, due to               c
c       convection, and the rate, scnvct, at which thermal energy      c
c       is transferred to the matter due to viscous dissipation        c
c       of turbulent eddies.                                           c
c                                                                      c
c    Subprograms called:                                               c
c      dtgvnpsye                                                       c
c                                                                      c
c    Input arguments:                                                  c
c  jmin       : inner radial zone of region for which convective       c
c                stability is to be tested.                            c
c  jmax       : outer radial zone of region for which convective       c
c                stability is to be tested.                            c
c  i_ray      : index denoting a specific radial ray                   c
c                                                                      c
c    Output arguments:                                                 c
c        none                                                          c
c                                                                      c
c    Include files:                                                    c
c      array_module, numerical_module, physcnst_module                 c
c      convect_module, eos_snc_x_module, mdl_cnfg_module,                c
c      prb_cntl_module, shock_module                                   c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvchy(jmin,jmax,i_ray)

      USE array_module
      USE numerical_module
      USE physcnst_module

      USE convect_module
      USE eos_snc_x_module
      USE mdl_cnfg_module
      USE prb_cntl_module
      USE shock_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,i_ray,jmax,jmin,jp1
      double precision rhojmhin,rhojmhout,rhojphin,
     *                 rhojphout,scnvctj,tguess,wrkjmh,wrkjph
c                                                                      c
c                                                                      c
      dimension scnvctj(nx)
c----------------------------------------------------------------------c
c        Initialize.                                                   c
c----------------------------------------------------------------------c
      do 1000 j = 1,nx
       scnvctj(j)  = zero
       scnvct(j)   = zero
       pcnvct(j)   = zero
 1000 continue
c----------------------------------------------------------------------c
c         if ipcnvct ne 0 compute the convective turbulent pressure.   c
c                                                                      c
c        pcnvct(j): turbulent pressure of zone j (dynes/cm3).          c
c----------------------------------------------------------------------c
      if ( ipcnvct .eq. 0 ) go to 2500
      do 2000 j = jmin,jmax
       pcnvct(j)  = rho(j)*half*( ulcnvct(j)**2 + ulcnvct(j-1)**2 )
 2000 continue
c----------------------------------------------------------------------c
c        If iscnvct ne 0 compute scnvctj(j). the energy released (per  c
c         unit mass of mass zones j and j+1) by the convective         c
c         interchange ot the mass dmcnvct(j) between mass zones j and  c
c         j+1. This energy will be put back into zones j and j+1 in    c
c         form of heat.                                                c
c----------------------------------------------------------------------c
 2500 continue
      if ( iscnvct .eq. 0 ) go to 4500
      do 3000 j = jmin,jmax-1
       if ( dmcnvct(j) .ne. zero ) then
c                                                                      c
        jp1        = j + 1
c----------------------------------------------------------------------c
c        Compute the work done on interchanging the mass dmcnvct(j)    c
c         at constant entropy, electron fraction, and in pressure      c
c         equilibrium with the background between mass zones j and     c
c         j+1.                                                         c
c----------------------------------------------------------------------c
        rhojmhout  = rho(j)
        rhojmhin   = rho(j)
        tguess     = t(j)
        call dtgvnpsye(j,aesv(j,1,i_ray),aesv(jp1,3,i_ray),ye(jp1),
     *                 rhojmhin,tguess)
        rhojphout  = rho(jp1)
        rhojphin   = rho(jp1)
        tguess     = t(jp1)
        call dtgvnpsye(j,aesv(jp1,1,i_ray),aesv(j,3,i_ray),ye(j),
     *                 rhojphin,tguess)
        scnvctj(j) = half * ( rhojphout/rhojphin - rhojmhout/rhojmhin )
     *             * ( g * rstmss(j)/r(j)**2 ) * dr(j)
     *             * dmcnvct(j)/( dmrst(j) + dmrst(jp1) )
c----------------------------------------------------------------------c
c        Compute the extra work done on the mass dmcnvct(j), when      c
c         introduced into mass zone j, to adjust its density to that   c
c         of mass zone j. Do the same for the mass introduced into     c
c         mass zone j+1. These energies must then be removed from      c
c         these respective mass zones.                                 c
c----------------------------------------------------------------------c
        wrkjmh     = ( aesv(j  ,1,i_ray)/rho(j  )**2 )*( rho(j  )
     *             - rhojmhin ) * dmcnvct(j)/dmrst(j)
        wrkjph     = ( aesv(jp1,1,i_ray)/rho(jp1)**2 )*( rho(jp1)
     *             - rhojphin ) * dmcnvct(j)/dmrst(jp1)
        scnvct(j  )= scnvct(j  ) - wrkjmh
        scnvct(jp1)= scnvct(jp1) - wrkjph
c                                                                      c
       end if ! dmcnvct(j) ne 0
 3000 continue
c----------------------------------------------------------------------c
c        Add contributions of scnvctj(j) to scnvct(j).                 c
c----------------------------------------------------------------------c
      do 4000 j = jmin,jmax-1
       scnvct(j)   = scnvct(j) + scnvctj(j-1) + scnvctj(j)
 4000 continue
c                                                                      c
c                                                                      c
 4500 continue
      return
      end