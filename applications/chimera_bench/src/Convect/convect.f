c----------------------------------------------------------------------c
c                                                                      c
c    File:         convect                                             c
c    Module:       convect                                             c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/1/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To direct the testing for convective nstability and the         c
c       resulting computation of the convective mixing.                c
c                                                                      c
c    Subprograms called:                                               c
c      cnvclstb, cnvcmixl, cnvclvel, cnvcmixm, cvmixye, cvmixe,        c
c      cvmixpsi, cnvchy                                                c
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
c      array_module, numerical_module                                  c
c      convect_module, incrmnt_module, prb_cntl_module                 c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine convect(jmin,jmax,i_ray)

      USE array_module
      USE numerical_module

      USE convect_module
      USE incrmnt_module
      USE prb_cntl_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,jmax,jmin
c----------------------------------------------------------------------c
c        Return if ilcnvct = 0.                                        c
c----------------------------------------------------------------------c
      if ( ilcnvct .eq. 0 ) then
       do 1000 j = 1,nx
        aledoux(j)  =  zero
        psclht(j)   =  zero
        ulmixl(j)   =  zero
        ulcnvct(j)  =  zero
        dmcnvct(j)  =  zero
        pcnvct(j)   =  zero
        scnvct(j)   =  zero
        dyecnvt(j,i_ray) =  zero
 1000  continue
       return
      end if
c----------------------------------------------------------------------c
c        Test instability for Ledoux convection.                       c
c----------------------------------------------------------------------c
      call cnvclstb(jmin,jmax)
c----------------------------------------------------------------------c
c        Compute the mixing lengths.                                   c
c----------------------------------------------------------------------c
      call cnvcmixl(jmin,jmax)
c----------------------------------------------------------------------c
c        Compute mixing length velocities and convective velocities.   c
c----------------------------------------------------------------------c
      call cnvclvel(jmin,jmax)
c----------------------------------------------------------------------c
c        Compute the mass exchanged between radial zones by            c
c         convection.                                                  c
c----------------------------------------------------------------------c
      call cnvcmixm(jmin,jmax)
c----------------------------------------------------------------------c
c        Compute the change in the electron fraction due to            c
c         convective mixing.                                           c
c----------------------------------------------------------------------c
      call cvmixye(jmin,jmax,i_ray)
c----------------------------------------------------------------------c
c        Compute the change in the internal energy and, thence, the    c
c         temperature due to convective mixing.                        c
c----------------------------------------------------------------------c
      call cvmixe(jmin,jmax,i_ray)
c----------------------------------------------------------------------c
c        Compute the change in the zero moments of the neutrino        c
c         occupation numbers due to convective mixing.                 c
c----------------------------------------------------------------------c
      call cvmixpsi(jmin,jmax)
c----------------------------------------------------------------------c
c        Compute the turbulent pressure and the rate at which energy   c
c         is transferred to the matter by viscous dissipation.         c
c----------------------------------------------------------------------c
      call cnvchy(jmin,jmax)
c                                                                      c
c                                                                      c
      return
      end