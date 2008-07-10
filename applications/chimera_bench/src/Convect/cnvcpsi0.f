c----------------------------------------------------------------------c
c                                                                      c
c    File:         cnvcpsi0                                            c
c    Module:       cnvcpsi0                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/9/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To update the n-type neutrino occupation number due to          c
c       neutrino advection in a convective blob which evolves          c
c       from radial zone ji at density rhoi to radial zone jf at       c
c       density rhof.                                                  c
c                                                                      c
c    Subprograms called:                                               c
c      cvnucmve                                                        c
c                                                                      c
c    Input arguments:                                                  c
c  n        : neutrino type.                                           c
c  ji       : radial zone from which blob originates.                  c
c  jf       : radial zone at which blob terminates.                    c
c  rhoi     : density of radial zone ji.                               c
c  rhof     : density of radial zone jf.                               c
c                                                                      c
c    Output arguments:                                                 c
c  psi0aa   : neutrino occupation probability in blob after an         c
c              adiabatic compression or expansion from rhoi to rhof.   c
c                                                                      c
c    Include files:                                                    c
c      array_module, numerical_module                                  c
c      convect_module, nu_dist_module, nu_energy_grid_module,          c
c      prb_cntl_module                                                 c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvcpsi0(n,ji,jf,rhoi,rhof,psi0aa)

      USE array_module
      USE numerical_module

      USE convect_module
      USE nu_dist_module
      USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
      USE prb_cntl_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      logical first
      integer j1,jf,ji,k,k0,n
      double precision cvccfp,dpsk,dpskm,dpskp1,psi0a,
     *                 psi0aa,psi0m,rhof,rhoi,rncofa,vdriftk
c                                                                      c
c                                                                      c
      dimension rncofa(30),psi0aa(nez)
c                                                                      c
c                                                                      c
      data first   /.true.  /
c----------------------------------------------------------------------c
c        Return if ilcnvct = 0.                                        c
c----------------------------------------------------------------------c
          if ( ilcnvct .eq. 0 ) return
c----------------------------------------------------------------------c
c        Initialize.                                                   c
c----------------------------------------------------------------------c
      if ( first ) then
       do 1000 k = 1,nnugpmx
        rncofa(k)  = unui(k) * unui(k) * dunui(k)
 1000  continue
       first       = .false.
      end if ! first = true
c----------------------------------------------------------------------c
c     ***Change in psi0 due to advection***                            c
c----------------------------------------------------------------------c
      dpskm        = zero
      dpskp1       = zero
c----------------------------------------------------------------------c
c  Determine k0, the minimum value of k for which                      c
c                                                                      c        
c            vdriftk < adrftcv*ulcnvct                                 c
c                                                                      c        
c   where vdriftk is the drift velocity of neutrinos of energy zone k, c
c   ulcnvct is the convective velocity, and adrftcv is a arbitrary     c
c   parameter of order unity; k0 must be less than or equal to         c
c   nnugp(n) - 2. Neutrinos in zones for which                         c
c                                                                      c        
c            vdriftk > adrftcv*ulcnvct                                 c
c                                                                      c 
c   are not advected with the matter. In this case psi0aa(k) is set    c
c   equal to psi0(jf,k,n) so the the neutrinos of energy k delivered   c
c   to the zone are equal to those removed, resulting in no net        c
c   advection.                                                         c
c----------------------------------------------------------------------c
      k0           = 1
      do 2000 k = 1,nnugp(n)-2
       j1          = min( ji, jf )
       vdriftk     = psi1(j1,k,n)
     *             / ( half * ( psi0(j1,k,n) + psi0(j1+1,k,n) )
     *             + epsilon )
       if ( vdriftk .le. adrftcv * ulcnvct(j1) ) then
        k0         = k
        go to 2500
       end if ! vdriftk le adrftcv * ulcnvct(j1)
       psi0aa(k)   = psi0(jf,k,n)
 2000 continue
      return
c----------------------------------------------------------------------c
c  Calculate change in n-type neutrino occupation number due to        c
c   neutrino advection by the convective blobs.                        c
c                                                                      c
c  cvccfp    - the change in psi0(ji,k,n) due to neutrino advection    c
c               (computed in subroutine 'cvnucmve').                   c
c  dpsk      - amount by which psi0a exceeds 1.                        c
c  dpskp1:   - amount by which psi0(ji,k+1,n) must increase to avoid   c
c               overfilling psi0(ji,k,n).                              c
c  dpskm:    - amount by which psi0(ji,k,n) had to be increased to     c
c               avoid overfilling psi0(ji,k-1,n).                      c
c  psi0a:    - updated value of psi0(ji,k,n).                          c
c----------------------------------------------------------------------c
 2500 continue
      do 3000 k = k0,nnugp(n)
       call cvnucmve(n,ji,k,k0,rhoi,rhof,cvccfp,dpskm,dpskp1)
       psi0a       = psi0(ji,k,n) + cvccfp + dpskp1
c----------------------------------------------------------------------c
c  If psi0a < 0.0, set psia0 = 0.0 and reduce the occupation number    c
c   of the smaller adjacent energy zone so that the total n-neutrino   c
c   number is conserved.                                               c
c----------------------------------------------------------------------c
       if ( psi0a .lt. zero  .and.  k .ne. 1 ) then
        psi0m      = rncofa(k) * psi0a/rncofa(k-1)
        psi0aa(k-1)= psi0aa(k-1) + psi0m
        psi0a      = zero
       end if ! psi0a < 0  .and.  k ne 1
       dpskm       = dpskp1
c----------------------------------------------------------------------c
c  If psi0a > 1.0 and k < nnugp(n), set psi0a = 1.0 and shift          c
c   remaining n-neutrinos to the adjacent energy zone of the higher    c
c   energy so that the total n-neutrino number is conserved.           c
c----------------------------------------------------------------------c
       if ( psi0a .gt. one  .and.  k .lt. nnugp(n) ) then
        dpsk       = psi0a - one
        psi0a      = one
        dpskp1     = rncofa(k) * dpsk/rncofa(k+1)
       else
        dpskp1     = zero
       end if ! psi0a > 1  and  k < nnugp(n)
c----------------------------------------------------------------------c
c        Store psi0a in psi0aa(k).                                     c
c----------------------------------------------------------------------c
       psi0aa(k)   = psi0a
c                                                                      c
 3000 continue
c                                                                      c
      return
      end