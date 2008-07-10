c----------------------------------------------------------------------c
c                                                                      c
c    File:         cvmixpsi                                            c
c    Module:       cvmixpsi                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/9/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the change in the neutrino distribution              c
c       in each radial zone due to convective mixing.                  c
c                                                                      c
c    Subprograms called:                                               c
c      cnvcpsi0                                                        c
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
c      convect_module, mdl_cnfg_module, nu_dist_module,                c
c      nu_energy_grid_module, prb_cntl_module                          c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cvmixpsi(jmin,jmax)

      USE array_module
      USE numerical_module
      USE physcnst_module

      USE convect_module
      USE mdl_cnfg_module
      USE nu_dist_module
      USE nu_energy_grid_module, ONLY : nnugp
      USE prb_cntl_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,jm1,jmax,jmaxm1,jmin,jminp1,jp1,k,n
      double precision dpsi0m,psi0aa,psi0injjm1,psi0injjp1,psi0outjjm1,
     *                 psi0outjjp1
c                                                                      c
c                                                                      c
      dimension psi0aa(nez),dpsi0m(nx,nez)
      dimension psi0outjjp1(nx,nez),psi0injjp1(nx,nez),
     *          psi0outjjm1(nx,nez),psi0injjm1(nx,nez)
c----------------------------------------------------------------------c
c        Return if ipsimix = 0.                                        c
c----------------------------------------------------------------------c
      if ( ipsimix .eq. 0 ) return
c----------------------------------------------------------------------c
c        Sum over all neutrino types. Omit convective mixing of        c
c         n-neutrinos if nnugp(n) = 0.                                 c
c----------------------------------------------------------------------c
      do 6000 n = 1,nnu
       if ( nnugp(n) .eq. 0 ) go to 6000
c----------------------------------------------------------------------c
c  Compute the neutrino occupation numbers of a blob after it moves    c
c   outward or inward between radial zones j and j+1 undergoing a      c
c   density change equal to the density difference between the two     c
c   zones.                                                             c
c                                                                      c
c    psi0outjjp1(j,k) : zeroth moment of the neutrino occupation       c
c           number in a blob as it leaves radial zone j on its way     c
c           to radial zone j+1.                                        c
c    psi0injjp1(j,k)  : zeroth moment of the neutrino occupation       c
c           number in a blob as it arrives at radial zone j from       c
c           radial zone j+1.                                           c
c    psi0outjjm1(j,k) : zeroth moment of the neutrino occupation       c
c           number in a blob as it leaves radial zone j on its way     c
c           to radial zone j-1.                                        c
c    psi0injjm1(j,k)  : zeroth moment of the neutrino occupation       c
c           number in a blob as it arrives at radial zone j from       c
c           radial zone j-1.                                           c
c                                                                      c
c  Note that the innermost zone (jmin) exchanges matter only with the  c
c   zone adjacent to it on the outside, and the outermost zone (jmax)  c
c   exchanges matter only with the zone adjacent to it on the inside.  c
c                                                                      c
c  Subroutine cnvcpsi0(n,ji,jf,rhoi,rhof,psi0aa) computes the zeroth   c
c   moments of the neutrino occupation numbers having initial values   c
c   equal to the zeroth moments of radial zone ji when the density     c
c   is changed adiabatically from rhoi to rhof. The final values of    c
c   the zeroth moments are returned in array psi0aa.                   c
c                                                                      c
c        Initialize.                                                   c
c----------------------------------------------------------------------c
       do 1000 k = 1,nnugp(n)
        do 100 j = jmin,jmax
         psi0outjjp1(j,k)
     *             = zero
         psi0injjp1(j,k)
     *             = zero
         psi0outjjm1(j,k)
     *             = zero
         psi0injjm1(j,k)
     *             = zero
  100   continue
 1000  continue
c----------------------------------------------------------------------c
c        Determine psi0outjjp1(jmin,k) and psi0injjp1(jmin,k) for the  c
c         innermost radial zone if the convective velocities between   c
c         it and the outwardly adjacent radial zone are nonzero.       c
c----------------------------------------------------------------------c
       if ( ulcnvct(jmin) .ne. zero ) then
        do 1200 k = 1,nnugp(n)
         psi0outjjp1(jmin,k)
     *             = psi0(jmin,k,n)
 1200   continue
        jminp1     = jmin + 1
        call cnvcpsi0(n,jminp1,jmin,rho(jminp1),rho(jmin),psi0aa)
        do 1400 k = 1,nnugp(n)
         psi0injjp1(jmin,k)
     *             = psi0aa(k)
 1400   continue
       end if ! ulcnvct(jmin) ne zero
c----------------------------------------------------------------------c
c        Determine psi0outjjm1(jmax,k) and psi0injjm1(jmax,k) for the  c
c         outermost radial zone if the convective velocities between   c
c         it and the inwardly adjacent radial zone are nonzero.        c
c----------------------------------------------------------------------c
       if ( ulcnvct(jmax-1) .ne. zero ) then
        do 1600 k = 1,nnugp(n)
         psi0outjjm1(jmax,k)
     *             = psi0(jmax,k,n)
 1600   continue
        jmaxm1        = jmax - 1
        call cnvcpsi0(n,jmaxm1,jmax,rho(jmaxm1),rho(jmax),psi0aa)
        do 1800 k = 1,nnugp(n)
         psi0injjm1(jmax,k)
     *             = psi0aa(k)
 1800   continue
       end if ! ulcnvct(jmax-1) ne zero
c----------------------------------------------------------------------c
c        Determine psi0outjjp1(j,k) and psi0injjp1(j,k) for radial     c
c         zone j if the convective velocities between it and the       c
c         outwardly adjacent radial zone are nonzero.                  c
c----------------------------------------------------------------------c
       do 3000 j = jmin+1,jmax-1
c                                                                      c
        if ( ulcnvct(j) .ne. zero ) then
         do 2200 k = 1,nnugp(n)
          psi0outjjp1(j,k)
     *             = psi0(j,k,n)
 2200    continue
         jp1       = j + 1
         call cnvcpsi0(n,jp1,j,rho(jp1),rho(j),psi0aa)
         do 2400 k = 1,nnugp(n)
          psi0injjp1(j,k)
     *             = psi0aa(k)
 2400    continue
        end if ! ulcnvct(j) ne zero
c----------------------------------------------------------------------c
c        Determine  psi0outjjm1(j,k) and psi0injjm1(j,k) for radial    c
c         zone j if the convective velocities between it and the       c
c         inwardly adjacent radial zone are nonzero.                   c
c----------------------------------------------------------------------c
        if ( ulcnvct(j-1) .ne. zero ) then
         do 2600 k = 1,nnugp(n)
          psi0outjjm1(j,k)
     *             = psi0(j,k,n)
 2600    continue
         jm1       = j - 1
         call cnvcpsi0(n,jm1,j,rho(jm1),rho(j),psi0aa)
         do 2800 k = 1,nnugp(n)
          psi0injjm1(j,k)
     *             = psi0aa(k)
 2800    continue
        end if ! ulcnvct(j-1) ne zero
c                                                                      c
 3000  continue
c----------------------------------------------------------------------c
c  The change in the zero moments of the neutrino occupation numbers   c  
c   in a zone will be given by the NET zero moments brought into the   c
c   zone by convection between the zone of interest and its two        c
c   adjacent zones. For the innermost zone, there is no convective     c
c   source or sink of neutrinos below it. Consequently, the NET        c
c   change in the neutrino zero moment of this zone is given by        c
c   what neutrinos it exchanges with the zone just above it. For the   c
c   outermost zone, its NET change in the neutrino zero moment is      c
c   given by what neutrinos it exchanges with the zone below it.       c
c   Below, we first compute the net change in the mass weighted        c
c   neutrino zero moment (dpsi0m) of a zone j, then its new neutrino   c
c   zero moment (psi0).                                                c
c----------------------------------------------------------------------c
       do 4000 k = 1,nnugp(n)
        dpsi0m(jmin,k)
     *             = ( psi0injjp1(jmin,k) - psi0outjjp1(jmin,k) )
     *             * dmcnvct(jmin  )
        dpsi0m(jmax,k)
     *             = ( psi0injjm1(jmax,k) - psi0outjjm1(jmax,k) )
     *             * dmcnvct(jmax-1)
c                                                                      c
        do 400 j=jmin+1,jmax-1
         dpsi0m(j,k)
     *             = ( psi0injjm1(j,k) - psi0outjjm1(j,k) )
     *             * dmcnvct(j-1)
     *             + ( psi0injjp1(j,k) - psi0outjjp1(j,k) )
     *             * dmcnvct(j  )
  400   continue
 4000  continue
c                                                                      c
       do 5000 k = 1,nnugp(n)
        do 500 j=jmin,jmax
         psi0(j,k,n)
     *             = psi0(j,k,n) + dpsi0m(j,k)/dmrst(j)
  500   continue
 5000  continue
c                                                                      c
 6000 continue
c                                                                      c
c                                                                      c
      return
      end