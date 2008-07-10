c----------------------------------------------------------------------c
c                                                                      c
c    File:         cvnucmve                                            c
c    Module:       cvnucmve                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         2/17/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the change in the occupation numbers of the          c
c       neutrinos comoving with convecting blobs.                      c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c  j          : radial zone for which neutrino occupation adjustment   c
c                is to be computed.                                    c
c  k          : neutrino energy zone for which neutrino occupation     c
c                adjustment is to be computed.                         c
c  n          : neutrino type for which neutrino occupation            c
c                adjustment is to be computed.                         c
c  rhoi       : initial density of convective blob.                    c
c  rhof       : final density of convective blob.                      c
c  cvclim(k)  : fraction of neutrinos advected with the matter.        c
c  dpsk       : amount by which psi0(j,k+1,n) must increase to avoid   c
c                overfilling psi0(j,k,n).                              c
c  dpskm      : amount by which psi0(j,k,n) had to be increased to     c
c                avoid overfilling psi0(j,k-1,n).                      c
c                                                                      c
c    Output arguments:                                                 c
c  cvccfp     : the change in psi0(j,k,n) due to neutrino advection.   c
c                                                                      c
c    Include files:                                                    c
c      numerical_module                                                c
c      numerical_module, nu_energy_grid_module                         c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cvnucmve(n,j,k,k0,rhoi,rhof,cvccfp,dpskm1,dpsk)

      USE numerical_module

      USE nu_dist_module
      USE nu_energy_grid_module, ONLY : nnugp
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer j,k,k0,km1,kp1,n
      double precision cvccf2,cvccf3,cvccfk,cvccfm,cvccfp,cvccft,
     *                 dpsk,dpskm1,drat,drat43,erat,eratm1,
     *                 rhof,rhoi
c----------------------------------------------------------------------c
c        drat - ratio of final density to initial density of           c
c                convective blob.                                      c
c----------------------------------------------------------------------c
      drat         = rhof/rhoi
      drat43       = drat**(4.d+00/3.d+00)
c                                                                      c
      if ( drat .ge. one ) then 
c----------------------------------------------------------------------c
c                 *drat > 1.*                                          c
c                  k = k0                                              c
c----------------------------------------------------------------------c
       if ( k .eq. k0 ) then
        erat       = unu(j,k0)/unu(j,k0+1)
        cvccfp     = ( drat - one - erat * ( drat43 - one ) )
     *             * psi0(j,k,n) / ( one - erat )
        return
       end if ! k = k0
c----------------------------------------------------------------------c
c                  k0 < k < nnugp(n) - 1                               c
c----------------------------------------------------------------------c
       if ( k .gt. k0  .and.  k .lt. nnugp(n) - 1 ) then
        km1        = k - 1
        kp1        = k + 1
        eratm1     = unu(j,km1)/unu(j,k)
        erat       = unu(j,k)/unu(j,kp1)
        cvccfk     = ( ( drat - one - erat * ( drat43 - one ) )
     *             * psi0(j,k,n)
     *             / ( one - erat ) )
        if ( dpsk .ne. zero ) then
         cvccfk    = cvccfk + dpsk * ( unu(j,k) - unu(j,km1) )
     *                  / ( unu(j,kp1) - unu(j,k) )
        end if ! dpsk ne 0.
        cvccfm     = ( ( drat43 - drat ) * eratm1 * eratm1 * eratm1
     *             * ( dunu(j,km1)/dunu(j,k) ) * psi0(j,km1,n)
     *             / ( one - eratm1 ) )
        if ( dpskm1 .ne. zero ) then
         cvccfm    = cvccfm - dpskm1 * eratm1 * eratm1
     *             * dunu(j,km1) * ( unu(j,km1) - unu(j,k-2) )
     *             / ( dunu(j,k) * ( unu(j,k) - unu(j,km1) ) )
        end if ! dpskm1 ne 0.
        cvccfp     = cvccfk + cvccfm
        return
       end if ! k > k0 and k < nnugp(n)-1
c----------------------------------------------------------------------c
c                  k = nnugp(n) - 1                                    c
c----------------------------------------------------------------------c
       if ( k .eq. nnugp(n) - 1 ) then
        km1        = k - 1
        kp1        = k + 1
        eratm1     = unu(j,km1)/unu(j,k)
        erat       = unu(j,k)/unu(j,kp1)
        cvccfk     = ( ( drat - one - erat * ( drat43 - one ) )
     *             * psi0(j,k,n)/( one - erat ) )
     *             - ( ( drat43 - drat ) * ( dunu(j,kp1)
     *             / ( dunu(j,k) * erat * erat ) )
     *             * psi0(j,kp1,n)/( one - erat ) )
        if ( dpsk .ne. zero ) then
         cvccfk    = cvccfk + dpsk * ( unu(j,k) - unu(j,km1) )
     *             / ( unu(j,kp1) - unu(j,k) )
        end if ! dpsk ne 0.
        cvccfm     = ( ( drat43 - drat ) * eratm1 * eratm1 * eratm1
     *             * ( dunu(j,km1)/dunu(j,k) ) * psi0(j,km1,n)
     *             / ( one - eratm1 ) )
        if ( dpskm1 .ne. zero ) then
         cvccfm    = cvccfm - dpskm1 * eratm1 * eratm1
     *             * dunu(j,km1) * ( unu(j,km1) - unu(j,k-2) )
     *             / ( dunu(j,k) * ( unu(j,k) - unu(j,km1) ) )
        end if ! dpskm1 ne 0.
        cvccfp     = cvccfk + cvccfm
        return
       end if ! k = nnugp(n) - 1
c----------------------------------------------------------------------c
c                  k = nnugp(n)                                        c
c----------------------------------------------------------------------c
       if ( k .eq. nnugp(n) ) then
        km1        = k - 1
        eratm1     = unu(j,km1)/unu(j,k)
        cvccft     = ( ( drat43 - drat ) * eratm1 * eratm1 * eratm1
     *             * ( dunu(j,km1)/dunu(j,k) ) * psi0(j,km1,n)
     *             / ( one - eratm1 ) )
     *             + ( ( drat43 - drat )
     *             / ( one - eratm1 ) + ( drat - one ) ) * psi0(j,k,n)
        if ( dpskm1 .ne. zero ) then
         cvccft    = cvccft - dpskm1 * eratm1 * eratm1
     *             * dunu(j,km1) * ( unu(j,km1) - unu(j,k-2) )
     *             / ( dunu(j,k) * ( unu(j,k) - unu(j,km1) ) )
        end if ! dpskm1 ne 0.
        cvccfp     = cvccft
        return
       end if ! k = nnugp(n)
c                                                                      c
      else ! drat < one
c----------------------------------------------------------------------c
c                 *drat < 1.*                                          c
c                  k = k0                                              c
c----------------------------------------------------------------------c
       if ( k .eq. k0 ) then
        erat       = unu(j,k0+1)/unu(j,k0)
        cvccfp     = ( ( drat43 - drat ) * erat * erat * erat
     *             * ( dunu(j,k0+1)/dunu(j,k0) ) * psi0(j,k0+1,n)
     *             / ( one - erat ) )
     *             + ( ( drat43 - drat )
     *             / ( one - erat ) + ( drat - one ) ) * psi0(j,k0,n)
        return
       end if ! k = k0
c----------------------------------------------------------------------c
c                  k = k0 + 1                                          c
c----------------------------------------------------------------------c
       if ( k .eq. k0 + 1 ) then
        eratm1     = unu(j,k0+1)/unu(j,k0)
        erat       = unu(j,k0+2)/unu(j,k0+1)
        cvccf2     = ( ( drat - one - eratm1 * ( drat43 - one ) )
     *             * psi0(j,k0+1,n)
     *             / ( one - eratm1 ) )
     *             - ( ( drat43 - drat ) * ( dunu(j,k0)
     *             / ( dunu(j,k0+1) * eratm1 * eratm1 ) )
     *             * psi0(j,k0,n)/( one - eratm1 ) )
        cvccf3     = ( ( drat43 - drat ) * erat * erat * erat
     *             * ( dunu(j,k0+2)/dunu(j,k0+1) ) * psi0(j,k0+2,n)
     *             / ( one - erat ) )
        if ( dpsk .ne. zero ) then
         cvccf3    = cvccf3 + dpsk * ( unu(j,k0+1) - unu(j,k0) )
     *             / ( unu(j,k0+2) - unu(j,k0+1) )
        end if ! dpsk ne 0.
        cvccfp     = cvccf2 + cvccf3
        return
       end if ! k = k0 + 1
c----------------------------------------------------------------------c
c                  k0 + 1 < k < nnugp(n)                               c
c----------------------------------------------------------------------c
       if ( k .gt. k0 + 1  .and.  k .lt. nnugp(n) ) then
        km1        = k - 1
        kp1        = k + 1
        eratm1     = unu(j,k)/unu(j,km1)
        erat       = unu(j,kp1)/unu(j,k)
        cvccfk     = ( ( drat - one - eratm1 * ( drat43 - one ) )
     *             * psi0(j,k,n)
     *             / ( one - eratm1 ) )
        if ( dpskm1 .ne. zero ) then
         cvccfk    = cvccfk - dpskm1
     *             * dunu(j,km1) * ( unu(j,km1) - unu(j,k-2) )
     *             / ( eratm1 * eratm1 * dunu(j,k) * ( unu(j,k)
     *             - unu(j,km1) ) )
        end if ! dpskm1 ne 0.
        cvccft     = ( ( drat43 - drat ) * erat * erat * erat
     *             * ( dunu(j,kp1)/dunu(j,k) ) * psi0(j,kp1,n)
     *             / ( one - erat ) )
        if ( dpsk .ne. zero ) then
         cvccft    = cvccft + dpsk * ( unu(j,k) - unu(j,km1) )
     *             / ( unu(j,kp1) - unu(j,k) )
        end if ! dpsk ne 0.
        cvccfp     = cvccfk + cvccft
        return
       end if ! k > 2  and  k < nnugp(n)
c----------------------------------------------------------------------c
c                  k = nnugp(n)                                        c
c----------------------------------------------------------------------c
       if ( k .eq. nnugp(n) ) then
        km1        = k - 1
        eratm1     = unu(j,k)/unu(j,km1)
        cvccft     = ( ( drat - one - eratm1 * ( drat43 - one ) )
     *             * psi0(j,k,n)
     *             / ( one - eratm1 ) )
        if ( dpskm1 .ne. zero ) then
         cvccft    = cvccft - dpskm1
     *             * dunu(j,km1) * ( unu(j,km1) - unu(j,k-2) )
     *             / ( eratm1 * eratm1 * dunu(j,k) * ( unu(j,k)
     *             - unu(j,km1) ) )
        end if ! dpskm1 ne 0.
        cvccfp     = cvccft
       end if ! k = nnugp(n)
c                                                                      c
      end if ! drat >/< 0
c                                                                      c
c                                                                      c
      return
      end