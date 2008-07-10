c----------------------------------------------------------------------c
c                                                                      c
c    File:         cnvccmpe                                            c
c    Module:       cnvccmpe                                            c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         12/7/95                                             c
c                                                                      c
c    Purpose:                                                          c
c      To compute the internal energy e given rho, entropy, and ye.    c
c       Iteration is initially Newton-Rhapson, but defaults to         c 
c       bisection if there is a convergence failure. The quantity t    c
c       is inputed as the initial guess at the temperature and         c
c       outputed as the computed temperature.                          c
c                                                                      c
c    Subprograms called:                                               c
c       eqstt_x, esrgn_x                                                   c
c                                                                      c
c    Input arguments:                                                  c
c  j        : radial zone whose internal energy and change in          c
c              temperature is to be computed.                          c
c  rho      : density of radial zone j.                                c
c  s        : entropy of radial zone j.                                c
c  ye       : electron fraction of radial zone j.                      c
c  t        : temperature of radial zone j (inputed as an initial      c
c              guess.                                                  c
c    Output arguments:                                                 c
c  e        : internal energy of radial zone j.                        c
c  t        : temperature of radial zone j.                            c
c                                                                      c
c    Include files:                                                    c
c        numerical_module                                              c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvccmpe(j,rho,s,ye,t,e)

      USE numerical_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer it,itmax,itrgn,ivare,ivars,j,jj,nprint
      double precision dedd,dedt,dedy,dsdd,dsdt,dsdy,dt,e,
     *                 rho,s,stest,t,tmax,tmin,tol,tp,ye
c                                                                      c
c                                                                      c
      data nprint  /6           /  ! unit number for output
      data itmax   /40          /  ! maximum number of iterations
      data itrgn   /20          /  ! esrgn_x not called if it > itrgn
      data tol     /1.d-6       /  ! convergence tolerence
c                                                                      c
c                                                                      c
 1001 format (' s will not converge in subroutine cnvccmpe - NR')
 1003 format (' j=',i3,' s=',1pe14.7,' stest=',1pe14.7,' tp=',1pe14.7,
     * ' dt=',1pe14.7,' rho=',1pe10.3,' ye=',1pe10.3)
 2001 format (' s will not converge in subroutine cnvccmpe - Bisection')
 2003 format (' j=',i3,' s=',1pe14.7,' stest=',1pe14.7,' tmin=',1pe14.7,
     * ' tmax=',1pe14.7,' rho=',1pe10.3,' ye=',1pe10.3)
c----------------------------------------------------------------------c
c        Initialize                                                    c
c----------------------------------------------------------------------c
      jj           = 1
      ivare        = 2
      ivars        = 3
c----------------------------------------------------------------------c
c        Newton-Rhapson iteration; use jj = 1 in eos calls to avoid    c
c         reindexing the equation of state table for zone j.           c
c----------------------------------------------------------------------c
      tp           = t
      do 1000 it = 1,itmax
       if ( it .le. itrgn ) call esrgn_x(jj, i_ray, rho, tp, ye) 
       call eqstt_x( ivars, jj, i_ray, rho, tp, ye, stest, dsdd, dsdt,
     * dsdy)
       if ( dabs( s - stest ) .le. tol * s ) then
        t          = tp
        call eqstt_x( ivare, jj, i_ray, rho, tp, ye, e, dedd, dedt, dedy )
        return
       end if ! ( dabs( s - stest ) le tol*s
       dt          = ( s - stest )/( dsdt + epsilon )
       tp          = tp + dt
       tp          = dmax1( 0.7d+00*t, dmin1( 1.3d+00*t, tp ) )
 1000 continue
      write (nprint,1001)
      write (nprint,1003) j,s,stest,tp,dt,rho,ye
c----------------------------------------------------------------------c
c        Initialize for bisection                                      c
c----------------------------------------------------------------------c
      tmin         = 1.d+7    ! minimum temperature (guess)
      tmax         = 1.d+12   ! maximum temperature (guess)
c----------------------------------------------------------------------c
c        Bisection iteration; use jj = 1 in eos calls to avoid         c
c         reindexing the equation of state table for zone j.           c
c----------------------------------------------------------------------c
      do 2000 it = 1,itmax
       tp          = half*( tmin + tmax )
       if ( it .le. itrgn ) call esrgn_x( jj, i_ray, rho, tp, ye )
       call eqstt_x( ivars, jj, i_ray, rho, tp, ye, stest, dsdd, dsdt, 
     * dsdy )
       if ( dabs(stest - s) .le. tol * s ) then
        t           = tp
        call eqstt_x( ivare, jj, i_ray, rho, tp, ye, e, dedd, dedt, dedy )
        return
       end if ! dabs(stest - s) le tol*s
       if ( stest .le. s ) then
        tmin       = tp
       else
        tmax       = tp
       end if ! stest le s
 2000 continue
      write (nprint,2001)
      write (nprint,2003) j,s,stest,tmin,tmax,rho,ye
c                                                                      c
c                                                                      c
      return
      end