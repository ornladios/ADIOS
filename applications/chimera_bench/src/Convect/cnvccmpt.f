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
c      To compute the temperature given rho, internal energy e,        c
c       and ye. Iteration is initially Newton-Rhapson, but defaults    c 
c       to bisection if there is a convergence failure.                c
c       The quantity t is inputed as the initial guess at the          c
c       temperature and outputed as the computed temperature.          c
c                                                                      c
c    Subprograms called:                                               c
c      eqstt_x, esrgn_x                                                    c
c                                                                      c
c    Input arguments:                                                  c
c  j        : radial zone whose internal energy and change in          c
c              temperature is to be computed.                          c
c  rho      : density of radial zone j.                                c
c  e        : internal energy of radial zone j.                        c
c  ye       : electron fraction of radial zone j.                      c
c  t        : temperature of radial zone j (inputed as an initial      c
c              guess.                                                  c
c    Output arguments:                                                 c
c  t        : temperature of radial zone j.                            c
c                                                                      c
c    Include files:                                                    c
c        numerical_module                                              c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine cnvccmpt(j,rho,e,ye,t)

      USE numerical_module
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer it,itmax,itrgn,ivar,j,jj,nprint
      double precision dedd,dedt,dedy,dt,e,etest,
     *                 rho,t,tbmax,tbmin,tmax,tmin,tol,tp,ye
c                                                                      c
c                                                                      c
      data nprint  /6           /  ! unit number for output
      data itmax   /40          /  ! maximum number of iterations
      data itrgn   /0           /  ! esrgn_x not called if it > itrgn
      data tol     /1.d-6       /  ! convergence tolerence
      data tbmin   /1.d+7       /  ! minimum acceptible temperature (K)
      data tbmax   /1.d+12      /  ! maximum acceptible temperature (K)
c                                                                      c
c                                                                      c
 1001 format (' e will not converge in subroutine cnvccmpt - NR')
 1003 format (' j=',i3,' e=',1pe14.7,' etest=',1pe14.7,' tp=',1pe14.7,
     * ' dt=',1pe14.7,' rho=',1pe10.3,' ye=',1pe10.3)
 2001 format (' e will not converge in subroutine cnvccmpt - Bisection')
 2003 format (' j=',i3,' e=',1pe14.7,' etest=',1pe14.7,' tmin=',1pe14.7,
     * ' tmax=',1pe14.7,' rho=',1pe10.3,' ye=',1pe10.3)
c----------------------------------------------------------------------c
c        Initialize                                                    c
c----------------------------------------------------------------------c
      jj           = j
      ivar         = 2
c----------------------------------------------------------------------c
c        Newton-Rhapson iteration; use jj = 1 in eos calls to avoid    c
c         reindexing the equation of state table for zone j.           c
c----------------------------------------------------------------------c
      tp           = t
      do 1000 it = 1,itmax
       if ( tp .lt. tbmin  .or.  tp .gt. tbmax ) go to 1500
       if ( it .le. itrgn ) call esrgn_x( jj, i_ray, rho, tp, ye )
       call eqstt_x( ivar, jj, i_ray, rho, tp, ye, etest, dedd, dedt, 
     * dedy )
       if ( dabs( e - etest ) .le. tol * e ) then
        t          = tp
        return
       end if ! dabs( e - etest ) le tol*e
       dt          = ( e - etest )/( dedt + epsilon )
       tp          = tp + dt
 1000 continue
 1500 continue
      write (nprint,1001)
      write (nprint,1003) j,e,etest,tp,dt,rho,ye
c----------------------------------------------------------------------c
c        Initialize for bisection                                      c
c----------------------------------------------------------------------c
      tmin         = tbmin    ! minimum temperature (guess)
      tmax         = tbmax    ! maximum temperature (guess)
c----------------------------------------------------------------------c
c        Bisection iteration; use jj = 1 in eos calls to avoid         c
c         reindexing the equation of state table for zone j.           c
c----------------------------------------------------------------------c
      do 2000 it = 1,itmax
       tp          = half * ( tmin + tmax )
       if ( it .le. itrgn ) call esrgn_x( jj, i_ray, rho, tp, ye )
       call eqstt_x( ivar, jj, i_ray, rho, tp, ye, etest, dedd, dedt,
     * dedy )
       if ( dabs(etest - e) .le. tol * e ) then
        t          = tp
        return
       end if ! dabs(etest - e) le tol*e
       if ( etest .le. e ) then
        tmin       = tp
       else
        tmax       = tp
       end if ! etest le e
 2000 continue
      write (nprint,2001)
      write (nprint,2003) j,e,etest,tmin,tmax,rho,ye
      stop
c                                                                      c
c                                                                      c
      return
      end