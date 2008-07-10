c----------------------------------------------------------------------c
c                                                                      c
c    File:         mprove                                              c
c    Module:       mprove                                              c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         8/22/00                                             c
c                                                                      c
c    Purpose:                                                          c
c      Improves a solution vector x of the linear set of equations     c
c       ax = b. The matrix a, and the vectors b and x are input, as    c
c       is alud, the LU decomposition of a as returned by ludcmp,      c
c       and the vector indx also returned by that routine. On output,  c
c       only x is modified, to an improved set of values.              c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c   a         : matrix of coefficients                                 c
c   alud      : nxn LU decomposed matrix of coefficients               c
c   n         : dimension of a                                         c
c   np        : physical dimension of a                                c
c   indx      : record of row permutation                              c
c   b         : constant vector                                        c
c   x         : solution vector to be improved                         c
c                                                                      c
c    Output arguments:                                                 c
c   x         : improved solution vector                               c
c                                                                      c
c    Variables that must be passed through common:                     c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      none                                                            c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine mprove(a,alud,n,np,indx,b,x)
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer i,indx,j,n,nmax,np
      double precision a,alud,b,r,sdp,x,zero
c                                                                      c
c                                                                      c
      parameter (nmax = 100)
c                                                                      c
c                                                                      c
      dimension a(np,np),alud(np,np),indx(n),b(n),x(n),r(nmax)
c                                                                      c
c                                                                      c
      data zero /0.0d+00/
c----------------------------------------------------------------------c
c        Calculate the right-hand side, accumulating the residual.     c
c----------------------------------------------------------------------c
      do 1000 i = 1,n
       sdp              = -b(i)
       do 100 j = 1,n
        sdp             = sdp + a(i,j) * x(j)
  100  continue
       r(i)             = sdp
 1000 continue
c----------------------------------------------------------------------c
c        Solve for the error term, and subtract it from the old        c
c         solution.                                                    c
c----------------------------------------------------------------------c
      call lubksb(alud,n,np,indx,r)
      do 2000 i = 1,n
       x(i)              = x(i) - r(i)
 2000 continue
c                                                                      c
c                                                                      c
      return
      end