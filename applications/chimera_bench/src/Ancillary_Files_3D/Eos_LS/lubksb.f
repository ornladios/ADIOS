c----------------------------------------------------------------------c
c                                                                      c
c    File:         lubksb                                              c
c    Module:       lubksb                                              c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         8/22/00                                             c
c                                                                      c
c    Purpose:                                                          c
c      Solves the set of of n linear equations ax = b. Here a is       c
c       input, not as the matrix a but rather as its LU decomposition, c
c       determined by subroutine ludcmp. indx is input as the          c
c       permutation vector returned by ludcmp. b is input as the       c
c       right-hand side vector, and returns with the solution vector   c
c       x. a, n, np, and indx are not modified by this rotine and      c
c       can be left in place for successive calls with different       c
c       right-hand sides b. This routine takes into account the        c
c       possibility that b will begin with many zero elements, so it   c
c       is efficient for use in matrix inversion.                      c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c   a         : nxn LU decomposed matrix of coefficients               c
c   n         : dimension of a                                         c
c   np        : physical dimension of a                                c
c   indx      : record of row permutation                              c
c   b         : constant vector                                        c
c                                                                      c
c    Output arguments:                                                 c
c   b         : solution vector                                        c
c                                                                      c
c    Variables that must be passed through common:                     c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      none                                                            c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine lubksb(a,n,np,indx,b)
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer i,ii,indx,j,ll,n,np
      double precision a,b,sum,zero
c                                                                      c
c                                                                      c
      dimension a(np,np),indx(n),b(n)
c                                                                      c
c                                                                      c
      data zero /0.0d+00/
c----------------------------------------------------------------------c
c        When ii is set to a positive value, it will become the index  c
c         of the first nonvanishing element of b. We now do the        c
c         forward substitution. The only new wrinkle is to unscramble  c
c         the permutation as we go.                                    c
c----------------------------------------------------------------------c
      ii                = 0
      do 1000 i = 1,n
       ll               = indx(i)
       sum              = b(ll)
       b(ll)            = b(i)
       if ( ii .ne. 0 ) then
        do 100 j = ii,i-1
         sum            = sum - a(i,j) * b(j)
  100   continue
       else if ( sum .ne. zero ) then
        ii              = i
       end if ! ii ne 0
       b(i)             = sum
 1000 continue
c----------------------------------------------------------------------c
c        Now we do the back substitution.                              c
c----------------------------------------------------------------------c
      do 2000 i = n,1,-1
       sum              = b(i)
       if ( i .lt. n ) then
       do 200 j = i+1,n
        sum             = sum - a(i,j) * b(j)
  200  continue
       end if
       b(i)             = sum/a(i,i)
 2000 continue
c                                                                      c
c                                                                      c
      return
      end