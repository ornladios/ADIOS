c----------------------------------------------------------------------c
c                                                                      c
c    File:         ludcmp                                              c
c    Module:       ludcmp                                              c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         8/21/00                                             c
c                                                                      c
c    Purpose:                                                          c
c      Given an nxn matrix a, with physical dimension np, this         c
c       replaces it by the LU decomposition of a rowwise permutation   c
c       of itself. A is output arranged as                             c
c                                                                      c
c              b11  b12  b13  b14                                      c
c              a21  b22  b23  b24                                      c
c              a31  a32  b33  b34                                      c
c              a41  a42  a43  b44                                      c
c                                                                      c
c       indx is an output vector which records the row permutation     c
c       effected by the partial pivoting; d is output + or - 1         c
c       depending on whether the number of row interchanges was        c
c       even or odd, respectively. This routine is used in combina-    c
c       tion with lubksb to solve linear equations or invert a         c
c       matrix.                                                        c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c   a         : nxn matrix                                             c
c   n         : dimension of a                                         c
c   np        : physical dimension of a                                c
c                                                                      c
c    Output arguments:                                                 c
c   a         : LU decomposition of a                                  c
c   indx      : record of row permutation                              c
c   d         : odd or even permutation number index                   c
c                                                                      c
c    Variables that must be passed through common:                     c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      none                                                            c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine ludcmp( a, n, np, indx, d )

c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer i,imax,indx,j,j1,j2,k,n,np,nmax
      double precision a,aamax,d,dum,one,sum,tiny,vv,zero
c                                                                      c
c                                                                      c
      parameter (nmax = 100,tiny = 1.d-20)
c                                                                      c
c                                                                      c
      dimension a(np,np),indx(n),vv(nmax)
c                                                                      c
c                                                                      c
      data zero /0.0d+00/
      data one  /1.0d+00/

 1001 FORMAT (81(1pe11.3))

c----------------------------------------------------------------------c
c        Loop over rows to get the implicit scaling information.       c
c----------------------------------------------------------------------c
      d             = one
      do 1000 i = 1,n
       aamax        = zero
       do 100 j = 1,n
        if ( dabs(a(i,j)) .gt. aamax ) aamax = dabs(a(i,j))
  100  continue
       IF ( aamax .eq. zero ) THEN
         WRITE (6,*) 'Singular matrix.'
         DO j1 = 1,n
           WRITE (6,1001) (a(j1,j2),j2 = 1,n)
         END DO
         STOP
       END IF
       vv(i)        = one/aamax
 1000 continue
c----------------------------------------------------------------------c
c        Loop over columns of Crout's method.                          c
c----------------------------------------------------------------------c
      do 2000 j = 1,n
       if ( j .gt. 1 ) then
        do 200 i = 1,j-1
         sum        = a(i,j)
         if ( i .gt. 1 ) then
          do 20 k = 1,i-1
           sum      = sum - a(i,k) * a(k,j)
   20     continue
          a(i,j)    = sum
         end if ! i > 1
  200   continue
       end if ! j > 1
c----------------------------------------------------------------------c
c        Initialize the search for the largest pivot.                  c
c----------------------------------------------------------------------c
       aamax        = zero
       do 210 i = j,n
        sum         = a(i,j)
        if ( j .gt. 1 ) then
         do 21 k = 1,j-1
          sum       = sum - a(i,k) * a(k,j)
   21    continue
         a(i,j)     = sum
        end if ! j > 1
        dum         = vv(i) * dabs(sum)
        if ( dum .ge. aamax ) then
         imax       = i
         aamax      = dum
        end if ! dum ge aamax
  210  continue
c----------------------------------------------------------------------c
c        Do we need to interchange rows?                               c
c----------------------------------------------------------------------c
       if ( j .ne. imax ) then
        do 220 k = 1,n
         dum        = a(imax,k)
         a(imax,k)  = a(j,k)
         a(j,k)     = dum
  220   continue
        d           = -d
        vv(imax)    = vv(j)
       end if ! j ne imax
       indx(j)      = imax
c----------------------------------------------------------------------c
c        If the pivot element is zero the matrix is singular (at       c
c         least to the precision of the algorithm). For some applica-  c
c         tions on singular matrices, it is desirable to substitute    c
c         tiny for zero.                                               c
c----------------------------------------------------------------------c
       if ( a(j,j) .eq. zero ) a(j,j) = tiny
c----------------------------------------------------------------------c
c        Now, finally, divide by the pivot element.                    c
c----------------------------------------------------------------------c
       if ( j .ne. n ) then
        dum         = one/a(j,j)
        do 230 i = j+1,n
         a(i,j)     = a(i,j) * dum
  230   continue
       end if ! j ne n
c----------------------------------------------------------------------c
c        Go back for the next column in the reduction.                 c
c----------------------------------------------------------------------c
 2000 continue
c                                                                      c
c                                                                      c
      return
      end