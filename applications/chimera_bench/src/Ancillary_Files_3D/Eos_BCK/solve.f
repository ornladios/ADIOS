c***********************************************************************
c      comments of september 7,1984
c       solves the matrix equation a * x = b where
c     a is an n by n matrix, b and x n-component vectors
c
c
c***********************************************************************
c***********************************************************************
      subroutine solve(n,a,bb,x,fmbad)
c***********************************************************************
      implicit double precision (a-h,o-z)
      logical fmbad
      save
      parameter (nsize = 4)
      dimension a(nsize,nsize),bb(nsize),x(nsize)
      dimension xl(nsize,nsize),u(nsize,nsize),y(nsize)
      do 10 m=1,n
10    xl(m,m)=1.
        do 100 k=1,n
        do 40 j=k,n
        u(k,j)=a(k,j)
        if(k.eq.1)go to 40
        do 310 m=1,k-1
310      u(k,j)=u(k,j)-xl(k,m)*u(m,j)
40    continue
41      if(k.eq.n)go to 100
        do 50 i=k+1,n
        xl(i,k)=a(i,k)
        if(k.eq.1) go to 50
        do 45 m=1,k-1
45      xl(i,k)=xl(i,k)-xl(i,m)*u(m,k)
      if(u(k,k) .eq. 0.) then
      fmbad = .true.
      return
      endif
50      xl(i,k)=xl(i,k)/u(k,k)
100     continue
      do 200 i=1,n
      y(i)=bb(i)
      if(i.eq.1) go to 200
      do 150 k=1,i-1
150   y(i)=y(i)-xl(i,k)*y(k)
200   continue
      do 300 ii=1,n
      i=n+1-ii
      x(i)=y(i)
      if(i.eq.n)go to 300
      do 250 k=i+1,n
250   x(i)=x(i)-u(i,k)*x(k)
300   x(i)=x(i)/u(i,i)
      return
      end
c***********************************************************************
