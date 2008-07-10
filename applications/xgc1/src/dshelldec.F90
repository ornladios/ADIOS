      subroutine dshelldec(n,dx,iperm)
      integer n
      real*8 dx(n)
      real*8 dtemp,dxert
      integer iperm(n)

      integer j,k,m,maxmn
      integer itemp
      

      do i=1,n
        iperm(i) = i
      enddo



      m = n
23000 if((1.eq.1))then
      m = m/2
      if(m.eq.0)then
      return
      endif
      maxmn = n - m

      do 23004 j = 1,maxmn
      do 23006 k = j,1,-m
      dxert = dx(iperm(k+m))
      if(dxert.lt.dx(iperm(k)))then
      goto 23007
      endif

      itemp = iperm(k+m)
      iperm(k+m) = iperm(k)
      iperm(k) = itemp

!      dtemp = dx(k+m)
!      dx(k+m) = dx(k)
!      dx(k) = dtemp

23006 continue
23007 continue

! end do k
23004 continue
23005 continue

! end do j
      goto 23000
      endif
23001 continue

! end while
      end
