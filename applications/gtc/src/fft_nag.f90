! wrapper to call NAG FFT routines

subroutine fftr1d(isign,irank,scale,x,y,icount)
  use precision
  implicit none
  integer, intent(in) :: isign,irank,icount
  real(wp), intent(in) :: scale
  real(wp), intent(inout), dimension(0:irank-1) :: x
  complex(wp), intent(inout), dimension(0:irank/2) :: y
  interface
     subroutine scfftnag(isign,n,scale,x,y)
       use precision
       implicit none
       integer, intent(in) :: isign,n
       real(wp), intent(in) :: scale
       real(wp), intent(in), dimension(0:n-1) :: x
       complex(wp), intent(out), dimension(0:n/2) :: y
     end subroutine scfftnag
     subroutine csfftnag(isign,n,scale,x,y)
       use precision
       implicit none
       integer, intent(in) :: isign,n
       real(wp), intent(in) :: scale
       complex(wp), intent(in), dimension(0:n/2) :: x
       real(wp), intent(out), dimension(0:n-1) :: y
     end subroutine csfftnag
  end interface
  if (icount.gt.0.and.icount.le.3) then
     if(isign==1)then
        call scfftnag(1,irank,scale,x,y)
     else
        call csfftnag(-1,irank,scale,y,x)
     endif
  endif
end subroutine fftr1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! only for parallel direction
subroutine fftc1d(isign,irank,scale,x)
  use precision
  implicit none
  integer, intent(in) :: isign,irank
  real(wp), intent(in) :: scale
  complex(wp), intent(inout), dimension(0:irank-1) :: x
  complex(wp),  dimension(0:irank-1) :: y
  interface
     subroutine ccfftnag(isign,n,scale,x,y)
       use precision
       implicit none
       integer, intent(in) :: isign,n
       real(wp), intent(in) :: scale
       complex(wp), intent(in), dimension(0:n-1) :: x
       complex(wp), intent(out), dimension(0:n-1) :: y
     end subroutine ccfftnag
  end interface
  call ccfftnag(isign,irank,scale,x,y)
  x=y
end subroutine fftc1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! scfftnag, csfftnag, ccfftnag are work-alikes to the corresponding
! Cray routines, except that they omit the last three arguments
! (table, work, isys)

subroutine scfftnag(isign,n,scale,x,y)
  use precision
  implicit none

  integer, intent(in) :: isign,n
  real(wp), intent(in) :: scale
  real(wp), intent(in), dimension(0:n-1) :: x
  complex(wp), intent(out), dimension(0:n/2) :: y

  real(wp) :: rscale, iscale
  integer :: ifail,j,n2,nj
  real(doubleprec), dimension(0:n-1) :: dz,dwork

  if (abs(isign).ne.1) return

  dz=real(x,doubleprec)
  ifail=0
  call c06faf(dz,n,dwork,ifail)
  rscale=sqrt(real(n,wp))*scale
  iscale=isign*sqrt(real(n,wp))*scale
  y(0)=rscale*real(dz(0),wp)
  n2=(n-1)/2
  do j=1,n2
     nj=n-j
     y(j)=cmplx(rscale*real(dz(j),wp),iscale*real(dz(nj),wp),wp)
  end do
  if (mod(n,2).eq.0) then
     y(n2+1)=rscale*real(dz(n2+1),wp)
  end if

  return
end subroutine scfftnag

subroutine csfftnag(isign,n,scale,x,y)
  use precision
  implicit none

  integer, intent(in) :: isign,n
  real(wp), intent(in) :: scale
  complex(wp), intent(in), dimension(0:n/2) :: x
  real(wp), intent(out), dimension(0:n-1) :: y

  real(wp) :: rscale, iscale
  integer :: ifail,j,n2,nj
  real(doubleprec), dimension(0:n-1) :: dz,dwork

  if (abs(isign).ne.1) return

  rscale=sqrt(real(n,wp))*scale
  iscale=isign*sqrt(real(n,wp))*scale
  dz(0)=real(rscale*real(x(0),wp),doubleprec)
  n2=(n-1)/2
  do j=1,n2
     nj=n-j
     dz(j)=real(rscale*real(x(j),wp),doubleprec)
     dz(nj)=real(iscale*aimag(x(j)),doubleprec)
  end do
  if (mod(n,2).eq.0) then
     dz(n2+1)=real(rscale*real(x(n2+1),wp),doubleprec)
  end if
  call c06fbf(dz,n,dwork,ifail)
  y=real(dz,wp)
  return
end subroutine csfftnag

subroutine ccfftnag(isign,n,scale,x,y)
  use precision
  implicit none

  integer, intent(in) :: isign,n
  real(wp), intent(in) :: scale
  complex(wp), intent(in), dimension(0:n-1) :: x
  complex(wp), intent(out), dimension(0:n-1) :: y

  real(wp) :: rscale, iscale
  integer :: ifail,j
  real(doubleprec), dimension(0:n-1) :: dr,di,dwork

  if (abs(isign).ne.1) return

  rscale=sqrt(real(n,wp))*scale
  iscale=isign*sqrt(real(n,wp))*scale
  do j=0,n-1
     dr(j)=real(real(x(j),wp),doubleprec)
     di(j)=real(isign*aimag(x(j)),doubleprec)
  end do
  call c06fcf(dr,di,n,dwork,ifail)
  do j=0,n-1
     y(j)=cmplx(rscale*real(dr(j),wp),iscale*real(di(j),wp),wp)
  end do
  return
end subroutine ccfftnag
