subroutine checkoverlap( txy, xlo,xhi, ylo, yhi, has_overlap)
implicit none
real(kind=8) :: txy(2,0:2), xlo, xhi, ylo, yhi
logical :: has_overlap
!
! check whether a triangle and a rectangle has any overlap
!

integer, parameter ::  npoints = 2
real(kind=8), parameter :: one = 1.0d0
real(kind=8), parameter :: zero = 0.0d0
real(kind=8) :: lam1, lam2, lam3, xc, yc
integer :: ip, iq, ipp1, iqp1

integer, parameter :: lda = 2
integer, parameter :: ldb = 2
integer, parameter :: nrhs = 1
integer, parameter :: lwork = 100

real(kind=8) :: work(lwork)
real(kind=8) :: A(lda,2), B(ldb,nrhs)

real(kind=8) :: qxy(2,0:3)
real(kind=8) :: x1,y1, x2,y2, x3,y3, x4,y4

integer :: info, m, n
character :: trans
!debug-begin
integer, save :: ncalls = 0

ncalls = ncalls + 1
if (ncalls .ge. 5000) then
  stop 'in checkoverlap '
endif
!debug-end

  has_overlap = .false.

! --------------------------------------
! check whether triangle is in rectangle
! --------------------------------------
  do ip=0,npoints
  do iq=0,npoints
    lam1 = dble(ip)/dble(npoints)
    lam2 = dble(iq)/dble(npoints)
    lam3 = one - lam1 - lam2

    xc = txy(1,0)*lam1 + txy(1,1)*lam2 + txy(1,2)*lam3
    yc = txy(2,0)*lam1 + txy(2,1)*lam2 + txy(2,2)*lam3

    has_overlap = (xlo.le.xc).and.(xc.le.xhi).and.                  &
                  (ylo.le.yc).and.(yc.le.yhi)
    if (has_overlap) then
      return
    endif
  enddo
  enddo


! --------------------------------------
! check whether rectangle is in triangle
! --------------------------------------
  do ip=0,npoints
  do iq=0,npoints
    lam1 = dble(ip)/dble(npoints)
    lam2 = dble(iq)/dble(npoints)
    xc = lam1 * xlo + (one-lam1)*xhi
    yc = lam2 * ylo + (one-lam2)*yhi

    ! ---------------------------------
    ! solve for barycentric coordinates
    ! lam1 * v0 + lam2 * v1 + (1-lam1-lam2)*v2 = vc
    ! [v0-v2 | v1 - v2] [lam1,lam2] = vc - v2
    ! ---------------------------------
    A(1,1) = txy(1,0) - txy(1,2)
    A(1,2) = txy(1,1) - txy(1,2)
    A(2,1) = txy(2,0) - txy(2,2)
    A(2,2) = txy(2,1) - txy(2,2)
    B(1,1) = xc - txy(1,2)
    B(2,1) = yc - txy(2,2)

    trans = 'N'
    m = 2
    n = 2
    call dgels( trans, m,n, nrhs, A, lda, B, ldb, work, lwork, info)
    lam1 = B(1,1)
    lam2 = B(2,1)
    has_overlap = (lam1 .ge. zero).and.(lam2 .ge.zero).and. &
                  (one .ge. lam1 + lam2)
    if ((info.eq.0).and.(has_overlap)) then
       return
    endif
   enddo
   enddo


! ---------------------------
! check intersection of sides
! ---------------------------
  qxy(1,0) = xlo
  qxy(2,0) = ylo

  qxy(1,1) = xhi
  qxy(2,1) = ylo

  qxy(1,2) = xhi
  qxy(2,2) = yhi

  qxy(1,3) = xlo
  qxy(2,3) = yhi

  do ip=0,3
  do iq=0,2
     ipp1 = mod(ip+1,3)
     iqp1 = mod(iq+1,2)

     x1 = txy(1,iq)
     y1 = txy(2,iq)
     x2 = txy(1,iqp1)
     y2 = txy(2,iqp1)

     x3 = qxy(1,ip)
     y3 = qxy(2,ip)
     x4 = qxy(1,ipp1)
     y4 = qxy(2,ipp1)

!    -----------------------------------------------------------------
!    check intersection
!    v1*lam1 + (1-lam1)*v2 = lam2*v3 + (1-lam2)*v4
!    (v1-v2)*lam1 + v2  = (v3-v4)*lam2 + v4
!   [v1-v2 | v4-v3] (lam1,lam2) = v4 - v2
!    -----------------------------------------------------------------
     A(1,1) = x1-x2
     A(2,1) = y1-y2
     A(1,2) = x4-x3
     A(2,2) = y4-y3

     B(1,1) = x4-x2
     B(2,1) = y4-y2
    
     trans = 'N'
     m = 2
     n = 2
     call dgels( trans, m,n,nrhs, A, lda, B, ldb, work, lwork, info)
     lam1 = B(1,1)
     lam2 = B(2,1)
     has_overlap = (zero.le.lam1).and.(lam1.le.one).and.    &
                   (zero.le.lam2).and.(lam2.le.one)
     if (has_overlap) then
       return
     endif

   enddo
   enddo

   has_overlap = .false.
   return
end subroutine checkoverlap
