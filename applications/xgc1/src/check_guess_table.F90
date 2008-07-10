subroutine check_guess_table(grid)
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  real (kind=8) :: x(2),p(3)
  integer :: itr,k(2),i,j,itr_guess,count
  real (kind=8) :: cnum,csum
  real (kind=8) :: psi_interpol, psi
  real (kind=8) :: dcount, cmin,cavg,cmax, cvar, csum2
  integer :: n1,n2

  real (kind=8), parameter :: one = 1.0d0
  real (kind=8), parameter :: zero = 0.0d0

  integer, allocatable,  dimension(:) :: lcount


  integer, parameter :: idebug = 1
  integer, parameter :: iodev = 17

! ------------------------------------------------
! compute average and std variation  of triangles per box cell
! ------------------------------------------------

  n1 = grid%guess_n(1)
  n2 = grid%guess_n(2)
  count = maxval( grid%guess_count(1:n1,1:n2) )
  allocate( lcount(0:count) )
  lcount(:) = 0
  do j=1,n2
  do i=1,n1
    count = grid%guess_count(i,j)
    lcount(count) = lcount(count) + 1
  enddo
  enddo



  csum = zero
  do i=1,ubound(lcount,1)
    csum = csum + dble(i)*dble(lcount(i))
  enddo
  cavg = csum/dble(n1*n2-lcount(0))




  print*,'total number of boxes,n1,n2,n1*n2 ',n1,n2,n1*n2
  print*,'number of nonzero boxes ', lcount(0)
  print*,'average number ', dble(sum(lcount))/(n1*n2-lcount(0))
  print*,'max number ', ubound(lcount,1)

  do i=0,ubound(lcount,1)
    print*,'i,lcount(i) ',i,lcount(i)
  enddo

  deallocate( lcount )

!debug-begin
! write out coordinates of triangular mesh
  if (idebug.ge.2) then 
  open(iodev,file='node.txt',access='sequential',form='formatted')
  rewind(iodev)
  do i=1,ubound(grid%x,2)
     write(iodev,*) grid%x(1,i), grid%x(2,i)
  enddo
  close(iodev)
  endif
!debug-end



  return



end subroutine check_guess_table
