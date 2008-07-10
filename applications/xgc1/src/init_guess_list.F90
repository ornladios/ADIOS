subroutine init_guess_list( grid )
  use grid_class
  implicit none
  type(grid_type) :: grid



  real(kind=8), parameter :: eps = 1.0d-10
  real(kind=8), dimension(2) :: xy, xy1,xy2,xy3,  xymin,xymax
  real(kind=8) :: p1,p2,p3, pmin, pmax, dp
  real(kind=8), dimension(3) :: p
  real(kind=8) :: xmin, xmax, ymin, ymax
  real(kind=8) :: xlo,xhi, ylo, yhi

  integer :: n(2), ilo,ihi,jlo,jhi
  integer, pointer, dimension(:,:) :: guess_count, guess_xtable
  integer, pointer, dimension(:) :: guess_list, nhits 
  real(kind=8), pointer, dimension(:) :: dhits
  integer :: ierr

  logical, parameter :: use_dhits = .true.
  logical :: need_swap

  integer, dimension(2) :: ijmin, ijmax, ij
  integer :: ii,jj, i,j,k,kk,   i1,i2,j1,j2, ntotal,nmax
  integer :: ifree, ip, istart,iend, itemp, ipos, npoint
  integer :: itr

  integer :: icount
  integer, allocatable, dimension(:) :: ilist
  real(kind=8) :: txy(1:2,0:2)
  logical :: has_overlap

  logical :: is_found, is_ok
  logical :: need_check_overlap = .false.

  real(kind=8) :: dtemp, darea
  real(kind=8) :: triangle_area, x1,y1,x2,y2,x3,y3
  triangle_area(x1,y1, x2,y2, x3,y3) = ((x2*y3-x3*y2) -                   &
               (x1*y3-x3*y1) + (x1*y1-x2*y2) )/2.0d0


  n = grid%guess_n
  grid%guess_d = (grid%guess_max - grid%guess_min)/n
  grid%inv_guess_d = 1D0/grid%guess_d


  allocate( guess_xtable( n(1), n(2) ) )
  allocate( guess_count(  n(1), n(2) ) )

  ilo = 1
  ihi = n(1)
  jlo = 1
  jhi = n(2)

  guess_count( :, :) = 0

  !       ---------------------------------------------------
  !       count number of (potential) triangles in each cell
  !       ---------------------------------------------------
  do itr=1,grid%ntriangle

     xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
     xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
     xy3(1:2)= grid%x(1:2,grid%nd(3,itr))

     !       -----------------------------------------
     !       determine the bottom left and upper right
     !       corners of enclosing rectangle
     !       -----------------------------------------
     xmin = min( xy1(1), min(xy2(1), xy3(1)) )
     xmax = max( xy1(1), max(xy2(1), xy3(1)) )
     ymin = min( xy1(2), min(xy2(2), xy3(2)) )
     ymax = max( xy1(2), max(xy2(2), xy3(2)) )


     xymin(1) = xmin
     xymin(2) = ymin
     ijmin(1:2)= (xymin(1:2)-grid%guess_min(1:2))*grid%inv_guess_d +1

     xymax(1) = xmax
     xymax(2) = ymax
     ijmax(1:2)= (xymax(1:2)-grid%guess_min(1:2))*grid%inv_guess_d +1


     i1 = max(ilo,min(ihi, ijmin(1)))
     i2 = max(ilo,min(ihi, ijmax(1)))
     j1 = max(jlo,min(jhi, ijmin(2)))
     j2 = max(jlo,min(jhi, ijmax(2)))


     do j=j1,j2
        do i=i1,i2
           guess_count(i,j) = guess_count(i,j) + 1
        enddo
     enddo

  enddo

  !       ---------------------
  !       prepare guess_list(:)
  !       ---------------------
  ntotal = 0
  nmax = 0
  do j=jlo,jhi
     do i=ilo,ihi
        ntotal = ntotal + guess_count(i,j)
        nmax = max( nmax, guess_count(i,j) )
     enddo
  enddo

  allocate( guess_list( ntotal ) )
  ifree = 1

  do j=jlo,jhi
     do i=ilo,ihi
        guess_xtable( i,j ) = ifree
        ifree = ifree + guess_count(i,j)
     enddo
  enddo

  !       ----------------------------
  !       copy the triangles into list
  !       ----------------------------

  guess_count(:,:) = 0

  do itr=1,grid%ntriangle

     xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
     xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
     xy3(1:2)= grid%x(1:2,grid%nd(3,itr))

     !       -----------------------------------------
     !       determine the bottom left and upper right
     !       corners of enclosing rectangle
     !       -----------------------------------------
     xmin = min( xy1(1), min(xy2(1), xy3(1)) )
     xmax = max( xy1(1), max(xy2(1), xy3(1)) )
     ymin = min( xy1(2), min(xy2(2), xy3(2)) )
     ymax = max( xy1(2), max(xy2(2), xy3(2)) )


     xymin(1) = xmin
     xymin(2) = ymin
     ijmin(1:2)= (xymin(1:2)-grid%guess_min(1:2))*grid%inv_guess_d +1

     xymax(1) = xmax
     xymax(2) = ymax
     ijmax(1:2)= (xymax(1:2)-grid%guess_min(1:2))*grid%inv_guess_d +1


     i1 = max(ilo,min(ihi, ijmin(1)))
     i2 = max(ilo,min(ihi, ijmax(1)))
     j1 = max(jlo,min(jhi, ijmin(2)))
     j2 = max(jlo,min(jhi, ijmax(2)))


     do j=j1,j2
        do i=i1,i2
           ip = guess_xtable(i,j) + guess_count(i,j)
           guess_list(ip) = itr

           guess_count(i,j) = guess_count(i,j) + 1
        enddo
     enddo

  enddo

  !debug
!  write(*,*) 'ntotal, nmax, navg ', ntotal,nmax,                          &
!       &            dble(ntotal)/dble(n(1)*n(2))

  !       --------------------------------
  !       sort triangles by number of hits
  !       --------------------------------
  allocate( nhits( ntotal ), dhits(ntotal),stat=ierr )
  call assert( ierr.eq.0,'allocate nhits, ntotal=',ntotal)
  nhits(:) = 0
  dhits(:) = 0.0d0

  icount = maxval( guess_count(ilo:ihi,jlo:jhi) )
  allocate( ilist( icount + 1 ) )

  npoint = 4

    do itr=1,grid%ntriangle

     xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
     xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
     xy3(1:2)= grid%x(1:2,grid%nd(3,itr))

     darea = triangle_area(xy1(1),xy1(2),xy2(1),xy2(2),xy3(1),xy3(2))
     dp = 1.0d0/npoint

     do ii=0,npoint
     do jj=0,npoint-ii
           p1 = dble(ii)*dp
           p2 = dble(jj)*dp
           p3 = 1.0d0 - p(1) - p(2)

           xy(1) = p1*xy1(1) + p2*xy2(1) + p3*xy3(1)
           xy(2) = p1*xy1(2) + p2*xy2(2) + p3*xy3(2)



           ij = (xy-grid%guess_min)*grid%inv_guess_d + 1
           i = max( ilo, min(ihi, ij(1) ) )
           j = max( jlo, min(jhi, ij(2) ) )

        

           istart = guess_xtable(i,j)
           iend = istart + guess_count(i,j)-1
           do k=istart,iend
              is_found = (itr.eq.guess_list(k))
              if (is_found) then
                 nhits(k) = nhits(k) + 1
                 dhits(k) = dhits(k) + darea
                 exit
              endif
           enddo

        enddo ! do ii
        enddo ! do jj

  enddo ! do itrig

  !       ------------------------------------------
  !       move triangles with more hits to the front
  !       ------------------------------------------
  do j=jlo,jhi
  do i=ilo,ihi

        if (guess_count(i,j).le.0) cycle


        xlo = grid%guess_min(1) + (i-1)*grid%inv_guess_d(1)
        xhi = grid%guess_min(1) + (i)*grid%inv_guess_d(1)

        ylo = grid%guess_min(2) + (j-1)*grid%inv_guess_d(2)
        yhi = grid%guess_min(2) + (j)*grid%inv_guess_d(2)

        istart = guess_xtable(i,j)
        iend = istart + guess_count(i,j)-1

        !       -------------------------------------
        !       simple bubble sort, ok for short list
        !       -------------------------------------
        do k=istart,iend
           do kk=k+1,iend
            if (use_dhits) then
              need_swap =  dhits(kk).gt.dhits(k)
            else
              need_swap = nhits(kk).gt.nhits(k)
            endif

            if(need_swap) then
                 ! --------------
                 ! swap data
                 ! --------------
                 itemp = nhits( k )
                 nhits(k) = nhits(kk)
                 nhits(kk) = itemp

                 dtemp = dhits(k)
                 dhits(k) = dhits(kk)
                 dhits(kk) = dtemp
           

                 itemp = guess_list( k )
                 guess_list(k) = guess_list(kk)
                 guess_list(kk) = itemp

              endif
           enddo
        enddo

        !        ---------------------------
        !        double check list is sorted
        !        ---------------------------

         do k=istart,iend-1
          if (use_dhits) then
           is_ok = dhits(k) .ge. dhits(k+1)
          else
           is_ok = nhits(k) .ge. nhits(k+1)
          endif

          if (.not.is_ok) then
              write(*,*) 'guess_list(:) ',guess_list(istart:iend)
              write(*,*) 'dhits(:) ',dhits(istart:iend)
              write(*,*) 'nhits(:) ',nhits(istart:iend)
              stop '** error in init_guess_list ** '
           endif
         enddo


!      --------------------------------------------
!      scrutinize and remove triangles with no overlap
!      --------------------------------------------
       if (need_check_overlap) then

       icount = 0
       do k=istart,iend
          itr =  guess_list(k)

          if (nhits(k).ge.1) then
             icount = icount + 1
             ilist(icount) = guess_list(k)
          else

           xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
           xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
           xy3(1:2)= grid%x(1:2,grid%nd(3,itr))

           txy(1:2,0) = xy1(1:2)
           txy(1:2,1) = xy2(1:2)
           txy(1:2,2) = xy3(1:2)

            call checkoverlap( txy, xlo, xhi, ylo, yhi, has_overlap) 
            if (has_overlap) then
               icount = icount + 1
               ilist(icount) = guess_list(k)
            endif
          endif
       enddo

!      -----------------
!      copy results back
!      -----------------
       guess_list(istart:(istart+icount-1)) = ilist(1:icount)
       guess_count(i,j) = icount

      endif

     enddo
  enddo


  deallocate( ilist )
  deallocate( nhits, dhits )

  !       --------------------------------
  !       assign back to grid data structure
  !       --------------------------------
  grid%guess_list => guess_list
  grid%guess_xtable => guess_xtable
  grid%guess_count => guess_count


!debug-begin
!   call check_guess_table(grid)
!debug-end

  return
end subroutine  init_guess_list
