!       compute the triangle and barycentric coordinates for a list of particles
!
subroutine search_vtr( grid, num, xy, tr_save, p_save )
  use grid_class
  implicit none
  type(grid_type) :: grid
  integer :: num
  real(kind=8) :: xy(2,num)
  integer :: tr_save(num)
  real(kind=8) :: p_save(3,num)



  logical :: use_preprocess 
  integer :: i
  integer, parameter :: iodev = 19

  !        ---------------------------------------------
  !        if there are not too many particles,
  !        no precessing may be faster since the high
  !        cost associated with preprocessing cannot be
  !        amortized by many particles
  !        ---------------------------------------------

  use_preprocess = .false.
  if (use_preprocess) then
     call search_vtr_preprocess(grid,num,xy,tr_save,p_save)
  else
     call search_vtr_nopreprocess(grid,num,xy,tr_save,p_save)
  endif

  return
end subroutine search_vtr



subroutine search_vtr_nopreprocess( grid, num, xy, tr_save, p_save )
  use grid_class
  implicit none
  type(grid_type) :: grid
  integer :: num
  real(kind=8) :: xy(2,num)
  integer :: tr_save(num)
  real(kind=8) :: p_save(3,num)
  !
  integer :: ilo,ihi,jlo,jhi,ipart,inode,ncount
  integer :: ii,jj,i,j,it,itr,istart,iend
  real(kind=8) :: dx,dy, xx,yy, xx3,yy3
  real(kind=8) :: map11,map12,map21,map22
  real(kind=8) :: p(1:3), p1,p2,p3, t1,t2
  real(kind=8) :: guess_min(2), inv_guess_d(2) 
  logical :: is_found, is_ok

  real(kind=8), parameter :: one = 1.0d0
  real(kind=8), parameter :: eps = 20.0d0*epsilon(one)
  real(kind=8), parameter :: tol = 1.0d-6

  integer, parameter :: idebug = 2

  ncount = 0
  if (idebug.ge.1) then
     call cpu_time(t1)
  endif

  ilo = lbound( grid%guess_table, 1 )
  ihi = ubound( grid%guess_table, 1 )

  jlo = lbound( grid%guess_table, 2 )
  jhi = ubound( grid%guess_table, 2 )

  guess_min(1) = grid%guess_min(1)
  guess_min(2) = grid%guess_min(2)
  inv_guess_d(1) = grid%inv_guess_d(1)
  inv_guess_d(2) = grid%inv_guess_d(2)

  do ipart=1,num
     xx = xy(1,ipart)
     yy = xy(2,ipart)
     ii = (xx - guess_min(1))*inv_guess_d(1) + 1
     jj = (yy - guess_min(2))*inv_guess_d(2) + 1
     i = max(ilo, min(ihi, ii ) )
     j = max(jlo, min(jhi, jj ) )

     tr_save(ipart) = -1


     istart = grid%guess_xtable(i,j)
     iend = istart + grid%guess_count(i,j) - 1
     do it=istart,iend

        if (idebug.ge.1) then
           ncount = ncount + 1
        endif

        itr = grid%guess_list(it)

        map11 = grid%mapping(1,1,itr)
        map12 = grid%mapping(1,2,itr)
        map21 = grid%mapping(2,1,itr)
        map22 = grid%mapping(2,2,itr)

        xx3 = grid%mapping(1,3,itr)
        yy3 = grid%mapping(2,3,itr)

        dx = xx - xx3
        dy = yy - yy3

        p1 = map11 * dx + map12 * dy
        p2 = map21 * dx + map22 * dy
        p3 = one - p1 - p2

        is_found = (p1 .ge. -eps).and.                                  &
             &                  (p2 .ge. -eps).and.                                  &
             &                  (p3 .ge. -eps)
        if (is_found) then
           tr_save(ipart) = itr
           p_save(1,ipart) = p1
           p_save(2,ipart) = p2
           p_save(3,ipart) = p3
           exit
        endif

     enddo ! end do it

  enddo ! end do ipart

  if (idebug.ge.1) then
     call cpu_time(t2)
     write(*,*) 'search_vtr_nopre: num, ncount,time ',num,ncount,t2-t1
  endif

  if (idebug.ge.2) then
     call cpu_time(t1)
     do ipart=1,num
        call search_tr(grid, xy(1:2,ipart),itr,p)
        is_ok = (itr .eq. tr_save(ipart))
        is_ok = is_ok .and.                                               &
             &          (maxval( abs(p(1:3)-p_save(1:3,ipart)) ).lt.tol)
        if (.not.is_ok) then
           write(*,*) 'search_vtr_nopre: itr,tr_save',itr,tr_save(ipart)
           write(*,*) 'p(1:3) ',p(1:3)
           write(*,*) 'p_save ',p_save(1:3,ipart)
           stop '** error in search_vtr '
        endif
     enddo
     call cpu_time(t2)
     write(*,*) 'search_vtr_nopr: old time ',t2-t1
  endif


  return
end subroutine search_vtr_nopreprocess

subroutine search_vtr_preprocess( grid, num, xy, tr_save, p_save )
  use grid_class
  implicit none
  type(grid_type) :: grid
  integer :: num
  real(kind=8) :: xy(2,num)
  integer :: tr_save(num)
  real(kind=8) :: p_save(3,num)

  integer :: ilo,ihi, jlo, jhi
  integer, dimension(:,:), allocatable :: icount, xtable
  integer, dimension(:), allocatable :: table

  integer :: i,j,k, ip, ii, jj, j2, i2, jstart,jend,jsize
  integer :: istart,iend, ipart
  integer :: it, imax, itr,  inode, ierr


  real(kind=8) :: map(1:2,1:2), xy3(1:2)
  integer, dimension(:), allocatable :: jp, jsum

  real(kind=8) :: guess_min(2), inv_guess_d(2) 
  real(kind=8) :: xy123(2,3),p123(3)

  logical :: is_in, is_ok

  logical :: use_matmul = .false.
  integer, parameter :: idebug = 2
  integer :: itriangle,jneed
  integer, allocatable, dimension(:) :: jfound
  real(kind=8) :: p1,p2,p3,  dx1, dx2, p(3)
  logical :: is_computed

  integer :: ncount
  real(kind=8) :: t1,t2

  real(kind=8), parameter :: one = 1.0d0
  real(kind=8), parameter :: eps = 20.0d0*epsilon(one)
  real(kind=8), parameter :: tol = 1.0d-6

  !       --------------------------------------------------
  !       1st pass to count number of particles in each box
  !       --------------------------------------------------
  if (idebug.ge.1) then
     ncount = 0
     call cpu_time(t1)
  endif


  do ipart=1,num
     tr_save(ipart) = -1
  enddo

  ilo = lbound( grid%guess_table, 1 )
  ihi = ubound( grid%guess_table, 1 )

  jlo = lbound( grid%guess_table, 2 )
  jhi = ubound( grid%guess_table, 2 )

  allocate( icount(ilo:ihi,jlo:jhi),                                  &
       &       xtable(ilo:ihi,jlo:jhi), table(num),                           &
       &       jp(jlo:jhi), jsum(jlo:jhi), stat=ierr )
  call assert( ierr.eq.0,'search_vtr: allocate(icount)',ierr)

  if (idebug.ge.1) then
     write(*,*) 'search_vtr:num,isize,jsize ',                         &
          &               num,(ihi-ilo+1),(jhi-jlo+1)
  endif

  do j=jlo,jhi
     do i=ilo,ihi
        icount(i,j) = 0
     enddo
  enddo

  guess_min(1) = grid%guess_min(1)
  guess_min(2) = grid%guess_min(2)
  inv_guess_d(1) = grid%inv_guess_d(1)
  inv_guess_d(2) = grid%inv_guess_d(2)

  do ipart=1,num
     ii = (xy(1,ipart) - guess_min(1))*inv_guess_d(1) + 1
     jj = (xy(2,ipart) - guess_min(2))*inv_guess_d(2) + 1
     i = max(ilo, min(ihi, ii ) )
     j = max(jlo, min(jhi, jj ) )

     icount(i,j) = icount(i,j) + 1
  enddo

  do j=jlo,jhi
     jsum(j) = 0
  enddo

  imax = 0
  do j=jlo,jhi
     do i=ilo,ihi
        imax = max( icount(i,j), imax )
        jsum(j) = jsum(j) + icount(i,j)
     enddo
  enddo

  allocate( jfound(max(1,imax)), stat=ierr)
  call assert(ierr.eq.0,'search_vtr: allocate(jfound)',ierr)


  !       --------------
  !       setup pointers
  !       --------------
  ip = 1
  do j=jlo,jhi
     jp(j) = ip
     ip = ip + jsum(j)
  enddo


  do j=jlo,jhi
     do i=ilo,ihi
        xtable(i,j) = jp(j)
        jp(j) = jp(j) + icount(i,j)
     enddo
  enddo

  !       ----------------------------------------
  !       2nd pass to deposit particles into boxes
  !       ----------------------------------------
  do j=jlo,jhi
     do i=ilo,ihi
        icount(i,j) = 0
     enddo
  enddo

  do k=1,num

     ii = (xy(1,k) - guess_min(1))*inv_guess_d(1) + 1
     jj = (xy(2,k) - guess_min(2))*inv_guess_d(2) + 1
     i = max(ilo, min(ihi, ii ) )
     j = max(jlo, min(jhi, jj ) )

     ip = xtable(i,j)  + icount(i,j)
     table(ip) = k
     icount(i,j) = icount(i,j) + 1
  enddo

  !       -------------------------------------------
  !       for each box, compute the triangle number 
  !       and barycentric coordinates
  !       -------------------------------------------
  do j=jlo,jhi
     do i=ilo,ihi
        if (icount(i,j).le.0) cycle

        jneed = icount(i,j)

        istart = grid%guess_xtable(i,j)
        iend = istart + grid%guess_count(i,j) - 1

        jstart = xtable(i,j)
        jend = jstart + icount(i,j)-1
        jsize = jend - jstart + 1


        jfound(1:jsize) = -1

        do it=istart,iend
           if (jneed.eq.0) exit

           itr = grid%guess_list(it)

           map(1,1) = grid%mapping(1,1,itr)
           map(2,1) = grid%mapping(2,1,itr)
           map(1,2) = grid%mapping(1,2,itr)
           map(2,2) = grid%mapping(2,2,itr)

           xy3(1) = grid%mapping(1,3,itr)
           xy3(2) = grid%mapping(2,3,itr)




           !  -------------------------------------------------
           !  want (0 <= p1 <= 1) & (0 <= p2 <= 1) & (0 <= p3 <= 1)
           !  where p3 = 1 - p1  - p2
           !  note that
           !  if p1 >= 0 and p2 >= 0 and 1 >= p1 + p2
           !  then it implies 1 >= p1 and 1 >= p2 as well
           !  1 >= p1 + p2 is equivalent to 1 - p1 - p2 >= 0
           !  or p3 >= 0
           !  1 >= p3 is equivalent to 1 >= 1 - p1 - p2
           !  or 0 >= -p1 - p2 or 0 <= p1 + p2, which is true 
           !  since we assumed p1 >= 0 and p2 >= 0 
           !
           !  thus we need only 3 tests instead of 6 tests
           !  -------------------------------------------------
           do jj=1,jsize
              is_computed = (jfound(jj).eq.itr)
              if (is_computed) cycle

              if (idebug.ge.1) then
                 ncount = ncount + 1
              endif

              ipart = table( (jstart-1) + jj)
              dx1 = xy(1, ipart) - xy3(1)
              dx2 = xy(2, ipart) - xy3(2)


              p1 = map(1,1)*dx1 + map(1,2)*dx2
              p2 = map(2,1)*dx1 + map(2,2)*dx2
              p3 = one - p1  - p2


              is_in = (p1.ge.-eps).and.                              &
                   &                   (p2.ge.-eps).and.                              &
                   &                   (p3.ge.-eps)
              if (is_in) then


                 tr_save(ipart) =  itr
                 p_save(1,ipart) = p1
                 p_save(2,ipart) = p2
                 p_save(3,ipart) = p3

                 jfound(jj) = itr
                 jneed = jneed - 1

                 if (idebug.ge.3) then
                    call search_tr(grid,xy(1:2,ipart),itriangle,p)
                    is_ok = tr_save(ipart) .eq.itriangle
                    is_ok = is_ok.and.                                     &
                         &                  maxval(abs(p_save(1:3,ipart)-p(1:3))).le.tol 
                    if (.not.is_ok) then
                       write(*,*) 'ipart,itr, itriangle ',itr,itriangle 
                       write(*,*) 'p(:) ',p(1:3)
                       write(*,*) 'p_save(:) ',p_save(1:3,ipart)
                       call t_coeff(grid,itr,xy(1:2,ipart),p)
                       write(*,*) 't_coeff,itr ',p(1:3)
                       call t_coeff(grid,itriangle,xy(1:2,ipart),p)
                       write(*,*) 't_coeff,itriangle ',p(1:3)
                       stop 'stop for debugging '
                    endif
                 endif

              endif
           enddo ! do jj

        enddo ! do it

        !       ------------------------------------------------
        !       double check, each point should be in a triangle
        !       ------------------------------------------------
        if (idebug.ge.3) then
           do jj=1,jsize
              ipart = table( (jstart-1) + jj)
              itr = tr_save(ipart)
              is_ok = (itr .ge.  1)
              if (is_ok) then
                 call t_coeff(grid,itr,xy(1:2,ipart),p)
                 is_ok = (p(1) .ge. -eps).and.                                   &
                      (p(2) .ge. -eps).and.                                   &
                      (p(3) .ge. -eps)
              endif

              if (.not. is_ok) then
                 write(*,*) 'ipart, itr ',ipart,itr
                 write(*,*) 'xy(:) ', xy(1,ipart),xy(2,ipart)

                 do it=istart,iend
                    itr  = grid%guess_list(it)

                    do i2=1,3
                       inode = grid%nd(i2,itr)
                       xy123(1:2,i2) = grid%x(1:2,inode)
                    enddo

                    call t_coeff( grid, itr, xy(1:2,ipart), p123)
                    write(*,9010) itr,p123
9010                format(1x,'itr ',i7,'p123(:) ',3(1pe20.10,1x))

                 enddo
              endif
              call assert(is_ok,'search_vtr: invalid itr, ipart=',ipart)
           enddo
        endif

     enddo ! do i=ilo,ihi
  enddo ! do j=jlo,jhi


  deallocate( table, xtable, icount, jsum, jp, stat=ierr )
  call assert( ierr.eq.0,'search_vtr: deallocate(table)',ierr)


  deallocate( jfound, stat=ierr )
  call assert( ierr.eq.0,'search_vtr: deallocat(jfound)',ierr)

  if (idebug.ge.1) then
     call cpu_time(t2)
     write(*,*) 'search_vtr:num,ncount,time ',num,ncount,t2-t1
  endif



  if (idebug.ge.2) then
     call cpu_time(t1)
     do ipart=1,num
        call search_tr(grid, xy(1:2,ipart),itr,p)
        is_ok = (itr .eq. tr_save(ipart))
        is_ok = is_ok .and.                                               &
             &          (maxval( abs(p(1:3)-p_save(1:3,ipart)) ).lt.tol)
        if (.not.is_ok) then
           write(*,*) 'search_vtr: itr,tr_save',itr,tr_save(ipart)
           write(*,*) 'p(1:3) ',p(1:3)
           write(*,*) 'p_save ',p_save(1:3,ipart)
           stop '** error in search_vtr '
        endif
     enddo
     call cpu_time(t2)
     write(*,*) 'search_vtr_pre: old time ',t2-t1
  endif



  return
end subroutine search_vtr_preprocess


subroutine assert(lcond,msg,ival)
  implicit none
  logical lcond
  character(len=*) msg
  integer ival

  if (.not.lcond) then
     write(*,*) msg,ival
     stop '** assertion error ** '
  endif

  return
end subroutine assert
