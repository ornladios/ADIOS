module grid_class
  use boundary_class
  type grid_type
     integer :: nnode      !! number of nodes in a grid system
     integer :: ntriangle  !! number of trianble in a grid system

     !for nodes
     integer, pointer :: nn(:,:)  !! next node with (1) B-field direction and (2) anti B-field direction
     integer, pointer :: p(:)     !! property value (boundary, ... etc)
     integer, pointer :: rgn(:)   !! region value for each node point
     integer, pointer :: it(:)  !! poloidal angle indext     
     real (kind=8), pointer :: psi(:)  !! psi value of node point (mabye need to change to a integer value)     
     real (kind=8), pointer :: x(:,:)  !! R-Z Positon of each node point

     !for triangle
     integer, pointer :: nd(:,:)  !! 3 node numbers that each triangle has
     integer, pointer :: adj(:,:) !! 3 adjacent triangle

     real (kind=8), pointer :: mapping(:,:,:)  !! shape function coefficient

     !for node->triangle structure
     integer, pointer :: tr_node(:,:)  !! triangle index that has a given node point 
     integer, pointer :: num_t_node(:) !! number of trangle for each node, triangle index for each node

     !guess table    
     integer :: guess_n(2)             !! number of array size that
     real (kind=8) :: guess_min(2), guess_max(2), guess_d(2), inv_guess_d(2)  !! min, max, delta for guess table in RZ space
     integer, pointer :: guess_table(:,:)        !! guess value for triangle search

     !phi direction
     integer :: iphi_offset  !! toroidal angle indexing offset
     integer :: nphi      !! number of toroidal cross section of one PE
     real (kind=8) :: delta_phi,phimin,phimax       !! Delta toroidal angle 

     
     ! simple boundary sepcification
!     integer :: inner_bd(4),outer_bd(8)
     type(boundary_type) :: bd

     !psi average
     integer :: npsi
     integer, pointer ::ntheta(:), itheta0(:)

#ifdef CIRCULAR_SPECIAL
     ! cc find index
     real (kind=8), pointer :: radial(:) ! nfl
     real (kind=8), pointer :: theta(:) ! nnode
#endif

     !node volume 
     real (kind=8),pointer ::  surf_vol(:), inv_node_vol(:),node_area(:),node_vol(:)
     !triangle volume
     real (kind=8),pointer ::  tr_vol(:),tr_area(:)
     !temporory variables for calculations - number of node point
     real (kind=8), pointer  :: rtmp1(:),rtmp2(:)
#ifdef XGC_DEBUG10
     real (kind=8), pointer  :: rtmp3(:),rtmp4(:),rtmp5(:),rtmp6(:),rtmp7(:),rtmp8(:)
#endif
     
     integer, pointer :: old2new(:)
!
     integer, pointer, dimension(:) :: guess_list
     integer, pointer, dimension(:,:) :: guess_xtable, guess_count

     real (kind=8) , pointer :: bfield(:,:)
!     integer, pointer :: sort(:), inv_sort(:)   ! sort(original)=sorted number , inv_sort(sorted_number) = original. sort(axis)=1, inv_sort(1)=axis
  end type grid_type

contains
  !! read grid information from file
  subroutine init_grid(grid,filename1,filename2)
    use sml_module
    implicit none
    include 'mpif.h'
    type(grid_type) :: grid
    character (len=65), intent(in) :: filename1,filename2 !! (1) node file, (2) triangle element file
    integer, parameter :: nodefile=2001,trifile=2002
    integer :: i,j,k,dum,n_n,n_t,ierr
    integer :: mype,total_pe
    integer, pointer :: ipsi(:),bd_flag(:),old_node_num(:)
    real (kind=8) , external :: psi_interpol
    real (kind=8) :: psi_val
    integer :: old_node_max, si




    if(sml_mype==0) then
       !read node -----------------------------------------------------------------
       !open file
       open(nodefile, file=filename1,status='old',form='formatted')
       open(trifile, file=filename2,status='old',form='formatted')
       
       !read
       ! read number of nodes
       read(nodefile,*) grid%nnode, dum,dum,dum
    endif

    call MPI_BCAST(grid%nnode, 1,  MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    
    n_n=grid%nnode
    ! read nodes
    allocate(grid%x(2,n_n), grid%nn(2,n_n),grid%p(n_n),grid%rgn(n_n),grid%it(n_n),grid%psi(n_n)) 
    allocate(ipsi(n_n),bd_flag(n_n),old_node_num(n_n))
    allocate(grid%node_area(n_n),grid%inv_node_vol(n_n),grid%node_vol(n_n))
    allocate(grid%bfield(4,n_n))
!    allocate(grid%sort(n_n),grid%inv_sort(n_n))
    if(sml_mype==0) then
       do i=1, n_n
          read(nodefile,*) dum, grid%x(1,i),grid%x(2,i),ipsi(i),grid%p(i),grid%nn(1,i),&
               grid%nn(2,i),grid%rgn(i),grid%it(i),old_node_num(i),bd_flag(i)
          !       print *, dum,grid%x(1,i),grid%x(2,i),ipsi(i),grid%p(i),grid%nn(1,i),grid%nn(2,i),grid%rgn(i),grid%it(i)
       enddo       
    endif
    !Broad Casting data
    call MPI_BCAST(grid%x, n_n*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ipsi, n_n, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(grid%p, n_n, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(grid%nn, n_n*2, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(grid%rgn, n_n, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(grid%it, n_n, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(old_node_num, n_n, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(bd_flag, n_n, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    ! end of reading and broadcasting -------------------------------------------------
           
    old_node_max=maxval(old_node_num)
    if(sml_mype==0) print *, 'old_node_max =', old_node_max, 'grid%nnode =', grid%nnode
    if(old_node_max>=2*n_n) print *, 'Warning : Too much difference of old_node index. Check node file.'
    
    !    if(minval(old_node_num)<=0) then 
    !    print *, 'Error : wrong old_node index', minval(old_node_num),minloc(old_node_num),&
    !         old_node_max,maxloc(old_node_num)
    !    stop
    !    endif
    allocate(grid%old2new(old_node_max))
    grid%old2new=0
    do i=1, n_n
       grid%old2new(old_node_num(i))=i
    enddo
    
    do i=1, n_n
       do j=1, 2
          if(grid%nn(j,i)>0) grid%nn(j,i)=grid%old2new(grid%nn(j,i))
       enddo
    enddo
    
    allocate(grid%rtmp1(n_n),grid%rtmp2(n_n))
#ifdef XGC_DEBUG10
    allocate(grid%rtmp3(n_n),grid%rtmp4(n_n))
    allocate(grid%rtmp5(n_n),grid%rtmp6(n_n))
    allocate(grid%rtmp7(n_n),grid%rtmp8(n_n))
!    rtmp3=0D0
!    rtmp4=0D0
!    rtmp5=0D0
!    rtmp6=0D0
!    rtmp7=0D0
    grid%rtmp8=0D0
#endif    
    
    ! get sorted list

!    do i=1, n_n
!       grid%inv_sort(i)=i
!    enddo
    !**
!    call qsort(1,n_n) ! returns sorted list grid%sort
    
!    do i=1, n_n
!       grid%sort(grid%inv_sort(i))=i
!    enddo

    !set psi
!    si=1  !sorted i
    i=1 !grid%inv_sort(si)
    grid%npsi=1
    psi_val=psi_interpol(grid%x(1,i),grid%x(2,i),0,0)
    grid%psi(i)=psi_val
    do si=2, n_n
       i=si !grid%inv_sort(si)
       if(ipsi(i)/=ipsi(si-1))  then
          psi_val=psi_interpol(grid%x(1,i),grid%x(2,i),0,0)
          grid%npsi=grid%npsi+1
       endif
       if(ipsi(i)==10000) then !set psi value for wall boundary
          psi_val=psi_interpol(grid%x(1,i),grid%x(2,i),0,0)
       endif
       grid%psi(i)=psi_val
    enddo


    !set ntheta, itheta0
    allocate(grid%ntheta(grid%npsi),grid%itheta0(grid%npsi))

    k=1
    grid%ntheta(1)=1
    grid%itheta0(1)=1
    do si=2, n_n
       i=si !grid%inv_sort(si)
       if(ipsi(i)/=ipsi(si-1)) then
          k=k+1
          grid%ntheta(k)=1
          grid%itheta0(k)=si
       else
          grid%ntheta(k)=grid%ntheta(k)+1
       endif
    enddo


    ! read number of triangles - element file
    if(sml_mype==0) then
       read(trifile,*) grid%ntriangle,dum,dum
    endif
    !Broad Casting data
    call MPI_BCAST(grid%ntriangle, 1,  MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

    n_t=grid%ntriangle
    ! number of triangle
!    --------------------------------------------------------
!    use the column grid%mapping(1:2,3,itr) in mapping array to hold the 
!    coordinate for 3rd vertex in triangle itr
!    --------------------------------------------------------
    allocate( grid%nd(3,n_t), grid%adj(3,n_t),grid%mapping(2,3,n_t))
    allocate( grid%tr_vol(n_t),grid%tr_area(n_t) )
    !read from file
    if(sml_mype==0) then
       do i=1, grid%ntriangle
          read(trifile,*) dum,grid%nd(1,i),grid%nd(2,i),grid%nd(3,i)
       enddo
    endif
    
    call MPI_BCAST(grid%nd, n_t*3, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    ! end of reading and broadcasting ----------------------------------------------

    call init_triangle(grid)

    call get_pe_info(mype,total_pe)
    call get_nphi(grid%nphi)
    grid%delta_phi=atan(1.D0)*8D0/real(sml_nphi_total)
    grid%iphi_offset=grid%nphi* (mype/sml_pe_per_plane)
    grid%phimin=real(grid%iphi_offset)*grid%delta_phi
    grid%phimax=grid%phimin + grid%delta_phi*real(grid%nphi)
    !    print *, mype, grid%iphi_offset, grid%phimin/atan(1D0), grid%phimax/atan(1D0)

    i=1
    do while(bd_flag(i)==1)
       i=i+1
    enddo
!    grid%inner_bd(1)=1    
!    grid%inner_bd(2)=i-1
    grid%bd%in%start=1
    grid%bd%in%end=i-1

    i=n_n
    do while(bd_flag(i)==1)
       i=i-1
    enddo
!    grid%outer_bd(1)=i+1
!    grid%outer_bd(2)=n_n
!    grid%outer_bd(3)=n_n
!    grid%outer_bd(4)=n_n
    grid%bd%out1%start=i+1
    grid%bd%out1%end  =n_n
    grid%bd%out2%start=n_n
    grid%bd%out2%end  =n_n



!    grid%inner_bd(3)=grid%inner_bd(1)
!    grid%inner_bd(4)=grid%inner_bd(2)
!    grid%outer_bd(5)=grid%outer_bd(1)
!    grid%outer_bd(6)=grid%outer_bd(2)
!    grid%outer_bd(7)=grid%outer_bd(3)
!    grid%outer_bd(8)=grid%outer_bd(4)

    deallocate(ipsi,bd_flag,old_node_num)

    call get_node_vol(grid)


    do i=1, n_n 
       call bvec_interpol(grid%x(1,i),grid%x(2,i),0D0,grid%bfield(1,i),grid%bfield(2,i),grid%bfield(3,i))
       grid%bfield(4,i)=sqrt(sum(grid%bfield(:,i)**2))
    enddo

#ifdef CIRCULAR_SPECIAL
    allocate(grid%radial(grid%npsi), grid%theta(n_n))
    call get_cc_values(grid)
#endif
!!$  contains
!!$    recursive subroutine qsort(left,right)
!!$      implicit none
!!$      integer :: left,right
!!$      integer :: i,j,p, store
!!$      !initialize
!!$      p=(left+right)/2
!!$      if(right<=left) return
!!$      
!!$      ! partition
!!$      call sort_swap(right,p) ! move pivot to right position
!!$      p=right                 ! indicate pivot
!!$      store=left-1            ! set to zero
!!$      do i=left, right-1
!!$         if(compare(i,p)<0) then
!!$            store=store+1
!!$            call sort_swap(store,i)
!!$         endif
!!$      enddo
!!$      call sort_swap(p,store+1)
!!$      p=store+1
!!$      !//resursion
!!$      call qsort(left,p-1)
!!$      call qsort(p+1,right)
!!$    end subroutine qsort
!!$    subroutine sort_swap(i,j)
!!$      implicit none
!!$      integer :: i,j
!!$      integer :: tmp
!!$      
!!$      if(i==j) return
!!$      tmp=grid%inv_sort(i)
!!$      grid%inv_sort(i)=grid%inv_sort(j)
!!$      grid%inv_sort(j)=tmp
!!$    end subroutine sort_swap
!!$    integer function compare(iin,jin)
!!$      implicit none
!!$      integer, intent(in) :: iin,jin
!!$      integer :: i,j
!!$      integer :: i1,i3,j1,j3
!!$      real (kind=8) :: i2,j2
!!$      
!!$      i=grid%inv_sort(iin)
!!$      j=grid%inv_sort(jin)
!!$
!!$      i1=grid%rgn(i)
!!$      i2=ipsi(i)
!!$      i3=grid%it(i)
!!$      j1=grid%rgn(j)
!!$      j2=ipsi(j)
!!$      j3=grid%it(j)
!!$      
!!$      if( i1<j1 .or. i1==j1 .and. i2<j2 .or. i1==j1 .and. i2==j2 .and. i3<j3) then
!!$         compare=-1
!!$      else if(i1==j1 .and. i2==j2 .and. i3==j3 ) then
!!$         compare=0
!!$      else
!!$         compare=1
!!$      endif
!!$    end function compare
  end subroutine init_grid
  
  !! initialize adjacent triangle
  subroutine init_triangle(grid)      
    implicit none
    type(grid_type) :: grid
    integer :: i,j,k,l,kp,lp      
    real (kind=8) :: dx1(2),dx2(2),det
    integer :: maxnum,nd
!
    integer :: ndsize, jj, njj, iv, kk, jtrig
    integer, allocatable, dimension(:) :: jlist
    logical, allocatable, dimension(:) :: is_seen
    integer, dimension(4) :: ilist4, jlist4
    logical :: isfound

    ! mapping matrix init
    do i=1, grid%ntriangle

       dx1=grid%x(:,grid%nd(1,i)) - grid%x(:,grid%nd(3,i))
       dx2=grid%x(:,grid%nd(2,i)) - grid%x(:,grid%nd(3,i))

       det=1./( dx1(1)*dx2(2) - dx2(1)*dx1(2))

       grid%mapping(1,1,i)=dx2(2)*det
       grid%mapping(1,2,i)=-dx2(1)*det
       grid%mapping(2,1,i)=-dx1(2)*det
       grid%mapping(2,2,i)=dx1(1)*det


       grid%mapping(:,3,i)=grid%x(:,grid%nd(3,i))
    enddo

    ! set node->triangle array
    allocate(grid%num_t_node(grid%nnode))
    do i=1,2
       grid%num_t_node(:)=0
       do j=1, grid%ntriangle
          do k=1,3
             nd=grid%nd(k,j)
             grid%num_t_node( nd )= grid%num_t_node( nd ) + 1
             if(i==2) grid%tr_node(grid%num_t_node(nd),nd)=j
          enddo
       enddo
       if(i==1) then
          maxnum=maxval(grid%num_t_node(:))
          allocate( grid%tr_node(maxnum,grid%nnode) )
       endif
    enddo
    
!   ----------------------------------------------
!   use temporary storage jlist(:) and is_seen(:)
!   to speed up search for adjacent triangle
!   ----------------------------------------------

    allocate( is_seen(grid%ntriangle), jlist( 3*maxnum + 1 ) )
    is_seen(:) = .false.
    jlist(:) = 0

    grid%adj=0
    ! find adjacent triangle  ->To do :  make it efficient using tr_node later
    do i=1, grid%ntriangle
       ! tr1=>grid%triangle(i)


!      -------------------------------------------
!      store potential list of neighbor triangles
!      in "jlist". Use is_seen(:) to avoid duplicates
!      -------------------------------------------
       is_seen(i) = .true.

       njj = 0
       do k=1,3
         iv = grid%nd(k,i)
         ndsize = grid%num_t_node(iv)

         do kk=1,ndsize
          jtrig = grid%tr_node(kk,iv)
          if (.not.is_seen( jtrig )) then
            njj = njj + 1
            jlist(njj) = jtrig
            is_seen(jtrig) = .true.
          endif
         enddo
       enddo

!      ---------------------------------------
!      reset boolean vector for next iteration
!      ---------------------------------------
       is_seen(i) = .false.
       is_seen( jlist(1:njj) ) = .false.
       
       ilist4(1:3) = grid%nd(1:3,i)
       ilist4(4) = ilist4(1)

       do jj=1,njj
          j = jlist(jj)

          jlist4(1:3) = grid%nd(1:3,j)
          jlist4(4) = jlist4(1)

          if(i/=j) then
             do k=1,3
                kp= k+1
                do l=1,3
                   lp=l+1
                   if ((ilist4(k)+ilist4(kp)).ne.                              &
     &                 (jlist4(l)+jlist4(lp))) then
                      cycle
                   endif
 
                   isfound = ((ilist4(k)  .eq.jlist4(l)).and.                 &
     &                        (ilist4(kp).eq.jlist4(lp))) .or.                &
     &                       ((ilist4(k)  .eq.jlist4(lp)).and.                &
     &                        (ilist4(kp).eq.jlist4(l)))
                   if (isfound) then
                      grid%adj(mod(kp,3)+1,i)=j
                      exit
                   endif
                enddo
             enddo
          end if
       enddo
    enddo

  deallocate( is_seen, jlist )

  end subroutine init_triangle

  !! Coefficient calculation
  subroutine t_coeff(grid,itr,x,p)
    implicit none
    type(grid_type) , intent(in) :: grid
    integer, intent(in) :: itr
    real (kind=8), intent(in) :: x(2)
    real (kind=8), intent(out) :: p(3)

    integer :: nd
    real (kind=8) :: dx(2)

    ! dx=x - grid%x(:,grid%nd(3,itr))
    dx(1:2) = x(1:2) - grid%mapping(1:2,3,itr)
    p(1:2)= grid%mapping(1:2,1,itr)*dx(1) + grid%mapping(1:2,2,itr)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)


  end subroutine t_coeff

  !! search triangle with guess table
  subroutine search_tr(grid,x,itr,p)
    implicit none
    type (grid_type) , intent(in) :: grid  ! grid system
    real (kind=8), intent(in) :: x(2)    ! input position
    integer, intent(out) :: itr !triangle number which contains x
    real (kind=8), intent(out) :: p(3)   !weiting value
    integer :: init,count

    integer :: itr2
    real(kind=8), dimension(3) :: p2

    ! initial guess
    call guess(grid,x,init)

    ! search 
    call search_tr_with_guess(grid,x,init,itr,p,count)
    !eorror message for debug only
    !    if(itr<0) then
    !       print *,"Search failed", x,init,p
    !    endif


  end subroutine search_tr

  !!return initial guess value for triangle search
  subroutine guess(grid,x,init)
    implicit none
    type(grid_type), intent(in) :: grid
    real (kind=8), intent(in) :: x(2)
    integer, intent(out) :: init
    integer :: i(2)

    i= (x-grid%guess_min)*grid%inv_guess_d +1
    !error message for debug only
!    if(i(1)<=0 .or. i(1)>grid%guess_n(1)) then
!       print *, 'Invaild number for guess table- R',i(1),x(1),grid%guess_min(1),grid%guess_max(1)
!    endif
!    if(i(2)<=0 .or. i(2)>grid%guess_n(2)) then
!       print *, 'Invaild number for guess table- Z',i(2),x(2),grid%guess_min(2),grid%guess_max(2)
!    endif
    
    i=min(max(i,1),grid%guess_n)
    init=grid%guess_table(i(1),i(2)) 

  end subroutine guess

  !! initialize guess table
  subroutine init_guess_table(grid)
    implicit none
    type(grid_type) :: grid
    integer :: n(2),i,j,k(2),l,init,itr,count
    real (kind=8) :: x(2),p(3)
    real (kind=8) :: dist, min_dist
    integer ::  min_dist_node,nd,itr2,find
    ! grid%guess%n, min, max  should be initialized

    logical, parameter :: use_original = .false.
    integer, parameter :: itr_init = 1
    logical :: is_found
    integer, dimension(2) :: ij,ijmin,ijmax
    real (kind=8), dimension(2) :: xy,xy1,xy2,xy3,xyc,xymin,xymax
    real (kind=8) :: xmin,xmax,ymin,ymax
    integer :: ilo,ihi,jlo,jhi,i1,j1,i2,j2
    logical, parameter :: do_refine = .false.

    logical :: isok
    integer :: it,ierr
    real(kind=8), allocatable, dimension(:) :: parea
    integer, allocatable, dimension(:) :: iperm

    real(kind=8), parameter :: zero = 0.0d0



    n = grid%guess_n
    grid%guess_d = (grid%guess_max - grid%guess_min)/n
    grid%inv_guess_d = 1D0/grid%guess_d

    allocate( grid%guess_table(n(1),n(2)) )
    grid%guess_table=0 ! added for safety 2006/10/11

    call init_guess_list(grid)

    if (use_original) then

    ! initialize
    init=1

    do i=1, n(1)
       do j=1, n(2)
          ! find most probable triangle index
          k(1)=i
          k(2)= mod(i,2)*j + mod(i+1,2)*(n(2)-j+1)
          x=grid%guess_d*(real(k)-0.5)+grid%guess_min
          
          
          !find nearest node point
          min_dist=1D50
          do nd=1, grid%nnode
             dist=(grid%x(1,nd)-x(1))**2 + (grid%x(2,nd)-x(2))**2
             if(min_dist > dist) then
                min_dist=dist
                min_dist_node=nd
             endif
          enddo

          !search the triangles that has this node point
          find=0
          do l=1, grid%num_t_node(min_dist_node)
             itr=grid%tr_node(l,min_dist_node)
             call t_coeff(grid,itr,x,p)
             ! if( minval(p) >= 0D0 .and. maxval(p) <= 1D0 ) then
             if (minval(p) >= zero) then
                find=1                
                exit
             endif
          enddo
          
          if(find/=1) then
             itr2=itr
             call search_tr_with_guess(grid,x,init,itr2,p,count)
          endif

          if(itr2<0) then
             grid%guess_table(k(1),k(2))=itr
          else
             grid%guess_table(k(1),k(2))=itr2
          endif
             
          
          !search all triangle
          if(itr<0) then
             do l=1, grid%ntriangle
                !get mapping value
                call t_coeff(grid,l,x,p)
                !check inside
                !if( minval(p) >= 0. .and. maxval(p) <=1 ) then
                if (minval(p) >= zero) then
                   itr=l
                endif
             enddo
          endif

          grid%guess_table(k(1),k(2))=itr

!          if(itr>0) init=itr ! set new initial triangle value to old finding one
       enddo
    enddo

  else


!   -----------------------------------  
!   The guess_table only need a close enough
!   starting triangle. This should affect
!   efficiency but not correctness
!   --------------------------------------  

!   -----------------------------------
!   default value for starting triangle
!   -----------------------------------  
    ilo = lbound(grid%guess_table,1)
    ihi = ubound(grid%guess_table,1)
    jlo = lbound(grid%guess_table,2)
    jhi = ubound(grid%guess_table,2)

    itr = itr_init
    grid%guess_table(:,:) = itr




!   --------------------------------
!   simple but sloppy initialization
!   --------------------------------

!   -------------------------------------
!   sort the triangles by the patch size
!   handle the smaller patches last
!   -------------------------------------
    allocate( parea(grid%ntriangle),iperm(grid%ntriangle),stat=ierr )
    call assert(ierr.eq.0,'allocate(parea),ntriangle ',grid%ntriangle)
    
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
     
        parea(itr) = (ymax-ymin)*(xmax-xmin)
     enddo

    call dshelldec(grid%ntriangle,parea,iperm)
!   --------------------------------------
!   check parea(iperm(:)) decreasing order
!   --------------------------------------
    do itr=1,grid%ntriangle-1
       isok = parea(iperm(itr)).ge.parea(iperm(itr+1))
       call assert( isok,'parea not in sorted order',itr)
    enddo


    do it=1,grid%ntriangle
        itr = iperm(it)

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
        grid%guess_table( i1:i2, j1:j2 ) = itr

    enddo

  deallocate( parea, iperm, stat=ierr)
  call assert(ierr.eq.0,'deallocate(parea)',ierr)


!   -----------------------------------------
!   refine triangle assignment  in guess_table
!   -----------------------------------------
  
    do itr=1,grid%ntriangle
        xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
        xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
        xy3(1:2)= grid%x(1:2,grid%nd(3,itr))


        xy(1:2) = xy1(1:2)
        ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
        i1 = max(ilo,min(ihi, ij(1)) )
        j1 = max(jlo,min(jhi, ij(2)) )
        grid%guess_table( i1,j1 ) = itr

        xy(1:2) = xy2(1:2)
        ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
        i1 = max(ilo,min(ihi, ij(1)) )
        j1 = max(jlo,min(jhi, ij(2)) )
        grid%guess_table( i1,j1 ) = itr

        xy(1:2) = xy3(1:2)
        ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
        i1 = max(ilo,min(ihi, ij(1)) )
        j1 = max(jlo,min(jhi, ij(2)) )
        grid%guess_table( i1,j1 ) = itr

!       --------------------
!       centroid of triangle
!       --------------------
        xy(1:2) = (xy1(1:2) + xy2(1:2) + xy3(1:2) )/3.0d0
        ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
        i1 = max(ilo,min(ihi, ij(1)) )
        j1 = max(jlo,min(jhi, ij(2)) )
        grid%guess_table( i1,j1 ) = itr

    enddo


!   ------------------------------------------
!   use triangle determined by init_guess_list
!   ------------------------------------------
    do j=1,ubound(grid%guess_count,2)
    do i=1,ubound(grid%guess_count,1)
       if (grid%guess_count(i,j) .ge. 1) then
          itr = grid%guess_list( grid%guess_xtable(i,j) )
          grid%guess_table(i,j) = itr
       endif
    enddo
    enddo





  endif


!debug-begin
!   print*,'grid%ntriangle ', grid%ntriangle
!   print*,'end of init_guess_table'
!debug-end

  return
  end subroutine init_guess_table



 !! search triangle with initial guess
  subroutine search_tr_with_guess(grid,x,init,itr,p,count)
    implicit none
    type(grid_type),target :: grid
    real (kind=8), intent(in) :: x(2)
    integer, intent(in) :: init
    integer, intent(out) :: itr
    real (kind=8), intent(out) :: p(3)
    integer :: find,min_index(1),next_itr,current_itr
    !debug
    integer :: tmp,i,count,mloc
    real (kind=8) :: mval
    real (kind=8), parameter :: zero = 0.0d0
    real (kind=8), parameter :: eps = 1.0d-10


    find=0
    count=0

    itr=init
    if(init<=0 .or. init>grid%ntriangle) then
       print *, 'invalid guess. init=1 is used instead',init
       itr=1
    endif

    
    do while(find/=1 .and. count<grid%nnode)
       count=count+1
       call t_coeff(grid,itr,x,p)
       ! if( minval(p) >= 0D0 .and. maxval(p) <=1D0 ) then
       if (minval(p) >= zero) then
          find=1
       else
          min_index=minloc(p)
          next_itr=grid%adj(min_index(1),itr)
          if(next_itr>0) then
             itr=next_itr
          else
             itr=-1
             return
          endif
       endif
    enddo
    if(count>=grid%nnode) then 
       call t_coeff(grid,itr,x,p)
       ! if( minval(p) > -1D-10 .and. maxval(p) <1+1D-10) then
       if (minval(p) > -eps) then
          find=1
          mval=-1D0
          do i=1, 3
             if(mval < p(i)) then
                mval=p(i)
                mloc=i
             endif
          enddo
          p(:)=0D0
          p(mloc)=1D0          
       else
          itr=-1
          print *, 'search error : p=',p
!       print *, 'Too large loop number in search tr'
       endif
    endif

  end subroutine search_tr_with_guess
  

end module grid_class

include 'init_guess_list.F90'
include 'check_guess_table.F90'
include 'dshelldec.F90'
include 'checkoverlap.F90'
include 'search_vtr.F90'


!! Return MPI processor information
subroutine get_pe_info(mype,totalpe)
  use sml_module
  implicit none
  integer :: mype, totalpe

  mype=sml_mype
  totalpe=sml_totalpe

end subroutine get_pe_info

!! return number toroidal cross section 
subroutine get_nphi(nphi_each)
  use sml_module
  implicit none
  integer :: nphi_each
  
  nphi_each=sml_plane_per_pe
end subroutine get_nphi
  
subroutine get_node_vol(grid)
#ifdef CIRCULAR_SPECIAL
  use eq_module
#endif
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
!  real (kind=8) , allocatable :: tr_area(:)  !,tr_vol(:)
  real (kind=8) :: vol_tmp,dx1(2),dx2(2),com_r,area,area_tmp
  integer :: i,j,tr,nd(3)
#ifdef CIRCULAR_SPECIAL
  integer :: n, n1, n2
  real (kind=8) :: dr,r1, r2,factor,sumv
#endif

  !  allocate(tr_vol(grid%ntriangle),tr_area(grid%ntriangle)
!  allocate(tr_area(grid%ntriangle))
  

  do i=1, grid%ntriangle
     ! find triangle area - 1/2 * vec1 x vec2
     nd(:)=grid%nd(:,i)
     dx1(:)=grid%x(:,nd(1))- grid%x(:,nd(3))
     dx2(:)=grid%x(:,nd(2))- grid%x(:,nd(3))
     area= 0.5D0 * abs( dx1(1)*dx2(2) - dx1(2)*dx2(1))
     ! find triangle center of mass - Radius
     com_r = 1D0/3D0 * ( grid%x(1,nd(1)) + grid%x(1,nd(2)) +grid% x(1,nd(3)) )
     !get volume - area * 2pi * Radius / nphi
     grid%tr_vol(i)= area * sml_2pi * com_r / real(sml_nphi_total)
     grid%tr_area(i)=area
  enddo
  do i=1, grid%nnode
!     vol_tmp =0D0
     area_tmp=0D0
     do j=1, grid%num_t_node(i)
!        vol_tmp =vol_tmp +tr_vol( grid%tr_node(j,i))/3D0
        area_tmp=area_tmp+grid%tr_area(grid%tr_node(j,i))/3D0
     enddo
     grid%node_area(i)=area_tmp
     grid%node_vol(i)=area_tmp*sml_2pi*grid%x(1,i)/real(sml_nphi_total)
     grid%inv_node_vol(i)=1D0/grid%node_vol(i)
  enddo

  allocate(grid%surf_vol(grid%npsi))
  do i=1, grid%npsi
     vol_tmp=0D0
     do j=1, grid%ntheta(i)
        vol_tmp=vol_tmp + 1D0/grid%inv_node_vol( grid%itheta0(i) + j - 1 )
     enddo
     grid%surf_vol(i)=vol_tmp*real(sml_nphi_total)
  enddo

  do i=1, grid%ntriangle
     
  enddo

#ifdef CIRCULAR_SPECIAL
  do i=1, grid%npsi
     if(i==1) then
        n1=grid%itheta0(i)
        n2=grid%itheta0(i+1)
        factor=1.
     else if(i==grid%npsi) then
        n1=grid%itheta0(i-1)
        n2=grid%itheta0(i)
        factor=1.
     else       
        n1=grid%itheta0(i-1)
        n2=grid%itheta0(i+1)
        factor=0.5
     endif
     r1=sqrt((grid%x(1,n1)-eq_axis_r)**2 + (grid%x(2,n1)-eq_axis_z)**2)
     r2=sqrt((grid%x(1,n2)-eq_axis_r)**2 + (grid%x(2,n2)-eq_axis_z)**2)
     dr=factor*(r2-r1)     
     grid%surf_vol(i)=(r2**2-r1**2)*sml_pi*factor*sml_2pi*eq_axis_r/real(sml_nphi_total)
      
     sumv=0D0
     do j=1,grid%ntheta(i)
        n=grid%itheta0(i)+j-1
        grid%node_vol(n)=grid%x(1,n)**2 
        sumv=sumv+grid%x(1,n)**2
     enddo     
     do j=1, grid%ntheta(i)
        n=grid%itheta0(i)+j-1
        grid%node_vol(n)=grid%surf_vol(i)/sumv*grid%node_vol(n)
        grid%inv_node_vol(n)=1D0/grid%node_vol(n)
     enddo
  enddo
  grid%surf_vol=real(sml_nphi_total)*grid%surf_vol
#endif
  !debug
  if(sml_mype==0) print *, 'sum(tr_vol)',sum(grid%tr_vol)
  if(sml_mype==0) print *, 'sum(node_area)',sum(grid%node_area)
  if(sml_mype==0) print *, 'sum(node_vol)',sum(grid%node_vol)
  if(sml_mype==0) print *, 'sum(surf_vol)',sum(grid%surf_vol)

!  deallocate(tr_area) !,tr_vol)
end subroutine get_node_vol


subroutine search_tr2( grid, xy, itr, p )
  use grid_class
  implicit none
  type(grid_type) :: grid
  real(kind=8) :: xy(2)
  integer :: itr
  real(kind=8) :: p(3)

  real(kind=8), parameter :: zero = 0.0d0
  real(kind=8), parameter :: eps = 10.0d0*epsilon(zero)
  integer :: ij(2), istart,iend, k, itrig
  integer :: i,j,  ilo,ihi,  jlo,jhi
  real(kind=8) ::  dx(2), pmin, pmax, dp
  logical :: is_found


  ilo = lbound( grid%guess_table, 1 )
  ihi = ubound( grid%guess_table, 1 )

  jlo = lbound( grid%guess_table, 2 )
  jhi = ubound( grid%guess_table, 2 )

  ij = (xy - grid%guess_min)*grid%inv_guess_d + 1
  i = max(ilo, min(ihi, ij(1)) )
  j = max(jlo, min(jhi, ij(2)) )


  istart = grid%guess_xtable(i,j)
  iend = istart + grid%guess_count(i,j) - 1


  itr = -1
  do k=istart,iend
     itrig = grid%guess_list(k)
     ! call t_coeff( grid, itrig, xy, p )

    dx(1:2) = xy(1:2) - grid%mapping(1:2,3,itrig)
    p(1:2)= grid%mapping(1:2,1,itrig)*dx(1) +                          &
            grid%mapping(1:2,2,itrig)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)

     if (minval(p) .ge. -eps) then
        itr = itrig
        exit
     endif
  enddo

  return
end subroutine search_tr2





