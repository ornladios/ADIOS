	module adi_module
        use prof_mod
	use mpi_module
	implicit none

!
! unsymmetric block tridiagonal band matrix
!
! [  A1  C1               ]
! [  B2  A2   C2          ]
! [      B3   A3   C3 ... ]
!
!
! diagonal blocks A's are msize x msize fully dense
!
! subdiagonal blocks B's are msize x msize diagonal
!
! superdiagonal blocks C's are msize x msize diagonal
!


!       -----------------------------
!       default everything is private 
!       -----------------------------
        private

!	------------------------------------------------------
!	determine whether to use MPI_Gatherv and MPI_Scatterv
!	or use explicit MPI_Isend and MPI_Irecv
!	------------------------------------------------------
        logical, parameter :: use_gatherv = .false.
        logical, parameter :: use_scatterv = .false.

	integer, parameter :: idebug = 1

!	--------------------------------------------------------
!	this may change depending on whether there is automatic
!	promotion of real to real*8
!	--------------------------------------------------------
	integer, parameter :: mpi_datatype = MPI_DOUBLE_PRECISION

	integer, parameter :: root = 0
	real, parameter :: zero = 0.0
	real, parameter :: one = 1.0

	logical :: is_mapping_initialized = .false.

	integer :: mypid = 0
	integer :: nproc = 1

	real, allocatable, dimension(:,:,:) :: ALU
	real, allocatable, dimension(:,:) :: d,dl,du,du2 
	integer, allocatable, dimension(:,:) :: ALU_ipiv, ipiv

	integer, allocatable,dimension(:) ::  ibstart,ibend,ibsize
	integer :: gnblocks
!	----------------------
!	interface block to MPI
!	----------------------
	interface

	subroutine MPI_ISEND( buf, count, datatype, dest, 
     &		tag, comm, request, ierror)
	integer count,datatype,dest,tag,comm,request,ierror
	real*8 buf(*)
	end subroutine MPI_ISEND

	subroutine MPI_IRECV( buf, count, datatype, source, 
     &		tag, comm, request, ierror )
	integer count,datatype,source,tag,comm,request,ierror
	real*8 buf(*)
	end subroutine MPI_IRECV


	subroutine MPI_COMM_RANK( comm, rank, ierror )
	integer comm,rank,ierror
	end subroutine MPI_COMM_RANK

	subroutine MPI_COMM_SIZE( comm, size, ierror )
	integer comm, size, ierror
	end subroutine MPI_COMM_SIZE

	end interface


!	-------------------------
!	interface block to LAPACK
!	-------------------------
	interface

!	------------------------------
!	LU factorization with pivoting
!	------------------------------
	subroutine dgetrs( trans, n, nrhs,
     &			A, lda,ipiv,B,ldb, info )
	character trans
	integer n,nrhs,lda,ldb,info
	real*8 A(lda,*),B(ldb,nrhs)
	integer ipiv(*)
	end subroutine dgetrs

	subroutine dgetrf( m,n, A, lda, ipiv, info )
	integer m,n,lda,info
	integer ipiv(*)
	real*8 A(lda,*)
	end subroutine dgetrf

!	-------------------
!	tridiagonal solvers
!	-------------------

	subroutine dgttrf( m, dl, d, du, du2, ipiv, info )
	integer m,info,ipiv(*)
	real*8 dl(*),d(*),du(*),du2(*)
	end subroutine dgttrf

	subroutine dgttrs(trans,m,nrhs,dl,d,du,du2,ipiv,B,ldb,info)
	character trans
	integer m,nrhs,ldB,info
	real*8 dl(*),d(*),du(*),du2(*),B(ldb,nrhs)
	integer ipiv(*)
	end subroutine dgttrs


	end interface


!       --------------------------
!       export only the following
!       --------------------------
        public :: adi_matvec, adi_dotprod,adi_vnorm, adi_dotprod_safe
        public :: adi_bicgstab, adi_gmres, adi_fixpoint, msgtag

	public :: adi_tridiag, adi_init
	public :: adi_gather,adi_scatter, adi_cleanup


	contains

	subroutine assert( lcond, message, ivalue )
	implicit none
	logical, intent(in) :: lcond
	character(len=*) :: message
	integer, intent(in) :: ivalue

	if (.not.lcond) then
	   write(*,*) message, ivalue
	   stop '** assertioin error ** '
	endif
	return
	end subroutine assert

	integer function msgtag(src,dest)
	implicit none
	integer, intent(in) :: src,dest
	integer, parameter :: msgbase = 999

	msgtag = msgbase + src + dest*nproc
	return
	end function msgtag


	subroutine adi_mapping( nblocks, context )
	implicit none
	integer, intent(in) :: nblocks, context
!	---------------------------------------------------
!	adi_mapping setup internal data structure
!	on the global decomposition on a 1-D processor grid
!	---------------------------------------------------


!	---------------
!	local variables
!	---------------
	integer, allocatable, dimension(:) :: ibuf
	integer :: ierror,iproc, count,datatype,op


        call MPI_COMM_SIZE( context, nproc, ierror )
        if (ierror.ne.MPI_SUCCESS) then
          write(*,*) 'adi_mapping: mpi_comm_size returns ierror',ierror
          stop '** error in adi_mapping '
        endif

        call MPI_COMM_RANK(context, mypid, ierror )
        if (ierror.ne.MPI_SUCCESS) then
          write(*,*) 
     &        'adi_mapping: mpi_comm_rank returns ierror ', ierror
          stop '** error in adi_mapping '
        endif

	allocate(
     &		ibstart(0:(nproc-1)),
     &		ibend(0:(nproc-1)),
     &		ibsize(0:(nproc-1)) )
	

!	-------------------
!	setup decomposition
!	-------------------
	allocate( ibuf(0:(nproc-1)) )
	do iproc=0,(nproc-1)
	  ibuf(iproc) = 0
	enddo
	ibuf(mypid) = nblocks

	count = nproc
	datatype = MPI_INTEGER
	op = MPI_SUM
	call MPI_ALLREDUCE( ibuf, ibsize, count, datatype, op, 
     &			context, ierror )
	call assert( ierror .eq. MPI_SUCCESS,
     &		'adi_init: mpi_allreduce returns ierror ', ierror )


	ibstart(0) = 1
	do iproc=0,(nproc-2)
	   ibend(iproc) = ibstart(iproc) + ibsize(iproc)-1
	   ibstart(iproc+1) = ibend(iproc) + 1
	enddo
	iproc = nproc-1
	ibend( iproc ) = ibstart(iproc) + ibsize(iproc)-1

	call assert( ibsize(mypid) .eq. nblocks,
     &		'adi_mapping: error in ibsize(mypid) ', ibsize(mypid))
	call assert( ibend(mypid) - ibstart(mypid) + 1 .eq. nblocks,
     &		'adi_mapping: error in ibend or ibstart, nblocks = ',nblocks)

	gnblocks = 0
	do iproc=0,(nproc-1)
	  gnblocks = gnblocks + ibsize(iproc)
	enddo

	call assert( gnblocks .eq. ibend( (nproc-1) ),
     &		'adi_mapping: invalid ibend( (nproc-1) ) ', 
     &		ibend( (nproc-1)) )

	
	deallocate( ibuf )

	return
	end subroutine adi_mapping


        subroutine adi_scatter( msize, nblocks,  gnblocks, 
     &		Alocal, Aglobal, context )
	implicit none
	integer, intent(in) :: msize,nblocks, gnblocks,context
	real, intent(inout), dimension(msize,nblocks) :: Alocal
	real, intent(in), dimension(gnblocks,msize) :: Aglobal


!	---------------------------------------
!	Note: gbuf and buf start indexing from 0
!	to be compatible with C array indexing
!	---------------------------------------
	real, dimension(0:msize*gnblocks) :: gbuf
	real, dimension(0:msize*nblocks) :: buf
	integer :: iproc,ipos,i,iblock
	integer ::  sendtype, recvtype, datatype, ierror,recvcount
	integer, dimension(0:(nproc-1)) :: sendcounts, displs

        integer, dimension(MPI_STATUS_SIZE,0:nproc) :: status
        integer, dimension(0:nproc) :: request
        integer :: isrc,idest,tag


c        call profstart('adi_scatter')

	if (mypid.eq.root) then
	    ipos = 0
	    do iblock=1,gnblocks
	    do i=1,msize
	       gbuf(ipos) = Aglobal(iblock,i)
	       ipos = ipos + 1
	    enddo
	    enddo

	    if (idebug.ge.2) then
		print*,'in scatter, sum(Aglobal) ', sum( Aglobal)
	    endif

	    do iproc=0,(nproc-1)
	       sendcounts(iproc) = ibsize(iproc) * msize
	    enddo

	    do iproc=0,(nproc-1)
		ipos = (ibstart(iproc)-1)*msize 
		displs(iproc) = ipos
	    enddo
	endif


	sendtype = mpi_datatype
	recvtype = mpi_datatype
	recvcount = nblocks * msize
	
        if (use_scatterv) then
!	----------------------------------------------
!	use MPI_Scatterv
!	may have unexplained slow down on some machine
!	with over 17 MPI tasks
!	----------------------------------------------
c	   call profstart('adi_scatter:mpi_scatterv')

	   call MPI_SCATTERV( gbuf, sendcounts, displs, sendtype,
     &		buf, recvcount, recvtype, 
     &		root, context, ierror)

c	   call profend('adi_scatter:mpi_scatterv')
           call assert( ierror.eq.MPI_SUCCESS,                                &
     &       'adi_scatter: mpi_scatterv returns ierror ', ierror )
        else
!	---------------
!	use isend/irecv
!	---------------
c          call profstart('adi_scatter:msg')
          isrc = root
          idest = mypid
          tag = msgtag( isrc,idest )

          call MPI_Irecv( buf, recvcount, recvtype,                           &
     &            isrc, tag, context, request(nproc), ierror )
          call assert( ierror.eq.MPI_SUCCESS,                                 &
     &            'adi_scatter: mpi_irecv returns ierror ',ierror)

          if (mypid.ne.root) then
             call MPI_Wait( request(nproc), status, ierror )
             call assert( ierror.eq.MPI_SUCCESS,                              &
     &           'adi_scatter:mpi_wait returns ierror ', ierror )

          else
             do iproc=0,(nproc-1)
               isrc = root
               idest = iproc
               tag = msgtag( isrc,idest )
               ipos = displs(iproc)

               call MPI_Isend( gbuf(ipos), sendcounts(iproc),                 &
     &           sendtype, idest, tag, context, request(iproc), ierror )
               call assert( ierror.eq.MPI_SUCCESS,                            &
     &            'adi_scatter:mpi_isend returns ierror ', ierror )
     
	     enddo

             call MPI_Waitall( nproc+1, request, status, ierror )
             call assert( ierror.eq.MPI_SUCCESS,                              &
     &             'adi_scatter:mpi_waitall returns ierror ', ierror )
          endif


c          call profend('adi_scatter:msg')

        endif



!	--------
!	copy out
!	--------
	ipos = 0
	do iblock=1,nblocks
	do i=1,msize
	   Alocal(i,iblock) = buf(ipos)
	   ipos = ipos + 1
	enddo
	enddo

c        call profend('adi_scatter')

        if (idebug.ge.2) then
           if (mypid.eq.root) then
             print*,'in scatter: sum(Alocal) ', sum(Alocal)
           endif
        endif
              

	return
	end subroutine adi_scatter


        subroutine adi_gather( msize, nblocks,  gnblocks, 
     &		Alocal, Aglobal, context )
	implicit none
	integer, intent(in) :: msize,nblocks, gnblocks,context
	real, intent(in), dimension(msize,nblocks) :: Alocal
	real, intent(inout), dimension(gnblocks,msize) :: Aglobal


!	---------------------------------------
!	Note: gbuf and buf start indexing from 0
!	to be compatible with C array indexing
!	---------------------------------------
	real, dimension(0:msize*gnblocks) :: gbuf
	real, dimension(0:msize*nblocks) :: buf
	integer :: iproc,ipos,i,iblock
	integer ::  sendtype, recvtype, datatype, ierror,sendcount
	integer, dimension(0:(nproc-1)) :: recvcounts, displs

	integer,dimension(0:nproc) :: request
	integer,dimension(MPI_STATUS_SIZE,0:nproc) :: status
        integer :: tag,isrc,idest


        if (idebug.ge.2) then
           if (mypid.eq.root) then
              print*,'in gather: (sum(Alocal)) ', (sum(Alocal))
           endif
        endif


c        call profstart('adi_gather')

	ipos = 0
	do iblock=1,nblocks
	do i=1,msize
	   buf(ipos) = Alocal(i,iblock)
	   ipos = ipos + 1
	enddo
	enddo

	sendcount = msize*nblocks
	sendtype = mpi_datatype

	do iproc=0,(nproc-1)
	   recvcounts(iproc) = ibsize(iproc) * msize
	enddo

	do iproc=0,(nproc-1)
	    ipos = (ibstart(iproc)-1)*msize 
	    displs(iproc) = ipos
	enddo

	recvtype = mpi_datatype

	
        if (use_gatherv) then
!	---------------
!	use MPI_Gatherv
!	may have unexplained slow down on some machine
!	with over 17 MPI tasks
!	---------------
c	   call profstart('adi_gather:mpi_gatherv')
	   call MPI_GATHERV( buf, sendcount, sendtype, 
     &		gbuf, recvcounts, displs, recvtype, 
     &		root, context, ierror)

           call assert( ierror.eq.MPI_SUCCESS,                                &
     &       'adi_gather: mpi_gatherv returns ierror ',ierror)

c	   call profend('adi_gather:mpi_gatherv')
        else
!          ---------------
!          use isend/irecv
!          ---------------
c	  call profstart('adi_gather:msg')

           if (mypid.ne.root) then

             isrc = mypid
	     idest = root
             tag = msgtag( isrc, idest )

             call MPI_Isend( buf, sendcount, sendtype, idest,                 &
     &          tag, context, request(nproc), ierror )
             call assert( ierror.eq.MPI_SUCCESS,                              &
     &       'adi_gather:mpi_isend returns ierror ', ierror )

             call MPI_Wait( request(nproc), status, ierror )
             call assert( ierror.eq.MPI_SUCCESS,                              &
     &          'adi_gather:mpi_iwait returns ierror ', ierror )
           else

!	-----------------
!	post irecvs first
!	-----------------
             do iproc=0,(nproc-1)
               ipos = displs(iproc)
	       isrc = iproc
               idest = root
               tag = msgtag( isrc, idest )
               call MPI_Irecv( gbuf(ipos),                                    &
     &              recvcounts(iproc), recvtype,                              &
     &              isrc, tag, context, request(iproc), ierror )
               call assert( ierror.eq.MPI_SUCCESS,                            &
     &            'adi_gather:mpi_irecv returns ierror ', ierror)
	     enddo

             isrc = mypid
	     idest = root
             tag = msgtag( isrc, idest )

             call MPI_Isend( buf, sendcount, sendtype, idest,                 &
     &          tag, context, request(nproc), ierror )
             call assert( ierror.eq.MPI_SUCCESS,                              &
     &       'adi_gather:mpi_isend returns ierror ', ierror )
	     
	     call MPI_Waitall( nproc+1, request, status, ierror )
             call assert( ierror.eq.MPI_SUCCESS,                              &
     &          'adi_gather:mpi_waitall returns ierror ', ierror)

            endif
c	  call profend('adi_gather:msg')

        endif



	if (mypid.eq.root) then

	    ipos = 0
	    do iblock=1,gnblocks
	    do i=1,msize
		Aglobal(iblock,i) = gbuf( ipos )
		ipos = ipos + 1
	    enddo
	    enddo

	    if (idebug.ge.2) then
		print*,'in gather, sum(Aglobal) = ', sum( Aglobal)
	    endif

	endif

c        call profend('adi_gather')

	return
	end subroutine adi_gather




	subroutine adi_init( msize, nblocks, A, B, C, context )
	implicit none

	integer, intent(in) :: msize, nblocks
	integer, intent(in) :: context
	real, intent(in),dimension(msize,msize,nblocks) :: A
	real, intent(in),dimension(msize,nblocks) :: B, C

!	---------------
!	local variables
!	---------------

	integer :: i,j,iblock,ierror,lda
	integer, dimension(nblocks) :: linfo_lu
	integer, dimension(msize) :: linfo
	integer :: mm,nn
	real, dimension(msize,nblocks) :: Adiag

	if (.not. is_mapping_initialized) then
	    is_mapping_initialized = .true.
	    call adi_mapping( nblocks, context )
	endif



!$omp   parallel do private(i,iblock)
        do iblock=1,nblocks
	do i=1,msize
	   Adiag(i,iblock) = A(i,i,iblock)
	enddo
	enddo

!	--------------------------
!	collect off-diagonal parts
!	--------------------------



        allocate(
     &		dl(gnblocks,msize),
     &		d(gnblocks,msize),
     &		du(gnblocks,msize),
     &		du2(gnblocks,msize),
     &          ipiv(gnblocks,msize))
   

!$omp   parallel do private(i,iblock)
        do i=1,msize
        do iblock=1,gnblocks

	dl(iblock,i) = zero
	d(iblock,i) = zero
	du(iblock,i) = zero
	du2(iblock,i) = zero
	ipiv(iblock,i) = 0

        enddo
        enddo

        call adi_gather( msize, nblocks,  gnblocks,
     &          Adiag, d, context )
        call adi_gather( msize, nblocks,  gnblocks,
     &          B, dl, context )
        call adi_gather( msize, nblocks,  gnblocks,
     &          C, du, context )



! ---------------------------------------
! factor off-diagonal tridiagonal systems
! ---------------------------------------

	if (mypid.eq.root) then


!$omp      parallel do private(i)
	   do i=1,msize
	       call dgttrf( gnblocks, dl(2,i), 
     &			d(1,i), du(1,i), du2(1,i),
     &			ipiv(1,i), linfo(i) )
	   enddo

	   do i=1,msize
	   if (linfo(i).ne.0) then
	      write(*,*) 
     &		'adi_init:dgttrf returns info= ',linfo(i)
	      stop '** error in adi_init '
	   endif
	   enddo

	endif

! ----------------------------------------
! perform LU factorization of dense blocks
! ----------------------------------------
	allocate( 
     &		ALU(msize,msize,nblocks),
     &		ALU_ipiv(msize,nblocks) )


!$omp   parallel do private(iblock,i,j)
        do iblock=1,nblocks
          do j=1,msize
          do i=1,msize
	   ALU(i,j,iblock) =  A(i,j,iblock)
          enddo
          enddo

          do i=1,msize
           ALU_ipiv(i,iblock) = 0
          enddo
        enddo

	lda = msize
	mm = msize
	nn = msize

!$omp   parallel do private(iblock)
	do iblock=1,nblocks
	   call dgetrf( mm,nn, ALU(1,1,iblock), lda,
     &		ALU_ipiv(1,iblock), linfo_lu(iblock) )
	enddo

	do iblock=1,nblocks
	   if (linfo_lu(iblock).ne.0) then
	      write(*,*) 'adi_init: dgetrf returns info = ', 
     &            linfo_lu(iblock)
	      stop '** error in adi_init '
	   endif
	enddo




	return
	end subroutine adi_init

	subroutine adi_cleanup()
	implicit none
	logical, parameter :: deallocate_mapping = .false.

	if (deallocate_mapping) then

		deallocate(ibstart)
		deallocate(ibend)
		deallocate(ibsize)


	endif

 
	deallocate(ALU)
	deallocate(ALU_ipiv)




	deallocate(dl)
	deallocate(du)
	deallocate(du2)
	deallocate(d)
	deallocate(ipiv)




	return
	end subroutine adi_cleanup


	  

	subroutine adi_dblock( msize,nblocks, A,B,C, rhs,  context )
	implicit none

	integer, intent(in) :: msize,nblocks, context
	real, intent(in),dimension(msize,msize,nblocks) :: A
	real, intent(in),dimension(msize,nblocks) :: B,C
	real, intent(inout),dimension(msize,nblocks)  ::  rhs

!	---------------
!	local variables
!	---------------
	real, dimension(msize,nblocks)  ::  x
	integer :: i,iblock,nrhs,ldalu,ldb,mm
	integer, dimension(nblocks) :: linfo_lu
	character :: trans


!	-------------------------
!	no communication required
!	-------------------------

!	-------
!	copy in
!	-------
	x(:,:) = rhs(:,:)


	trans = 'N'
	nrhs = 1
        ldalu = msize
        ldb = msize
	mm = msize



!$omp   parallel do private(iblock)
	do iblock=1,nblocks

	call dgetrs( trans, mm, nrhs, 
     &        ALU(1,1,iblock), ldalu, ALU_ipiv(1,iblock),
     &        x(1,iblock), ldb, linfo_lu(iblock) )

	enddo


	do iblock=1,nblocks
	  if (linfo_lu(iblock).ne.0) then
	    write(*,*) 'adi_dblock: dgetrs returns ',linfo_lu(iblock)
	    stop '** error in adi_dblock '
	  endif
	enddo

!	--------
!	copy out
!	--------
	rhs(:,:) = x(:,:)


	return
	end subroutine adi_dblock


	subroutine adi_tridiag( msize,nblocks, A,B,C, rhs,  context )
	implicit none

        integer, intent(in) :: msize,nblocks, context 
        real, intent(in),dimension(msize,msize,nblocks) :: A
        real, intent(in),dimension(msize,nblocks) :: B,C
        real, intent(inout),dimension(msize,nblocks)  ::  rhs


!	----------------
!	local variables
!	----------------
	integer :: i,iblock,  ierror,nrhs,ldb
	real, dimension(gnblocks,msize) :: grhs
	integer, dimension(msize) :: linfo




c        call profstart('gath_scatt')
!	-------------------------
!	collect to root processor
!	-------------------------
	call adi_gather(msize,nblocks,gnblocks, rhs, grhs, context )


!	--------------------------
!	perform tridiagonal solves
!	--------------------------

        if (mypid.eq.root) then

	nrhs = 1
	ldb = gnblocks
	do i=1,msize
           call dgttrs( 'Notrans', gnblocks, nrhs,
     &                  dl(2,i),d(1,i),du(1,i),
     &                  du2(1,i),ipiv(1,i),grhs(1,i),ldb,linfo(i) )
	enddo

	do i=1,msize
	   if (linfo(i).ne.0) then
		write(*,*) 'adi_tridiag: dgttrs returns info ', linfo(i)
		stop '** error in adi_tridiag '
	   endif
	enddo

        endif
		


!	----------------
!	send results out
!	----------------

	call adi_scatter(msize,nblocks,gnblocks, rhs, grhs, context )

c        call profend('gath_scatt')
c        call profstat(mypid)
	end subroutine adi_tridiag



	subroutine adi_matvec( msize,nblocks, A, B, C, x, Ax, 
     &		context, cal_dblock)
	implicit none

!	--------------------------------------
!	if cal_dblock is true, the matvec
!	involves the dense blocks
!	if cal_dblock is false, the matvec
!	involves only the off diagonal blocks
!	--------------------------------------

	integer, intent(in) :: msize,nblocks,context
	logical, intent(in) :: cal_dblock
	real, intent(in),dimension(msize,msize,nblocks) :: A
	real, intent(in),dimension(msize,nblocks) :: B, C
	real, intent(in),dimension(msize,nblocks) :: x
	real, intent(inout),dimension(msize,nblocks) :: Ax


!	---------------
!	local variables
!	---------------
	integer :: i,iblock
	integer :: src,dest,datatype,tag,count
	integer :: ierror, incx,incy,lda
	integer, dimension(MPI_STATUS_SIZE) :: status
	integer :: rleft_request, sleft_request
	integer :: rright_request, sright_request
	integer :: leftpid,rightpid
	logical :: need_recv_left,need_recv_right
	integer :: mm,nn

	real, dimension(msize) :: xleft,xright
	real :: alpha, beta

        logical, parameter :: use_gemv = .true.

        if (.not. is_mapping_initialized) then
            is_mapping_initialized = .true.
            call adi_mapping( nblocks, context )
        endif



!	-----------------
!	send out messages
!	-----------------
	need_recv_left =  (ibstart(mypid) .gt. 1)
	need_recv_right = (ibend(mypid) .lt. gnblocks)

	if (need_recv_left) then
	   leftpid  = max(0,min( (nproc-1), mypid-1))

	   count = msize
	   datatype = mpi_datatype
	   src  = leftpid
	   dest = mypid
	   tag = msgtag( src, dest )
	   call MPI_IRECV( xleft, count, datatype, src, tag, 
     &			context, rleft_request, ierror )
	   if (ierror.ne.MPI_SUCCESS) then
		write(*,*) 'adi_matvec: mpi_irecv returns ', ierror
		stop '** error in adi_matvec '
	   endif


	   iblock = 1
	   count = msize
	   datatype = mpi_datatype
	   src = mypid
	   dest = leftpid
	   tag = msgtag( src,dest )
	   call MPI_ISEND( x(1,iblock), count, datatype,
     &		dest, tag, context, sleft_request, ierror )
	   if (ierror.ne.MPI_SUCCESS) then
		write(*,*) 'adi_matvec: mpi_isend returns ', ierror
		stop '** error in adi_matvec '
	   endif
	endif


	if (need_recv_right) then
	   rightpid = max(0,min( (nproc-1), mypid+1))

	   count = msize
	   datatype = mpi_datatype
	   src = rightpid
	   dest = mypid
	   tag = msgtag( src, dest )
	   call MPI_IRECV( xright, count, datatype, src, tag, 
     &			context, rright_request, ierror )
	   if (ierror.ne.MPI_SUCCESS) then
		write(*,*) 'adi_matvec: mpi_irecv returns ', ierror
		stop '** error in adi_matvec '
	   endif


	   iblock = nblocks
	   count = msize
	   datatype = mpi_datatype
	   src = mypid
	   dest = rightpid
	   tag = msgtag( src,dest )
	   call MPI_ISEND( x(1,iblock), count, datatype,
     &		dest, tag, context, sright_request, ierror )
	   if (ierror.ne.MPI_SUCCESS) then
		write(*,*) 'adi_matvec: mpi_isend returns ', ierror
		stop '** error in adi_matvec '
	   endif
	endif



!	----------------------------------------
! 	perform computation in the interior part
!	----------------------------------------
	if (cal_dblock) then

!	-----------------------------------
!	multiply with dense diagonal blocks
!	-----------------------------------
	alpha = one
	beta = zero
	incx = 1
	incy = 1
	lda = msize
	mm = msize
	nn = msize

!$omp	parallel do private(iblock)
	do iblock=1,nblocks
	  if (use_gemv) then
	   call dgemv( 'N', mm,nn,
     &		alpha, A(1,1,iblock),lda,
     &			x(1,iblock),incx,
     &		beta,   Ax(1,iblock),incy )
          else
	    Ax(1:msize,iblock) = matmul( 
     &         A(1:msize,1:msize,iblock), x(1:msize,iblock))
          endif
	
	enddo

	else

!$omp   parallel do private(i,iblock)
	do iblock=1,nblocks
	do i=1,msize
	  Ax(i,iblock) = zero
	enddo
	enddo

	endif

!	-------------------
!	lower diagonal part
!	-------------------

!$omp   parallel do private(iblock,i)
	do iblock=2,nblocks
	   do i=1,msize
		Ax(i,iblock) = Ax(i,iblock) + 
     &			B(i,iblock) * x(i,iblock-1)
	   enddo
	enddo

!	-------------------
!	upper diagonal part
!	-------------------
!$omp   parallel do private(iblock,i)
	do iblock=1,(nblocks-1)
	   do i=1,msize
		Ax(i,iblock) =  Ax(i,iblock) + 
     &			C(i,iblock) * x(i,iblock+1)
	   enddo
	enddo
	

	if (need_recv_left) then
	   call MPI_WAIT( sleft_request, status, ierror )
	   call assert( ierror .eq. MPI_SUCCESS,
     &		'adi_matvec: mpi_wait of sleft_request returns ', ierror)

	   call MPI_WAIT( rleft_request, status, ierror )
	   call assert( ierror .eq. MPI_SUCCESS,
     &		'adi_matvec: mpi_wait of rleft_request returns ', ierror)

	   iblock = 1
	   do i=1,msize
	 	Ax(i,iblock) = Ax(i,iblock) +
     &			B(i,iblock) * xleft(i)
	   enddo
        endif

	if (need_recv_right) then
	   call MPI_WAIT( sright_request, status, ierror )
	   call assert( ierror .eq. MPI_SUCCESS,
     &		'adi_matvec: mpi_wait of sright_request returns ', ierror)

	   call MPI_WAIT( rright_request, status, ierror )
	   call assert( ierror .eq. MPI_SUCCESS,
     &		'adi_matvec: mpi_wait of rright_request returns ', ierror)


	   iblock = nblocks
	   do i=1,msize
		Ax(i,iblock) = Ax(i,iblock) +
     &			C(i,iblock) * xright(i)
	   enddo
	 endif

	return
	end subroutine adi_matvec
		

	   

        subroutine adi_precond( msize,nblocks, A,B,C, rhs,  context )
        implicit none
 
        integer, intent(in) :: msize,nblocks, context
        real, intent(in),dimension(msize,msize,nblocks) :: A
        real, intent(in),dimension(msize,nblocks) :: B,C
        real, intent(inout),dimension(msize,nblocks)  ::  rhs               


!	---------------
!	local variables
!	---------------
	logical, parameter :: cal_dblock = .true.
	integer :: i,iblock
	real, dimension(msize,nblocks) :: rvec, x1
	logical :: need_print
!	-------------------------------------
!	istep = 0 for no preconditioning
!	istep = 1 for preconditioning with just diagonal block
!	istep = 2 for ADI like precondition
!	-------------------------------------
	integer, parameter :: istep = 2

!	----------------------------------
!	first phase, apply diagonal blocks
!	----------------------------------
	need_print = (mypid.eq.root).and.(idebug.ge.1)

	if (istep .le. 0) then
		return
	endif

	do iblock=1,nblocks
	do i=1,msize
	   x1(i,iblock) = rhs(i,iblock)
	enddo
	enddo

	call adi_dblock( msize,nblocks, A,B,C, x1, context )

	if (istep .le. 1) then
!$omp      parallel do private(i,iblock)
	   do iblock=1,nblocks
	   do i=1,msize
	     rhs(i,iblock) = x1(i,iblock)
	   enddo
	   enddo

	   return
	endif

!	----------------
!	compute residual
!	----------------
        call adi_matvec( msize,nblocks, A, B, C, x1, rvec,
     &          context, cal_dblock)                                            
	if (cal_dblock) then
!$omp      parallel do private(i,iblock)
	   do iblock=1,nblocks
	   do i=1,msize
		rvec(i,iblock) = rhs(i,iblock) - rvec(i,iblock)
	   enddo
	   enddo
	else
!$omp      parallel do private(i,iblock)
	   do iblock=1,nblocks
	   do i=1,msize
	       rvec(i,iblock) = -rvec(i,iblock)
	   enddo
	   enddo
	endif

!	------------------------------------
!	second correction
!	solve (diag(A) + B + C) * dx = rvec
!	------------------------------------

  	call adi_tridiag( msize,nblocks, A,B,C, rvec,  context )  

!	---------------
!	add correction
!	---------------
!$omp   parallel do private(i,iblock)
	do iblock=1,nblocks
	do i=1,msize
	   rhs(i,iblock) = x1(i,iblock) + rvec(i,iblock)
	enddo
	enddo

	return
	end subroutine adi_precond



        real function adi_dotprod_safe( msize,nblocks,
     &          x,y,context )
	implicit none
	integer, intent(in) :: msize,nblocks,context
	real, intent(in), dimension(msize,nblocks) :: x,y

!	---------------
!	local variables
!	---------------
	real, parameter :: zero = 0.0
	real, allocatable, dimension(:) :: inx,outx,iny,outy
	real :: ddotprod,dvalue
	integer :: i,iblock, op,count,datatype,ierror,
     &			total_size,isize,ioff


	isize = msize*nblocks
	count = 1

	datatype = MPI_INTEGER
        op = MPI_SUM
        call MPI_ALLREDUCE( isize, total_size, count,
     &                  datatype, op, context, ierror)

	allocate( inx(total_size), outx( total_size ), 
     &		iny(total_size), outy(total_size) )
	inx(:) = zero
	outx(:) = zero
	iny(:) = zero
	outy(:) = zero

	do iblock=1,nblocks
	ioff = mypid * (msize*nblocks) + (iblock-1)*msize
	inx( (ioff+1):(ioff+msize) ) =  x(1:msize,iblock)
	iny( (ioff+1):(ioff+msize) ) =  y(1:msize,iblock)
	enddo
	
	datatype = mpi_datatype
        op = MPI_SUM
	count = total_size
        call MPI_ALLREDUCE( inx, outx, count,
     &                  datatype, op, context, ierror)
        call MPI_ALLREDUCE( iny, outy, count,
     &                  datatype, op, context, ierror)

	
	ddotprod = zero
	do i=1,total_size
	  ddotprod = ddotprod + outx(i)*outy(i)
	enddo


	deallocate( outy, iny, outx, inx )
	
	adi_dotprod_safe = ddotprod

	return
	end function adi_dotprod_safe
     

	real function adi_dotprod( msize,nblocks, x, y, context )
	implicit none
	integer, intent(in) :: msize,nblocks,context
	real, intent(in), dimension(msize,nblocks) :: x,y


!	---------------
!	local variables
!	---------------
	real :: ddotprod,dvalue
	integer :: i,iblock, op,count,datatype,ierror
	logical, parameter :: use_safe = .false.

	if (use_safe) then
	  ddotprod = adi_dotprod_safe( msize, nblocks,
     &			x,y,context )
          adi_dotProd = ddotprod
	else 

	ddotprod = zero
!$omp   parallel do private(iblock,i) reduction(+:ddotprod)
	do iblock=1,nblocks
	do i=1,msize
	   ddotprod = ddotprod + x(i,iblock) * y(i,iblock)
	enddo
	enddo


	dvalue = -999.9
	count = 1
	datatype = mpi_datatype
	op = MPI_SUM
	call MPI_ALLREDUCE( ddotprod, dvalue, count, 
     &			datatype, op, context, ierror)
	call assert( ierror .eq. MPI_SUCCESS,
     &		'adi_dotprod: mpi_allreduce returns ', ierror)


	adi_dotprod = dvalue

	endif
	return
	end function adi_dotprod






	subroutine adi_vnorm( msize,nblocks, x, vnorm, context )
	implicit none
        integer, intent(in) :: msize,nblocks,context
        real, intent(in), dimension(msize,nblocks) :: x
	real, intent(inout) :: vnorm
	

	vnorm = adi_dotprod( msize, nblocks, x,x, context )
	vnorm = sqrt( abs(vnorm) )
	return
	end subroutine adi_vnorm


	subroutine adi_fixpoint( msize,nblocks, A,B,C,  
     &		rhs, niter, tol, context)
	implicit none
	integer, intent(in) :: msize,nblocks
	integer, intent(in) :: niter,context
	real, intent(in) :: tol
	real, intent(in), dimension(msize,msize,nblocks) :: A
	real, intent(in), dimension(msize,nblocks) :: B,C
	real, intent(inout), dimension(msize,nblocks) :: rhs

!	------------------------------------------------
!	if tol is exact zero, then use only niter
!	as termination, otherwise, do both niter and tol
!	------------------------------------------------


!	---------------
!	local variables
!	---------------
	logical, parameter :: cal_dblock = .true.
	real, dimension(msize,nblocks) :: bvec, x
	logical :: need_norm, isconverged
	real, dimension(0:niter) :: resid_norm
	real :: vnorm
	integer :: i,iblock,k
	logical :: need_print




!	-----------------------------
!	perform initial factorization
!	-----------------------------
        call adi_init( msize, nblocks, A, B, C, context )
	need_print = (mypid.eq.root).and.(idebug.ge.1)


	need_norm = (tol .gt. zero)

!$omp   parallel do private(i,iblock)
	do iblock=1,nblocks
	do i=1,msize
	  bvec(i,iblock) = rhs(i,iblock)
	enddo
	enddo

	k = 0
        if (need_norm) then
                call adi_vnorm( msize,nblocks, bvec, vnorm, context )
                resid_norm(k) = vnorm

                if (need_print) then
                  write(*,*) 'adi_fixpoint: k, resid_norm(k) ',
     &                  k,resid_norm(k)                                         
		endif
        endif                                                                 

        call adi_precond( msize,nblocks, A,B,C, bvec,  context )      

!$omp parallel do private(i,iblock)
	do iblock=1,nblocks
	do i=1,msize
	   x(i,iblock) = bvec(i,iblock)
	enddo
	enddo

	do k=1,niter

!	  -----------------
!	  compute  residual
!	  -----------------
          call adi_matvec( msize,nblocks, A, B, C, x, bvec,
     &          context, cal_dblock)                                            

!$omp parallel do private(i,iblock)
	  do iblock=1,nblocks
	  do i=1,msize
	    bvec(i,iblock) = rhs(i,iblock) - bvec(i,iblock)
	  enddo
	  enddo

	  if (need_norm) then
        	call adi_vnorm( msize,nblocks, bvec, vnorm, context )         
	        resid_norm(k) = vnorm

	        if (need_print) then
		  write(*,*) 'adi_fixpoint: k, resid_norm(k) ',
     &			k,resid_norm(k)
		endif

		isconverged = resid_norm(k) .lt. tol* resid_norm(0)
		if (isconverged) exit
	  endif

!	  --------------------
!	  apply preconditioner
!	  perform correction
!	  --------------------
          call adi_precond( msize,nblocks, A,B,C, bvec,  context )      

!$omp parallel do private(i,iblock)
	  do iblock=1,nblocks
	  do i=1,msize
	      x(i,iblock) = x(i,iblock) + bvec(i,iblock)
	  enddo
	  enddo
	enddo

!$omp parallel do private(i,iblock)
	do iblock=1,nblocks
	do i=1,msize
	  rhs(i,iblock) = x(i,iblock)
	enddo
	enddo


	call adi_cleanup()


	return
	end subroutine adi_fixpoint

	

        subroutine adi_bicgstab( msize,nblocks, A,B,C,
     &          rhs, niter, tol, context)                                       

        implicit none
        integer, intent(in) :: msize,nblocks
        integer, intent(in) :: niter,context
        real, intent(in) :: tol
        real, intent(in), dimension(msize,msize,nblocks) :: A
        real, intent(in), dimension(msize,nblocks) :: B,C
        real, intent(inout), dimension(msize,nblocks) :: rhs
 
!       ---------------
!       local variables
!       ---------------
        logical, parameter :: cal_dblock = .true.
        real, dimension(msize,nblocks) ::  x,svec,shat,pvec,phat
        real, dimension(msize,nblocks) ::  rvec,rtilde,vvec,tvec
	real, dimension(0:niter) :: rho,beta,alpha,omega
        logical :: is_converged, is_small_enough
        real, dimension(0:niter) :: resid_norm
        real :: rtilde_times_vvec,tvec_times_svec,tvec_times_tvec
	real :: rvec_norm,svec_norm
        integer :: i,iblock,k
	logical :: need_print

	real :: dtiny  
 


!	-----------------------------
!	perform initial factorization
!	-----------------------------
        call adi_init( msize, nblocks, A, B, C, context )
	dtiny = 10.0*tiny(one)
	need_print = (mypid.eq.root).and.(idebug.ge.1)



	resid_norm(:) = zero
	rho(:) = zero
	alpha(:) = zero
	beta(:) = zero
	omega(:) = zero
 
!$omp   parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize

	svec(i,iblock) = zero
	shat(i,iblock) = zero
	pvec(i,iblock) = zero
	phat(i,iblock) = zero
	rvec(i,iblock) = zero
	rtilde(i,iblock) = zero
	vvec(i,iblock) = zero
	tvec(i,iblock) = zero

        enddo
        enddo
 

!	---------------------------------------------------
!	Compute r_0 = b - A*x_0, for some initial guess x_0
!	Choose rtilde (for example rtilde = r_0)
!	---------------------------------------------------

!$omp   parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize

	x(i,iblock) = zero
	rvec(i,iblock) = rhs(i,iblock)
	rtilde(i,iblock) = rvec(i,iblock)

        enddo
        enddo


!	------------------------
!	compute initial residual
!	------------------------
	k = 0
	call adi_vnorm( msize,nblocks, rvec, rvec_norm, context )         
	resid_norm(k) = rvec_norm

	if (need_print) then
		write(*,*) 'adi_bicgstab: k,resid_norm(k) ',
     &			k,resid_norm(k)
	endif


	do k=1,niter

	  rho(k-1) = adi_dotprod( msize,nblocks, 
     &		rtilde, rvec, context )
	  if (abs(rho(k-1)) .lt. dtiny ) then
		write(*,*) 
     &		   'adi_bicgstab: rho(k-1) too small ', 
     &			rho(k-1), ' method fails'
		exit
	  endif

	  if (k.eq.1) then

!$omp         parallel do private(i,iblock)
              do iblock=1,nblocks
              do i=1,msize
		pvec(i,iblock) = rvec(i,iblock)
              enddo
              enddo
	  else
		beta(k-1) = ( rho(k-1)/rho(k-2) ) * 
     &			(alpha(k-1)/omega(k-1))

!$omp        parallel do private(i,iblock)
             do iblock=1,nblocks 
             do i=1,msize
		pvec(i,iblock) = rvec(i,iblock) + 
     &		  beta(k-1)*
     &             (pvec(i,iblock) - omega(k-1)*vvec(i,iblock))
             enddo
             enddo
	  endif

!	----------------------
!	solve M phat = p_(i)
!	----------------------
!$omp    parallel do private(i,iblock)
         do iblock=1,nblocks
         do i=1,msize
	  phat(i,iblock) = pvec(i,iblock)
         enddo
         enddo

          call adi_precond( msize,nblocks, A,B,C, phat,  context )              


!	------------------
!	v_(i) = A * phat
!	------------------
          call adi_matvec( msize,nblocks, A, B, C, phat, vvec,
     &          context, cal_dblock)                                            


	rtilde_times_vvec = adi_dotprod(msize,nblocks,
     &			rtilde,vvec,context)
	alpha(k) = rho(k-1) /  rtilde_times_vvec

!$omp   parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize
	  svec(i,iblock) = rvec(i,iblock) - 
     &                   alpha(k) * vvec(i,iblock)
        enddo
        enddo

!	-----------------------------------------------------------------
!	check norm of svec; 
!	if small enough: set x_(k) = x_(i-1) + alpha_i * phat and stop
!	-----------------------------------------------------------------

	call adi_vnorm( msize,nblocks, svec, svec_norm, context )         
	is_small_enough = svec_norm .lt. dtiny
	if (is_small_enough) then

!$omp      parallel do private(i,iblock)
           do iblock=1,nblocks
           do i=1,msize
	     x(i,iblock) = x(i,iblock) + alpha(k) * phat(i,iblock)
           enddo
           enddo

	   if (need_print) then
		write(*,*) 'svec is small enough ', svec_norm
	   endif

	   exit
	endif


!	-------------------
!	solve M shat = svec
!	-------------------
!$omp   parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize
	   shat(i,iblock) = svec(i,iblock)
        enddo
        enddo

        call adi_precond( msize,nblocks, A,B,C, shat,  context )              

!	----------------
!	tvec  = A * shat
!	----------------
          call adi_matvec( msize,nblocks, A, B, C, shat, tvec,
     &          context, cal_dblock)                                            


	tvec_times_svec = adi_dotprod( msize,nblocks, 
     &			tvec, svec, context )
	tvec_times_tvec = adi_dotprod( msize,nblocks,
     &			tvec, tvec, context )
	omega(k) = tvec_times_svec / tvec_times_tvec


!$omp   parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize
	   x(i,iblock) = x(i,iblock) + 
     &          alpha(k) * phat(i,iblock) + 
     &		omega(k) * shat(i,iblock)
	   rvec(i,iblock) = svec(i,iblock) - 
     &          omega(k) * tvec(i,iblock)
        enddo
        enddo

!	--------------------------------------------
!	check for convergence; continue if necessary
!	--------------------------------------------
	call adi_vnorm( msize,nblocks, rvec, rvec_norm, context )         
	resid_norm(k) = rvec_norm

	if (need_print) then
		write(*,*) 'adi_bicgstab: k,resid_norm(k) ',
     &			k,resid_norm(k)
	endif

	is_converged = (resid_norm(k) .lt. tol * resid_norm(0))
	if (is_converged) then
	   if (need_print) then
	      write(*,*) 'convergence, k,resid_norm(k) ', 
     &			k, resid_norm(k)
	   endif
	   exit
	endif

	if (abs(omega(k)).lt.dtiny) then
	   if (need_print) then
		write(*,*) 'omega(k) too small, k,omega(k) ',
     &				k,omega(k)
	   endif
	   exit
	endif

	enddo


!	---------------
!	copy result out
!	---------------
!$omp   parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize
	   rhs(i,iblock) = x(i,iblock)
        enddo
        enddo

	call adi_cleanup()

	return
	end subroutine adi_bicgstab


	subroutine adi_rotmat( a,b,  c,s )
!
! ----------------------------------------------------------
! Compute the Givens rotation matrix parameters for a and b.
! ----------------------------------------------------------
	implicit none
	real, intent(in) :: a,b
	real, intent(inout) :: c,s

	real :: temp
	real :: dtiny

        dtiny = 10.0*tiny(b)

	if (abs(b).le.dtiny) then
	   c = one
	   s = zero
	else 
          if ( abs(b) .gt. abs(a) ) then
	   temp = a / b
	   s = one / sqrt( one + temp * temp )
	   c = temp * s
	  else
	   temp = b / a
	   c = one / sqrt( one + temp * temp )
	   s = temp * c
          endif
         endif

	 return
	 end subroutine adi_rotmat

           

        subroutine adi_gmres( msize,nblocks, A,B,C,
     &          rhs, tol, restrt, max_it, context)

        implicit none
        integer, intent(in) :: msize,nblocks
        integer, intent(in) :: context
	integer, intent(in) :: restrt, max_it
        real, intent(in) :: tol
        real, intent(in), dimension(msize,msize,nblocks) :: A
        real, intent(in), dimension(msize,nblocks) :: B,C
        real, intent(inout), dimension(msize,nblocks) :: rhs

!	---------------
!	local variables
!	---------------
	logical, parameter ::  cal_dblock = .true.
	real :: dtiny

	logical :: is_converged, need_print
	integer :: iter,i,j,k,m,n,i_save,iblock
	real, dimension(restrt+1) :: e1,svec
	real :: error,bnrm2, rvec_norm, wvec_norm
	real, dimension(restrt) :: cs,sn

	real, dimension(msize,nblocks,restrt+1) :: V
	real, dimension(msize,nblocks) :: 
     &               rvec,bvec, x,tmpvec,wvec
	real :: temp,wvec_times_V_k

	real, dimension(restrt+1,restrt+1) :: H
	real, dimension(restrt,restrt) :: H_LU
	real, dimension(restrt) :: yvec
	integer, dimension(restrt) :: H_ipiv
	integer :: info,nrhs,mm,nn,lda,ldb
	

	call assert( restrt .ge. 1,
     &		'adi_gmres: invalid restrt ', restrt )
	call assert( max_it .ge. 1,
     &		'adi_gmres: invalid max_it ', max_it )

!       -----------------------------
!       perform initial factorization
!       -----------------------------
        call adi_init( msize, nblocks, A, B, C, context )

	dtiny = 10.0*tiny(one)
        need_print = (mypid.eq.root).and.(idebug.ge.1)

	iter = 0

!$omp parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize
	  bvec(i,iblock) = rhs(i,iblock)
        enddo
        enddo

        call adi_vnorm( msize,nblocks, bvec, bnrm2, context )
        if (need_print) then
           write(*,*) 'adi_gmres: bnrm2, dtiny ', bnrm2,dtiny
	endif
	if (bnrm2 .lt. dtiny) then
	   bnrm2 = one;
	endif
	




!	--------------------------
!	rvec = M \ (bvec - A * xvec )
!	error = norm( rvec ) / bnrm2
!	if ( error < tol ) return; end
!	--------------------------


!$omp parallel do private(i,iblock)
        do iblock=1,nblocks
        do i=1,msize
	   x(i,iblock) = zero
	   rvec(i,iblock)  = bvec(i,iblock)
        enddo
        enddo


        call adi_precond( msize,nblocks, A,B,C, rvec,  context )

        call adi_vnorm( msize,nblocks, rvec, rvec_norm, context )
	error = rvec_norm / bnrm2
        if (need_print) then
	   write(*,*) 'adi_gmres: initial error ', error
	endif

	if (error .lt. tol) then

!$omp parallel do private(i,iblock)
            do iblock=1,nblocks
            do i=1,msize
	       x(i,iblock) = rvec(i,iblock)
            enddo
            enddo

	    if (need_print) then
		write(*,*) 'adi_gmres:early return '
		write(*,*) 'rvec_norm, bnrm2 ', rvec_norm, bnrm2
	    endif
	    goto 999
	endif
	

!    --------------------------------
!    initialize workspace
!
!    [n,n] = size(A);                                  
!    m = restrt;
!    V(1:n,1:m+1) = zeros(n,m+1);
!    H(1:m+1,1:m) = zeros(m+1,m);
!    cs(1:m) = zeros(m,1);
!    sn(1:m) = zeros(m,1);
!    e1    = zeros(n,1); <-- perhaps should be zeros(m,1)
!    e1(1) = 1.0;
!    --------------------------------
        n = msize * nblocks
        m = restrt
   
!$omp   parallel do private(k,j,iblock)
        do iblock=1,nblocks
        do j=1,(m+1)
        do k=1,msize
          V(k,iblock,j) = zero
        enddo
        enddo
        enddo

        H(1:(m+1),1:m) = zero
        cs(1:m) = zero
        sn(1:m) = zero
        e1(:) = zero
        e1(1) = one


!	----------------
!	main outer loop
!	----------------
	do iter = 1,max_it
     

!	---------------------
!      r = M \ ( b-A*x );
!      V(:,1) = r / norm( r );
!      s = norm( r )*e1;
!	---------------------

          call adi_matvec( msize,nblocks, A, B, C, x, rvec,
     &          context, cal_dblock)


!$omp parallel do private(i,iblock)
          do iblock=1,nblocks
          do i=1,msize
	     rvec(i,iblock) = rhs(i,iblock) - rvec(i,iblock)
          enddo
          enddo

          call adi_precond( msize,nblocks, A,B,C, rvec,  context )

          call adi_vnorm( msize,nblocks, rvec, rvec_norm, context )

!$omp     parallel do private(i,iblock)
          do iblock=1,nblocks
          do i=1,msize
	    V(i,iblock,1) = rvec(i,iblock) / rvec_norm
          enddo
          enddo

           
	  svec(:) = rvec_norm * e1(:)

!	  ---------------------------------------------
!	  contruct orthogonal basis using  Gram-Schmidt
!	  ---------------------------------------------
	  do i=1,m
	    i_save = i
!	  ------------------------
!         w = M \ (A*V(:,i));                         
!         for k = 1:i,
!           H(k,i)= w'*V(:,k);
!           w = w - H(k,i)*V(:,k);
!         end
!	  ------------------------

!$omp parallel do private(j,iblock)
              do iblock=1,nblocks
              do j=1,msize
	        tmpvec(j,iblock) = V(j,iblock,i)
              enddo
              enddo

	      call adi_matvec( msize,nblocks, A,B,C, 
     &		tmpvec, wvec, context, cal_dblock)

	      call adi_precond( msize,nblocks, A,B,C, wvec, context )

	      do k=1,i

!$omp           parallel do private(j,iblock)
                do iblock=1,nblocks
                do j=1,msize
		  tmpvec(j,iblock) = V(j,iblock,k)
                enddo
                enddo

		wvec_times_V_k = adi_dotprod( msize,nblocks, 
     &			wvec, tmpvec, context )
		H(k,i) = wvec_times_V_k

!$omp           parallel do private(j,iblock)
                do iblock=1,nblocks
                do j=1,msize
		  wvec(j,iblock) = wvec(j,iblock) - 
     &                             H(k,i) * V(j,iblock,k)
                enddo
                enddo

	      enddo

!	  ----------------------------------
!         H(i+1,i) = norm( w );
!         V(:,i+1) = w / H(i+1,i);
!         for k = 1:i-1,                              % apply Givens rotation
!            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
!            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
!            H(k,i)   = temp;
!         end
!	  ----------------------------------
	        call adi_vnorm( msize,nblocks, wvec, wvec_norm, context )
		H(i+1,i) = wvec_norm

!$omp           parallel do private(j,iblock)
                do iblock=1,nblocks
                do j=1,msize
		  V(j,iblock,i+1) = wvec(j,iblock) / H(i+1,i)
                enddo
                enddo

		do k=1,(i-1)
		   temp     =  cs(k) * H(k,i) + sn(k) * H(k+1,i)
		   H(k+1,i) = -sn(k) * H(k,i) + cs(k) * H(k+1,i)
		   H(k,i)   = temp
		enddo


!	  -------------------------------
!         % form i-th rotation matrix
!         [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); 
!
!         % approximate residual norm
!         temp   = cs(i)*s(i);                        
!         s(i+1) = -sn(i)*s(i);
!         s(i)   = temp;
!         H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
!         H(i+1,i) = 0.0;
!	  -------------------------------

	  call adi_rotmat( H(i,i), H(i+1,i),  cs(i),sn(i) )

	  temp   =  cs(i) * svec(i)
	  svec(i+1) = -sn(i)*svec(i)
	  svec(i)   = temp
	  H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i)
	  H(i+1,i) = zero


!	  -----------------------------
!         error  = abs(s(i+1)) / bnrm2;
!         if ( error <= tol ),                        % update approximation
!            y = H(1:i,1:i) \ s(1:i);                 % and exit
!            x = x + V(:,1:i)*y;
!            break;
!         end
!	  -----------------------------

	  error = abs( svec(i+1) )/ bnrm2
	  if (need_print) then
	    write(*,*) 'adi_gmres:iter,i,error ',iter,i,error
	  endif

	  if (error .le. tol) then
               H_LU(1:i,1:i) = H(1:i,1:i)
	       H_ipiv(1:i) = 0
	       mm = i
	       nn = i
	       lda = m
	       info = 0
               call dgetrf( mm,nn, H_LU, lda,
     &		   H_ipiv, info )
	       call assert( info.eq.0,
     &		'adi_gmres: dgetrf returns info ', info )


	       yvec(1:mm) = svec(1:mm)
	       nrhs = 1
 	       ldb = m
	       call dgetrs( 'N', mm, nrhs,
     &		    H_LU, lda, H_ipiv,
     &		    yvec, ldb, info )
	       call assert( info.eq.0,
     &		 'adi_gmres: dgetrs returns info ', info )

	   	do j=1,i

!$omp             parallel do private(k,iblock)
                  do iblock=1,nblocks
                  do k=1,msize
		     x(k,iblock) = x(k,iblock) + 
     &                    V(k,iblock,j) * yvec(j)
                  enddo
                  enddo

		enddo


                if (need_print .and. (idebug .ge.2)) then
                   write(*,*) 'adi_gmres: yvec '
                   do j=1,i
                     write(*,*) yvec(j)
                   enddo
                endif

		exit
	    endif
         enddo         ! end do i

!	------------------------------------------------------
!	 fix possible indexing error in gmres.m templates code
!	 set i=i_save to some reasonable value
!	 outside of do loop, i can be undefined
!	------------------------------------------------------
         i = i_save


!	-------------------------------
!      if ( error <= tol ), break, end
!      y = H(1:m,1:m) \ s(1:m);
!      x = x + V(:,1:m)*y;                            % update approximation
!	-------------------------------

	if (error .le. tol ) then
           if (need_print) then
	       write(*,*) 'adi_gmres: early exit, i,iter ',i,iter
	       write(*,*) 'error, bnrm2 ', error,bnrm2
	   endif

	   exit
	endif

	H_LU(1:m,1:m) = H(1:m,1:m)
	H_ipiv(1:m) = 0
	mm = m
	nn = m
	lda = m
	info = 0
	call dgetrf( mm,nn,H_LU,lda,H_ipiv,info)
	call assert(info.eq.0,
     &		'adi_gmres: dgetrf return info = ', info )

	nrhs = 1
	ldb = m
	yvec(1:m) = svec(1:m)
	call dgetrs( 'N', mm,nrhs, H_LU,lda,H_ipiv,
     &		yvec, ldb, info )
	call assert( info.eq.0,
     &		'adi_gmres: dgetrs return info = ', info )

	do j=1,m

!$omp     parallel do private(k,iblock)
          do iblock=1,nblocks
          do k=1,msize
	     x(k,iblock) = x(k,iblock) + 
     &             V(k,iblock,j) * yvec(j)
          enddo
          enddo

	enddo

        if (need_print .and. (idebug .ge. 2)) then
           write(*,*) 'adi_gmres: yvec '
	   do j=1,m
	     write(*,*) yvec(j)
	   enddo
	endif


!	---------------
!      r = M \ ( b-A*x )                              % compute residual
!      s(i+1) = norm(r);
!	---------------

	
	  call adi_matvec( msize,nblocks, A,B,C,
     &          x, rvec, context, cal_dblock)


!$omp     parallel do private(k,iblock)
          do iblock=1,nblocks
          do k=1,msize
	    rvec(k,iblock)  = rhs(k,iblock) - rvec(k,iblock)
          enddo
          enddo

	  call adi_precond( msize,nblocks, A,B,C, rvec, context )

	  call adi_vnorm( msize,nblocks, rvec, rvec_norm, context )
	  
          if (i+1.gt.ubound(svec,1)) then
	     write(*,*) 'out of bound error '
	     write(*,*) 'i,m,restrt ',i,m,restrt
	  endif
	  svec(i+1) = rvec_norm

!	----------------------------
!      error = s(i+1) / bnrm2;                        % check convergence
!      if ( error <= tol ), break, end;
!	----------------------------

	error = svec(i+1) / bnrm2
        if (need_print) then
	   write(*,*) 'adi_gmres: iter,i,error ',iter,i,error
	endif



	if ( error .le. tol) then
	   if (need_print) then
              write(*,*) 'adi_gmres: early exit in iter loop '
	      write(*,*) 'i,iter,error,bnrm2 ',i,iter,error,bnrm2
	   endif
	      
	   exit
	endif

      enddo   


      is_converged = (error .le. tol)
      if (.not.is_converged) then
	if (mypid.eq.root) then
	  write(*,*) 
     &	    'adi_gmres: no convergence, error = ', error
          write(*,*) 'restrt,max_it = ', restrt,max_it
        endif
      endif


 999   continue
!	-------------
!       copy result out
!	-------------

!$omp parallel do private(k,iblock)
      do iblock=1,nblocks
      do k=1,msize
         rhs(k,iblock) = x(k,iblock)
      enddo
      enddo

      call adi_cleanup()
      return
      end subroutine adi_gmres


	end module  adi_module
