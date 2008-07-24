subroutine shifti
  use global_parameters
  use particle_array
  use particle_decomp
  use particle_tracking
  implicit none
  
  integer i,m,msendleft(2),msendright(2),mrecvleft(2),mrecvright(2),mtop,&
       kzi(mimax),m0,msend,mrecv,idest,isource,isendtag,irecvtag,nzion,&
       iright(mimax),ileft(mimax),isendcount,irecvcount,&
       istatus(MPI_STATUS_SIZE),ierror,iteration,lasth
  real(wp),dimension(:,:),allocatable :: recvbuf,sendleft,sendright
  real(wp) zetaright,zetaleft,pi_inv
  character(len=8) cdum

  real(wp) aright,aleft,alpha,beta,kappa,pi2sq_inv,pi2,pi4_inv
#ifdef _SX
  integer nflag,imm,msend_r,msend_l
#else
  integer nflag,imm,msend_r(4),msend_l(4)
  integer kzi_r(mimax/4,4),kzi_l(mimax/4,4)
#endif

  if(track_particles > 0)then   ! Keep track of particles
     nzion=14
  else
     nzion=12
  endif
  if(nzion /= 2*nparam)then
  ! nzion and 2*nparam have to be equal so we do a quick check.
  ! We use nzion instead of nparam to let the vectorizing compiler
  ! know what this value is for optimization purposes.
    write(0,*)'Error in shifti: nzion not equal to nparam!'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif
  pi_inv=1.0_wp/pi
  m0=1
  iteration=0

! Quantities used in the vectorized version of the main loop
  pi2=2.0_wp*pi
  pi2sq_inv=1.0_wp/(pi2*pi2)
  pi4_inv=1.0_wp/(4.0_wp*pi)
  
100 iteration=iteration+1
  if(iteration>ntoroidal)then
     write(0,*)'endless particle sorting loop at PE=',mype
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  msend=0
  msendright=0
  msendleft=0

  if(m0 <= mi)then
     !do m=m0,mi
     do m=1,mi
        kzi(m)=0
     enddo

#ifdef _SX
! ******************************************************************
! Section for the NEC SX type of vector processor. It should work on any
! vector machine.
! ******************************************************************
     msend_r=0
     msend_l=0

  CALL FTRACE_REGION_BEGIN("SHIFTI_1")

    !!! do m=m0,mi
     do m=1,mi
        zetaright=min(2.0*pi,zion(3,m))-zetamax
        zetaleft=zetamin-zion(3,m)

        alpha=pi2*aint(1.0-pi4_inv*zetaleft)
        beta=pi2*aint(1.0-pi4_inv*zetaright)
        kappa=pi2*aint(1.0+zetaleft*zetaright*pi2sq_inv)

        aright=(alpha+zetaleft) - (beta+zetaright) - kappa
        aleft=(alpha+zetaleft) - (beta+zetaright) + kappa

        if( aright > 0.0 )then
! # of particle to move right
           msend_r=msend_r+1
           iright(msend_r)=m
        endif
     enddo

     msendright(1)=msend_r

     do m=m0,mi
        zetaright=min(2.0*pi,zion(3,m))-zetamax
        zetaleft=zetamin-zion(3,m)

        alpha=pi2*aint(1.0-pi4_inv*zetaleft)
        beta=pi2*aint(1.0-pi4_inv*zetaright)
        kappa=pi2*aint(1.0+zetaleft*zetaright*pi2sq_inv)

        aright=(alpha+zetaleft) - (beta+zetaright) - kappa
        aleft=(alpha+zetaleft) - (beta+zetaright) + kappa

        if( aleft < 0.0 )then
! # of particle to move left
           msend_l=msend_l+1
           ileft(msend_l)=m
        endif
     enddo

     msendleft(1)=msend_l

  CALL FTRACE_REGION_END("SHIFTI_1")

     if(nhybrid == 0) then
       nflag = 0
       do imm = 1,msendright(1)
         m = iright(imm)
         if(m == ntracer)then
           msendright(2) = m
           nflag = 1
         endif
       enddo
       do imm = 1,msendleft(1)
         m = ileft(imm)
         if(m == ntracer)then
           msendleft(2) = m
           nflag = 1
         endif
       enddo
       if(nflag == 1) ntracer = 0
     endif

     msend=0
     i=1
     do m=1,msend_r
        do while (ileft(i) < iright(m))
           msend=msend+1
           kzi(msend)=ileft(i)
           i=i+1
           if(i > msend_l)then
           ! Set kzi_l(i,imm) to value > than kzi_r(msend_r(imm),imm)
             ileft(i)=iright(msend_r)+1000
           endif
        enddo
        msend=msend+1
        kzi(msend)=iright(m)
     enddo
   ! If ileft(i) > iright(msend_r) and i < msend_l
     do m=i,msend_l
        msend=msend+1
        kzi(msend)=ileft(m)
     enddo

#else
! ******************************************************************
! Section specific to the CRAY-X1 in MSP mode (for multi-streaming)
! ******************************************************************
     msend_r=0
     msend_l=0

!dir$ preferstream
     do imm=1,4
!dir$ prefervector
        do m=(imm-1)*mi/4+1,imm*mi/4
           zetaright=min(2.0*pi,zion(3,m))-zetamax
           zetaleft=zetamin-zion(3,m)

           alpha=pi2*aint(1.0-pi4_inv*zetaleft)
           beta=pi2*aint(1.0-pi4_inv*zetaright)
           kappa=pi2*aint(1.0+zetaleft*zetaright*pi2sq_inv)

           aright=(alpha+zetaleft) - (beta+zetaright) - kappa
           aleft=(alpha+zetaleft) - (beta+zetaright) + kappa

           if( aright > 0.0 )then
! # of particle to move right
              msend_r(imm)=msend_r(imm)+1
              kzi_r(msend_r(imm),imm)=m
           endif

           if( aleft < 0.0 )then
! # of particle to move left
              msend_l(imm)=msend_l(imm)+1
              kzi_l(msend_l(imm),imm)=m
           endif
        enddo
     enddo

     msendright=0
     msendleft=0
     do imm=1,4
        do m=1,msend_r(imm)
           iright(msendright(1)+m)=kzi_r(m,imm)
        enddo
        msendright(1)=msendright(1)+msend_r(imm)
        do m=1,msend_l(imm)
           ileft(msendleft(1)+m)=kzi_l(m,imm)
        enddo
        msendleft(1)=msendleft(1)+msend_l(imm)
     enddo

     msend=0
     do imm=1,4
        i=1
        do m=1,msend_r(imm)
           do while (kzi_l(i,imm) < kzi_r(m,imm))
              msend=msend+1
              kzi(msend)=kzi_l(i,imm)
              i=i+1
              if(i > msend_l(imm))then
              ! Set kzi_l(i,imm) to value > than kzi_r(msend_r(imm),imm)
                kzi_l(i,imm)=kzi_r(msend_r(imm),imm)+1000
              endif
           enddo
           msend=msend+1
           kzi(msend)=kzi_r(m,imm)
        enddo
      ! If kzi_l(i,imm) > kzi_r(msend_r(imm),imm) and i < msend_l(imm)
        do m=i,msend_l(imm)
           msend=msend+1
           kzi(msend)=kzi_l(m,imm)
        enddo
     enddo

     if(nhybrid == 0) then
       nflag = 0
       do m = 1,msendright(1)
         if(iright(m) == ntracer)then
           msendright(2) = m
           nflag = 1
         endif
       enddo
       do m = 1,msendleft(1)
         if(ileft(m) == ntracer)then
           msendleft(2) = m
           nflag = 1
         endif
       enddo
       if(nflag == 1) ntracer = 0
     endif

#endif

  endif  ! if(m0 <= mi)then

  if(iteration>1)then
! total # of particles to be shifted
     mrecv=0
     msend=msendleft(1)+msendright(1)

     call MPI_ALLREDUCE(msend,mrecv,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)

! no particle to be shifted, return
     if ( mrecv == 0 ) then
        return
     endif
  endif

! an extra space to prevent zero size when msendright(1)=msendleft(1)=0
  allocate(sendright(nzion,max(msendright(1),1)),sendleft(nzion,max(msendleft(1),1)))

! pack particle to move right
     if(track_particles > 0)then
        do m=1,msendright(1)
           sendright(1,m)=zion(1,iright(m))
           sendright(2,m)=zion(2,iright(m))
           sendright(3,m)=zion(3,iright(m))
           sendright(4,m)=zion(4,iright(m))
           sendright(5,m)=zion(5,iright(m))
           sendright(6,m)=zion(6,iright(m))
           sendright(7,m)=zion(7,iright(m))

           sendright(8,m)=zion0(1,iright(m))
           sendright(9,m)=zion0(2,iright(m))
           sendright(10,m)=zion0(3,iright(m))
           sendright(11,m)=zion0(4,iright(m))
           sendright(12,m)=zion0(5,iright(m))
           sendright(13,m)=zion0(6,iright(m))
           sendright(14,m)=zion0(7,iright(m))
        enddo
     else
        do m=1,msendright(1)
           sendright(1,m)=zion(1,iright(m))
           sendright(2,m)=zion(2,iright(m))
           sendright(3,m)=zion(3,iright(m))
           sendright(4,m)=zion(4,iright(m))
           sendright(5,m)=zion(5,iright(m))
           sendright(6,m)=zion(6,iright(m))

           sendright(7,m)=zion0(1,iright(m))
           sendright(8,m)=zion0(2,iright(m))
           sendright(9,m)=zion0(3,iright(m))
           sendright(10,m)=zion0(4,iright(m))
           sendright(11,m)=zion0(5,iright(m))
           sendright(12,m)=zion0(6,iright(m))
        enddo
     endif

! pack particle to move left
     if(track_particles > 0)then
        do m=1,msendleft(1)    
           sendleft(1,m)=zion(1,ileft(m))
           sendleft(2,m)=zion(2,ileft(m))
           sendleft(3,m)=zion(3,ileft(m))
           sendleft(4,m)=zion(4,ileft(m))
           sendleft(5,m)=zion(5,ileft(m))
           sendleft(6,m)=zion(6,ileft(m))
           sendleft(7,m)=zion(7,ileft(m))

           sendleft(8,m)=zion0(1,ileft(m))
           sendleft(9,m)=zion0(2,ileft(m))
           sendleft(10,m)=zion0(3,ileft(m))
           sendleft(11,m)=zion0(4,ileft(m))
           sendleft(12,m)=zion0(5,ileft(m))
           sendleft(13,m)=zion0(6,ileft(m))
           sendleft(14,m)=zion0(7,ileft(m))
        enddo
     else
        do m=1,msendleft(1)    
           sendleft(1,m)=zion(1,ileft(m))
           sendleft(2,m)=zion(2,ileft(m))
           sendleft(3,m)=zion(3,ileft(m))
           sendleft(4,m)=zion(4,ileft(m))
           sendleft(5,m)=zion(5,ileft(m))
           sendleft(6,m)=zion(6,ileft(m))

           sendleft(7,m)=zion0(1,ileft(m))
           sendleft(8,m)=zion0(2,ileft(m))
           sendleft(9,m)=zion0(3,ileft(m))
           sendleft(10,m)=zion0(4,ileft(m))
           sendleft(11,m)=zion0(5,ileft(m))
           sendleft(12,m)=zion0(6,ileft(m))
        enddo
     endif

!  endif

! send # of particle to move right to neighboring PEs of same particle
! domain.
  mrecvleft=0
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendright,2,MPI_INTEGER,idest,isendtag,&
       mrecvleft,2,MPI_INTEGER,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! send # of particle to move left
  mrecvright=0
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendleft,2,MPI_INTEGER,idest,isendtag,&
       mrecvright,2,MPI_INTEGER,isource,irecvtag,toroidal_comm,istatus,ierror)
 
! Total number of particles to receive
  mrecv=mrecvleft(1)+mrecvright(1)

! need extra particle array
  if(mi-msend+mrecv > mimax)then
     write(0,*)"need bigger particle array",mimax,mi+mrecvleft(1)+mrecvright(1)
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
! We assume that there is not enough memory available to this process to
! allow the allocation of a temporary array to hold the particle data.
! We use a disk file instead although this slows down the code.
! open disk file
!     write(cdum,'("TEMP.",i5.5)')mype
!     open(111,file=cdum,status='replace',form='unformatted')
!    
! record particle information to file
!     write(111)zion(1:nparam,1:mi)
!     write(111)zion0(1:nparam,1:mi)
!    
! make bigger array
!     deallocate(zion,zion0)
!     mimax=2*(mi+mrecvleft(1)+mrecvright(1))-mimax
!     allocate(zion(nparam,mimax),zion0(nparam,mimax))
!    
! read back particle information
!     rewind(111)
!     read(111)zion(1:nparam,1:mi)
!     read(111)zion0(1:nparam,1:mi)
!     close(111)
  endif

  !  write(0,*)mype,' msend=',msend,'  mrecv=',mrecv

! Allocate receive buffer
  allocate(recvbuf(nzion,max(mrecv+1,1)))
 
! send particle to right and receive from left
  recvbuf=0.0_wp
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  isendcount=msendright(1)*nzion
  irecvcount=mrecvleft(1)*nzion
  call MPI_SENDRECV(sendright,isendcount,mpi_Rsize,idest,isendtag,recvbuf,&
       irecvcount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! send particle to left and receive from right
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  isendcount=msendleft(1)*nzion
  irecvcount=mrecvright(1)*nzion
  call MPI_SENDRECV(sendleft,isendcount,mpi_Rsize,idest,isendtag,&
       recvbuf(1,mrecvleft(1)+1),irecvcount,mpi_Rsize,isource,irecvtag,&
       toroidal_comm,istatus,ierror)
 
   !write(0,*)mype,' istep=',istep,'  irk=',irk,'  After receiving buffer'

! The particle array zion() now has "msend" holes at different positions
! in the array. We want to fill these holes with the received particles.
! We need to check for the possibility of having more holes than received
! particles or vice versa. The positions of the holes have been saved in
! the kzi() vector.
  if(track_particles > 0)then
!dir$ concurrent
     do m=1,min(msend,mrecv)
        zion(1,kzi(m))=recvbuf(1,m)
        zion(2,kzi(m))=recvbuf(2,m)
        zion(3,kzi(m))=recvbuf(3,m)
        zion(4,kzi(m))=recvbuf(4,m)
        zion(5,kzi(m))=recvbuf(5,m)
        zion(6,kzi(m))=recvbuf(6,m)
        zion(7,kzi(m))=recvbuf(7,m)

        zion0(1,kzi(m))=recvbuf(8,m)
        zion0(2,kzi(m))=recvbuf(9,m)
        zion0(3,kzi(m))=recvbuf(10,m)
        zion0(4,kzi(m))=recvbuf(11,m)
        zion0(5,kzi(m))=recvbuf(12,m)
        zion0(6,kzi(m))=recvbuf(13,m)
        zion0(7,kzi(m))=recvbuf(14,m)
     enddo
  else
!dir$ concurrent
     do m=1,min(msend,mrecv)
        zion(1,kzi(m))=recvbuf(1,m)
        zion(2,kzi(m))=recvbuf(2,m)
        zion(3,kzi(m))=recvbuf(3,m)
        zion(4,kzi(m))=recvbuf(4,m)
        zion(5,kzi(m))=recvbuf(5,m)
        zion(6,kzi(m))=recvbuf(6,m)

        zion0(1,kzi(m))=recvbuf(7,m)
        zion0(2,kzi(m))=recvbuf(8,m)
        zion0(3,kzi(m))=recvbuf(9,m)
        zion0(4,kzi(m))=recvbuf(10,m)
        zion0(5,kzi(m))=recvbuf(11,m)
        zion0(6,kzi(m))=recvbuf(12,m)
     enddo
  endif

! tracer particle
  if( nhybrid==0 )then
     if( mrecvleft(2) > 0 .and. mrecvleft(2) <= min(msend,mrecv))then
        ntracer=kzi(mrecvleft(2))
     elseif( mrecvright(2) > 0 .and. (mrecvleft(1)+mrecvright(2)) <= min(msend,mrecv))then
        ntracer=kzi(mrecvleft(1)+mrecvright(2))
     endif
  endif

! If mrecv > msend, we have some leftover received particles after having
! filled all the holes. We place these extra particles in the empty space
! at the top of the particle array.
  if(mrecv > msend)then
     if(track_particles > 0)then
        do m=1,(mrecv-msend)
           zion(1,mi+m)=recvbuf(1,msend+m)
           zion(2,mi+m)=recvbuf(2,msend+m)
           zion(3,mi+m)=recvbuf(3,msend+m)
           zion(4,mi+m)=recvbuf(4,msend+m)
           zion(5,mi+m)=recvbuf(5,msend+m)
           zion(6,mi+m)=recvbuf(6,msend+m)
           zion(7,mi+m)=recvbuf(7,msend+m)

           zion0(1,mi+m)=recvbuf(8,msend+m)
           zion0(2,mi+m)=recvbuf(9,msend+m)
           zion0(3,mi+m)=recvbuf(10,msend+m)
           zion0(4,mi+m)=recvbuf(11,msend+m)
           zion0(5,mi+m)=recvbuf(12,msend+m)
           zion0(6,mi+m)=recvbuf(13,msend+m)
           zion0(7,mi+m)=recvbuf(14,msend+m)
        enddo
     else
        do m=1,(mrecv-msend)
           zion(1,mi+m)=recvbuf(1,msend+m)
           zion(2,mi+m)=recvbuf(2,msend+m)
           zion(3,mi+m)=recvbuf(3,msend+m)
           zion(4,mi+m)=recvbuf(4,msend+m)
           zion(5,mi+m)=recvbuf(5,msend+m)
           zion(6,mi+m)=recvbuf(6,msend+m)

           zion0(1,mi+m)=recvbuf(7,msend+m)
           zion0(2,mi+m)=recvbuf(8,msend+m)
           zion0(3,mi+m)=recvbuf(9,msend+m)
           zion0(4,mi+m)=recvbuf(10,msend+m)
           zion0(5,mi+m)=recvbuf(11,msend+m)
           zion0(6,mi+m)=recvbuf(12,msend+m)
        enddo
     endif
   ! tracer particle
     if( nhybrid==0 )then
        if( mrecvleft(2) > 0 .and. mrecvleft(2) > msend )then
           ntracer=mi+mrecvleft(2)-msend
        elseif( mrecvright(2) > 0 .and. (mrecvleft(1)+mrecvright(2)) > msend )then
           ntracer=mi+(mrecvleft(1)+mrecvright(2)-msend)
        endif
     endif

   ! Update value of mi, which is the # of particles remaining on local PE
     mi=mi+(mrecv-msend)

! If mrecv < msend, we still have holes in the array and we fill them by
! taking particles from the top and moving them into the holes.
  else if(mrecv < msend)then
     mtop=mi
   ! # of particles remain on local PE
     mi=mi-msend+mrecv
   ! fill the remaining holes
     lasth=msend
     do i=mrecv+1,msend
        m=kzi(i)
        if (m > mi) exit  !Break out of the DO loop if m > mi
        do while(mtop == kzi(lasth))
           mtop=mtop-1
           lasth=lasth-1
        enddo
        if( nhybrid == 0 .and. mtop == ntracer )ntracer=m
        zion(1:nparam,m)=zion(1:nparam,mtop)
        zion0(1:nparam,m)=zion0(1:nparam,mtop)
        mtop=mtop-1
        if (mtop == mi) exit  !Break out of the DO loop
     enddo
  endif

  deallocate(sendleft,sendright,recvbuf)
  !!!m0=mi-mrecvright(1)-mrecvleft(1)+1
  goto 100
  
end subroutine shifti


