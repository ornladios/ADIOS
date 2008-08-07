subroutine restart_io(iop)
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use particle_decomp
  use data_type
  use particle_tracking
  implicit none
  integer merror, mrequest, mfmode
  integer i,j,k,subsize,startidx,endidx
  integer mquantity,mflx,n_mode,mstepfinal,noutputs
  real(wp) dum
  character(*),intent(in)::iop
  INTEGER(KIND=MPI_OFFSET_KIND) mype_filesize, sum_filesize

#if ADIOS
#define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
#define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0))
#define ADIOS_READ(a,b) call adios_read(a,'b'//char(0),b)
  integer*8 buf_id, group_id
  real(wp),dimension(:),allocatable::zion0_read,zelectron0_read
  character(len=50) :: restart_fname,dirstr
#endif

  if(iop/="read" .and. iop/="write")then
     write(*,*)'*** subroutine restart_io (iop <> "read" or "write")',iop
     call MPI_ABORT(MPI_COMM_WORLD,1,merror)
     return
  endif
  !call system("mkdir "// "restart" //char(0))

  if(mype==0)write(*,*)"MFLUX,MPSI,mstepall:",mflux,mpsi,mstepall
!  write(restart_fname,'(a,i5.5,"_",i5.5,".bp")')"restart_dir/restart_",myrank_toroidal, mstepall+istep
  write(restart_fname,'(a,i5.5,".bp")')"restart_dir/restart_",myrank_toroidal
!  #ifdef __TIMER
!  call MPI_BARRIER(MPI_COMM_WORLD,merror)
!  start_time=MPI_WTIME()
!  #endif

#if ADIOS
     ! setup the element path for this node
     write(dirstr,'("/node",i5.5,"/param")')mype
     dirstr=trim(dirstr)//char(0)
!     call MPI_BARRIER(MPI_COMM_WORLD,merror)
!     start_time = MPI_WTIME()
     call adios_get_group (group_id, "restart"//char(0))
     ! set the path for all vars in the type for proper sizing
     call adios_set_path (group_id,dirstr//char(0));
     restart_fname=trim(restart_fname)//char(0)
     if (iop=="read") then
        call adios_open_read (buf_id, group_id, restart_fname)
     else
        call adios_open (buf_id, group_id, restart_fname)
     endif
     ! write the sizing paramters for both reading and writing
     ADIOS_WRITE(buf_id,partd_comm)
     ADIOS_WRITE(buf_id,mflux)
     ADIOS_WRITE(buf_id,mpsi+1)
     ADIOS_WRITE(buf_id,mzeta+1)
     ADIOS_WRITE(buf_id,nparam)
     ADIOS_WRITE(buf_id,mimax)
     ADIOS_WRITE(buf_id,mgrid)
     ADIOS_WRITE(buf_id,memax)
     ADIOS_WRITE(buf_id,nhybrid)
     ADIOS_WRITE(buf_id,2*nhybrid)
  if(iop=="write")then
     ADIOS_WRITE(buf_id,mzeta)
     ADIOS_WRITE(buf_id,mi)
     ADIOS_WRITE(buf_id,me)
     ADIOS_WRITE(buf_id,ntracer)
     ADIOS_WRITE(buf_id,etracer)
     ADIOS_WRITE(buf_id,rdtemi)
     ADIOS_WRITE(buf_id,rdteme)
     ADIOS_WRITE(buf_id,ptracer)
     ADIOS_WRITE(buf_id,pfluxpsi)
     ADIOS_WRITE(buf_id,phi00)
     ADIOS_WRITE(buf_id,phip00)
     ADIOS_WRITE(buf_id,zonali)
     ADIOS_WRITE(buf_id,zonale)
     call adios_write(buf_id,"zion0"//char(0),zion0(6,:))
     ADIOS_WRITE(buf_id,phi)
     ADIOS_WRITE(buf_id,zion)
     if(nhybrid>0)then
        call adios_write(buf_id,"zelectron0"//char(0),zelectron0(6,:))
        ADIOS_WRITE(buf_id,zelectron)
        ADIOS_WRITE(buf_id,phisave)
     endif
!!     write(*,*) "Chens debugging shit",mype
     call adios_close (buf_id)
!   #ifdef __TIMER
!     call MPI_BARRIER(MPI_COMM_WORLD,merror)
!     end_time=MPI_WTIME()
!     call MPI_REDUCE(mype_filesize,sum_filesize,1,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,merror)
!     if(mype==0)then
!        write(stdout,'("NtoM Time(s), MBytes, MB/s ",a)')iop
!        write(stdout,*)sum_filesize/1024/1024,mype_filesize/1024/1024
!        write(stdout,*)end_time-start_time,mype_filesize*numberpe/1024/1024,mype_filesize*numberpe/((end_time-start_time)*1024*1024)
!     endif
!   #endif

! S.Ethier 01/30/04 Save a copy of history.out and sheareb.out for restart
     if(mype==0 .and. istep<=mstep)then
        open(777,file='history_restart.out',status='replace')
        rewind(ihistory)
        read(ihistory,101)j
        write(777,101)j
        read(ihistory,101)mquantity
        write(777,101)mquantity
        read(ihistory,101)mflx
        write(777,101)mflx
        read(ihistory,101)n_mode
        write(777,101)n_mode
        read(ihistory,101)mstepfinal
        noutputs=mstepfinal-mstep/ndiag+istep/ndiag
        write(777,101)noutputs
        if(mype==0)write(*,*)"restart history:",noutputs  
        do i=0,(mquantity+mflx+4*n_mode)*noutputs
           read(ihistory,102)dum
           write(777,102)dum
        enddo
        close(777)

        ! Now do sheareb.out
        open(777,file='sheareb_restart.out',status='replace')
        if(istep==mstep)open(444,file='sheareb.out',status='old')
        rewind(444)
        read(444,101)j
        write(777,101)j
        read(444,101)mflx
        write(777,101)mflx

        do i=1,mpsi*noutputs*j
           read(444,102)dum
           write(777,102)dum
        enddo
        close(777)
        if(istep==mstep)close(444)
     endif

101  format(i6)
102  format(e12.6)

  else
     allocate(zion0_read(mimax))
     ADIOS_READ(buf_id,mzeta)
     ADIOS_READ(buf_id,mi)
     ADIOS_READ(buf_id,me)
     
     ADIOS_READ(buf_id,ntracer)
     ADIOS_READ(buf_id,etracer)
     ADIOS_READ(buf_id,rdtemi)
     ADIOS_READ(buf_id,rdteme)
     ADIOS_READ(buf_id,ptracer)
     ADIOS_READ(buf_id,pfluxpsi)
     ADIOS_READ(buf_id,phi00)
     ADIOS_READ(buf_id,phip00)
     ADIOS_READ(buf_id,zonali)
     ADIOS_READ(buf_id,zonale)
     call adios_read(buf_id,"zion0"//char(0),zion0_read)
     ADIOS_READ(buf_id,phi)
     ADIOS_READ(buf_id,zion)
     if(nhybrid>0)then
      !  allocate(zelectron0_read(memax))
        call adios_read(buf_id,"zelectron0"//char(0),zion0_read)
        ADIOS_READ(buf_id,zelectron)
        ADIOS_READ(buf_id,phisave)
     endif
     call adios_close (buf_id)
     if(mype==0)write(*,*)"param end::",size(phisave,1),size(phisave,2)
     zion0(6,:)=zion0_read 
     deallocate(zion0_read)
     if(nhybrid>0)then
        zelectron0(6,:)=zion0_read
        deallocate(zelectron0_read)
     endif
     !write(*,*)"mype=",mype," mi=",mi
  endif  ! end of read
  if(mype==1)write(*,*)"READ FILE CLOSED"
#endif
end subroutine restart_io
