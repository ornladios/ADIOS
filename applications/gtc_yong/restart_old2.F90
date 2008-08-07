module mpiio_params
  use precision
  integer mfhandle,merror, mstatus(MPI_STATUS_SIZE),mfmode, mrequest,mpi_offset_type
  INTEGER(KIND=MPI_OFFSET_KIND) mype_filesize, sum_filesize
  DOUBLE PRECISION start_time,end_time
  character(len=32) restart_fname
  integer:: rank,dims(4)
  integer:: dtype,group_id(1)
  character(len=100)aname,dirstr,aval
end module mpiio_params

subroutine restart_io(iop)
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use particle_decomp
  use mpiio_params
  use data_type
  use particle_tracking
  implicit none

  integer i,j,k,subsize,startidx,endidx
  integer mquantity,mflx,n_mode,mstepfinal,noutputs
  real(wp) dum
  character(*),intent(in)::iop

#if ADIOS
#define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
#define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0))
#define ADIOS_READ(a,b) call adios_read(a,'b'//char(0),b)
  integer*8 file_handle, io_type
  real(wp),dimension(:),allocatable::zion0_read
#endif

  if(iop/="read" .and. iop/="write")then
     write(0,*)'*** subroutine restart_io (iop <> "read" or "write")',iop
     call MPI_ABORT(MPI_COMM_WORLD,1,merror)
     return
  endif
  mpi_offset_type=MPI_INTEGER4 
  !if(mype==1)write(stdout,*)"MFLUX,MPSI:",mflux,mpsi
  write(restart_fname,'(a,i5.5,"_",i5.5,".bp")')"restart_",myrank_toroidal,mstepall+istep
  #ifdef __TIMER
  call MPI_BARRIER(MPI_COMM_WORLD,merror)
  start_time=MPI_WTIME()
  #endif

#if ADIOS
     ! setup the element path for this node
     write(dirstr,'("/node",i5.5,"/param")')mype
     dirstr=trim(dirstr)//char(0)
     call MPI_BARRIER(MPI_COMM_WORLD,merror)
     start_time = MPI_WTIME()
     call adios_get_group (io_type, "restart"//char(0))
     ! set the path for all vars in the type for proper sizing
     call adios_set_path (io_type,dirstr//char(0));
     restart_fname=trim(restart_fname)//char(0)
     if (iop=="read") then
        call adios_open_read (file_handle, io_type, restart_fname)
     else
        call adios_open (file_handle, io_type, restart_fname)
     endif
     ! write the sizing paramters for both reading and writing
     ADIOS_WRITE(file_handle,mflux)
     ADIOS_WRITE(file_handle,mpsi+1)
     ADIOS_WRITE(file_handle,mzeta+1)
     ADIOS_WRITE(file_handle,nparam)
     ADIOS_WRITE(file_handle,mimax)
     ADIOS_WRITE(file_handle,mgrid)
     ADIOS_WRITE(file_handle,partd_comm)
  if(iop=="write")then
     ADIOS_WRITE(file_handle,mzeta)
     ADIOS_WRITE(file_handle,mi)
     ADIOS_WRITE(file_handle,me)
     ADIOS_WRITE(file_handle,ntracer)
     ADIOS_WRITE(file_handle,etracer)
     ADIOS_WRITE(file_handle,rdtemi)
     ADIOS_WRITE(file_handle,rdteme)
     ADIOS_WRITE(file_handle,ptracer)
     ADIOS_WRITE(file_handle,pfluxpsi)
     ADIOS_WRITE(file_handle,phi00)
     ADIOS_WRITE(file_handle,phip00)
     ADIOS_WRITE(file_handle,zonali)
     ADIOS_WRITE(file_handle,zonale)
     call adios_write(file_handle,"zion0",zion0(7,:))
     ADIOS_WRITE(file_handle,phi)
     ADIOS_WRITE(file_handle,zion)
     call adios_get_data_size (file_handle, mype_filesize)
     call adios_close (file_handle)
   #ifdef __TIMER
     call MPI_BARRIER(MPI_COMM_WORLD,merror)
     end_time=MPI_WTIME()
     call MPI_REDUCE(mype_filesize,sum_filesize,1,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,merror)
     if(mype==0)then
        write(stdout,'("NtoM Time(s), MBytes, MB/s ",a)')iop
        write(stdout,*)sum_filesize/1024/1024,mype_filesize/1024/1024
        write(stdout,*)end_time-start_time,mype_filesize*numberpe/1024/1024,mype_filesize*numberpe/((end_time-start_time)*1024*1024)
     endif
   #endif

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
     if(mype==0)write(*,*)"param::",mpsi,mzeta
     allocate(zion0_read(mimax))
     ADIOS_READ(file_handle,mzeta)
     ADIOS_READ(file_handle,mi)
     ADIOS_READ(file_handle,me)
     ADIOS_READ(file_handle,ntracer)
     ADIOS_READ(file_handle,etracer)
     ADIOS_READ(file_handle,rdtemi)
     ADIOS_READ(file_handle,rdteme)
     ADIOS_READ(file_handle,ptracer)
     ADIOS_READ(file_handle,pfluxpsi)
     ADIOS_READ(file_handle,phi00)
     ADIOS_READ(file_handle,phip00)
     ADIOS_READ(file_handle,zonali)
     ADIOS_READ(file_handle,zonale)
     call adios_read(file_handle,"zion0"//char(0),zion0_read)
     zion0(7,:)=zion0_read 
     ADIOS_READ(file_handle,phi)
     ADIOS_READ(file_handle,zion)
     call adios_get_data_size (file_handle, mype_filesize)
     deallocate(zion0_read)
     call adios_close (file_handle)
     if(myrank_partd<nproc_partd-1)call MPI_Wait(mrequest,mstatus,merror)
  endif  ! end of read
#endif
     if(mype==0)write(stdout,*)"FILE CLOSED"
end subroutine restart_io
