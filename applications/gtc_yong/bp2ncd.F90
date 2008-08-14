module precision
  use mpi
  integer, parameter :: doubleprec=selected_real_kind(12),&
       singleprec=selected_real_kind(6),&
       defaultprec=kind(0.0)

#ifdef DOUBLE_PRECISION
  integer, parameter :: wp=doubleprec,mpi_Rsize=MPI_DOUBLE_PRECISION,&
                        mpi_Csize=MPI_DOUBLE_COMPLEX
#else
  integer, parameter :: wp=singleprec,mpi_Rsize=MPI_REAL,&
                        mpi_Csize=MPI_COMPLEX
#endif

end module precision

program bp2ncd
!!reading .bp 3d fluid data -> netcdf
  use precision
  use mpi  
  implicit none
  include 'netcdf.inc'
  integer :: istatus,ncid,dataid,dataid1(2),dataid2(3),dimid(2)
  integer :: mzeta,mpoloidal,mgrid,mgrid1,kp,mzetamax,mgrid_tor
  integer:: mz,m,n,inqdim(10),inqlen(10),i,jt,mpsi,ii,nstart,nend,ndel,ip,k
!!mpi parameters:
  integer :: mype,numberpe,tag
  integer :: ntoroidal,npartdom
  integer :: toroidal_comm,partd_comm
  integer :: nproc_toroidal,myrank_toroidal
  integer  :: toroidal_domain_location,particle_domain_location
  integer :: ierror
  character(len=5) cdum
  character(60) vname,cpdum,namef
  character(len=80) namein,dirpath
  integer,dimension(:),allocatable :: mtheta,igrid,itran,radial,eachdatatmp
  real(wp),dimension(:,:),allocatable :: eachdata,alldata
  integer :: ncdout
  logical file_exist
  character(len=4) tdum
  character(len=9) vdum
  character(len=20) fdum
  character(len=100):: input_fname,dirstr,nameout

  namelist /bp2ncd_in_params/ dirpath,nstart,nend,ndel
!ADIOS
  integer*8 group_handle,type_id
  integer next,previous,current
   
  #if ADIOS 
     #define ADIOS_READ_PATH(a,b,c) call adios_read_path(a,'b'//char(0),b,c//char(0))
     #define ADIOS_READ(a,b) call adios_read(a,'b'//char(0),b)
  #endif

! MPI initialize
  call mpi_init(ierror)
  call mpi_comm_size(mpi_comm_world,numberpe,ierror)
  call mpi_comm_rank(mpi_comm_world,mype,ierror)
  
tag=50  !!mpi tag: any integer number
!!XY open bp2ncd.in file for input parameters
! Test if the input file ncdpost.in exists
  if(mype==0)then
     ncdout=9
     open(ncdout,file='bp2ncd.out',status='replace')
     inquire(file='bp2ncd.in',exist=file_exist)
     if (file_exist) then
        open(55,file='bp2ncd.in',status='old')
        read(55,nml=bp2ncd_in_params)
        close(55)
        write(ncdout,nml=bp2ncd_in_params)

     else
        write(ncdout,*)'******************************************'
        write(ncdout,*)'*** NOTE!!! Cannot find file bp2ncd.in !!!'
        write(ncdout,*)'*** run exit...'
        write(ncdout,*)'******************************************'
     endif

  endif


 call mpi_bcast(nstart,1,mpi_integer,0,mpi_comm_world,ierror)
 call mpi_bcast(nend,1,mpi_integer,0,mpi_comm_world,ierror)
 call mpi_bcast(ndel,1,mpi_integer,0,mpi_comm_world,ierror)


#if ADIOS
   CALL adios_init ("config_bp2ncd.xml"//char(0), MPI_COMM_WORLD, MPI_COMM_SELF, MPI_INFO_NULL)
! open file RUNcoord.bp for grid dimensions and other parameters
   fdum='RUNcoord.bp'//char(0)
   call adios_get_group(type_id,"output3d.0"//char(0))
   call adios_open_read(group_handle,type_id,fdum)
!!read dimension
   ADIOS_READ(group_handle,mpsi)
   ADIOS_READ(group_handle,mgrid)
   ADIOS_READ(group_handle,mzeta)
   ADIOS_READ(group_handle,kp)
   ADIOS_READ(group_handle,myrank_toroidal)
   ADIOS_READ(group_handle,nproc_toroidal) 

     !!debug output
        if(mype==0)then
           open(127,file='readdiag.out',status='replace')
           write(127,*)'mpsi= ',mpsi
           write(127,*)'mzeta= ',mzeta
           write(127,*)'kp= ',kp
           write(127,*)'nproc_toroidal= ',nproc_toroidal
           write(127,*)'myrank_toroidal= ',myrank_toroidal
           write(127,*)'mgrid= ',mgrid
           close(127)
         endif

   if(allocated(radial))deallocate(radial)
   allocate(radial(1:mpsi+1))
   if(allocated(itran))deallocate(itran)
   allocate(itran(0:mpsi))
   if(allocated(mtheta))deallocate(mtheta)
   allocate(mtheta(0:mpsi))
   if(allocated(igrid))deallocate(igrid)
   allocate(igrid(0:mpsi))
!!initialization
   radial=0.0
   itran=0
   igrid=0
   mtheta=0

 
   ADIOS_READ(group_handle,radial)
   call adios_read(group_handle,"mtheta"//char(0),mtheta)
   mtheta=mtheta-1
   ADIOS_READ(group_handle,itran)
 
   mgrid1=0
   do i=1,mpsi 
      igrid(i)=igrid(i-1)+mtheta(i-1)+1
      mgrid1=mgrid1+mtheta(i-1)
   enddo
   mgrid1=mgrid1+mtheta(mpsi)
   mpoloidal=mgrid1
   
  if(nproc_toroidal .ne. numberpe)then
     write(*,*)'*** incorrect # of  processors, should be = ',nproc_toroidal
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
     return
  endif

!!write Rundimen.ncd: 8k dataset 
! write run-dimension
  if(mype==0)then
        write(ncdout,*)"read mgrid = ",mgrid
        write(ncdout,*)"calc mgrid = ",mgrid1
        fdum='RUNdimen.ncd'
! open netcdf data file
        istatus=nf_create(fdum,nf_clobber,ncid)

! define data array dimension
        istatus=nf_def_dim(ncid,'scalar',1,dimid(1))
        istatus=nf_def_dim(ncid,'flux-surfaces',mpsi+1,dimid(2))
! define data array id
        istatus=nf_def_var(ncid,'PE-number',nf_int,1,dimid(1),dataid1(1))
        istatus=nf_def_var(ncid,'flux-surface-number',nf_int,1,dimid(1),dataid1(2))
        istatus=nf_def_var(ncid,'radial-grids',nf_real,1,dimid(2),dataid2(1))
        istatus=nf_def_var(ncid,'poloidal-grids',nf_int,1,dimid(2),dataid2(2))
        istatus=nf_def_var(ncid,'index-shift',nf_int,1,dimid(2),dataid2(3))
! end of netcdf definition
        istatus=nf_enddef(ncid)

! write data
        istatus=nf_put_var_int(ncid,dataid1(1),nproc_toroidal)
        istatus=nf_put_var_int(ncid,dataid1(2),mpsi+1)
        istatus=nf_put_var_real(ncid,dataid2(1),radial)
        istatus=nf_put_var_int(ncid,dataid2(2),mtheta+1)
        istatus=nf_put_var_int(ncid,dataid2(3),itran)
! check error
        if (istatus .ne. nf_noerr) then
           print *, nf_strerror(istatus)
        endif
! close netcdf data file
        istatus=nf_close(ncid)
  endif


!!  if(mype.eq.0)then
!!     kp=1   !add the plane zeta=0 to mype=0 (zetamin=0 for mype=0)
!!  else
!!     kp=0   !no changes for mype > 0
!!  endif
  
  do n=nstart,nend,ndel
    
     write(input_fname,'("PHI_",i5.5,".bp")')n
     input_fname=trim(input_fname)//char(0)
     if(mype==0)then
        write(ncdout,*)input_fname," started"
     endif
     call adios_get_group (type_id, "output3d.1"//char(0))
     call adios_open_read (group_handle, type_id, input_fname)
     ADIOS_READ(group_handle,mzeta)
     ADIOS_READ(group_handle,kp)
     !!!call adios_read(group_handle,"mzeta+kp"//char(0),mzeta)
     !!!mzeta=mzeta-kp
     mzetamax=mzeta*numberpe+1
     mgrid_tor=mpoloidal*mzeta

     if(n==nstart .and. mype==0)then
        write(ncdout,*)'mzeta= ',mzeta, 'mzetamax= ',mzetamax
     endif

     if(allocated(eachdata))deallocate(eachdata)
     allocate (eachdata(mpoloidal,mzeta+kp))
     eachdata=0.0
     if(allocated(eachdatatmp))deallocate(eachdatatmp)
     allocate (eachdatatmp(mgrid_tor))
     eachdatatmp=0.0
     if(allocated(alldata))deallocate(alldata)
     allocate (alldata(mpoloidal,mzetamax))
     alldata=0.0
     call adios_read(group_handle,"phi"//char(0),eachdata)
     if(mype==0)then
        do k=1,mzeta+kp
           alldata(:,k)=eachdata(:,k)
        enddo
     else
	 eachdatatmp=reshape(eachdata,(/mgrid_tor/))
     endif
	 
     do ip=1,numberpe-1
        if(mype==0)then
           call mpi_recv(eachdatatmp,mgrid_tor,mpi_Rsize,ip,tag,mpi_comm_world,istatus,ierror)
           do k=1,mzeta
              alldata(1:mpoloidal,ip*mzeta+1+k)=eachdatatmp(mpoloidal*(k-1)+1:mpoloidal*k)
           enddo
        else
           call mpi_send(eachdatatmp,mgrid_tor,mpi_Rsize,0,tag,mpi_comm_world,istatus,ierror)
        endif
		call MPI_BARRIER(MPI_COMM_WORLD,ierror)
     enddo
     
     !!write alldata to ncdcdf file with
     if(mype==0)then
        write(ncdout,*)input_fname, " reading in"
        write(nameout,'("PHI_",i0,".ncd")')n
        vname='Potential'
        vname=trim(vname)
        istatus=nf_create(trim(nameout),nf_clobber,ncid)
        istatus=nf_def_dim(ncid,'poloidal_grids',mpoloidal,dimid(1))
        istatus=nf_def_dim(ncid,'toroidal_grids',mzetamax,dimid(2))
        istatus=nf_def_var(ncid,vname,nf_real,2,dimid,dataid)
        istatus=nf_enddef(ncid)
        istatus=nf_put_var_real(ncid,dataid,alldata)

     ! check error
        if (istatus .ne. nf_noerr) then
           print *, nf_strerror(istatus)
        endif
        istatus=nf_close(ncid)
     
        write(ncdout,*)n,nameout," done"

     endif

  enddo


#else

    if(mype==0)then
       write(ncdout,*)"no adios flags and then no adios compiled in the code"
    end if

#endif

  if(mype==0)then
     close(ncdout)
  endif

  CALL adios_finalize (mype)


! MPI finalize
  call mpi_finalize(ierror)


end program bp2ncd
