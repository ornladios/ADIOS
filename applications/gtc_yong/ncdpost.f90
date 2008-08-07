program ncdpost
  implicit none
  include 'netcdf.inc'
  integer istatus,ncid,dataid,dimid(2),rundimen,mzeta,mpoloidal,numberpe,&
       mz,m,n,inqdim(10),inqlen(10),i,jt,mpsi,ii,nstart,nend,ndel
  character(len=5) cdum
  character(40) vname,nameout,cpdum,namef
  character(len=80) namein,dirpath
  integer,dimension(:),allocatable :: mtheta,igrid,itran
  real,dimension(:,:),allocatable :: eachdata,alldata
  integer :: ncdout
  logical file_exist
  namelist /ncd_in_params/ dirpath,nstart,nend,ndel

 
! open file for grid dimensions and other parameters
  ncdout=9
  open(ncdout,file='ncdout.out',status='replace')

!!XY open ncdpost.in file
! Test if the input file ncdpost.in exists
    inquire(file='ncdpost.in',exist=file_exist)
    if (file_exist) then
       open(55,file='ncdpost.in',status='old')
       read(55,nml=ncd_in_params)
       close(55)
       write(ncdout,nml=ncd_in_params)

    else
       write(ncdout,*)'******************************************'
       write(ncdout,*)'*** NOTE!!! Cannot find file ncdpost.in !!!'
       write(ncdout,*)'*** Using default run parameters...'
       write(ncdout,*)'******************************************'
    endif



  namein='RUNdimen.ncd'

  write(ncdout,*)'... opening ',namein
  istatus=nf_open(namein,nf_nowrite,ncid)

! read data
  istatus=nf_inq_varid(ncid,'PE-number',dataid)
  istatus=nf_get_var_int(ncid,dataid,numberpe)
!  write(0,*)' numberpe =',numberpe
!  write(*,'(" Number of processors involved in I/O = ",$)')
!  read(*,*)numberpe
  numberpe=32

  istatus=nf_inq_varid(ncid,'flux-surface-number',dataid)
  istatus=nf_get_var_int(ncid,dataid,mpsi)
  mpsi=mpsi-1
  write(0,*)' mpsi =',mpsi
  write(ncdout,*)' mpsi =',mpsi
  allocate(mtheta(0:mpsi),igrid(0:mpsi),itran(0:mpsi))

  istatus=nf_inq_varid(ncid,'poloidal-grids',dataid)
  istatus=nf_get_var_int(ncid,dataid,mtheta)
  mtheta=mtheta-1

  istatus=nf_inq_varid(ncid,'index-shift',dataid)
  istatus=nf_get_var_int(ncid,dataid,itran)

! close netcdf data file
  istatus=nf_close(ncid)

! starting point for a poloidal grid
  igrid(0)=1
  do i=1,mpsi
     igrid(i)=igrid(i-1)+mtheta(i-1)+1
  enddo

! Directory where to find the files
!  write(*,'(" Where are the files? :",$)')
!  read(*,'(a80)')dirpath
!!********************
!  dirpath='/tmp/work/hardes/tem82'
!!*************************
 
  write(ncdout,*)trim(dirpath)



! open data file
!!*************************
!nstart=5
!nend=2500
!ndel=5
!!**************************


  write(ncdout,*)'start_time, end_time, time_interval'
!  read(*,*)nstart,nend,ndel
  write(ncdout,*)nstart,nend,ndel
!  nstart=nstart-3

  do n=nstart,nend,ndel

     write(ncdout,*)n," started"
     
! coordinates data
     if(n<=0 .and. n==nstart)then
        write(cdum,'("Xcoor")')
        vname='Xcoordina'
     elseif(n<=0 .and. n==nstart+ndel)then
        write(cdum,'("Ycoor")')
        vname='Ycoordina'
     elseif(n<=0 .and. n==nstart+2*ndel)then
        write(cdum,'("Zcoor")')
        vname='Zcoordina'
! potential data
     else
        vname='Potential'
     endif

! read each PE file
     do m=0,numberpe-1
        if(n <= 0)then
           if(m < 10)then
              write(namef,'("NCD",a5,".00",i1)')cdum,m
           elseif(m < 100)then
              write(namef,'("NCD",a5,".0",i2)')cdum,m
           else
              write(namef,'("NCD",a5,".",i3)')cdum,m
           endif
        else
           write(namef,'("PHI_",i0,"_",i0,".ncd")')n,m
        endif
        write(ncdout,*)'namef = ',trim(namef)
        namein=trim(dirpath)//'/'//trim(namef)
        write(ncdout,*)'   ... namein =',trim(namein)

! open NCD data file
        istatus=nf_open(namein,nf_nowrite,ncid)

        if(m==0 .and. n==nstart)then
! real data dimension
           if(n <= 0)then
              istatus=nf_inq_varid(ncid,'Xcoordina',dataid)
           else
              istatus=nf_inq_varid(ncid,'Potential',dataid)
           endif
           istatus=nf_inq_vardimid(ncid,dataid,inqdim)
           istatus=nf_inq_dimlen(ncid,inqdim(1),inqlen(1))
           istatus=nf_inq_dimlen(ncid,inqdim(2),inqlen(2))
           mpoloidal=inqlen(1)
           mz=inqlen(2)
! allocate memory
           mzeta=mz*numberpe+1
           write(0,*)'allocating alldata = ',mpoloidal,mzeta,mz
           write(ncdout,*)'allocating alldata = ',mpoloidal,mzeta,mz
           allocate (eachdata(mpoloidal,mz),alldata(mpoloidal,mzeta))
        endif

        istatus=nf_inq_varid(ncid,vname,dataid)
        istatus=nf_get_var_real(ncid,dataid,eachdata)
        istatus=nf_close(ncid)

! check error
        if (istatus .ne. nf_noerr) then
           print *, nf_strerror(istatus)
        endif

        alldata(:,mz*m+2:mz*(m+1)+1)=eachdata
     enddo

! toroidal BC
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        alldata(ii+1:ii+jt,1)=cshift(alldata(ii+1:ii+jt,mzeta),-itran(i))
        alldata(ii,1)=alldata(ii+jt,1)
     enddo

! write data to a single file
     if(n <= 0)then
        write(nameout,'(a5,".ncd")')cdum
     else
        write(nameout,'("PHI_",i0,".ncd")')n
     endif
     istatus=nf_create(nameout,nf_clobber,ncid)
     istatus=nf_def_dim(ncid,'poloidal_grids',mpoloidal,dimid(1))
     istatus=nf_def_dim(ncid,'toroidal_grids',mzeta,dimid(2))
     istatus=nf_def_var(ncid,vname,nf_real,2,dimid,dataid)
     istatus=nf_enddef(ncid)
     istatus=nf_put_var_real(ncid,dataid,alldata)

! check error
        if (istatus .ne. nf_noerr) then
           print *, nf_strerror(istatus)
        endif
     istatus=nf_close(ncid)
     
     write(ncdout,*)n,nameout," done"

  enddo

  close(ncdout)
end program ncdpost
