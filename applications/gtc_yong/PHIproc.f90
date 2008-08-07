program PHIproc

implicit none
include 'netcdf.inc'

integer :: istatus,ncid,dataid,dimid(2),rundimen,mzeta,mpoloidal,,numberpe,&
       mz,m,n,inqdim(10),inqlen(10),i,jt,mpsi,ii,nstart,nend,ndel

character(len=5)cdum
character(len=40):: vname,nameout,cpdum,namef
character(len=80):: namein,dirpath

integer,dimension(:,:),allocatable::mtheta,igrid,itran
integer::procout
logical :: file_exist

namelist/dataproc_in_params/dirpath,nstart,nend,ndel

procout=9
open(procout,file='dataproc.out',status='replace')

!test if exist input file

inquire(file='PHIproc.in',exist=file_exist)

if(file_exist)then
   open(55,file='PHIproc.in',status='old')
   read(55,nml=dataproc_in_params)
   close(55)
   write(procout,nml=dataproc_in_params)
else
   write(procout,*)'******************************'
   write(procout,*)'*** NOTE!!! Cannot find file PHIproc.in !!!'
   write(procout,*)'*** Please set up input deck...'
   write(procout,*)'******************************************'
   stop
   return
endif

!reading data from 'RUNdimen'
write(procout,*)'...opening', namein
istatus=nf_open(namein,nf_nowrite,ncid)

istatus=nf_inq_varid(ncid,'PE-number',dataid)
istatus=nf_get_var_int(ncid,dataid,numberpe)

istatus=nf_inq_varid(ncid,'flux-surface-number',dataid)
istatus=nf_get_var_int(ncid,dataid,mpsi)
mpsi=mpsi-1
write(procout,*)'mpsi=', mpsi
write(procout,*)'numberpe=', numberpe

if(allocated(igrid))deallocate(igrid)
if(allocated(mtheta))deallocate(mtheta)
if(allocated(itran))deallocate(itran)

allocate(mtheta(0:mpsi),igrid(0:mpsi),itran(0:mpsi)

istatus=nf_inq_varid(ncid,'poloidal-grids',dataid)
istatus=nf_get_var_int(ncid,dataid,mtheta)
mtheta=mtheta-1

istatus=nf_inq_varid(ncid,'index-shift',dataid)
istatus=nf_get_var_int(ncid,dataid,itran)

!close netcdf data file
istatus=nf_close(ncid)

print,'mtheta=',mtheta


end program PHIproc
