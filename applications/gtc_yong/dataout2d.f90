Subroutine dataout2d
  use global_parameters
  use field_array
  include 'netcdf.inc'

  integer i,j,k,jt,kk,mzbig,nfile,icount,ierror,istatus,dataid,ncid,dimid(2)
  real(singleprec) tdum,wt,wz,zdum,datatmp(0:mzeta,mgrid),eachdata(mpsi,mtdiag/ntoroidal),alldata(mpsi,mtdiag),delt(0:mpsi)      
  character(len=12) cdum
  character(len=4) ddum(3)

  write(ddum(1),'("dpot")')
  write(ddum(2),'("dtem")')
  write(ddum(3),'("heat")')

  do kk=1,1

     if (kk==1)then
        datatmp=phi/(gyroradius*gyroradius)
     elseif(kk==2)then
        do k=1,mzeta
           datatmp(k,:)=dtemper(:,k)
        enddo
     else
        do k=1,mzeta
           datatmp(k,:)=heatflux(:,k)
        enddo
     endif

! find the value at theta=0 by interpolating from fieldline coordinates to magnetic coordinates assuming mzeta=1
     delt=real(mtheta)
     mzbig=mtdiag/ntoroidal
     do k=1,mzbig
        wz=real(k)/real(mzbig)
        
! toroidal position (0,1]
        zdum=(real(mrank_toroidal)+wz)/real(ntoroidal)

        do i=1,mpsi

! poloidal location (0,1]
           tdum=10.0-zdum*qtinv(i)
           tdum=(tdum-aint(tdum))*delt(i)
           jt=max(0,min(mtheta(i)-1,int(tdum)))
           tdum=tdum-real(jt)
           jt=igrid(i)+jt

! linear interpolation along field line
           eachdata(i,k)=((1.0-wt)*datatmp(1,jt)+wt*datatmp(1,jt+1))*wz+&
                (1.0-wz)*((1.0-wt)*datatmp(0,jt)+wt*datatmp(0,jt+1))
        enddo
     enddo

! gather data to PE0
     icount=mpsi*mtdiag/ntoroidal
     call MPI_GATHER(eachdata,icount,MPI_REAL,alldata,icount,MPI_REAL,0,MPI_COMM_WORLD,ierror)

     if(mype==0)then
! write to data file in NetCDF format
        nfile=(mstepall+istep)/ndiag
        if(nfile<10)then
           write(cdum,'(a4,"000",i1,".ncd")')ddum(kk),nfile
        elseif(nfile<100)then
           write(cdum,'(a4,"00",i2,".ncd")')ddum(kk),nfile
        elseif(nfile<1000)then
           write(cdum,'(a4,"0",i3,".ncd")')ddum(kk),nfile
        else
           write(cdum,'(a4,i4,".ncd")')ddum(kk),nfile
        endif
        
! open netcdf data file
        istatus=nf_create(cdum,nf_clobber,ncid)
! define data array dimension
        istatus=nf_def_dim(ncid,'radial_grid',mpsi,dimid(1))
        istatus=nf_def_dim(ncid,'toroidal_grid',mtdiag,dimid(2))
! define data array id
        istatus=nf_def_var(ncid,ddum(kk),nf_real,2,dimid,dataid)
! end of netcdf definition
        istatus=nf_enddef(ncid)
! write data
        istatus=nf_put_var_real(ncid,dataid,alldata)
        istatus=nf_close(ncid)
     endif

  enddo

end subroutine dataout2d


