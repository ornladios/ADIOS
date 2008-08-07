Subroutine dataout3d
  use global_parameters
  use field_array
  use particle_decomp
  use particle_tracking
  include 'netcdf.inc'

  integer i,j,k,n,ij,istatus,ncid,dataid2(3),dimid(2),dataid1(2)
  real(singleprec) data3d(mgrid,mzeta),radial(mpsi+1)
  real(singleprec) r,theta,zeta,theta0
  character(len=1) cdum(3)
  character(len=4) tdum
  character(len=9) vdum
  character(len=20) fdum



  if(istep==ndiag)then
     if(mype==0)then
     !!!   write(*,*)'istep=',istep
        do i=0,mpsi
           radial(i+1)=(a0+deltar*real(i))/gyroradius
        enddo
        write(*,*)'write RUNdimen.ncd'
! write run-dimension
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
        istatus=nf_put_var_int(ncid,dataid1(1),numberpe)
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

! grid coordinates,each PE writes to a separate file 
     cdum(1)='X'
     cdum(2)='Y'
     cdum(3)='Z'
     do n=1,3
        do k=1,mzeta
           zeta=zetamin+deltaz*real(k)
           do i=0,mpsi
              r=a0+deltar*real(i)
              do j=0,mtheta(i)
                 ij=igrid(i)+j
                 theta=deltat(i)*real(j)+zeta*qtinv(i)
!                 theta0=theta+r*sin(theta)
                 theta0=theta
! grid coordinates (x,y,z), use geometry center as the origin 
                 if(n==1)then
                    data3d(ij,k)=cos(zeta)*(1.0+r*cos(theta0))
                 elseif(n==2)then
                    data3d(ij,k)=-sin(zeta)*(1.0+r*cos(theta0))
                 else
                    data3d(ij,k)=r*sin(theta0)
                 endif
              enddo
           enddo
        enddo
        
! coordinate data file name
        if(myrank_toroidal < 10)then
           write(fdum,'("NCD",a1,"coor.00",i1)')cdum(n),myrank_toroidal
        elseif(myrank_toroidal < 100)then
           write(fdum,'("NCD",a1,"coor.0",i2)')cdum(n),myrank_toroidal
        else
           write(fdum,'("NCD",a1,"coor.",i3)')cdum(n),myrank_toroidal
        endif

! variable name
        write(vdum,'(a1,"coordina")')cdum(n)

! open netcdf data file
        istatus=nf_create(fdum,nf_clobber,ncid)
! define data array dimension
        istatus=nf_def_dim(ncid,'poloidal_grid',mgrid,dimid(1))
        istatus=nf_def_dim(ncid,'toroidal_grid',mzeta,dimid(2))
! define data array id
        istatus=nf_def_var(ncid,vdum,nf_real,2,dimid,dataid1(1))
! end of netcdf definition
        istatus=nf_enddef(ncid)
! write data
        istatus=nf_put_var_real(ncid,dataid1(1),data3d)
        istatus=nf_close(ncid)        
     enddo
  endif

! potential data file
  if(myrank_partd==0)then
     do k=1,mzeta
        do i=0,mpsi
           do j=0,mtheta(i)
              ij=igrid(i)+j
              data3d(ij,k)=phi(k,ij)
           enddo
        enddo
     enddo
  
     write(fdum,'("PHI_",i0,"_",i0,".ncd")')(mstepall+istep),myrank_toroidal

   ! variable name
     write(vdum,'("Potential")')

   ! open netcdf data file
     istatus=nf_create(fdum,nf_clobber,ncid)
   ! define data array dimension
     istatus=nf_def_dim(ncid,'poloidal_grid',mgrid,dimid(1))
     istatus=nf_def_dim(ncid,'toroidal_grid',mzeta,dimid(2))
   ! define data array id
     istatus=nf_def_var(ncid,vdum,nf_real,2,dimid,dataid1(1))
   ! end of netcdf definition
     istatus=nf_enddef(ncid)
   ! write data
     istatus=nf_put_var_real(ncid,dataid1(1),data3d)
     istatus=nf_close(ncid)
  endif
        
end subroutine dataout3d

