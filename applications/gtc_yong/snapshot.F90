subroutine snapshot
  use global_parameters
  use particle_array
  use particle_decomp
  use field_array
  use diagnosis_array
  implicit none
  
  integer,parameter :: mbin_psi=5, mbin_u=41
  integer i,ibin,ip,iu,j,jj,jm,jt,k,kz,m,nsnap,ierror,icount
  real(wp) eright,uright,uleft,pdum,tdum,b,upara,energy,pitch,q,dela,delu,&
       qdum,r,vthi_inv,delr,delz,phisn(mzeta,mtheta(mpsi/2)),&
       eachpsi(mtheta(mpsi/2)*mzetamax),vthe_inv,aion_inv,ae_inv,wt,&
       delt(0:mpsi),weight,adum(0:mpsi),bdum(mbin_u,mbin_psi)
  real(wp) ubin_max,pbin_max,ebin_max,xtmp,ytmp,ztmp
  real(wp),dimension(0:mpsi) :: marker,fflows,dflows,ftem,dtem
  real(wp),dimension(mbin_u,mbin_psi) :: ubin,dubin,pbin,dpbin,ebin,debin
!  real(wp),dimension(mzeta,0:mpsi,0:mthetamax) :: ftotal,tpara,tperp
  real(wp),external :: boozer2x,boozer2z,r2psi
  real(wp),dimension(:,:,:),allocatable::booz_out,flux3darray
  real(wp),dimension(:,:),allocatable::dataout_u,dataout_r
  character(len=11) date,time
  character(len=13) cdum

#ifdef ADIOS
  character(len=30) bpdum
  integer*8 filehandle, groupfilehandle
  integer*8 :: group_proffilehandle
  integer :: adios_err
  integer ::total_size 
!  #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b,adios_err)
#endif

! write particle data for re-run
!!  call restart_write
  call restart_io("write")
! number of poloidal grid
  jm=mtheta(mpsi/2)

! zion space grids: parallel velocity, kinetic energy, pitch angle
  uright=umax
  uleft=-uright  
  dela=real(mbin_psi)/(a1-a0)
  delu=1.0/(uright-uleft)
  eright=1.0/(umax*umax/4.0)
  vthi_inv=aion/(gyroradius*abs(qion))
  vthe_inv=vthi_inv*sqrt(tite*aelectron/aion)
  aion_inv=1.0/aion
  ae_inv=1.0/aelectron
  delr=1.0/deltar
  delt=1.0/deltat
  delz=1.0/deltaz

! radial and 3-d volume data
  ubin=0.0
  pbin=0.0
  ebin=0.0
  dubin=0.0
  dpbin=0.0
  debin=0.0
  marker=0.0
  fflows=0.0
  dflows=0.0
  ftem=0.0
  dtem=0.0
!  ftotal=0.0
!  tpara=0.0
!  tperp=0.0

  if(nhybrid==0)then

     do m=1,mi 
        weight=zion0(6,m)

        r=sqrt(2.0*zion(1,m))
        ibin=max(1,min(mbin_psi,1+int((r-a0)*dela)))
        ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
        jt=max(0,min(mtheta(ip),int(zion(2,m)*delt(ip)+0.5)))
        kz=max(1,min(mzeta,1+int((zion(3,m)-zetamin)*delz)))
     
        b=1.0/(1.0+r*cos(zion(2,m)))
        upara=zion(4,m)*b*qion*aion_inv*vthi_inv !normalized by v_thi
        energy=max(1.0e-20,0.5*upara*upara+&
             zion(6,m)*zion(6,m)*b*aion_inv*vthi_inv*vthi_inv) !normalized by T_i
        pitch=upara/sqrt(2.0*energy)
     
        !iu=max(1,min(mbin_u,1+int(real(mbin_u)*(upara-uleft)*delu)))
        iu=1+int(real(mbin_u-1)*(upara-uleft)*delu)
        if(iu >= 1 .and. iu <= mbin_u)then
           ubin(iu,ibin)=ubin(iu,ibin)+weight
           dubin(iu,ibin)=dubin(iu,ibin)+zion(5,m)
        endif
     
        !iu=max(1,min(mbin_u,1+int(real(mbin_u)*energy*eright)))
        iu=1+int(real(mbin_u-1)*energy*eright)
        if(iu >= 1 .and. iu <= mbin_u)then
           ebin(iu,ibin)=ebin(iu,ibin)+weight
           debin(iu,ibin)=debin(iu,ibin)+zion(5,m)
        endif
        
        !iu=max(1,min(mbin_u,1+int(real(mbin_u)*0.5*(pitch+1.0))))
        iu=1+int(real(mbin_u-1)*0.5*(pitch+1.0))
        if(iu >= 1 .and. iu <= mbin_u)then
           pbin(iu,ibin)=pbin(iu,ibin)+weight
           dpbin(iu,ibin)=dpbin(iu,ibin)+zion(5,m)
        endif
        
! radial profile of marker, v_para, and temperature
        marker(ip)=marker(ip)+weight
        fflows(ip)=fflows(ip)+upara*weight
        dflows(ip)=dflows(ip)+zion(5,m)*upara
        ftem(ip)=ftem(ip)+energy*weight
        dtem(ip)=dtem(ip)+zion(5,m)*energy
     
! parallel and perpendicular pressure perturbation
!     ftotal(kz,ip,jt)=ftotal(kz,ip,jt)+weight
!     tpara(kz,ip,jt)=tpara(kz,ip,jt)+zion(5,m)*0.5*upara*upara
!     tperp(kz,ip,jt)=tperp(kz,ip,jt)+zion(5,m)*energy
     enddo
  
  else

     do m=1,me 
        weight=zelectron0(6,m)

        r=sqrt(2.0*zelectron(1,m))
        ibin=max(1,min(mbin_psi,1+int((r-a0)*dela)))
        ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
        jt=max(0,min(mtheta(ip),int(zelectron(2,m)*delt(ip)+0.5)))
        kz=max(1,min(mzeta,1+int((zelectron(3,m)-zetamin)*delz)))
     
        b=1.0/(1.0+r*cos(zelectron(2,m)))
        upara=zelectron(4,m)*b*qelectron*ae_inv*vthe_inv
        energy=max(1.0e-20,0.5*upara*upara+&
             zelectron(6,m)*zelectron(6,m)*b*ae_inv*vthe_inv*vthe_inv)
        pitch=upara/sqrt(2.0*energy)
     
        iu=max(1,min(mbin_u,1+int(real(mbin_u)*(upara-uleft)*delu)))
        ubin(iu,ibin)=ubin(iu,ibin)+weight
        dubin(iu,ibin)=dubin(iu,ibin)+zelectron(5,m)
     
        iu=max(1,min(mbin_u,1+int(real(mbin_u)*energy*eright)))
        ebin(iu,ibin)=ebin(iu,ibin)+weight
        debin(iu,ibin)=debin(iu,ibin)+zelectron(5,m)
        
        iu=max(1,min(mbin_u,1+int(real(mbin_u)*0.5*(pitch+1.0))))
        pbin(iu,ibin)=pbin(iu,ibin)+weight
        dpbin(iu,ibin)=dpbin(iu,ibin)+zelectron(5,m)
        
! radial profile of marker, v_para, and temperature
        marker(ip)=marker(ip)+weight
        fflows(ip)=fflows(ip)+upara*weight
        dflows(ip)=dflows(ip)+zelectron(5,m)*upara
        ftem(ip)=ftem(ip)+energy*weight
        dtem(ip)=dtem(ip)+zelectron(5,m)*energy
    
! parallel and perpendicular pressure perturbation
!     ftotal(kz,ip,jt)=ftotal(kz,ip,jt)+weight
!     tpara(kz,ip,jt)=tpara(kz,ip,jt)+zelectron(5,m)*0.5*upara*upara
!     tperp(kz,ip,jt)=tperp(kz,ip,jt)+zelectron(5,m)*energy
     enddo

  endif

!  do i=0,mpsi
!     tpara(:,i,mtheta(i))=tpara(:,i,mtheta(i))+tpara(:,i,0)
!     tperp(:,i,mtheta(i))=tperp(:,i,mtheta(i))+tperp(:,i,0)
!     ftotal(:,i,mtheta(i))=ftotal(:,i,mtheta(i))+ftotal(:,i,0)
!  enddo

! do i=0,mpsi
!    ftotal(:,i,1:mtheta(i))=max(1.0e-6,ftotal(:,i,1:mtheta(i)))
!     tpara(:,i,1:mtheta(i))=tpara(:,i,1:mtheta(i))/ftotal(:,i,1:mtheta(i))
!     tperp(:,i,1:mtheta(i))=tperp(:,i,1:mtheta(i))/ftotal(:,i,1:mtheta(i))
!     ftotal(:,i,1:mtheta(i))=ftotal(:,i,1:mtheta(i))*markeri(:,i,1:mtheta(i))
!     ftotal(:,i,0)=ftotal(:,i,mtheta(i))
!  enddo

  icount=mbin_psi*mbin_u
  call MPI_REDUCE(ubin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)ubin=bdum
  call MPI_REDUCE(dubin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)dubin=bdum
  call MPI_REDUCE(pbin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)pbin=bdum
  call MPI_REDUCE(dpbin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)dpbin=bdum
  call MPI_REDUCE(ebin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)ebin=bdum
  call MPI_REDUCE(debin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)debin=bdum
  call MPI_REDUCE(marker,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)marker=adum
  call MPI_REDUCE(fflows,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)fflows=adum
  call MPI_REDUCE(dflows,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)dflows=adum
  call MPI_REDUCE(ftem,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)ftem=adum
  call MPI_REDUCE(dtem,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0)dtem=adum

  where(marker .lt. 1.0)marker=1.0
  fflows=fflows/marker
  dflows=dflows/marker
  ftem=ftem/marker
  dtem=dtem/marker
  marker(1:mpsi)=marker(1:mpsi)*pmarki(1:mpsi)

! gather flux surface profile
  eachpsi=0.0
  do j=1,jm
     phisn(:,j)=phi(1:mzeta,igrid(mpsi/2)+j)
  enddo
  call MPI_GATHER(phisn,jm*mzeta,mpi_Rsize,eachpsi,&
       jm*mzeta,mpi_Rsize,0,toroidal_comm,ierror)
  
  if(mype == 0)then
! normalization
     !dubin=dubin/max(1.0,ubin)
     !dpbin=dpbin/max(1.0,pbin)
     !debin=debin/max(1.0,ebin)

!!******
     allocate(booz_out(mpsi,0:jm,3))
     allocate(dataout_r(1:mpsi,9))
     allocate(dataout_u(1:mbin_u,mbin_psi*6+3))
     allocate(flux3darray(mzeta,ntoroidal,jm))

     ubin_max=maxval(ubin)
     pbin_max=maxval(pbin)
     ebin_max=maxval(ebin)
     dubin=dubin/ubin_max
     dpbin=dpbin/pbin_max
     debin=debin/ebin_max
         
! open snapshot output file
     !nsnap=int(real(mstepall+istep)*tstep/1000.0)
     nsnap=mstepall+istep
     write(cdum,'("snap",i5.5,".out")')nsnap
     open(snapout,file=cdum,status='replace')
         
! time step and grid dimension,write every other point in poloidal direction
     write(snapout,102)real(mstepall+istep)*tstep,q0+0.5*q1+0.25*q2
     write(snapout,101)mbin_u,mbin_psi,mpsi,jm,mzetamax,mpsi/2
     
! write velocity space grids: v_para, pitch, energy to data file
! number of quantities
     i=9
     write(snapout,101)i
     do i=1,mbin_u
        xtmp=uleft+(uright-uleft)*real(i-1)/real(mbin_u-1)
        dataout_u(i,1)=xtmp
        write(snapout,102)xtmp
     enddo
     do i=1,mbin_u
        xtmp=-1.0+2.0*real(i-1)/real(mbin_u-1)
        dataout_u(i,2)=xtmp
        write(snapout,102)xtmp
     enddo
     do i=1,mbin_u
        xtmp=1.0/eright*real(i-1)/real(mbin_u-1)
        dataout_u(i,3)=xtmp
        write(snapout,102)xtmp
     enddo


! velocity space distribution function
     do j=1,mbin_psi
        write(snapout,102)(ubin(i,j),i=1,mbin_u)
        write(snapout,102)(dubin(i,j),i=1,mbin_u)
        write(snapout,102)(pbin(i,j),i=1,mbin_u)
        write(snapout,102)(dpbin(i,j),i=1,mbin_u)
        write(snapout,102)(ebin(i,j),i=1,mbin_u)
        write(snapout,102)(debin(i,j),i=1,mbin_u)

        dataout_u(1:mbin_u,(j-1)*6+4)=ubin(1:mbin_u,j)
        dataout_u(1:mbin_u,(j-1)*6+5)=dubin(1:mbin_u,j)
        dataout_u(1:mbin_u,(j-1)*6+6)=pbin(1:mbin_u,j)
        dataout_u(1:mbin_u,(j-1)*6+7)=dpbin(1:mbin_u,j)
        dataout_u(1:mbin_u,(j-1)*6+8)=ebin(1:mbin_u,j)
        dataout_u(1:mbin_u,(j-1)*6+9)=debin(1:mbin_u,j)

     enddo

! radial grid
     i=9
     write(snapout,101)i
     do i=1,mpsi
        write(snapout,102)(a0+deltar*real(i))/gyroradius
     enddo
     write(snapout,102)(zonali(i),i=1,mpsi)
     write(snapout,102)(zonale(i),i=1,mpsi)
     write(snapout,102)(phip00(i)/gyroradius,i=1,mpsi)
!     write(snapout,102)(phi00(i)/gyroradius**2,i=1,mpsi)
     write(snapout,102)(marker(i),i=1,mpsi)
     write(snapout,102)(fflows(i),i=1,mpsi)
     write(snapout,102)(dflows(i),i=1,mpsi)
     write(snapout,102)(ftem(i),i=1,mpsi)
     write(snapout,102)(dtem(i),i=1,mpsi)

     dataout_r(1:mpsi,2)=zonali(1:mpsi)
     dataout_r(1:mpsi,3)=zonale(1:mpsi)
     dataout_r(1:mpsi,4)=phip00(1:mpsi)
     dataout_r(1:mpsi,5)=marker(1:mpsi)
     dataout_r(1:mpsi,6)=fflows(1:mpsi)
     dataout_r(1:mpsi,7)=dflows(1:mpsi)
     dataout_r(1:mpsi,8)=ftem(1:mpsi)
     dataout_r(1:mpsi,9)=dtem(1:mpsi)

! poloidal cross section: marker,density,potential,parallel flow and temperature
     i=3
     write(snapout,101)i
     do j=0,jm
        jj=j
        if(j>jm/2)jj=j-jm
        do i=1,mpsi
           pdum=r2psi(a0+deltar*real(i))
           tdum=pi+2.0*pi*real(jj)/real(jm)
           xtmp=boozer2x(pdum,tdum)
           booz_out(i,j,1)=xtmp
           write(snapout,102)xtmp          
        enddo
     enddo
     do j=0,jm
        jj=j
        if(j>jm/2)jj=j-jm
        do i=1,mpsi
           pdum=r2psi(a0+deltar*real(i))
           tdum=pi+2.0*pi*real(jj)/real(jm)
           xtmp=boozer2z(pdum,tdum)
           booz_out(i,j,2)=xtmp
           write(snapout,102)xtmp
        enddo
     enddo

     do j=0,jm
        jj=j
        if(j>jm/2)jj=j-jm
        do i=1,mpsi
           tdum=pi+2.0*pi*real(jj)/real(jm)
           jt=max(1,min(mtheta(i),1+int(tdum/deltat(i))))
           wt=tdum/deltat(i)-real(jt-1)
           xtmp=(wt*phi(0,igrid(i)+jt)+(1.0-wt)*phi(0,igrid(i)+jt-1))/&
                gyroradius**2
           booz_out(i,j,3)=xtmp
           write(snapout,102)xtmp
!!           write(snapout,102)wt*densityi(0,igrid(i)+jt)+(1.0-wt)*densityi(0,igrid(i)+jt-1)
        enddo
     enddo

! flux surface profile for marker, density, potential, parallel flow and temperature
     i=1
     write(snapout,101)i
     do j=1,jm
        do kz=1,ntoroidal
           do k=1,mzeta
              xtmp=eachpsi((kz-1)*mzeta*jm+(j-1)*mzeta+k)/gyroradius**2
              flux3darray(k,kz,j)=xtmp
              write(snapout,102)xtmp
           enddo
        enddo
     enddo

! eigenmode structure
     write(snapout,101)num_mode,m_poloidal,nmode
     do i=1,mpsi
        do kz=1,num_mode
           do j=1,m_poloidal
              write(snapout,102)eigenmode(j,kz,i)
           enddo
        enddo
     enddo

! close snapshot file
     close(snapout)

!!write snapshot.adio
#ifdef ADIOS
     write(bpdum,'("snap_",i5.5,".bp")')nsnap
     bpdum=trim(bpdum)//char(0)

!!! modified by zf2 for timing
!     CALL open_start_for_group(group_proffilehandle, "snapshot"//char(0),istep)

     !!!call adios_get_group (groupfilehandle, "snapshot"//char(0))
     call adios_open (filehandle, "snapshot"//char(0), bpdum, "w"//char(0), adios_err)

!     CALL open_end_for_group(group_proffilehandle,istep)

!     CALL write_start_for_group(group_proffilehandle,istep)

call adios_group_size(filehandle, 4 + &
4 + &
4 + &
4 + &
4 + &
4 + &
4 + &
4 + &
4 + &
4 + &
4*(mbin_u)*(3+mbin_psi*6) + &
4*(mpsi)*(9) + &
4*(mpsi)*(jm+1)*(3) + &
4 + &
4 + &
4*(mzeta)*(ntoroidal)*(jm) + &
4 + &
4*(num_mode) + &
4 + &
4*(mpsi)*(num_mode)*(m_poloidal), &
total_size, &
0, &
adios_err &
)


call adios_write (filehandle, "timestep"//char(0), (mstepall+istep),adios_err)
call adios_write (filehandle, "time"//char(0), q0+0.5*q1+0.25*q2,adios_err)
call adios_write(filehandle,"3+mbin_psi*6"//char(0),3+mbin_psi*6,adios_err)
call adios_write(filehandle,"mbin_u"//char(0),mbin_u,adios_err)
call adios_write(filehandle,"mbin_psi"//char(0),mbin_psi,adios_err)
call adios_write(filehandle,"mpsi"//char(0),mpsi,adios_err)
call adios_write(filehandle,"jm"//char(0),jm,adios_err)
call adios_write(filehandle,"jm+1"//char(0),jm+1,adios_err)
call adios_write(filehandle,"mzetamax"//char(0),mzetamax,adios_err)
call adios_write(filehandle,"num_of_quantities"//char(0),i, adios_err)
call adios_write(filehandle,"dataout_u"//char(0),dataout_u,adios_err)
call adios_write(filehandle,"dataout_r"//char(0),dataout_r,adios_err)
call adios_write (filehandle, "dataout_p"//char(0), booz_out, adios_err)
call adios_write(filehandle,"ntoroidal"//char(0),ntoroidal,adios_err)
call adios_write(filehandle,"mzeta"//char(0),mzeta,adios_err)
call adios_write (filehandle, "dataout_f"//char(0), flux3darray, adios_err)
call adios_write(filehandle,"num_mode"//char(0),num_mode,adios_err)
call adios_write (filehandle, "modeeigen"//char(0), nmode, adios_err)
call adios_write(filehandle,"m_poloidal"//char(0),m_poloidal,adios_err)
call adios_write (filehandle, "dataeigen"//char(0), eigenmode, adios_err)


!     call adios_write (filehandle, "timestep"//char(0), (mstepall+istep),adios_err)
!     call adios_write (filehandle, "time"//char(0), q0+0.5*q1+0.25*q2,adios_err)
     !!!!if(mype==0)write(*,*)"filename,groupfilehandle,filehandle",bpdum,groupfilehandle,filehandle
!     call adios_write (filehandle, "3+mbin_psi*6"//char(0), 3+mbin_psi*6,adios_err)
!     call adios_write (filehandle, "jm+1"//char(0), jm+1, adios_err)
!     ADIOS_WRITE(filehandle,mbin_u)
!     ADIOS_WRITE(filehandle,mbin_psi)
!     ADIOS_WRITE(filehandle,mpsi)
!     ADIOS_WRITE(filehandle,jm)
!     ADIOS_WRITE(filehandle,mzetamax)
!     call adios_write(filehandle,"num_of_quantities"//char(0),i, adios_err)
     ! poloidal cross section: marker,density,potential,parallel flow and temperature
     !!call adios_write (filehandle, "jm+1"//char(0), jm+1)
!     ADIOS_WRITE(filehandle,jm)
!     ADIOS_WRITE(filehandle,ntoroidal)
!     ADIOS_WRITE(filehandle,mzeta)
!     ADIOS_WRITE(filehandle,num_mode)
!     ADIOS_WRITE(filehandle,m_poloidal)
!     call adios_write (filehandle, "modeeigen"//char(0), nmode, adios_err)
!     ADIOS_WRITE(filehandle,dataout_u)
!     ADIOS_WRITE(filehandle,dataout_r)
!     call adios_write (filehandle, "dataout_p"//char(0), booz_out, adios_err)
!     call adios_write (filehandle, "dataout_f"//char(0), flux3darray, adios_err)
!     call adios_write (filehandle, "dataeigen"//char(0), eigenmode, adios_err)

!     CALL write_end_for_group(group_proffilehandle,istep)

!     CALL close_start_for_group(group_proffilehandle,istep)

     call adios_close(filehandle, adios_err)

!     CALL close_end_for_group(group_proffilehandle,istep)

#endif
     deallocate(booz_out)  
     deallocate(dataout_u)  
     deallocate(dataout_r)  
     deallocate(flux3darray)  


  endif !!if(mype==0)

! record program end time
  if(istep .eq. mstep)then
     call date_and_time(date,time)
     if(mype == 0) then
        if(stdout /= 6 .and. stdout /= 0)open(stdout,file='stdout.out',status='old',position='append')
        write(stdout,*) 'Program ends at DATE=', date, 'TIME=', time
        if(stdout /= 6 .and. stdout /= 0)close(stdout)
     end if
  endif
  
101 format(i6)
102 format(e10.4)
  
end subroutine snapshot

