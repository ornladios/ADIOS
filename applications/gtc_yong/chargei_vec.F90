subroutine chargei
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use particle_decomp
  implicit none

! "vlen" represent the vector length on the vector machine.
  integer,parameter:: vlen=256,vlenp1=vlen+1

  integer m,i,im,j,k,ip,jt,kk,ii,j11,j10,j01,j00,ierror,larmor,ij,ipjt,&
       icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) weight,wemark,rdum,tdum,delr,delt(0:mpsi),delz,smu_inv,r,wz1,wz0,&
       wp1,wp0,wt11,wt10,wt01,wt00,adum(0:mpsi),tflr,damping,pi2_inv,&
       psitmp,thetatmp,zetatmp,rhoi,deltheta
  real(wp) sendl(mgrid),recvr(mgrid)
  real(wp) dnitmp(0:mzeta,mgrid),partial(vlenp1,0:mzeta,mgrid)

#ifdef _SX
!cdir duplicate_noinit(delt,256)
#endif

  integer jstat,imm

  delr=1.0/deltar 
  delt=2.0*pi/deltat
!cdir du_update(delt)

  delz=1.0/deltaz
  smu_inv=sqrt(aion)/(abs(qion)*gyroradius)
  pi2_inv=0.5/pi
  densityi=0.0

#ifdef _SX
! SX-6 profiling routines
  CALL FTRACE_REGION_BEGIN("CHARGEI_R1")
#endif

!$omp parallel do private(m,larmor,psitmp,thetatmp,zetatmp,rhoi,r,ip,jt,ipjt,&
!$omp& wz1,kk,rdum,ii,wp1,tflr,im,tdum,j00,j01)
  do m=1,mi
     psitmp=zion(1,m)
     thetatmp=zion(2,m)
     zetatmp=zion(3,m)
     rhoi=zion(6,m)*smu_inv
	
     r=sqrt(2.0*psitmp)
     ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
#ifdef _SX
     jt=max(0,min(mmtheta(ip),int(thetatmp*pi2_inv*delt(ip)+0.5)))
     ipjt=iigrid(ip)+jt   ! Trick to reduce bank conflicts
#else
     jt=max(0,min(mtheta(ip),int(thetatmp*pi2_inv*delt(ip)+0.5)))
     ipjt=igrid(ip)+jt
#endif
     wz1=(zetatmp-zetamin)*delz
     kk=max(0,min(mzeta-1,int(wz1)))
     kzion(m)=kk
     wzion(m)=wz1-real(kk)

!cdir expand=5
     do larmor=1,4
        rdum=delr*max(0.0,min(a1-a0,r+rhoi*pgyro(larmor,ipjt)-a0))
        ii=max(0,min(mpsi-1,int(rdum)))
        wp1=rdum-real(ii)
        wpion(m,larmor)=wp1

! particle position in theta
        tflr=thetatmp+rhoi*tgyro(larmor,ipjt)

! inner flux surface
        im=ii
#ifdef _SX
        tdum=pi2_inv*(tflr-zetatmp*qqtinv(im))+10.0  ! Reduces bank conflicts
        tdum=(tdum-aint(tdum))*delt(im)
        j00=max(0,min(mmtheta(im)-1,int(tdum))) ! Reduces bank conflicts
        jtion0(m,larmor)=iigrid(im)+j00 ! Reduces bank conflicts
        wtion0(m,larmor)=tdum-real(j00)
#else
        tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0
        tdum=(tdum-aint(tdum))*delt(im)
        j00=max(0,min(mtheta(im)-1,int(tdum)))
        jtion0(m,larmor)=igrid(im)+j00
        wtion0(m,larmor)=tdum-real(j00)
#endif

! outer flux surface
        im=ii+1
#ifdef _SX
        tdum=pi2_inv*(tflr-zetatmp*qqtinv(im))+10.0 ! Reduces bank conflicts
        tdum=(tdum-aint(tdum))*delt(im)
        j01=max(0,min(mmtheta(im)-1,int(tdum))) ! Reduces bank conflicts
        jtion1(m,larmor)=iigrid(im)+j01 ! Reduces bank conflicts
        wtion1(m,larmor)=tdum-real(j01)
#else
        tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0
        tdum=(tdum-aint(tdum))*delt(im)
        j01=max(0,min(mtheta(im)-1,int(tdum)))
        jtion1(m,larmor)=igrid(im)+j01
        wtion1(m,larmor)=tdum-real(j01)
#endif
     enddo		
  enddo	

#ifdef _SX
  CALL FTRACE_REGION_END("CHARGEI_R1")
#endif

  if(istep==0)return

  partial=0.   ! Set array elements to zero

#ifdef _SX
  CALL FTRACE_REGION_BEGIN("CHARGEI_R2")
#endif

  do m=1,mi,vlen
     do imm=1,min(vlen,mi-m+1)
        wz1=zion(5,m+imm-1)*wzion(m+imm-1)  !wz1=weight*wzion
        wz0=zion(5,m+imm-1)-wz1
        kk=kzion(m+imm-1)
        !!!do larmor=1,4
        !----------------------------------------------------
        !larmor=1
           wp1=wpion(m+imm-1,1)
           wp0=1.0-wp1

           wt10=wp0*wtion0(m+imm-1,1)
           wt00=wp0-wt10

           wt11=wp1*wtion1(m+imm-1,1)
           wt01=wp1-wt11

           ij=jtion0(m+imm-1,1)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt00
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt00

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt10
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt10

           ij=jtion1(m+imm-1,1)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt01
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt01

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt11
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt11

        !----------------------------------------------------
        !larmor=2
           wp1=wpion(m+imm-1,2)
           wp0=1.0-wp1

           wt10=wp0*wtion0(m+imm-1,2)
           wt00=wp0-wt10

           wt11=wp1*wtion1(m+imm-1,2)
           wt01=wp1-wt11

           ij=jtion0(m+imm-1,2)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt00
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt00

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt10
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt10

           ij=jtion1(m+imm-1,2)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt01
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt01

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt11
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt11

        !----------------------------------------------------
        !larmor=3
           wp1=wpion(m+imm-1,3)
           wp0=1.0-wp1

           wt10=wp0*wtion0(m+imm-1,3)
           wt00=wp0-wt10

           wt11=wp1*wtion1(m+imm-1,3)
           wt01=wp1-wt11

           ij=jtion0(m+imm-1,3)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt00
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt00

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt10
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt10

           ij=jtion1(m+imm-1,3)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt01
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt01

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt11
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt11

        !----------------------------------------------------
        !larmor=4
           wp1=wpion(m+imm-1,4)
           wp0=1.0-wp1

           wt10=wp0*wtion0(m+imm-1,4)
           wt00=wp0-wt10

           wt11=wp1*wtion1(m+imm-1,4)
           wt01=wp1-wt11

           ij=jtion0(m+imm-1,4)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt00
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt00

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt10
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt10

           ij=jtion1(m+imm-1,4)
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt01
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt01

           ij=ij+1
           partial(imm,kk,ij) = partial(imm,kk,ij) + wz0*wt11
           partial(imm,kk+1,ij)   = partial(imm,kk+1,ij)   + wz1*wt11

        !!!!!enddo
     enddo
  enddo

#ifdef _SX
  CALL FTRACE_REGION_END("CHARGEI_R2")

  CALL FTRACE_REGION_BEGIN("CHARGEI_R3")
#endif

! Merge into main array "densityi" the results stored into the temporary
! array "partial". For OpenMP, the loop is enclosed in a critical
! section so that one thread at a time updates densityi() since each
! thread has its own private copy of the array "partial".
!$omp critical
!cdir noloopchg
  do kk=0,mzeta
    do ij=1,mgrid
      do imm=1,vlen
        densityi(kk,ij) = densityi(kk,ij) + partial(imm,kk,ij)
      enddo
    enddo
  enddo

#ifdef _SX
  CALL FTRACE_REGION_END("CHARGEI_R3")
#endif

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij,kk)
    do ij=1,mgrid
       do kk=0,mzeta
          dnitmp(kk,ij)=densityi(kk,ij)
          densityi(kk,ij)=0.
       enddo
    enddo
    call MPI_ALLREDUCE(dnitmp,densityi,(mgrid*(mzeta+1)),mpi_Rsize,&
                       MPI_SUM,partd_comm,ierror)
  endif

! poloidal end cell, discard ghost cell j=0
  do i=0,mpsi
     densityi(:,igrid(i)+mtheta(i))=densityi(:,igrid(i)+mtheta(i))+densityi(:,igrid(i))
  enddo

! toroidal end cell
  sendl=densityi(0,:)
  recvr=0.0
  icount=mgrid
  !!!idest=mod(myrank_toroidal-1+ntoroidal,ntoroidal)
  idest=left_pe
  !!!isource=mod(myrank_toroidal+1,ntoroidal)
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource

! send densityi to left and receive from right
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
       recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
  if(myrank_toroidal == ntoroidal-1)then
! B.C. at zeta=2*pi is shifted
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        densityi(mzeta,ii+1:ii+jt)=densityi(mzeta,ii+1:ii+jt)+&
             cshift(recvr(ii+1:ii+jt),itran(i))
     enddo
  else
! B.C. at zeta<2*pi is continuous
     densityi(mzeta,:)=densityi(mzeta,:)+recvr
  endif
  
! zero out charge in radial boundary cell
  do i=0,nbound-1
     densityi(:,igrid(i):igrid(i)+mtheta(i))=&
	densityi(:,igrid(i):igrid(i)+mtheta(i))*real(i)/real(nbound)
     densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=&
	densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*real(i)/real(nbound)
  enddo

! flux surface average and normalization  
  zonali=0.0
!$omp parallel do private(i,j,k,ij)
  do i=0,mpsi
!cdir noloopchg
     do k=1,mzeta
        do j=1,mtheta(i)
           ij=igrid(i)+j
           zonali(i)=zonali(i)+0.25*densityi(k,ij)
           densityi(k,ij)=0.25*densityi(k,ij)*markeri(k,ij)
        enddo        
     enddo
  enddo

! global sum of phi00, broadcast to every toroidal PE
  call MPI_ALLREDUCE(zonali,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  zonali=adum*pmarki

! densityi subtracted (0,0) mode
!$omp parallel do private(i,j,k,ij)
  do i=0,mpsi
!cdir noloopchg
     do k=1,mzeta
        do j=1,mtheta(i)
           ij=igrid(i)+j
           densityi(k,ij)=densityi(k,ij)-zonali(i)
        enddo
     enddo
! poloidal BC condition
     densityi(1:mzeta,igrid(i))=densityi(1:mzeta,igrid(i)+mtheta(i))
  enddo
  
! enforce charge conservation for zonal flow mode
  rdum=0.0
  tdum=0.0
  do i=1,mpsi-1
     r=a0+deltar*real(i)
     rdum=rdum+r
     tdum=tdum+r*zonali(i)
  enddo
  tdum=tdum/rdum
  ddeni=tdum !for diagnostic
  do i=1,mpsi-1
     zonali(i)=zonali(i)-tdum
  enddo

end subroutine chargei

