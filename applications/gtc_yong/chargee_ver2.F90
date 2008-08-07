subroutine chargee
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use particle_decomp
  implicit none

  integer m,i,im,j,k,ip,jt,kk,ii,j11,j10,j01,j00,ierror,larmor,ij,ipjt,&
       icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) weight,wemark,rdum,tdum,delr,delt(0:mpsi),delz,r,wz1,wz0,wp1,&
       wp0,wt11,wt10,wt01,wt00,adum(0:mpsi),tflr,damping,pi2_inv,&
       psitmp,thetatmp,zetatmp,deltheta
  real(wp) sendl(mgrid),recvr(mgrid)
  real(wp) dnitmp(0:mzeta,mgrid)

#ifdef _OPENMP
! We need the following temporary array only when using OpenMP
!!!!  real dnitmp(0:mzeta,mgrid)
#endif

  delr=1.0/deltar 
  delt=2.0*pi/deltat
  delz=1.0/deltaz
  pi2_inv=0.5/pi
  densitye=0.0

!$omp parallel do private(m,larmor,psitmp,thetatmp,zetatmp,r,ip,jt,ipjt,wz1,&
!$omp& kk,rdum,ii,wp1,tflr,im,tdum,j00,j01)
  do m=1,me
     psitmp=zelectron(1,m)
     thetatmp=zelectron(2,m)
     zetatmp=zelectron(3,m)

     r=sqrt(2.0*psitmp)
     ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
     jt=max(0,min(mtheta(ip),int(thetatmp*pi2_inv*delt(ip)+0.5)))
     ipjt=igrid(ip)+jt
     wz1=(zetatmp-zetamin)*delz
     kk=max(0,min(mzeta-1,int(wz1)))
     kzelectron(m)=kk
     wzelectron(m)=wz1-real(kk)

     rdum=delr*max(0.0,min(a1-a0,r-a0))
     ii=max(0,min(mpsi-1,int(rdum)))
     wp1=rdum-real(ii)
     wpelectron(m)=wp1

! particle position in theta
     tflr=thetatmp

! inner flux surface
     im=ii
     tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0
     tdum=(tdum-aint(tdum))*delt(im)
     j00=max(0,min(mtheta(im)-1,int(tdum)))
     jtelectron0(m)=igrid(im)+j00
     wtelectron0(m)=tdum-real(j00)

! outer flux surface
     im=ii+1
     tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0
     tdum=(tdum-aint(tdum))*delt(im)
     j01=max(0,min(mtheta(im)-1,int(tdum)))
     jtelectron1(m)=igrid(im)+j01
     wtelectron1(m)=tdum-real(j01)

  enddo	

#ifdef _OPENMP
! The following lines are OpenMP directives for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
! When compiling with the OpenMP option (-qsmp=omp for the XLF compiler on
! AIX and -mp for the MIPSpro compiler on IRIX), the preprocessor variable
! _OPENMP is automatically defined and we use it here to add pieces of codes
! needed to avoid contentions between threads. In the following parallel
! loop each thread has its own private array "dnitmp()" which is used to
! prevent the thread from writing into the shared array "densitye()"
!
!$omp parallel private(dnitmp)
  dnitmp=0.   ! Set array elements to zero
#endif
!$omp do private(m,larmor,weight,kk,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)
  do m=1,me
     weight=zelectron(5,m)

     kk=kzelectron(m)
     wz1=weight*wzelectron(m)
     wz0=weight-wz1     

     wp1=wpelectron(m)
     wp0=1.0-wp1

     wt10=wp0*wtelectron0(m)
     wt00=wp0-wt10

     wt11=wp1*wtelectron1(m)
     wt01=wp1-wt11

#ifdef _OPENMP
! Use thread-private temporary array dnitmp to store the results
     ij=jtelectron0(m)
     dnitmp(kk,ij) = dnitmp(kk,ij) + wz0*wt00
     dnitmp(kk+1,ij)   = dnitmp(kk+1,ij)   + wz1*wt00
        
     ij=ij+1
     dnitmp(kk,ij) = dnitmp(kk,ij) + wz0*wt10
     dnitmp(kk+1,ij)   = dnitmp(kk+1,ij)   + wz1*wt10

     ij=jtelectron1(m)
     dnitmp(kk,ij) = dnitmp(kk,ij) + wz0*wt01
     dnitmp(kk+1,ij)   = dnitmp(kk+1,ij)   + wz1*wt01
        
     ij=ij+1
     dnitmp(kk,ij) = dnitmp(kk,ij) + wz0*wt11
     dnitmp(kk+1,ij)   = dnitmp(kk+1,ij)   + wz1*wt11
#else
! If no loop-level parallelism, use original algorithm (write directly
! into array "densitye()".
     ij=jtelectron0(m)
     densitye(kk,ij) = densitye(kk,ij) + wz0*wt00
     densitye(kk+1,ij)   = densitye(kk+1,ij)   + wz1*wt00

     ij=ij+1
     densitye(kk,ij) = densitye(kk,ij) + wz0*wt10
     densitye(kk+1,ij)   = densitye(kk+1,ij)   + wz1*wt10

     ij=jtelectron1(m)
     densitye(kk,ij) = densitye(kk,ij) + wz0*wt01
     densitye(kk+1,ij)   = densitye(kk+1,ij)   + wz1*wt01

     ij=ij+1
     densitye(kk,ij) = densitye(kk,ij) + wz0*wt11
     densitye(kk+1,ij)   = densitye(kk+1,ij)   + wz1*wt11
#endif

  enddo

#ifdef _OPENMP
! For OpenMP, we now write the results accumulated in each thread-private
! array dnitmp() back into the shared array densitye(). The loop is enclosed
! in a critical section so that one thread at a time updates densitye(). 
!$omp critical
  do ij=1,mgrid
     do kk=0,mzeta
        densitye(kk,ij) = densitye(kk,ij) + dnitmp(kk,ij)
     enddo
  enddo
!$omp end critical
#endif
!$omp end parallel

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij,kk)
    do ij=1,mgrid
       do kk=0,mzeta
          dnitmp(kk,ij)=densitye(kk,ij)
          densitye(kk,ij)=0.
       enddo
    enddo
    call MPI_ALLREDUCE(dnitmp,densitye,(mgrid*(mzeta+1)),mpi_Rsize,&
                       MPI_SUM,partd_comm,ierror)
  endif


! poloidal end cell
  do i=0,mpsi
     densitye(:,igrid(i)+mtheta(i))=densitye(:,igrid(i)+mtheta(i))+densitye(:,igrid(i))
  enddo

! toroidal end cell
  sendl=densitye(0,:)
  recvr=0.0
  icount=mgrid
  !!idest=mod(mype-1+numberpe,numberpe)
  idest=left_pe
  !!isource=mod(mype+1,numberpe)
  isource=right_pe
  !!isendtag=mype
  isendtag=myrank_toroidal
  irecvtag=isource

! send densitye to left and receive from right
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
       recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
  if(myrank_toroidal == ntoroidal-1)then
! B.C. at zeta=2*pi is shifted
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        densitye(mzeta,ii+1:ii+jt)=densitye(mzeta,ii+1:ii+jt)+&
             cshift(recvr(ii+1:ii+jt),itran(i))
     enddo
  else
! B.C. at zeta<2*pi is continuous
     densitye(mzeta,:)=densitye(mzeta,:)+recvr
  endif
  
! zero out charge in radial boundary cell
  do i=0,nbound-1
     densitye(:,igrid(i):igrid(i)+mtheta(i))=&
        densitye(:,igrid(i):igrid(i)+mtheta(i))*real(i)/real(nbound)
     densitye(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=&
        densitye(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*real(i)/real(nbound)
  enddo

! flux surface average and normalization  
  zonale=0.0
!$omp parallel do private(i,j,k,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        do k=1,mzeta
           ij=igrid(i)+j
           zonale(i)=zonale(i)+densitye(k,ij)
           densitye(k,ij)=densitye(k,ij)*markere(k,ij)
        enddo        
     enddo
  enddo

! toroidal sum of phi00, broadcast to every PE
  call MPI_ALLREDUCE(zonale,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  zonale=adum*pmarke

! densitye subtracted (0,0) mode
!$omp parallel do private(i,j,k,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        do k=1,mzeta
           ij=igrid(i)+j
           densitye(k,ij)=densitye(k,ij)-zonale(i)
        enddo
     enddo
! poloidal BC condition
     densitye(1:mzeta,igrid(i))=densitye(1:mzeta,igrid(i)+mtheta(i))
  enddo
  
! enforce charge conservation for zonal flow mode
  rdum=0.0
  tdum=0.0
  do i=1,mpsi-1
     r=a0+deltar*real(i)
     rdum=rdum+r
     tdum=tdum+r*zonale(i)
  enddo
  tdum=tdum/rdum

 !! if(mype==0)then
    !!print *,'densitye=',densitye
 !!   print *,'zonale=',zonale
 !!   print *,'tdum=',tdum
 !! endif
  ddene=tdum !for diagnostic
  do i=1,mpsi-1
     zonale(i)=zonale(i)-tdum
  enddo

end subroutine chargee
