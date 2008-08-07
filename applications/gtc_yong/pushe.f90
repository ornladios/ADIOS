subroutine pushe(icycle,irke)
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use particle_tracking,only :track_particles
  implicit none

  integer i,im,ip,j,jt,kk,m,ii,j11,j10,j01,j00,jtheta,num,ierror,ij,&
       ipjt,icycle,irke
  real(wp) cost0,sint0,cost,sint,rdum,tdum,dtime,q,r,b,g,gp,ri,rip,dbdp,dbdt,&
       dedb,deni,upara,energy,kappa,dptdp,dptdt,dptdz,vdr,pdot,tdot,zdot,rdot,&
       wdot,wdrift,vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,pi2,uleft,uright,&
       adum(0:mpsi),pstat,xdum,ydum,xdot,ydot,g_inv,wz1,wz0,wp1,wp0,wt11,&
       wt10,wt01,wt00,temp(0:mpsi),dtemp(0:mpsi),denp(0:mpsi),g2,tem,rfac,&
       sbound,qdum,tflr,wpara,psimax,psimin,paxis,pdum,wdrive,ainv,rinv,qinv,&
       psitmp,thetatmp,zetatmp,rhoi,e1,e2,e3,e4,cmratio,cinv,pptpt,dzonal,&
       tem_inv,dtem(0:mpsi),dden(0:mpsi),ddum(0:mpsi),ddentmp(0:mpsi),&
       hfluxpsetmp(0:mpsi)
  real(doubleprec) d_efluxe,d_pfluxe,d_dflowe,d_entropye

  delr=1.0/deltar
  delt=2.0*pi/deltat
  delz=1.0/deltaz
  pi2=2.0*pi
  pi2_inv=0.5/pi
  g_inv=1.0/gyroradius
  sbound=1.0
  if(nbound==0)sbound=0.0
  psimax=0.5*a1*a1
  psimin=0.5*a0*a0
  paxis=0.5*(8.0*gyroradius)**2
  cmratio=qelectron/aelectron
  cinv=1.0/qelectron
  tem_inv=1.0/(gyroradius*gyroradius)

  if(irke==1)then
! 1st step of Runge-Kutta method
     dtime=0.25*tstep/real(ncycle)

     if(icycle==1)then
! save electron info
        if(irk*ihybrid==1)then
           ntracer1=ntracer
           me1=me
!$omp parallel do private(m)
           do m=1,me
              zelectron1(1:nparam,m)=zelectron(1:nparam,m)
           enddo
        else

! electron start from initial position
           ntracer=ntracer1
           me=me1
!$omp parallel do private(m)
           do m=1,me
              zelectron(1:nparam,m)=zelectron1(1:nparam,m)
           enddo
        endif
     endif

!$omp parallel do private(m)
     do m=1,me
        zelectron0(1:5,m)=zelectron(1:5,m)
     enddo

     if(track_particles == 1)then
        do m=1,me
          zelectron0(7:nparam,m)=zelectron(7:nparam,m)
        enddo
     endif

! 2nd step of Runge-Kutta method
  else
     dtime=0.5*tstep/real(ncycle)
     uright=umax*gyroradius*sqrt(aelectron)
     uleft=-uright
  endif

!$omp parallel do private(m,psitmp,thetatmp,zetatmp,r,ip,jt,ipjt,wz1,&
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

! gather e_field using 4-point gyro-averaging, sorting in poloidal angle
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
!$omp parallel do private(m,e1,e2,e3,e4,kk,wz1,wz0,ij,wp0,wt00,&
!$omp& wt10,wp1,wt01,wt11)
  do m=1,me
     e1=0.0
     e2=0.0
     e3=0.0
     e4=0.0
     kk=kzelectron(m)
     wz1=wzelectron(m)
     wz0=1.0-wz1

     ij=jtelectron0(m)
     wp0=1.0-wpelectron(m)
     wt00=1.0-wtelectron0(m)
     e1=e1+wp0*wt00*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
     e2=e2+wp0*wt00*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
     e3=e3+wp0*wt00*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     e4=e4+wp0*wt00*(wz0*phit(kk,ij)+wz1*phit(kk+1,ij))
        
     ij=ij+1
     wt10=1.0-wt00
     e1=e1+wp0*wt10*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
     e2=e2+wp0*wt10*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
     e3=e3+wp0*wt10*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     e4=e4+wp0*wt10*(wz0*phit(kk,ij)+wz1*phit(kk+1,ij))
        
     ij=jtelectron1(m)
     wp1=1.0-wp0
     wt01=1.0-wtelectron1(m)
     e1=e1+wp1*wt01*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
     e2=e2+wp1*wt01*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
     e3=e3+wp1*wt01*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     e4=e4+wp1*wt01*(wz0*phit(kk,ij)+wz1*phit(kk+1,ij))
        
     ij=ij+1
     wt11=1.0-wt01
     e1=e1+wp1*wt11*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
     e2=e2+wp1*wt11*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
     e3=e3+wp1*wt11*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     e4=e4+wp1*wt11*(wz0*phit(kk,ij)+wz1*phit(kk+1,ij))

     wpelectron(m)=e1
     wtelectron0(m)=e2
     wzelectron(m)=e3
     wtelectron1(m)=e4
  enddo

! marker temperature
!  denp=0.0
!  temp=0.0
!  dtemp=0.0
!  do m=1,mi
!     sint=sin(zelectron(2,m))
!     cost0=cos(zelectron(2,m)+r*sint)
!     b=1.0/(1.0+r*cost0)
!     upara=zelectron(4,m)*b
!     denp(ip)=denp(ip)+1.0
!     temp(ip)=temp(ip)+0.5*upara*upara+zelectron(6,m)*zelectron(6,m)*b
!     dtemp(ip)=dtemp(ip)+upara
!  enddo
!  call MPI_ALLREDUCE(denp,adum,mpsi,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
!  denp=adum
!  call MPI_ALLREDUCE(temp,adum,mpsi,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
!  temp=adum
!  call MPI_ALLREDUCE(dtemp,adum,mpsi,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
!  dtemp=adum
!  temp=g_inv*g_inv*temp/(1.5*denp)
!  dtemp=dtemp/denp
  temp=1.0
  dtemp=0.0

  temp=1.0/(temp*rteme)
  ainv=1.0/a
! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private(m,r,rinv,ii,wp0,wp1,tem,dzonal,q,qinv,&
!$omp& cost,sint,cost0,sint0,&
!$omp& b,g,gp,ri,rip,dbdp,dbdt,dedb,deni,upara,energy,rfac,kappa,dptdp,&
!$omp& dptdt,dptdz,pptpt,wpara,wdrift,vdr,wdrive,wdot,pdot,tdot,zdot,rdot)
  do m=1,me
     r=sqrt(2.0*zelectron(1,m))
     rinv=1.0/r
     ii=max(0,min(mpsi-1,int((r-a0)*delr)))
     wp0=real(ii+1)-(r-a0)*delr
     wp1=1.0-wp0
     tem=g_inv*g_inv*(wp0*temp(ii)+wp1*temp(ii+1))
     dzonal=wp0*phip00(ii)+wp1*phip00(ii+1)
     q=q0+q1*r*ainv+q2*r*r*ainv*ainv
     qinv=1.0/q
     cost=cos(zelectron(2,m))
     sint=sin(zelectron(2,m))
!     cost0=cos(zelectron(2,m)+r*sint)
!     sint0=sin(zelectron(2,m)+r*sint)
     b=1.0/(1.0+r*cost)
     g=1.0
     gp=0.0
!     ri=r*r*qinv
!     rip=(2.0*q0+q1*r*ainv)*qinv*qinv
     ri=0.0
     rip=0.0
     dbdp=-b*b*cost*rinv
     dbdt=b*b*r*sint
     dedb=cinv*(zelectron(4,m)*zelectron(4,m)*qelectron*b*cmratio+zelectron(6,m)*&
          zelectron(6,m))
     deni=1.0/(g*q + ri + zelectron(4,m)*(g*rip-ri*gp))
     upara=zelectron(4,m)*b*cmratio
     energy=0.5*aelectron*upara*upara+zelectron(6,m)*zelectron(6,m)*b
     rfac=rw*(r-rc)
     rfac=rfac*rfac
     rfac=rfac*rfac*rfac
     rfac=exp(-rfac)
     kappa=1.0-sbound+sbound*rfac
!     kappa=((min(umax*umax,energy*tem)-1.5)*kappate+kappan)*kappa*rinv
     kappa=((energy*tem-1.5)*kappate+kappan)*kappa*rinv

! perturbed quantities
     dptdp=wpelectron(m)
     dptdt=wtelectron0(m)
     dptdz=wzelectron(m)-wtelectron0(m)*qinv
     pptpt=wtelectron1(m)

! ExB drift in radial direction for w-dot and flux diagnostics
     vdr=q*(ri*dptdz-g*dptdt)*deni
     wdrive=vdr*kappa
     wpara=pptpt*qelectron*tem
     wdrift=q*(g*dedb*dbdt + g*dptdt - ri*dptdz)*deni*dzonal*qelectron*tem*rinv  !!XY
     wdot=(zelectron0(6,m)-paranl*zelectron(5,m))*(wdrive+wpara+wdrift)

! self-consistent and external electric field for marker orbits
     dptdp=dptdp*nonlinear+gyroradius*(flow0+flow1*r*ainv+flow2*r*r*ainv*ainv)
     dptdt=dptdt*nonlinear
     dptdz=dptdz*nonlinear

! particle velocity
     pdot = q*(-g*dedb*dbdt - g*dptdt + ri*dptdz)*deni
     tdot = (upara*b*(1.0-q*gp*zelectron(4,m)) + q*g*(dedb*dbdp + dptdp))*deni
     zdot = (upara*b*q*(1.0+rip*zelectron(4,m)) - q*ri*(dedb*dbdp + dptdp))*deni
     rdot = ((gp*zelectron(4,m)-1.0)*(dedb*dbdt + paranl*dptdt)-&
          paranl*q*(1.0+rip*zelectron(4,m))*dptdz)*deni 
         
! update particle position
!     if(zelectron0(1,m) < paxis)then
! particles close to axis use (x,y) coordinates
!        pdum=sqrt(zelectron(1,m))
!        xdum   = pdum*cost  
!        ydum   = pdum*sint
!        pdum=1.0/zelectron(1,m)
!        xdot   = 0.5*pdot*xdum*pdum-ydum*tdot
!        ydot   = 0.5*pdot*ydum*pdum+xdum*tdot
!        pdum=sqrt(zelectron0(1,m))
!        xdum   = pdum*cos(zelectron0(2,m)) + dtime*xdot
!        ydum   = pdum*sin(zelectron0(2,m)) + dtime*ydot
!        zelectron(1,m) = max(1.0e-8*psimax,xdum*xdum+ydum*ydum)
!        zelectron(2,m) = sign(1.0,ydum)*acos(max(-1.0,min(1.0,xdum/sqrt(zelectron(1,m)))))
!     else
     zelectron(1,m) = max(1.0e-8*psimax,zelectron0(1,m)+dtime*pdot)
     zelectron(2,m) = zelectron0(2,m)+dtime*tdot
!     endif

     zelectron(3,m) = zelectron0(3,m)+dtime*zdot
     zelectron(4,m) = zelectron0(4,m)+dtime*rdot
     zelectron(5,m) = zelectron0(5,m)+dtime*wdot

! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
! S.Ethier 04/13/2004  The hand-coded procedure gives wrong answers in some
! cases so it is better to use the modulo() function.

!!!     zelectron(2,m)=zelectron(2,m)*pi2_inv+10.0 !period of 1
!!!     zelectron(2,m)=2.0*pi*(zelectron(2,m)-aint(zelectron(2,m))) ![0,2*pi)
!!!     zelectron(3,m)=zelectron(3,m)*pi2_inv+10.0
!!!     zelectron(3,m)=2.0*pi*(zelectron(3,m)-aint(zelectron(3,m)))

     zelectron(2,m)=modulo(zelectron(2,m),pi2)
     zelectron(3,m)=modulo(zelectron(3,m),pi2)

! 02/20/2004  The modulo function seems to prevent streaming on the X1
! 02/20/2004 mod() does the same thing as modulo but it streams.
! 02/23/2004 need to do mod((pi2+zelectron),pi2) instead of
! mod(zelectron,pi2) in order to catch the negative values of zelectron.
!     zelectron(2,m)=mod((pi2+zelectron(2,m)),pi2)
!     zelectron(3,m)=mod((pi2+zelectron(3,m)),pi2)


! store GC information for flux measurements
     wpelectron(m)=vdr*rinv
     wtelectron0(m)=energy
     wzelectron(m)=b
 enddo
      
 if(irke==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,me
       if(zelectron(1,m) > psimax)then
          zelectron(1,m)=zelectron0(1,m)
          zelectron(2,m)=2.0*pi-zelectron0(2,m)
          zelectron(3,m)=zelectron0(3,m)
          zelectron(4,m)=zelectron0(4,m)
          zelectron(5,m)=zelectron0(5,m)
          
       elseif(zelectron(1,m) < psimin)then
          zelectron(1,m)=zelectron0(1,m)
          zelectron(2,m)=2.0*pi-zelectron0(2,m)
          zelectron(3,m)=zelectron0(3,m)
          zelectron(4,m)=zelectron0(4,m)
          zelectron(5,m)=zelectron0(5,m)
       endif
    enddo

! restore temperature profile
    if(nonlinear > 0.5)then
       if(mod(istep,ndiag)==0 .and. ihybrid==nhybrid .and. icycle==2*ncycle)then
          dtem=0.0
          dden=0.0
          do m=1,me
             r=sqrt(2.0*zelectron(1,m))
             ii=max(0,min(mpsi,int((r-a0)*delr+0.5)))
             dtem(ii)=dtem(ii)+wtelectron0(m)*zelectron(5,m)
             dden(ii)=dden(ii)+1.0
          enddo

          call MPI_ALLREDUCE(dtem,ddum,mpsi+1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
          dtem=ddum
          call MPI_ALLREDUCE(dden,ddum,mpsi+1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
          dden=ddum
          dtem=dtem*tem_inv/max(1.0,dden) !perturbed thermal energy
          do i=1,mpsi/2
             rdum=dtem(1)
             tdum=dtem(mpsi-1)
             dtem(1:mpsi-1)=0.5*dtem(1:mpsi-1)+0.25*(dtem(0:mpsi-2)+dtem(2:mpsi))
             dtem(0)=0.5*(dtem(0)+rdum)
             dtem(mpsi)=0.5*(dtem(mpsi)+tdum)
          enddo
          tdum=0.01*real(ndiag)
          rdteme=(1.0-tdum)*rdteme+tdum*dtem

!$omp parallel do private(m,r,ip)
          do m=1,me
             r=sqrt(2.0*zelectron(1,m))
             ii=max(0,min(mpsi,int((r-a0)*delr+0.5)))
             zelectron(5,m)=zelectron(5,m)-(wtelectron0(m)*tem_inv-1.5)*rdteme(ii)
          enddo
       endif
    endif

 endif

!!!!XY Y. Xiao 6/20/07 for electron heat flux radial profile
       
 if(idiag==0 .and. ihybrid*icycle*irke==1)then
! fluxes diagnose
    d_efluxe=0.0
    d_pfluxe=0.0
    d_dflowe=0.0
    d_entropye=0.0
    hfluxpse=0.0  !!!!XY
    dden=0.0  !!!!XY
    do m=1,me
       r=sqrt(2.0*zelectron(1,m))  !!!!XY
       ii=max(0,min(mpsi,int((r-a0)*delr+0.5))) !!!!XY
       dden(ii)=dden(ii)+1.0   !!!!XY

       vdrenergy=wpelectron(m)*(wtelectron0(m)-1.5*gyroradius*gyroradius)*&
            zelectron0(5,m) !assume T_e=1
       hfluxpse(ii)=hfluxpse(ii)+vdrenergy !!!!XY
       d_efluxe=d_efluxe+vdrenergy
       d_pfluxe=d_pfluxe+wpelectron(m)*zelectron0(5,m)
       d_dflowe=d_dflowe+wzelectron(m)*zelectron0(4,m)*zelectron0(5,m)
       d_entropye=d_entropye+zelectron0(5,m)*zelectron0(5,m)
    enddo
    efluxe=d_efluxe
    pfluxe=d_pfluxe
    dflowe=d_dflowe
    entropye=d_entropye
    
    hfluxpsetmp=0.0  !!!!XY
    call MPI_ALLREDUCE(hfluxpse,hfluxpsetmp,mpsi+1,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror) !!!!XY
    hfluxpse=hfluxpsetmp !!!!XY
    ddentmp=0.0   !!!!XY
    call MPI_ALLREDUCE(dden,ddentmp,mpsi+1,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror) !!!!XY
    dden=ddentmp !!!!XY
    hfluxpse=hfluxpse/dden !!!!XY


 endif
 
end subroutine pushe
