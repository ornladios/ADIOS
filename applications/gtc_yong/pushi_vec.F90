subroutine pushi
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  implicit none

  integer i,im,ip,j,jt,k,kz,kk,m,ii,j11,j10,j01,j00,jtheta,num,ierror,&
       larmor,ij,ipjt
  real(wp) cost0,sint0,cost,sint,rdum,tdum,dtime,q,r,b,g,gp,ri,rip,dbdp,dbdt,&
       dedb,deni,upara,energy,kappa,dptdp,dptdt,dptdz,vdr,pdot,tdot,zdot,rdot,&
       wdot,wdrift,vdrenergy,delr,pi2,pi2_inv,uleft,uright,adum(0:mpsi),&
       dmark(mflux),pstat,xdum,ydum,xdot,ydot,tem_inv,wz1,wz0,wp1,wp0,wt11,&
       wt10,wt01,wt00,temp(0:mpsi),dtemp(0:mpsi),denp(0:mpsi),g2,tem,rfac,&
       sbound,qdum,tflr,d_inv,epara,wpara,psimax,psimin,paxis,pdum,wdrive,&
       ainv,rinv,qinv,psitmp,thetatmp,zetatmp,rhoi,e1,e2,e3,cmratio,cinv,&
       wpi(3,mi),rdel,vthi,dtem(mflux),dden(mflux),ddentmp(mflux),&
       dtemtmp(mflux),vdrtmp(mflux),dmarktmp(mflux),hfluxpsitmp(0:mpsi)


  delr=1.0/deltar
  pi2=2.0*pi
  !pi2_inv=0.5/pi
  sbound=1.0
  if(nbound==0)sbound=0.0
  psimax=0.5*a1*a1
  psimin=0.5*a0*a0
  paxis=0.5*(8.0*gyroradius)**2
  cmratio=qion/aion
  cinv=1.0/qion
  vthi=gyroradius*abs(qion)/aion
  tem_inv=1.0/(aion*vthi*vthi)
  d_inv=real(mflux)/(a1-a0)

  if(irk==1)then
! 1st step of Runge-Kutta method
     dtime=0.5*tstep
!$omp parallel do private(m)
     do m=1,mi
        zion0(1:5,m)=zion(1:5,m)
     enddo

     vdrtmp=0.0

! 2nd step of Runge-Kutta method
  else
     dtime=tstep
     uright=umax*vthi
     uleft=-uright

     if(nonlinear<0.5)vdrtmp=0.0
     if(nonlinear>0.5)vdrtmp=pfluxpsi
  endif

#ifdef _SX
  CALL FTRACE_REGION_BEGIN("PUSHI_R1")
#endif

! gather e_field using 4-point gyro-averaging, sorting in poloidal angle
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
!$omp parallel do private(m,e1,e2,e3,kk,wz1,wz0,larmor,ij,wp0,wt00,&
!$omp& wt10,wp1,wt01,wt11)
  do m=1,mi
     e1=0.0
     e2=0.0
     e3=0.0
     kk=kzion(m)
     wz1=wzion(m)
     wz0=1.0-wz1

     !!!!do larmor=1,4
     !----------------------------------------------------------------
     !larmor=1

        ij=jtion0(m,1)
        wp0=1.0-wpion(m,1)
        wt00=1.0-wtion0(m,1)
        e1=e1+wp0*wt00*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt00*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt00*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt10=1.0-wt00
        e1=e1+wp0*wt10*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt10*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt10*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=jtion1(m,1)
        wp1=1.0-wp0
        wt01=1.0-wtion1(m,1)
        e1=e1+wp1*wt01*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt01*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt01*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt11=1.0-wt01
        e1=e1+wp1*wt11*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt11*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt11*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     !----------------------------------------------------------------
     !larmor=2

        ij=jtion0(m,2)
        wp0=1.0-wpion(m,2)
        wt00=1.0-wtion0(m,2)
        e1=e1+wp0*wt00*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt00*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt00*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt10=1.0-wt00
        e1=e1+wp0*wt10*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt10*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt10*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=jtion1(m,2)
        wp1=1.0-wp0
        wt01=1.0-wtion1(m,2)
        e1=e1+wp1*wt01*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt01*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt01*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt11=1.0-wt01
        e1=e1+wp1*wt11*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt11*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt11*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     !----------------------------------------------------------------
     !larmor=3

        ij=jtion0(m,3)
        wp0=1.0-wpion(m,3)
        wt00=1.0-wtion0(m,3)
        e1=e1+wp0*wt00*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt00*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt00*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt10=1.0-wt00
        e1=e1+wp0*wt10*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt10*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt10*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=jtion1(m,3)
        wp1=1.0-wp0
        wt01=1.0-wtion1(m,3)
        e1=e1+wp1*wt01*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt01*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt01*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt11=1.0-wt01
        e1=e1+wp1*wt11*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt11*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt11*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     !----------------------------------------------------------------
     !larmor=4

        ij=jtion0(m,4)
        wp0=1.0-wpion(m,4)
        wt00=1.0-wtion0(m,4)
        e1=e1+wp0*wt00*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt00*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt00*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt10=1.0-wt00
        e1=e1+wp0*wt10*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp0*wt10*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp0*wt10*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=jtion1(m,4)
        wp1=1.0-wp0
        wt01=1.0-wtion1(m,4)
        e1=e1+wp1*wt01*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt01*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt01*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
        
        ij=ij+1
        wt11=1.0-wt01
        e1=e1+wp1*wt11*(wz0*evector(1,kk,ij)+wz1*evector(1,kk+1,ij))
        e2=e2+wp1*wt11*(wz0*evector(2,kk,ij)+wz1*evector(2,kk+1,ij))
        e3=e3+wp1*wt11*(wz0*evector(3,kk,ij)+wz1*evector(3,kk+1,ij))
     !----------------------------------------------------------------

     !!!!enddo

     wpi(1,m)=0.25*e1
     wpi(2,m)=0.25*e2
     wpi(3,m)=0.25*e3

  enddo

#ifdef _SX
  CALL FTRACE_REGION_END("PUSHI_R1")
#endif

! primary ion marker temperature and parallel flow velocity
  temp=1.0
  dtemp=0.0
  temp=1.0/(temp*rtemi*aion*vthi*vthi) !inverse local temperature
  ainv=1.0/a

#ifdef _SX
 CALL FTRACE_REGION_BEGIN("PUSHI_R2")
#endif

! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private(m,r,rinv,ii,ip,wp0,wp1,tem,q,qinv,cost,sint,cost0,&
!$omp& sint0,b,g,gp,ri,rip,dbdp,dbdt,dedb,deni,upara,energy,rfac,kappa,dptdp,&
!$omp& dptdt,dptdz,epara,vdr,wdrive,wpara,wdrift,wdot,pdot,tdot,zdot,rdot)
  do m=1,mi
     r=sqrt(2.0*zion(1,m))
     rinv=1.0/r
     ii=max(0,min(mpsi-1,int((r-a0)*delr)))
     ip=max(1,min(mflux,1+int((r-a0)*d_inv)))
     wp0=real(ii+1)-(r-a0)*delr
     wp1=1.0-wp0
     tem=wp0*temp(ii)+wp1*temp(ii+1)
     q=q0+q1*r*ainv+q2*r*r*ainv*ainv
     qinv=1.0/q
     cost=cos(zion(2,m))
     sint=sin(zion(2,m))
!     cost0=cos(zion(2,m)+r*sint)
!     sint0=sin(zion(2,m)+r*sint)
     b=1.0/(1.0+r*cost)
     g=1.0
     gp=0.0
!     ri=r*r*qinv
!     rip=(2.0*q0+q1*r*ainv)*qinv*qinv
     ri=0.0
     rip=0.0
     dbdp=-b*b*cost*rinv
     dbdt=b*b*r*sint
     dedb=cinv*(zion(4,m)*zion(4,m)*qion*b*cmratio+zion(6,m)*zion(6,m))
     deni=1.0/(g*q + ri + zion(4,m)*(g*rip-ri*gp))
     upara=zion(4,m)*b*cmratio
     energy=0.5*aion*upara*upara+zion(6,m)*zion(6,m)*b
     rfac=rw*(r-rc)
     rfac=rfac*rfac
     rfac=rfac*rfac*rfac
     rfac=exp(-rfac)
     kappa=1.0-sbound+sbound*rfac
!     kappa=((min(umax*umax,energy*tem)-1.5)*kappati+kappan)*kappa*rinv
     kappa=((energy*tem-1.5)*kappati+kappan)*kappa*rinv

! perturbed quantities
     dptdp=wpi(1,m)
     dptdt=wpi(2,m)
     dptdz=wpi(3,m)-wpi(2,m)*qinv
     epara=-wpi(3,m)*b*q*deni

! subtract net particle flow
     dptdt=dptdt+vdrtmp(ip)

! ExB drift in radial direction for w-dot and flux diagnostics
     vdr=q*(ri*dptdz-g*dptdt)*deni
     wdrive=vdr*kappa
     wpara=epara*(upara-dtemp(ii))*qion*tem
     wdrift=q*(g*dbdt*dptdp-g*dbdp*dptdt+ri*dbdp*dptdz)*deni*dedb*qion*tem
     wdot=(zion0(6,m)-paranl*zion(5,m))*(wdrive+wpara+wdrift)

! self-consistent and external electric field for marker orbits
     dptdp=dptdp*nonlinear+gyroradius*(flow0+flow1*r*ainv+flow2*r*r*ainv*ainv)
     dptdt=dptdt*nonlinear
     dptdz=dptdz*nonlinear

! particle velocity
     pdot = q*(-g*dedb*dbdt - g*dptdt + ri*dptdz)*deni
     tdot = (upara*b*(1.0-q*gp*zion(4,m)) + q*g*(dedb*dbdp + dptdp))*deni
     zdot = (upara*b*q*(1.0+rip*zion(4,m)) - q*ri*(dedb*dbdp + dptdp))*deni
     rdot = ((gp*zion(4,m)-1.0)*(dedb*dbdt + paranl*dptdt)-&
          paranl*q*(1.0+rip*zion(4,m))*dptdz)*deni 
         
! update particle position
!     if(zion0(1,m) < paxis)then
! particles close to axis use (x,y) coordinates
!        pdum=sqrt(zion(1,m))
!        xdum   = pdum*cost  
!        ydum   = pdum*sint
!        pdum=1.0/zion(1,m)
!        xdot   = 0.5*pdot*xdum*pdum-ydum*tdot
!        ydot   = 0.5*pdot*ydum*pdum+xdum*tdot
!        pdum=sqrt(zion0(1,m))
!        xdum   = pdum*cos(zion0(2,m)) + dtime*xdot
!        ydum   = pdum*sin(zion0(2,m)) + dtime*ydot
!        zion(1,m) = max(1.0e-8*psimax,xdum*xdum+ydum*ydum)
!        zion(2,m) = sign(1.0,ydum)*acos(max(-1.0,min(1.0,xdum/sqrt(zion(1,m)))))
!     else
     zion(1,m) = max(1.0e-8*psimax,zion0(1,m)+dtime*pdot)
     zion(2,m) = zion0(2,m)+dtime*tdot
!     endif

     zion(3,m) = zion0(3,m)+dtime*zdot
     zion(4,m) = zion0(4,m)+dtime*rdot
     zion(5,m) = zion0(5,m)+dtime*wdot
     
! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!!!     zion(2,m)=zion(2,m)*pi2_inv+10.0 !period of 1
!!!     zion(2,m)=2.0*pi*(zion(2,m)-aint(zion(2,m))) ![0,2*pi)
!!!     zion(3,m)=zion(3,m)*pi2_inv+10.0
!!!     zion(3,m)=2.0*pi*(zion(3,m)-aint(zion(3,m)))

!     zion(2,m)=modulo(zion(2,m),pi2)
!     zion(3,m)=modulo(zion(3,m),pi2)

! 02/20/2004  The modulo function seems to prevent streaming on the X1
! 02/20/2004 mod() does the same thing as modulo but it streams.
! 02/23/2004 need to do mod((pi2+zion),pi2) instead of mod(zion,pi2) in order
!  to catch the negative values of zion.
     zion(2,m)=mod((pi2+zion(2,m)),pi2)
     zion(3,m)=mod((pi2+zion(3,m)),pi2)

! store GC information for flux measurements
     wpi(1,m)=vdr*rinv
     wpi(2,m)=energy
     wpi(3,m)=b
 enddo 
     
#ifdef _SX
  CALL FTRACE_REGION_END("PUSHI_R2")
#endif

 if(irk==2)then

#ifdef _SX
  CALL FTRACE_REGION_BEGIN("PUSHI_R3")
#endif

! out of boundary particle
!$omp parallel do private(m)
    do m=1,mi
       if(zion(1,m) > psimax)then
          zion(1,m)=zion0(1,m)
          zion(2,m)=2.0*pi-zion0(2,m)
          zion(3,m)=zion0(3,m)
          zion(4,m)=zion0(4,m)
          zion(5,m)=zion0(5,m)
          
       elseif(zion(1,m) < psimin)then
          zion(1,m)=zion0(1,m)
          zion(2,m)=2.0*pi-zion0(2,m)
          zion(3,m)=zion0(3,m)
          zion(4,m)=zion0(4,m)
          zion(5,m)=zion0(5,m)
       endif
    enddo

#ifdef _SX
  CALL FTRACE_REGION_END("PUSHI_R3")
#endif

! Restore temperature profile when running a nonlinear calculation
! (nonlinear=1.0) without the velocity space nonlinearity (paranl=0.0).
    if(nonlinear > 0.5 .and. paranl < 0.5)then
       if(mod(istep,ndiag)==0)then
          dtem=0.0
          dden=0.0
!$omp parallel do private(m,cost,b,upara)
          do m=1,mi
             wpi(1,m)=sqrt(2.0*zion(1,m))
             cost=cos(zion(2,m))
             b=1.0/(1.0+wpi(1,m)*cost)
             upara=zion(4,m)*b*cmratio
             wpi(2,m)=0.5*aion*upara*upara+zion(6,m)*zion(6,m)*b
          enddo
          do m=1,mi
             ip=max(1,min(mflux,1+int((wpi(1,m)-a0)*d_inv)))
             dtem(ip)=dtem(ip)+wpi(2,m)*zion(5,m)
             dden(ip)=dden(ip)+1.0
          enddo
          dtemtmp=0.
          ddentmp=0.
! S.Ethier 06/04/03  According to the MPI standard, the send and receive
! buffers cannot be the same for MPI_Reduce or MPI_Allreduce. It generates
! an error on the Linux platform. We thus make sure that the 2 buffers
! are different.
          call MPI_ALLREDUCE(dtem,dtemtmp,mflux,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
          call MPI_ALLREDUCE(dden,ddentmp,mflux,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
          dtem=dtemtmp*tem_inv/max(1.0,ddentmp) !perturbed temperature
          tdum=0.01*real(ndiag)
          rdtemi=(1.0-tdum)*rdtemi+tdum*dtem

!$omp parallel do private(m,ip)
          do m=1,mi
             ip=max(1,min(mflux,1+int((wpi(1,m)-a0)*d_inv)))
             zion(5,m)=zion(5,m)-(wpi(2,m)*tem_inv-1.5)*rdtemi(ip)
         ! Correction to include paranl effect according to Weixing
         ! 08/27/04 The correction may not be necessary according to Jerome
         !          because the restoration uses delta_f/f0
         !!!    zion(5,m)=zion(5,m)-(1.0-paranl*zion(5,m))*&
         !!!                 (wpi(2,m)*tem_inv-1.5)*rdtemi(ip)
          enddo
       endif
    endif
 endif

 if(idiag==0)then
! fluxes diagnose at irk=1
    rmarker=0.0
    eflux=0.0
    pfluxi=0.0
    dflowi=0.0
    entropyi=0.0
    hfluxpsi=0.0
    dmark=0.0
    dden=0.0
    particles_energy=0.0
    do m=1,mi
       r=sqrt(2.0*zion0(1,m))

! radial location in diagnostic bin
       ip=max(1,min(mflux,1+int((r-a0)*d_inv)))
       ii=max(0,min(mpsi,int((r-a0)*delr+0.5)))
       vdrenergy=wpi(1,m)*(wpi(2,m)-1.5*aion*vthi*vthi*rtemi(ii))*zion0(5,m)

! radial profile of heat flux
       hfluxpsi(ii)=hfluxpsi(ii)+vdrenergy ! energy flux profile       

! marker,energy,particle,momentum fluxes,parallel flows,entropy and kinetic energy
       rmarker(ip)=rmarker(ip)+zion0(6,m)
       eflux(ip)=eflux(ip)+vdrenergy
       pfluxi=pfluxi+wpi(1,m)*zion0(5,m)
       dflowi=dflowi+wpi(3,m)*zion0(4,m)*zion0(5,m)
       entropyi=entropyi+zion0(5,m)*zion0(5,m)
     ! S.Ethier 9/21/04  We want to keep track of the total energy in the
     ! particles, so we add up the contributions of energy*weight for each
     ! particle. This kinetic energy was saved in wpi(2,m) in the irk=1
     ! section above.
       particles_energy(1)=particles_energy(1)+wpi(2,m)*zion0(5,m)
       particles_energy(2)=particles_energy(2)+wpi(2,m)

       dmark(ip)=dmark(ip)+wpi(1,m)*r
       dden(ip)=dden(ip)+1.0
    enddo
    hfluxpsitmp=0.
    call MPI_ALLREDUCE(hfluxpsi,hfluxpsitmp,mpsi+1,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
    hfluxpsi=hfluxpsitmp*pmarki

    dmarktmp=0.
    call MPI_ALLREDUCE(dmark,dmarktmp,mflux,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
    ddentmp=0.
    call MPI_ALLREDUCE(dden,ddentmp,mflux,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
    dmark=dmarktmp/max(1.0,ddentmp)
    tdum=0.01*real(ndiag)
    pfluxpsi=(1.0-tdum)*pfluxpsi+tdum*dmark
 endif
 
end subroutine pushi
