subroutine diagnosis
  use global_parameters
  use particle_array
  use particle_decomp
  use field_array
  use diagnosis_array
  use data_type
  use mpiio_params
  implicit none

  integer,parameter :: mquantity=20
  integer i,j,k,ij,m,ierror,istatus(MPI_STATUS_SIZE),file_hd,startidx,endidx
  real(wp) r,b,g,q,kappa,fdum(21),adum(21),ddum(mpsi),ddumtmp(mpsi),tracer(6),&
       tmarker,xnormal(mflux),rfac,dum,epotential,vthi,vthe,tem_inv
  real(wp),external :: boozer2x,boozer2z,bfield,psi2r
  character(len=20) cdum
  character(len=30) hdfdum
  real dataout(57)
  integer*8 file_handle, io_type
  integer fh, err
  integer*8 offset ! since on first write we prepend, we need offset for second
#define ADIOS_WRITE(a,b) call adios_write (a,'b'//char(0),b)
  save xnormal
#define __NEWC
#ifdef __NERSC
#define FLUSH flush_
#else
#define FLUSH flush
#endif

! send tracer information to PE=0
! Make sure that PE 0 does not send to itself if it already holds the tracer.
! The processor holding the tracer particle has ntracer>0 while all the
! other processors have ntracer=0.
  tracer=0.0
  if(ntracer>0)then
     if(nhybrid==0)then
        tracer(1:5)=zion0(1:5,ntracer)
        tracer(6)=zion(6,ntracer)
     else
        tracer(1:5)=zelectron0(1:5,ntracer)
        tracer(6)=zelectron(6,ntracer)
     endif
     if(mype/=0)call MPI_SEND(tracer,6,mpi_Rsize,0,1,MPI_COMM_WORLD,ierror)
  endif
  if(mype==0 .and. ntracer==0)call MPI_RECV(tracer,6,mpi_Rsize,MPI_ANY_SOURCE,&
                  1,MPI_COMM_WORLD,istatus,ierror)

! initialize output file


  if(istep==ndiag)then
     if(mype==0)then
        if(irun==0)then
             offset = 0
     open(666,FILE="indextable",position='append',status='replace')
     write(666,101)istep,offset
     close(666)
             ! first time run 
             mstepall=0
             mstepdiag=0
             !write(hdfdum, '("DIAGNOSIS_",i5.5,".bp")')mstepall+istep
             !hdfdum=trim(hdfdum)//char(0)
#ifdef ADIOS
             hdfdum="diagnosis.bp"//char(0)
             call adios_get_group (io_type, "diagnosis.0"//char(0))
             call adios_open (file_handle, io_type, hdfdum)
             !call adios_get_data_size(file_handle,offset) ! use below to append
             ADIOS_WRITE(file_handle,irun)
             ADIOS_WRITE(file_handle,mquantity)
             ADIOS_WRITE(file_handle,mflux)
             ADIOS_WRITE(file_handle,num_mode)
             call adios_write (file_handle,"mstep_ndiag"//char(0), mstep/real(ndiag))
             call adios_write (file_handle,"tstep_ndiag"//char(0), tstep*real(ndiag))
             call adios_get_data_size (file_handle,offset) ! use below to append
             call adios_close (file_handle)
#endif
             ! E X B history
! tracer particle energy and angular momentum
! S.Ethier 4/15/04  The calculation of "etracer" and "ptracer" is now done
! only for irun=0 since we save the values of these variables in the
! restart files. This ensures continuity of the tracer particle energy
! and momentum plots.
             b=bfield(tracer(1),tracer(2))
             g=1.0
             r=sqrt(2.0*tracer(1))
             q=q0+q1*r/a+q2*r*r/(a*a)
!            epotential=gyroradius*r*r*(0.5*flow0+flow1*r/(3.0*a)+0.25*flow2*r*r/(a*a))
             epotential=0.0
             if(nhybrid==0)etracer=tracer(6)*tracer(6)*b &
                  +0.5*(tracer(4)*b*qion)**2/aion+epotential
             if(nhybrid>0)etracer=tracer(6)*tracer(6)*b &
                  +0.5*(tracer(4)*b*qelectron)**2/aelectron+epotential
             ptracer(1)=tracer(1)
             ptracer(2)=tracer(4)*g
             ptracer(3)=q
             ptracer(4)=0.0
             write(stdout,*)'etracer =',etracer
        else
#ifdef ADIOS
             !write(hdfdum, '("DIAGNOSIS_",i5.5,".bp")')mstepdiag
             !hdfdum=trim(hdfdum)//char(0)
             hdfdum="diagnosis.bp"//char(0)
             call adios_get_group (io_type, "diagnosis.1"//char(0))
             call adios_open (file_handle, io_type, hdfdum)

             ADIOS_WRITE(file_handle,irun)
             call adios_write (file_handle, "mstep_ndiag"//char(0), (mstep/real(ndiag)))
             call adios_close (file_handle)
#endif
! E X B history
! restart, copy old data to backup file to protect previous run data
! # of run

! if everything goes OK, copy old data to history.out
        endif

! flux normalization
        do i=1,mflux
             r=a0+(a1-a0)/real(mflux)*(real(i)-0.5)

! divided by local Kappa_T to obtain chi_i=c_i rho_i
             rfac=rw*(r-rc)
             rfac=rfac*rfac
             rfac=rfac*rfac*rfac
             rfac=max(0.1,exp(-rfac))
           
             kappa=1.0
             if(nbound==0)kappa=0.0
             kappa=1.0-kappa+kappa*rfac
           
             xnormal(i)=1.0/(kappa*kappati*gyroradius)
             write(stdout,*)"kappa_T at radial_bin=",i,kappa*kappati
        enddo
        do i=1,mpsi
             r=a0+(a1-a0)/real(mpsi)*(real(i)-0.5)

! divided by local Kappa_T to obtain chi_i=c_i rho_i
             rfac=rw*(r-rc)
             rfac=rfac*rfac
             rfac=rfac*rfac*rfac
             rfac=max(0.1,exp(-rfac))

             kappa=1.0
             if(nbound==0)kappa=0.0
             kappa=1.0-kappa+kappa*rfac
           
             gradt(i)=1.0/(kappa*kappati*gyroradius)
        enddo
     endif
     call MPI_BCAST(mstepall,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  endif
! global sum of fluxes
! All these quantities come from summing up contributions from the particles
! so the MPI_REDUCE involves all the MPI processes.
  vthi=gyroradius*abs(qion)/aion
  vthe=vthi*sqrt(aion/(aelectron*tite))
  tem_inv=1.0_wp/(aion*vthi*vthi)
  fdum(1)=efield
  fdum(2)=entropyi
  fdum(3)=entropye
  fdum(4)=dflowi/vthi
  fdum(5)=dflowe/vthe
  fdum(6)=pfluxi/vthi
  fdum(7)=pfluxe/vthi
  fdum(8)=efluxi*tem_inv/vthi
  fdum(9)=efluxe*tem_inv/vthi
  fdum(10:14)=eflux*tem_inv/vthi
  fdum(15:19)=rmarker
! S.Ethier 9/21/04 Add the total energy in the particles to the list
! of variables to reduce.
  fdum(20)=particles_energy(1)
  fdum(21)=particles_energy(2)
  call MPI_REDUCE(fdum,adum,21,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)

! radial envelope of fluctuation intensity
  ddum=0.0_wp
  do i=1,mpsi
           ddum(i)=ddum(i)+sum((phi(1:mzeta,igrid(i)+1:igrid(i)+mtheta(i)))**2)
  enddo
  call MPI_REDUCE(ddum,ddumtmp,mpsi,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  if(mype==0 .and. (mstepall+istep)>mstepdiag)then
     fdum=adum
     do i=1,mpsi
        ddum(i)=tem_inv*sqrt(ddumtmp(i)/real(mzetamax*mtheta(i)))
     enddo
! record tracer information
     dataout(1) = boozer2x(tracer(1),tracer(2))-1.0
     dataout(2) = boozer2z(tracer(1),tracer(2))
     dataout(3) = tracer(2)
     dataout(4) = tracer(3)
     dataout(5) = tracer(5) ! tracer weight
     b=bfield(tracer(1),tracer(2))
     g=1.0
     r=sqrt(2.0*tracer(1))
     q=q0+q1*r/a+q2*r*r/(a*a)
     
     if(nhybrid==0)dataout(6) = tracer(4)*b*qion/(aion*vthi)
     if(nhybrid>0)dataout(6) = tracer(4)*b*qelectron/(aelectron*vthe)
! energy error
     epotential=0.0

     if(nhybrid==0)dataout(7) = (tracer(6)*tracer(6)*b &
          +0.5*(tracer(4)*b*qion)**2/aion+epotential)/etracer-1.0
     if(nhybrid>0)dataout(7) = (tracer(6)*tracer(6)*b &
          +0.5*(tracer(4)*b*qelectron)**2/aelectron+epotential)/etracer-1.0

     dum=(tracer(1)-ptracer(1))/(0.5*(ptracer(3)+q))-(tracer(4)*g-ptracer(2))
     ptracer(1)=tracer(1)
     ptracer(2)=tracer(4)*g
     ptracer(3)=q
     ptracer(4)=ptracer(4)+dum
     dataout(8) = ptracer(4) ! momentum error
! normalization
     fdum(15:19)=max(1.0,fdum(15:19))
   !!  fdum(1)=fdum(1)/real(numberpe)
   ! At this point, fdum(1)/real(numberpe) is the volume average of phi**2.
   ! Since sqrt(y1)+sqrt(y2) is not equal to sqrt(y1+y2), we need to take
   ! the square root after the MPI_Reduce operation and final volume average.
   ! Modifications to the calculation of efield were also made in smooth.
     fdum(1)=sqrt(fdum(1)/real(numberpe))
     tmarker=sum(fdum(15:19))
     fdum(2:9)=fdum(2:9)/tmarker
     fdum(10:14)=fdum(10:14)*xnormal(1:5)/fdum(15:19)
! write fluxes to history file
     dataout(9) = ddeni
     dataout(10) = ddene
     dataout(11) = eradial
     dataout(12:25) = fdum(1:14)
     do i=1,2
        do j=1,num_mode
           do k=1,2
              dataout(25+((i-1)*num_mode+j-1)*2+k) = amp_mode(k,j,i)
           enddo
        enddo
     enddo
     open(666,file="diag_tstep.dat",status='unknown')
     write(666,101)mstepall+(istep/msnap)*msnap,mstepall+(istep/ndiag)*ndiag
     close(666)
#ifdef ADIOS
     call adios_get_group (io_type, "diagnosis.2"//char(0))
     hdfdum="diagnosis.bp"//char(0)
     call adios_open_append (file_handle, io_type, hdfdum, offset)
     call adios_write (file_handle, "index_timestep"//char(0), (mstepall+istep)/ndiag)
     call adios_write (file_handle, "axisx"//char(0), real(boozer2x(tracer(1),tracer(2))-1.0))
     call adios_write (file_handle, "axisz"//char(0), real(boozer2z(tracer(1),tracer(2))))
     call adios_write (file_handle, "tracer2"//char(0), real (tracer(2)))
     call adios_write (file_handle, "tracer3"//char(0), tracer(3))
     call adios_write (file_handle, "tracer5"//char(0), tracer(5))
     if(nhybrid==0) then
       call adios_write (file_handle, "tracer4"//char(0), tracer(4)*b*qion/(aion*vthi))
     endif 
     if(nhybrid>0) then
       call adios_write (file_handle, "tracer4"//char(0), tracer(4)*b*qelectron/(aelectron*vthe))
     endif
     if(nhybrid==0) then
       call adios_write (file_handle, "tracer6"//char(0), &
        (tracer(6)*tracer(6)*b+0.5*(tracer(4)*b*qion)**2/aion+epotential)/etracer-1.0)
     endif
     if(nhybrid>0) then
       call adios_write (file_handle, "tracer6"//char(0), &
        (tracer(6)*tracer(6)*b+0.5*(tracer(4)*b*qelectron)**2/aelectron+epotential)/etracer-1.0)
     endif
     call adios_write (file_handle, "ptracer4"//char(0), ptracer(4))
     ADIOS_WRITE(file_handle,ddeni)
     ADIOS_WRITE(file_handle,ddene)
     ADIOS_WRITE(file_handle,eradial)
     ADIOS_WRITE(file_handle,fdum)
     call adios_write (file_handle, "/history/dataset"//char(0), dataout)
     ADIOS_WRITE(file_handle,num_mode)
     ADIOS_WRITE(file_handle,amp_mode)
! radial-time 2D data
     ADIOS_WRITE(file_handle,mpsi)
     call adios_write (file_handle, "phip00"//char(0), phip00(1:mpsi)/vthi)
     ADIOS_WRITE(file_handle,ddum)
     ADIOS_WRITE(file_handle,hfluxpsi)
     call adios_close (file_handle)
#endif
     write(stdout,*)istep+mstepall,fdum(1),eradial,fdum(12),fdum(20)
endif
101 format(i16)

end subroutine diagnosis
