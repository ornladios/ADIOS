Subroutine output3d
  use global_parameters
  use field_array
  use particle_decomp
  use particle_tracking

  integer i,j,k,kp,n,ij
  real(wp) dataout(mgrid,mzeta),radial(0:mpsi)
  real(wp) r,theta,zeta,theta0
  character(len=1) cdum(3)
  character(len=20) vdum
  character(len=20) fdum

!  real(doubleprec) tt0,tt1,rtc

  real(wp):: curr_time


! Write out the grid and its cartesian mapping only once at the beginning.
  if((mstepall+istep)==isnap .and. myrank_partd==0)then

     do i=0,mpsi
        radial(i)=(a0+deltar*real(i))/gyroradius
     enddo

   ! Filename for grid quantities
     !fdum='RUNcoord.h5'
     write(fdum,'("RUNcoord_",i0,".dat")')myrank_toroidal

     open(888,file=fdum,status='replace',form='unformatted')

     write(888)mpsi,mgrid,mzeta
     write(888)mtheta(0:mspi),igrid(0:mpsi),qtinv(0:mpsi),radial(0:mpsi)

!  ! Mapping of magnetic coordinates onto cartesian grid and write into file
   ! We explicitely add the zeta=0 poloidal plane since it will be used by
   ! the vizualization software to connect the torus. mype=0 takes care of
   ! this, having one extra plane to write into the output file (kp=1). 
     cdum(1)='X'
     cdum(2)='Y'
     cdum(3)='Z'
     dataout=0.0
     if(myrank_toroidal.eq.0)then
       !!!!kp=1   !add the plane zeta=0 to mype=0 (zetamin=0 for mype=0)
       kp=0  !Don't need to add copy of the plane for serial version
     else
       kp=0   !no changes for mype > 0
     endif
     do n=1,3
        do k=1,mzeta+kp
           zeta=zetamin+deltaz*real(k-kp)
           do i=0,mpsi
              r=a0+deltar*real(i)
              !!!!do j=0,mtheta(i)
              do j=1,mtheta(i) !Remove duplicated point at 0
                 ij=igrid(i)+j
                 theta=deltat(i)*real(j)+zeta*qtinv(i)
                 theta0=theta+r*sin(theta)
! grid coordinates (x,y,z), use geometry center as the origin 
                 if(n==1)then
                    dataout(ij,k)=cos(zeta)*(1.0+r*cos(theta0))
                 elseif(n==2)then
                    dataout(ij,k)=-sin(zeta)*(1.0+r*cos(theta0))
                 else
                    dataout(ij,k)=r*sin(theta0)
                 endif
              enddo
           enddo
        enddo

      ! Write X (n=1), Y (n=2), and Z (n=3)
        write(888)dataout

     enddo

!    Close the file.
     close(888)
  endif

! potential data file
  dataout=0.0
  if(myrank_toroidal.eq.0)then
    !!!!!kp=1   !add the plane zeta=0 to mype=0 (zetamin=0 for mype=0)
    kp=0  !Don't need to add copy of the plane for serial version
  else
    kp=0   !no changes for mype > 0
  endif
  do k=1,mzeta+kp
     do i=0,mpsi
        !!!!do j=0,mtheta(i)
        do j=1,mtheta(i) !Remove duplicated point at 0
           ij=igrid(i)+j
           dataout(ij,k)=phi(k-kp,ij)
        enddo
     enddo
  enddo

  if(myrank_partd==0)then

   ! Filename for grid quantities
     write(fdum,'("PHI_",i0,"_",i0,".dat")')(mstepall+istep),myrank_toroidal

   ! Create the file.
     open(889,file=fdum,status='replace',form='unformatted')

   ! Write current time step and simulation time
     curr_time=real(mstepall+istep)*tstep  ! Attribute data
     write(889)(mstepall+istep),curr_time

   ! Write the dataset
     write(889)dataout

   ! Close the file.
     close(889)

  endif

end

