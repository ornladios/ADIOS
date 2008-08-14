Subroutine dataout3d
!! using adio
  use global_parameters
  use particle_decomp
  use field_array
  use particle_tracking

  integer i,j,k,kp,n,ij,istatus,ncid,dataid2(3),dimid(2),dataid1(2)
  real(wp) theta_start(mpsi+1,mzeta+1),radial(mpsi+1)
  real(wp) dataout(mgrid,mzeta+1)
  real(wp) X(mgrid,mzeta+1)
  real(wp) Y(mgrid,mzeta+1)
  real(wp) Z(mgrid,mzeta+1)
  real(wp) r,theta,zeta,theta0
  integer :: rank
  character(len=1) cdum(3)
  character(len=60),SAVE :: fdum
  character(len=100) output_fname
  character(len=100) dirstr 
  real(doubleprec) tt0,tt1,rtc
! HDF5 declarations
  integer :: ierr                 ! Return error flag
  real(wp):: curr_time
  integer :: merror               ! MPI error flag

#if ADIOS
  integer :: adios_err
  #define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0),adios_err)
  #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b,adios_err)
  integer*8 group_handle,type_id
  integer next,previous,current 
  integer*8 :: group_prof_handle 
#endif

!!first of all, write down the dimension information

   
! Since we have several MPI processes on each plane, only one of these will
! participate in writing out the data for that plane. Let's pick the processes
! with myrank_partd=0 on each plane.
if(myrank_partd==0)then
  if(myrank_toroidal==0)then
     previous=-1
  else
     previous=myrank_toroidal-1
  endif
  current = myrank_toroidal
  if(myrank_toroidal==(nproc_toroidal-1))then
    next=-1
  else
    next=myrank_toroidal+1
  endif
  if(istep==ndiag)then
     write(fdum,'("phi_dir/RUNcoord_",i5.5,".bp")')myrank_toroidal
     fdum=trim(fdum)//char(0)
!!!     fdum='phi_dir/RUNcoord.bp'//char(0)
!!! modified by zf2
!!!     write(dirstr,'("/Coordinate_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
     write(dirstr,'("/Coordinate")') !!!(mstepall+istep)/ndiag
     dirstr=trim(dirstr)//char(0)

!!! modified by zf2 for timing
     CALL MPI_BARRIER(toroidal_comm,merror)

     CALL open_start_for_group(group_prof_handle, "output3d.0"//char(0),istep)

     !!!call adios_get_group(type_id,"output3d.0"//char(0))
     call adios_open(group_handle,"output3d.0"//char(0),fdum,"w"//char(0),adios_err)
     call adios_set_path (group_handle,dirstr,adios_err)

    CALL open_end_for_group(group_prof_handle,istep)

  !! call adios_group_by(type_id,"mype"//char(0),toroidal_comm,previous,current,next)
   !!  if(myrank_toroidal==0)then
        do i=0,mpsi
           radial(i+1)=(a0+deltar*real(i))/gyroradius
        enddo
   !!  endif
!  ! Mapping of magnetic coordinates onto cartesian grid and write into file
   ! We explicitely add the zeta=0 poloidal plane since it will be used by
   ! the vizualization software to connect the torus. myrank_toroidal=0 takes care of
   ! this, having one extra plane to write into the output file (kp=1). 
     cdum(1)='X'
     cdum(2)='Y'
     cdum(3)='Z'
    !! if(myrank_toroidal.eq.0)then
    !!   kp=1   !add the plane zeta=0 to myrank_toroidal=0 (zetamin=0 for myrank_toroidal=0)
    !! else
    !!   kp=0   !no changes for myrank_toroidal > 0
    !! endif
     kp=0
     !do n=1,3
        do k=1,mzeta+kp
           zeta=zetamin+deltaz*real(k-kp)
           do i=0,mpsi
              r=a0+deltar*real(i)
              theta_start(i+1,k)=zeta*qtinv(i)
              do j=0,mtheta(i)
                 ij=igrid(i)+j
                 theta=deltat(i)*real(j)+zeta*qtinv(i)
                !!! theta0=theta+r*sin(theta)
                 theta0=theta
                 !if(n==1)then
                    X(ij,k)=cos(zeta)*(1.0+r*cos(theta0))
                 !elseif(n==2)then
                    Y(ij,k)=-sin(zeta)*(1.0+r*cos(theta0))
                 !else
                    Z(ij,k)=r*sin(theta0)
                 !endif
              enddo
           enddo
        enddo
!      If n=1, write out the theta_start(mpsi+1,mzeta) array
       !if(n.eq.1)then
         !!dimsf(1)=mpsi+1     ! First dimension of theta_start
         !!dimsf(2)=mzeta*ntoroidal+1  ! = mzetamax+1 (for the extra plane)
       !endif
!!write dimension

    CALL write_start_for_group(group_prof_handle,istep)

        ADIOS_WRITE(group_handle,toroidal_comm)
        ADIOS_WRITE(group_handle,mpsi+1)
        ADIOS_WRITE(group_handle,mzeta+kp)
!!write data
        ADIOS_WRITE(group_handle,myrank_toroidal)
        ADIOS_WRITE(group_handle,mpsi)
        ADIOS_WRITE(group_handle,mzeta)
        ADIOS_WRITE(group_handle,kp)
        ADIOS_WRITE(group_handle,nproc_toroidal)

!!! added by zf2
        if(myrank_toroidal==0) then 
          ADIOS_WRITE(group_handle,radial)
          call adios_write(group_handle,"mtheta"//char(0),(mtheta+1),adios_err)
          ADIOS_WRITE(group_handle,itran)
        endif
      !!!   call adios_write(group_handle,"mpsi_dim"//char(0),mpsi)

!!! added by zf2
!!!     ADIOS_WRITE(group_handle,(mzeta+kp)*nproc_toroidal)
!!!     ADIOS_WRITE(group_handle,(mzeta+kp)*myrank_toroidal)

      !!!   ADIOS_WRITE(group_handle,theta_start)

        ADIOS_WRITE(group_handle,mgrid)

!!! added by zf2
        ADIOS_WRITE(group_handle,mgrid*nproc_toroidal)
        ADIOS_WRITE(group_handle,mgrid*myrank_toroidal)
    
        ADIOS_WRITE(group_handle,X)
        ADIOS_WRITE(group_handle,Y)
        ADIOS_WRITE(group_handle,Z)

    CALL write_end_for_group(group_prof_handle,istep)

    CALL close_start_for_group(group_prof_handle,istep)

       !! ADIOS_WRITE(group_handle,mype)
        call adios_close(group_handle,adios_err)

    CALL close_end_for_group(group_prof_handle,istep)

    !! enddo  !! loop of n
!!    Close the file.
     !!debug output
       !! if(myrank_toroidal==0)then
        !!   open(127,file='newdiag.out',status='replace')
        !!   write(127,*)'mpsi= ',mpsi
        !!   write(127,*)'mzeta= ',mzeta
        !!   write(127,*)'kp= ',kp
        !!   write(127,*)'nproc_toroidal= ',nproc_toroidal
        !!   write(127,*)'mgrid= ',mgrid
        !!   close(127)
        !! endif
  endif !! end if if(istep==ndiag)

! potential data file
!! if(myrank_toroidal.eq.0)then
!!    kp=1   !add the plane zeta=0 to myrank_toroidal=0 (zetamin=0 for myrank_toroidal=0)
!!  else
!!    kp=0   !no changes for myrank_toroidal > 0
!!  endif
  kp=0
  do k=1,mzeta+kp
     do i=0,mpsi
        do j=0,mtheta(i)
           ij=igrid(i)+j
           dataout(ij,k)=phi(k-kp,ij)
        enddo
     enddo
  enddo

!!!  write(output_fname,'("phi_dir/PHI_",i5.5,"_",i5.5,".bp")')mstepall+istep,myrank_toroidal !!(1+(mstepall+istep-ndiag)/isnap)
  write(output_fname,'("phi_dir/PHI_",i5.5,".bp")')mstepall+istep !!(1+(mstepall+istep-ndiag)/isnap)
  output_fname=trim(output_fname)//char(0)
  !!dirstr="Potential"//char(0)

!!! modified by zf2
!!!  write(dirstr,'("/Potential_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
!!!  dirstr=trim(dirstr)//char(0)
  dirstr="Potential"//char(0)

!!!  write(dirstr,'("/Potential_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
!!!  dirstr=trim(dirstr)//char(0)
  !!!if(mype==0)write(stdout,*)"DIRSTR: ",dirstr

!!! modified by zf2 for timing
  CALL MPI_BARRIER(toroidal_comm,merror)

  CALL open_start_for_group(group_prof_handle, "output3d.1"//char(0),istep)

  !!!call adios_get_group (type_id, "output3d.1"//char(0))
!!  call adios_group_by (type_id,"mype"//char(0), toroidal_comm, previous, current, next)
  call adios_open (group_handle, "output3d.1"//char(0), output_fname,"w"//char(0),adios_err)
  call adios_set_path (group_handle,dirstr,adios_err)

  CALL open_end_for_group(group_prof_handle,istep)

!  CALL write_start_for_group(group_prof_handle,istep)

!!  call adios_write_offset (group_handle, phi_offset)
  ADIOS_WRITE(group_handle,toroidal_comm)
  ADIOS_WRITE(group_handle,mzeta+kp)  ! TODO check this for correctness

  ADIOS_WRITE(group_handle,myrank_toroidal)
  ADIOS_WRITE(group_handle,nproc_toroidal)
  ADIOS_WRITE(group_handle,mgrid)
  ADIOS_WRITE(group_handle,mzeta)
  ADIOS_WRITE(group_handle,kp)
  ADIOS_WRITE(group_handle,mpsi)
  ADIOS_WRITE(group_handle,mpsi+1)
  !!ADIOS_WRITE(group_handle,mype)

!!! added by zf2
  ADIOS_WRITE(group_handle,mgrid*nproc_toroidal)
  ADIOS_WRITE(group_handle,mgrid*myrank_toroidal)

  if(myrank_toroidal==0) then
    call adios_write(group_handle,"mtheta"//char(0),(mtheta+1),adios_err)
  endif
  call adios_write (group_handle, "phi"//char(0), dataout,adios_err)

  CALL write_end_for_group(group_prof_handle,istep)

  CALL close_start_for_group(group_prof_handle,istep)

  call adios_close (group_handle)

  CALL close_end_for_group(group_prof_handle,istep)

endif !!end myrank_partd   

! if(myrank_partd==0)then
! The processes not participating in writing to file wait at the following
! barrier. We use the partd_comm communicator since it does not require a
! fully global synchronization.
call MPI_BARRIER(partd_comm,ierr)

end

