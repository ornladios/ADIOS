!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine outputs the potential data for a single poloidal plane.
! We pick the first plane held by mype=0.
!
! NOTE: THIS IS A VERSION OF OUTPUT3D THAT WRITES OUT A POLOIDAL
!       GRID OF CONSTANT DTHETA INSTEAD OF THE REAL GRID WHICH
!       IS OF CONSTANT R*DTHETA. LINEAR INTERPOLATION IS USED TO
!       GET THE EXTRA POINTS BETWEEN 2 REAL GRID POINTS.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine output2d
  use global_parameters
  use field_array
  use particle_decomp
  use particle_tracking
  implicit none

  integer i,j,k,kp,n,ij
  real(wp) dataout(mthetamax+1,mpsi+1),radial(mpsi+1)
  real(wp) r,theta,zeta,theta0,dtheta
  integer :: rank
  integer, dimension(2) :: nval,dimsf
  integer, dimension(2) :: offset
  character(len=2) cdum(3)
  character(len=20) vdum
  character(len=20),SAVE :: fdum
  character(len=60) :: dataname

  real(doubleprec) tt0,tt1,rtc

  integer :: ierr                 ! Return error flag
  integer :: info,iunit
  real(wp):: curr_time

  info = MPI_INFO_NULL

! Write out the grid and its cartesian mapping only once at the beginning.
  if((mstepall+istep)==isnap .and. mype==0)then

   ! Angular step between to points on a flux surface
     dtheta=2.0_wp*pi/real(mthetamax)

   ! Filename for grid quantities
     fdum='XYcoord.h5'//char(0)

   ! Create the file collectively.
     iunit=1
     call open_hdf5_file(iunit,filename)

   ! Write datasets to file. Only the master process writes the datasets.
     do i=0,mpsi
        radial(i+1)=(a0+deltar*real(i))/gyroradius
     enddo

   ! Write scalars first
     rank=0   ! scalar has rank=0
     dimsf(1)=0
   ! Write "mpsi+1"
     dataname='number_of_flux_surfaces'//char(0)
     call write_hdf5_f(iunit,(mpsi+1),'INT',rank,dims,trim(dataname))

   ! Write the number of poloidal angles on each flux surface (mthetamax+1)
     dataname='number_of_poloidal_angles'//char(0)
     call write_hdf5_f(iunit,(mthetamax+1),'INT',rank,dimsf,trim(dataname))

   ! Write 1D arrays
     rank=1
     dimsf(1)=mpsi+1

     dataname='radial_grid'//char(0)
     call write_hdf5_f(iunit,radial,'FLOAT',rank,dimsf,trim(dataname))

!  ! Mapping of magnetic coordinates onto cartesian grid and write into file
   ! We explicitely add the zeta=0 poloidal plane since it will be used by
   ! the vizualization software to connect the torus. mype=0 takes care of
   ! this, having one extra plane to write into the output file (kp=1). 
     cdum(1)='X'//char(0)
     cdum(2)='Y'//char(0)
     cdum(3)='Z'//char(0)
     do n=1,2
        k=0
        zeta=zetamin+deltaz*real(k)
           do i=0,mpsi
              r=a0+deltar*real(i)
              do j=0,mthetamax
                 theta0=dtheta*real(j)
                 if(n==1)then  ! X
                    dataout(j,i)=r*cos(theta0)
                 elseif(n==2)then ! Y
                    dataout(j,i)=r*sin(theta0)
                 endif
              enddo
           enddo
        enddo

     ! Create the data space for the 2D dataset.
       rank=2
       dimsf(1)=mthetamax+1    ! First dimension of dataout
       dimsf(2)=mpsi+1         ! Second dimension (flux surfaces)

     ! Determine size and offset of dataset for each process writing to file
       nval(1)=mgrid     ! First dimension of dataout
       offset(1)=0
       if(myrank_toroidal==0)then
         nval(2)=mzeta+1  !2nd dim of dataout (extra plane myrank_toroidal=0)
         offset(2)=0
       else
         nval(2)=mzeta    !Second dimension of dataout (=mzetamax/ntoroidal)
         offset(2)=myrank_toroidal*mzeta+1  !shifted by +1 because of extra plane myrank_toroidal=0
       endif

       tt0=rtc()
     ! Write the dataset in parallel in HDF5 file.
       call write_parhdf5_f(iunit,dataout,'FLOAT',rank,dimsf,nval,&
                            offset,toroidal_comm,cdum(n))
       tt1=rtc()

     !!  write(*,*)mype,'  time for collective HDF5 write =',(tt1-tt0)
     enddo

   ! Close the file.
     call close_hdf5_file(iunit)
  endif

! potential data file
  if(myrank_toroidal==0)then
    kp=1   !add the plane zeta=0 to mype=0 (zetamin=0 for mype=0)
  else
    kp=0   !no changes for mype > 0
  endif
  do k=1,mzeta+kp
     do i=0,mpsi
        do j=0,mtheta(i)
           ij=igrid(i)+j
           dataout(ij,k)=phi(k-kp,ij)
        enddo
     enddo
  enddo

  if(myrank_partd==0)then
   ! Filename for grid quantities
     write(fdum,'("PHI_",i0,".h5",a1)')(1+(mstepall+istep-ndiag)/isnap),char(0)

   ! Create the file collectively.
     iunit=1
     call open_parhdf5_file(iunit,toroidal_comm,fdum)

   ! We now write the current calculation time as an attribute to the
   ! root group. We first need to get the root group id and then create the
   ! attribute data space, data type, etc. Only one processor needs to write
   ! the attribute.

   ! Create the data space for the dataset.
     rank=2
     dimsf(1)=mgrid     ! First dimension of dataout
     dimsf(2)=mzeta*ntoroidal+1  ! = mzetamax

   ! Write dataset name into variable
     write(vdum,'("Potential_t",i0,a1)')(mstepall+istep)/ndiag,char(0)

   ! Each process defines dataset in memory and writes it to the hyperslab
   ! in the file.
     nval(1)=mgrid     ! First dimension of dataout
     offset(1)=0
     if(myrank_toroidal==0)then
       nval(2)=mzeta+1  !2nd dim of dataout (extra plane for myrank_toroidal=0)
       offset(2)=0
     else
       nval(2)=mzeta    !Second dimension of dataout (=mzetamax/numberpe)
       offset(2)=myrank_toroidal*mzeta+1  !shifted by +1 because of extra plane myrank_toroidal=0
     endif

     tt0=rtc()
   ! Write the dataset
     call write_parhdf5_f(iunit,dataout,'FLOAT',rank,dimsf,nval,&
                          offset,toroidal_comm,trim(vdum))
     tt1=rtc()

   ! Close the file.
     call close_hdf5_file(iunit)

   endif
end

