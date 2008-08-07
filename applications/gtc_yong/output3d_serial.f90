Subroutine output3d
  use global_parameters
  use field_array
  use particle_decomp
  use particle_tracking
!  use h5global
  use hdf5

  integer i,j,k,kp,n,ij,istatus,ncid,dataid2(3),dimid(2),dataid1(2)
  real(wp) dataout(mgrid,mzeta),radial(mpsi+1)
  real(wp) r,theta,zeta,theta0
  integer :: rank,intval
  character(len=1) cdum(3)
  character(len=20) vdum
  character(len=20),SAVE :: fdum

!  real(doubleprec) tt0,tt1,rtc

! HDF5 declarations
  integer(HSIZE_T), dimension(7) :: dimsfi = (/0,0,0,0,0,0,0/)
  integer(HSIZE_T), dimension(2) :: count,dimsf
  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in memory
  integer(HID_T) :: plist_id      ! Property list identifier
  integer(HID_T) :: grp_id        ! Group identifier 
  integer(HID_T) :: time_id       ! Attribute identifier 
  integer(HID_T) :: mpsi_id       ! Attribute identifier 
  integer(HID_T) :: npe_id        ! Attribute identifier 
  integer(HID_T) :: mtheta_id     ! Attribute identifier 
  integer(HID_T) :: radial_id     ! Attribute identifier 
  integer(HID_T) :: itran_id      ! Attribute identifier 
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier 
  integer(HID_T) :: a1Dspace_id   ! Attribute Dataspace identifier 
  integer :: ierr                 ! Return error flag
  integer :: comm,info
  real(wp):: curr_time

  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

! Initializes the HDF5 library and the Fortran90 interface
  call h5open_f(ierr)

! Write out the grid and its cartesian mapping only once at the beginning.
  if(istep==isnap .and. myrank_partd==0)then

   ! Filename for grid quantities
     !fdum='RUNcoord.h5'
     write(fdum,'("RUNcoord_",i0,".h5")')myrank_toroidal

   ! Create the file.
     call h5fcreate_f(fdum,H5F_ACC_TRUNC_F,file_id,ierr)

   ! Create scalar data space for scalar attributes. 
     call h5screate_f(H5S_SCALAR_F,aspace_id,ierr)

   ! Create 1D data space for 1D attributes. 
     rank=1
     dimsf(1)=mpsi+1   ! Dimension of 1D attributes
     call h5screate_simple_f(rank,dimsf,a1Dspace_id,ierr)

   ! Create group datasets.
     call h5dcreate_f(file_id,"number_of_flux_surfaces",H5T_NATIVE_INTEGER,&
                      aspace_id,mpsi_id,ierr)
     call h5dcreate_f(file_id,"number_of_PE",H5T_NATIVE_INTEGER,aspace_id,&
                      npe_id,ierr)
     call h5dcreate_f(file_id,"radial_grid",H5T_NATIVE_REAL,a1Dspace_id,&
                      radial_id,ierr)
     call h5dcreate_f(file_id,"poloidal_grid",H5T_NATIVE_INTEGER,a1Dspace_id,&
                      mtheta_id,ierr)
     call h5dcreate_f(file_id,"index_shift",H5T_NATIVE_INTEGER,a1Dspace_id,&
                      itran_id,ierr)
     
   ! Write datasets to file. Only the master process writes the datasets.
     !!!if(mype==0)then
        do i=0,mpsi
           radial(i+1)=(a0+deltar*real(i))/gyroradius
        enddo

        dimsfi(1)=0
        intval=mpsi+1
        call h5dwrite_f(mpsi_id,H5T_NATIVE_INTEGER,intval,dimsfi,ierr)
        call h5dwrite_f(npe_id,H5T_NATIVE_INTEGER,numberpe,dimsfi,ierr)
        dimsfi(1)=dimsf(1)
        call h5dwrite_f(radial_id,H5T_NATIVE_REAL,radial,dimsfi,ierr)
        call h5dwrite_f(mtheta_id,H5T_NATIVE_INTEGER,(mtheta+1),dimsfi,ierr)
        call h5dwrite_f(itran_id,H5T_NATIVE_INTEGER,itran,dimsfi,ierr)
     !!!endif

   ! Release datasets ids, datasets data spaces, and group id.
     call h5dclose_f(mpsi_id,ierr)
     call h5dclose_f(npe_id,ierr)
     call h5dclose_f(mtheta_id,ierr)
     call h5dclose_f(radial_id,ierr)
     call h5dclose_f(itran_id,ierr)
     call h5sclose_f(aspace_id,ierr)
     call h5sclose_f(a1Dspace_id,ierr)

!  ! Mapping of magnetic coordinates onto cartesian grid and write into file
   ! We explicitely add the zeta=0 poloidal plane since it will be used by
   ! the vizualization software to connect the torus. mype=0 takes care of
   ! this, having one extra plane to write into the output file (kp=1). 
     cdum(1)='X'
     cdum(2)='Y'
     cdum(3)='Z'
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

!      Create the data space for the 2D dataset.
       rank=2
       dimsf(1)=mgrid     ! First dimension of dataout
       dimsf(2)=mzeta+kp  ! = mzetamax+1 (for the extra plane)
       dimsfi(1)=dimsf(1)
       dimsfi(2)=dimsf(2)
       call h5screate_simple_f(rank,dimsf,filespace,ierr)

!      Create the dataset with default properties.
       call h5dcreate_f(file_id,cdum(n),H5T_NATIVE_REAL,filespace,&
                        dset_id,ierr)

       !tt0=rtc()
!      Write the dataset collectively.
       call h5dwrite_f(dset_id,H5T_NATIVE_REAL,dataout, dimsfi,ierr)
       !tt1=rtc()

!      Close dataspaces, the dataset and property list.
       call h5dclose_f(dset_id,ierr)
       call h5sclose_f(filespace,ierr)

     !!  write(*,*)mype,'  time for collective HDF5 write =',(tt1-tt0)

     enddo

!    Close the file.
     call h5fclose_f(file_id,ierr)
  endif

! potential data file
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

  !!!!!if(mod((istep-ndiag),isnap)==0)then
  if(myrank_partd==0)then
   ! We open a new hdf5 file at the first diagnostic call after an output
   ! to the snapshot files.

   ! Filename for grid quantities
     write(fdum,'("PHI_",i0,"_",i0,".h5")')(mstepall+istep),myrank_toroidal

   ! Create the file.
     call h5fcreate_f(fdum,H5F_ACC_TRUNC_F,file_id,ierr)

! We now write the current calculation time as an attribute to the
! root group. We first need to get the root group id and then create the
! attribute data space, data type, etc. Only one processor needs to write
! the attribute.

   ! Create the data space for the dataset.
     rank=2
     dimsf(1)=mgrid     ! First dimension of dataout
     dimsf(2)=mzeta+kp  ! = mzetamax
     dimsfi(1)=dimsf(1)
     dimsfi(2)=dimsf(2)
     call h5screate_simple_f(rank,dimsf,filespace,ierr)

! Write dataset name into variable
     write(vdum,'("Potential_t",i0)')(mstepall+istep)

! Create the dataset with default properties.
     call h5dcreate_f(file_id,vdum,H5T_NATIVE_REAL,filespace,&
                   dset_id,ierr)

! Create scalar data space for the attribute. 
     call h5screate_f(H5S_SCALAR_F,aspace_id,ierr)

! Create dataset attribute.
     call h5acreate_f(dset_id,"time",H5T_NATIVE_REAL,aspace_id,&
                   time_id,ierr,H5P_DEFAULT_F)
     
     curr_time=real(mstepall+istep)*tstep  ! Attribute data

! Write the time dataset data. Here, dimsfi is just a dummy argument (ignored).
       call h5awrite_f(time_id,H5T_NATIVE_REAL,curr_time,dimsfi,ierr)
     
! Release attribute id, attribute data space, and group id.
     call h5aclose_f(time_id,ierr)
     call h5sclose_f(aspace_id,ierr)

     !tt0=rtc()
! Write the dataset
     call h5dwrite_f(dset_id,H5T_NATIVE_REAL,dataout, dimsfi,ierr)
     !tt1=rtc()

! Close dataspaces, the dataset and property list.
     call h5dclose_f(dset_id,ierr)
     call h5sclose_f(filespace,ierr)

! Close the file.
     call h5fclose_f(file_id,ierr)

! Close the HDF5 library and the Fortran90 interface
     call h5close_f(ierr)
  !!   write(*,*)mype,'  time for collective HDF5 Potential write =',(tt1-tt0)
  endif

end

