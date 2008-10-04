!=========================================================================================
  module param_m
!=========================================================================================
! module for param variables

  implicit none
!-----------------------------------------------------------------------------------------
! integers

  integer mode        !mode in which to run

  integer numdim      !number of dimensions of problem

  integer nx          !local number of grid points in x-direction
  integer ny          !local number of grid points in y-direction
  integer nz          !local number of grid points in z-direction

  integer nx_g        !global number of grid points in x-direction
  integer ny_g        !global number of grid points in y-direction
  integer nz_g        !global number of grid points in z-direction

  integer nxyz        !miscelaneous for compatibility checks
  integer nxyz2       !miscelaneous for compatibility checks
  integer nxm         !miscelaneous for compatibility checks
  integer nym         !miscelaneous for compatibility checks
  integer nzm         !miscelaneous for compatibility checks

  integer nxyz_g      !miscelaneous for compatibility checks
  integer nxyz2_g     !miscelaneous for compatibility checks
  integer nxm_g       !miscelaneous for compatibility checks
  integer nym_g       !miscelaneous for compatibility checks
  integer nzm_g       !miscelaneous for compatibility checks

  integer npx         !number of processors in x-direction
  integer npy         !number of processors in y-direction
  integer npz         !number of processors in z-direction

  integer nsc         !number of chemical species (excluding N2)

  integer nvar_tot    !total number of variables (equations) = 5 + nsc

  integer n_elem      !number of elements in chemical mechanism
  integer n_spec      !number of chemical species including N2
  integer n_reac      !total number of reactions

  integer ntr         !number of reactions with third-body efficiencies  

  integer n_reg       !number of registers used in runge-kutta integration

  integer periodic_x  !switch for periodicity in x-direction
  integer periodic_y  !switch for periodicity in y-direction
  integer periodic_z  !switch for periodicity in z-direction

  integer vary_in_x   !switch to turn on x-direction
  integer vary_in_y   !switch to turn on y-direction
  integer vary_in_z   !switch to turn on z-direction

  integer iorder      !order of derivatives
  integer iforder     !order of filter

  character*8  dat_1, dat_2     !for start and end date (wall clock)
  character*10 tim_1, tim_2     !for start and end time (wall clock)

!-----------------------------------------------------------------------------------------
  contains
!========================================================================================
  subroutine initialize_param(io,myid,ierr,gcomm)
!========================================================================================
! routine sets various parameters
!----------------------------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------------------------
  include 'mpif.h'
!----------------------------------------------------------------------------------------
! declarations passed in

  integer io
  integer myid, ierr, gcomm

! local declarations

!!$  integer nxnynz(3)
  integer i, iflag
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    write(io,*) 'initializing param module...'
    write(io,*)
  endif
!----------------------------------------------------------------------------------------
! set iflag

  iflag=0
!----------------------------------------------------------------------------------------
  if(myid.eq.0) then

!   error trapping for number of grid points in x-direction

    if(mod(nx_g,npx).ne.0 ) then

      if( myid.eq.0 ) then
        write(io,*) ' Grid Pts in X dimension ',  nx_g,  &
                    ' are not exactly divisible by ', npx ,  &
                    ' number of PEs in the X dimension '
      endif

      iflag=1

    endif

!   error trapping for number of grid points in x-direction

    if(mod(ny_g,npy).ne.0 ) then

      if(myid.eq.0) then
        write(io,*) ' Grid Pts in Y dimension ',  ny_g,  &
                    ' are not exactly divisible by ', npy ,  &
                    ' number of PEs in the Y dimension '
      endif

      iflag=1

    endif

!   error trapping for number of grid points in x-direction

    if(mod(nz_g,npz).ne.0 ) then

      if( myid.eq.0 ) then
        write(io,*) ' Grid Pts in Z dimension ',  nz_g,  &
                    ' are not exactly divisible by ', npz ,  &
                    ' number of PEs in the Z dimension '
      endif

      iflag=1

    endif

!   check for domain size conflict with stencils size in x-direction

    if((vary_in_x.eq.1).and.(nx_g.lt.9)) then

      if(myid.eq.0) then
        write(io,*) 'input error: nx_g < 9 and vary_in_x = 1'
      endif

      iflag=1

    endif

!   check for domain size conflict with stencils size in y-direction

    if((vary_in_y.eq.1).and.(ny_g.lt.9)) then

      if(myid.eq.0) then
        write(io,*) 'input error: ny_g < 9 and vary_in_y = 1'
      endif

      iflag=1

    endif

!   check for domain size conflict with stencils size in z-direction

    if((vary_in_z.eq.1).and.(nz_g.lt.9)) then

      if(myid.eq.0) then
        write(io,*) 'input error: nz_g < 9 and vary_in_z = 1'
      endif

      iflag=1

    endif

!   set local number of grid points

    nx = nx_g / npx
    ny = ny_g / npy
    nz = nz_g / npz

!   set some other stuff

    nxyz = max(nx, ny, nz);      nxyz_g = max(nx_g, ny_g, nz_g)
    nxyz2= 2*nxyz         ;      nxyz2_g= 2*nxyz_g      
    nxm  = max(1,nx - 1)  ;      nxm_g  = max(1,nx_g - 1)   
    nym  = max(1,ny - 1)  ;      nym_g  = max(1,ny_g - 1)   
    nzm  = max(1,nz - 1)  ;      nzm_g  = max(1,nz_g - 1)   

!   set chemistry parameters based on chemkin initialization

    call set_number_elem_spec_reac(n_elem,n_spec,n_reac,myid,io)
    call set_number_third_body_reactions(ntr,io)

    n_reac=n_reac*2   !convert n_reac for DNS/getrates purposes (2x)
    nsc=n_spec-1      !set number of chemical species for DNS purposes
    nvar_tot=nsc+5    !number of species + 5 (density, energy, momentum)

  endif
!----------------------------------------------------------------------------------------
! check status of error

  call MPI_Bcast(iflag,1,MPI_INTEGER,0,gcomm,ierr)

  if(iflag.eq.1) call terminate_run(io,0)  !must be called by all processors
!----------------------------------------------------------------------------------------
! broadcast parameters
!----------------------------------------------------------------------------------------
! broadcast local grid parameters

  call MPI_Bcast(nx       ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(ny       ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nz       ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nxyz     ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nxyz2    ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nxm      ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nym      ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nzm      ,1, MPI_INTEGER, 0, gcomm, ierr)

! broadcast global grid parameters

  call MPI_Bcast(nxyz_g  ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nxyz2_g ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nxm_g   ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nym_g   ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nzm_g   ,1, MPI_INTEGER, 0, gcomm, ierr)

! broadcast chemistry parameters

  call MPI_Bcast(nsc     ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nvar_tot,1, MPI_INTEGER, 0, gcomm, ierr)

  call MPI_Bcast(n_elem  ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(n_spec  ,1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(n_reac  ,1, MPI_INTEGER, 0, gcomm, ierr)

  call MPI_Bcast(ntr     ,1, MPI_INTEGER, 0, gcomm, ierr)

! broadcast runge-kutta parameters

  call MPI_Bcast(n_reg   ,1, MPI_INTEGER, 0, gcomm, ierr)

! sync processors

  call MPI_Barrier( gcomm,ierr )
!----------------------------------------------------------------------------------------
! miscellaneous initializations to be set on all processors
!----------------------------------------------------------------------------------------
! set number of dimensions

  numdim = 0
  if(vary_in_x.eq.1) numdim = numdim + 1
  if(vary_in_y.eq.1) numdim = numdim + 1
  if(vary_in_z.eq.1) numdim = numdim + 1

!!$  nxnynz(1)=nx
!!$  nxnynz(2)=ny
!!$  nxnynz(3)=nz
!!$
!!$! set number of waves (nw) and nh
!!$
!!$  if (numdim == 3) then
!!$    nw=int(min(nx_g,ny_g,nz_g)/2)-1
!!$  elseif(numdim == 2) then
!!$    if( vary_in_x==1 .and. vary_in_y==1 ) then
!!$      nw = int(min(nx_g,ny_g)/2)-1
!!$    elseif (vary_in_x==1 .and. vary_in_z==1) then
!!$      nw = int(min(nx_g,nz_g)/2)-1
!!$    elseif (vary_in_y==1 .and. vary_in_z==1) then
!!$      nw = int(min(ny_g,nz_g)/2)-1
!!$    endif
!!$  endif
!!$  nh = 1 + nw / 2

!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif
!----------------------------------------------------------------------------------------
  return
  end subroutine initialize_param
!----------------------------------------------------------------------------------------
  end module param_m
