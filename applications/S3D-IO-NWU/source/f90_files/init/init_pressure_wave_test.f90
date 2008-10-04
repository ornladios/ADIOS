!========================================================================================
  subroutine initialize_pressure_wave_test( io, yspecies, temp, pressure, u )
!========================================================================================
! This subroutine initializes the following quantities:                        !
!                                                                              !
!   yspecies - Species mass fractions                                          !
!   temp     - Temperature                                                     !
!   pressure - pressure field                                                  !
!   u        - Mean velocity field                                             !
!                                                                              !
! It is currently set up to do a gaussian temperature spike in an otherwise    ! 
! quiescent field to initialize a pressure wave.                               !
!----------------------------------------------------------------------------------------
  use param_m
  use topology_m
  use grid_m, only : x, y, z, xmin, xmax, ymin, ymax, zmin, zmax
  use reference_m
  use chemkin_m, only : species_name

!  use feedTurb_m, only : vel_bar, meanVel
  use work_m, only : gauss => work1_1

  implicit none
!----------------------------------------------------------------------------------------
! declarations passed in

  real, intent(inout), dimension(nx,ny,nz) :: temp, pressure
  real, intent(inout), dimension(nx,ny,nz,n_spec) :: yspecies
  real, intent(inout), dimension(nx,ny,nz,3) :: u

  integer, intent(in) :: io

! local declarations

  integer :: i,j,k,L

  real :: maxTemp = 1000,       ambientTemp = 300
  real :: variance

  real :: maxPt(3)
  real :: junk, facx, facy, facz, max, max_g

  real, parameter :: pi = 3.14159265358979
  character*2 ext

!----------------------------------------------------------------------------------------
! write header

  if(myid == 0) then
    write(io,*) 'initializing pressure wave test...'
    write(io,*)
  endif
!-------------------------------------------------------------------------------------
! set pressure

  pressure = pres_atm / p_ref  !set to one atmosphere
!-------------------------------------------------------------------------------------
! set species for air

  yspecies = 0.0
  do L=1,n_spec,1
    if(trim(species_name(L)).eq.'O2') yspecies(:,:,:,L)=0.233
    if(trim(species_name(L)).eq.'N2') yspecies(:,:,:,L)=0.767
  enddo
!-------------------------------------------------------------------------------------
! set mean velocities

  u(:,:,:,1) = 0.10               ! mean u-velocity set to zero
  u(:,:,:,2) = 0.0                ! mean v-velocity set to zero
  u(:,:,:,3) = 0.0                ! mean w-velocity set to zero

!-------------------------------------------------------------------------------------
! set mean convective velocity for turbulence feeding. See "feedTurb_m.f90"
! these variables MUST be set at RESTARTS

!  meanVel = maxval(u(:,:,:,1))          ! scanning velocity
!  allocate(vel_bar(ny,nz))
!  vel_bar = u(1,:,:,1)                  ! mean velocity at inlet plane
  ! ideally, vel_bar is uniform and meanVel==vel_bar

!-------------------------------------------------------------------------------------
! set temperature

  ambientTemp = ambientTemp / t_ref
  maxTemp = maxTemp / t_ref

! set center of spike to center of domain

  maxPt = 0.0
  if (vary_in_x == 1)   maxPt(1) = (xmax-xmin) / 2.
  if (vary_in_y == 1)   maxPt(2) = (ymax-ymin) / 2.
  if (vary_in_z == 1)   maxPt(3) = (zmax-zmin) / 2.

! set the variance

  variance = 0.0
  if (vary_in_x == 1)   variance = variance + (xmax-maxPt(1))**2
  if (vary_in_y == 1)   variance = variance + (ymax-maxPt(2))**2
  if (vary_in_z == 1)   variance = variance + (zmax-maxPt(3))**2
  variance = 0.01*variance

! compute temperature using a gaussian distribution centered at the middle of the domain

  facx = 0.0
  facy = 0.0
  facz = 0.0

  do k=1,nz
     do j=1,ny
        do i=1,nx

           if (vary_in_x == 1)   facx = x(i)-maxPt(1)
           if (vary_in_y == 1)   facy = y(j)-maxPt(2)
           if (vary_in_z == 1)   facz = z(k)-maxPt(3)

           junk = sqrt( facx**2 + facy**2 + facz**2 )

           gauss(i,j,k) = 1/sqrt(2.0*pi*variance) * exp(-junk**2 / (2.0*variance))

        enddo
     enddo
  enddo

  max = maxval(gauss)

  call MPI_Allreduce(max,max_g,1,MPI_REAL8,MPI_MAX,gcomm,ierr)

  gauss = gauss / max_g

  temp = ambientTemp + (maxTemp-ambientTemp) * gauss
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif

!-------------------------------------------------------------------------------------
  end subroutine initialize_pressure_wave_test
