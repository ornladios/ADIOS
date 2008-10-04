!========================================================================================
subroutine initialize_field(io)
!========================================================================================
! initializes primitive and solution variables amoung other things
!----------------------------------------------------------------------------------------
  use topology_m

  use runtime_m, only : run_title, time_save, time_save_inc, i_restart, time, i_time
  use runtime_m, only : tstep, time_restart

  use rk_m, only : tstep_init, cont_switch, tstep_min

  use param_m, only : nsc, n_elem, n_spec, mode, nx, ny, nz, nx_g, ny_g, nz_g
  use param_m, only : vary_in_x, vary_in_y, vary_in_z

  use reference_m

  !use turbulence_m

  use thermchem_m

  use variables_m

  use bc_m, only : store_hard_bc, nrf_x0, nrf_xl, nrf_y0, nrf_yl, nrf_z0, nrf_zl

  use grid_m, only : xmin, xmax, ymin, ymax, zmin, zmax, x

!  use feedTurb_m, only : setupTurbFeed, dumpVelField

  implicit none
!----------------------------------------------------------------------------------------
! declarations passed in

  integer, intent(in) :: io

! local declarations

  real :: a_l             !local speed of sound
  real :: a_max           !maximum global speed of sound
  real :: a_min           !minimum global speed of sound
  real :: Lxyz(3)         !domain lengths in each direction
  real :: L_max           !longest distance in domain
  real :: a_time_max      !time for max wave to travel the longest distance (sec)
  real :: a_time_min      !time for min wave to travel the longest distance (sec)
  integer :: dir, i1,i2,i3
  real :: L1,L2,L3
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    write(io,*) 'initializing variables...'
    write(io,*)
  endif
!----------------------------------------------------------------------------------------
! allocate arrays

  call allocate_variables_arrays(1)
!----------------------------------------------------------------------------------------
! initialize primative variables (yspecies, temp, pressure, and u only!)
! with user-specified routines.
!
! Note: Upon restart (i_restart=1), user-specified initialization routines should only 
!       initialize any necessary modules and not the primative variables. Any primative
!       variables initialized with i_restart=1 will be overwritten when the restart
!       files are read.
!----------------------------------------------------------------------------------------
  call MPI_Barrier(gcomm, ierr)

  i_time=0
  time=0.0
  time_save=time_save_inc
  time_restart=0.0

  select case (trim(run_title))

! add initialization cases here
!  case ('hcci')
!     call initialize_hcci(io,yspecies,temp,pressure,u)
!  case ('bomb')
!     call initialize_bomb(io,yspecies,temp,pressure,u)
!  case ('ignition_test')
!     call initialize_ignition_test(io,yspecies,temp,pressure,u)
!  case ('1d_flame_test')
!     call initialize_1d_flame_test(io,yspecies,temp,pressure,u)
  case ('pressure_wave_test')
     call initialize_pressure_wave_test(io,yspecies,temp,pressure,u)
!  case ('turbulent_air_test')
!     call initialize_turbulent_air_test(io,yspecies,temp,pressure,u)
!  case ('opthinrad_test')
!     call initialize_opthinrad_test(io,yspecies,temp,pressure,u)
!  case ('compheat_test')
!     call initialize_compheat_test(io,yspecies,temp,pressure,u)
  case default
     if(myid==0) then
        write(io,*) 'improper initialization of primative variables'
        write(io,*) 'a initialization routine does not exist for run_title: ',  &
             trim(run_title)
     endif
     call terminate_run(io,0)  !must be called by all processors
  end select
!----------------------------------------------------------------------------------------
! Restart code from prevous data files.
! Restarting will overwrite any initialization done in the initialization routines above!
!----------------------------------------------------------------------------------------

  if( (i_restart==1) .and. (mode==0) )  then

    if(myid == 0) write(io,*) 'restarting from saved files...'

!   read restart files

    call read_savefile(io)

!   set initial timestep for controller to a conservative fraction of
!   previous timestep

    if(cont_switch==1) tstep_init=(tstep*time_ref)/3.0  !factor of 3.0 is a conservative measure

!   set time_restart (seconds)

    time_restart=time*time_ref
    if(myid == 0) write(io,'(a16,1pe9.3,a6)') ' restart time = ',time_restart,' (sec)'

!   set save time (bug fix on restarts by evatt)
    if(time==0.0) then 
      time_save=time_save_inc
      tstep_init=max(tstep_min,tstep_init)
      tstep=max(tstep_min,tstep_init)
    end if

  endif
!----------------------------------------------------------------------------------------
! write header

  call MPI_Barrier(gcomm, ierr)

  if(myid.eq.0) then
    write(io,*) 'continuing with initialization...'
    write(io,*)
  endif
!----------------------------------------------------------------------------------------
! set derived variables from initialized primative variables

!  call set_variables(io)
!----------------------------------------------------------------------------------------
! store hard boundary conditions (volum must be set)

!  call store_hard_bc(io,q(1,1,1,1,1),temp,pressure,volum,i_restart)
!----------------------------------------------------------------------------------------
! sync processors

  call MPI_Barrier(gcomm,ierr)
!----------------------------------------------------------------------------------------
! initialize turbulence and/or vortex pairs

!  if(i_restart==0) then
!     select case (i_turbulence)
!     case(1)  ! initialize isotropic turbulence
!        call initialize_turbulence(io,q,xmin,xmax,ymin,ymax,zmin,zmax,x,re)
!     case(2)  ! initialize vortex pairs
!        call initialize_vortex( q, u, io )
!     end select
!  elseif(i_turbulence==1)then   
!    call get_turbulence_params(re,io,1)
!  endif


  call MPI_Barrier(gcomm,ierr)
!----------------------------------------------------------------------------------------
! store hard boundary conditions (volum must be set)

!  call store_hard_bc(io,q(1,1,1,1,1),temp,pressure,volum,i_restart)
!----------------------------------------------------------------------------------------

  call MPI_Barrier(gcomm,ierr)

! setup turbulence feed boundary conditions
! really only set up for feeding in x-direction...

!  if(i_turbulence==1) then
!     dir = 0
!     if(nrf_x0<=0 .or. nrf_xl<=0) then
!        dir = 1
!        i1 = nx;       i2 = ny;         i3 = nz;
!        L1 = xmax-xmin;  L2 = ymax-ymin;  L3 = zmax-zmin;
!     elseif(nrf_y0<=0 .or. nrf_yl<=0) then
!        dir=2
!        i1 = ny;       i2 = nx;         i3 = nz;
!        L1 = ymax-ymin;  L2 = xmax-xmin;  L3 = zmax-zmin;
!     elseif(nrf_z0<=0 .or. nrf_zl<=0) then
!        dir=3
!        i1 = nz;       i2 = nx;         i3 = ny;
!        L1 = zmax-zmin;  L2 = xmax-xmin;  L3 = ymax-ymin;
!     endif
!
!     if(dir/=0) then
!        if(myid==0) write(io,'(\A\)')'Setting up turbulent inflow conditions.'
!        call setupTurbFeed( io, i1,i2,i3, L1, L2, L3, dir )
!        call dumpVelField      ! dumps inflow fluctuations to an ASCII file.
!     endif
!  endif

!----------------------------------------------------------------------------------------
! sync processors

  call MPI_Barrier(gcomm,ierr)
!----------------------------------------------------------------------------------------
! update primative variables to be absolutely consistent with derived variables

!  call get_mass_frac(q(1,1,1,1,1),volum,yspecies)
!  call get_velocity_vec(u,q(1,1,1,1,1),volum)
!  call calc_temp(temp,q(:,:,:,5,1)*volum, u, yspecies )   !set T, Cp_mix
!  call calc_gamma( gamma, cpmix, mixMW )                  !set gamma
!  call calc_press( pressure, q(:,:,:,4,1), temp, mixMW )  !set pressure
!------------------------------------------------------------------------------
! calculate maximum speed of sound for given initial conditions

!  a_l=maxval(sqrt(gamma*pressure*volum))*a_ref
!  call MPI_Allreduce(a_l,a_max,1,MPI_REAL8,MPI_MAX,gcomm,ierr)

!  a_l=minval(sqrt(gamma*pressure*volum))*a_ref
!  call MPI_Allreduce(a_l,a_min,1,MPI_REAL8,MPI_MIN,gcomm,ierr)

!  if(myid.eq.0) then
!    write(io,100) 'maximum speed of sound = ',a_max, ' m/s'
!    write(io,100) 'mimimum speed of sound = ',a_min, ' m/s'
!    write(io,*)
!  endif

! set longest domain distance

!  Lxyz(:)=0.0
!  if(vary_in_x==1) Lxyz(1)=(xmax-xmin)*l_ref
!  if(vary_in_y==1) Lxyz(2)=(ymax-ymin)*l_ref
!  if(vary_in_z==1) Lxyz(3)=(zmax-zmin)*l_ref

!  L_max=sqrt(Lxyz(1)**2 + Lxyz(2)**2 + Lxyz(3)**2)

! calculate time for min and max acoustic waves to traverse domain
!
!  a_time_max=L_max/a_max
!  a_time_min=L_max/a_min

!  if((nx/=1).and.(ny/=1).and.(nz/=1)) then
!    if(myid.eq.0) then
!      write(io,101) 'maximum acoustic time =',a_time_max, ' sec'
!      write(io,101) 'minimum acoustic time =',a_time_min, ' sec'
!      write(io,*)
!    endif
!  endif
!----------------------------------------------------------------------------------------
! write closing header

  if(myid.eq.0) then
    write(io,102) 'initialization completed for time = ', time*time_ref
    call write_header(io,'-')
  endif
!----------------------------------------------------------------------------------------
! format statements

  100 format(1x,a,f7.2,a)
  101 format(1x,a,1pe9.2,a)
  102 format(1x,a,1pe9.2)
!----------------------------------------------------------------------------------------
  return
end subroutine initialize_field
