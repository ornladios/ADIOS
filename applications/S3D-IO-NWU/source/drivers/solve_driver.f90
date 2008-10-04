!$Id: solve_driver.f90,v 1.5.2.2 2006/04/04 18:25:21 rsankar Exp $
!========================================================================================
  subroutine solve_driver(io)
!========================================================================================
! routine drives the solution of the governing equations
!----------------------------------------------------------------------------------------
  use topology_m, only : gcomm, ierr, myid

  use param_m, only : dat_1

  use reference_m, only : time_ref
  use reference_m, only : initialize_reference  !routine reference

  use runtime_m, only : i_restart, i_time_res, i_time, i_time_end
  use runtime_m, only : i_time_mon, i_time_tec, i_time_save
  use runtime_m, only : time, tstep, time_save, time_save_inc

  !use filter_m, only : initialize_filter, allocate_filter_arrays  !routine references

  use rk_m, only : q_err_max
  use rk_m, only : i_time_cont
!  use rk_m, only : initialize_rk, allocate_rk_arrays  !routine references

  use thermchem_m, only : initialize_thermchem        !routine reference
  use thermchem_m, only : allocate_thermchem_arrays   !routine reference

  use variables_m, only : temp, pressure
  use variables_m, only : allocate_variables_arrays  !routine reference

  use work_m, only : allocate_work_arrays  !routine reference

  use grid_m, only : initialize_grid  !routine reference
  use grid_m, only : allocate_grid_arrays  !routine reference

  use bc_m, only : initialize_bc  !routine reference
  use bc_m, only : allocate_bc_arrays  !routine reference

!  use transport_m, only : pr

!  use transport_m, only : initialize_transport  !routine reference
!  use transport_m, only : allocate_transport_arrays  !routine reference

  use derivative_m, only : initialize_derivative  !routine reference
  use derivative_m, only : allocate_derivative_arrays  !routine reference

  ! Added by Ramanan 05/04/2005.
  ! Tabulate dG/T reaction-wise.
!  use gibbs_rxn_table_m, only: init_gibbstable

! wkliao: read/write files through MPI I/O
  use mpi_io_m
  use adios_m

  implicit none
!----------------------------------------------------------------------------------------
! declarations of variables passed in

  integer, intent(in) :: io

! declarations of local variables

  integer :: L, i, j
  real :: temp_max, pres_max  !maximum temperature and pressure in domain
  character*10 :: tim_start, tim_stop     !for start and end time (wall clock)
  integer :: tec_del
!----------------------------------------------------------------------------------------
! set variable tec_del (0=save all tecplot files, 1=delete previous file upon new write)

  tec_del=1
!----------------------------------------------------------------------------------------
!-- allocate work space
  call allocate_work_arrays(1)

!-- initialize reference module
!  call initialize_reference(io,pr)
  call initialize_reference(io,0.708)

!-- initialize rk module
!  call initialize_rk(io)

!-- initialize transport module
!  call initialize_transport(io)

!-- initialize filter module
!  call initialize_filter(io)

!-- initialize derivative module
  call initialize_derivative(io)

!-- initialize thermo and chemistry module (must be done after reference module)
  call initialize_thermchem(io)

!-- initialize grid module (must be done after reference module)
  call initialize_grid(io)

!-- initialize boundary conditions module
  call initialize_bc(io)

! wkliao: define MPI I/O file type to be used to define file view
  if (io_method .EQ. 1) then
      call mpi_io_set_filetype(0)
  elseif (io_method .EQ. 4) then
      call start_adios
  endif

!-- initialize field
  call initialize_field(io)

!-- generate active file (must be done last)
  call generate_active_file(io)

!-- tar input files for record
#ifndef SYSTEMCALLWONTWORK
!  call tar_input_files(io)
#endif

!-- initialize Rxn-wise Gibbs tabulation
#ifdef VECTORVERSION
!  call init_gibbstable
#endif

!----------------------------------------------------------------------------------------
! write files

  if(i_restart.eq.0) then
    call write_savefile(io)
!    call write_tecplot_skip(io,1)
!    call write_tecplot_animate_file(io,0,2)
!    call write_basic_tecplot_file(io,tec_del)
!    call write_1D_gnuplot_files(io)
  endif
!----------------------------------------------------------------------------------------
! check spatial resolution requirement of initial field

  if(i_time_res.ge.0) then
!    call check_resolution(io)
  endif
!----------------------------------------------------------------------------------------
! start loop over all time steps
!----------------------------------------------------------------------------------------
  i_time=0

  call date_and_time(dat_1,tim_start)

  do while (i_time < i_time_end)
!----------------------------------------------------------------------------------------
!   increment i_time

    i_time = i_time + 1

!   integrate in time

    !call integrate(io)
!----------------------------------------------------------------------
! wkliao:
!    if(myid .eq. 0) write(io, *) 'Sleeping for 10 seconds'
!    call sleep(10)

!   advance time step

    tstep = 1e-6/time_ref
    time=time+tstep

!   sync processors
!   Commented out by Ramanan - 04/18/05
!    call MPI_Barrier(gcomm,ierr)

!   monitor min/max values and look at active file for updates

    if(i_time_mon.gt.0) then

       if((i_time.le.10).or.(mod(i_time,i_time_mon).eq.0)) then

!        call monitor(io,temp_max,pres_max)  !max temp and pres passed out for write below

!       write timestep and error to screen

        if(myid.eq.0) then
          write(io,100) i_time,           &  !integer
                        time *(time_ref), &  !seconds
                        tstep*(time_ref), &  !seconds
                        q_err_max,        &  !dimensionless
                        temp_max,         &  !K
                        pres_max             !atm
        endif

      endif

!     read active file

      if(mod(i_time,i_time_mon).eq.0) then
!        call read_active_file(io)
      endif

    endif

!   write save file

    if(i_time_save.gt.0) then
      if((mod(i_time,i_time_save).eq.0).or.(time*time_ref.gt.time_save)) then
        call write_savefile(io)
        if(time*time_ref.gt.time_save) then  !increment save time
          time_save=time_save+time_save_inc
        endif
      endif
    endif

!   write tecplot file

    if((i_time_tec.gt.0).and.(mod(i_time,i_time_tec).eq.0)) then
!      call write_tecplot_skip(io,1)
!      call write_basic_tecplot_file(io,tec_del)
!      call write_1D_gnuplot_files(io)
    endif

!   check spatial resolution by analyzing the energy spectrum

    if(i_time_res.gt.0) then
      if(mod(i_time,i_time_res).eq.0) then
!        call check_resolution(io)
      endif
    endif

!   flush io

#   if PC || SGI
#   else
      if(myid.eq.0) call flush(io)
#   endif

!   always write save files for last timestep (if not already written)

    if((i_time.ge.i_time_end).and.(mod(i_time,i_time_save).ne.0)) then
!      call write_savefile(io)
    endif

!   always write tecplot files for last timestep (if not already written)

    if((i_time.ge.i_time_end).and.(mod(i_time,i_time_tec).ne.0)) then
!      call write_tecplot_skip(io,1)
!      call write_basic_tecplot_file(io,tec_del)
!      call write_1D_gnuplot_files(io)
    endif
!----------------------------------------------------------------------------------------
! end of loop over all time steps


  enddo

! wkliao: free up MPI I/O file type
  if (io_method .EQ. 1) then
      call mpi_io_set_filetype(-1)
  elseif (io_method .EQ. 4) then
      call end_adios
  endif

  call date_and_time(dat_1,tim_stop)
  if(myid==0) then
     call write_header(io,'=')
     write(io,*)'Timestepping began at: '
     call write_date_and_time(dat_1,tim_start,io)
     write(io,*)'Timestepping ended at: '
     call write_date_and_time(dat_1,tim_stop,io)
     call write_header(io,'=')
  endif

!----------------------------------------------------------------------------------------
! deallocate arrays

  call allocate_derivative_arrays(-1)
  !call allocate_filter_arrays(-1)
  call allocate_thermchem_arrays(-1)
!  call allocate_transport_arrays(-1)
!  call allocate_rk_arrays(-1)
  call allocate_variables_arrays(-1)
  call allocate_work_arrays(-1)
  call allocate_grid_arrays(-1)
  call allocate_bc_arrays(-1)
!----------------------------------------------------------------------------------------
! format statements

  100 format(i8,2x,1pe9.3,2x,1pe9.3,2x,1pe9.3,2x,0pf7.2,2x,0pf6.2)
!----------------------------------------------------------------------------------------
  return
  end subroutine solve_driver
