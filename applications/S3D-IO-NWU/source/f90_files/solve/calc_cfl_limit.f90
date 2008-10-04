!=========================================================================================
  subroutine calc_cfl_limit(tstep_cfl,u,volum,gamma,pressure,  &
                            time_ref,delx,dely,delz,cfl_no)
!=========================================================================================
! calculates maximum timestep corresponding to inviscid CFL condition
! 
! BUG FIX 25-NOV-2004 Evatt Hawkes
! This routine was hard-wired for uniform grid.  Now using scale_1x to figure a local dx.
! Passing through the module to avoid changing the args for the subroutine.
!
!-----------------------------------------------------------------------------------------
  use topology_m
  use param_m, only : nx, ny, nz, vary_in_x, vary_in_y, vary_in_z 
  use work_m, only : a => work1_1  !alias speed of sound to work array to save storage

! for nonuniform mesh
  use param_m, only : nx_g, ny_g, nz_g
  use grid_m, only : unif_grid_x,unif_grid_y,unif_grid_z,scale_1x, scale_1y, scale_1z

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  real tstep_cfl
  real u(nx,ny,nz,3)
  real gamma(nx,ny,nz)
  real volum(nx,ny,nz)
  real pressure(nx,ny,nz)
  real time_ref
  real delx, dely, delz, cfl_no

! local declarations

  real tstep_cfl_l(3)     !local cfl tstep for each coordinate direction
  real tstep_cfl_g(3)     !global cfl tstep for each coordinate direction

  real large              !a really large number
  integer i               !counter
!-----------------------------------------------------------------------------------------
! return if zero dimensions

  tstep_cfl=0.0
  if((vary_in_x.eq.0).and.(vary_in_y.eq.0).and.(vary_in_z.eq.0)) return
!-----------------------------------------------------------------------------------------
! set large values

  large=huge(tstep_cfl)   !for setting large values

  tstep_cfl=large
  tstep_cfl_l=large
!-----------------------------------------------------------------------------------------
! calculate speed of sound

  a=sqrt(gamma*pressure*volum)
!-----------------------------------------------------------------------------------------
! calculate local minimum time step for x-direction

! 25-NOV-2004 Evatt Hawkes: nonuniform mesh fix

  if(vary_in_x.eq.1) then
    if(unif_grid_x==1)then
      tstep_cfl_l(1)=minval(delx/(abs(u(:,:,:,1))+a(:,:,:)))
    else
      do i = 1,nx
        tstep_cfl_l(1)=min(tstep_cfl_l(1),&
                           minval(real(nx_g-1)/(scale_1x(i)* (abs(u(i,:,:,1))+a(i,:,:)) )))
      enddo
    endif
  endif

  if(vary_in_y.eq.1) then
    if(unif_grid_y==1)then
      tstep_cfl_l(2)=minval(dely/(abs(u(:,:,:,2))+a(:,:,:)))
    else
      do i = 1,ny
        tstep_cfl_l(2)=min(tstep_cfl_l(2),&
                           minval(real(ny_g-1)/(scale_1y(i)*abs(u(:,i,:,2))+a(:,i,:))))
      enddo
    endif
  endif

  if(vary_in_z.eq.1) then
    if(unif_grid_z==1)then
      tstep_cfl_l(3)=minval(delz/(abs(u(:,:,:,3))+a(:,:,:)))
    else
      do i = 1,nz
        tstep_cfl_l(3)=min(tstep_cfl_l(3),&
                           minval(real(nz_g-1)/(scale_1z(i)*abs(u(:,:,i,3))+a(:,:,i))))
      enddo
    endif
  endif
!-----------------------------------------------------------------------------------------
! find global minimum

  call MPI_Allreduce(tstep_cfl_l,tstep_cfl_g,3,MPI_REAL8,MPI_MIN,gcomm,ierr)
!-----------------------------------------------------------------------------------------
! find direction minimum and multiply by time_ref and cfl number to get cfl time step

  tstep_cfl=cfl_no*minval(tstep_cfl_g)*time_ref
!-----------------------------------------------------------------------------------------
  return
  end subroutine calc_cfl_limit
