!=========================================================================================
module bc_m
!=========================================================================================
! module for bc variables

! conventions:
!  nrf_x0  controls the BC type on the X boundary.
!
!      = (1)  ->  Nonreflecting BC's with appropriate diffusive BC's
!                 See files:  nscbc.f90, compute_d.f90, compute_L.f90, update_L.f90, bc_flux.f90
!
!      = (0)  ->  Standard hard inflow.  Velocity components, T, Yi, and the normal stress
!                 are specified on the boundary.  See files: impose_hard_bc.f90, bc_flux.f90
!
!      = (-2) ->  Hard inflow with a Characteristic correction on density.  No condition on any
!                 flux quantities.  See files: inflowA_BC.f90 as well as impose_hard_bc.f90
!
!      = (-5) ->  No-slip, adiabatic wall.
!                 See files: impose_hard_bc.f90, bc_flux.f90
!
!      = (-6) ->  No-slip, isothermal wall.
!                 See files: impose_hard_bc.f90, bc_flux.f90
!
  implicit none
  public
  
  integer nrf_x0    !non-reflecting/hard BC at x=0 boundary
  integer nrf_xl    !non-reflecting/hard BC at x=Lx boundary
  integer nrf_y0    !non-reflecting/hard BC at y=0 boundary
  integer nrf_yl    !non-reflecting/hard BC at y=Ly boundary
  integer nrf_z0    !non-reflecting/hard BC at z=0 boundary
  integer nrf_zl    !non-reflecting/hard BC at z=Lz boundary

  real relax_ct     !relaxation parameter for non-reflecting boundary conditions
  real pout         !initial pressure (saved) for non-reflecting boundary conditions
  real :: beta=15.0 !relaxation for soft inflow boundary-normal velocity component

  real, allocatable :: qx_bc(:,:,:,:)  !hard-inflow storage array for x-direction
  real, allocatable :: qy_bc(:,:,:,:)  !hard-inflow storage array for y-direction
  real, allocatable :: qz_bc(:,:,:,:)  !hard-inflow storage array for z-direction
!-----------------------------------------------------------------------------------------

contains

!=========================================================================================

  subroutine allocate_bc_arrays(flag)
    use param_m, only : nx, ny, nz, nvar_tot, nsc
    implicit none
    integer, intent(in) :: flag
    if(flag.eq.1) then

       if(nrf_x0<=0 .or. nrf_xl<=0) then
          allocate(qx_bc(2,ny,nz,nvar_tot))
          qx_bc=0.0
       endif

       if(nrf_y0<=0 .or. nrf_yl<=0) then
          allocate(qy_bc(nx,2,nz,nvar_tot))
          qy_bc=0.0
       endif

       if(nrf_z0<=0 .or. nrf_zl<=0) then
          allocate(qz_bc(nx,ny,2,nvar_tot))
          qz_bc=0.0
       endif

    elseif(flag.eq.-1) then
       if(nrf_x0<=0 .or. nrf_xl<=0)  deallocate(qx_bc)
       if(nrf_y0<=0 .or. nrf_yl<=0)  deallocate(qy_bc)
       if(nrf_z0<=0 .or. nrf_zl<=0)  deallocate(qz_bc)
    endif
    return
  end subroutine allocate_bc_arrays

!=========================================================================================

  subroutine initialize_bc(io)
    use topology_m, only : myid
    implicit none
    integer, intent(in) :: io

    if(myid.eq.0) then
       write(io,*) 'initializing boundary condition module...'
       write(io,*)
    endif
    call allocate_bc_arrays(1)
    if(myid.eq.0) then
       call write_header(io,'-')
    endif
    return
  end subroutine initialize_bc

!=========================================================================================

  subroutine store_hard_bc(io,q,temp,pressure,volum,i_restart)
    ! inititialize bc arrays
    use topology_m, only : myid
    use param_m, only : nx, ny, nz, nvar_tot, n_reg
    implicit none

    integer, intent(in) :: io
    real, intent(in) ::  q(nx,ny,nz,nvar_tot,n_reg)
    real, intent(in), dimension(nx,ny,nz) :: temp, pressure, volum    !1/rho
    integer, intent(in) :: i_restart

    integer i,j,k,L

    if(myid.eq.0) then
       write(io,*) 'storing hard inflow conditions...'
       write(io,*)
    endif
  !-----------------------------------------------------------------------------------------
  ! store the primative values of variables at boundaries
  ! note the division by rho everywhere
  ! L = 1       -> u-component of velocity
  ! L = 2       -> v-component of velocity
  ! L = 3       -> w-component of velocity
  ! L = 4       -> unity (density is never imposed at boundaries)
  ! L = 5       -> total energy
  ! L = 5 + nsc -> species
  !
  ! save boundary temperature for hard inlet condition (over-write previous energy values)
  ! L = 5 -> temperature

    if(nrf_x0<=0 .or. nrf_xl<=0) then
       do L=1,nvar_tot
          qx_bc(1,:,:,L) = q( 1,:,:,L,1)*volum( 1,:,:)
          qx_bc(2,:,:,L) = q(nx,:,:,L,1)*volum(nx,:,:)
       enddo
       qx_bc(1,:,:,5) = temp( 1,:,:)
       qx_bc(2,:,:,5) = temp(nx,:,:)
    endif

    if(nrf_y0<=0 .or. nrf_yl<=0) then
       do L=1,nvar_tot
          qy_bc(:,1,:,L) = q(:, 1,:,L,1)*volum(:, 1,:)
          qy_bc(:,2,:,L) = q(:,ny,:,L,1)*volum(:,ny,:)
       enddo
       qy_bc(:,1,:,5) = temp(:, 1,:)
       qy_bc(:,2,:,5) = temp(:,ny,:)
    endif

    if(nrf_z0<=0 .or. nrf_zl<=0) then
       do L=1,nvar_tot
          qz_bc(:,:,1,L) = q(:,:, 1,L,1)*volum(:,:, 1)
          qz_bc(:,:,2,L) = q(:,:,nz,L,1)*volum(:,:,nz)
       enddo
       qz_bc(:,:,1,5) = temp(:,:, 1)
       qz_bc(:,:,2,5) = temp(:,:,nz)
    endif

    ! set pout

    if(i_restart.ne.1) then
       pout=pressure(1,1,1)    !only valid for uniform initial pressure, see Scott
    endif

    return
  end subroutine store_hard_bc

!=========================================================================================

end module bc_m
