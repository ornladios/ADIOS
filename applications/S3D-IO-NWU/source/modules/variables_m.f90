!=========================================================================================
  module variables_m
!=========================================================================================
! module for variables variables

  implicit none
!-----------------------------------------------------------------------------------------
! derived variables

!  real, allocatable :: q(:,:,:,:,:) !solution vector (rho*u,rho*v,rho*w,rho,rho*e,rho*Y_i)

! primative variables

  real, allocatable :: yspecies(:,:,:,:)    !mass fractions for ALL species
  real, allocatable :: u(:,:,:,:)           !velocity vector (non-dimensional)
  real, allocatable :: volum(:,:,:)         !inverse of density (non-dimensional)
  real, allocatable :: pressure(:,:,:)      !pressure (non-dimensional)
  real, allocatable :: temp(:,:,:)          !temprature (non-dimensional)
!-----------------------------------------------------------------------------------------
  contains
!=========================================================================================
  subroutine allocate_variables_arrays(flag)
!=========================================================================================
! allocate variables arrays
!-----------------------------------------------------------------------------------------
  use param_m, only : nx, ny, nz, nsc, n_reg, nvar_tot

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed

  integer flag
!-----------------------------------------------------------------------------------------
! variables arrays

  if(flag.eq.1) then

!    allocate(q(nx,ny,nz,nvar_tot,n_reg));   q=0.0

    allocate(yspecies(nx,ny,nz,nsc+1));     yspecies=0.0
    allocate(u(nx,ny,nz,3));                u=0.0
    allocate(volum(nx,ny,nz));              volum=0.0
    allocate(pressure(nx,ny,nz));           pressure=0.0
    allocate(temp(nx,ny,nz));               temp=0.0    

  elseif(flag.eq.-1) then

!    deallocate(q)

    deallocate(yspecies)
    deallocate(u)
    deallocate(volum)
    deallocate(pressure)
    deallocate(temp)

  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine allocate_variables_arrays
!=========================================================================================
  subroutine get_mass_frac(q,volum,yspecies)
!=========================================================================================
! computes species mass fractions and 1/rho by extracting them from the U-vector
!
! q         - U-vector (rho*u, rho*v, rho*w, rho, rho*e_0, rho*Y_i)
! volum     - 1/rho
! yspecies  - species mass fraction -> Y_i
!             note that the array has an extra element
!             the (nsc+1)th element ultimately holds the degree
!             to which the sum of the mass fractions does NOT equal unity
!-----------------------------------------------------------------------------------------
  use param_m

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  real volum(nx,ny,nz)
  real q(nx,ny,nz,nvar_tot)
  real yspecies(nx,ny,nz,nsc+1)

! local declarations

  integer m, L
!-----------------------------------------------------------------------------------------
! calculate 1/rho

  volum(:,:,:) = 1.0/q(:,:,:,4)

! set nsc+1 species to unity

  yspecies(:,:,:,nsc+1) = 1.0

! calculate yspecies

  do m = 1, nsc
    yspecies(:,:,:,m) = volum(:,:,:) * q(:,:,:,5+m)
    yspecies(:,:,:,nsc+1) = yspecies(:,:,:,nsc+1) - yspecies(:,:,:,m)
  enddo
!-----------------------------------------------------------------------------------------
  return
  end subroutine get_mass_frac
!-----------------------------------------------------------------------------------------
  end module variables_m
!=========================================================================================
  subroutine get_velocity_vec(u,q,volum)
!=========================================================================================
! computes velocity vector extracting it from the U-vector
!-----------------------------------------------------------------------------------------
  use param_m, only : nx, ny, nz, nvar_tot

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  real u(nx,ny,nz,3)
  real q(nx,ny,nz,nvar_tot)
  real volum(nx,ny,nz)
!-----------------------------------------------------------------------------------------
! fill u-vector

  u(:,:,:,1) = q(:,:,:,1) * volum(:,:,:)
  u(:,:,:,2) = q(:,:,:,2) * volum(:,:,:)
  u(:,:,:,3) = q(:,:,:,3) * volum(:,:,:)
!-----------------------------------------------------------------------------------------
  return
  end subroutine get_velocity_vec
