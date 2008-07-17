! turbulence diffusion
! first create 2002/11/13
! last modified 2003/10/05 by CS Chang
! last modified 2003/1025 - time varying turbulence diffusion

subroutine tbl_diffusion(istep,ptl)
  use ptl_module
  use tbl_module
  use sml_module
  use eq_module
  implicit none
  type(ptl_type) :: ptl
  type(species_type), pointer :: sp
  integer :: i,istep,k
  real (kind=8) :: r,z,tbl_d
  real (kind=8) :: dx, dt,psi
  real (kind=8), external :: psi_interpol, ranx
  target ptl
  real (kind=8) :: slope
  integer :: do_tbl
  integer :: tbl_slope_on, tbl_slope_side
  real (kind=8) :: tbl_slope_up, tbl_slope_down
  
  
  tbl_slope_on=1
  tbl_slope_up= 0.3638 ! + 20 degree
  tbl_slope_down = -0.3638 ! -20 degree

  if(istep > tbl_stop_time) return ! no turbulence diff after the first profile evaluation

  sp=>ptl%ion
   
  if(mod(istep,tbl_period)==0) then
     if(mod(istep,tbl_period*10)==0 .and. sml_mype==0) then
!       call tbl_coeff_out
     endif
     
      k=istep/tbl_mult_period+1   !2003/10/23
     k=min(k,tbl_mult_max)
     
     dt=sml_dt*tbl_period  ! time interval 
     ! Use a reflection condition at the boundaries.
    
     do i=1, sp%num
        if(sp%gid(i)>0) then
           r=sp%phase(1,i)
           z=sp%phase(2,i)
           ! diffusion
 !!          tbl_d = tbl_mult(k)*tbl_d_coeff
           if(tbl_slope_on==1) then  ! slope check
              do_tbl=0
              if(r>eq_axis_r) then ! right side only
                 slope= (z-eq_axis_z) /(r-eq_axis_r)
                 if( tbl_slope_down < slope .and. slope < tbl_slope_up) then
                    do_tbl=1
                 endif
              endif
           else
              do_tbl=1
           endif
              
           if(do_tbl == 1) then
              tbl_d= tbl_d_coeff
              dx=sqrt(0.5*dt*tbl_d)   ! diffusion step size
              r=r + dsign(1.0D0, 0.5D0-ranx())*dx 
              z=z + dsign(1.0D0, 0.5D0-ranx())*dx
              
              ! radial current
              !           r=r+ tbl_mu*er*dt
              !           z=z+ tbl_mu*ez*dt
              
              ! prevent diffusing out from simulation range
              psi=psi_interpol(r,z,0,0)
              
              !  Reflection condition at the boundaries
              if(psi > tbl_inpsi .AND. psi<sml_outpsi) then
                 sp%phase(1,i)=r
                 sp%phase(2,i)=z
                 sp%phase(9,i)=psi
              endif
           endif
        endif
     100 continue
     enddo
     sml_angle_stored=0
  endif


end subroutine tbl_diffusion

