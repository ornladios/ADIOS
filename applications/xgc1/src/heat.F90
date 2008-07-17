subroutine heating(istep,ptl)
  use heat_module
  use ptl_module
  use sml_module
  use perf_monitor
  implicit none
  type(ptl_type) :: ptl
  type(species_type), pointer :: sp
  real (kind=8) :: En_sum, En, psi,B, alpha, sqrt_alpha, sum(1), dum(1)
  real (kind=8) , external :: B_interpol, psi_interpol
  integer :: i,istep,k
  target ptl
  
  sp=>ptl%ion
  k=istep/heat_mult_period+1
  k=min(k,heat_mult_max)
  En_sum=0D0
  do i= 1,sp%num
     if(sp%gid(i) > 0) then
        psi= sp%phase(9,i)
        if(heat_inpsi < psi .AND. psi < heat_outpsi) then
           B=B_interpol(sp%phase(1,i), sp%phase(2,i),sp%phase(3,i)) 
           En=sp%phase(5,i)*B + 0.5D0*ptl_c2_2m(1)*(sp%phase(4,i)*B)**2
           En_sum=En_sum+En*sp%phase(8,i)
        endif
     endif
  enddo
  
  ! 
  sum(1)=En_sum
  call monitor_start (HEAT_RED_)
  call my_mpi_allreduce(sum, dum, 1)
  call monitor_stop (HEAT_RED_)
  En_sum=dum(1)

  alpha= (En_sum + heat_mult(k)* heat_power*sml_dt*real(heat_period))/En_sum
  sqrt_alpha=sqrt(alpha)
  

  do i=1, sp%num
     if(sp%gid(i) > 0) then
        psi=sp%phase(9,i)
        if(heat_inpsi < psi .AND. psi < heat_outpsi) then
           sp%phase(5,i)=sp%phase(5,i)*alpha
           sp%phase(4,i)=sp%phase(4,i)*sqrt_alpha
        endif
     endif
  enddo
end subroutine heating
     
