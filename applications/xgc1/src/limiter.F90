!***********************************************
! Limiter check routine 
! first created 2003/11/03
! last modifed 2003/03/07
! bug fix - 2003/03/27
! bug fix - 2003/11/03
!************************************************
!!$subroutine limiter_check(i)
!!$  use ptl_module
!!$  use lim_module
!!$  use eq_module, only : eq_x_z
!!$  use sml_module, only : sml_norm_r, sml_deltaf
!!$  implicit none
!!$  integer, intent(in) :: i
!!$  real (kind=8):: r,z,phi,psi,b,weight
!!$  integer :: j,rl,upper,lower,mid
!!$  real (kind=8), external :: psi_interpol, b_interpol
!!$  
!!$!  sp=ptl_species(i)
!!$  r=ptl_phase(i,1)
!!$  z=ptl_phase(i,2)
!!$  phi=ptl_phase(i,3)
!!$  psi=ptl_phase(i,9)
!!$    
!!$  if(z < eq_x_z .OR. psi> lim_psi_min) then
!!$     !trivial case
!!$     if(z>=lim_zmax) then
!!$        !store it to the top position
!!$        ptl_gid(i)=-ptl_gid(i)
!!$        if(r>lim_r0_up) then
!!$           rl=1
!!$        else
!!$           rl=2
!!$        endif
!!$        b=b_interpol(r,z,phi)
!!$        weight=ptl_phase(i,8)
!!$        if(sml_deltaf==1) then
!!$           weight=weight*ptl_phase(i,6)
!!$        endif
!!$        lim_en(lim_store_mz,rl)=lim_en(lim_store_mz,rl) &
!!$             + (0.5D0 * (ptl_phase(i,4)*B)**2 + ptl_phase(i,5)*B)*weight
!!$        lim_weight(lim_store_mz,rl)=lim_weight(lim_store_mz,rl) + weight
!!$
!!$     elseif(z<=lim_zmin) then
!!$        !store it to the bottom position
!!$        ptl_gid(i)=-ptl_gid(i)
!!$        if(r>lim_r0_down) then
!!$           rl=1
!!$        else
!!$           rl=2
!!$        endif
!!$        B=B_interpol(r,z,phi)
!!$        weight=ptl_phase(i,8)
!!$        if(sml_deltaf==1) then
!!$           weight=weight*ptl_phase(i,6)
!!$        endif
!!$        lim_en(1,rl)=lim_en(1,rl) &
!!$             + (0.5D0 * (ptl_phase(i,4)*B)**2 + ptl_phase(i,5)*B)*weight
!!$        lim_weight(1,rl)=lim_weight(1,rl) + weight
!!$
!!$     else
!!$        if( r > lim_r0_down .and. z<0 .OR. r>lim_r0_up .and. z>0) then !right
!!$           rl=1  !RIGHT
!!$        else
!!$           rl=2  !LEFT
!!$        endif
!!$        
!!$        lower=1
!!$        upper=lim_zindex(rl)
!!$        do while ( upper-lower > 1 )
!!$           mid=(lower+upper)/2
!!$           if(  z < lim_z(mid,rl) ) then
!!$              upper=mid
!!$           else
!!$              lower=mid
!!$           endif
!!$        enddo
!!$        ! check condition
!!$        if( (3.-real(2*rl))*  & 
!!$             (r- lim_r(lower,rl))*(lim_z(upper,rl)-lim_z(lower,rl)) &
!!$             > (3.-real(2*rl))*(lim_r(upper,rl)-lim_r(lower,rl))*(z-lim_z(lower,rl))  ) then
!!$           !outside
!!$           ptl_gid(i)=-ptl_gid(i)
!!$           !store it
!!$           j=(z-lim_zmin)/lim_dz + 1   ! dz=(zmax-zmin)/store_mz
!!$           j=max(1,min(lim_store_mz,j))
!!$           b=b_interpol(r,z,phi)
!!$           weight=ptl_phase(i,8)
!!$           if(sml_deltaf==1) then
!!$              weight=weight*ptl_phase(i,6)
!!$           endif
!!$           lim_en(j,rl)= lim_en(j,rl) &                ! energy sum
!!$                + (0.5D0 * (ptl_phase(i,4)*B)**2 + ptl_phase(i,5)*B)*weight
!!$           lim_weight(j,rl)=lim_weight(j,rl)+weight ! weight sum
!!$        endif
!!$     endif
!!$  endif
!!$end subroutine

subroutine limiter_setup
  use lim_module
  use eq_module, only : eq_x_z,eq_x_psi,eq_axis_r,eq_axis_z
  use sml_module, only : sml_2pi,&
       sml_bd_min_r, sml_bd_max_r, sml_bd_min_z, sml_bd_max_z
  implicit none
  real (kind=8) :: zmin, zmax,r,z,psi,psi_interpol,dr,old_r
  real (kind=8), allocatable :: tmp_r(:,:), tmp_z(:,:)
  integer :: i,j,izmin,izmax,offset,zindex(2),rl,size,index
  integer, parameter :: unit=92
  if(lim_filename/='') then
     call limiter_read
  endif
  allocate(tmp_r(lim_mdata,2), tmp_z(lim_mdata,2))
  tmp_r=eq_axis_r
  tmp_z=eq_axis_z
  zmin=0D0
  zmax=0D0

  do  i=1, lim_mdata
     if(zmin > lim_org_z(i)) then
        zmin=lim_org_z(i)
        izmin=i
     elseif(zmax < lim_org_z(i)) then
        zmax=lim_org_z(i)
        izmax=i
     endif
  enddo

  lim_zmin=zmin
  lim_zmax=zmax
  lim_dz=(lim_zmax-lim_zmin)/real(lim_store_mz)
  lim_r0_up=lim_org_r(izmax)
  lim_r0_down=lim_org_r(izmin)

  if( izmax<izmin ) then
     offset=lim_mdata
  else
     offset=0
  endif
!  print *, izmax,izmin,lim_zmin,lim_zmax
! stop
!------- One direction ---------------------------
  rl=1
  
  tmp_z(1,rl)=lim_org_z(izmin)
  tmp_r(1,rl)=lim_org_r(izmin)

  j=2
  do i= izmin+1, izmax+offset
     index=mod(i-1,lim_mdata)+1
     if(lim_org_z(index) > tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index)
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     else if ( lim_org_z(index) == tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index) + 1E-5
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     endif
  enddo
  zindex(rl)=j-1

!------- the other direction ---------------------------
  rl=2
  
  tmp_z(1,rl)=lim_org_z(izmin)
  tmp_r(1,rl)=lim_org_r(izmin)

  j=2
  do i= izmin-1, izmax+offset-lim_mdata , -1
     index=mod(i+lim_mdata-1,lim_mdata)+1
     if(lim_org_z(index) > tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index)
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     else if ( lim_org_z(index) == tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index) + 1E-5
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     endif
  enddo
  zindex(rl)=j-1

!----------- setting up ---------------------------------
  size=maxval(zindex)
  allocate( lim_r(size,2), lim_z(size,2) )
  allocate( lim_weight(lim_store_mz,2), lim_en(lim_store_mz,2),lim_ds(lim_store_mz,2))
  lim_weight=0D0
  lim_en=0D0

  if( maxval(tmp_r(:,1)) > maxval(tmp_r(:,2)) ) then
     lim_r=tmp_r(1:size,:)
     lim_z=tmp_z(1:size,:)
     lim_zindex=zindex
  else
     lim_r(:,1)=tmp_r(1:size,2)
     lim_r(:,2)=tmp_r(1:size,1)
     lim_z(:,1)=tmp_z(1:size,2)
     lim_z(:,2)=tmp_z(1:size,1)
     lim_zindex(1)=zindex(2)
     lim_zindex(2)=zindex(1)
  endif
  deallocate(tmp_r,tmp_z)
  ! finding smallest psi value -----------------------------
  lim_psi_min=1D30
  do rl=1, 2
     do i=1,lim_zindex(rl)
        r=lim_r(i,rl)
        z=lim_z(i,rl)
        !print *,r, z
        !if inside of interpolation range
        if(r > sml_bd_min_r .and. r < sml_bd_max_r &
             .and. z > sml_bd_min_z .and. z < sml_bd_max_z ) then
           psi=psi_interpol(r,z,0,0)
           if(z>eq_x_z .AND. psi < lim_psi_min) then
              lim_psi_min=psi
           endif
        endif
     enddo
  enddo

  !writing out the limiter data for plot and debug------------
  write(unit, *)'# limiter setup output'
  write(unit, *)'# num =', lim_zindex
  write(unit, *)'# psi_min = ', lim_psi_min/eq_x_psi
  write(unit, *)'# (r,z) values of limiter position'
  do rl=1,2
     do i=1, lim_zindex(rl)
        write(unit, *) lim_r(i,rl),lim_z(i,rl)
     enddo
     write(unit,*) ' '
  enddo


  !store the ds values

  do rl=1, 2
     j=1 
     old_r=lim_r0_down
     do i=1, lim_store_mz
        z=lim_zmin + real(i)*lim_dz
        z=min(lim_zmax,z)
        ! find proper j value
        do while( z<lim_z(j,rl) .OR. lim_z(j+1,rl) <z ) 
           j=j+1
        enddo
        r=  (lim_r(j+1,rl)-lim_r(j,rl)) * (z-lim_z(j,rl)) / (lim_z(j+1,rl)-lim_z(j,rl))&
             + lim_r(j,rl)
        dr=r-old_r
        lim_ds(i,rl)=sqrt(lim_dz**2+dr**2)*sml_2pi*( 0.5*(r+old_r))
        old_r=r

     enddo
  enddo

end subroutine limiter_setup

subroutine limiter_read
  use lim_module
  implicit none
  integer, parameter :: fileunit=114
  integer :: i
  character :: char(10)

  open(unit=fileunit,file=lim_filename,status='old',action='read')
!  read(fileunit,100) char
  read(fileunit,101) lim_mdata
  allocate(lim_org_r(lim_mdata),lim_org_z(lim_mdata))

  do i=1, lim_mdata
     read(fileunit,102) lim_org_r(i),lim_org_z(i)
  !   print *, lim_org_r(i), lim_org_z(i)
  enddo


100 format (a2)
101 format (i8)
102 format (e19.13,1x,e19.13)
  close(fileunit)

end subroutine

subroutine limiter_write
  use lim_module
  use sml_module, only :  sml_mype
  implicit none
  integer rl,i
  real (kind=8) :: tmp(lim_store_mz,2),dum(lim_store_mz,2) &
       ,sum_w(lim_store_mz,2),sum_e(lim_store_mz,2)
  real (kind=8) :: z
  


! mpi_sum - weight
  tmp=lim_weight
  call my_mpi_reduce(tmp,dum,lim_store_mz*2)
  sum_w=dum
!  tmp=lim_weight(:,2)
!  call my_mpi_reduce(tmp,dum,lim_store_mz)
!  sum_w(:,2)=dum
! mpi_sum - energy
  tmp=lim_en
  call my_mpi_reduce(tmp,dum,lim_store_mz*2)
  sum_e=dum
!  tmp=lim_en(:,2)
!  call my_mpi_reduce(tmp,dum,lim_store_mz)
!  sum_e(:,2)=dum


!debug
!  sum_w=lim_weight
!  sum_e=lim_en
  
  if(sml_mype==0) then
     write(35,*) '# z    outer_num   inner_num  outer_energy  inner_energy '
     do i=1,lim_store_mz
        z= (real(i)-0.5)*lim_dz + lim_zmin
        write(35,*) z, sum_w(i,1)/(lim_ds(i,1)), &
             sum_w(i,2)/(lim_ds(i,2)), &
             sum_e(i,1)/(lim_ds(i,1)),&
             sum_e(i,2)/(lim_ds(i,2)),&
             sum_w(i,1),sum_w(i,2),&
             sum_e(i,1),sum_e(i,2)
        
        ! lim_ds is the area sqrt(dz**2+dr**2)*2piR
        ! normalization factor need
     enddo
     write(35,*) ' '
  endif
end subroutine

     

