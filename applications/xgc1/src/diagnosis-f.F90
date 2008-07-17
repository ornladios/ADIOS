!!$subroutine diagnosis_f(ptl,istep)
!!$  use sml_module
!!$  use diag_module  ! diag_f_nE, diag_f_npitch, diag_f_Emax(), diag_f_count
!!$  use ptl_module
!!$  use eq_module
!!$  use perf_monitor
!!$  implicit none
!!$  type(ptl_type), intent(in) :: ptl
!!$  type(species_type), pointer :: sp
!!$  integer, intent(in) :: istep
!!$  real (kind=8) , dimension(diag_f_nE, diag_f_npitch,diag_f_npsi,2) :: dumf, dumf_c
!!$  real (kind=8) :: r, z, psi, phi, rho, mu, w1, w3, psi_c
!!$
!!$
!!$  save dumf, dumf_c
!!$
!!$  sp=>ptl%ion
!!$  diag_f_count=diag_f_count+1
!!$  
!!$  do i=1, sp%num
!!$     if(sp%gid(i)>0) then
!!$        r=sp%phase(1,iptl)
!!$        z=sp%phase(2,iptl)
!!$        if( r> eq_axis_r .AND. z < diag_f_slope1*(r-eq_axis_r) .AND. z > diag_f_slope2*(r-eq_axis_r) .and. z>eq_x_z ) then
!!$
!!$           ! get E, pitch, psi, psi_c
!!$           phi=sp%phase(3,iptl)
!!$           rho=sp%phase(4,iptl)
!!$           mu=sp%phase(5,iptl)
!!$           w1 = sp%phase(6,iptl)
!!$           w3 = sp%phase(8,iptl)           
!!$           psi=sp%phase(9,iptl)
!!$           call bphi_interpol_rzpsi(r,z,phi,psi,bphi)
!!$           if(sml_minusB/=1) then
!!$              psi_c=psi + rho*bphi*r
!!$           else
!!$              psi_c=psi - rho*bphi*r
!!$           endif
!!$           
!!$           ! 8-point (3D) linear interpolation         
!!$           ipsi(1)=int((psi-diag_f_pin)/diag_f_dpsi) + 1
!!$           ipsi_c(1)=int((psi_c-diag_f_pin)/diag_f_dpsi) + 1
!!$           iE(1)
!!$           ipitch(1)
!!$
!!$           do jp=1,2
!!$              do je=1,2
!!$                 do jk=1,2
!!$                
!!$                    
!!$
!!$    
!!$                    
!!$                 enddo
!!$              enddo
!!$           enddo
!!$           
!!$           
!!$        endif
!!$     endif
!!$  enddo
!!$  
!!$  
!!$
!!$
!!$
!!$end subroutine diagnosis_f

subroutine diagnosis_f(ptl,istep)
  use sml_module
  use diag_module
  use ptl_module
  use eq_module
  use rf_module
  use perf_monitor
  implicit none
  type(ptl_type), intent(in) :: ptl
  type(species_type), pointer :: sp
  integer,intent(in) :: istep
  real (kind=8),dimension(diag_f_nv,-diag_f_nv:diag_f_nv,diag_f_npsi) :: dumf,dumf0,dumf_c,dumf0_c
  real (kind=8) :: dv(diag_f_npsi)
  real (kind=8) :: r, z, psi,phi, rho, mu, w1, w3, b, v_para, v_perp,psi_c,bphi
  real (kind=8) :: wp, wl, f0, ff , jc_v, jc_x, jc_inv,sum1,sum2,sum1_c,sum2_c
  integer :: p1, p2, l1, l2
  integer :: i,j,k, iptl,ipsi,ipsi_c
  real (kind=8),external :: b_interpol
  character (len=20):: filename
  target ptl
  save dumf, dumf0,dumf_c,dumf0_c

  sp=>ptl%ion

  diag_f_count=diag_f_count+1

  do i = 1, diag_f_npsi
     dv(i) = diag_f_vmax*diag_f_vt(i)/real(diag_f_nv)
  enddo

  do iptl=1, sp%num
     if(sp%gid(iptl)>0) then
        r=sp%phase(1,iptl)
        z=sp%phase(2,iptl)
        phi=sp%phase(3,iptl)
        rho=sp%phase(4,iptl)
        mu=sp%phase(5,iptl)
        w1 = sp%phase(6,iptl)
        w3 = sp%phase(8,iptl)

        if( r> eq_axis_r .AND. z < diag_f_slope1*(r-eq_axis_r) .AND. z > diag_f_slope2*(r-eq_axis_r) .and. z>eq_x_z ) then
!        if(1==1) then
           psi=sp%phase(9,iptl)
           call bphi_interpol_rzpsi(r,z,phi,psi,bphi)
           if(sml_minusB/=1) then
              psi_c=psi + rho*bphi*r
           else
              psi_c=psi - rho*bphi*r
           endif
     
           ipsi= int((psi-diag_f_pin)/diag_f_dpsi) + 1
           ipsi_c=int((psi_c-diag_f_pin)/diag_f_dpsi) + 1
           if(  ipsi_c < 1 .or. ipsi_c > diag_f_npsi) ipsi_c=0
           if( 1 <= ipsi .and.  ipsi <= diag_f_npsi) then
              b=b_interpol(r,z,phi)

              ! different normalization for each species ???
              v_para = ptl_c_m(1)*rho*b
              v_perp = sqrt(2D0*B*mu/ptl_mass(1))
              ! --------------------------------------------

              if(v_para >= 0) then
                 l1 = int(v_para/dv(ipsi)) 
              else 
                 l1 = int(v_para/dv(ipsi))  - 1
              endif
              l2 = l1 + 1
              wl = v_para/dv(ipsi) - l1


              p1 = int(v_perp/dv(ipsi)) + 1
              p2 = p1 + 1
              !wp = v_perp/dv(ipsi, sp) - p1 +1
              wp = 0  ! prevent singularity -- nearest neighbor

              jc_x= diag_f_n0(ipsi)*diag_f_dvol(ipsi)/sml_sqrt2pi
              jc_v = dv(ipsi)**3*(2*p1-1)/&
                   (2.0*diag_f_vt(ipsi)*diag_f_vt(ipsi)*diag_f_vt(ipsi))
              jc_inv = 1.0/(jc_x*jc_v)
!               jc_inv=1
    
              ff = w3*w1 + w3
              f0 = w3
              !dubug check
!              if( jc_inv<0.)  print *, 'minus jc_inv', jc_inv
!              if( wl <0. .or. wl > 1.) print *, 'worng wl', wl, l1, l2, v_para, dv(ipsi,sp)
!              if( w3 <0. ) print *, 'wrong w3' , w3
!              if( dv(ipsi,sp) < 0. ) print *, 'wring dv', dv(ipsi,sp)

              if(p1 >= 1 .and. p1 <= diag_f_nv .and. l1 >= -diag_f_nv .and. l1 <= diag_f_nv) then
                 diag_f(p1,l1,ipsi) = diag_f(p1,l1,ipsi) &
                      + (1.0-wl)*(1.0-wp)*jc_inv*ff
                 diag_f0(p1,l1,ipsi) = diag_f0(p1,l1,ipsi) &
                      + (1.0-wl)*(1.0-wp)*jc_inv*f0
                 if(ipsi_c>0) then
                     diag_f_c(p1,l1,ipsi_c) = diag_f(p1,l1,ipsi_c) &
                          + (1.0-wl)*(1.0-wp)*jc_inv*ff
                     diag_f0_C(p1,l1,ipsi_c) = diag_f0(p1,l1,ipsi_c) &
                          + (1.0-wl)*(1.0-wp)*jc_inv*f0
                  endif                    
              endif

              if(p1 >= 1 .and. p1 <= diag_f_nv .and. l2 >= -diag_f_nv .and. l2 <= diag_f_nv) then
                 diag_f(p1,l2,ipsi) = diag_f(p1,l2,ipsi) &
                      + wl*(1.0-wp)*jc_inv*ff
                 diag_f0(p1,l2,ipsi) = diag_f0(p1,l2,ipsi) &
                      + wl*(1.0-wp)*jc_inv*f0
                 if(ipsi_c>0) then
                     diag_f_c(p1,l2,ipsi_c) = diag_f(p1,l2,ipsi_c) &
                          + wl*(1.0-wp)*jc_inv*ff
                     diag_f0_c(p1,l2,ipsi_c) = diag_f0(p1,l2,ipsi_c) &
                          + wl*(1.0-wp)*jc_inv*f0
                  endif
              endif

              if(p2 >= 1 .and. p2 <= diag_f_nv .and. l2 >= -diag_f_nv .and. l2 <= diag_f_nv) then
                 diag_f(p2,l2,ipsi) = diag_f(p2,l2,ipsi) &
                      + wl*wp*jc_inv*ff
                 diag_f0(p2,l2,ipsi) = diag_f0(p2,l2,ipsi) &
                      + wl*wp*jc_inv*f0
                 if(ipsi_c>0) then
                     diag_f_c(p2,l2,ipsi_c) = diag_f(p2,l2,ipsi_c) &
                          + wl*wp*jc_inv*ff
                     diag_f0_c(p2,l2,ipsi_c) = diag_f0(p2,l2,ipsi_c) &
                          + wl*wp*jc_inv*f0
                  endif
              endif

              if(p2 >= 1 .and. p2 <= diag_f_nv .and. l1 >= -diag_f_nv .and. l1 <= diag_f_nv) then
                 diag_f(p2,l1,ipsi) = diag_f(p2,l1,ipsi) &
                      + (1.0-wl)*wp*jc_inv*ff
                 diag_f0(p2,l1,ipsi) = diag_f0(p2,l1,ipsi) &
                      + (1.0-wl)*wp*jc_inv*f0
                 if(ipsi_c>0) then
                    diag_f_c(p2,l1,ipsi_c) = diag_f(p2,l1,ipsi_c) &
                         + (1.0-wl)*wp*jc_inv*ff
                    diag_f0_c(p2,l1,ipsi_c) = diag_f0(p2,l1,ipsi_c) &
                         + (1.0-wl)*wp*jc_inv*f0
                 endif
              endif
           endif
        endif
     endif
  end do

  if(diag_f_count >= diag_f_mod ) then
     call monitor_start (DIAG_F_RED_)
     call my_mpi_reduce(diag_f,dumf,diag_f_nv*(2*diag_f_nv+1)*diag_f_npsi)
     call my_mpi_reduce(diag_f0,dumf0,diag_f_nv*(2*diag_f_nv+1)*diag_f_npsi)
     call my_mpi_reduce(diag_f_c,dumf_c,diag_f_nv*(2*diag_f_nv+1)*diag_f_npsi)
     call my_mpi_reduce(diag_f0_c,dumf0_c,diag_f_nv*(2*diag_f_nv+1)*diag_f_npsi)
     call monitor_stop (DIAG_F_RED_)

     if(sml_mype==0) then
        do i=1, diag_f_npsi     
           write(filename,'("fort.f",".",i2.2,".",i3.3)') i, istep/(diag_f_skip*diag_f_mod)
           open(333, file=filename, status='replace', form='formatted') 
           write(333,*) '# time=', sml_time, 'psi=', diag_f_pin+(i-0.5)*diag_f_dpsi
           do j=1, diag_f_nv
              do k=-diag_f_nv, diag_f_nv
                 write(333,1000) diag_f_vmax*real(k)/real(diag_f_nv),&
                      diag_f_vmax*(j-0.5)/real(diag_f_nv), &
                      dumf(j,k,i)/real(diag_f_count),&
                      dumf0(j,k,i)/real(diag_f_count),&
                      dumf_c(j,k,i)/real(diag_f_count),&
                      dumf0_c(j,k,i)/real(diag_f_count)
              enddo
              write(333,*)' '
           enddo
           close(333)
           ! canonical distribution function
           !write(filename,'("fort.fc",".",i2.2,".",i3.3)') i, istep/(diag_f_skip*diag_f_mod)
           !open(333, file=filename, status='replace', form='formatted') 
           !write(333,*) '# time=', sml_time, 'psi_c=', diag_f_pin+(i-0.5)*diag_f_dpsi
           !do j=1, diag_f_nv
           !   do k=-diag_f_nv, diag_f_nv
           !      write(333,1000) diag_f_vmax*real(k)/real(diag_f_nv),&
           !           diag_f_vmax*(j-0.5)/real(diag_f_nv), &
           !           dumf_c(j,k,i)/real(diag_f_count),dumf0_c(j,k,i)/real(diag_f_count)
           !   enddo
           !   write(333,*)' '
           !enddo
           !close(333)
        enddo
        
        do i=1, diag_f_npsi
           ! v_perp dist
           write(filename,'("fort.fperp",".",i2.2,".",i3.3)') i, istep/(diag_f_skip*diag_f_mod)
           open(333, file=filename, status='replace', form='formatted') 
           do j=1, diag_f_nv
              sum1=sum(dumf(j,:,i))
              sum2=sum(dumf0(j,:,i))
              sum1_c=sum(dumf_c(j,:,i))
              sum2_c=sum(dumf0_c(j,:,i))
              write(333,1001) diag_f_vmax*(j-0.5)/real(diag_f_nv), &
                   sum1/real(diag_f_count),sum2/real(diag_f_count),&
                   sum1_c/real(diag_f_count),sum2_c/real(diag_f_count)
           enddo
           close(333)
           
           ! v_para dist
           write(filename,'("fort.fpara",".",i2.2,".",i3.3)') i, istep/(diag_f_skip*diag_f_mod)
           open(333, file=filename, status='replace', form='formatted') 
           do k=-diag_f_nv, diag_f_nv
              sum1=sum(dumf(:,k,i))
              sum2=sum(dumf0(:,k,i))
              sum1_c=sum(dumf_c(:,k,i))
              sum2_c=sum(dumf_c(:,k,i))
              write(333,1001) diag_f_vmax*real(k)/real(diag_f_nv), &
                   sum1/real(diag_f_count),sum2/real(diag_f_count),&
                   sum1_c/real(diag_f_count),sum2_c/real(diag_f_count)
           enddo
           close(333)
        enddo
     
1000    format(6(e19.13,' '))
1001    format(5(e19.13,' '))
     endif
     diag_f_count=0
     diag_f=0D0
     diag_f0=0D0
  endif
end subroutine diagnosis_f
