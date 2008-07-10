subroutine output_bfield
  use eq_module
  use sml_module
  implicit none
  real (kind=8) :: psi,br,bz,bphi
  real (kind=8) :: r,z,rstart,rend,dr
  real (kind=8), external:: psi_interpol
  integer :: i,N=1000
#ifdef ADIOS
#define ADIOS_WRITE(a,b) call adios_write (a,'b'//char(0),b)
#define ADIOS_WRITE_VAR(a,b,c) call adios_write (a,'b'//char(0),c)
  integer*8 :: buf_id,grp_id
  real (kind=8), allocatable :: rvec(:),psi_eq_x_psi(:),dpdr(:),abs_bphi(:)
  allocate(rvec(N),psi_eq_x_psi(N),dpdr(N),abs_bphi(N))
#endif
  z=0
  rstart=sml_bd_min_r ! general r start --2002/02/01
  rend=sml_bd_max_r ! general r end -- 2002/02/01
  dr=(rend-rstart)/real(N)
#ifndef ADIOS
  write(60,*) '#R   Psi/eq_x_psi  dpdr  dpdz   d2pdr2 d2pdrdz d2pdz2'
#endif
  do i=1, N
     r=rstart+dr*real(i)
     psi=psi_interpol(r,z,0,0)
     call bvec_interpol(r,z,3.141592/2.,br,bz,bphi)
#ifdef ADIOS
     rvec(i)=r
     psi_eq_x_psi(i)=psi/eq_x_psi;
     dpdr(i)=dsqrt(br**2+bz**2)
     abs_bphi(i)=abs(bphi)
#else
     write(60,1000) (r),psi/eq_x_psi, dsqrt(br**2+bz**2), abs(bphi), r
!          if(psi > NC_x_psi*1.1) then
!             return
!          endif
#endif 
  enddo
#ifdef ADIOS
!  call adios_open(buf_id,grp_id,"fort.60"//char(0),"fort.bp"//char(0),"a"//char(0))
  call adios_open(buf_id,"fort.60"//char(0),"fort.bp"//char(0),"a"//char(0))
!  call adios_set_path(grp_id,"fort.60"//char(0))
  call adios_set_path(buf_id,"fort.60"//char(0))
  call adios_gwrite(buf_id,"fort.60")
  call adios_close(buf_id)
  deallocate(rvec,psi_eq_x_psi,dpdr,abs_bphi)
#endif
1000 format (5(e19.13,' '))
end subroutine output_bfield

subroutine output_initial
  use eq_module
  use efld_module
  use sml_module
  implicit none
  real (kind=8) :: r,z,rstart,rend,dr
  real (kind=8) :: er,ez,ephi,epot,psi,psi_interpol,br,bz,bphi
  real (kind=8), external ::init_tempi_ev,tempe_ev,init_den,den_imp
  integer :: i,N=1000
#ifdef ADIOS
  integer*8 :: buf_id,grp_id
  real (kind=8), allocatable :: r_arr(:),Epoten(:),psi_eq_x_psi(:),Er_minor_radial(:),er_arr(:),ez_arr(:),ephi_arr(:),ion_temp(:),e_temp(:),back_den(:),impurity_den(:)
  allocate(r_arr(N),psi_eq_x_psi(N),Epoten(N),Er_minor_radial(N),er_arr(N),ez_arr(N),ephi_arr(N),ion_temp(N),e_temp(N),back_den(N),impurity_den(N))
#endif
  z=0
  rstart=eq_axis_r ! to start from magnetic axis
  !  rstart=0.25D0  ! starting position of r
  !  rend=1.6       ! end postion of r -- 1.6 is NSTX
   rend=sml_bd_max_r ! general end position of r --2002/02/01
  dr=(rend-rstart)/real(N)
#ifndef ADIOS
  write(61,*) '#Psi     Epoten   Er(minor-radial) er  ez  ephi   r   ion-temp e-temp back-den impurity-den'
#endif
! Initialize efield values to zero
  er=0.
  ez=0.
  ephi=0.
  epot=0.
  do i=1, N
     r=rstart+dr*real(i)
     psi=psi_interpol(r,z,0,0)
     call bvec_interpol(r,z,3.141592/2.,br,bz,bphi)
!     call efield_static(r,z,er,ez,ephi,epot)
#ifndef ADIOS
     write(61,1000) psi/eq_x_psi,epot, &
         dsqrt(er**2+ez**2), &
         er,ez,ephi,(r-eq_axis_r),init_tempi_ev(psi),tempe_ev(psi), &
         init_den(psi), den_imp(psi)
#else
!     if(psi > NC_x_psi*1.1) then
!        return
!     endif
     psi_eq_x_psi(i)=psi/eq_x_psi
     Epoten(i)=epot
     Er_minor_radial(i)=dsqrt(er**2+ez**2)
     er_arr(i)=er
     ez_arr(i)=ez
     ephi_arr(i)=ephi
     r_arr(i)=r-eq_axis_r
     ion_temp(i)=init_tempi_ev(psi)
     e_temp(i)=tempe_ev(psi)
     back_den(i)=init_den(psi)
     impurity_den(i)=den_imp(psi)
#endif
  enddo
#ifdef ADIOS
  call adios_open(buf_id,"fort.61"//char(0),"fort.bp"//char(0),"a"//char(0))
  !call adios_open(buf_id,grp_id,"fort.61"//char(0),"fort.bp"//char(0),"a"//char(0))
  !call adios_set_path(grp_id,"fort.61"//char(0))
  call adios_set_path(buf_id,"fort.61"//char(0))
  call adios_gwrite(buf_id,"fort.61")
  call adios_close(buf_id)
  deallocate(psi_eq_x_psi,Epoten,Er_minor_radial,er_arr,ez_arr,ephi_arr,r_arr,ion_temp,e_temp,back_den,impurity_den)
#else
  close(61)
#endif
1000 format (11(e19.13,' '))
!1001 format (2(e19.13,' '))
end subroutine

! routines to check interpolation vailidty
subroutine output_psi_der
  use eq_module
  use sml_module
  implicit none
  real (kind=8) :: psi,br,bz,bphi
  real (kind=8) :: r,z,rstart,rend,dr,zstart,zend,dz
  real (kind=8), external:: psi_interpol
  integer :: i,Nr,Nz,j
  Nz=10
  Nr=1000
  zstart=sml_bd_min_z !
  zend=sml_bd_max_z 
  dz=(zend-zstart)/real(Nz)

  rstart=sml_bd_min_r! general r start --2002/02/01
  rend=sml_bd_max_r ! general r end -- 2002/02/01
  dr=(rend-rstart)/real(Nr)
  write(60,*) '#minor_r   Psi/eq_x_psi  bp  bt   R'

  do j=1, Nz

  do i=1, Nr
     r=rstart+dr*real(i)
     z=zstart+dz*real(j)
     psi=psi_interpol(r,z,0,0)
     call bvec_interpol(r,z,3.141592/2.,br,bz,bphi)
     write(1060,1000) (r),psi/eq_x_psi, psi_interpol(r,z,1,0),psi_interpol(r,z,0,1),&
          psi_interpol(r,z,2,0),psi_interpol(r,z,1,1),psi_interpol(r,z,0,2)
     !     if(psi > NC_x_psi*1.1) then
     !        return
     !     endif
  enddo
     write(1060,*) ' '
  enddo



1000 format (7(e19.13,' '))
end subroutine output_psi_der
