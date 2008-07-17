!****************************************************************************
! retrun I
!
! first created : 2000/08/26
! last modified : 2002/02/08
!adopted from xorbit
!
!****************************************************************************
real (kind=8)  function I_interpol(in_psi, ideriv,region)
#if defined(PSPLINE)
    use EZspline_obj
    use EZspline
    use itp_module, only : spl_psi,itp_min_psi,itp_max_psi
#else 
! IMSL or MY_IMSL wrappers or portable pppack routines
    use itp_module, only : itp_I_cscoef,itp_I_break,itp_mpsi,itp_min_psi,itp_max_psi
#endif
    use eq_module, only : eq_x_psi
    use sml_module, only : sml_bt_sign

    implicit none
    integer , intent(in) :: ideriv,region
    real (kind=8) , intent(in) :: in_psi
    real (kind=8)  :: psi
#if defined(PSPLINE)
    integer :: ier
    real (kind=8) :: r8value
#elif defined(MY_IMSL)
    real (kind=8) , external :: my_dppder
#elif defined(IMSL)
    real (kind=8) , external :: dcsder
#else
! portable pppack routines
    real (kind=8) , external :: my_dppder
#endif


!    sign=-1D0 !-1D0 : cocurrent, 1D0 :counter current 
!    sml_bt_sign can be changed in setup.f90 2002/02/08
    
    if(region == 3 ) then

       psi=min(eq_x_psi,itp_max_psi) ! for itp_max_psi < eq_x_psi case 2002/01/22
       if(ideriv == 0) then
#if defined(PSPLINE)
          call EZspline_interp(spl_psi,psi,r8value,ier)
          call EZspline_error(ier)
          I_interpol=sml_bt_sign*r8value
#elif defined(MY_IMSL)
          I_interpol=sml_bt_sign*my_dppder(ideriv,psi,4,itp_mpsi-1,itp_I_break,itp_I_cscoef)
#elif defined(IMSL)
          I_interpol=sml_bt_sign*dcsder(ideriv,psi,itp_mpsi-1,itp_I_break,itp_I_cscoef)
#else
          I_interpol=sml_bt_sign*my_dppder(ideriv,psi,4,itp_mpsi-1,itp_I_break,itp_I_cscoef)
#endif
       else
          I_interpol=0
       endif

    else

       psi = in_psi
       if(psi < itp_min_psi) then
          if(psi < itp_min_psi - 1D-4) then
             print *, 'psi range exceeded',psi
             call err_count
          endif
	  psi=itp_min_psi
       else if(psi > itp_max_psi) then
	  psi=itp_max_psi ! I is constant outside of edge
       endif
#if defined(PSPLINE)
       call EZspline_derivative(spl_psi,ideriv,psi,r8value,ier)
       call EZspline_error(ier)
       I_interpol=sml_bt_sign*r8value
#elif defined(MY_IMSL)
       I_interpol=sml_bt_sign*my_dppder(ideriv,psi,4,itp_mpsi-1,itp_I_break,itp_I_cscoef)
#elif defined(IMSL)
       I_interpol=sml_bt_sign*dcsder(ideriv,psi,itp_mpsi-1,itp_I_break,itp_I_cscoef)
#else
       I_interpol=sml_bt_sign*my_dppder(ideriv,psi,4,itp_mpsi-1,itp_I_break,itp_I_cscoef)
#endif

    endif
end function
!**************************************************************
! B-field interpolation
!**************************************************************

real (kind=8) function b_interpol(r,z,phi)
  use eq_module
  use sml_module, only : sml_minusB
  use rpl_module
  implicit none
  real (kind=8) , intent(in) :: r,z,phi
  real (kind=8)              :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=8)              :: ripp,dripp_dr,dripp_dz
  real (kind=8) , external   :: I_interpol,psi_interpol
  psi     = psi_interpol(r,z,0,0)
  dpsi_dr = psi_interpol(r,z,1,0)
  dpsi_dz = psi_interpol(r,z,0,1)
  
  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol(psi,0,3)
  else
     fi=I_interpol(psi,0,1)
  endif
  
  br=- dpsi_dz / r
  bz= dpsi_dr / r
  bphi=fi / r
  if(rpl_mode==1) then
     call get_ripple(r,z,ripp,dripp_dr,dripp_dz)
     bphi=bphi*(1D0 + dcos(rpl_N_coil*phi)*ripp) 
  endif

  b_interpol= sqrt(br**2+bz**2+bphi**2)

!  if(sml_minusB==1) then
!     br=-br
!     bz=-bz
!     bphi=-bphi
!  endif

end function b_interpol

! return B_phi
subroutine bphi_interpol_rzpsi(r,z,phi,psi,bphi)
  use eq_module 
  use sml_module, only : sml_minusB
  use rpl_module
  implicit none
  real (kind=8) , intent(in) :: r,z,phi,psi
  real (kind=8)              :: bphi,fi
  real (kind=8)              :: ripp, dripp_dr, dripp_dz
  real (kind=8) , external   :: I_interpol
 
  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol(psi,0,3)
  else
     fi=I_interpol(psi,0,1)
  endif
  
  bphi=fi / r
  if(rpl_mode==1) then
     call get_ripple(r,z,ripp,dripp_dr,dripp_dz)
     bphi=bphi*(1D0 + dcos(rpl_N_coil*phi)*ripp) 
  endif
  
  if(sml_minusB==1) then
     bphi=-bphi
  endif

end subroutine


subroutine bvec_interpol(r,z,phi,br,bz,bphi)
  use eq_module 
  use sml_module, only : sml_minusB
  use rpl_module
  implicit none
  real (kind=8) , intent(in) :: r,z,phi
  real (kind=8)              :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=8)              :: ripp, dripp_dr, dripp_dz
  real (kind=8) , external   :: I_interpol,psi_interpol
  psi     = psi_interpol(r,z,0,0)
  dpsi_dr = psi_interpol(r,z,1,0)
  dpsi_dz = psi_interpol(r,z,0,1)
  
  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol(psi,0,3)
  else
     fi=I_interpol(psi,0,1)
  endif
  
  br=- dpsi_dz / r
  bz= dpsi_dr / r
  bphi=fi / r
  if(rpl_mode==1) then
     call get_ripple(r,z,ripp,dripp_dr,dripp_dz)
     bphi=bphi*(1D0 + dcos(rpl_N_coil*phi)*ripp) 
  endif

  if(sml_minusB==1) then
     br=-br
     bz=-bz
     bphi=-bphi
  endif
end subroutine


#if defined(PSPLINE)

real (kind=8) function psi_interpol_pspline(r,z,r_der,z_der)
    use ITP_module
    use EZspline_obj
    use EZspline

    implicit none
    real (kind=8), intent(in) :: r, z
    integer , intent(in) :: r_der, z_der

    integer :: ier
    real (kind=8) :: r8value
  

#if defined(USE_EXPLICIT)
    call EZspline_interp(spl(r_der,z_der),r,z,r8value,ier)
#else
    call EZspline_derivative(spl(0,0),r_der,z_der,r,z,r8value,ier)
#endif
    call EZspline_error(ier)

    psi_interpol_pspline = r8value
end function

#endif


real (kind=8) function psi_interpol_raw(r,z,r_der,z_der)
    use ITP_module
#if !defined(IMSL)
    use bspline
#endif

    implicit none
    real (kind=8), intent(in) :: r, z
    integer , intent(in) :: r_der, z_der
#if defined(MY_IMSL)
! MY_IMSL wrapper routine
    real (kind=8) ,external:: my_db22dr
#elif defined(IMSL)
    real (kind=8) ,external:: dbs2dr
#endif

#if defined(MY_IMSL)
    psi_interpol_raw=my_db22dr(r_der, z_der,r,z,ITP_korder_rz,ITP_korder_rz,ITP_r_KNOT, &
         ITP_z_KNOT, itp_mr, itp_mz, ITP_psi_BSCOEF)
#else
    psi_interpol_raw=dbs2dr(r_der, z_der,r,z,ITP_korder_rz,ITP_korder_rz,ITP_r_KNOT, &
         ITP_z_KNOT, itp_mr, itp_mz, ITP_psi_BSCOEF)
#endif
end function


!****************************************************************************
! retrun psi and dpsi as a function of r, z
!
! first created : 2000/10/19
! last modified : 2000/10/19
!****************************************************************************

real (kind=8) function psi_interpol(r_in,z_in,r_der,z_der)
  use itp_module
  implicit none
  real (kind=8), intent(in) :: r_in, z_in
  integer , intent(in) :: r_der, z_der
  real (kind=8) :: r,z
#if defined(PSPLINE)
  real (kind=8) ,external:: psi_interpol_pspline
#else
  real (kind=8) ,external:: psi_interpol_raw
#endif


! check r,z values
  r=r_in
  z=z_in
  if( r < itp_min_r .OR. r> itp_max_r) then
     print *, 'improper r in psi_interpol (r,z,min_r,max_r):' , r,z, itp_min_r, itp_max_r
     call err_count
     r=max(itp_min_r,min(itp_max_r,r))
  endif
  if( z< itp_min_z .OR. z > itp_max_z) then
     print *, 'improper z in psi_interpol (r,z,min_z,max_z):' , r,z, itp_min_z, itp_max_z
     call err_count
     z=max(itp_min_z,min(itp_max_z,z))
  endif


#if defined(PSPLINE)
  psi_interpol=psi_interpol_pspline(r,z,r_der,z_der)
#else
  psi_interpol=psi_interpol_raw(r,z,r_der,z_der)
#endif

  
end function psi_interpol
  



!****************************************************************************
! initialization of interpolation routines
!
! first created : 2000/08/26
! last modified : 2003/01/28
! adopted from xorbit
!****************************************************************************
subroutine init_interpolation
    use eq_module
    use itp_module

#if defined(PSPLINE)
    use EZspline_obj
    use EZspline
#endif

#if !defined(IMSL)
    use bspline
#endif

    !use numerical_libraries
    implicit none
    integer :: i,j, max_mpsi, max_nearx_mr , max_mr ,region
    real (kind=8)  :: psi,theta,r,z

#if defined(PSPLINE)
    integer, parameter :: idebug = 0
    integer, dimension(2) :: BCS1, BCS2
    integer :: ier,  ir,iz,  r_der, z_der
    real (kind=8), dimension(eq_mr,eq_mz) :: rzgrid
    real (kind=8) :: rr,zz, abserr,relerr 
    real (kind=8) :: r8value, maxrelerr, maxabserr
    real (kind=8), parameter  ::  one = 1.0d0, zero = 0.0d0
    real (kind=8), external :: psi_interpol_raw
#endif

    max_mpsi=eq_mpsi
    max_mr=eq_mr

    !set valid psi, theta range for interpolation
    itp_min_psi=eq_psi_grid(1)
    itp_max_psi=eq_psi_grid(itp_mpsi)
        
    !------------------------------------------------------------

    !initializing interpolation routines
    
    ! for psi grid 
    !common - psi
!    print *,eq_psi_grid(1:2)
!    print *, 'psi_grid'
!    do i=1, eq_mpsi
!       write(111,*) i, eq_psi_grid(i)
!    enddo
!    close(111)
    
#if defined(PSPLINE)
    BCS1 = 0 ! not a knot
    BCS2 = 0 ! not a knot
    call EZspline_init(spl_psi,eq_mpsi,BCS1,ier)
    call EZspline_error(ier)

    spl_psi%x1 = eq_psi_grid

    call EZspline_setup(spl_psi,eq_I,ier)
    call EZspline_error(ier)
#elif defined(MY_IMSL)
    CALL my_db2nak(eq_mpsi,eq_psi_grid,itp_korder,itp_psi_knot)
    call dcsint(itp_mpsi,eq_psi_grid,eq_I,itp_I_break,itp_I_cscoef)
#elif defined(IMSL)
    CALL dbsnak(eq_mpsi,eq_psi_grid,itp_korder,itp_psi_knot)
    call dcsint(itp_mpsi,eq_psi_grid,eq_I,itp_I_break,itp_I_cscoef)
#else
! portable routines
    CALL dbsnak(eq_mpsi,eq_psi_grid,itp_korder,itp_psi_knot)
    call my_dcsint(itp_mpsi,eq_psi_grid,eq_I,itp_I_break,itp_I_cscoef)
#endif

    !all rz mode interpolation

    ! for r-z grid  -- additional term added for numerical safety.
!    print *, 'rz grid'
    itp_min_r=eq_rgrid(1)! *(1D0-1D-5) + 1D-5*eq_rgrid(2)
    itp_max_r=eq_rgrid(eq_mr)! *(1D0-1D-5)+ 1D-5 *eq_rgrid(eq_mr-1)
    itp_min_z=eq_zgrid(1) !*(1D0-1D-5) + 1D-5*eq_zgrid(2)
    itp_max_z=eq_zgrid(eq_mz) !*(1D0-1D-5) + 1D-5*eq_zgrid(eq_mz-1)

#if defined(MY_IMSL)
    CALL my_db2nak(eq_mr,eq_rgrid,itp_korder_rz,itp_r_knot)
    CALL my_db2nak(eq_mz,eq_zgrid,itp_korder_rz,itp_z_knot)
#else
! standard IMSL or portable bspline routine
    CALL dbsnak(eq_mr,eq_rgrid,itp_korder_rz,itp_r_knot)
    CALL dbsnak(eq_mz,eq_zgrid,itp_korder_rz,itp_z_knot)
#endif

#if defined(PSPLINE)

    BCS1 = 0 ! not a knot
    BCS2 = 0 ! not a knot
#if defined(USE_EXPLICIT)
    do r_der=0,2
    do z_der=0,2-r_der
      call EZspline_init(spl(r_der,z_der),eq_mr,eq_mz,BCS1,BCS2,ier)
      call EZspline_error(ier)
      spl(r_der,z_der)%x1 = eq_rgrid
      spl(r_der,z_der)%x2 = eq_zgrid
    enddo
    enddo
#else
     call EZspline_init(spl(0,0),eq_mr,eq_mz,BCS1,BCS2,ier)
     call EZspline_error(ier)
     spl(0,0)%x1 = eq_rgrid
     spl(0,0)%x2 = eq_zgrid
#endif

#endif


    !debug
    !do i=1, eq_mr
    !   do j=1, eq_mz
    !      write(800,800) eq_rgrid(i),eq_zgrid(j),eq_psi_rz(i,j)
    !   enddo
    !   write(800,*) ' '
    !enddo
    !do i=1, eq_mr+itp_korder_rz
    !   write(800,900) itp_r_KNOT(i)
    !enddo
    !do i=1, eq_mz+itp_korder_rz
    !   write(800,900) itp_z_KNOT(i)
    !enddo
900 format (e19.13)    
800 format (e19.13,' ',e19.13,' ',e19.13)

    !for psi

! standard IMSL or portable bspline routine
! also needed to test initialization of pspline interpolation
    CALL dbs2IN (eq_mr, eq_rgrid, eq_mz, eq_zgrid, eq_psi_rz, max_mr,  &
         itp_korder_rz, itp_korder_rz, itp_r_KNOT, itp_z_KNOT, itp_psi_BSCOEF)

#if defined(PSPLINE)

    if (idebug.ge.1) then
    write(*,*) 'mindr ',minval(eq_rgrid(2:eq_mr)-eq_rgrid(1:(eq_mr-1)))
    write(*,*) 'maxdr ',maxval(eq_rgrid(2:eq_mr)-eq_rgrid(1:(eq_mr-1)))
    write(*,*) 'mindz ',minval(eq_zgrid(2:eq_mz)-eq_zgrid(1:(eq_mz-1)))
    write(*,*) 'maxdz ',maxval(eq_zgrid(2:eq_mz)-eq_zgrid(1:(eq_mz-1)))
    endif

#if defined(USE_EXPLICIT)

    do r_der=0,2
    do z_der=0,2-r_der

      do iz=1,eq_mz
      do ir=1,eq_mr
       rr = eq_rgrid(ir)
       zz = eq_zgrid(iz)
       rzgrid(ir,iz) = psi_interpol_raw(rr,zz,r_der,z_der)
      enddo
      enddo


      call EZspline_setup(spl(r_der,z_der),rzgrid,ier)
      call EZspline_error(ier)

      !     double check

      if (idebug.ge.1) then
      maxabserr = zero
      maxrelerr = zero

      do iz=1,eq_mz-1
      do ir=1,eq_mr-1
         rr = (eq_rgrid(ir) + eq_rgrid(ir+1))/2.0d0
         zz = (eq_zgrid(iz) + eq_zgrid(iz+1))/2.0d0

         call EZspline_interp(spl(r_der,z_der),rr,zz, r8value,ier)
         call EZspline_error(ier)

         abserr = abs(r8value - psi_interpol_raw(rr,zz,r_der,z_der))
         relerr = abserr / max( one, abs(rzgrid(ir,iz)) )

         maxabserr = max( abserr, maxabserr )
         maxrelerr = max( relerr, maxrelerr )
      enddo
      enddo
      write(*,*) 'r_der,z_der ',r_der,z_der
      write(*,*) 'use_explicit:maxabserr,maxrelerr ', maxabserr,maxrelerr
      endif

    enddo
    enddo
#else
    call EZspline_setup(spl(0,0),eq_psi_rz,ier)
    call EZspline_error(ier)
#endif
    
      !   check spline evaluation of derivative

      if (idebug.ge.1) then

      do r_der=0,2
      do z_der=0,2-r_der
      maxabserr = zero
      maxrelerr = zero

      do iz=1,eq_mz-1
      do ir=1,eq_mr-1
         rr = (eq_rgrid(ir) + eq_rgrid(ir+1))/2.0d0
         zz = (eq_zgrid(iz) + eq_zgrid(iz+1))/2.0d0

         call EZspline_derivative(spl(0,0),r_der,z_der,rr,zz,r8value,ier)
         call EZspline_error(ier)

         abserr = abs(r8value - psi_interpol_raw(rr,zz,r_der,z_der))
         relerr = abserr / max( one, abs(rzgrid(ir,iz)) )

         maxabserr = max( abserr, maxabserr )
         maxrelerr = max( relerr, maxrelerr )
      enddo
      enddo
      write(*,*) 'r_der,z_der ',r_der,z_der
      write(*,*) 'pspline: maxabserr,maxrelerr ', maxabserr,maxrelerr
      enddo
      enddo

      endif
#endif
end subroutine

