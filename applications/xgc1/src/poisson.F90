subroutine poisson(grid,psn,iflag)
  use psn_class
  use grid_class
  use sml_module
  use efld_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: iflag

!  if(efld_mode==-1) then
!     return
!  endif
   
  if(.true.) then
     call poisson_two_solvers(grid,psn,iflag)
  else !(.false.) then
     !     call poisson_iter_solvers(grid,psn,iflag)
  endif
  
end subroutine poisson

subroutine init_poisson(grid,psn,iflag)
  use psn_class
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: iflag
  integer :: i
  real (kind=8) :: psi
  real (kind=8), external :: init_tempi_ev, tempe_ev
  
  ! tempe and tite is not used.
  do i=1, grid%nnode
     psi=grid%psi(i)
     psn%tempi_ev(i)=init_tempi_ev(psi)
     psn%tite(i)=psn%tempi_ev(i)/tempe_ev(psi)
     psn%tempe(i)=tempe_ev(psi)*sml_ev2j
  enddo
  
  
  if(.true.) then
     call init_poisson_two_solvers(grid,psn,iflag)
  else
     
  endif
  

  if(sml_add_pot0>0) then
     allocate( psn%add_pot0(grid%nnode))
     
     ! prepare additional 0 - potential 
     if(sml_add_pot0==1) then
        call read_add_pot0(grid,psn)
     elseif(sml_add_pot0==2) then
        call neo_pot0_simple(grid,psn)
     else
        psn%add_pot0=0D0
     endif
  endif


end subroutine init_poisson


subroutine read_add_pot0(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: num, i
  
  open(unit=18,file=sml_add_pot0_file,action='read')
  read(18,*) num
  if(num/=grid%nnode) then
     if(sml_mype==0) print *, 'Error : Grid number missmatch in read_add_pot0',num, grid%nnode
     stop
  endif
  do i=1, num
     read(18,*) psn%add_pot0(i)
  enddo

  close(18)
end subroutine read_add_pot0


subroutine poisson_two_solvers(grid,psn,iflag)
  use psn_class
  use grid_class
  use sml_module
  use perf_monitor
  use smooth_module
  use eq_module, only : eq_den_out
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: iflag
  integer :: k, loop_num, ierr
  real (kind=8) , pointer :: idensity_ad(:), phitmp(:)
  real (kind=8), allocatable :: dentmp(:), dentmp2(:), dentmp3(:), dentmp4(:)
  real (kind=8), allocatable :: sendl(:),recvr(:)
  real (kind=8) :: t1, t2, t3
  save dentmp, dentmp2, dentmp3, dentmp4, phitmp
  save sendl,recvr



  if(iflag==0) then
     allocate(dentmp(grid%nnode),dentmp2(grid%nnode),dentmp3(grid%nnode),dentmp4(grid%nnode))
     allocate(sendl(grid%nnode),recvr(grid%nnode))
     dentmp=0D0 !for safty
     dentmp2=0D0
     dentmp3=0D0
     dentmp4=0D0
     if(sml_mype==0) print *, '**************** Poisson_two_solvers is used ******************'
     return
  endif
  

  ! get n0**************************************************************************************************
  call monitor_start (PHASE0_)
!  psn%density0_full=psn%edensity0
  
!  call smooth_pol0(grid,psn%density0_full,smooth00)  
!  if(sml_deltaf==1) then
!     psn%density0_full(:)=psn%density0_full(:) + psn%f0_density0(:)
!  endif
!  psn%density0_full=max(psn%density0_full,eq_den_out)

  ! set boundary *********************************************************************************************
  call set_boundary_values(psn%idensity0,0D0,psn%cbd0)
  call set_boundary_values(psn%edensity0,0D0,psn%cbd0)
  call monitor_stop (PHASE0_)

  call monitor_start (PHASE1_)
  ! RHS density for 00 mode
  
  dentmp=psn%idensity0-psn%edensity0 ! n=0 - mode
  call smooth_pol0(grid,dentmp,smooth00)  ! making 00 mode
  if(sml_zero_out_total_charge) then  ! need opimization
     call zero_out_total_charge(grid,psn,dentmp,dentmp4) ! zero out totoal charge - dentmp4= dentmp - int(dentmp)/total_vol
  else
     dentmp4=dentmp
  endif
  call set_boundary_values(dentmp4,0D0,psn%cbd0)
  
  
#ifdef XGC_DEBUG10
        grid%rtmp4=dentmp4
#endif
  ! solve 00 (neoclassical) field****************************************************************************
  call monitor_start (SOLVER1_)
  if(sml_use_simple00) then
     call simple00(grid,psn,dentmp4,psn%pot0)
  else
     dentmp3=dentmp4 !/psn%density0_full   ! need to be optimized
     if(sml_sheath_mode/=0) call apply_wall_boundary_condition(grid,psn,dentmp3)
     call petsc_solve(grid%nnode, dentmp3, psn%pot0, psn%solver00%comm, psn%solver00, ierr)
  endif
  call monitor_stop (SOLVER1_)
  call smooth_pol0(grid,psn%pot0,smooth00)  ! ????????????????????????????? keep it generally?
  call monitor_stop (PHASE1_)

#ifdef TAVG_POT  
  t1=1D0
  t2=0.01D0
  t3=(t1-t2)*(1D0 - real(sml_istep-1)/real(10D0)) +t2
  t3=max(t3,t2)
  psn%pot0=t3*psn%pot0+(1D0-t3)*grid%rtmp8
  if(sml_ipc==2) then
     grid%rtmp8=psn%pot0
  endif
#endif

  ! solve turbulence field**********************************************************************************
  if(sml_turb_efield==1) then
     call monitor_start (PHASE3_)
     
     if(sml_turb_efield==1) then
        loop_num=grid%nphi
     else
        loop_num=1
     endif
     
     ! For all toroidal section
     ! the case when phitmp=>psn%pot0 can be eliminated. No function now.
     do k=1, loop_num
        if(sml_turb_efield==1) then
           idensity_ad=>psn%idensity(:,k)
           phitmp=>psn%dpot(:,k)
        else
           idensity_ad=>psn%idensity0(:)
           phitmp=>psn%pot0
        endif
        
        dentmp3(:)= (idensity_ad(:) - psn%edensity0(:) - dentmp(:)) !/psn%density0_full(:)  ! ni - ne - <ni - ne>
        
        call set_boundary_values(dentmp3,0D0,psn%cbdH)
        phitmp=0D0

#ifdef XGC_DEBUG10
        grid%rtmp3=dentmp3
#endif
        call monitor_start (SOLVER3_)        
        call petsc_solve(grid%nnode, dentmp3, phitmp, psn%solverH%comm, psn%solverH, ierr)     
        call monitor_stop (SOLVER3_)

        call smooth_pol0(grid,phitmp,smoothH)
!        call smooth_r(phitmp,grid%rtmp1,smooth_r1,grid)
        
     enddo

     !************************* extract 00 mode **********************************************************
     call extract_00mode(grid,psn%dpot)
                        
     !**************************** mode selection*********************************************************
     if(sml_mode_select_on==1) then
        call mode_selection(sml_mode_select_n,grid,psn)
     endif
          
     call monitor_stop (PHASE3_)

  endif
  
  
  !********************** Phase 4 *********************************************
  !send and receive potential value to and from adjacent PE
  !***************************************************************************
  if(sml_turb_efield==1) then
     call send_recv_potential( psn%dpot, recvr,sendl,grid%nnode,grid%nphi)
  endif
  
  !******************** Phase 5 ***********************************************
  ! get space derivative of potential
  ! psn%E_para(node,1:nphi )  psn%E_perp(2 vec, ntriangle, 0:nphi)
  !***************************************************************************
  if(sml_add_pot0/=0) then
     psn%pot0=psn%pot0+psn%add_pot0
  endif
#if !defined(CIRCULAR_SPECIAL)
  call get_potential_grad(grid,psn)
#else
  call get_potential_grad_cc(grid,psn)
#endif

end subroutine poisson_two_solvers

subroutine init_poisson_two_solvers(grid,psn,iflag)
  use psn_class
  use grid_class
  use sml_module
  use ptl_module, only : ptl_mass
  use eq_module, only : eq_axis_r, eq_axis_z
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: iflag
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscpc.h"
  type(mat_type) :: matl, matr, matr2
  integer :: max_mat_width
  real (kind=8),allocatable :: alpha(:), beta(:), tiev(:), teev(:), den(:),b2(:) !rhoi2_tiev(:)
  integer :: itr,nd(3),ierr,in_bd_node
  real (kind=8) :: x_center(2),psi
  real (kind=8) :: factor, psi_in_poisson_bd, psi_out_poisson_bd
  real (kind=8), external :: gyro2_tev,psi_interpol, init_tempi_ev_wz, tempe_ev_wz, b_interpol, init_den_wz
  real (kind=8), parameter :: offset=2D0
  allocate( alpha(grid%ntriangle), beta(grid%ntriangle),tiev(grid%ntriangle),teev(grid%ntriangle),den(grid%ntriangle) ,b2(grid%ntriangle) )  ! ,rhoi2_tiev(grid%ntriangle) )
  max_mat_width=sml_max_mat_width
  call new_mat(matl,grid%nnode,max_mat_width)
  call new_mat(matr,grid%nnode,max_mat_width)
  call new_mat(matr2,grid%nnode,max_mat_width)
  

  ! setup basic values
  do itr=1, grid%ntriangle
     nd=grid%nd(:,itr)
     x_center=( grid%x(:,nd(1)) + grid%x(:,nd(2)) + grid%x(:,nd(3)) )/3D0
#ifdef FLAT_GYRORADIUS
     x_center(1)=eq_axis_r
     x_center(2)=eq_axis_z
#endif
!     rhoi2_tiev(itr)=gyro2_tev(x_center)
     psi=psi_interpol(x_center(1),x_center(2),0,0)
     tiev(itr)=init_tempi_ev_wz(psi,x_center(2))
     teev(itr)=tempe_ev_wz(psi,x_center(2))
     b2(itr)=b_interpol(x_center(1),x_center(2),0D0)**2
     den(itr)=init_den_wz(psi,x_center(2))
  enddo
     
  !Cold electron on the poisson boundary
  if(sml_bd_Te_mode/=0) then
     !get boundary psi value of both boundary
     if(sml_bd_Te_mode==1 .or. sml_bd_Te_mode==2) then
        in_bd_node=max(1,psn%pbd0%in%end )
!        x_center=grid%x(:,grid%inv_sort(in_bd_node))
        x_center=grid%x(:,in_bd_node)
        psi_in_poisson_bd = psi_interpol(x_center(1),x_center(2),0,0)
     endif
     if(sml_bd_Te_mode==2 .or. sml_bd_Te_mode==3) then
!        x_center=grid%x(:,grid%inv_sort(grid%outer_bd(1)))
        x_center=grid%x(:,psn%pbd0%out1%start )
        psi_out_poisson_bd= psi_interpol(x_center(1),x_center(2),0,0)
     endif

     do itr=1, grid%ntriangle
        nd=grid%nd(:,itr)
        x_center=( grid%x(:,nd(1)) + grid%x(:,nd(2)) + grid%x(:,nd(3)) )/3D0        
        psi=psi_interpol(x_center(1),x_center(2),0,0)
        
        factor=1D0
        ! inside
        if(sml_bd_Te_mode==1 .or. sml_bd_Te_mode==2) then
           factor=0.5D0*(1D0+tanh( (psi - psi_in_poisson_bd)/sml_bd_Te_width*4D0 - offset ))
        endif
        !outside
        if(sml_bd_Te_mode==2 .or. sml_bd_Te_mode==3) then
           factor=factor*0.5D0*(1D0+tanh( -(psi - psi_out_poisson_bd)/sml_bd_Te_width*4D0 - offset ))
        endif
        factor=max(factor,0.5D0*(1D0+tanh(-offset)))
        teev(itr)=teev(itr)*factor
     enddo
  endif



  ! 00 solver
  if(.NOT. sml_use_simple00 ) then
     if(sml_fem_matrix==1) then
        ! LHS
        ! Alpha and Beta
        alpha=-ptl_mass(1)/sml_e_charge*den/b2
        beta=0D0
        call form_fem_matrix(matl, grid%ntriangle, grid%nd,grid%nnode, grid%x, alpha, beta, ierr)
        ! RHS
        psn%solver00%n_rhs_mat=1
        if(sml_use_pade) then
           print *, 'Warning Pade approximation is not working now'           
           alpha=0D0
        else
           alpha=0D0  ! no del term on RHS
        endif
        beta=1D0  ! diagonal (identity)
        call form_fem_matrix(matr, grid%ntriangle, grid%nd,grid%nnode, grid%x, alpha, beta, ierr)
        !
     else ! integral method
        !
        psn%solver00%n_rhs_mat=0
!        call form_integral_matrix(grid,psn,matl)        
        
     endif
     
     ! **************************************
     ! initialize petsc_solver from my matrix
     ! **************************************
     psn%solver00%comm = petsc_comm_world
     psn%solver00%mype=sml_mype
     psn%solver00%totalpe=sml_totalpe
     psn%solver00%prefix=1
     call initialize_solver(grid,psn,psn%solver00,matl,matr,matr2)
     call del_mat(matl)
     call del_mat(matr)
     call del_mat(matr2)
  endif
  
  ! turb solver
  if(sml_turb_efield==1) then
     max_mat_width=sml_max_mat_width
     call new_mat(matl,grid%nnode,max_mat_width)
     call new_mat(matr,grid%nnode,max_mat_width)
     call new_mat(matr2,grid%nnode,max_mat_width)

     if(sml_fem_matrix==1) then 
        ! LHS
        if(sml_use_pade) then
           print *, 'Warning Pade approximation is not working now (turb sovler)'
!           alpha=-(tiev/teev+1D0)*rhoi2_tiev
!           beta=1D0/teev
           alpha=-ptl_mass(1)/sml_e_charge*den/b2
           beta=den/teev
        else
           alpha=-ptl_mass(1)/sml_e_charge*den/b2
           beta=den/teev
!           alpha=-rhoi2_tiev
!           beta=1D0/teev
        endif
        call form_fem_matrix(matl, grid%ntriangle, grid%nd,grid%nnode, grid%x, alpha, beta, ierr)
        
        psn%solverH%n_rhs_mat=1     
        if(sml_use_pade) then
!           alpha=-rhoi2_tiev
           alpha=0D0
        else
           alpha=0D0  ! no del term on RHS
        endif
        beta=1D0  ! diagonal (identity)
        call form_fem_matrix(matr, grid%ntriangle, grid%nd,grid%nnode, grid%x, alpha, beta, ierr)        
        
!        if(psn%solverH%n_rhs_mat>=2) then
!           alpha=
!           beta=
!           call form_fem_matrix(matr2,grid%ntriangle, grid%nd,grid%nnode, grid%x, alpha, beta, ierr)
!        endif
     else ! integral method
        psn%solverH%n_rhs_mat=0
        !        call form_int_matrix(matl,        
     endif
     ! **************************************
     ! initialize petsc_solver from my matrix
     ! **************************************
     psn%solverH%comm = sml_plane_comm
     psn%solverH%mype=sml_plane_mype
     psn%solverH%totalpe=sml_plane_totalpe
     psn%solverH%prefix=2
     call initialize_solver(grid,psn,psn%solverH,matl,matr,matr2)

     call del_mat(matl)
     call del_mat(matr)
     call del_mat(matr2)
  endif
  
  deallocate(alpha,beta,tiev,teev,den,b2)

end subroutine init_poisson_two_solvers

subroutine form_fem_matrix( Amat, iemax, indices, npts, & ! iemax = ntriangle, indices = nd(3,nt), npts=nnode
     points, alpha, beta, ierr )              ! points=x(2,nn), 
  !========================================================================
  !c
  !c
  !c     input:
  !c     iemax: number of elements
  !c     iet1,2,3[iemax]: nodes of triangles of mesh
  !c     npts: number of mesh points
  !c     points[2,npts]: coordinates of points
  !c
  !c     output:
  !c     Amat: matrix to add tangent to
  !c     ierr: error code
  !c
  use mat_class
  implicit none
  integer indices(3,iemax),iemax,ierr,npts
  real (kind=8) ::  points(2,npts),factor(npts)
  type(mat_type) :: Amat
  real (kind=8) :: alpha(iemax),beta(iemax)
  real (kind=8) ::  c1,c2,xlt(2,3),arr(3,3),arr2(3*3),ul(3,4), rhoi_th, x_center(2),den,psi
  integer ie,i,j,k,kk,ind(3),nel,nst
  integer :: itmp,jtmp

  nel = 3                   ! points per element
  nst = 3                   ! size of elem matrix
  ul  = 0D0                   ! displacement, velocit and acceloration (not used)
  !c2    = 0.d0             ! mass (beta)

!  call MatZeroEntries(Amat,ierr)

  do ie=1,iemax
     ind(:) = indices(:,ie)
     if(ind(1) .le. 0) stop 'ind(1) == 0'
     
     xlt(1,1) = points(1,ind(1))
     xlt(1,2) = points(1,ind(2))
     xlt(1,3) = points(1,ind(3))
     xlt(2,1) = points(2,ind(1))
     xlt(2,2) = points(2,ind(2))
     xlt(2,3) = points(2,ind(3))
     
     !obtain gyro_radius
     !taus= 1.
     !     c1=1
     c1 = -alpha(ie)  ! del^2
     c2 = beta(ie)   ! mass
     arr = 0D0
     call therm2d(c1,c2,ul,xlt,arr,nel,nst,2,1)
     ! zero based
     ind = ind - 1

!     call MatSetValues(Amat,nel,ind,nel,ind,arr,ADD_VALUES,ierr)  ! petsc matrix
     do itmp=1, nel
        do jtmp=1, nel
           call set_value(Amat,ind(itmp)+1,ind(jtmp)+1,arr(itmp,jtmp),1)  !????? check i,j order of arr
        enddo
     enddo
  enddo

end subroutine form_fem_matrix

subroutine initialize_solver(grid,psn,solver, matl, matr, matr2)
  use grid_class
  use psn_class
  use mat_class
  use xgc_solver_module
  use sml_module, only : sml_mype, sml_boundary_diagonal, sml_sheath_mode
!  use boundary_class 
  implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"  
  type(grid_type) :: grid
  type(psn_type):: psn
  type(mat_type) :: matl, matr, matr2
  type(xgc_solver) :: solver
  real (kind=8), parameter :: epsilon=0D0
  real (kind=8) :: diagonal
  integer :: maxn1, maxn2, maxn3, ierr
  integer :: i,j,k
  type(boundary_type) ::  pbd,cbd
!  logical, external :: is_inside
  
  diagonal = sml_boundary_diagonal
  if(solver%prefix==1) then ! solver0
     pbd=psn%pbd0
     cbd=psn%cbd0
  else ! solverH
     pbd=psn%pbdH
     cbd=psn%cbdH
  endif 


  call get_max_width(matl , maxn1)
  call get_max_width(matr , maxn2)
  call get_max_width(matr2 , maxn3)
  maxn1=max(maxn1,1)
  maxn2=max(maxn2,1)
  maxn3=max(maxn3,1)
  call new_solver( solver, grid%nnode, maxn1, maxn2, maxn3, ierr )

  ! 00 solver  - apply outer boundary condition only
  do i=1, matl%n
     do j=1, matl%nelement(i)
        k=matl%eindex(j,i)
        if( is_inside(i,pbd) .and. (is_inside(k,pbd) .or. sml_sheath_mode/=0) ) then
           call MatSetValues(solver%mat,1,i-1,1,k-1,matl%value(j,i),ADD_VALUES,ierr)
           CHKERRQ(ierr)
           call MatSetValues(solver%mat,1,k-1,1,i-1,epsilon,ADD_VALUES,ierr)
           CHKERRQ(ierr)
        endif
     enddo
     if( .NOT. is_inside(i,pbd) ) then
        call MatSetValues(solver%mat,1,i-1,1,i-1,diagonal,ADD_VALUES,ierr)
        CHKERRQ(ierr)
     endif
  enddo

  if(solver%n_rhs_mat >=1) then
     do i=1, matr%n
        do j=1, matr%nelement(i)
           k=matr%eindex(j,i)
           if( is_inside(i,cbd) .and. is_inside(k,cbd) ) then
              call MatSetValues(solver%rhs_mat,1,i-1,1,k-1,matr%value(j,i),ADD_VALUES,ierr)
              CHKERRQ(ierr)
              call MatSetValues(solver%rhs_mat,1,k-1,1,i-1,epsilon,ADD_VALUES,ierr)
              CHKERRQ(ierr)
           endif
        enddo
        if( .NOT. is_inside(i,cbd) ) then
           call MatSetValues(solver%rhs_mat,1,i-1,1,i-1,diagonal,ADD_VALUES,ierr)
           CHKERRQ(ierr)
        endif
     enddo
  endif

  if(solver%n_rhs_mat >=2) then
     do i=1, matr2%n
        do j=1, matr2%nelement(i)
           k=matr2%eindex(j,i)
           if( is_inside(i,cbd) .and. is_inside(k,cbd) ) then
              call MatSetValues(solver%rhs2_mat,1,i-1,1,k-1,matr2%value(j,i),ADD_VALUES,ierr)
              CHKERRQ(ierr)
              call MatSetValues(solver%rhs2_mat,1,k-1,1,i-1,epsilon,ADD_VALUES,ierr)
              CHKERRQ(ierr)
           endif
        enddo
!        if( .NOT. is_inside(i,grid,2) ) then
!           call MatSetValues(solver%rhs2_mat,1,i-1,1,i-1,diagonal,ADD_VALUES,ierr)
!           CHKERRQ(ierr)
!        endif
     enddo
  endif

 ! Matrix assembly
  call MatAssemblyBegin( solver%mat, MAT_FINAL_ASSEMBLY, ierr )
  CHKERRQ(ierr)
  call MatAssemblyEnd( solver%mat, MAT_FINAL_ASSEMBLY, ierr )
  CHKERRQ(ierr)
  if(solver%n_rhs_mat>=1) then
     call MatAssemblyBegin( solver%rhs_mat, MAT_FINAL_ASSEMBLY, ierr )
     CHKERRQ(ierr)
     call MatAssemblyEnd( solver%rhs_mat, MAT_FINAL_ASSEMBLY, ierr )
     CHKERRQ(ierr)
  endif
  if(solver%n_rhs_mat>=2) then
     call MatAssemblyBegin( solver%rhs2_mat, MAT_FINAL_ASSEMBLY, ierr )
     CHKERRQ(ierr)
     call MatAssemblyEnd( solver%rhs2_mat, MAT_FINAL_ASSEMBLY, ierr )
     CHKERRQ(ierr)  
  endif


  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call parallelize_solver(solver, grid%nnode, solver%comm,&
       solver%mype, sml_mype, solver%totalpe, ierr )


!  if(sml_turb_efield==1 .or. sml_smtrc_ad_elec==1) then
!     if(sml_mype==0) print *, 'Parallelize solverH'
!     call parallelize_solver(psn%solverH, grid%nnode, sml_plane_comm,&
!          sml_plane_mype, sml_mype, sml_plane_totalpe, ierr )
!  endif

  if(sml_mype==0) print *, 'End of init_poisson'

end subroutine initialize_solver

subroutine parallelize_solver( solver, neq, comm, mype_loc, mype, npe_loc, ierr )
  use xgc_solver_module
  use sml_module, only : sml_totalpe
  implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscviewer.h"
  integer          neq, comm, ierr, mype_loc, mype, npe_loc
  type(xgc_solver) :: solver
  Mat              Bmat,Amat
  Vec              vec
  double precision tol,value(17),zeros(17)
  integer          i,flg,j,col(17),inz,n_loc,start,kk
  PetscViewer      viewer
  MatPartitioning  mpart
  IS               mis,isn,id_is,xgc_petsc_is
  PetscInt         cnt(sml_totalpe) ! the size should be sml_total_pe in general

  !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Amat = solver%mat
  if( mype / npe_loc .eq. -1) then
     call PetscViewerASCIIOpen(comm, 'A.m',viewer,ierr)
     call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
     if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)
     call MatView( Amat, viewer, ierr )
     if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)
     call PetscViewerDestroy(viewer,ierr)
     !         call PetscViewerDrawOpen( comm, PETSC_NULL_CHARACTER,PETSC_NULL_CHARACTER,0,0,700,700,viewer,ierr )
     !         tol = 1.
     !         call VecSet( solver%xVec, tol, ierr )
     !     call VecSetRandom( solver%xVec, PETSC_NULL, ierr )
     !         call MatMult( Amat, solver%xVec, solver%bVec, ierr )
     !         call VecView( solver%bVec, viewer, ierr )
     !         call PetscSleep( 5, ierr )
     !         call PetscViewerDestroy(viewer,ierr)
  endif

  call check_point(1)
  call MatPartitioningCreate( comm, mpart, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(2)
  call MatPartitioningSetAdjacency( mpart, Amat, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(3)
  call MatPartitioningSetFromOptions( mpart, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(4)
  call MatPartitioningApply( mpart, mis, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(5)
  call MatPartitioningDestroy( mpart, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(6)
  call ISPartitioningToNumbering( mis, xgc_petsc_is, ierr )
  call check_point(7)
  call ISPartitioningCount( mis, cnt, ierr )
  call check_point(8)
  call ISDestroy(mis,ierr)
  call check_point(9)
  !        side effect
  call check_point(10)
  call ISInvertPermutation( xgc_petsc_is, cnt(mype_loc+1),&
       solver%petsc_xgc_is, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(11)
  call ISSort( solver%petsc_xgc_is, ierr )
  call check_point(12)
  call ISDestroy( xgc_petsc_is, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(13)
  call ISAllGather( solver%petsc_xgc_is, isn, ierr )
  call check_point(14)
  call MatGetSubMatrix(Amat,solver%petsc_xgc_is,isn,PETSC_DECIDE,&
       MAT_INITIAL_MATRIX,Bmat,ierr)
  call check_point(15)
  call MatDestroy( Amat, ierr )
  solver%mat = Bmat
  Amat = Bmat

  ! RHS
  call check_point(16)
  if(solver%n_rhs_mat>=1) then
     Bmat = solver%rhs_mat
     call Matgetsubmatrix(Bmat,solver%petsc_xgc_is,isn,PETSC_DECIDE,MAT_INITIAL_MATRIX,solver%rhs_mat,ierr)
     call MatDestroy( Bmat, ierr )
  endif
  call check_point(17)
  call ISDestroy( isn, ierr )
  call check_point(18)
  ! allocate: solver%scat solver%xVec solver%bVec
  call MatGetOwnershipRange( Amat, start, kk, ierr )
  n_loc = kk - start
  call check_point(19)
  call VecCreateMPI( comm, n_loc, neq, solver%xVec, ierr )
  call check_point(20)
  call VecSetFromOptions( solver%xVec, ierr )
  call check_point(21)
  call VecDuplicate( solver%xVec, solver%bVec, ierr )
  call check_point(22)

  call ISCreateStride(PETSC_COMM_SELF,n_loc,start,1,id_is,ierr)
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  !     scatter object
  call check_point(23)
  call VecCreateSeq( PETSC_COMM_SELF, neq, vec, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  ! side effect
  call check_point(24)
  call VecScatterCreate(vec,solver%petsc_xgc_is,solver%bVec,&
       id_is,solver%scat,ierr)
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )

  call check_point(25)
  call ISDestroy( id_is, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(26)
  call VecDestroy( vec, ierr )
  if(ierr .ne. 0) call MPI_ABORT( comm, 1, ierr )
  call check_point(27)
  call KSPSetOperators(solver%ksp, Amat, Amat, SAME_NONZERO_PATTERN, ierr)
  CHKERRQ(ierr)
  
end subroutine parallelize_solver

subroutine check_point(n)
  use sml_module
  implicit none
  integer :: n,ierr
  include 'mpif.h'

  call mpi_barrier(MPI_COMM_WORLD,ierr)  
  if(sml_mype==0) print *, 'C',n
end subroutine check_point


!cccccccccccccccccsubroutine petsc_solve cccccccccccccccccccccccccccccc
subroutine petsc_solve( nn, bb, xx, comm, solver, ierr )
!c     solve linear system.
!c     Solver (KSP) and matrix are in linear_solve.h
!c
!c     input:
!c       nn: global number of equations in system
!c       bb: rhs[nn]
!c       xx: solution[nn]
!c       comm: communicator for sysetem (ie, a poloidal plane)
!c
!c     output:
!c       bb: RHS[nn]
!c     ierr: error code
!c
  use xgc_solver_module
  implicit none
  integer nn, comm, ierr
  PetscScalar xx(nn),bb(nn),yy(nn)
  type(xgc_solver) :: solver
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscpc.h"
  integer          iplot,i,its,j
#include "include/finclude/petscviewer.h"
  Vec              vec_yy   ! line 789
  !PetscViewer  viewer

  if(solver%mat.eq.PETSC_NULL_OBJECT) call MPI_ABORT(comm,1,ierr)

  !     scatter rhs into PETSc vector
  call VecCreateSeqWithArray(PETSC_COMM_SELF,nn,bb,vec_yy,ierr)
  if(ierr .ne. 0) then
     print *, 'Error in petsc_solve - VecCreateSeqWithArray :',ierr
     call MPI_ABORT(comm,1,ierr)
  endif

  !     call VecSet( solver%bVec, 1.e30, ierr ) ! debug
#ifdef PETSC_OLD
  call VecScatterBegin( vec_yy, solver%bVec, INSERT_VALUES, SCATTER_FORWARD, solver%scat, ierr )
#else
  call VecScatterBegin( solver%scat, vec_yy, solver%bVec, INSERT_VALUES, SCATTER_FORWARD, ierr )
#endif
  if(ierr .ne. 0) then
     print *, 'Error in petsc_solve - VecScatterBegin :',ierr
     call MPI_ABORT(PETSC_COMM_WORLD,1,ierr)
  endif
#ifdef PETSC_OLD
  call VecScatterEnd( vec_yy, solver%bVec, INSERT_VALUES, SCATTER_FORWARD, solver%scat, ierr )
#else
  call VecScatterEnd( solver%scat, vec_yy, solver%bVec, INSERT_VALUES, SCATTER_FORWARD, ierr )
#endif
  if(ierr .ne. 0) then
     print *, 'Error in petsc_solve - VecScatterEnd :',ierr
     call MPI_ABORT( comm, 1, ierr )
  endif


  !debug
  !call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",viewer,ierr);
  !  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'mat.output',1,viewer,ierr);
  !call MatView(solver%mat, viewer,ierr)
  !call VecView(solver%bVec, viewer,ierr)
  !call PetscViewerDestroy(viewer,ierr)
  !  stop

!  if( solver%use_mass_matrix ) then
  if( solver%n_rhs_mat>=1 ) then
     ! bb <- (I + A) * bb or I * bb
     call MatMult( solver%rhs_mat, solver%bVec, solver%xVec, ierr )
     if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)
     call VecCopy( solver%xVec, solver%bVec, ierr )
     if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)
  end if

  if( solver%n_rhs_mat>=2 ) then
     !Not ready 
     ! add one more vector to solve
  endif

  !     Solve
  call KSPSolve( solver%ksp, solver%bVec, solver%xVec, ierr )
  if(ierr .ne. 0) then
     print *, 'Error in petsc_solve - KSPSolve :',ierr
     call MPI_ABORT(comm,1,ierr)
  endif

  !     scatter solution to xgc array
  call VecDestroy( vec_yy, ierr )
  call VecCreateSeqWithArray(PETSC_COMM_SELF,nn,yy,vec_yy,ierr)
  if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)

  call VecSet( vec_yy, 0.d0, ierr )
  if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)
#ifdef PETSC_OLD
  call VecScatterBegin( solver%xVec, vec_yy, INSERT_VALUES, SCATTER_REVERSE, solver%scat, ierr )
#else
  call VecScatterBegin( solver%scat, solver%xVec, vec_yy, INSERT_VALUES, SCATTER_REVERSE, ierr )
#endif
  if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)
#ifdef PETSC_OLD
  call VecScatterEnd( solver%xVec, vec_yy, INSERT_VALUES, SCATTER_REVERSE, solver%scat, ierr )
#else
  call VecScatterEnd( solver%scat, solver%xVec, vec_yy, INSERT_VALUES, SCATTER_REVERSE, ierr )
#endif
  if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)

  !     add all to all -- not super efficient
  call MPI_ALLREDUCE(yy,xx,nn,MPI_DOUBLE_PRECISION, MPI_SUM,comm,ierr)
  ! clean up
  call VecDestroy( vec_yy, ierr )
  if(ierr .ne. 0) call MPI_ABORT(comm,1,ierr)

end subroutine petsc_solve



! send and receive other plane potential
subroutine send_recv_potential( dpot, recvr,sendl,nnode,nphi)
  use sml_module
  use perf_monitor
  implicit none
  include 'mpif.h'
  integer , intent(in) :: nnode, nphi
  real (kind=8) :: recvr(nnode), sendl(nnode), dpot(nnode,-1:nphi+1)  
  integer :: icount, idest,isource,isendtag,irecvtag,ierr
  integer :: istatus

  call monitor_start (PHASE4_)
  
  ! receive 0
  sendl=dpot(:,nphi)
  recvr=0D0
  icount=nnode
  isource= modulo(sml_mype-sml_pe_per_plane,sml_totalpe)
  idest  = modulo(sml_mype+sml_pe_per_plane,sml_totalpe)
  isendtag=sml_mype
  irecvtag=isource
  
  call mpi_sendrecv(sendl,icount,MPI_DOUBLE_PRECISION,idest,isendtag,&
       recvr,icount,MPI_DOUBLE_PRECISION,isource,irecvtag,MPI_COMM_WORLD,istatus,ierr)

  dpot(:,0) = recvr(:)

  ! receive -1
  sendl=dpot(:,nphi-1)
  recvr=0D0
  icount=nnode
  isource= modulo(sml_mype-sml_pe_per_plane,sml_totalpe)
  idest  = modulo(sml_mype+sml_pe_per_plane,sml_totalpe)

  isendtag=sml_mype
  irecvtag=isource

  call mpi_sendrecv(sendl,icount,MPI_DOUBLE_PRECISION,idest,isendtag,&
       recvr,icount,MPI_DOUBLE_PRECISION,isource,irecvtag,MPI_COMM_WORLD,istatus,ierr)

  dpot(:,-1) = recvr(:)

  ! receive nphi+1
  sendl=dpot(:,1)
  recvr=0D0
  icount=nnode
  isource= modulo(sml_mype+sml_pe_per_plane,sml_totalpe)
  idest  = modulo(sml_mype-sml_pe_per_plane,sml_totalpe)
  isendtag=sml_mype
  irecvtag=isource
  call mpi_sendrecv(sendl,icount,MPI_DOUBLE_PRECISION,idest,isendtag,&
       recvr,icount,MPI_DOUBLE_PRECISION,isource,irecvtag,MPI_COMM_WORLD,istatus,ierr)

  dpot(:,nphi+1) = recvr(:)
  call monitor_stop (PHASE4_)


end subroutine send_recv_potential

subroutine get_potential_grad(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  use smooth_module
  use perf_monitor
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,iphi,itr,nodes(3),j
  real (kind=8) :: dp1,dp2
  integer :: iphi_p, iphi_m
  real (kind=8) :: dpot_p, dpot_m
  integer :: nodes_p(3),nodes_m(3),nd,sgn
  real (kind=8) :: area_sum, E0(2), E(2,0:grid%nphi)

  call monitor_start (PHASE5_)
  ! store grad-2d phi0  in E_perp0
  do itr=1, grid%ntriangle
     nodes(:)=grid%nd(:,itr)
     dp1=psn%pot0(nodes(1))-psn%pot0(nodes(3)) 
     dp2=psn%pot0(nodes(2))-psn%pot0(nodes(3)) 
     psn%E_perp0_tr(:,itr)= -( dp1*grid%mapping(1,1:2,itr) + dp2*grid%mapping(2,1:2,itr) )        
  enddo
  ! store grad-2d phi0 + phi in E_perp

  if(sml_turb_efield==1) then
     do iphi=0, grid%nphi
        do itr=1, grid%ntriangle
           nodes(:)=grid%nd(:,itr)
           if(sml_no_00_efield==0) then
              dp1=psn%pot0(nodes(1))-psn%pot0(nodes(3)) & 
                   + psn%dpot(nodes(1),iphi)- psn%dpot(nodes(3),iphi)
              dp2=psn%pot0(nodes(2))-psn%pot0(nodes(3)) &
                   + psn%dpot(nodes(2),iphi)- psn%dpot(nodes(3),iphi)
           else
              dp1= psn%dpot(nodes(1),iphi)- psn%dpot(nodes(3),iphi)
              dp2= psn%dpot(nodes(2),iphi)- psn%dpot(nodes(3),iphi)
           endif
           psn%E_perp_tr(:,itr,iphi)= -( dp1*grid%mapping(1,1:2,itr) + dp2*grid%mapping(2,1:2,itr) )        
        enddo
     enddo
  endif

  ! Get area averge of perp E-field
  if(sml_turb_efield==1) then
     do i=1, grid%nnode
        area_sum=0D0     
        E=0D0
        do j=1, grid%num_t_node(i)
           itr=grid%tr_node(j,i)
           area_sum=area_sum+grid%tr_area(itr)
           E(:,:)=E(:,:)+ psn%E_perp_tr(:,itr,0:grid%nphi) * grid%tr_area(itr)
        enddo
        psn%E_perp_node(:,i,:)=E(:,:)/area_sum
     enddo
  else
     do i=1, grid%nnode
        area_sum=0D0     
        E0=0D0
        do j=1, grid%num_t_node(i)
           itr=grid%tr_node(j,i)
           area_sum=area_sum+grid%tr_area(itr)
           E0(:)=E0(:)  + psn%E_perp0_tr(:,itr) * grid%tr_area(itr)
        enddo
        psn%E_perp0_node(:,i)=E0(:)/area_sum
     enddo
  endif
  
  
  if(sml_bt_sign>0) then
     sgn=1
  else
     sgn=-1
  endif

  ! store grad_para (phi0 + phi) in E_para
  if(sml_turb_efield==1) then
#if !defined(LINEAR_E_PARA)
     do i=1, grid%nnode
        do iphi=0, grid%nphi
           iphi_p=iphi+sgn
           iphi_m=iphi-sgn
           nodes_p=grid%nd(:,psn%bfollow_tr(1,i))
           nodes_m=grid%nd(:,psn%bfollow_tr(2,i))
           dpot_p=0D0
           dpot_m=0D0
           do nd=1,3
              dpot_p=dpot_p + psn%dpot(nodes_p(nd),iphi_p)*psn%bfollow_p(nd,1,i)
              dpot_m=dpot_m + psn%dpot(nodes_m(nd),iphi_m)*psn%bfollow_p(nd,2,i)
           enddo
           psn%E_para(i,iphi) = -(dpot_p-dpot_m)*psn%bfollow_1_dx(i)    !(psn%pot0(i)
        enddo
     enddo
#else
     ! particles are between plane 0 and plane 1 when nphi==1
     ! sign = -1 case only
     if(sgn/=-1) then
        print *, 'Not implimented yet - get_potential_grad.'
        stop
     endif
     ! E-field between plane 0 and 1 is assigned E_para(:,1)  
     ! E_para = - (dpot_p - dpot_current)     
     do i=1, grid%nnode
        do iphi=1, grid%nphi
           iphi_p=iphi+sgn
           !           iphi_m=iphi-sgn
           nodes_p=grid%nd(:,psn%bfollow_tr(1,i))
           !           nodes_m=grid%nd(:,psn%bfollow_tr(2,i))
           dpot_p=0D0
           dpot_m=psn%dpot(i,iphi)
           !           dpot_m=0D0
           do nd=1,3
              dpot_p=dpot_p + psn%dpot(nodes_p(nd),iphi_p)*psn%bfollow_p(nd,1,i)
              !              dpot_m=dpot_m + psn%dpot(nodes_m(nd),iphi_m)*psn%bfollow_p(nd,2,i)
           enddo
           psn%E_para(i,iphi) = -(dpot_p-dpot_m)*psn%bfollow_1_dx(i)*2D0    !(psn%pot0(i)
        enddo
     enddo
#endif
  endif
  call monitor_stop (PHASE5_)
  
end subroutine get_potential_grad
subroutine zero_out_total_charge(grid,psn,den_org,den_zero)
  use grid_class
  use psn_class
  implicit none
!  integer :: flag
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8) :: den_org(grid%nnode), den_zero(grid%nnode), total,area
  integer :: ibd, obd
  
  ibd=psn%cbd0%in%end+1
  obd=psn%cbd0%out1%start-1
  
  ! Volume?? Area??
!  total=sum( den_org(grid%inv_sort(ibd:obd))*grid%node_area(grid%inv_sort(ibd:obd)) )
  total=sum( den_org(ibd:obd)*grid%node_area(ibd:obd) )
!  area=sum(grid%node_area(grid%inv_sort(ibd:obd)))
  area=sum(grid%node_area(ibd:obd))
  
  den_zero=den_org - total/area

end subroutine zero_out_total_charge

subroutine set_decaying_boundary(den,grid)
  use grid_class
  use sml_module
  implicit none
!  integer :: flag
  type(grid_type) :: grid
  real (kind=8) :: den(grid%nnode)
  real (kind=8) :: sml_f0_psi_c, sml_f0_1_psi_w, alpha
  integer :: i

  sml_f0_psi_c=0.5*(sml_inpsi+sml_outpsi)
  sml_f0_1_psi_w=1D0/( 0.35*(sml_outpsi-sml_inpsi) )

  do i=1, grid%nnode
    alpha= exp( - ((grid%psi(i) - sml_f0_psi_c)*sml_f0_1_psi_w )**6 )
    den(i)=den(i)* alpha
  enddo

end subroutine


!! return gyro radius with a give mu and position
real (kind=8) function gyro_radius(x,mu)
  use ptl_module
  implicit none
  real (kind=8),intent(in) :: x(2),mu
  real (kind=8) :: B
  real (kind=8),external :: b_interpol

  B=b_interpol(x(1),x(2),0D0)
  gyro_radius=sqrt(mu/B/ptl_c2_2m(1))


end function gyro_radius


!!return gyro radius of a thermal ion at a given position
real (kind=8) function gyro_radius2(x)
  use sml_module
  use ptl_module
  implicit none
  real (kind=8),intent(in) :: x(2)
  real (kind=8) :: B,psi,tev
  real (kind=8),external :: b_interpol,psi_interpol
  real (kind=8),external :: init_tempi_ev

  B=b_interpol(x(1),x(2),0D0)
  psi=psi_interpol(x(1),x(2),0,0)
  tev=init_tempi_ev(psi)

  gyro_radius2=sqrt(Tev*sml_ev2j/ptl_c2_2m(1)/2D0)/B
  

end function gyro_radius2

real (kind=8) function gyro2_tev(x)
  use sml_module
  use ptl_module
  implicit none
  real (kind=8),intent(in) :: x(2)
  real (kind=8) :: B,psi,tev
  real (kind=8),external :: b_interpol,psi_interpol
  real (kind=8),external :: init_tempi_ev

  B=b_interpol(x(1),x(2),0D0)
  psi=psi_interpol(x(1),x(2),0,0)
!  tev=init_tempi_ev(psi)

!  gyro2_tev=sml_ev2j/ptl_c2_2m(1)/2D0/B**2
  gyro2_tev=sml_ev2j/(ptl_c2_2m(1)*2D0*B**2)

end function gyro2_tev

!! Write out matrix element for poisson equation
!subroutine output_matrix(grid,psn)
!  use psn_class
!  use grid_class
!  implicit none
!  integer :: i,j
!  type(psn_type) :: psn
!  type(grid_type) :: grid
!  
!  open(unit=201,file='poisson.matrix',status='replace')
!  do i=1, grid%nnode
!     do j=1, psn%mat%nelement(i)
!        write(201,*) i,psn%mat%eindex(j,i),psn%mat%value(j,i),psn%mat%value(j,i)/grid%node_area(i)
!     enddo
!  enddo
!  close(201)
!
!end subroutine output_matrix

! return .true. if ith node point is inside of simulation boundary (extended)  
! type == 1 check only outer boundary
! type == 2 check both inner and outer boundary 
!!$logical function is_inside(usi, grid, type) 
!!$  use grid_class  
!!$  use sml_module, only : sml_zero_inner_bd
!!$  implicit none
!!$  integer, intent(in) :: usi
!!$  integer, intent(in) :: type
!!$  type(grid_type) :: grid
!!$  integer :: i
!!$  
!!$  i=grid%sort(usi)
!!$
!!$  if(type==1 .and. sml_zero_inner_bd==0 ) then
!!$     if( (grid%outer_bd(1) <= i .and. i<= grid%outer_bd(2)) .or. &
!!$          (grid%outer_bd(3) <= i .and. i<= grid%outer_bd(4)) ) then
!!$        is_inside=.false.
!!$     else
!!$        is_inside=.true.
!!$     endif
!!$
!!$  else !type 2
!!$     if( (grid%inner_bd(1)<= i .and. i<= grid%inner_bd(2)) .or. & 
!!$          (grid%outer_bd(1) <= i .and. i<= grid%outer_bd(2)) .or. &
!!$          (grid%outer_bd(3) <= i .and. i<= grid%outer_bd(4)) ) then 
!!$        is_inside=.false.
!!$     else
!!$        is_inside=.true.
!!$     endif
!!$  endif
!!$  
!!$end function is_inside
!!$
!!$subroutine set_boundary_values(arr,value,grid,type)
!!$  use grid_class
!!$  use sml_module, only : sml_zero_inner_bd
!!$  implicit none
!!$  type(grid_type) :: grid
!!$  real (kind=8) :: arr(grid%nnode), value
!!$  integer, intent(in) :: type
!!$  integer :: i
!!$
!!$  if(type == 1 .and. sml_zero_inner_bd==0) then  ! only outer boundary - potential bd
!!$     arr(grid%inv_sort(grid%outer_bd(1):grid%outer_bd(2)))=value
!!$     arr(grid%inv_sort(grid%outer_bd(3):grid%outer_bd(4)))=value
!!$  elseif(type==2) then                           ! both inner (if exist) and outer boundary - pot bd
!!$     if(grid%inner_bd(2)>0) then
!!$        arr(grid%inv_sort(grid%inner_bd(1):grid%inner_bd(2)))=value
!!$     endif
!!$     arr(grid%inv_sort(grid%outer_bd(1):grid%outer_bd(2)))=value
!!$     arr(grid%inv_sort(grid%outer_bd(3):grid%outer_bd(4)))=value
!!$  elseif(type==3 .and. sml_zero_inner_bd==0 ) then     ! only outer bd - charge bd
!!$     arr(grid%inv_sort(grid%outer_bd(5):grid%outer_bd(6)))=value
!!$     arr(grid%inv_sort(grid%outer_bd(7):grid%outer_bd(8)))=value
!!$  else ! type ==4               both inner (if exist) and outer boundary - charge bd
!!$     if(grid%inner_bd(4)>0) then
!!$        arr(grid%inv_sort(grid%inner_bd(3):grid%inner_bd(4)))=value
!!$     endif
!!$     arr(grid%inv_sort(grid%outer_bd(5):grid%outer_bd(6)))=value
!!$     arr(grid%inv_sort(grid%outer_bd(7):grid%outer_bd(8)))=value     
!!$  endif
!!$
!!$end subroutine set_boundary_values



! mode selection - works only the nodes between boundaries
! refine ngrid setting for generalization
subroutine mode_selection(nmode,grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  include 'mpif.h'
  integer, intent(in) :: nmode
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: ngrid,ngrid2,grid_offset,ngrid_per_pe,nphi_total,nplane_per_pe,npe_per_plane 
  integer :: iphi_loc,iphi_offset, iseg,idest,isource,isendtag,irecvtag,istart
  integer :: i,inode,itmp1,itmp2,itmp3,istatus,ierr
  real (kind=8) ,allocatable :: buff1(:,:), buff2(:,:),buff3(:,:)
  complex (kind=8), allocatable :: ftc(:),iftc(:)
  complex (kind=8) :: cn

  ! 1. initialization
  nphi_total=sml_nphi_total
  nplane_per_pe=sml_plane_per_pe
  npe_per_plane=sml_pe_per_plane
  iphi_offset=grid%iphi_offset

!  grid_offset=grid%inner_bd(2)+1  ! hard to exclude boundary node because it is not ordered now
  grid_offset=1
!  ngrid= min(grid%outer_bd(1),grid%outer_bd(3)) - grid_offset ! ngrid is the number of grid point inside boundary condition
  ngrid=grid%nnode
 
  ! Parallelizaing - ngrid_per_pe : number of grid (inside boundary) per cpu 
  ! ngrid2 - number of grid total computed.  Node number which is  > ngrid will be ignored later.
  ngrid_per_pe=ngrid/sml_totalpe 
  if(ngrid==sml_totalpe*ngrid_per_pe) then   
     ngrid2=ngrid  
  else
     ngrid_per_pe=ngrid_per_pe+1
     ngrid2=sml_totalpe*ngrid_per_pe
  endif

  allocate( buff1(ngrid2,nplane_per_pe),buff3(ngrid2,nplane_per_pe),buff2(ngrid_per_pe,nphi_total) )
  allocate( ftc(nphi_total), iftc(nphi_total) )

  ! preparing (inverse) Fourier transform coefficient.
  do i=1, nphi_total
!     ftc(i)=cexp(cmplx(0D0,sml_2pi*real((i-1)*nmode,8)/real(nphi_total,8),8))
     ftc(i)=cexp(cmplx(0D0,sml_2pi*real((i-1)*nmode)/real(nphi_total)))
     iftc(i)=conjg(ftc(i))/real(nphi_total,8)
  enddo

  ! store delta potential value to the buff1 - with grid offset.
  do iphi_loc=1, nplane_per_pe
     buff1(1:ngrid,iphi_loc)=psn%dpot(grid_offset:grid_offset-1+ngrid,iphi_loc)
     buff1(ngrid+1:ngrid2,iphi_loc)=0D0
  enddo


  ! 2. send - receive field data
  ! send buff1 data (given plane field data) to buff2 data (poloidal angle (or grid number) localized data)
  do iseg=1, nphi_total
     iphi_loc=1+(iseg-1)/sml_totalpe
     idest  =mod(sml_mype+npe_per_plane*(iseg-1),sml_totalpe)
     isource=mod(sml_mype-npe_per_plane*(iseg-1),sml_totalpe)  
     isource=mod(isource+sml_totalpe,sml_totalpe) ! to preven minus
     isendtag=sml_mype
     irecvtag=isource
!     istart= 1 + ngrid_per_pe*mod( sml_mype+iseg-1 , sml_totalpe )
     istart= 1 + ngrid_per_pe*idest
     !sending buff1(istart:istart+ngrid_per_pe,iphi)
     !receiving buff2(1:ngrid_per_pe,1+modulo(iphi_offset+iseg-1,nphi_total))
     itmp1=istart
     itmp2=istart+ngrid_per_pe-1
!     itmp3=mod(  iphi_offset + nplane_per_pe*(1-iseg) + iphi_loc-1 , nphi_total   )
!     itmp3= 1 + mod(itmp3 + nphi_total, nphi_total)
     itmp3=iphi_loc + isource*nplane_per_pe/npe_per_plane
     ! send and receive data
     call MPI_SENDRECV( buff1(itmp1:itmp2,iphi_loc), ngrid_per_pe, MPI_DOUBLE_PRECISION,idest,isendtag,&
          buff2(1:ngrid_per_pe,itmp3), ngrid_per_pe, MPI_DOUBLE_PRECISION,isource,irecvtag,&
          MPI_COMM_WORLD,istatus,ierr)
  enddo


  do inode=1, ngrid_per_pe
     
     ! 3. Fourier transform. Single mode direct cal  
     cn=0D0
     do i=1, nphi_total     
        cn= cn + buff2(inode,i)*ftc(i)
     enddo

     ! 4. inverse fouirer transform. Single mode direct cal
     do i=1, nphi_total
        buff2(inode,i)=real(cn*iftc(i))*2D0 ! 2D0 for -nmode summation : nmode and (-nmode) are complex conjugate and should be summed.
     enddo

  enddo

  ! 5. send - receive filtered data
  buff1=0D0
  do iseg=1, nphi_total
     iphi_loc=1+(iseg-1)/sml_totalpe
     idest  =mod(sml_mype-npe_per_plane*(iseg-1),sml_totalpe)
     idest  =mod(idest+sml_totalpe,sml_totalpe) ! to preven minus
     isource=mod(sml_mype+npe_per_plane*(iseg-1),sml_totalpe)  
     isendtag=sml_mype
     irecvtag=isource
!     istart= 1 + ngrid_per_pe*mod( sml_mype+iseg-1 , sml_totalpe )
     istart= 1 + ngrid_per_pe*isource
     !sending buff1(istart:istart+ngrid_per_pe,iphi)
     !receiving buff2(1:ngrid_per_pe,1+modulo(iphi_offset+iseg-1,nphi_total))
     itmp1=istart
     itmp2=istart+ngrid_per_pe-1
!     itmp3=mod(  iphi_offset + nplane_per_pe*(1-iseg) + iphi_loc-1 , nphi_total   )
!     itmp3= 1 + mod(itmp3 + nphi_total, nphi_total)
     itmp3=iphi_loc + idest*nplane_per_pe/npe_per_plane

     ! send and receive data
     call MPI_SENDRECV( buff2(1:ngrid_per_pe,itmp3), ngrid_per_pe,MPI_DOUBLE_PRECISION,idest,isendtag,&
          buff1(itmp1:itmp2,iphi_loc), ngrid_per_pe, MPI_DOUBLE_PRECISION,isource,irecvtag,&
          MPI_COMM_WORLD,istatus,ierr)
  enddo

  ! 6. broad cast between group (for multiple cpu per plane)
  if(npe_per_plane > 1 ) then
     call MPI_ALLREDUCE(buff1,buff3,ngrid2,MPI_DOUBLE_PRECISION, MPI_SUM,sml_plane_comm,ierr)
  else
     buff3=buff1
  endif
  do iphi_loc=1, nplane_per_pe
     psn%dpot(grid_offset:grid_offset-1+ngrid,iphi_loc)=buff3(1:ngrid,iphi_loc)
  enddo
  ! 7. finalize
  deallocate(buff1,buff2,buff3,ftc,iftc)

end subroutine mode_selection


subroutine simple00(grid,psn,den,pot)
  use grid_class
  use psn_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8) :: den(grid%nnode), pot(grid%nnode)
  real (kind=8) :: den00(grid%npsi), phip00(grid%npsi),phi00(grid%npsi), r(grid%npsi)
#ifdef ADIOS
  integer*8::buf
  real (kind=8) :: temppsi(grid%npsi)
#endif
  real (kind=8) :: deltar, get_mid_r, gyroradius, gyro_radius2, x(2)
  integer :: i,nd1,nd2,bd1,bd2,itmp, icount=0, ismooth
  save icount
  icount=icount+1

#ifdef NO_NZERO_MODE
  pot=0D0
  return
#endif

  ! basic parameters
  x(1)=eq_axis_r
  x(2)=eq_axis_z
  gyroradius = gyro_radius2(x)

  bd1=1
  bd2=grid%npsi

  do i=2, grid%npsi-1
!     if( grid%itheta0(i) <= grid%inner_bd(4) .and. grid%inner_bd(4) < grid%itheta0(i+1) ) bd1=i
!     if( grid%itheta0(i) <= grid%outer_bd(5) .and. grid%outer_bd(5) < grid%itheta0(i+1) ) bd2=i
      if( grid%itheta0(i) <= psn%pbd0%in%end .and. psn%pbd0%in%end < grid%itheta0(i+1) ) bd1=i
      if( grid%itheta0(i) <= psn%pbd0%out1%start .and. psn%pbd0%out1%start < grid%itheta0(i+1) ) bd2=i
  enddo

  ! get r(minor-radius)
  do i=1, grid%npsi
     !r= get_mid_r(grid%psi(grid%itheta0(i))) - eq_axis_r   
!     itmp = grid%inv_sort(grid%itheta0(i))
     itmp = grid%itheta0(i)
     r(i)= sqrt( (grid%x(1,itmp)-eq_axis_r)**2  + (grid%x(2,itmp)-eq_axis_z)**2 )
  enddo
  ! 
  deltar=(r(grid%npsi-1)-r(2))/real(grid%npsi-3)  ! for special case only - equal radius


  ! grid -> 1D reduce
  do i=1, grid%npsi
!     den00(i)= den(grid%inv_sort(grid%itheta0(i)))   ! +1 / +0 ???
     den00(i)= den(grid%itheta0(i))   ! +1 / +0 ???
  enddo

  ! (-0.0625 0.25 0.625 0.25 -0.0625) radial smoothing of (0,0) mode density den00
  do ismooth=1,sml_simple00_nsmth
     phip00(bd1)=den00(bd1)
     phip00(bd2)=den00(bd2)
     phip00(bd1+1)=den00(bd1+3)
     phip00(bd2-1)=den00(bd2-3)
     phip00(bd1+2:bd2-2)=den00(bd1:bd2-4)+den00(bd1+4:bd2)
     phip00(bd1+1:bd2-1)=0.625*den00(bd1+1:bd2-1)+0.25*(den00(bd1:bd2-2)+den00(bd1+2:bd2))-&
          0.0625*phip00(bd1+1:bd2-1)
     den00(bd1:bd2)=phip00(bd1:bd2)
  enddo

  phip00=0D0
  do i=bd1+1, bd2
     phip00(i)=phip00(i-1) + 0.5*deltar*(r(i-1)*den00(i-1)+r(i)*den00(i))
  enddo
  
  do i=1, grid%npsi
     phip00(i)=-phip00(i)/(r(i)+1D-50)
  enddo


  ! FLR
  phi00 = den00 * gyroradius**2
  do i=bd1+1, bd2 -1 
     phip00(i) = phip00(i) + 0.5D0 *(phi00(i+1)-phi00(i-1))/deltar
  enddo

  ! 0,0 mode potential
  phi00=0.0
  do i=1+1, grid%npsi
     phi00(i)=phi00(i-1) + 0.5D0 * deltar *(phip00(i-1)+phip00(i))
  enddo

  phi00 = phi00 - phi00(grid%npsi)


  !enforce phi00=0 at radial boundary
!  do i=bd1+1, bd2
!     phi00(i)=phi00(i)-phi00(bd2)*real(i-bd1)/real(bd2-bd1) ! -fomula from original GTC
     ! log correction for more accurate calculation
!     phi00(i)=phi00(i)-phi00(bd2)*log(r(i)/r(bd1))/log(r(bd2)/r(bd1))
!  enddo

  !Dimension factor
!  phi00(:)=phi00(:) * psn%tempi_ev(grid%nnode/2)/psn%density0_full(grid%nnode/2) /gyroradius**2
  phi00(:)=phi00(:) * psn%tempi_ev(grid%nnode/2)/eq_den_edge /gyroradius**2

  ! 1D -> grid 
  do i=1, grid%npsi
     nd1=grid%itheta0(i)
     nd2=grid%itheta0(i)+grid%ntheta(i)-1
!     pot(grid%inv_sort(nd1:nd2))=phi00(i)
     pot(nd1:nd2)=phi00(i)
  enddo

  if(sml_mype==0 .and. mod(icount,20)==0) then
#ifndef ADIOS
     open(unit=987, file='fort.pot_zero', position='append')
#endif
     do i=1, grid%npsi
#ifdef ADIOS
        temppsi(i)=grid%psi(grid%itheta0(i))
#else
        write(987,1000) r(i), phi00(i), den00(i), phip00(i),grid%psi(grid%itheta0(i))
#endif
     enddo
#ifndef ADIOS
     write(987,*) ' '
     close(987)
#else 
    if(sml_mype==0)write(*,*)"i am writing fort.pot_zero to fort.bp"
    call adios_open(buf,"fort.pot_zero"//char(0),"fort.bp"//char(0),"a"//char(0))
    call adios_set_path(buf,"fort.pot_zero"//char(0))
    call adios_gwrite(buf,"fort.pot_zero")
    call adios_close(buf)
#endif
  endif
  
1000 format(200(e19.13,1x))

end subroutine simple00

subroutine neo_pot0_simple(grid,psn)
  use sml_module
  use eq_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8), allocatable :: pot(:), dpot(:), psi(:), p(:),tempi_ev(:)
  integer :: i,  np, nd1, nd2
#if ADIOS
  integer*8 :: buf
  real(kind=8), allocatable :: temp(:)
#endif
  real (kind=8), external :: init_den, init_tempi_ev
  np=grid%npsi

  allocate( pot(np), dpot(np-1), psi(np), p(np), tempi_ev(np-1) )
#if ADIOS
  allocate(temp(np))
#endif
  ! get psi from ordered list
  do i=1, np
     psi(i)=grid%psi(grid%itheta0(i))  ! wrap with reordered 
!     psi(i)=grid%psi( grid%to_ordered_number(grid%itheta0(i)) )
  enddo

  ! get pressure from initial profile
  do i=1, np
     p(i)=init_den(psi(i))*init_tempi_ev(psi(i))     
  enddo

  ! get temperature 
  do i=1, np-1
     tempi_ev(i)=init_tempi_ev(0.5D0*(psi(i)+psi(i+1)))
  enddo

  ! get dpot
  dpot(1:np-1)=-tempi_ev(1:np-1)*(log(p(2:np))-log(p(1:np-1)))/(psi(2:np)-psi(1:np-1))
  
  
  ! get pot from integration
  pot(np)=0D0
  do i=np-1, 1, -1
     pot(i)=pot(i+1)-dpot(i)*(psi(i+1)-psi(i))     
  enddo

  if(sml_mype==0) then
     do i=1, np
#ifdef ADIOS
        temp(i)=psi(i)/eq_x_psi
#else
        write(700,*) psi(i)/eq_x_psi, pot(i)     
#endif
     enddo
#ifdef ADIOS
     call adios_open(buf,"fort.700"//char(0),"fort.bp"//char(0),"a"//char(0))
     call adios_set_path(buf,"fort.700"//char(0))
     call adios_gwrite(buf,"fort.700")
     call adios_close(buf)
#else
     close(700)
#endif
  endif
  ! assign to grid

  do i=1, np
     nd1=grid%itheta0(i)
     nd2=grid%itheta0(i)+grid%ntheta(i)-1
!     psn%add_pot0(grid%inv_sort(nd1:nd2))=pot(i)
     psn%add_pot0(nd1:nd2)=pot(i)
  enddo
  deallocate( pot, dpot, psi, p, tempi_ev)
#ifdef ADIOS 
  deallocate(temp)
#endif
  

end subroutine

subroutine get_add_pot0_psi(psi_in,pot,grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8), intent(in) :: psi_in
  real (kind=8), intent(out) :: pot
  real (kind=8) :: psi,x(2), p(3)
  integer :: itr, nodes(3), j
  logical, parameter :: USE_SEARCH_TR2 = .true.
  real (kind=8) , external :: get_mid_r
  
  ! exit with 0
  if(sml_add_pot0<=0) then
     pot=0D0
     return
  endif
  ! set psi
  psi=min(max(psi_in,sml_inpsi),sml_outpsi)

  ! get r,z=0 from psi
  x(1)=get_mid_R(psi)
  x(2)=eq_axis_z

  ! find triangle
        
  ! find position and save it.
  if (USE_SEARCH_TR2) then
     call search_tr2(grid,x,itr,p)
  else
     call search_tr(grid,x,itr,p)
  endif
  

  ! get potential value
  pot=0D0
  if( itr > 0 ) then
     
     nodes=grid%nd(:,itr)
     
     do j=1, 3
        pot= pot + p(j) * psn%add_pot0(nodes(j))
     enddo
  else
     print *, 'Warning : itr <0 in get_add_pot0_psi', psi, psi_in
     call err_count
  endif


end subroutine get_add_pot0_psi


subroutine extract_00mode(grid,phi)
  use grid_class
  use sml_module
  use smooth_module
  implicit none
  type(grid_type) :: grid
  real (kind=8) :: phi(grid%nnode,-1:grid%nphi+1)
  real (kind=8) :: tmp(grid%nnode)
  integer :: i

  do i=1, grid%nnode
     tmp(i)=sum(phi(i,1:grid%nphi))
  enddo
  ! sum-up
  call my_mpi_allreduce(tmp,grid%rtmp1,grid%nnode) !rtmp1 is used for a temporory purpose. variable
  ! get density from charge summation
  call smooth_pol0(grid, grid%rtmp1, smooth00)

  do i=1, grid%nphi
     phi(:,i)=phi(:,i) - grid%rtmp1(:)/(real(grid%nphi*sml_totalpe))
  enddo
  

end subroutine extract_00mode

subroutine apply_wall_boundary_condition(grid,psn,den)
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8) :: den(grid%nnode)
  integer :: i

  do i=1, psn%nwall
     den(psn%wall_start+i-1)=psn%sheath_pot(i)
  enddo
  
end subroutine apply_wall_boundary_condition

subroutine  sheath_pot_init(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i
  real (kind=8) :: x(2), Te, psi
  real (kind=8) :: tempe_ev_wz, psi_interpol
! apply sheat potential according to temperature


   do i=1, psn%nwall
      x=grid%x(:,psn%wall_start+i-1)
      psi=psi_interpol(x(1),x(2),0,0)
      ! get initial temperature
      Te=tempe_ev_wz(psi,x(2))
      ! apply alpha *  T
      psn%sheath_pot(i)=sml_sheath_init_pot_factor*Te
   enddo 

end subroutine
!******************************************************************************************************
!********  Funtions for debugging ********************************************************************
!******************************************************************************************************
!subroutine enforce_modeled_charge(gird,psn)



!end subroutine enforce_modeled_charge



!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!***********************************************************************************************************
!---------------------------- Olds functions ----------------------------------------------------------------

