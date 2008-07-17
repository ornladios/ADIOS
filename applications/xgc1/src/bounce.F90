!first created 2002/07/03
!last modified 2004/01/20
! 2002/08/16 minimum search -> motion following
! let the ptl bounce in at the psi boundary
! 2003/02/04 adopted from xorbit , conentric part added
! 2004/01/20 motion following -> binary search. 1st order correction.
! 2008/01/26 implememnt 2D newton method for more accurate bounce operation   

subroutine bounce(new_phase,old_phase,rtn,nphase)
  use eq_module , only : eq_axis_r,eq_x_psi
  use sml_module
  use bnc_module
  implicit none
  integer, intent(in) :: nphase
  real (kind=8),intent(inout) :: new_phase(nphase),old_phase(nphase)
  integer, intent(out) :: rtn
  integer :: inner
  real (kind=8) :: r, z, psi,sign
  real (kind=8) :: b , br, bz, bphi,mid,z1,z2,z_pmin
  integer :: i,j,count
  integer, parameter :: JMAX=5
  real (kind=8), external :: psi_interpol,z_psi_min,I_interpol,B_interpol
  !  real (kind=8), parameter :: ZTOL=1D-6, NPSI_TOL=1D-5
  !  real (kind=8), parameter :: ZTOL=1D-10, NPSI_TOL=1D-5
  real (kind=8), parameter :: ZTOL=1D-6, NPSI_TOL=1D-4, BTOL=1D-8, PTOL=1D-8
  real (kind=8) :: psitol
  real (kind=8) :: delz,deltab2,deltap2
  real (kind=8) :: RBt,oldB,oldP,newB,newP,deltab,deltap,sign2,deltar,deltaz,zh,rh,dz1,dz2,dz3
  real (kind=8) :: dpsi_dr,dpsi_dz,d2psi_d2r,d2psi_drdz,d2psi_d2z,db_dr,db_dz,denom
  logical :: use_current, diverge
  psitol=NPSI_TOL*eq_x_psi
  rtn=0

  ! reverse angle --   z -> -z, phase variable -> old phase variable except z
  if(sml_concentric==1) then
     psi=new_phase(9)
!     if(psi<=sml_inpsi .or. sml_bounce==2 .and. psi>= sml_outpsi) then
     if( (sml_bounce==1 .or. sml_bounce==2) .and. psi < sml_inpsi .or. (sml_bounce==2 .or. sml_bounce==3) .and. psi > sml_outpsi ) then
        new_phase(1:nphase)=old_phase(1:nphase)
        new_phase(2)=-old_phase(2)
        new_phase(9)=psi_interpol(new_phase(1),new_phase(2),0,0)
        if(sml_bounce_zero_weight==1) new_phase(6:7)=0D0
     endif
  else 
     ! bounce at inner (& outer )boundary
     ! Finding z value which has same psi and r of old particle position
     ! This position gives same mu, P_phi.
     ! similar energy if B is similar.
     ! for exact calculation, 1st order correction is needed.
     psi=new_phase(9)  ! psi value of new phase
     r=old_phase(1)    ! radius of old position
     !     if((psi <= sml_inpsi .or. sml_bounce==2 .and. psi >= sml_outpsi) .and. &
     !        r < bnc_max_r .and. r > bnc_min_r .and. psi < eq_x_psi ) then                
     if(( (sml_bounce==1 .or. sml_bounce==2) .and. psi < sml_inpsi .or. (sml_bounce==2 .or. sml_bounce==3) .and. psi > sml_outpsi) .and. &
          r < bnc_max_r .and. r > bnc_min_r .and. psi < eq_x_psi ) then                
        rtn=1 ! a big jump will be happen
        if(psi<=sml_inpsi) then
           inner=1
        else
           inner=0
        endif
        psi=old_phase(9)  !  psi is old position value
        !        if(psi > (sml_inpsi-psitol) .and. psi < (sml_outpsi+psitol)) then  !old psi is inside simulation range
        if(psi >= (sml_inpsi-psitol) .and. psi <=(sml_outpsi+psitol)) then  !old psi is inside simulation range
           z=old_phase(2)
           z_pmin=z_psi_min(r)

           ! bug fixed. sign should be -1 
           !                 if(r>eq_axis_r) then
           !                    sign=1
           !                 else
           sign=-1
           !                 endif

           ! set boundary of binary search
           z1=z_pmin              
           if(z>z_pmin) then                                   
              !             z2=max(z_pmin - (z-z_pmin+0.01)*2. , sml_bd_min_z) ! heuristic formula
              if( z_pmin- (z-z_pmin+0.01)*2. > sml_bd_min_z ) then
                 z2= z_pmin- (z-z_pmin+0.01)*2.
              else
                 z2=sml_bd_min_z
              endif
           else
              !                z2=min(z_pmin - (z-z_pmin+0.01/sml_norm_r)*2.  , sml_bd_max_z) ! heuristic formula
              if(z_pmin - (z-z_pmin+0.01)*2. < sml_bd_max_z) then
                 z2=z_pmin - (z-z_pmin+0.01)*2
              else
                 z2=sml_bd_max_z
              endif
           endif

           ! find z value using binary search.-------------------------
           do while (abs(z1-z2) > ZTOL)
              mid=0.5D0*(z1+z2)
              if(sign *(psi_interpol(r,mid,0,0)-psi) < 0 ) then
                 z2=mid
              else
                 z1=mid
              endif
              !                    write(400,*) mid,z1,z2,z1-z2,psi,psi_interpol(r,mid,0,0)
           enddo

           !           z1=0.5D0*(z1+z2)
           if(inner==1) then  ! z1 gives larger psi for inner bd.
              z1=z2
           endif
           !------------------------------------------------------------

           new_phase(1:nphase)=old_phase(1:nphase)
           new_phase(2)=z1              

           ! 1st order correction  
           !           delz=(psi-new_phase(9))/psi_interpol(new_phase(1),z1,0,1)
           !           new_phase(2)=new_phase(2)+delz
           !           new_phase(9)=psi_interpol(new_phase(1),new_phase(2),0,0)

           ! (dR,dZ) correction for attaining position with same B and psi values as old position  
           oldB=B_interpol(old_phase(1),old_phase(2))
           oldP=old_phase(9)
           deltab2=oldB-B_interpol(new_phase(1),new_phase(2))
           deltap2=oldP-psi_interpol(new_phase(1),z1,0,0)

           !init
           use_current=.false.
           diverge=.false.
           do j=1, JMAX
              ! Newton-Raphson procedure is iterated 
              psi=psi_interpol(new_phase(1),new_phase(2),0,0)
              dpsi_dr=psi_interpol(new_phase(1),new_phase(2),1,0)
              dpsi_dz=psi_interpol(new_phase(1),new_phase(2),0,1)
              RBt=I_interpol(psi,0,1)
              newB=sqrt(RBt**2+dpsi_dr**2+dpsi_dz**2)/new_phase(1)
              deltab=oldB-newB
              deltap=oldP-psi
              if(((dabs(deltab)/oldB)<BTOL).and.((dabs(deltap)/oldP)<PTOL)) then 
                 use_current=.true.
                 exit
              endif

              d2psi_d2r=psi_interpol(new_phase(1),new_phase(2),2,0)
              d2psi_drdz=psi_interpol(new_phase(1),new_phase(2),1,1)
              d2psi_d2z=psi_interpol(new_phase(1),new_phase(2),0,2)

              db_dr=-newB/new_phase(1)+1D0/(newB*(new_phase(1)**2))*(dpsi_dr*d2psi_d2r+dpsi_dz*d2psi_drdz)
              db_dz=1D0/(newB*(new_phase(1)**2))*(dpsi_dr*d2psi_drdz+dpsi_dz*d2psi_d2z)
              
              denom=dpsi_dr*db_dz-dpsi_dz*db_dr
              deltar=(db_dz*deltap-dpsi_dz*deltab)/denom
              deltaz=(-db_dr*deltap+dpsi_dr*deltab)/denom
              
              ! move to new (r,z) 
              new_phase(1)=new_phase(1)+deltar
              new_phase(2)=new_phase(2)+deltaz
              
              ! check diverge
              if(new_phase(1) < sml_bd_min_r .or. new_phase(1) > sml_bd_max_r &
                   .or. z < sml_bd_min_z .or. z > sml_bd_max_z ) then
                 diverge=.true.
                 use_current=.false. ! for safety
                 exit
              endif                 
           enddo

           if(.NOT. use_current .and. .NOT. diverge ) then
              ! loop ended not satisfying the TOL codition, nor diverge
              psi=psi_interpol(new_phase(1),new_phase(2),0,0)
              dpsi_dr=psi_interpol(new_phase(1),new_phase(2),1,0)
              dpsi_dz=psi_interpol(new_phase(1),new_phase(2),0,1)
              RBt=I_interpol(psi,0,1)
              newB=sqrt(RBt**2+dpsi_dr**2+dpsi_dz**2)/new_phase(1)
              deltab=oldB-newB
              deltap=oldP-psi
              ! use original (binary search position)
              if(((deltab2/oldB)**2+(deltap2/oldP)**2)<((deltab/oldB)**2+(deltap/oldP)**2)) then
                 use_current=.false.
              endif
           endif
           
           if(.NOT. use_current) then  
              new_phase(1:nphase)=old_phase(1:nphase)
              new_phase(2)=z1
           endif
           
           new_phase(9)=psi_interpol(new_phase(1),new_phase(2),0,0)
           ! end of second (dR,dZ) correction for same B

           if(new_phase(9) < (sml_inpsi-psitol) .or. new_phase(9) > (sml_outpsi+psitol)) then
              !           if(new_phase(9) <= sml_inpsi-psitol .or. new_phase(9) >= sml_outpsi+psitol) then
              print *, 'Fail finding proper psi', new_phase(9)/eq_x_psi, old_phase(9)/eq_x_psi, psi/eq_x_psi
              print *, 'oldB, newB', oldB, newB
              rtn=-1
           endif
           if(sml_bounce_zero_weight==1) new_phase(6:7)=0D0           
        else ! old position was outside of boundary
           rtn=-1
           !           print *, 'ptl_elimination', i, sml_mype
        endif
     endif
  endif

end subroutine bounce

real (kind=8) function z_psi_min(r)
  use sml_module
  use bnc_module
  implicit none
  real (kind=8),intent(in) :: r
  integer :: i
  real (kind=8) :: aa,bb
  
  i=int((r-bnc_min_r)/bnc_dr) +1
  i=min(bnc_nr-1,max(1,i))
  bb=(r-bnc_min_r)/bnc_dr + 1D0 -i
  aa=1D0-bb
!  print *, 'zaa', r, i, aa,bb,bnc_z_psi_min(i)
  z_psi_min=bnc_z_psi_min(i)*aa + bnc_z_psi_min(i+1)*bb

end function z_psi_min

!finding bnc_z_psi_min(r_index)
! z_psi_min is z value of the position in whicth psi is minimum at given R
! This line is nearly midplane and z_psi_min is nearly 0.
  
subroutine bounce_setup
  use sml_module
  use bnc_module
  implicit none
  integer :: i
  real (kind=8) ::  xguess,errrel,gtol,a,b,x,fx,gx
  integer :: maxfn
  real (kind=8), external :: bnc_psi_imsl,bnc_dpsi_imsl
#if !defined(IMSL)
  real (kind=8), external :: fmin
#endif
#ifdef ADIOS
#define ADIOS_WRITE(a,b) call adios_write (a,'b'//char(0),b)
  integer*8 :: grp_id,buf_id
  real (kind=8), allocatable :: bnc_r_dr(:)
  allocate(bnc_r_dr(bnc_nr))
#endif
  xguess=0D0
  errrel=1D-7
  gtol=1D-7
  a=sml_bd_min_z * 0.3  ! heuristic formula
  b=sml_bd_max_z * 0.3
  
  do i=1, bnc_nr
     bnc_arg_pass_r=bnc_min_r + bnc_dr*real(i-1)
#if defined(IMSL)
     call duvmid(bnc_psi_imsl,bnc_dpsi_imsl,xguess,errrel,gtol,maxfn,a,b,x,fx,gx)
     bnc_z_psi_min(i)=x
#else
     bnc_z_psi_min(i)=fmin(a,b,bnc_psi_imsl,gtol)
#endif
  enddo

  if(sml_mype==0) then
     do i=1, bnc_nr
#ifdef ADIOS   
        bnc_r_dr(i)=bnc_min_r+bnc_dr*real(i-1)
#else   
        write(103,*) (bnc_min_r+bnc_dr*real(i-1)), bnc_z_psi_min(i)
#endif
     enddo
#ifdef ADIOS
     call adios_open(buf_id,grp_id,"fort.103"//char(0),"fort.bp"//char(0),"a"//char(0))
     call adios_set_path(grp_id,"fort.103"//char(0))
     call adios_gwrite(buf_id,"fort.103")
     call adios_close(buf_id)
     deallocate(bnc_r_dr)
#else 
     close(100)
#endif
  endif

end subroutine bounce_setup

real (kind=8) function bnc_psi_imsl(z)
  use bnc_module
  implicit none
  real (kind=8) :: psi_interpol,z

  bnc_psi_imsl=psi_interpol(bnc_arg_pass_r,z,0,0)
end function bnc_psi_imsl

real (kind=8) function bnc_dpsi_imsl(z)
  use bnc_module
  implicit none
  real (kind=8) :: psi_interpol,z

  bnc_dpsi_imsl=psi_interpol(bnc_arg_pass_r,z,0,1)
end function bnc_dpsi_imsl
