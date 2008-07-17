!!particle (ion) pushing routine using R-K 2nd order method
subroutine push(istep,ipc,grid,psn,sp)
  use sml_module
  use ptl_module
  use fld_module
  use grid_class
  use psn_class
  use perf_monitor
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  type(species_type) :: sp
  integer, intent(in) :: ipc, istep !! RK2 index
  real (kind=8) :: y(sp%nphase),yprime(sp%nphase),dt_now,time_now,new_phase(sp%nphase),old_phase(sp%nphase)
  integer :: i,rtn,j
  real (kind=8) , external :: psi_interpol
  character (len=5) :: err_str(2)
  logical, parameter :: USE_SEARCH_TR2 = .true.
  interface
     subroutine reset_derivs_without_e(new_phase,sp,i,time_now)
       use sml_module
       use ptl_module
       type(species_type) :: sp
       integer :: i
       real (kind=8) :: new_phase(sp%nphase),time_now
     end subroutine reset_derivs_without_e
  end interface
  err_str(1)='ion'
  err_str(2)='elec'

  if(sml_push_mode==1 .or. sp%stored_derivs<sml_pc_order) then  !RK2
     if(ipc==1) then
        dt_now=0.5D0*sml_dt
        time_now=sml_time
        sp%phase0(:,1:sp%num)=sp%phase(:,1:sp%num)  !store initial phase variable
     else
        dt_now=sml_dt
        time_now=sml_time+0.5D0*sml_dt
     endif
  else  ! P-C
     dt_now=sml_dt  ! dt is the same for PC
     if(ipc==1) then 
        time_now=sml_time
        sp%phase0(:,1:sp%num)=sp%phase(:,1:sp%num)  !store initial phase variable
     else
        time_now=sml_time+sml_dt
     endif
  endif


  call monitor_start (PUSH_ALL_PART_)
  do i=1, sp%num
     y(:)=sp%phase(:,i)
     if(sp%gid(i)>0) then
        if(sml_push_mode==1 .or. sp%stored_derivs < sml_pc_order) then  ! RK2
           new_phase(:)=sp%phase0(:,i) + dt_now*sp%derivs(:,sp%pc_index(1),i)
        else  !PC
           if(sml_pc_order==3) then
              new_phase(:)=sp%phase0(:,i) + dt_now* &
                   ( sml_pc_coeff(1,ipc)* sp%derivs(:,sp%pc_index(1),i) &
                   + sml_pc_coeff(2,ipc)* sp%derivs(:,sp%pc_index(2),i) &
                   + sml_pc_coeff(3,ipc)* sp%derivs(:,sp%pc_index(3),i)  )
           elseif(sml_pc_order==4) then
              new_phase(:)=sp%phase0(:,i) + dt_now* &
                   ( sml_pc_coeff(1,ipc)* sp%derivs(:,sp%pc_index(1),i) &
                   + sml_pc_coeff(2,ipc)* sp%derivs(:,sp%pc_index(2),i) &
                   + sml_pc_coeff(3,ipc)* sp%derivs(:,sp%pc_index(3),i) &
                   + sml_pc_coeff(4,ipc)* sp%derivs(:,sp%pc_index(4),i) )
           elseif(sml_pc_order==6) then
              new_phase=0D0
              do j=1, sml_pc_order
                 new_phase(:)=new_phase(:) + sml_pc_coeff(j,ipc) * sp%derivs(:,sp%pc_index(j),i)
              enddo
              new_phase(:)=sp%phase0(:,i) + dt_now*new_phase(:)
           else
              print *, 'Error : Unprepared predictor-corrector order', sml_pc_order
              !! generalized expression will be implimented later
              stop
           endif
        endif

        new_phase(9)=psi_interpol(new_phase(1),new_phase(2),0,0)



        ! bounce 
        if(ipc==2 .and. sml_bounce/=0) then
           old_phase(:)=sp%phase0(:,i)
           call bounce(new_phase,old_phase,rtn,sp%nphase)
           if(rtn<0)  then
              call remove_particle(sp,i,-2)
!              print *, 'particle eliminated due to bounce error', i, sml_mype, err_str(sp%type),sp%gid(i)
           elseif (rtn==1 .and. sml_push_mode/=1) then
              call reset_derivs_without_e(new_phase,sp,i,time_now)
!              sp%reset_derivs(i)=1
           endif
        endif

        ! time advance one step
        sp%phase(:,i)= new_phase(:)
!     else
!        call remove_particle(sp,i,-1)
!        ptl_rz_outside=0
!        print *, 'particle eliminated due to rz_outside', i, sml_mype, err_str(sp%type), sp%gid(i)
     endif
  enddo
  call monitor_stop (PUSH_ALL_PART_)

  ! update pc_index. 
  ! This operation should be done at the end of push(ipc=1), 
  ! predictor step requires another space and does not use oldest derivs data 
  if(ipc==1) then
     call shift_pc_index(sp%pc_index,sml_pc_order)
  else
     ! update stored_derivs
     ! RK2 initial steps for PC
     sp%stored_derivs=min(sp%stored_derivs+1,sml_pc_order) !update stored number after pushing step is over
  endif
end subroutine push

subroutine phase_derivs(istep,ipc,grid,psn,sp)
  use sml_module
  use ptl_module
  use fld_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  type(species_type) :: sp
  integer, intent(in) :: ipc, istep !! RK2 index
  real (kind=8) :: y(sp%nphase),yprime(sp%nphase),dt_now,time_now,new_phase(sp%nphase),old_phase(sp%nphase)
  integer :: i,rtn,j
  real (kind=8) , external :: psi_interpol
  character (len=5) :: err_str(2)
  real (kind=8) :: dum

  err_str(1)='ion'
  err_str(2)='elec'

  if(sml_push_mode==1 .or. sp%stored_derivs<sml_pc_order-1) then  !RK2
     if(ipc==1) then
        time_now=sml_time
     else
        time_now=sml_time+0.5D0*sml_dt
     endif
  else  ! P-C
     dt_now=sml_dt  ! dt is the same for PC
     if(ipc==1) then 
        time_now=sml_time
     else
        time_now=sml_time+sml_dt
     endif
  endif


  do i=1, sp%num
     y(:)=sp%phase(:,i)
     if(sp%gid(i)>0) then
        ! Save space information
        fld%r=y(1)
        fld%z=y(2)
        fld%phi=y(3)
        ptl_rz_outside=0
        ! obtain B-field information 
        call field(fld,time_now)
        if(ptl_rz_outside==0) then
           ! obtain gyro-averaged E-field for each particle : use position information from 'charge'
           call efield(grid,psn,sp,i,fld,time_now)
           if(sml_deltaf==1) call f0_info(grid,i,sp,fld)
           ! get time derivatives
           call derivs_sp(fld,time_now,y,yprime,sp%type)
                      
           ! Fill empty derivative space
!           if(sml_push_mode==1 .or. sp%stored_derivs < sml_pc_order) then ! RK2 mode
              !Do nohting
!           else 
!              if(sp%reset_derivs(i)==1 .and. ipc==1) then
!                 do j=1, sml_pc_order
!                    sp%derivs(:,j,i)=yprime(:)
!                 enddo
!                 sp%reset_derivs(i)=0
!              endif
!           endif
           ! above routine is not needed - derivative values are filled always except rk2 mode

           
           !*****************************************************
           !               Save derivatives
           !******************************************************
           sp%derivs(:,sp%pc_index(1),i)=yprime(:) 

#ifdef COREVV2
           ! for heat flux diagnosis of core turbulence verification
           ! get spatial coordinate derivative only from ExB drift
                   
           
           !over_B2=  1D0/( fld%br**2 + fld%bz**2 + fld%bphi **2 )
           dum=   1D0/( fld%br**2 + fld%bz**2 + fld%bphi **2 )
           ptl_exb(1,i)=  (fld%bz * fld%Ephi - fld%Bphi * fld%Ez) * dum
           ptl_exb(2,i)=  (fld%bphi * fld%Er - fld%br * fld%Ephi) * dum


#endif           
           
        else
           call remove_particle(sp,i,-1)
           ptl_rz_outside=0
           print *, 'particle eliminated due to rz_outside', i, sml_mype, err_str(sp%type), sp%gid(i)
        endif
     endif
  enddo


end subroutine phase_derivs

subroutine reset_derivs_without_e(new_phase,sp,i,time_now)
  use sml_module
  use ptl_module
  implicit none
  type(species_type):: sp
  integer :: i
  real (kind=8) :: new_phase(sp%nphase),time_now
  real (kind=8) :: y(sp%nphase), yprime(sp%nphase),t,dt
  integer :: k
  interface
     subroutine single_phase_derivs_without_E(y,yprime,sp,i,time)
       use ptl_module
       use fld_module
       type(species_type) :: sp
       real (kind=8) :: y(sp%nphase), yprime(sp%nphase),time
     end subroutine single_phase_derivs_without_E
     subroutine rk4_single_without_e(y,yprime,sp,i,t,dt)
       use ptl_module
       use fld_module
       type(species_type) :: sp
       real (kind=8) :: y(sp%nphase), yprime(sp%nphase),t,dt
       integer ,intent(in):: i
     end subroutine rk4_single_without_e
  end interface
  y=new_phase
  t=time_now
  dt=-sml_dt
  ! get derivative
  call single_phase_derivs_without_E(y,yprime,sp,i,t)  
  !store derivative - it doesn't need because ipc=2 derivative is meaningless. just for safety.
  sp%derivs(:,sp%pc_index(1),i)=yprime

  do k=2,sml_pc_order
     if(sp%gid(i)>0) then
        ! advance backward
        call rk4_single_without_e(y,yprime,sp,i,t,dt)
        ! get derivative
        call single_phase_derivs_without_E(y,yprime,sp,i,t) 
        sp%derivs(:,sp%pc_index(k),i)=yprime        
     else
        return
     endif
  enddo
  
end subroutine reset_derivs_without_e
subroutine single_phase_derivs_without_E(y,yprime,sp,i,time)
  use ptl_module
  use fld_module
  use sml_module
  implicit none
  type(species_type) :: sp
  real (kind=8) :: y(sp%nphase), yprime(sp%nphase),time
  integer :: i
  type(fld_type) :: fld

  fld%r=y(1)
  fld%z=y(2)
  fld%phi=y(3)
  ptl_rz_outside=0
  call field(fld,time)
  if(ptl_rz_outside==0) then
     fld%Er=0D0
     fld%Ez=0D0
     fld%Ephi=0D0
     if(sml_deltaf==1) then 
        !call f0_info(grid,i,sp,fld) - cannot use this
        fld%f0_den=1D0  ! arbitary non zero value
        fld%f0_temp=1D0 ! arbitraty non zero value
        fld%f0_gradn(:)=0D0    ! zero gradient is assumed
        fld%f0_gradT(:)=0D0    ! zero gradient
     endif
     call derivs_sp(fld,time,y,yprime,sp%type)
  else
     call remove_particle(sp,i,-1)
     return
  endif

end subroutine single_phase_derivs_without_E

subroutine rk4_single_without_e(y,yprime,sp,i,t,dt)
  use ptl_module
  use fld_module
  implicit none
  type(species_type) :: sp
  real (kind=8) :: y(sp%nphase), yprime(sp%nphase),t,dt
  integer ,intent(in):: i
  real (kind=8) :: h6, hh, th,dym(sp%nphase), dyt(sp%nphase),yt(sp%nphase)
  interface
     subroutine single_phase_derivs_without_E(y,yprime,sp,i,time)
       use ptl_module
       use fld_module
       type(species_type) :: sp
       real (kind=8) :: y(sp%nphase), yprime(sp%nphase),time
     end subroutine single_phase_derivs_without_E
  end interface
  hh=dt*0.5D0
  h6=dt/6D0
  th=t+HH

  !write(1004,*) sml_istep, sml_time, sml_dt, x,h,xh
  !added to get dydx in this routine
!  call derivs(sp,x,y,dydx)
  !diagnosis -- save it for fast diagnosis
!     diag_derivs(i,:)=dydx

  yt=y + hh* yprime

!  CALL DERIVS(fld,XH,YT,DYT,sp_type,grid)
  call single_phase_derivs_without_E(yt,dyt,sp,i,th)
  yt=y+hh*dyt

!  CALL DERIVS(fld,XH,YT,DYM,sp_type,grid)
  call single_phase_derivs_without_E(yt,dym,sp,i,th)
  yt=y+dt*dym
  dym=dyt + dym

!  CALL DERIVS(fld,X+H,YT,DYT,sp_type,grid)
  call single_phase_derivs_without_E(yt,dyt,sp,i,t+dt)
  y=y + h6 * ( yprime+dyt+2D0*dym )

end subroutine rk4_single_without_e





! Code from gtc
subroutine shift_ie(grid,ptl,init)
  use sml_module
  use ptl_module
  use grid_class
  use diag_module  
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(ptl_type) :: ptl
  integer, intent(in) :: init
  type(species_type), pointer :: sp
  integer i,m,msendleft(2),msendright(2),mrecvleft(2),mrecvright(2),mtop,&
       m0,msend,mrecv,idest,isource,isendtag,irecvtag,&
       isendcount,irecvcount,istatus(MPI_STATUS_SIZE),ierror,iteration
  integer,dimension(:),allocatable :: ihole
  real (kind=8),dimension(:,:),allocatable :: recvleft,recvright,sendleft,sendright
#ifdef SMALL_GID
  integer,dimension(:),allocatable :: recvleftid,recvrightid,sendleftid,sendrightid
#else
  integer (kind=8),dimension(:),allocatable :: recvleftid,recvrightid,sendleftid,sendrightid
#endif
  real (kind=8) :: zetaright,zetaleft,zetamid,pi_inv
  character(len=8) cdum
  integer :: sp_type
  integer :: st,last
  integer :: ptl_shiftmax
  integer, pointer :: ptl_rshift(:), ptl_lshift(:),ptl_hole(:)
  integer :: sendsize,arr_start, arr_end, pc_index
  save ptl_rshift,ptl_lshift,ptl_hole
  target ptl

  if(init==0) then
     ptl_shiftmax=ptl%ion%maxnum / 3
     allocate(ptl_rshift(ptl_shiftmax),ptl_lshift(ptl_shiftmax),ptl_hole(2*ptl_shiftmax))
     return
  endif

  do sp_type=1, 1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif

     ! Modulo operation
     do i=1, sp%num
        if(sp%phase(3,i)>= sml_2pi .or. sp%phase(3,i)< 0D0 ) then
           sp%phase(3,i)=modulo(sp%phase(3,i),sml_2pi)
        endif
     enddo

     ! Shift operation
     pi_inv=1.D0/sml_pi
     m0=1
     iteration=0
     mrecv=1
     
     INFLOOP: do while (1==1)

        iteration=iteration+1
        if(iteration>sml_totalpe)then
           write(sml_stdout,*)'endless particle sorting loop at PE=',sml_mype
           stop
        endif
        
        msend=0
        msendleft=0
        msendright=0
!        print *, 'iter=',iteration,sml_mype

        !************* Phase 1 *******************************
        ! Finding particle to be shifted
        !*****************************************************
        if(m0 <= sp%num)then
           
           do m=m0,sp%num
              zetaright=min(sml_2pi,sp%phase(3,m))-grid%phimax
              zetaleft=sp%phase(3,m)-grid%phimin

              if( zetaright*zetaleft > 0 )then
!                 print *, '[',sml_mype,']',msend,phase(m,3)*pi_inv,grid%phimax*pi_inv, grid%phimin*pi_inv
                 msend=msend+1
                 ptl_hole(msend)=m
                 zetaright=zetaright*0.5D0*pi_inv
                 zetaright=zetaright-real(floor(zetaright))

                 if( zetaright < 0.5 )then
                    ! # of particle to move right
                    msendright(1)=msendright(1)+1
                    ptl_rshift(msendright(1))=m
                    ! keep track of tracer
                    if( m == diag_tracer_n )then
                       msendright(2)=msendright(1)
                       diag_tracer_n=0
                    endif
                    ! # of particle to move left
                 else
                    msendleft(1)=msendleft(1)+1
                    ptl_lshift(msendleft(1))=m
                    if( m == diag_tracer_n )then
                       msendleft(2)=msendleft(1)
                       diag_tracer_n=0
                    endif
                 endif
              endif
           enddo
        endif

        ! total # of particles to be shifted for whole CPUs
        call monitor_start (SHIFT_IE_RED_)
        mrecv=0
        msend=msendright(1)+msendleft(1)
        call MPI_ALLREDUCE(msend,mrecv,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
        call monitor_stop (SHIFT_IE_RED_)

!        print *, '[',sml_mype,']','msend,right,left,recv',msend,msendright(1),msendleft(1),mrecv
        ! no particle to be shifted, free memory and return

        if ( mrecv == 0 ) then
           if(msend/=0)  then
              print *,'Error in shfit', msendright(1), msendleft(1), sml_mype
              stop
           endif
!           print *, sml_mype,'interation # :', iteration
           !     if(sml_mype==5) print *, m0,ptl_num
!           print *, pnum,sml_mype
          
           call shift_check(grid,sp) ! debug only
        
           exit INFLOOP
        endif
           
        ! an extra space to prevent zero size when msendright(1)=msendleft(1)=0
!        print *, 'mem allocation',sml_mype !debug
        !allocate memory for particle sending
        sendsize=sp%nphase2
        if(sml_push_mode/=1) then !P-C
           sendsize=sendsize + sp%nphase * sml_pc_order
        endif
        allocate(sendright(sendsize,max(msendright(1),1)),sendleft(sendsize,max(msendleft(1),1)))
        allocate(sendrightid(max(msendright(1),1)),sendleftid(max(msendleft(1),1)))
!        sendright=0D0 !debug
!        sendleft=0D0  !debug
!        sendrightid=0 !debug
!        sendleftid=0  !debug

        ! pack particle to move right
        do m=1, msendright(1)
           sendright(1:sp%nphase,m)=sp%phase(1:sp%nphase,ptl_rshift(m))
           sendright(sp%nphase+1:sp%nphase2,m)=sp%phase0(1:sp%nphase,ptl_rshift(m)) ! "
           if(sml_push_mode/=1) then
              do pc_index=1, sml_pc_order
                 arr_end=sp%nphase2 + sp%nphase*pc_index
                 arr_start=arr_end -sp%nphase +1
                 sendright(arr_start:arr_end,m)=sp%derivs(:,pc_index,ptl_rshift(m))
              enddo
           endif
           sendrightid(m)=sp%gid(ptl_rshift(m))
        enddo
        do m=1, msendleft(1)
           sendleft(1:sp%nphase,m)=sp%phase(1:sp%nphase,ptl_lshift(m))
           sendleft(sp%nphase+1:sp%nphase2,m)=sp%phase0(1:sp%nphase,ptl_lshift(m))
           if(sml_push_mode/=1) then
              do pc_index=1, sml_pc_order
                 arr_end=sp%nphase2 + sp%nphase*pc_index
                 arr_start=arr_end -sp%nphase +1
                 sendleft(arr_start:arr_end,m)=sp%derivs(:,pc_index,ptl_lshift(m))
              enddo
           endif
           sendleftid(m)=sp%gid(ptl_lshift(m))
        enddo
!        print *, '[',sml_mype,']','sendright',pnum,mtop,last
                
        mtop=sp%num
        ! # of particles remain on local PE
        sp%num=sp%num-msendleft(1)-msendright(1)
        ! fill the hole
        last=msend
!        print *, '[',sml_mype,']','pnum,mtop,last',pnum,mtop,last
        FILLHOLE : do i=1, msend
           m=ptl_hole(i)
           if( m > sp%num )  exit FILLHOLE
           !when empty space in the end - possible only for i=1
           do while( mtop == ptl_hole(last) )
              mtop=mtop-1
              last=last-1
           enddo
           if( mtop == diag_tracer_n ) diag_tracer_n=m
           sp%phase(1:sp%nphase,m)  = sp%phase(1:sp%nphase,mtop)
           sp%phase0(1:sp%nphase,m) = sp%phase0(1:sp%nphase,mtop)
           sp%gid(m) = sp%gid(mtop)
           if(sml_push_mode/=1) then
              sp%derivs(:,:,m)=sp%derivs(:,:,mtop)
!              sp%reset_derivs(m)=sp%reset_derivs(mtop)
           endif
           mtop=mtop-1
           if( mtop == sp%num) exit FILLHOLE
        enddo FILLHOLE

        !*************** To Right ****************************
        !*****************************************************
        !        print *, 'send_num right',sml_mype
        ! send # of particle to move right
        call monitor_start (SHIFT_IE_SR_R_)
        mrecvleft=0
!        idest=mod(sml_mype+1,sml_totalpe)
!        isource=mod(sml_mype-1+sml_totalpe,sml_totalpe)
        idest=modulo(sml_mype+sml_pe_per_plane,sml_totalpe)
        isource=modulo(sml_mype-sml_pe_per_plane,sml_totalpe)
        isendtag=sml_mype
        irecvtag=isource
        call MPI_SENDRECV(msendright,2,MPI_INTEGER,idest,isendtag,&
             mrecvleft,2,MPI_INTEGER,isource,irecvtag,MPI_COMM_WORLD,istatus,ierror)

        allocate(recvleft(sendsize,max(mrecvleft(1),1)))
        recvleft=0D0 !debug
        allocate(recvleftid(max(mrecvleft(1),1)))
        recvleftid=0 !debug
        
!        print *, 'send_data right',sml_mype,mrecvleft(1),msendright(1)
        ! send particle to right and receive from left
        recvleft=0D0
        isendcount=msendright(1)*sendsize
        irecvcount=mrecvleft(1)*sendsize
        call MPI_SENDRECV(sendright,isendcount,MPI_REAL8,idest,&
             isendtag,recvleft,irecvcount,MPI_REAL8,isource,&
             irecvtag,MPI_COMM_WORLD,istatus,ierror)
#ifdef SMALL_GID        
        call MPI_SENDRECV(sendrightid,msendright(1),MPI_INTEGER,idest,&
             isendtag,recvleftid,mrecvleft(1),MPI_INTEGER,isource,&
             irecvtag,MPI_COMM_WORLD,istatus,ierror)
#else
        call MPI_SENDRECV(sendrightid,msendright(1),MPI_INTEGER8,idest,&
             isendtag,recvleftid,mrecvleft(1),MPI_INTEGER8,isource,&
             irecvtag,MPI_COMM_WORLD,istatus,ierror)
#endif
        call monitor_stop (SHIFT_IE_SR_R_)
        
        
!        call MPI_BARRIER(MPI_COMM_WORLD,ierror)



        !************** To Left *****************************
        !****************************************************
        !        print *, 'send_num left',sml_mype
        ! send # of particle to move left
        call monitor_start (SHIFT_IE_SR_L_)
        mrecvright=0
!        idest=mod(sml_mype-1+sml_totalpe,sml_totalpe)
!        isource=mod(sml_mype+1,sml_totalpe)
        idest=modulo(sml_mype-sml_pe_per_plane,sml_totalpe)
        isource=modulo(sml_mype+sml_pe_per_plane,sml_totalpe)
        isendtag=sml_mype
        irecvtag=isource
        call MPI_SENDRECV(msendleft,2,MPI_INTEGER,idest,isendtag,&
             mrecvright,2,MPI_INTEGER,isource,irecvtag,MPI_COMM_WORLD,istatus,ierror)

        allocate(recvright(sendsize,max(mrecvright(1),1)))
        recvright=0D0 !debug
        allocate(recvrightid(max(mrecvright(1),1)))
        recvrightid=0 !debug

        ! send particle to left and receive from right
        recvright=0.0
        isendcount=msendleft(1)*sendsize
        irecvcount=mrecvright(1)*sendsize
        call MPI_SENDRECV(sendleft,isendcount,MPI_REAL8,idest,&
             isendtag,recvright,irecvcount,MPI_REAL8,isource,&
             irecvtag,MPI_COMM_WORLD,istatus,ierror)

#ifdef SMALL_GID        
        call MPI_SENDRECV(sendleftid,msendleft(1),MPI_INTEGER,idest,&
             isendtag,recvrightid,mrecvright(1),MPI_INTEGER,isource,&
             irecvtag,MPI_COMM_WORLD,istatus,ierror)
#else        
        call MPI_SENDRECV(sendleftid,msendleft(1),MPI_INTEGER8,idest,&
             isendtag,recvrightid,mrecvright(1),MPI_INTEGER8,isource,&
             irecvtag,MPI_COMM_WORLD,istatus,ierror)
#endif

        call monitor_stop (SHIFT_IE_SR_L_)
        
!        call mpi_barrier(mpi_comm_world,ierror) !debug
!        print *, 'tracer',sml_mype !debug

!        print *, '[',sml_mype,']','mrecvl,mrevr',mrecvleft(1),mrecvright(1)

        ! tracer particle
        if( mrecvleft(2) > 0 )then
           diag_tracer_n=sp%num+mrecvleft(2)
        elseif( mrecvright(2) > 0 )then
           diag_tracer_n=sp%num+mrecvleft(1)+mrecvright(2)
        endif
        
        
        ! need extra particle array
        if(sp%num+mrecvleft(1)+mrecvright(1) > sp%maxnum)then
           print *, 'not enough ptl array size', sp%num,sp%num+mrecvleft(1)+mrecvright(1),sp%maxnum,sml_mype, sp%type
           call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
        endif
        
        !**************** Unpack *********************88
        !***********************************************
     
!        print *, 'particle  unpack', sml_mype
        ! unpack particle, particle moved from left
        
        do m=1,mrecvleft(1)
           sp%phase(1:sp%nphase,m+sp%num)=recvleft(1:sp%nphase,m)
           sp%phase0(1:sp%nphase,m+sp%num)=recvleft(sp%nphase+1:sp%nphase2,m)
           if(sml_push_mode/=1) then
              do pc_index=1, sml_pc_order
                 arr_end=sp%nphase2 + sp%nphase*pc_index
                 arr_start=arr_end -sp%nphase +1
                 sp%derivs(:,pc_index,m+sp%num) = recvleft(arr_start:arr_end,m)
              enddo
!              sp%reset_derivs(m+sp%num)=0
           endif
           sp%gid(m+sp%num)=recvleftid(m)
        enddo
        sp%num=sp%num+mrecvleft(1)

        ! particle moved from right
        do m=1,mrecvright(1)
           sp%phase(1:sp%nphase,m+sp%num)=recvright(1:sp%nphase,m)
           sp%phase0(1:sp%nphase,m+sp%num)=recvright(sp%nphase+1:sp%nphase2,m)
           if(sml_push_mode/=1) then
              do pc_index=1, sml_pc_order
                 arr_end=sp%nphase2 + sp%nphase*pc_index
                 arr_start=arr_end -sp%nphase +1
                 sp%derivs(:,pc_index,m+sp%num) = recvright(arr_start:arr_end,m)
              enddo
!              sp%reset_derivs(m+sp%num)=0
           endif
           sp%gid(m+sp%num)=recvrightid(m)
        enddo
        sp%num=sp%num+mrecvright(1)
        

        deallocate(sendleft,sendright,recvleft,recvright)
        deallocate(sendleftid,sendrightid,recvleftid,recvrightid)
        m0=sp%num-mrecvright(1)-mrecvleft(1)+1
!        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
       
     enddo INFLOOP
  enddo !species loop

end subroutine shift_ie

subroutine shift_check(grid,sp)
  use ptl_module
  use grid_class
  use sml_module
  implicit none
  integer :: i
  type(grid_type):: grid
  type(species_type) :: sp
  
  do i= 1, sp%num
     if(grid%phimin > sp%phase(3,i) .or. grid%phimax < sp%phase(3,i) ) then
        print *, 'err in shift check', sp%phase(3,i),grid%phimin, grid%phimax, sml_mype,i,sp%type
     endif
  enddo

end subroutine shift_check



subroutine shift_check2(grid,sp,message,flag)
  use ptl_module
  use grid_class
  use sml_module
  implicit none  
  character (len=30) :: message
  integer :: i,flag
  type(grid_type):: grid
  type(species_type) sp
  integer :: count
  

  count=0
  do i= 1, sp%num
     if(grid%phimin > sp%phase(3,i) .or. grid%phimax < sp%phase(3,i) ) then
        count=count+1
     endif
  enddo
  if(count>0) then
     print *, 'shift check :: ',message, count, sml_mype,sp%type,flag
  endif

end subroutine shift_check2




!!****************************************************************************
!!> derivatives of electron phase
!!
!!<***************************************************************************
subroutine derivs_sp(fld,T, Y,YPRIME,sp_type)
  use fld_module
  use efld_module, only : efld_npsi
  use sml_module
  use ptl_module
  use rpl_module
  implicit none
  type(fld_type), intent(in) :: fld
  integer, intent(in) :: sp_type 
  real (kind=8), intent(in) :: T, Y(ptl_nphase_max) 
  real (kind=8), intent(out) :: YPRIME(ptl_nphase_max)

  real (kind=8):: D, B, nb_curl_nb, dbdr, dbdz, fr,fp,fz,r,z,phi,rho,mu,dbdphi
  real (kind=8):: ripp,dripp_dr,dripp_dz
  integer :: i
  real (kind=8) :: mass, charge,c_m  !relative charge and mass to normalizing charge and mass for multi-species
  !! for weigth calculation
  real (kind=8) :: exb1,exb2,energy,gradn(2),&
       gradt(2),exbgradn, exbgradt, exbgradb
  !! variables for optimization. 1/B, 1/B^2, rho^2, mu*rho^2*B, B^2
  real (kind=8):: over_B, over_B2,cmrho2,cmrho,murho2b,murho2b_c,b2,inv_r
  real (kind=8) :: exb_on,exb_comp
  real (kind=8) :: energy_max, sml_wdot_energy_max, one_m_w
  real (kind=8) :: yprime_df(ptl_nphase_max)
  sml_wdot_energy_max=10D0
  
  mass=ptl_mass(sp_type) !
  charge=ptl_charge(sp_type) !-1.D0/ptl_charge
  exb_on=sml_exb_on(sp_type)
  exb_comp=1D0-exb_on
  c_m=charge/mass
  
  
  r=y(1)
  z=y(2)
  phi=y(3)
  rho=y(4)
  mu=y(5)    
  inv_r=1D0/r
  
  B = sqrt( fld%br**2 + fld%bz**2 + fld%bphi **2 )
  b2=b**2
  over_B=1/B
  over_B2=over_B**2
  
  ! normalized b dot curl of normalized b
  nb_curl_nb= -1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - (fld%bphi/r  + fld%dbpdr)*fld%bz )
  D=1.D0/ ( 1.D0 + rho * nb_curl_nb )
  
  dbdphi=0D0
  !    if(RPL_mode==1) then  ! don't use with delta-f , not implimented yet
  !       call get_ripple(r,z,ripp,dripp_dr,dripp_dz)
  !       dbdphi=B*rpl_N_coil*dsin(rpl_N_coil*phi)*ripp
  !    endif
  
  dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
  dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
  ! bug fix -- rho**2 B grad B term added 2002/1/15
  !! optimization variables
  cmrho2=c_m*rho**2
  cmrho =c_m*rho
  murho2b=(mu+charge*cmrho2 *B)
  murho2b_c=murho2b/charge

  ! F is grad H/q
  fr = exb_on*fld%Er - (murho2B_c) * dbdr ! fr is force over q . -gradient of  Hamiltonian
  fp = exb_on*fld%Ephi - (murho2B_c) * dbdphi/r  ! modified by shlee 5/30/2001 -- /r added
  fz = exb_on*fld%Ez - (murho2B_c) * dbdz
  
  yprime(1)= D*( (fld%bz*Fp - fld%Bphi * Fz) * over_B2         &
       +  cmrho * fld%br                       &
       +  cmrho2 * (fld%dbzdp*inv_r - fld%dbpdz ) )
  yprime(2)= D*( (fld%bphi * fr - fld%br * fp ) * over_B2      &
       +  cmrho * fld%bz                       &
       +  cmrho2 * (fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r) )
  yprime(3)= D*( (fld%br * fz - fld%bz * fr) * over_B2         &
       +  cmrho * fld%bphi                     &
       +  cmrho2 * ( fld%dbrdz - fld%dbzdr) )
  yprime(3)=yprime(3)*inv_r  ! modified by shlee 5/30/2001
  ! modified equation from lagrangian 2001/09  , 2002/01/23 -- 1/D corrected to D
  
  fr = fr + exb_comp*fld%Er
  fp = fp + exb_comp*fld%Ephi
  fz = fz + exb_comp*fld%Ez
  
  
  yprime(4)=D*over_B2 *( &
       fld%br*fr + fld%bz*fz + fld%bphi*fp &
       + rho*( fr*(fld%dbzdp*inv_r-fld%dbpdz) + fz*(fld%bphi*inv_r+fld%dbpdr-fld%dbrdp*inv_r) + fp*(fld%dbrdz-fld%dbzdr)) &
       )
  
  if( sml_deltaf==1 ) then
     yprime_df=yprime !default
     energy=(0.5D0*charge*cmrho2*b2 + mu*b)
     if( sml_canonical_maxwell==1) then  
        exb1= D*( (fld%bz*fld%Ephi - fld%Bphi * fld%Ez) * over_B2 )
        exb2= D*( (fld%bphi * fld%Er - fld%br * fld%Ephi ) * over_B2 ) 
        
        !             exb1= yprime(1) - D*(rho*fld%br)
        !             exb2= yprime(2) - D*(rho*fld%bz)
        !    exb3= D*( (fld%br * fz - fld%bz * fr) * over_B2	    
        !    exb3=exb3/r   ! no need for weight calculation of symmetric background

        gradn(:)=fld%f0_gradn(:) * (1D0 + rho*fld%dIdpsi)
        gradt(:)=fld%f0_gradt(:) * (1D0 + rho * fld%dIdpsi )
        !             dndrhon=ddendpsi *fld%Bphi*r /den  !no need
        !             dtdrhot=dtempdpsi*fld%Bphi*r /temp  !no need
     elseif( sml_dwdt_exb_only==1) then
        exb1= D*( (fld%bz*fld%Ephi - fld%Bphi * fld%Ez) * over_B2 )
        exb2= D*( (fld%bphi * fld%Er - fld%br * fld%Ephi ) * over_B2 ) 
        
        yprime_df(1)=exb1 !+ D* cmrho * fld%br
        yprime_df(2)=exb2 !+ D* cmrho * fld%bz
        yprime_df(3)= inv_r*D*( (fld%br * fld%Ez - fld%bz * fld%Er) * over_B2 )             
        yprime_df(4)=D*over_B2 *( &
             fld%br*fld%Er + fld%bz*fld%Ez + fld%bphi*fld%Ephi &
             + rho*( fld%Er*(fld%dbzdp*inv_r-fld%dbpdz) + &
             fld%Ez*(fld%bphi*inv_r+fld%dbpdr-fld%dbrdp*inv_r) +&
             fld%Ephi*(fld%dbrdz-fld%dbzdr)) &
             )
       
        gradn(:)=fld%f0_gradn(:)
        gradt(:)=fld%f0_gradt(:)

        if(sml_no_para_nonlinear==1) yprime(4)=yprime(4)-yprime_df(4)
     else
        !!add other drift effect in exb: gradb drift and curl drift
        exb1= yprime(1) - D*(cmrho*fld%br)
        exb2= yprime(2) - D*(cmrho*fld%bz) 
        
        gradn(:)=fld%f0_gradn(:)
        gradt(:)=fld%f0_gradt(:)
        !             dndrhon=0d0    !no need
        !             dtdrhot=0d0    !no need
        !             print *, ddendpsi,dtempdpsi,energy/temp-1.5D0!,energy/temp,energy,mass,rho2,mu,b
     endif 
     exbgradn= (exb1*gradn(1) + exb2*gradn(2))/fld%f0_den
     exbgradT= (exb1*gradT(1) + exb2*gradT(2))/fld%f0_Temp
     exbgradB= (exb1*dbdr   + exb2*dbdz)
     
     
     !          yprime(6)= (1./y(8)- y(7))* (&
     energy_max=min( fld%f0_temp*sml_wdot_energy_max, energy )
     
!     yprime(6)= (1.D0 - y(7))* (&
     if(sml_dwdt_fix_bg) then 
        one_m_w=1D0
     else
        one_m_w=1D0 - y(7)
     end if
     yprime(6)= one_m_w* (&
          - exbgradn &
          - exbgradt*(energy_max/fld%f0_temp-1.5D0) &
!          + (yprime(1)*fld%Er + yprime(2)*fld%Ez + yprime(3)*fld%r*fld%Ephi)/fld%f0_temp &          
          + (murho2B*(dbdr*yprime_df(1)+dbdz*yprime_df(2)) + yprime_df(4)*c_m*charge*rho*b2)/fld%f0_temp &
!          + ((yprime(1)-yprime_df(1))*fld%Er + (yprime(2)-yprime_df(2))*fld%Ez + (yprime(3)-yprime_df(3))*fld%r*fld%Ephi)/fld%f0_temp &          
          !      
          !          +  (yprime(4)*c_m*charge_rel*rho*b2 +murho2b*(dbdr* yprime(1)+ dbdz*yprime(2)))/fld%f0_temp &
          ! /
          !               eleminate below term 2003/11/11 energy conservation f0 equation
          !               + exbgradB*(murho2b*B/temp) &
          !               - yprime(4) *(dndrhon  + dtdrhot*(energy/temp -1.5D0) - rho*B*B/temp) &
          )
     
     if(sml_supress_weight_growth) then
        if(abs(y(7))> sml_weight_max) then
           if( y(7)*yprime(6) > 0 ) then
              yprime(6)=0D0
           endif
        endif
     endif

     yprime(7)= yprime(6)
  else
     yprime(6:7)=0.D0
  endif
  
  
  yprime(5)=0D0  ! constant magnetic moment
  yprime(9)=0D0
  yprime(8)=0D0
    
  if(sml_turb_efield==0) then
     yprime(3)=0D0  !debug
  endif
end subroutine derivs_sp



!!****************************************************************************
!!> calculate field of a given point using interpolation funtions
!!  adopted from xorbit
!!
!!  first created : 2000/10/19
!!  last modified : 2006/02/01
!!  B->-B routine added (commented out)
!!  time dependance of field is added (2002/6/19)
!!  2006/02/01 fld module update
!!  2002/09/10 code optimization for speed
!!  2002/11/18 code modification for gxc
!!****************************************************************************
subroutine field(fld,t)
    use fld_module
    use ptl_module, only : ptl_rz_outside
    use sml_module, only :sml_minusB
    use eq_module, only : eq_min_r,eq_max_r,eq_min_z,eq_max_z, eq_x_psi, eq_x_z
    implicit none
    type(fld_type) :: fld  !! Field information
    real (kind=8), intent(in) :: t  !! time
    real (kind=8) :: psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z
    real (kind=8) :: r,z
    real (kind=8) , external :: psi_interpol
    real (kind=8) , external :: I_interpol
    real (kind=8) :: r2, over_r,over_r2 !! variables for opimization


    r=fld%r
    z=fld%z
    r2=r**2
    over_r=1/r
    over_r2=over_r**2
    if(r<eq_min_r) then
	r=eq_min_r
	ptl_rz_outside=1
    else if (r>eq_max_r)then
	r=eq_max_r
	ptl_rz_outside=1
    endif
    if(z<eq_min_z) then
	z=eq_min_z
	ptl_rz_outside=1
    else if (z>eq_max_z)then
	z=eq_max_z
	ptl_rz_outside=1
    endif

    psi = psi_interpol(r,z,0,0)
    dpsi_dr = psi_interpol(r,z,1,0)
    dpsi_dz = psi_interpol(r,z,0,1)
    d2psi_d2r= psi_interpol(r,z,2,0)
    d2psi_drdz= psi_interpol(r,z,1,1)
    d2psi_d2z= psi_interpol(r,z,0,2)

    fld%psi=psi
    fld%dpsidr=dpsi_dr
    fld%dpsidz=dpsi_dz
    ! added 2001/06/01  - lower bound of psi
    if(psi<0D0) then
	psi=0D0
    endif

    !fld_q=q_interpol(psi,0)  --> no need
    
    if(psi<eq_x_psi .AND. z<eq_x_z) then
	fld%I=I_interpol(psi,0,3)
	fld%dIdpsi = I_interpol(psi,1,3)
    else
        fld%I=I_interpol(psi,0,1)
	fld%dIdpsi = I_interpol(psi,1,1)
    endif

    
    fld%br=- dpsi_dz *over_r
    fld%bz= dpsi_dr *over_r
    fld%bphi=fld%I *over_r
	
    !derivativs
    fld%dbrdr=dpsi_dz *over_r2 - d2psi_drdz *over_r
    fld%dbrdz=- d2psi_d2z *over_r
    fld%dbrdp=0D0
    
    fld%dbzdr= - dpsi_dr * over_r2 + d2psi_d2r *over_r
    fld%dbzdz= d2psi_drdz *over_r
    fld%dbzdp=0D0
           
    fld%dbpdr= dpsi_dr * fld%dIdpsi *over_r - fld%I *over_r2
    fld%dbpdz= fld%dIdpsi * dpsi_dz *over_r
    fld%dbpdp=0D0
!    call efield(r,z,fld_psi,t)  ! set fld_er, ez, ephi value for a give r,z value  -> move to derivs

    if (sml_minusB==1) then
       ! set minus B -----------
       fld%br=-fld%br
       fld%bz=-fld%bz
       fld%bphi=-fld%bphi
       
       fld%dbrdr=-fld%dbrdr
       fld%dbrdz=-fld%dbrdz
       fld%dbrdp=-fld%dbrdp
       
       fld%dbzdr=-fld%dbzdr
       fld%dbzdz=-fld%dbzdz
       fld%dbzdp=-fld%dbzdp
       
       fld%dbpdr=-fld%dbpdr
       fld%dbpdz=-fld%dbpdz
       fld%dbpdp=-fld%dbpdp
       !---------------------------
       
    endif


end subroutine



!-------------- obsolete routines----------------------------

!!*********************************************************
!!> routine from NRF
!!  modified a little. (program style only)
!!<*********************************************************
!!$SUBROUTINE single_RK4(y,dydx,X,H,YOUT,DERIVS)
!!$  use ptl_module, only : ptl_maxnum, ptl_nphase
!!$!  use efld_module, only : efld_npsi,efld_dpdp,efld_mode
!!$  !use sml_module
!!$  implicit none
!!$  external derivs
!!$  
!!$  real (kind=8), intent(inout) :: h,x, yout(ptl_nphase),y(ptl_nphase)
!!$  real (kind=8) :: dydx(ptl_nphase)
!!$  real (kind=8) :: h6, hh, xh,dym(ptl_nphase), dyt(ptl_nphase),yt(ptl_nphase)
!!$!  save dydx,dym,dyt,yt
!!$  
!!$  HH=H*0.5D0
!!$  H6=H/6D0
!!$  XH=X+HH
!!$  !write(1004,*) sml_istep, sml_time, sml_dt, x,h,xh
!!$  !added to get dydx in this routine
!!$!  call derivs(sp,x,y,dydx)  
!!$  !diagnosis -- save it for fast diagnosis
!!$!     diag_derivs(i,:)=dydx
!!$  
!!$  yt=y + hh* dydx
!!$
!!$  CALL DERIVS(XH,YT,DYT)
!!$  yt=y+hh*dyt
!!$
!!$  CALL DERIVS(XH,YT,DYM)
!!$  yt=y+h*dym
!!$  dym=dyt + dym
!!$  
!!$  CALL DERIVS(X+H,YT,DYT)
!!$  yout=y + h6 * ( dydx+dyt+2D0*dym )
!!$
!!$END subroutine
!!$
!! Evaluating ripple field strength.
subroutine get_ripple(r_in,z_in,ripp,drippdr,drippdz) 
  use rpl_module 
  implicit none
  real (kind=8), intent(in) :: r_in,z_in
  real (kind=8), intent(out) :: ripp,drippdr,drippdz
  real (kind=8) :: tau
  
  !real (kind=8)  :: rpl_N_coil, rpl_R0, rpl_tau0, rlp_ratio, rpl_elong
  !r0 : center of elipse
  !elong : elongation of elipse
  !tau0 : increase rate
  !ratio : relative strength to B at the center of elipse
  !ripple model
  ! ripple strenght = ratio*B*exp(tau/tau0)
  ! tau= (R-R0)**2 + elong*z**2
  ! ripp = ripple strength / B

  tau=sqrt( (r_in - rpl_R0)**2 + rpl_elong*(z_in)**2 )
  ripp= rpl_ratio*exp(tau/rpl_tau0)
  
  drippdr=0D0 !incomplete
  drippdz=0D0 !incomplete


end subroutine get_ripple
!!$
!!$!! Only for test the accuracy of Euyler method
!!$SUBROUTINE single_RK1(y,dydx,X,H,YOUT,DERIVS)
!!$  use ptl_module, only : ptl_maxnum, ptl_nphase
!!$!  use efld_module, only : efld_npsi,efld_dpdp,efld_mode
!!$  !use sml_module
!!$  implicit none
!!$  external derivs
!!$  
!!$  real (kind=8), intent(inout) :: h,x, yout(ptl_nphase),y(ptl_nphase)
!!$  real (kind=8) :: dydx(ptl_nphase)
!!$  real (kind=8) :: h6, hh, xh,dym(ptl_nphase), dyt(ptl_nphase),yt(ptl_nphase)
!!$!  save dydx,dym,dyt,yt
!!$  
!!$  HH=H*0.5D0
!!$  H6=H/6D0
!!$  XH=X+HH
!!$  
!!$  yout=y + h* dydx
!!$
!!$END subroutine
!!$
!!$!! Only for test the accuracy of RK2.
!!$SUBROUTINE single_RK2(y,dydx,X,H,YOUT,DERIVS)
!!$  use ptl_module, only : ptl_maxnum, ptl_nphase
!!$!  use efld_module, only : efld_npsi,efld_dpdp,efld_mode
!!$  !use sml_module
!!$  implicit none
!!$  external derivs
!!$  
!!$  real (kind=8), intent(inout) :: h,x, yout(ptl_nphase),y(ptl_nphase)
!!$  real (kind=8) :: dydx(ptl_nphase)
!!$  real (kind=8) :: h6, hh, xh,dym(ptl_nphase), dyt(ptl_nphase),yt(ptl_nphase)
!!$!  save dydx,dym,dyt,yt
!!$  
!!$  HH=H*0.5D0
!!$  H6=H/6D0
!!$  XH=X+HH
!!$  !write(1004,*) sml_istep, sml_time, sml_dt, x,h,xh
!!$  !added to get dydx in this routine
!!$!  call derivs(sp,x,y,dydx)  
!!$  !diagnosis -- save it for fast diagnosis
!!$!     diag_derivs(i,:)=dydx
!!$  
!!$  yt=y + hh* dydx
!!$
!!$  call derivs(xh,yt,dym)
!!$  
!!$  yout=y+ h*dym
!!$!  CALL DERIVS(XH,YT,DYT)
!!$!  yt=y+hh*dyt
!!$
!!$!  CALL DERIVS(XH,YT,DYM)
!!$!  yt=y+h*dym
!!$!  dym=dyt + dym
!!$  
!!$!  CALL DERIVS(X+H,YT,DYT)
!!$!  yout=y + h6 * ( dydx+dyt+2D0*dym )
!!$
!!$
!!$END subroutine


subroutine modulo_ie(grid,ptl)
  use sml_module
  use ptl_module
  use grid_class
  use diag_module
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(ptl_type) :: ptl
  type(species_type),pointer :: sp
  integer i,m
  target ptl

  if(sp%type==1) then
     sp=>ptl%ion
  else
     sp=>ptl%elec
  endif

  ! Modulo operation
  do i=1, sp%num
     if(sp%phase(3,i)>= sml_2pi .or. sp%phase(3,i)< 0D0 ) then
        sp%phase(3,i)=modulo(sp%phase(3,i),sml_2pi)
     endif
  enddo
  
end subroutine modulo_ie


subroutine remove_particle(sp,i,flag)
  use ptl_module
  use sml_module
  implicit none
  integer, intent(in) :: i,flag
  type(species_type) :: sp
  !error diagnosis
  real (kind=8) :: b_interpol, b,ekin,pitch
  
  if(sp%gid(i) < 0 ) return
  sp%gid(i)=-sp%gid(i)
  if(sp%lost_num < sp%lost_nummax) then
     sp%lost_num=sp%lost_num + 1
     sp%lost_index(sp%lost_num)=i
  endif  

  if(flag==-2 .and. sml_mype==0) then
     b=b_interpol(sp%phase(1,i),sp%phase(2,i),0D0)
     call rho_mu_to_ev_pitch2(sp%phase(4,i),sp%phase(5,i),b,ekin,pitch,sp%type)
     write(400+sp%type,*) sp%phase(1,i), sp%phase(2,i),&
          ekin, pitch,&          
          sp%phase0(1,i),sp%phase0(2,i) 
  endif

  if(flag==-1 .and. sml_mype==0) then
     write(450+sp%type,*) sp%phase(1,i), sp%phase(2,i),&
          ekin, pitch,&          
          sp%phase0(1,i),sp%phase0(2,i) 
  endif
end subroutine remove_particle

subroutine memory_cleaning(ptl)
  use ptl_module
  use neu_module
  use rf_module
  use sml_module, only :sml_mype,sml_angle_stored, sml_electron_on
  implicit none
  type(ptl_type) :: ptl
  type(species_type), pointer :: sp
  integer :: i,j,ptmp
  integer :: sp_type 
  real (kind=8) :: y(ptl_nphase_max)
  target ptl

  do sp_type=1, 1+sml_electron_on
     if(sp_type==1) then
        sp=>ptl%ion
     else
        sp=>ptl%elec
     endif
     
     if(sp%lost_num<=sp%lost_nummax) then
        ! algorithm 
        ! BASIC : for every empty space, move last particle to this 
        ! EXCEPTION : last particle can be lost particle
        !           : If so, eliminate the particle and find next non-empty space
        !           : When the index indicate erased space, just ignore it
        do j=1, sp%lost_num
           !find non-empty last particle
           do while( sp%gid(sp%num)<0 ) 
              y(:)=sp%phase(:,sp%num)              
              ptmp=sp%num
              call eliminate_one(ptmp,y,sp%gid(ptmp),sp_type)
              sp%num=sp%num-1
              if(sp%num < 1) exit
           enddo
           
           i=sp%lost_index(j)

           ! ignore erased space 
           ! i cannot be sp%num because sp%num has positive gid value
           if( i < sp%num ) then
              y(:)=sp%phase(:,i)
              call eliminate_one(i,y,sp%gid(i),sp_type)

              !move last particle to i-th memory
              sp%phase(:,i)=sp%phase(:,sp%num)
              sp%phase0(:,i)=sp%phase0(:,sp%num)
              sp%gid(i)=sp%gid(sp%num)
              sp%derivs(:,:,i)=sp%derivs(:,:,sp%num)
!              sp%reset_derivs(i)=sp%reset_derivs(sp%num)              
              if( sml_angle_stored==1 .and. sp_type==1) sp%angle(i)=sp%angle(sp%num)
              if (rf_on==1 .and. sp_type==1) rf_r_save(i)=rf_r_save(sp%num)
              sp%num=sp%num-1
           endif
        enddo
     else
        !sequencial search whole particle
        i=1
        do while ( i < sp%num)
           if(sp%gid(i) < 0) then
              !calculate weight summation of lost particle
              y(:)=sp%phase(:,i)
              call eliminate_one(i,y,sp%gid(i),sp_type)
              !eliminate i-th particle and move last particle to i-th memory
              sp%phase(:,i) =sp%phase(:,sp%num)
              sp%phase0(:,i)=sp%phase0(:,sp%num)
              sp%gid(i)=sp%gid(sp%num)
              if( sml_angle_stored==1 .and. sp_type==1) sp%angle(i)=sp%angle(sp%num)
              if (rf_on==1 .and. sp_type==1) rf_r_save(i)=rf_r_save(sp%num)        
              sp%num=sp%num-1
           else
              ! next particle
              i=i+1        
           endif
        enddo
        
        !last particle exception
        if(sp%gid(sp%num) < 0 ) then
           y(:)=sp%phase(:,sp%num)
           ptmp=sp%num
           call eliminate_one(ptmp,y,sp%gid(ptmp),sp_type)
           sp%num=sp%num-1
        endif
     endif

     ! debug - check gid again
!!$     do i=1, pnum
!!$        if(gid(i)<0) then
!!$           print *, 'error in memory cleaning',i,pnum, sml_mype
!!$           stop
!!$        endif
!!$     enddo
     sp%lost_num=0
  enddo

end subroutine memory_cleaning

! Particle elimination - side effect calculation
subroutine eliminate_one(i,y,gid,sp_type)
  use neu_module
  use ptl_module
  implicit none
  integer,intent(in) :: i,sp_type
#ifdef SMALL_GID
  integer, intent(in) :: gid
#else
  integer (kind=8), intent(in) :: gid
#endif
  real (kind=8) :: y(ptl_nphase_max)

  if(gid<0 .and. sp_type==1) then  ! for full particle simulation only
     neu_weight_sum_lost=neu_weight_sum_lost + y(8)
  endif


end subroutine eliminate_one

subroutine f0_info(grid,i,sp,fld)
  use sml_module
  use grid_class
  use fld_module
  use ptl_module
  use eq_module, only : eq_x_psi
  implicit none
  type(grid_type) :: grid
  type(fld_type) :: fld
  integer, intent(in) :: i
  type(species_type) :: sp
  integer :: ip,nodes(3),init,itr,count
  real (kind=8) :: gradn(2),p(3),dtempdpsi,ddendpsi,x(2),den,psi,z
  !! f0 ion temp, initial ion density, dT/dpsi, dn/dpsi
  real (kind=8), external :: f0_temp_ev, f0_den,df0_temp_ev_dpsi,df0_den_dpsi 
  real (kind=8) :: tmp
  real (kind=8) :: sml_f0_psi_c, sml_f0_1_psi_w, alpha

  logical, parameter :: USE_SEARCH_TR2 = .true.

  sml_f0_psi_c=0.5*(sqrt(sml_inpsi)+sqrt(sml_outpsi))
  sml_f0_1_psi_w=1D0/( 0.35*(sqrt(sml_outpsi)-sqrt(sml_inpsi)) )

  psi=fld%psi
  z=fld%z

  if(sml_deltaf_f0_mode==0 .or. fld%psi < eq_x_psi .and. sml_deltaf_f0_mode==1) then
     den=f0_den(psi,z)
     ddendpsi=df0_den_dpsi(psi,z)
     gradn(1)=ddendpsi * fld%dpsidr
     gradn(2)=ddendpsi * fld%dpsidz
  
     !temperature
     fld%f0_temp=f0_temp_ev(psi,z)* sml_ev2j
     dtempdpsi=df0_temp_ev_dpsi(psi,z)*sml_ev2j

  else if(sml_deltaf_f0_mode==-1) then
     den=f0_den(psi,z)
     tmp=1D0/sqrt(fld%dpsidr**2+fld%dpsidz**2)  ! drdpsi 
     alpha= exp( - ((sqrt(psi) - sml_f0_psi_c)*sml_f0_1_psi_w )**6 )
     ddendpsi=den*(-sml_f0_1_Ln)*tmp*alpha
     gradn(1)=ddendpsi * fld%dpsidr
     gradn(2)=ddendpsi * fld%dpsidz

     !temperature
     fld%f0_temp=f0_temp_ev(psi,z)* sml_ev2j
     dtempdpsi=fld%f0_temp*(-sml_f0_1_Lt)*tmp*alpha
          
  else
     den=0D0
     if(sp%type==1) then
        ! Ion -----------------------------------------------------------------
        init=sp%tr_save(1,1,i)
        if(init > 0 ) then
           x(1)=fld%r
           x(2)=fld%z

           if (USE_SEARCH_TR2) then
             call search_tr2(grid,x,itr,p)
           else
             call search_tr_with_guess(grid,x,init,itr,p,count)
           endif

           nodes=grid%nd(:,itr)
           do ip=1,3
              den=den + p(ip) * sml_f0_n(nodes(ip))
           enddo
           gradn(1:2)=sml_f0_gradn(1:2,itr)
        else
           print *, 'f0_info ion invalid tr'
        endif
     else
        !Electron --------------------------------------------------------------
        ! parallel field
        if(sp%tr_save(1,1,i)>0) then
           nodes=grid%nd(:,sp%tr_save(1,1,i))
           do ip=1,3
              den=den + sp%p_save(ip,1,1,i) *  sml_f0_n(nodes(ip))
           enddo
           gradn(1:2)=sml_f0_gradn(1:2,sp%tr_save(1,1,i))
        else
           print *, 'f0_info electron invalid tr', i,sp%tr_save(1,1,i)
        endif
     endif

     !temperature - flux function f0
     fld%f0_temp=f0_temp_ev(psi,z)* sml_ev2j
     dtempdpsi=df0_temp_ev_dpsi(psi,z)*sml_ev2j
  endif
  
  fld%f0_den=den
  fld%f0_gradn(:)=gradn(:)
  
  fld%f0_gradT(1)=dtempdpsi * fld%dpsidr
  fld%f0_gradT(2)=dtempdpsi * fld%dpsidz
end subroutine

real (kind=8) function f0_den2(grid,r,z,psi)
  use grid_class
  use sml_module
  use eq_module, only : eq_x_psi
  implicit none
  type(grid_type):: grid
  integer :: itr,ip,nodes(3)
  real (kind=8), intent(in) :: r,z,psi
  real (kind=8) :: x(2),p(3),den
  real (kind=8), external :: f0_den

  logical, parameter :: USE_SEARCH_TR2 = .true.

  if(sml_deltaf_f0_mode==0 .or. psi < eq_x_psi .and. sml_deltaf_f0_mode==1) then
     den=f0_den(psi,z)
  else
     den=0D0
     x(1)=r
     x(2)=z

     if (USE_SEARCH_TR2) then
       call search_tr2(grid,x,itr,p)
     else
       call search_tr(grid,x,itr,p)
     endif

     if(itr>0) then
        nodes=grid%nd(:,itr)
        do ip=1,3
           den=den + p(ip) * sml_f0_n(nodes(ip))
        enddo
     else
        den=f0_den(psi,z)
     endif
  endif
  
  f0_den2=den
end function f0_den2


! initialize f0 grid data
! it should be called after grid & smooth initialization 
! and before load phase
subroutine init_f0(grid)
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  integer :: i, nodes(3)
  real (kind=8) :: r,z,psi,dp1,dp2
  real (kind=8), external :: f0_den,psi_interpol
  ! memory allocation
  allocate(sml_f0_n(grid%nnode),sml_f0_gradn(2,grid%ntriangle))
  
  ! assign data
  do i=1, grid%nnode
     r=grid%x(1,i)
     z=grid%x(2,i)
     psi=psi_interpol(r,z,0,0)

     sml_f0_n(i)=f0_den(psi,z)
  enddo
  
  !density vanish smoothly on the wall
  do i=grid%bd%out1%start,grid%bd%out1%end
     sml_f0_n(i)=1D15
  enddo
  


  !calculate gradient
  do i=1, grid%ntriangle
     nodes(:)=grid%nd(:,i)
     dp1=sml_f0_n(nodes(1))-sml_f0_n(nodes(3))
     dp2=sml_f0_n(nodes(2))-sml_f0_n(nodes(3))
     sml_f0_gradn(:,i)= -( dp1*grid%mapping(1,1:2,i) + dp2*grid%mapping(2,1:2,i) )
  enddo
  
end subroutine init_f0


subroutine shift_pc_index(index,n)
  implicit none
  integer, intent(in) :: n
  integer :: index(n)
  integer :: tmp(n)

  tmp(2:n) = index(1:n-1)
  tmp(1)=index(n)
  
  index=tmp

end subroutine shift_pc_index

