!"fft_cray.f90": 1D FFT
subroutine fftr1d(isign,irank,scale,x,y,icount)

  use precision
  implicit none
  integer isign,irank,icount
  real(wp) x(0:irank-1),scale
  complex(wp) y(0:irank/2)

  integer,dimension(3),SAVE :: initial=(/0,0,0/)
  real(doubleprec),dimension(:),allocatable,SAVE :: table1,table2,table3
  real(singleprec),dimension(:),allocatable :: worksp
  real(doubleprec),dimension(:),allocatable :: workdp

  if(wp==singleprec)then
     allocate(worksp(4+4*irank))
  else
     allocate(workdp(4+4*irank))
  endif

  if(icount==1)then
! initialize cray fft
     if(initial(icount)/=1)then
        initial(icount)=1
        if(.not.ALLOCATED(table1))allocate(table1(100+4*irank))
        if(wp==singleprec)then
           call scfft(0,irank,scale,x,y,table1,worksp,0)
        else
           call dzfft(0,irank,scale,x,y,table1,workdp,0)
        endif
     endif
! forward and backward 1d real fft
     if(isign==1)then
        if(wp==singleprec)then
           call scfft(1,irank,scale,x,y,table1,worksp,0)
        else
           call dzfft(1,irank,scale,x,y,table1,workdp,0)
        endif
     else
        if(wp==singleprec)then
           call csfft(-1,irank,scale,y,x,table1,worksp,0)
        else
           call zdfft(-1,irank,scale,y,x,table1,workdp,0)
        endif
     endif

  elseif(icount==2)then
     if(initial(icount)/=1)then
        initial(icount)=1
        if(.not.ALLOCATED(table2))allocate(table2(100+4*irank))
        if(wp==singleprec)then
           call scfft(0,irank,scale,x,y,table2,worksp,0)
        else
           call dzfft(0,irank,scale,x,y,table2,workdp,0)
        endif
     endif
     if(isign==1)then
        if(wp==singleprec)then
           call scfft(1,irank,scale,x,y,table2,worksp,0)
        else
           call dzfft(1,irank,scale,x,y,table2,workdp,0)
        endif
     else
        if(wp==singleprec)then
           call csfft(-1,irank,scale,y,x,table2,worksp,0)
        else
           call zdfft(-1,irank,scale,y,x,table2,workdp,0)
        endif
     endif

  elseif(icount==3)then
     if(initial(icount)/=1)then
        initial(icount)=1
        if(.not.ALLOCATED(table3))allocate(table3(100+4*irank))
        if(wp==singleprec)then
           call scfft(0,irank,scale,x,y,table3,worksp,0)
        else
           call dzfft(0,irank,scale,x,y,table3,workdp,0)
        endif
     endif
     
     if(isign==1)then
        if(wp==singleprec)then
           call scfft(1,irank,scale,x,y,table3,worksp,0)
        else
           call dzfft(1,irank,scale,x,y,table3,workdp,0)
        endif
     else
        if(wp==singleprec)then
           call csfft(-1,irank,scale,y,x,table3,worksp,0)
        else
           call zdfft(-1,irank,scale,y,x,table3,workdp,0)
        endif
     endif
  endif
  
  if(wp==singleprec)then
     deallocate(worksp)
  else
     deallocate(workdp)
  endif

end subroutine fftr1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fftc1d(isign,irank,scale,x)
  
  use precision
  implicit none
  integer :: isign,irank
  real(wp) scale
  complex(wp) x(0:irank-1)

  integer,SAVE :: initial=0
  real(doubleprec),dimension(:),allocatable,SAVE :: table
  real(singleprec),dimension(:),allocatable :: worksp
  real(doubleprec),dimension(:),allocatable :: workdp

  if(initial/=1)then
     initial=1
     if(.not.ALLOCATED(table))allocate(table(100+8*irank))
     if(wp==singleprec)then
        allocate(worksp(1))
        call ccfft(0,irank,scale,x,x,table,worksp,0)
        deallocate(worksp)
     else
        allocate(workdp(1))
        call zzfft(0,irank,scale,x,x,table,workdp,0)
        deallocate(workdp)
     endif
  endif

  if(wp==singleprec)then
     allocate(worksp(8*irank))
     call ccfft(isign,irank,scale,x,x,table,worksp,0)
     deallocate(worksp)
  else
     allocate(workdp(8*irank))
     call zzfft(isign,irank,scale,x,x,table,workdp,0)
     deallocate(workdp)
  endif

end subroutine fftc1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
