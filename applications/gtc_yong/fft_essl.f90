!"fft_essl.f90": 1D FFT
!
! S. Ethier Updated on 08/07/2001
!    The fftr1d and fftc1d routines are used as interfaces between GTC
!    and IBM's ESSL library calls.
!  see:  http://hpcf.nersc.gov/software/ibm/essl/essl67.html
!
subroutine fftr1d(isign,irank,scale,x,y,icount)
  use precision
  implicit none
  integer isign,irank,icount,naux11,naux21
  real(doubleprec),dimension(:),allocatable :: aux11,aux21
  real(doubleprec) :: aux3(1)
  real(wp) x(irank),scale
  complex(wp) y(irank/2+1)

  if(icount>=1 .and. icount<=3)then
    if (irank<=16384) then
       naux11=25000
       naux21=20000
    else
       naux11=20000+0.82*irank  !ESSL requirement
       naux21=20000+0.57*irank
    endif
    allocate(aux11(naux11))
    allocate(aux21(naux21))
! forward and backward 1d real fft
    if(wp==singleprec)then
      if(isign==1)then
        call srcft(1,x,0,y,0,irank,1,1,scale,aux11,naux11,aux21,naux21,aux3,1)
        call srcft(0,x,0,y,0,irank,1,1,scale,aux11,naux11,aux21,naux21,aux3,1)
      else
        call scrft(1,y,0,x,0,irank,1,-1,scale,aux11,naux11,aux21,naux21,aux3,1)
        call scrft(0,y,0,x,0,irank,1,-1,scale,aux11,naux11,aux21,naux21,aux3,1)
      endif
    else   !Call double precision version of the FFT routines
      if(isign==1)then
        call drcft(1,x,0,y,0,irank,1,1,scale,aux11,naux11,aux21,naux21,aux3,1)
        call drcft(0,x,0,y,0,irank,1,1,scale,aux11,naux11,aux21,naux21,aux3,1)
      else
        call dcrft(1,y,0,x,0,irank,1,-1,scale,aux11,naux11,aux21,naux21,aux3,1)
        call dcrft(0,y,0,x,0,irank,1,-1,scale,aux11,naux11,aux21,naux21,aux3,1)
      endif
    endif
    deallocate(aux11,aux21)

  endif
  
end subroutine fftr1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fftc1d(isign,irank,scale,x)
  
  use precision
  implicit none
  integer initial,isign,irank,naux1,naux2,icheck
  real(doubleprec),dimension(:),allocatable :: aux1,aux2
  real(wp) scale
  complex(wp) x(irank)

  if (irank<=8192) then
     naux1=20000
     naux2=20000
  else
     naux1=20000+1.14*irank  !ESSL requirement
     naux2=20000+1.14*irank
  endif
  allocate(aux1(naux1))
  allocate(aux2(naux2))
  if(wp==singleprec)then
     call scft(1,x,1,0,x,1,0,irank,1,1,scale,aux1,naux1,aux2,naux2)
     call scft(0,x,1,0,x,1,0,irank,1,isign,scale,aux1,naux1,aux2,naux2)
  else
     call dcft(1,x,1,0,x,1,0,irank,1,1,scale,aux1,naux1,aux2,naux2)
     call dcft(0,x,1,0,x,1,0,irank,1,isign,scale,aux1,naux1,aux2,naux2)
  endif
  deallocate(aux1,aux2)

end subroutine fftc1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
