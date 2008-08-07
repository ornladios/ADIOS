program xyrestart
implicit none
integer:: ihistory,j,mquantity
ihistory
    open(777,file='history_restart.out',status='unknown')
     rewind(ihistory)
   !!!  read(ihistory,101)j
   !!!  write(777,101)j
   !!!  read(ihistory,101)mquantity
   !!!  write(777,101)mquantity
   !!!  read(ihistory,101)mflx
   !!!  write(777,101)mflx
   !!!  read(ihistory,101)n_mode
   !!!  write(777,101)n_mode
   !!!  read(ihistory,101)mstepfinal
   !!!  noutputs=mstepfinal-mstep/ndiag+istep/ndiag
   !!!  write(777,101)noutputs
   !!!  do i=0,(mquantity+mflx+4*n_mode)*noutputs
   !!!     read(ihistory,102)dum
   !!!     write(777,102)dum
   !!!  enddo
   !!!  close(777)

   ! Now do sheareb.out
   !!!  open(777,file='sheareb_restart.out',status='unknown')
   !!!  rewind(444)
   !!!  read(444,101)j
   !!!  write(777,101)j
   !!!  do i=1,mpsi*noutputs
   !!!     read(444,102)dum
   !!!     write(777,102)dum
   !!!   enddo
   !!!  close(777)
!  endif

101 format(i6)
102 format(e12.6)
