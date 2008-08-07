program xydiag

! E X B history
           open(444,file='sheareb.out',status='old',position='append')
           
! restart, copy old data to backup file to protect previous run data
           open(555,file='history.out',status='old')
! # of run
           read(555,101)j
           
           write(cdum,'("histry",i1,".bak")')j-10*(j/10)
           open(ihistory,file=cdum,status='replace')
           write(ihistory,101)j
           do i=1,3
              read(555,101)j
              write(ihistory,101)j
           enddo
! # of time step
           read(555,101)j
           write(ihistory,101)j
           do i=0,(mquantity+mflux+4*num_mode)*j !i=0 for tstep
              read(555,102)dum
              write(ihistory,102)dum
           enddo
           close(555)
           close(ihistory)

! if everything goes OK, copy old data to history.out
           open(666,file=cdum,status='old')
           open(ihistory,file='history.out',status='replace')            
           read(666,101)j
           irun=j+1
           write(ihistory,101)irun               
           do i=1,3
              read(666,101)j
              write(ihistory,101)j
           enddo
           read(666,101)mstepall
           write(ihistory,101)mstepall+mstep/ndiag
           do i=0,(mquantity+mflux+4*num_mode)*mstepall
              read(666,102)dum
           .   write(ihistory,102)dum
           enddo
           close(666)
           mstepall=mstepall*ndiag
        endif

end program
