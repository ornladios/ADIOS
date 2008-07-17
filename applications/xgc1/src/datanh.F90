      real (kind=8) function datanh (x)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      implicit none

! arguments
      real (kind=8) x

! local variables
      integer nterms
      real (kind=8) atnhcs(27), dxrel, sqeps, y

!
! series for atnh       on the interval  0.          to  2.50000e-01
!                                        with weighted error   6.86e-32
!                                         log weighted error  31.16
!                               significant figures required  30.00
!                                    decimal places required  31.88
!
      data atnhcs(  1) / +.9439510239319549230842892218633d-1      /
      data atnhcs(  2) / +.4919843705578615947200034576668d-1      /
      data atnhcs(  3) / +.2102593522455432763479327331752d-2      /
      data atnhcs(  4) / +.1073554449776116584640731045276d-3      /
      data atnhcs(  5) / +.5978267249293031478642787517872d-5      /
      data atnhcs(  6) / +.3505062030889134845966834886200d-6      /
      data atnhcs(  7) / +.2126374343765340350896219314431d-7      /
      data atnhcs(  8) / +.1321694535715527192129801723055d-8      /
      data atnhcs(  9) / +.8365875501178070364623604052959d-10     /
      data atnhcs( 10) / +.5370503749311002163881434587772d-11     /
      data atnhcs( 11) / +.3486659470157107922971245784290d-12     /
      data atnhcs( 12) / +.2284549509603433015524024119722d-13     /
      data atnhcs( 13) / +.1508407105944793044874229067558d-14     /
      data atnhcs( 14) / +.1002418816804109126136995722837d-15     /
      data atnhcs( 15) / +.6698674738165069539715526882986d-17     /
      data atnhcs( 16) / +.4497954546494931083083327624533d-18     /
      data atnhcs( 17) / +.3032954474279453541682367146666d-19     /
      data atnhcs( 18) / +.2052702064190936826463861418666d-20     /
      data atnhcs( 19) / +.1393848977053837713193014613333d-21     /
      data atnhcs( 20) / +.9492580637224576971958954666666d-23     /
      data atnhcs( 21) / +.6481915448242307604982442666666d-24     /
      data atnhcs( 22) / +.4436730205723615272632320000000d-25     /
      data atnhcs( 23) / +.3043465618543161638912000000000d-26     /
      data atnhcs( 24) / +.2091881298792393474047999999999d-27     /
      data atnhcs( 25) / +.1440445411234050561365333333333d-28     /
      data atnhcs( 26) / +.9935374683141640465066666666666d-30     /
      data atnhcs( 27) / +.6863462444358260053333333333333d-31     /
!
      data nterms, dxrel, sqeps / 0, 2*0.d0 /
!
! external functions
      integer initds
      real (kind=8) dcsevl, d1mach
      external d1mach, dcsevl, initds

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (nterms.ne.0) go to 10
      nterms = initds (atnhcs, 27, 0.1*sngl(d1mach(3)) )
      dxrel = sqrt (d1mach(4))
      sqeps = sqrt (3.0d0*d1mach(3))
!
 10      y = dabs(x)
      if (y.ge.1.d0) then
        write(0,*) "datanh  dabs(x) ge 1"
        stop
      endif
!
      if (1.d0-y.lt.dxrel) then
        write(0,*) "datanh  answer lt half precision because dabs(x) too near 1"
      endif
!
      datanh = x
      if (y.gt.sqeps .and. y.le.0.5d0) &
        datanh = x*(1.0d0 + dcsevl (8.d0*x*x-1.d0, atnhcs, nterms) )
      if (y.gt.0.5d0) datanh = 0.5d0*log ((1.0d0+x)/(1.0d0-x) )
!
      return
      end function datanh
    
