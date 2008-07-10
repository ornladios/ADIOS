      real (kind=8) function derf (x)
! july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
! last revised june 1983.  w. fullerton, imsl, houston.
      implicit none

! arguments
      real (kind=8) x

! local variables
      integer nterf
      real (kind=8) erfcs(21), sqeps, sqrtpi, xbig, y

!
! series for erf        on the interval  0.          to  1.00000e+00
!                                        with weighted error   1.28e-32
!                                         log weighted error  31.89
!                               significant figures required  31.05
!                                    decimal places required  32.55
!
      data erfcs(  1)  / -.49046121234691808039984544033376d-1     /
      data erfcs(  2)  / -.14226120510371364237824741899631d+0     /
      data erfcs(  3)  / +.10035582187599795575754676712933d-1     /
      data erfcs(  4)  / -.57687646997674847650827025509167d-3     /
      data erfcs(  5)  / +.27419931252196061034422160791471d-4     /
      data erfcs(  6)  / -.11043175507344507604135381295905d-5     /
      data erfcs(  7)  / +.38488755420345036949961311498174d-7     /
      data erfcs(  8)  / -.11808582533875466969631751801581d-8     /
      data erfcs(  9)  / +.32334215826050909646402930953354d-10    /
      data erfcs( 10)  / -.79910159470045487581607374708595d-12    /
      data erfcs( 11)  / +.17990725113961455611967245486634d-13    /
      data erfcs( 12)  / -.37186354878186926382316828209493d-15    /
      data erfcs( 13)  / +.71035990037142529711689908394666d-17    /
      data erfcs( 14)  / -.12612455119155225832495424853333d-18    /
      data erfcs( 15)  / +.20916406941769294369170500266666d-20    /
      data erfcs( 16)  / -.32539731029314072982364160000000d-22    /
      data erfcs( 17)  / +.47668672097976748332373333333333d-24    /
      data erfcs( 18)  / -.65980120782851343155199999999999d-26    /
      data erfcs( 19)  / +.86550114699637626197333333333333d-28    /
      data erfcs( 20)  / -.10788925177498064213333333333333d-29    /
      data erfcs( 21)  / +.12811883993017002666666666666666d-31    /
!
      data sqrtpi / 1.77245385090551602729816748334115d0 /
      data nterf, xbig, sqeps / 0, 2*0.d0 /

! external functions
      integer initds
      real (kind=8) d1mach, dcsevl, derfc
      external d1mach, dcsevl, derfc, initds

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     if (nterf.ne.0) go to 10
      nterf = initds (erfcs, 21, 0.1*sngl(d1mach(3)))
      xbig = sqrt (-log(sqrtpi*d1mach(3)))
      sqeps = sqrt (2.0d0*d1mach(3))
!
 10   y = dabs(x)
      if (y.gt.1.d0) go to 20
!
! erf(x) = 1.0 - erfc(x)  for  -1.0 .le. x .le. 1.0
!
      if (y.le.sqeps) derf = 2.0d0*x/sqrtpi
      if (y.gt.sqeps) &
        derf = x*(1.0d0 + dcsevl (2.d0*x*x-1.d0, erfcs, nterf))
      return
!
! erf(x) = 1.0 - erfc(x) for abs(x) .gt. 1.0
!
 20   if (y.le.xbig) derf = dsign (1.0d0-derfc(y), x)
      if (y.gt.xbig) derf = dsign (1.0d0, x)
!
      return
      end function derf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real (kind=8) function derfc (x)
! july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
      implicit none

! arguments
      real (kind=8) x

! local variables
      integer nterf, nterfc, nterc2
      real (kind=8) erfcs(21), erfccs(59), erc2cs(49), sqeps, sqrtpi
      real (kind=8) xmax, xsml, y
      real (kind=4) eta

!
! series for erf        on the interval  0.          to  1.00000e+00
!                                        with weighted error   1.28e-32
!                                         log weighted error  31.89
!                               significant figures required  31.05
!                                    decimal places required  32.55
!
      data erfcs(  1)  / -.49046121234691808039984544033376d-1     /
      data erfcs(  2)  / -.14226120510371364237824741899631d+0     /
      data erfcs(  3)  / +.10035582187599795575754676712933d-1     /
      data erfcs(  4)  / -.57687646997674847650827025509167d-3     /
      data erfcs(  5)  / +.27419931252196061034422160791471d-4     /
      data erfcs(  6)  / -.11043175507344507604135381295905d-5     /
      data erfcs(  7)  / +.38488755420345036949961311498174d-7     /
      data erfcs(  8)  / -.11808582533875466969631751801581d-8     /
      data erfcs(  9)  / +.32334215826050909646402930953354d-10    /
      data erfcs( 10)  / -.79910159470045487581607374708595d-12    /
      data erfcs( 11)  / +.17990725113961455611967245486634d-13    /
      data erfcs( 12)  / -.37186354878186926382316828209493d-15    /
      data erfcs( 13)  / +.71035990037142529711689908394666d-17    /
      data erfcs( 14)  / -.12612455119155225832495424853333d-18    /
      data erfcs( 15)  / +.20916406941769294369170500266666d-20    /
      data erfcs( 16)  / -.32539731029314072982364160000000d-22    /
      data erfcs( 17)  / +.47668672097976748332373333333333d-24    /
      data erfcs( 18)  / -.65980120782851343155199999999999d-26    /
      data erfcs( 19)  / +.86550114699637626197333333333333d-28    /
      data erfcs( 20)  / -.10788925177498064213333333333333d-29    /
      data erfcs( 21)  / +.12811883993017002666666666666666d-31    /
!
! series for erc2       on the interval  2.50000e-01 to  1.00000e+00
!                                        with weighted error   2.67e-32
!                                         log weighted error  31.57
!                               significant figures required  30.31
!                                    decimal places required  32.42
!
      data erc2cs(  1) / -.6960134660230950112739150826197d-1      /
      data erc2cs(  2) / -.4110133936262089348982212084666d-1      /
      data erc2cs(  3) / +.3914495866689626881561143705244d-2      /
      data erc2cs(  4) / -.4906395650548979161280935450774d-3      /
      data erc2cs(  5) / +.7157479001377036380760894141825d-4      /
      data erc2cs(  6) / -.1153071634131232833808232847912d-4      /
      data erc2cs(  7) / +.1994670590201997635052314867709d-5      /
      data erc2cs(  8) / -.3642666471599222873936118430711d-6      /
      data erc2cs(  9) / +.6944372610005012589931277214633d-7      /
      data erc2cs( 10) / -.1371220902104366019534605141210d-7      /
      data erc2cs( 11) / +.2788389661007137131963860348087d-8      /
      data erc2cs( 12) / -.5814164724331161551864791050316d-9      /
      data erc2cs( 13) / +.1238920491752753181180168817950d-9      /
      data erc2cs( 14) / -.2690639145306743432390424937889d-10     /
      data erc2cs( 15) / +.5942614350847910982444709683840d-11     /
      data erc2cs( 16) / -.1332386735758119579287754420570d-11     /
      data erc2cs( 17) / +.3028046806177132017173697243304d-12     /
      data erc2cs( 18) / -.6966648814941032588795867588954d-13     /
      data erc2cs( 19) / +.1620854541053922969812893227628d-13     /
      data erc2cs( 20) / -.3809934465250491999876913057729d-14     /
      data erc2cs( 21) / +.9040487815978831149368971012975d-15     /
      data erc2cs( 22) / -.2164006195089607347809812047003d-15     /
      data erc2cs( 23) / +.5222102233995854984607980244172d-16     /
      data erc2cs( 24) / -.1269729602364555336372415527780d-16     /
      data erc2cs( 25) / +.3109145504276197583836227412951d-17     /
      data erc2cs( 26) / -.7663762920320385524009566714811d-18     /
      data erc2cs( 27) / +.1900819251362745202536929733290d-18     /
      data erc2cs( 28) / -.4742207279069039545225655999965d-19     /
      data erc2cs( 29) / +.1189649200076528382880683078451d-19     /
      data erc2cs( 30) / -.3000035590325780256845271313066d-20     /
      data erc2cs( 31) / +.7602993453043246173019385277098d-21     /
      data erc2cs( 32) / -.1935909447606872881569811049130d-21     /
      data erc2cs( 33) / +.4951399124773337881000042386773d-22     /
      data erc2cs( 34) / -.1271807481336371879608621989888d-22     /
      data erc2cs( 35) / +.3280049600469513043315841652053d-23     /
      data erc2cs( 36) / -.8492320176822896568924792422399d-24     /
      data erc2cs( 37) / +.2206917892807560223519879987199d-24     /
      data erc2cs( 38) / -.5755617245696528498312819507199d-25     /
      data erc2cs( 39) / +.1506191533639234250354144051199d-25     /
      data erc2cs( 40) / -.3954502959018796953104285695999d-26     /
      data erc2cs( 41) / +.1041529704151500979984645051733d-26     /
      data erc2cs( 42) / -.2751487795278765079450178901333d-27     /
      data erc2cs( 43) / +.7290058205497557408997703680000d-28     /
      data erc2cs( 44) / -.1936939645915947804077501098666d-28     /
      data erc2cs( 45) / +.5160357112051487298370054826666d-29     /
      data erc2cs( 46) / -.1378419322193094099389644800000d-29     /
      data erc2cs( 47) / +.3691326793107069042251093333333d-30     /
      data erc2cs( 48) / -.9909389590624365420653226666666d-31     /
      data erc2cs( 49) / +.2666491705195388413323946666666d-31     /
!
! series for erfc       on the interval  0.          to  2.50000e-01
!                                        with weighted error   1.53e-31
!                                         log weighted error  30.82
!                               significant figures required  29.47
!                                    decimal places required  31.70
!
      data erfccs(  1) / +.715179310202924774503697709496d-1        /
      data erfccs(  2) / -.265324343376067157558893386681d-1        /
      data erfccs(  3) / +.171115397792085588332699194606d-2        /
      data erfccs(  4) / -.163751663458517884163746404749d-3        /
      data erfccs(  5) / +.198712935005520364995974806758d-4        /
      data erfccs(  6) / -.284371241276655508750175183152d-5        /
      data erfccs(  7) / +.460616130896313036969379968464d-6        /
      data erfccs(  8) / -.822775302587920842057766536366d-7        /
      data erfccs(  9) / +.159214187277090112989358340826d-7        /
      data erfccs( 10) / -.329507136225284321486631665072d-8        /
      data erfccs( 11) / +.722343976040055546581261153890d-9        /
      data erfccs( 12) / -.166485581339872959344695966886d-9        /
      data erfccs( 13) / +.401039258823766482077671768814d-10       /
      data erfccs( 14) / -.100481621442573113272170176283d-10       /
      data erfccs( 15) / +.260827591330033380859341009439d-11       /
      data erfccs( 16) / -.699111056040402486557697812476d-12       /
      data erfccs( 17) / +.192949233326170708624205749803d-12       /
      data erfccs( 18) / -.547013118875433106490125085271d-13       /
      data erfccs( 19) / +.158966330976269744839084032762d-13       /
      data erfccs( 20) / -.472689398019755483920369584290d-14       /
      data erfccs( 21) / +.143587337678498478672873997840d-14       /
      data erfccs( 22) / -.444951056181735839417250062829d-15       /
      data erfccs( 23) / +.140481088476823343737305537466d-15       /
      data erfccs( 24) / -.451381838776421089625963281623d-16       /
      data erfccs( 25) / +.147452154104513307787018713262d-16       /
      data erfccs( 26) / -.489262140694577615436841552532d-17       /
      data erfccs( 27) / +.164761214141064673895301522827d-17       /
      data erfccs( 28) / -.562681717632940809299928521323d-18       /
      data erfccs( 29) / +.194744338223207851429197867821d-18       /
      data erfccs( 30) / -.682630564294842072956664144723d-19       /
      data erfccs( 31) / +.242198888729864924018301125438d-19       /
      data erfccs( 32) / -.869341413350307042563800861857d-20       /
      data erfccs( 33) / +.315518034622808557122363401262d-20       /
      data erfccs( 34) / -.115737232404960874261239486742d-20       /
      data erfccs( 35) / +.428894716160565394623737097442d-21       /
      data erfccs( 36) / -.160503074205761685005737770964d-21       /
      data erfccs( 37) / +.606329875745380264495069923027d-22       /
      data erfccs( 38) / -.231140425169795849098840801367d-22       /
      data erfccs( 39) / +.888877854066188552554702955697d-23       /
      data erfccs( 40) / -.344726057665137652230718495566d-23       /
      data erfccs( 41) / +.134786546020696506827582774181d-23       /
      data erfccs( 42) / -.531179407112502173645873201807d-24       /
      data erfccs( 43) / +.210934105861978316828954734537d-24       /
      data erfccs( 44) / -.843836558792378911598133256738d-25       /
      data erfccs( 45) / +.339998252494520890627359576337d-25       /
      data erfccs( 46) / -.137945238807324209002238377110d-25       /
      data erfccs( 47) / +.563449031183325261513392634811d-26       /
      data erfccs( 48) / -.231649043447706544823427752700d-26       /
      data erfccs( 49) / +.958446284460181015263158381226d-27       /
      data erfccs( 50) / -.399072288033010972624224850193d-27       /
      data erfccs( 51) / +.167212922594447736017228709669d-27       /
      data erfccs( 52) / -.704599152276601385638803782587d-28       /
      data erfccs( 53) / +.297976840286420635412357989444d-28       /
      data erfccs( 54) / -.126252246646061929722422632994d-28       /
      data erfccs( 55) / +.539543870454248793985299653154d-29       /
      data erfccs( 56) / -.238099288253145918675346190062d-29       /
      data erfccs( 57) / +.109905283010276157359726683750d-29       /
      data erfccs( 58) / -.486771374164496572732518677435d-30       /
      data erfccs( 59) / +.152587726411035756763200828211d-30       /
!
      data sqrtpi / 1.77245385090551602729816748334115d0 /
      data nterf, nterfc, nterc2, xsml, xmax, sqeps / 3*0, 3*0.d0 /

! external functions
      integer initds
      real (kind=8) d1mach, dcsevl
      external d1mach, dcsevl, initds

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (nterf.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      nterf = initds (erfcs, 21, eta)
      nterfc = initds (erfccs, 59, eta)
      nterc2 = initds (erc2cs, 49, eta)
!
      xsml = -sqrt (-log(sqrtpi*d1mach(3)))
      xmax = sqrt (-log(sqrtpi*d1mach(1)) )
      xmax = xmax - 0.5d0*log(xmax)/xmax - 0.01d0
      sqeps = sqrt (2.0d0*d1mach(3))
!
 10   if (x.gt.xsml) go to 20
!
! erfc(x) = 1.0 - erf(x)  for  x .lt. xsml
!
      derfc = 2.0d0
      return
!
 20   if (x.gt.xmax) go to 40
      y = dabs(x)
      if (y.gt.1.0d0) go to 30
!
! erfc(x) = 1.0 - erf(x)  for abs(x) .le. 1.0
!
      if (y.lt.sqeps) derfc = 1.0d0 - 2.0d0*x/sqrtpi
      if (y.ge.sqeps) then
        derfc = 1.0d0 - x*(1.0d0 + dcsevl (2.d0*x*x-1.d0, erfcs, nterf))
      endif
      return
!
! erfc(x) = 1.0 - erf(x)  for  1.0 .lt. abs(x) .le. xmax
!
 30   y = y*y
      if (y.le.4.d0) derfc = exp(-y)/dabs(x) * (0.5d0 + dcsevl ( &
        (8.d0/y-5.d0)/3.d0, erc2cs, nterc2) )
      if (y.gt.4.d0) derfc = exp(-y)/dabs(x) * (0.5d0 + dcsevl ( &
        8.d0/y-1.d0, erfccs, nterfc) )
      if (x.lt.0.d0) derfc = 2.0d0 - derfc
      return
!
 40   write(0,*) "derfc   x so big erfc underflows"
      derfc = 0.d0
      return
!
      end function derfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function initds (dos, nos, eta)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
! initialize the real (kind=8) orthogonal series dos so that initds
! is the number of terms needed to insure the error is no larger than
! eta.  ordinarily eta will be chosen to be one-tenth machine precision.
!
!             input arguments --
! dos    dble prec array of nos coefficients in an orthogonal series.
! nos    number of coefficients in dos.
! eta    requested accuracy of series.
!
      implicit none

! arguments
      integer nos
      real (kind=8) dos(nos)
      real (kind=4) eta

! local variables
      integer i, ii
      real (kind=4) err

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (nos.lt.1) then
        write(0,*) "initds  number of coefficients lt 1"
        stop
      endif
!
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
!
 20   if (i.eq.nos) then
        write(0,*) "initds  eta may be too small"
        stop
      endif

      initds = i
!
      return
      end function initds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real (kind=8) function dcsevl (x, a, n)
!
! evaluate the n-term chebyshev series a at x.  adapted from
! r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
!
!             input arguments --
! x      dble prec value at which the series is to be evaluated.
! a      dble prec array of n terms of a chebyshev series.  in eval-
!        uating a, only half the first coef is summed.
! n      number of terms in array a.
!
      implicit none

! arguments
      integer n
      real (kind=8) x, a(n)

! local variables
      integer i, ni
      real (kind=8) twox, b0, b1, b2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (n.lt.1) then
        write(0,*) "dcsevl  number of terms le 0"
        stop
      endif
      if (n.gt.1000) then
        write(0,*) "dcsevl  number of terms gt 1000"
        stop
      endif
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) then
        write(0,*) "dcsevl  x outside (-1,+1)"
      endif
!
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
!
      dcsevl = 0.5d0 * (b0-b2)
!
      return
      end function dcsevl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a portable version of D1MACH that should work for any compiler
! that complies with the Fortran 90 specifications for the intrinsic
! EPSILON.

      real (kind=8) function d1mach(i) 
!                                                                       
!  double-precision machine constants                                   
!                                                                       
!  d1mach( 1) = b**(emin-1), the smallest positive magnitude.           
!                                                                       
!  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.           
!                                                                       
!  d1mach( 3) = b**(-t), the smallest relative spacing.                 
!                                                                       
!  d1mach( 4) = b**(1-t), the largest relative spacing.                 
!                                                                       
!  d1mach( 5) = log10(b)                                                
!                                                                       
      implicit none

! arguments
      integer i

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      select case (i) 
        case (1) 
           d1mach = tiny(1d0) 
        case (2) 
           d1mach = huge(1d0) 
        case (3) 
           d1mach = epsilon(1d0)/radix(1d0) 
        case (4) 
           d1mach = epsilon(1d0) 
        case (5) 
           d1mach = radix(1d0) 
           d1mach = log10(d1mach) 
        case default 
        write(0,'(a,i10)') ' d1mach - i out of bounds', i 
        stop 
      end select 
      end function d1mach
