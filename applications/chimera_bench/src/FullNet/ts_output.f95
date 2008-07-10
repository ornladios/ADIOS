SUBROUTINE ts_output( kstep, enuc, edot, kts, knr ) 
!-----------------------------------------------------------------------
!  If the flag itso is > 0, full_net calls this routine to handle stepwise 
!  output.  
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE abundances, ONLY : y
USE nuclear_data, ONLY : aa
USE controls, ONLY : itso
USE conditions, ONLY : t, tdel, t9t, rhot
!     use match_data
!     use flux_data

IMPLICIT none
SAVE

INTEGER           :: ii(7), kstep, kout, kts, knr
REAL(KIND=double) :: xg(7), enuc, edot
!-----
!  An abundance snapshot is written to the binary file on unit 24.
!     Write(24) kstep,t,t9t,rhot,tdel,(y(k),k=1,nnet)
!-----
!  For itso>=2, output important mass fractions to the ASCII file on unit 22
IF ( itso >= 2 ) THEN
  kout        = 10 * kts + knr
!         xg=0.0
!         Do i=1,nnet
!             If(zz(i)<=1) Then
!                 ig=1
!             Elseif(zz(i)<=2) Then
!                 ig=2
!             Elseif(zz(i)<=8) Then
!                 ig=3
!             Elseif(zz(i)<=10) Then
!                 ig=4
!             Elseif(zz(i)<=12) Then
!                 ig=5
!             Elseif(zz(i)<=14) Then
!                 ig=6
!             Else
!                 ig=7
!             Endif
!             xg(ig)=xg(ig)+aa(i)*y(i)
!         Enddo
  ii          = (/1,2,6,11,12,13,14/)
  xg          = aa(ii)*y(ii)
  xg(7)       = SUM(aa*y) - 1.0d0
!         Write(50,"(i6,1es15.8,2es10.3,2es10.2,8es9.2,i4)") &
!    &        kstep,t,t9t,rhot,edot,enuc,tdel,(xg(i),i=1,7),kout
!         If(kstep>0) Then
!             call flux
!         Else
!            flx_int=0.0
!         Endif
!-----
!  For itso>=3, output time and thermo evolution to the screen
  IF ( itso >= 3 ) WRITE (6,"(i5,4es12.4,2i3)") &
&        kstep,t,tdel,t9t,rhot,knr,kts
END IF

RETURN
END SUBROUTINE ts_output
