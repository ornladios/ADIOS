SUBROUTINE interp_ec( j, jec, temp, rho, Ye )  
!-----------------------------------------------------------------------
!
!    File:         interp_ec
!    Module:       interp_ec
!    Type:         subprogram
!    Author:       W.R. Hix, Physics Division
!                  Oak Ridge National Laboratory, Oak Ridge TN 37831
!                  E.J. Lentz, Department of Physics and Astronomy
!                  University of Tennessee, Knoxville TN 37919
!
!    Date:         7/17/07
!
!    Purpose:
!      To compute the inverse mean free paths for the absorption and emission
!       of n-type neutrinos using the NSE folded table of electron capture 
!       rates.
!
!    Variables that must be passed through common:
!        spec, energy, jectab
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  rho         : matter density (g/cm**3)
!  t           : matter temperature (K)
!  ye          : electron fraction
!  xh          : heavy nucleus mass fraction
!  ah          : heavy nucleus mass number
!  cmpn        : free neutron chemical potential (excluding rest mass) (MeV)
!  cmpp        : free proton chemical potential (excluding rest mass) (MeV)
!  cmpe        : electron chemical potential (including rest mass) (MeV)
!
!    Output arguments:
!  absrnc      :  absorption inverse mean free path (/cm)
!  emitnc      :  emission inverse mean free path (/cm)
!
!    Input arguments (common):
!
!  iaenct     : 0; inverse mean free paths set to zero.
!               1; inverse mean free paths computed.
!  roaenct    : density above which rates are set to zero.
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module,
!  nu_energy_grid_module, prb_cntl_module, ec_table_module
!
!-----------------------------------------------------------------------

Use kind_module
Use physcnst_module, ONLY: kmev
Use nu_energy_grid_module, ONLY : nnugp, nnugpmx, unui, dunui, unubi
Use nu_dist_module, ONLY : unu, unub, dunu
Use ec_table_module
use edit_module, ONLY : nprint, nlog

implicit none

INTEGER, INTENT(in)              :: j             ! radial zone index
Real(KIND=double), dimension(nnugpmx), intent(out) :: jec
Real(KIND=double) :: temp, rho, ye
Real(KIND=double), dimension(0:npts) :: jecfine
Integer :: kfmin(nnugpmx),kfmax(nnugpmx)
Real(KIND=double) :: deltaupper(nnugpmx),deltalower(nnugpmx)
integer :: i,k,l,n,m,ii,kk
Integer :: rhoup,rhodown,yeup,yedown,tup,tdown
Real(KIND=double) :: c1,c2,c3,c3base,loctot,tmev
Integer :: irho,iye,itemp
Integer :: kmax,ktop,yeoffset(2)
Logical :: askew    

kmax=nnugp(1)
            
!-----------------------------------------------------------------------
!  Find temperature in table.  Temperature tabulated in MeV.
!  For temperature outside of table, set rates to small value
!-----------------------------------------------------------------------
tmev=temp*kmev
itemp=int(tmev/deltaT_ec)
tdown=itemp
tup=itemp+1
    
If (tup > tmx .or. tdown < tmn ) Then
!  Write(nprint,'(a,3es10.3,i4)') "temperature out of range in ec table ",rho,temp,ye,itemp
  jec=10.0d-203
  return
Endif

!-----------------------------------------------------------------------
!  Find density in table.
!  For density outside of table, set rates to small value
!-----------------------------------------------------------------------
irho=int(log10(rho)/deltalrho_ec)
rhodown=irho
rhoup = irho+1

If (rhoup > rhomx .or. rhodown < rhomn) Then
!  Write(nprint,'(a,3es10.3,i4)') "Density out of range in ec table ",rho,temp,ye,irho
  jec=10.0d-203
  return
Endif

!-----------------------------------------------------------------------
!  Find ye in table.  Ye grid is density dependent
!  For ye outside of table, allow extrapolation
!-----------------------------------------------------------------------

iye=int(Ye/deltaYe_ec)
yeoffset(1)=iye-yefloor(rhodown)
yeoffset(2)=iye-yefloor(rhoup)

If (sum(yeoffset) < -2) Then
!  Write(nprint,'(a,3es10.3,3i4)') "Ye out of range in ec table ",rho,temp,ye,iye,yeoffset
! jec=10.0d-203
! Return
Endif
 
If (sum(yeoffset) > (2*nye-2)) Then
!  Write(nprint,'(a,2es10.3,3i4)') "Ye out of range in ec table:",rho,ye,iye,yeoffset
! jec=10.0d-203
! Return
Endif

!-----------------------------------------------------------------------
!  Because of density dependence, Ye may require askew interpolation.
!-----------------------------------------------------------------------

askew=.false.
Do ii=1,2
  If (yeoffset(ii).lt.0) Then
    askew=.true.
    yeoffset(ii) =0
  Endif
Enddo

Do ii=1,2
  If (yeoffset(ii).gt.(nye-2)) Then
    askew=.true.
    yeoffset(ii) =nye-2
  Endif
Enddo

!-----------------------------------------------------------------------
!  Compute interpolation forefactors.
!-----------------------------------------------------------------------
c1=(log10(rho)-dble(rhodown)*deltalrho_ec)
!  /deltalrho_ec
 
c2=(log10(tmev/(dble(itemp)*deltaT_ec)))
! c2=(log10(tmev/(dble(itemp)*deltaT_ec)))
! /(log10(dble(tup)/dble(tdown)))  !a deltaT_ec factor is ommitted from numerator and denominator
 
c3base=(ye-dble(iye)*deltaYe_ec)/deltaYe_ec
c3=c3base + &
& c1*dble(yefloor(rhoup)+yeoffset(2)-yefloor(rhodown)-yeoffset(1)) !correction for skew If present
!Write(nprint,*) "c3 ecinterp:",i,askew,c3base,c3
  
!-----------------------------------------------------------------------
! Interpolate rate on table energy grid
!-----------------------------------------------------------------------

jecfine(:) = c3 & !note on indices: yeoffset=1 for rhodown, 2 for rhoup
&           *( (1.0-c1)*(1.0-c2) * jectab(:,rhodown,tdown,yeoffset(1)+1) &  !   add one in upper block for old yeup
&           +  c1*(1.0-c2)       * jectab(:,rhoup,tdown,yeoffset(2)+1) &
&           +  c2*(1.0-c1)       * jectab(:,rhodown,tup,yeoffset(1)+1) &
&           +  c1*c2             * jectab(:,rhoup,tup,yeoffset(2)+1) ) &
&           + (1.0-c3) &
&           *( (1.0-c1)*(1.0-c2) * jectab(:,rhodown,tdown,yeoffset(1)) &
&           +  c1*(1.0-c2)       * jectab(:,rhoup,tdown,yeoffset(2)) &
&           +  c2*(1.0-c1)       * jectab(:,rhodown,tup,yeoffset(1)) &
&           +  c1*c2             * jectab(:,rhoup,tup,yeoffset(2)) )

!Write(nprint,*) "ecinterp: now call exponential", maxval(jecfine)
 
jecfine=10.0**(jecfine)  

!-----------------------------------------------------------------------
! Match Table grid to MGFLD energy grid
!-----------------------------------------------------------------------

ktop=kmax ! top energy bin that has a EC contribution
Do k=1,kmax
  kfmin(k)=int(unub(j,k)/deltaE)+1
  kfmax(k)=int(unub(j,k+1)/deltaE)
  If (kfmax(k).gt.npts) kfmax(k)=npts
  If (kfmin(k).gt.npts) kfmin(k)=npts
  If (unub(j,k+1) .gt. energy(npts)) ktop=min(k,ktop)
!  Write(nprint,*) "Enu:", unub(j,k), energy(kfmin(k)), energy(kfmax(k))
Enddo
!Write(nprint,*) "ktop:",ktop, unub(j,ktop),unub(j,ktop+1)

!-----------------------------------------------------------------------
! Calculate integration constants for Table -> MGFLD energy grid
!-----------------------------------------------------------------------

Do k=1,ktop
  deltalower(k)=(energy(kfmin(k))-unub(j,k))/deltaE ! assuming no "undershoot" problems on the bottom.
  deltaupper(k)=(unub(j,k+1)-energy(kfmax(k)))/deltaE
  If (deltalower(k) < 0.d0) Then
!    Write (nprint,*) "deltalower < 0.",k,deltalower(k)
    deltalower(k)=0.d0
  Endif
  If (deltaupper(k) < 0.d0) Then
!    Write (nprint,*) "deltaupper < 0.",k,deltaupper(k)
    deltaupper(k)=0.d0
  Endif
  If (deltalower(k) > 1.d0) Then
!    Write (nprint,*) "deltalower > 1.",k,deltalower(k)
    deltalower(k)=1.d0
  Endif
  If (deltaupper(k) > 1.d0) Then
!    Write (nprint,*) "deltaupper > 1.",k,deltaupper(k)
    deltaupper(k)=1.d0
  Endif
Enddo
If(unub(j,kmax+1).gt. energy(npts)) deltaupper(ktop)=0.d0 

!-----------------------------------------------------------------------
! Integrate onto MGFLD energy grid
!-----------------------------------------------------------------------

Do k=1,ktop

  If (kfmax(k) .lt. kfmin(k)) Then ! If MGFLD bin fits within a single table bin
    loctot = 0.5d0*dunu(j,k) &
&          * ( (deltaupper(k)+1.d0-deltalower(k))*jecfine(kfmin(k)) &
&          + (1.d0-deltaupper(k)+deltalower(k))*jecfine(kfmax(k)) )
  Else 
    loctot = deltalower(k)*deltaE*0.5d0 & 
&          * ( deltalower(k)*jecfine(kfmin(k)-1) &
&          + (2.0-deltalower(k))*jecfine(kfmin(k)) )
    Do kk=kfmin(k),kfmax(k)-1
       loctot = loctot+0.5*deltaE*(jecfine(kk)+jecfine(kk+1))
    Enddo
    If(k<ktop) Then
      loctot = loctot+deltaupper(k)*deltaE*0.5 &
&            * ( (deltaupper(k))*jecfine(kfmax(k)+1) &
&            + (2.d0-deltaupper(k))*jecfine(kfmax(k)) )
     Endif
   Endif
   jec(k)=loctot/(dunu(j,k)*unu(j,k)**2)

Enddo


!.... check completeness of chosen input

loctot=0.d0
do k=1,kmax
  loctot=loctot+dunu(j,k)*unu(j,k)**2*jecfine(k)
enddo
!write(nprint,*) rho,ye,temp,loctot,irho,iye,itemp
 
Return 
End 
