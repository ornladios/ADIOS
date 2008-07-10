Subroutine read_ec_table
!-----------------------------------------------------------------------
!
!    File:         Read_ec_table
!    Module:       Read_ec_table
!    Type:         subprogram
!    Author:       W. R. Hix, ORNL
!                  E. J. Lentz, UTK
!
!    Date:         7/19/07 
!
!    Purpose:
!      To read the tabular electron capture grid and data.  
!
!    Variables that must be passed through common:
!  nnugp(n)            : number of energy zones for neutrinos of type n
!  unubi(k)            : energy of energy zone k (MeV)
!  dunu(k)             : energy width of energy zone k (MeV)
!
!    Subprograms called:
!      none
!
!    Input arguments:
!      none
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!
!  iaenct              : 0, e-neutrino emission from nuclei
!                          (tabular rates) omitted
!                      : 1, e-neutrino emission from nuclei 
!                          (tabular rates) included
!  ntemp, nrho,nye     : number of grid points in respective directions.
!  tmn,tmx             : minimum and maximum temperature grid points
!  rhomn,rhomx         : minimum and maximum density grid points
!
!    Output arguments (common):
!  jectab              : nuclear electron capture emissivity table
!  energy              : table energy grid
!
!    Modules used:
!      physcnst_module
!      ec_table_module, edit_module, prb_cntl_module
!
!-----------------------------------------------------------------------

use physcnst_module, ONLY : cvel,pi,hbar

use ec_table_module
Use edit_module, ONLY : nlog, nprint
USE prb_cntl_module, ONLY : iaenct, roaenct

implicit none

real(KIND=double) :: spec(0:npts)
Integer, parameter                 :: nrdata = 31   ! unit number for accessing the data files
Integer                            :: istat         ! open-close file flag
integer :: i,j,k,kk,jj
integer :: itemp,irho,iye,it,iy,ir
character(80) :: tablefile 
character (LEN=120) :: desc(2)
integer :: izone,nline,npart
real(KIND=double) :: rho,ye,temp,rate,dum
logical :: ispart
integer, parameter :: nperline=7, outperline=6
real(KIND=double), parameter :: gigK=11.60451 !MeV
real(KIND=double), parameter :: eps=0.01 ! 1% tolerance on grid readin
real(KIND=double), parameter :: const= 2.0*(pi*cvel)**2*hbar**3

!-----------------------------------------------------------------------
! Return if tabular ec rates are not used.
!-----------------------------------------------------------------------
If (iaenct /=1 ) Then
  Return
Endif

!-----------------------------------------------------------------------
!  Determine spectrum input file geometry
!-----------------------------------------------------------------------

nline=(npts+1)/nperline
npart=(npts+1)-nline*nperline
ispart= npart.gt.0

!-----------------------------------------------------------------------
!  Initialize table energy grid data 
!-----------------------------------------------------------------------
Do k=0,npts
  energy(k)=dble(k)*deltaE
Enddo
jectab=0.d0

!-----------------------------------------------------------------------
!  Open file
!-----------------------------------------------------------------------

tablefile = "../../MGFLD/Abem/ec_table.d"
Open (nrdata,file=trim(tablefile),status='old',iostat=istat)
If ( istat /= 0 ) Then
  Write (nprint,"(a,a,a)") 'EC_table data file ',trim(tablefile),' not found.'
  STOP
Endif
     
!-----------------------------------------------------------------------
!  Read table data
!-----------------------------------------------------------------------
Read(nrdata,"(a)") desc(1)
Read(nrdata,"(a)") desc(2)

Do it=tmn,tmx ; Do ir=rhomn,rhomx 
  Do iy=(nye-1),0,-1  
    Read(nrdata,"(6x,i2,6x,f12.6,5x,f9.6,5x,f12.6,5x,f9.6,5x,f12.6)") & 
&     izone,rho,ye,temp,dum,rate
    irho=int(log10(rho)/deltalrho_ec+eps)
    itemp=int(temp/gigK/deltaT_ec + eps)
    iye=int(ye/deltaYe_ec + eps)
    If (irho /= ir) Then
      Write(nprint,*) 'EC table misread: density ',irho, ' /= ',ir
      Stop
    Elseif (itemp /= it) Then
      Write(nprint,*) 'EC table misread: temperature ',itemp, ' /= ',it
      Stop
    Endif
    Read(nrdata,*) 
    Do j=1,nline
     Read(nrdata,*) (spec(k),k=(j-1)*nperline,j*nperline-1)
!    Write(nprint,'(7es13.6)') (spec(k),k=(j-1)*nperline,j*nperline-1)
    Enddo
    If (ispart) Read(nrdata,*) (spec(k),k=nline,nperline,npts)
!   If (ispart) Write(6,'(7es13.6)') (spec(k),k=nline,nperline,npts)
    jectab(:,ir,it,iy)=rate*spec
  Enddo
  If(it==tmx) Then
   yefloor(ir)=iye
!   Write(nprint,'(a,es10.3,a,f5.3)') "yefloor for rho=",rho,': ',.02*yefloor(ir)
  Endif
Enddo ; Enddo

Close(nrdata)

!-----------------------------------------------------------------------
!  Multiply spectrum by rate and take log of table
!-----------------------------------------------------------------------

Where (jectab <= 0.d0)
  jectab = -200.d0
Elsewhere
  jectab = log10(const*jectab)
Endwhere

!-----------------------------------------------------------------------
!  Signal completion 
!-----------------------------------------------------------------------
Write(nlog,'(a,a)') ' Finished reading ',trim(desc(1))
Write(nlog,'(a,a)') ' Table grid:',trim(desc(2))

Return

End Subroutine read_ec_table
