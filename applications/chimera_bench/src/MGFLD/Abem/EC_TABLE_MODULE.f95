!-----------------------------------------------------------------------
!    Module:       ec_table_module
!    Author:       W. R. Hix
!    Date:         7/20/07
!
!    Specification of the tabular electron capture table grid and data
!-----------------------------------------------------------------------

module ec_table_module

Use kind_module

integer, parameter :: npts=200, ntemp=37, nrho=11, nye=8 ! Revised Table
! integer, parameter :: npts=200, ntemp=35, nrho=11, nye=7  ! Original Table
integer, parameter :: tmx=40, tmn=(tmx+1-ntemp), rhomx=26, rhomn=(rhomx+1-nrho)
real(KIND=double), parameter :: deltaE=0.5,deltalrho_ec=.5,deltaT_ec=.1,deltaYe_ec=.02
real(KIND=double) :: energy(0:npts)
real(KIND=double) :: jectab(0:npts,rhomn:rhomx,tmn:tmx,0:(nye-1))
integer :: yefloor(rhomn:rhomx)

end module ec_table_module

