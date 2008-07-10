SUBROUTINE abemtnc( j, n, rho, t, ye, xh, ah, cmpn, cmpp, cmpe, absrnc, emitnc)
!-----------------------------------------------------------------------
!
!    File:         abemtnc
!    Module:       abemtnc
!    Type:         Subprogram
!    Author:       W.R. Hix, Physics Division
!                  Oak Ridge National Laboratory, Oak Ridge TN 37831
!
!    Date:         7/17/07
!
!    Purpose:
!      To compute the inverse mean free paths for the absorption and emission
!       of n-type neutrinos using the NSE folded table of electron capture 
!       rates.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  interp_ec_table
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
!  nu_energy_grid_module, nu_dist_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, one
USE physcnst_module, ONLY : rmu, kmev, dmnp

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx, unui
USE nu_dist_module, ONLY : unu
USE prb_cntl_module, ONLY : iaenct, roaenct

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j               ! radial zone index
INTEGER, INTENT(in)              :: n               ! neutrino flavor index

REAL(KIND=double), INTENT(in)    :: rho             ! density (g/cm^3)
REAL(KIND=double), INTENT(in)    :: t               ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye              ! electron fraction
REAL(KIND=double), INTENT(in)    :: xh              ! heavy nuclei mass fraction
REAL(KIND=double), INTENT(in)    :: ah              ! heavy nuclei mass number
REAL(KIND=double), INTENT(in)    :: cmpn            ! neutron chemical porential
REAL(KIND=double), INTENT(in)    :: cmpp            ! proton chemical porential
REAL(KIND=double), INTENT(in)    :: cmpe            ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout) :: absrnc(nnugpmx) ! inverse mean free path for absorption on free nucleons
REAL(KIND=double), INTENT(inout) :: emitnc(nnugpmx) ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: je              ! do index

REAL(KIND=double)                :: tmev            ! temperature (MeV)
REAL(KIND=double)                :: xnuc            ! number density of nuclei (cm^{-3})
REAL(KIND=double)                :: fexp            ! exponential function

EXTERNAL fexp

!-----------------------------------------------------------------------
!  Set emitnc and absrnc to zero and return if
!               iaenct = 0
!  or
!               rho > roaenct
!  or
!               n /= 1
!-----------------------------------------------------------------------

IF ( iaenct == 0  .or.  rho  >  roaenct  .or.  n /= 1  .or.  ah < 40 ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

emitnc            = zero
absrnc            = zero
tmev              = kmev * t
xnuc              = ( xh/( rmu * ah ) ) * rho

!-----------------------------------------------------------------------
!
!       \\\\\ E-NEUTRINO-NUCLEUS  EMISSION  AND ABSORPTION/////
!
!-----------------------------------------------------------------------

IF ( n == 1 ) THEN

!-----------------------------------------------------------------------
!  Calculate values for rate/per heavy nucleus by interpolation
!-----------------------------------------------------------------------

  CALL interp_ec( j, emitnc, t, rho, ye ) 

!-----------------------------------------------------------------------
!  Multiply by heavy nucleus abundance
!-----------------------------------------------------------------------

  emitnc = emitnc * xnuc

!-----------------------------------------------------------------------
!  Compute the absorption rate using detailed balance
!-----------------------------------------------------------------------
  
!  absrnc = emitnc * fexpv((unu(j,:) + dmnp + cmpn - cmpp - cmpe)/tmev)
!  Forall (je=1:nnugpmx) absrnc(je) = emitnc(je) * fexp((unu(j,je) + dmnp + cmpn - cmpp - cmpe)/tmev)
   DO je = 1,nnugpmx 
     absrnc(je) = emitnc(je) * fexp((unu(j,je) + dmnp + cmpn - cmpp - cmpe)/tmev)
   END DO
END IF ! ( n == 1 )


RETURN
END SUBROUTINE abemtnc
