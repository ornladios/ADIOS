SUBROUTINE sctednncs( n, j, k, rho, t, ye, &
&  a0w, b0w, c0w, a1w, b1w, c1w, rmdnncs0, rmdnncs1, rmdnncs, &
&  arnncs, brnncs, crnncs, aenncs, benncs, cenncs, arscte, brscte, &
&  crscte, arsctu, brsctu, crsctu, arscti, brscti, crscti, &
&  arsctd, brsctd, crsctd, aesctu, besctu, cesctu, aescti, &
&  bescti, cescti, aesctd, besctd, cesctd )
!-----------------------------------------------------------------------
!
!    File:         sctarate
!    Module:       sctarate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To evaluate the neutrino-nucleon inelastic scattering rates.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!      kind_module, numerical_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: j             ! radial zone index
INTEGER, INTENT(in)               :: k             ! neutrino energy index

REAL(KIND=double), INTENT(in)     :: rho           ! density (g cm^{-3})
REAL(KIND=double), INTENT(in)     :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)     :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)    :: a0w
REAL(KIND=double), INTENT(out)    :: b0w
REAL(KIND=double), INTENT(out)    :: c0w
REAL(KIND=double), INTENT(out)    :: a1w
REAL(KIND=double), INTENT(out)    :: b1w
REAL(KIND=double), INTENT(out)    :: c1w

REAL(KIND=double), INTENT(out)    :: rmdnncs0
REAL(KIND=double), INTENT(out)    :: rmdnncs1
REAL(KIND=double), INTENT(out)    :: rmdnncs

REAL(KIND=double), INTENT(out)    :: arnncs
REAL(KIND=double), INTENT(out)    :: brnncs
REAL(KIND=double), INTENT(out)    :: crnncs
REAL(KIND=double), INTENT(out)    :: aenncs
REAL(KIND=double), INTENT(out)    :: benncs
REAL(KIND=double), INTENT(out)    :: cenncs

REAL(KIND=double), INTENT(out)    :: arscte
REAL(KIND=double), INTENT(out)    :: brscte
REAL(KIND=double), INTENT(out)    :: crscte

REAL(KIND=double), INTENT(out)    :: arsctu
REAL(KIND=double), INTENT(out)    :: brsctu
REAL(KIND=double), INTENT(out)    :: crsctu

REAL(KIND=double), INTENT(out)    :: arscti
REAL(KIND=double), INTENT(out)    :: brscti
REAL(KIND=double), INTENT(out)    :: crscti

REAL(KIND=double), INTENT(out)    :: arsctd
REAL(KIND=double), INTENT(out)    :: brsctd
REAL(KIND=double), INTENT(out)    :: crsctd

REAL(KIND=double), INTENT(out)    :: aesctu
REAL(KIND=double), INTENT(out)    :: besctu
REAL(KIND=double), INTENT(out)    :: cesctu

REAL(KIND=double), INTENT(out)    :: aescti
REAL(KIND=double), INTENT(out)    :: bescti
REAL(KIND=double), INTENT(out)    :: cescti

REAL(KIND=double), INTENT(out)    :: aesctd
REAL(KIND=double), INTENT(out)    :: besctd
REAL(KIND=double), INTENT(out)    :: cesctd

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

a0w            = zero
b0w            = zero
c0w            = zero
a1w            = zero
b1w            = zero
c1w            = zero
rmdnncs0       = zero
rmdnncs1       = zero
rmdnncs        = zero
arnncs         = zero
brnncs         = zero
crnncs         = zero
aenncs         = zero
benncs         = zero
cenncs         = zero
arscte         = zero
brscte         = zero
crscte         = zero
arsctu         = zero
brsctu         = zero
crsctu         = zero
arscti         = zero
brscti         = zero
crscti         = zero
arsctd         = zero
brsctd         = zero
crsctd         = zero
aesctu         = zero
besctu         = zero
cesctu         = zero
aescti         = zero
bescti         = zero
cescti         = zero
aesctd         = zero
besctd         = zero
cesctd         = zero

RETURN
END SUBROUTINE sctednncs
