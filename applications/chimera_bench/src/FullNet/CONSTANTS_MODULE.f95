MODULE constants
!===============================================================================
!  These are fundamental constants in CGS, except energies in MeV, Temp in GK
!===============================================================================

USE kind_module, ONLY : double

REAL(KIND=double), PARAMETER :: pi   =3.1415926536, hbar=6.582122d-22
REAL(KIND=double), PARAMETER :: amu  =1.036427d-18, bok =8.617385d-02
REAL(KIND=double), PARAMETER :: avn  =6.022137d+23, e2  =1.439830d-13
REAL(KIND=double), PARAMETER :: epmev=1.602177d-06

END MODULE constants
