!-----------------------------------------------------------------------
!    MODULE:       el_eos_modlue
!    Author:       S. W. Bruenn
!    Date:         10/26/02
!-----------------------------------------------------------------------

MODULE el_eos_module

USE kind_module
SAVE

!-----------------------------------------------------------------------
!  NSUBE, QSUBE - Electron density and concentration
!
!  NEPLUS - Positron density (Not set in our simple ultrarelativistic EOS)
!
!  ACOEF, BCOEF, CCOEF, DBDT - Temporary coefficients
!
!  MUSUBE - Electron chemical potential
!
!  EPRESS, EU, FSUBE, ES - Electron pressure, internal energy, Helmholtz
!   free energy, and entropy
!
!  PPRESS, PS, PU, PF - Photon pressure, internal energy, Helmholtz free
!   energy, and entropy
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: NSUBE, QSUBE
REAL(KIND=double)                               :: NEPLUS
REAL(KIND=double)                               :: ACOEF, BCOEF, CCOEF, DBDT
REAL(KIND=double)                               :: MUSUBE
REAL(KIND=double)                               :: EPRESS, EU, FSUBE, ES
REAL(KIND=double)                               :: PPRESS, PS, PU, PF
REAL(KIND=double)                               :: DEMUDT, DEMUDN, DEMUDY
REAL(KIND=double)                               :: DEPDT, DEPDN, DEPDY
REAL(KIND=double)                               :: DESDT, DESDN, DESDY
REAL(KIND=double)                               :: DEUDT, DEUDN, DEUDY
REAL(KIND=double)                               :: DPPDN, DPPDT, DPPDY
REAL(KIND=double)                               :: DPSDN, DPSDT, DPSDY
REAL(KIND=double)                               :: DPUDN, DPUDT, DPUDY


END MODULE el_eos_module
