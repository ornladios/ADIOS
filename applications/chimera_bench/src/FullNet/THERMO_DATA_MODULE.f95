MODULE thermo_data
!===============================================================================
!  This module contains the thermodynamic trajectory which the network follows
!===============================================================================

USE kind_module, ONLY : double

INTEGER, PARAMETER   :: nhmx= 500 ! The max number of thermo points
INTEGER              :: nh        ! The actual number of points

REAL(KIND=double)    :: tstart,tstop,th(nhmx),t9h(nhmx),rhoh(nhmx)

END module thermo_data
