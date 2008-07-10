SUBROUTINE nu_energy_agr_advct_inout_x( imin, imax, nx, nez, nnu, agr_c, &
& agra_c, psi0p )
!-----------------------------------------------------------------------
!
!    File:         nu_energy_agr_advct_inout_x
!    Module:       nu_energy_agr_advct_inout_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To receive variables from radial_ray_module, perform the neutrino
!       energy advection step due to changes in the lapse, and return to
!       radial_ray_module the updated variables.
!
!    Input arguments:
!
!  imin                  : minimum x-array index
!  imax                  : maximum x-array index
!  nx                    : x-array extent
!  nez                   : neutrino energy array extent
!  nnu                   : neutrino flavor array extent
!  agr_c                 : unshifted initial zone-centered lapse function
!  agra_c                : unshifted final zone-centered lapse function
!  psi0p                 : unshifted initial zero angular moments of the
!                           neutrino occupation number
!
!    Output arguments:
!
!  psi0p                 : unshifted updated zero angular moments of the
!                           neutrino occupation number
!
!    Subprograms called:
!  nu_energy_agr_advct_x : performs the x neutrino energy advection due
!   to changes in the lapse
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith

USE e_advct_module, ONLY : ivc_x, agrjmh, agrajmh, psi0, psi0_a
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin             ! minimum x-array index
INTEGER, INTENT(in)              :: imax             ! maximum x-array index
INTEGER, INTENT(in)              :: nx               ! lx-array extent
INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)       :: agr_c    ! unshifted initial zone-centered lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx)       :: agra_c   ! unshifted final zone-centered lapse function

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu) :: psi0p ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: jr_min          ! minimum radial zone index
INTEGER                          :: jr_max          ! maximum radial zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ SET UP FOR NEURINO ENERGY ADVECTION ///
!
!  Load variables received from radial_ray_module into e_advct_module
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize zone-center array size
!-----------------------------------------------------------------------

jr_min                       = imin + 1
jr_max                       = imax + 1

!-----------------------------------------------------------------------
!  Return if nnugpmx    = 0  or  ivc_x    = 0
!-----------------------------------------------------------------------

IF ( nnugpmx == 0  .or.  ivc_x == 0 ) RETURN

!-----------------------------------------------------------------------
!  Transfer psi0 to e_advct_module
!-----------------------------------------------------------------------

psi0(jr_min:jr_max,:,:)      = psi0p(imin:imax,:,:)

!-----------------------------------------------------------------------
!  Transfer state variables to shifted arrays in e_advct_module
!-----------------------------------------------------------------------

agrjmh (jr_min:jr_max)       = agr_c(imin:imax)
agrajmh(jr_min:jr_max)       = agr_c(imin:imax)

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute change in psi0 due to neutrino energy advection
!-----------------------------------------------------------------------

CALL nu_energy_agr_advct_x( jr_min, jr_max )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Updated neutrino distribution function
!-----------------------------------------------------------------------

psi0p(imin:imax,:,:)         = psi0_a(jr_min:jr_max,:,:)

RETURN
END SUBROUTINE nu_energy_agr_advct_inout_x
