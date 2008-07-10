SUBROUTINE luminosity_MD( imin, imax, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& nx, nez, nnu, r, psi1, unue, dunue, lum )
!-----------------------------------------------------------------------
!
!    File:         luminosity_MD
!    Module:       luminosity_MD
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To calculate the neutrino luminosities at each radial zone interface.
!
!    Subprograms called:
!  flux            : Computes the neutrino flux as a funciton of i and n
!
!    Input arguments:
!  imin            : inner unshifted radial zone number
!  imax            : outer unshifted radial zone number
!  ij_ray          : j-index of a radial ray
!  ik_ray          : k-index of a radial ray
!  ij_ray_dim      : number of y-zones on a processor before swapping
!  ik_ray_dim      : number of z-zones on a processor before swapping
!  nx              : x-array extent
!  nez             : neutrino energy array extent
!  nnu             : neutrino flavor array extent
!  r               : radius (cm)
!  psi1            : first moment of the neutrino distribution
!  unue            : radial zone-edged neutrino energy (MeV)
!  dunue           :radial zone-edged neutrino energy zone width (MeV)
!
!    Output arguments:
!  lum             : luminosity at each zone interface (foes)
!
!    Include files:
!  kind_module, array_module, physcnst_module
!  edit_module, mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpi, ecoef
USE physcnst_module, ONLY : ergfoe, cvel

USE edit_module, ONLY : nprint
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: imin          ! inner unshifted radial zone number
INTEGER, INTENT(in)               :: imax          ! outer unshifted radial zone number
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)             :: r        ! radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1  ! zero moment of f_neutrino
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim)     :: unue  ! radial zone-edged neutrino energy (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim)     :: dunue ! radial zone-edged neutrino energy zone width (MeV)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)    :: lum  ! neutrino luminosity

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: i             ! radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double), DIMENSION(nx)  :: area          ! area (cm^{2})

!-----------------------------------------------------------------------
!        Initialize.
!-----------------------------------------------------------------------

lum(imin:imax+1,:,ij_ray,ik_ray) = zero


!-----------------------------------------------------------------------
!        Return if no neutrinos.
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ NEUTRINO LUMINOSITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Compute the luminosities as a function of i and n.
!-----------------------------------------------------------------------

DO i = imin,imax+1
  area(i)                        = frpi * r(i) * r(i)
END DO

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO i = imin,imax+1
      lum(i,n,ij_ray,ik_ray)     = lum(i-imin+1,n,ij_ray,ik_ray) &
&      + area(i) * psi1(i,k,n,ij_ray,ik_ray) * cvel * ecoef * unue(i,k,ij_ray,ik_ray) * unue(i,k,ij_ray,ik_ray) &
&      * unue(i,k,ij_ray,ik_ray) * dunue(i,k,ij_ray,ik_ray) * ergfoe
    END DO ! i
  END DO ! k
END DO ! n

RETURN
END SUBROUTINE luminosity_MD
