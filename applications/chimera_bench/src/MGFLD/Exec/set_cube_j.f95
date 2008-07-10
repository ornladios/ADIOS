SUBROUTINE set_cube_j( j, rho, t, ye, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         set_cube_j
!    Module:       set_cube_j
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/11/02
!
!    Purpose:
!      To compute reset the equation of state and rate functions at the cube corners
!       for the jth radial zone.
!
!    Input variables passed through the calling statement
!  j         : radial zone index
!  rho       : shifted matter density array (g/cm**3).
!  t         : shifted matter matter temperature array (K).
!  ye        : shifted matter matter electron fraction array.
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!
!    Output variables passed through the calling statement
!      none
!
!    Subprograms called:
!  esrgn_x   : recomputes, if necessary, table of EOS quantities
!  eqstt_x   : recomputes interpolated EOS quantities
!  gammaj_x  : recomputes interpolated adiabatic gammas
!  abemset   : recomputes, if necessary, table of neutrino emission and
!               absorption rates
!  bremset   : recomputes, if necessary, table of nucleon-nucleus
!               bremsstrahlung rates
!  pairset   : recomputes, if necessary, table of pair annihilation rates
!  scataset  : recomputes, if necessary, table of neutrino-nucleus
!               inelastic scattering rates
!  scateset  : recomputes, if necessary, table of neutrino-electron
!               elastic scattering rates
!  scatiset  : recomputes, if necessary, table of neutrino-nucleon and
!               nucleus isoenergetic scattering rates
!  scatnnset : recomputes, if necessary, table of neutrino-nucleon
!               inelastic scattering rates
!  scatnAset : recomputes, if necessary, table of neutrino-nucleus
!               inelastic scattering rates
!  scatnset  : recomputes, if necessary, table of neutrino-nucleon
!               elastic scattering rates
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_x_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: nx
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY: aesv, aesvd, aesvt, aesvy, gam1, gam2, gam3, duesrc
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j                ! radial zone index
INTEGER, INTENT(in)              :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray           ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho  ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t    ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ye   ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: i                ! equation of state variable index
INTEGER                          :: istat            ! allocation status

REAL(KIND=double)                :: duddp            ! derivative with respect to density
REAL(KIND=double)                :: dudtp            ! derivative with respect to temperature
REAL(KIND=double)                :: dudyp            ! derivative with respect to electron fraction

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: umat ! internal energy (ergs/gm)

 1001 FORMAT (' Allocation problem for array ',a10,' in set_cube_j')
 2001 FORMAT (' Deallocation problem for array ',a10,' in set_cube_j')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (umat(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'umat      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

umat                    = zero

!-----------------------------------------------------------------------
!
!             \\\\\ RESET EOS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute internal energy using old cube corners
!-----------------------------------------------------------------------

CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), umat(j), duddp, &
& dudtp, dudyp )
  
!-----------------------------------------------------------------------
!  Recompute thermodynamic quantities
!-----------------------------------------------------------------------

CALL esrgn_x( j, ij_ray,ik_ray, rho(j), t(j), ye(j) )
DO i = 1,12
  CALL eqstt_x( i, j, ij_ray, ik_ray, rho(j), t(j), ye(j), aesv(j,i,ij_ray,ik_ray), &
& aesvd(j,i,ij_ray,ik_ray), aesvt(j,i,ij_ray,ik_ray), aesvy(j,i,ij_ray,ik_ray) )
END DO
CALL gammaj_x( j, rho(j), t(j), ij_ray, ik_ray, gam1(j,ij_ray,ik_ray), &
& gam2(j,ij_ray,ik_ray), gam3(j,ij_ray,ik_ray) )
  
!-----------------------------------------------------------------------
!  Compute internal energy discontinuities
!-----------------------------------------------------------------------

duesrc(j,ij_ray,ik_ray) = duesrc(j,ij_ray,ik_ray) + aesv(j,2,ij_ray,ik_ray) - umat(j)

!-----------------------------------------------------------------------
!
!             \\\\\ RESET NEITRONO RATES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Recompute neutrino rates
!-----------------------------------------------------------------------

IF ( nnugpmx > 0 ) THEN
  CALL abemset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL bremset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL pairset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL scataset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL scateset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL scatiset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL scatnset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL scatnAset( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL scatnnset( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
END IF

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (umat, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'umat      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE set_cube_j
