SUBROUTINE cnfg_energy( ij_ray, ik_ray, u_ge_tot, u_ie_tot, u_ke_tot, &
& u_ne_tot, u_ke_x_tot, u_ke_y_tot, u_ke_z_tot, radtot, elecn, totlpn, &
& nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         cnfg_energy
!    Module:       cnfg_energy
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/28/04
!
!    Purpose:
!      To compute the total energy of the core configuration.
!
!    Variables that must be passed through common:
!
!    Subprograms called:
!  eqstt_x    : interpolates EOS quantities
!
!    Input arguments:
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!  nx         : x_array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!  u_ge_tot   : total NT gravitational potential energy
!  u_ie_tot   : total internal energy
!  u_ke_tot   : total NT kinetic energy
!  u_ne_tot   : total NT neutrino energy
!  u_ke_x_tot : total NT x-kinetic energy of i_ray (ergs)
!  u_ke_y_tot : total NT y-kinetic energy of i_ray (ergs)
!  u_ke_z_tot : total NT z-kinetic energy of i_ray (ergs)
!  radtot     : total energy radiated by neutrinos
!  elecn      : total number of electrons - positrons residing in the core
!  totlpn     : net total number of leptons residing in the core
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, third, epsilon
USE physcnst_module, ONLY : cvel, g, ergfoe, rmu

USE edit_module, ONLY : nprint, g_pot
USE eos_snc_x_module, ONLY: duesrc
USE mdl_cnfg_module, ONLY : jr_min, jr_max, rho, t, ye, u, v, w, r, dmgrv, &
& dmrst, rstmss
USE nu_dist_module, ONLY : ncoefa, ecoefa, stwt, psi0, psi1, unuinfty, &
& unujinfty, vol, nnucr, unurad, elec_rad, nnurad, gamgr_nu, dc
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)    :: u_ge_tot      ! total gravitational potential energy
REAL(KIND=double), INTENT(out)    :: u_ie_tot      ! total internal energy
REAL(KIND=double), INTENT(out)    :: u_ke_tot      ! total NT kinetic energy
REAL(KIND=double), INTENT(out)    :: u_ne_tot      ! total NT neutrino energy
REAL(KIND=double), INTENT(out)    :: u_ke_x_tot    ! total NT x-kinetic energy
REAL(KIND=double), INTENT(out)    :: u_ke_y_tot    ! total NT y-kinetic energy
REAL(KIND=double), INTENT(out)    :: u_ke_z_tot    ! total NT z-kinetic energy
REAL(KIND=double), INTENT(out)    :: radtot        ! total energy radiated by neutrinos
REAL(KIND=double), INTENT(out)    :: elecn         ! total electron number
REAL(KIND=double), INTENT(out)    :: totlpn        ! total lepton number

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                           :: first = .true.

INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double)                 :: csqinv        ! 1/c^{2}
REAL(KIND=double)                 :: u_int         ! internal energy (ergs/g)
REAL(KIND=double)                 :: duddn         ! d(u_int)/dd
REAL(KIND=double)                 :: dudtn         ! d(u_int)/dt
REAL(KIND=double)                 :: dudyn         ! d(u_int)/dye
REAL(KIND=double)                 :: rnnutt        ! quantity for calculating the neutrino number

REAL(KIND=double), DIMENSION(nx)  :: u_ge_nt       ! gravitational energy 
REAL(KIND=double), DIMENSION(nx)  :: u_ie_nt       ! NR internal energy
REAL(KIND=double), DIMENSION(nx)  :: u_ke_nt       ! NR kinetic energy
REAL(KIND=double), DIMENSION(nx)  :: u_ne_nt       ! NR neutrino energy
REAL(KIND=double), DIMENSION(nx)  :: u_ke_x_nt     ! NR x-kinetic energy
REAL(KIND=double), DIMENSION(nx)  :: u_ke_y_nt     ! NR y-kinetic energy
REAL(KIND=double), DIMENSION(nx)  :: u_ke_z_nt     ! NR z-kinetic energy

REAL(KIND=double), DIMENSION(nx)  :: beta2         ! zone-centered mass
REAL(KIND=double), DIMENSION(nx)  :: gamma2        ! zone-centered mass

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  first                = .false.
  csqinv               = 1.d0/cvel**2
END IF

!-----------------------------------------------------------------------
!  Total energy radiated by neutrinos
!-----------------------------------------------------------------------

radtot                 = SUM( unurad(:,ij_ray,ik_ray) ) * ergfoe

!-----------------------------------------------------------------------
!
!                     \\\\\ TOTAL ENERGY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Gravitational potential energy
!-----------------------------------------------------------------------

u_ge_nt(jr_min:jr_max) = g_pot(jr_min:jr_max) * dmrst(jr_min:jr_max)

!........Get total

u_ge_tot               = SUM( u_ge_nt(jr_min:jr_max) ) * ergfoe

!-----------------------------------------------------------------------
!  Internal energy
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), u_int, duddn, &
&  dudtn, dudyn )
  u_ie_nt(j)           = ( u_int - duesrc(j,ij_ray,ik_ray) ) * dmrst(j) 
END DO ! j

!........Get total

u_ie_tot               = SUM( u_ie_nt(jr_min:jr_max) ) * ergfoe

!-----------------------------------------------------------------------
!  Newtonian kinetic energy
!-----------------------------------------------------------------------

u_ke_nt(jr_min:jr_max) = half * ( u(jr_min:jr_max) * u(jr_min:jr_max)   &
&                      +          v(jr_min:jr_max) * v(jr_min:jr_max)   &
&                      +          w(jr_min:jr_max) * w(jr_min:jr_max) )

u_ke_x_nt(jr_min:jr_max) = half * ( u(jr_min:jr_max) * u(jr_min:jr_max) )
u_ke_y_nt(jr_min:jr_max) = half * ( v(jr_min:jr_max) * v(jr_min:jr_max) )
u_ke_z_nt(jr_min:jr_max) = half * ( w(jr_min:jr_max) * w(jr_min:jr_max) )

!........Get total

u_ke_tot               = SUM( u_ke_nt  (jr_min:jr_max) * dmrst(jr_min:jr_max) ) * ergfoe
u_ke_x_tot             = SUM( u_ke_x_nt(jr_min:jr_max) * dmrst(jr_min:jr_max) ) * ergfoe
u_ke_y_tot             = SUM( u_ke_y_nt(jr_min:jr_max) * dmrst(jr_min:jr_max) ) * ergfoe
u_ke_z_tot             = SUM( u_ke_z_nt(jr_min:jr_max) * dmrst(jr_min:jr_max) ) * ergfoe

!-----------------------------------------------------------------------
!  Neutrino energy
!-----------------------------------------------------------------------

beta2 (jr_min:jr_max)  = u(jr_min:jr_max) * u(jr_min:jr_max)/( cvel * cvel )
gamma2(jr_min:jr_max)  = 1.d0/( 1.d0 - beta2(jr_min:jr_max) )

unujinfty(jr_min:jr_max,:) = zero
unuinfty(:)            = zero
nnucr(:,ij_ray,ik_ray) = zero

!-----------------------------------------------------------------------
!  Neutrino energy as a function of j and n (summed over k)
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max
      unujinfty(j,n)   = unujinfty(j,n) &
&                      + ecoefa(j,k) * psi0(j,k,n) * stwt(n) * vol(j) &
&                      + ecoefa(j,k) * ( psi1(j-1,k,n) + psi1(j,k,n) ) * stwt(n) * vol(j) * u(j)/cvel
    END DO ! j
  END DO ! k
END DO ! n

DO n = 1,nnu
  DO j = jr_min,jr_max
    unujinfty(j,n)     = gamma2(j) * ( unujinfty(j,n) + beta2(j) * unujinfty(j,n)/3.d0 )
  END DO ! j
END DO ! n

!-----------------------------------------------------------------------
!  Neutrino energy as a function of n (summed over j and k)
!-----------------------------------------------------------------------

DO n = 1,nnu
  unuinfty(n)          = SUM( unujinfty(jr_min:jr_max,n) )
END DO ! n

!-----------------------------------------------------------------------
!  Neutrino energy as a function of j (summed over k and n)
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  u_ne_nt(j)           = SUM( unujinfty(j,:) )
END DO ! j

!-----------------------------------------------------------------------
!  Neutrino number as a function of n
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      DO j = jr_min,jr_max
        rnnutt         = ncoefa(j,k) * psi0(j,k,n) * stwt(n) * vol(j)
        nnucr(n,ij_ray,ik_ray) = nnucr(n,ij_ray,ik_ray) + rnnutt
      END DO ! j
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!  Total neutrino energy
!-----------------------------------------------------------------------

u_ne_tot               = SUM( unuinfty(:) ) * ergfoe

!-----------------------------------------------------------------------
!
!             \\\\\ TOTAL LEPTON NUMBER ///
!
!-----------------------------------------------------------------------

elecn                  = SUM( vol(jr_min:jr_max) * rho(jr_min:jr_max) * ye(jr_min:jr_max) )/rmu

totlpn                 = elecn + elec_rad(ij_ray,ik_ray) + nnucr(1,ij_ray,ik_ray) &
&                      + nnurad(1,ij_ray,ik_ray) - nnucr(2,ij_ray,ik_ray)         &
&                      - nnurad(2,ij_ray,ik_ray)

RETURN
END SUBROUTINE cnfg_energy
