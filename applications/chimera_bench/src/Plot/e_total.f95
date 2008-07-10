SUBROUTINE e_total( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         e_total
!    Module:       e_total
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/28/99
!
!    Purpose:
!      To compute the total energy of the core configuration.
!
!    Variables that must be passed through common:
!
!    Subprograms called:
!  eqstt_x : interpolates quantities in the local EOS table
!
!    Input arguments:
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!        none
!
!    Output arguments (common):
!        none
!
!    Include files:
!  array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!  nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nz, nnu
USE numerical_module, ONLY : zero, half, third, epsilon
USE physcnst_module, ONLY : cvel, g, ergfoe, rmu

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, i_eplt, n_eplt, dt_eplot, pdv, data_path
USE eos_snc_x_module, ONLY : duesrc
USE mdl_cnfg_module, ONLY : jr_max, rho, t, ye, u, r, dmgrv, dmrst, rstmss
USE nu_dist_module, ONLY : ncoefa, ecoefa, stwt, psi0, psi1, unuinfty, &
& unujinfty, nnucr, vol, unurad, elec_rad, e_rad, nnurad, gamgr_nu
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY : dtnph, time

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)               :: ik_ray        ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                           :: first = .true.

INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index
INTEGER, PARAMETER                :: i_grav = 2    ! 1, zone-edged; 2, zone-centered

INTEGER                           :: itime         ! quantity for determining when to write to polt file
INTEGER                           :: itimeprev     ! quantity for determining when to write to polt file
INTEGER                           :: istat         ! open/close error flag

REAL(KIND=double)                 :: csqinv        ! 1/c^{2}
REAL(KIND=double)                 :: ugrstm        ! gravitational potential energy at j-1
REAL(KIND=double)                 :: ugrstp        ! gravitational potential energy at j
REAL(KIND=double)                 :: u_int         ! internal energy (ergs/g)
REAL(KIND=double)                 :: duddn         ! d(u_int)/dd
REAL(KIND=double)                 :: dudtn         ! d(u_int)/dt
REAL(KIND=double)                 :: dudyn         ! d(u_int)/dye
REAL(KIND=double)                 :: unucrtt       ! quantity for calculating the neutrino energy
REAL(KIND=double)                 :: rnnutt        ! quantity for calculating the neutrino number
REAL(KIND=double)                 :: radtot        ! total energy radiated by neutrinos (ergs)
REAL(KIND=double)                 :: radtot_foe    ! total energy radiated by neutrinos (foes)
REAL(KIND=double)                 :: utot_nt       ! total NT energy
REAL(KIND=double)                 :: utot_nts      ! total NT energy including energy glitches
REAL(KIND=double)                 :: utot_gr       ! total GR energy
REAL(KIND=double)                 :: utot_grs      ! total GR energy including energy glitches
REAL(KIND=double)                 :: duesrt        ! energy glitches

REAL(KIND=double)                 :: elecn         ! electron number
REAL(KIND=double)                 :: totlpn        ! lepton number

REAL(KIND=double)                 :: tmult         ! quantity for determining when to write to polt file

REAL(KIND=double), DIMENSION(nz)  :: u_ge_nt       ! NR gravitational energy 
REAL(KIND=double), DIMENSION(nz)  :: u_ie_nt       ! NR internal energy
REAL(KIND=double), DIMENSION(nz)  :: u_ke_nt       ! NR kinetic energy
REAL(KIND=double), DIMENSION(nz)  :: u_ne_nt       ! NR neutrino energy
REAL(KIND=double), DIMENSION(nz)  :: r_center      ! zone-centered radius
REAL(KIND=double), DIMENSION(nz)  :: m_center      ! zone-centered mass

 1001 format (' i_grav in e_total3 is neither 1 or 2')
 6001 format (' ncycle=',i7,' time=',1pe14.7,' utot_gr=',1pe10.3,' utot_nt=',1pe10.3,' utot_grs=',1pe10.3, &
& ' utot_nts=',1pe10.3,' n_l=',1pe10.3)
 6003 format (1pe15.7,4(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  csqinv           = 1.d0/cvel**2
END IF

!-----------------------------------------------------------------------
!  Total energy radiated by neutrinos
!-----------------------------------------------------------------------

radtot             = zero
DO n = 1,nnu
  radtot           = radtot + unurad(n,ij_ray,ik_ray)
END DO
radtot_foe         = radtot * ergfoe

!-----------------------------------------------------------------------
!  Total energy glitches
!-----------------------------------------------------------------------

duesrt             = zero

DO j = 2,jr_max
  duesrt           = duesrt + duesrc(j,ij_ray,ik_ray) * dmrst(j) * ergfoe
END DO

!-----------------------------------------------------------------------
!
!                    \\\\\ NEWTONIAN ENERGY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Newtonian gravitational potential energy
!-----------------------------------------------------------------------

ugrstp             = zero

IF ( i_grav == 1 ) THEN
  DO j = 2,jr_max
    ugrstm         = ugrstp
    ugrstp         = -rstmss(j) * g/r(j)
    u_ge_nt(j)     = half * ( ugrstm + ugrstp ) * dmrst(j)
  END DO
ELSE IF ( i_grav == 2 ) THEN
  DO j = 2,jr_max
    r_center(j)    = ( half * ( ( r(j) + epsilon )**3 + ( r(j-1) + epsilon )**3 ) )**third
    m_center(j)    = half * ( rstmss(j) + rstmss(j-1) )
  END DO
  DO j = 2,jr_max
    ugrstm         = ugrstp
    ugrstp         = -m_center(j) * g/r_center(j)
    u_ge_nt(j)     = ugrstp * dmrst(j)
  END DO
ELSE
  WRITE (nprint,1001)
  STOP
END IF

!-----------------------------------------------------------------------
!  Internal energy
!-----------------------------------------------------------------------

DO j = 2,jr_max
  CALL eqstt_x( 2, j, ij_ray,ik_ray, rho(j), t(j), ye(j), u_int, duddn, dudtn, dudyn )
  u_ie_nt(j)       = u_int * dmrst(j) 
END DO

!-----------------------------------------------------------------------
!  Newtonian kinetic energy
!-----------------------------------------------------------------------

DO j = 2,jr_max
  u_ke_nt(j)  = 0.25d0 * ( u(j) * u(j) + u(j-1) * u(j-1) ) * dmrst(j)
END DO

!-----------------------------------------------------------------------
!  Neutrino energy
!-----------------------------------------------------------------------

DO n = 1,nnu

  unuinfty(n)         = zero
  nnucr(n,ij_ray,ik_ray)      = zero

  DO j = 2,jr_max

    unujinfty(j,n)    = zero

    DO k = 1,nnugp(n)
      unucrtt         = ecoefa(j,k) * psi0(j,k,n) * stwt(n) * vol(j) * half * ( gamgr_nu(j) + gamgr_nu(j-1) ) &
&                     + ecoefa(j,k) * ( psi1(j-1,k,n) + psi1(j,k,n) ) * stwt(n) * vol(j) * u(j)/cvel
      unujinfty(j,n)  = unujinfty(j,n) + unucrtt
      unuinfty(n)     = unuinfty(n) + unucrtt
      rnnutt          = ncoefa(j,k) * psi0(j,k,n) * stwt(n) * vol(j)
      nnucr(n,ij_ray,ik_ray)  = nnucr(n,ij_ray,ik_ray) + rnnutt
    END DO ! k

  END DO ! j

END DO ! n

u_ne_nt               = zero
DO n = 1,nnu
  DO j = 2,jr_max
    u_ne_nt(j)        = u_ne_nt(j) + unujinfty(j,n)
  END DO ! j
END DO ! n

!-----------------------------------------------------------------------
!  Nt: totge + totie + totke + totne + rad + pdv
!-----------------------------------------------------------------------

utot_nt            = zero

DO j = 2,jr_max
  utot_nt          = utot_nt + ( u_ge_nt(j) + u_ie_nt(j) + u_ke_nt(j) + u_ne_nt(j) ) * ergfoe
END DO

utot_nt            = utot_nt + ( radtot + pdv(jr_max+1) + e_rad(ij_ray,ik_ray) ) * ergfoe

!-----------------------------------------------------------------------
!  Nt: totge + totie + totke + totne + rad + pdv - duesrc
!-----------------------------------------------------------------------

utot_nts           = utot_nt - duesrt

!-----------------------------------------------------------------------
!
!                     \\\\\ GR ENERGY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  GR: totge + totie + totke + totne + rad + pdv
!-----------------------------------------------------------------------

utot_gr            = zero

DO j = 2,jr_max
  utot_gr          = utot_gr + ( cvel**2 * ( dmgrv(j) ) - cvel**2 * dmrst(j) ) * ergfoe
END DO

utot_gr            = utot_gr + ( radtot + pdv(jr_max+1) ) * ergfoe

!-----------------------------------------------------------------------
!  GR: totge + totie + totke + totne + rad + pdv - duesrc
!-----------------------------------------------------------------------

utot_grs           = utot_gr - duesrt

!-----------------------------------------------------------------------
!
!                   \\\\\ TOTAL LEPTON NUMBER /////
!
!-----------------------------------------------------------------------

elecn                = zero

DO j = 2,jr_max
  elecn              = elecn + dmrst(j) * ye(j)/rmu
END DO

totlpn               = nnurad(1,ij_ray,ik_ray) + nnucr(1,ij_ray,ik_ray) + elecn - nnurad(2,ij_ray,ik_ray) &
&                    - nnucr(2,ij_ray,ik_ray) + elec_rad(ij_ray,ik_ray)

!-----------------------------------------------------------------------
!
!                 \\\\\ EDIT TOTAL ENERGY DATA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Print to nprint
!-----------------------------------------------------------------------

WRITE (nprint,6001) ncycle,time,utot_gr,utot_nt,utot_grs,utot_nts,totlpn 

!-----------------------------------------------------------------------
!  Return if i_eplt = 0
!-----------------------------------------------------------------------

IF ( i_eplt == 0 ) RETURN

!-----------------------------------------------------------------------
!  Evaluate criterion for writing entries to total energy plot files.
!-----------------------------------------------------------------------

tmult              = 1.d+3/dt_eplot
itime              = int( time * tmult )
itimeprev          = int( ( time - dtnph ) * tmult )

!-----------------------------------------------------------------------
!  Write to e_chk.d file if criteria satisfies
!-----------------------------------------------------------------------

IF ( itime > itimeprev  .or.  ncycle == 0 ) THEN

!-----------------------------------------------------------------------
!  Open e_chk.d file
!-----------------------------------------------------------------------

  OPEN (UNIT=n_eplt,FILE=TRIM(data_path)//'/e_chk.d',STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=n_eplt,FILE=TRIM(data_path)//'/e_chk.d', &
&  STATUS='old', POSITION='append')

!-----------------------------------------------------------------------
!  Write to e_chk.d file
!-----------------------------------------------------------------------

  WRITE (n_eplt,6003) time,utot_gr,utot_nt,utot_grs,utot_nts

!-----------------------------------------------------------------------
!  Close e_chk.d file
!-----------------------------------------------------------------------

  CLOSE (UNIT=n_eplt, STATUS='keep')

END IF ! itime > itimeprev or ncycle = 0

RETURN
END SUBROUTINE e_total
