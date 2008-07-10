SUBROUTINE radhyd_to_testl( i_ray_dim, nx, ny, nez, nnu, nnc, i_flag )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_testl
!    Module:       radhyd_to_testl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To load variables, call monitor, and print energy and lepton
!       number variables for monitoring energy and lepton conservation.
!
!    Subprograms called:
!  monitor    : computes quantities pertaining to the total energy and
!   lepton number 
!
!    Input arguments:
!  i_ray_dim  : number of rays to put on a processor
!  nx         : x_array extent
!  ny         : y_array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  nnc        : neutrino abundance array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, parallel_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpi

USE edit_module, ONLY: nprint, nlog, data_path
USE parallel_module, ONLY : myid
USE radial_ray_module, ONLY : imin, imax, jmin, jmax, rho_ci, rho_c, t_ci, &
& t_c, ye_ci, ye_c, u_c, x_ef, y_ei, psi0_c, psi1_e, xn_c, be_nuc_rep_c, &
& a_nuc_rep_c, z_nuc_rep_c, nse_c, dtnph, time, ncycle, t_bounce, rhobar, &
& d_omega
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                     :: i_ray_dim       ! number of rays to put on a processor
INTEGER, INTENT(in)                     :: nx              ! x-array extent
INTEGER, INTENT(in)                     :: ny              ! y-array extent
INTEGER, INTENT(in)                     :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)                     :: nnu             ! neutrino flavor array extent
INTEGER, INTENT(in)                     :: nnc             ! composition array extent
INTEGER, INTENT(in)                     :: i_flag          ! call location flag

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                                 :: first = .true.
LOGICAL                                 :: l_monitor       ! monitor write flag

INTEGER                                 :: j               ! y-array index
INTEGER                                 :: i_ray           ! index denoting a specific radial ray
INTEGER                                 :: istat           ! open and close file flag

REAL(KIND=double), DIMENSION(i_ray_dim) :: u_ge_tot        ! total NT gravitational potential energy of i_ray (ergs)
REAL(KIND=double), DIMENSION(i_ray_dim) :: u_ie_tot        ! total internal energy of i_ray (ergs)
REAL(KIND=double), DIMENSION(i_ray_dim) :: u_ke_tot        ! total NT kinetic energy of i_ray (ergs)
REAL(KIND=double), DIMENSION(i_ray_dim) :: u_ne_tot        ! total NT neutrino energy of i_ray (ergs)
REAL(KIND=double), DIMENSION(i_ray_dim) :: e_rad           ! the cumulative material energy entering i_ray grid (ergs)
REAL(KIND=double), DIMENSION(i_ray_dim) :: radtot          ! energy radiated by neutrinos in i_ray (ergs)
REAL(KIND=double), DIMENSION(i_ray_dim) :: u_tot           ! total energy of i_ray (ergs)
REAL(KIND=double), DIMENSION(i_ray_dim) :: elecn           ! total electron number of i_ray
REAL(KIND=double), DIMENSION(i_ray_dim) :: elec_rad        ! the net number of electrons advected into the I_ray grid
REAL(KIND=double), DIMENSION(i_ray_dim) :: nnucr_enu       ! the net number of e-type neutrinos currently residing in I_ray grid
REAL(KIND=double), DIMENSION(i_ray_dim) :: nnucr_enubar    ! the net number of e-type anatineutrinos currently residing in I_ray grid
REAL(KIND=double), DIMENSION(i_ray_dim) :: nnurad_enu      ! the cumulative number of e-type neutrinos emitted by the I_ray grid
REAL(KIND=double), DIMENSION(i_ray_dim) :: nnurad_enubar   ! the cumulative number of e-type antineutrinos emitted by the I_ray grid
REAL(KIND=double), DIMENSION(i_ray_dim) :: totlpn          ! total number of leptons emitted or residing in i_ray
REAL(KIND=double), DIMENSION(ny)        :: d_omega_n       ! solid angles subtended by radial rays normalized to unity

REAL(KIND=double)                       :: u_ge_tot_t      ! total NT gravitational potential energy of model (ergs)
REAL(KIND=double)                       :: u_ie_tot_t      ! total internal energy of model (ergs)
REAL(KIND=double)                       :: u_ke_tot_t      ! total NT kinetic energy of model (ergs)
REAL(KIND=double)                       :: u_ne_tot_t      ! total NT neutrino energy of model (ergs)
REAL(KIND=double)                       :: e_rad_t         ! the cumulative material energy entering model grid (ergs)
REAL(KIND=double)                       :: radtot_t        ! energy radiated by neutrinos in model (ergs)
REAL(KIND=double)                       :: u_tot_t         ! total energy of model (ergs)
REAL(KIND=double)                       :: elecn_t         ! total electron number of model
REAL(KIND=double)                       :: elec_rad_t      ! the net number of electrons advected into the model grid
REAL(KIND=double)                       :: nnucr_enu_t     ! the net number of e-type neutrinos currently residing in model grid
REAL(KIND=double)                       :: nnucr_enubar_t  ! the net number of e-type antineutrinos currently residing in model grid
REAL(KIND=double)                       :: nnurad_enu_t    ! the cumulative number of e-type neutrinos emitted by the model grid
REAL(KIND=double)                       :: nnurad_enubar_t ! the cumulative number of e-type antineutrinos emitted by the model grid
REAL(KIND=double)                       :: totlpn_t        ! total number of leptons emitted or residing in model

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' jmax=',i3,' > ny_max=',i3,' in subroutine radhyd_to_monitor')
 1001 FORMAT (' iflag    cycle       time        rhobar      u_ge_tot     u_ie_tot &
&    u_ke_tot     u_ne_tot      e_rad        radtot       u_tot'/)
 1003 FORMAT (2x,i5,i7,2x,9(es13.5))
 2001 FORMAT (' iflag    cycle       time        rhobar       elecn       elec_rad &
&     tot_enu    tot_enubar     rad_enu    rad_enubar     totlpn'/)
 2003 FORMAT (2x,i5,i7,2x,9(es13.5))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                     \\\\\ INITIALIZE /////
!
!-----------------------------------------------------------------------

l_monitor            = .false.

IF ( first ) THEN
  first              = .false.
  l_monitor          = .true.
END IF ! first

!-----------------------------------------------------------------------
!
!                 \\\\\ MONITOR WRITE CRITERIA /////
!
!-----------------------------------------------------------------------

IF ( t_bounce < 0.004d0 ) THEN
  IF ( MOD( ncycle, 1 ) == 0 ) l_monitor = .true.
ELSE
  IF ( MOD( ncycle, 1 ) == 0 ) l_monitor = .true.
END IF

IF ( .not. l_monitor ) RETURN

!-----------------------------------------------------------------------
!
!          \\\\\ TRANSFER VARIABLES TO AND FROM MONITOR /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!              ||||| Begin loop over radial rays |||||
!-----------------------------------------------------------------------

DO i_ray = 1,i_ray_dim

!-----------------------------------------------------------------------
!  Transfer variables to monitor
!-----------------------------------------------------------------------

  CALL monitor( imin, imax, i_ray, i_ray_dim, nx, nez, nnu, nnc, rho_ci, &
& rho_c, t_ci, t_c, ye_ci, ye_c, rhobar, x_el, u_c, psi0_c, psi1_e, dtnph, &
& time, ncycle, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, nse_c, &
& u_ge_tot(i_ray), u_ie_tot(i_ray), u_ke_tot(i_ray), u_ne_tot(i_ray), &
v
END DO

!-----------------------------------------------------------------------
!
!                \\\\\ SUM OVER RADIAL RAYS /////
!
!-----------------------------------------------------------------------

IF ( jmax == jmin ) THEN

  u_ge_tot_t           = u_ge_tot     (1)
  u_ie_tot_t           = u_ie_tot     (1)
  u_ke_tot_t           = u_ke_tot     (1)
  u_ne_tot_t           = u_ne_tot     (1)
  e_rad_t              = e_rad        (1)
  radtot_t             = radtot       (1)
  u_tot_t              = u_tot        (1)
  elecn_t              = elecn        (1)
  elec_rad_t           = elec_rad     (1)
  nnucr_enu_t          = nnucr_enu    (1)
  nnucr_enubar_t       = nnucr_enubar (1)
  nnurad_enu_t         = nnurad_enu   (1)
  nnurad_enubar_t      = nnurad_enubar(1)
  totlpn_t             = totlpn       (1)

ELSE ! jmax /= jmin

  u_ge_tot_t           = zero
  u_ie_tot_t           = zero
  u_ke_tot_t           = zero
  u_ne_tot_t           = zero
  e_rad_t              = zero
  radtot_t             = zero
  u_tot_t              = zero
  elecn_t              = zero
  elec_rad_t           = zero
  nnucr_enu_t          = zero
  nnucr_enubar_t       = zero
  nnurad_enu_t         = zero
  nnurad_enu_t         = zero
  nnurad_enubar_t      = zero
  totlpn_t             = zero

  d_omega_n(jmin:jmax) = d_omega(jmin:jmax)/frpi

  u_ge_tot_t           = SUM(u_ge_tot     (jmin:jmax) * d_omega_n(jmin:jmax))
  u_ie_tot_t           = SUM(u_ie_tot     (jmin:jmax) * d_omega_n(jmin:jmax))
  u_ke_tot_t           = SUM(u_ke_tot     (jmin:jmax) * d_omega_n(jmin:jmax))
  u_ne_tot_t           = SUM(u_ne_tot     (jmin:jmax) * d_omega_n(jmin:jmax))
  e_rad_t              = SUM(e_rad        (jmin:jmax) * d_omega_n(jmin:jmax))
  radtot_t             = SUM(radtot       (jmin:jmax) * d_omega_n(jmin:jmax))
  u_tot_t              = SUM(u_tot        (jmin:jmax) * d_omega_n(jmin:jmax))
  elecn_t              = SUM(elecn        (jmin:jmax) * d_omega_n(jmin:jmax))
  elec_rad_t           = SUM(elec_rad     (jmin:jmax) * d_omega_n(jmin:jmax))
  nnucr_enu_t          = SUM(nnucr_enu    (jmin:jmax) * d_omega_n(jmin:jmax))
  nnucr_enubar_t       = SUM(nnucr_enubar (jmin:jmax) * d_omega_n(jmin:jmax))
  nnurad_enu_t         = SUM(nnurad_enu   (jmin:jmax) * d_omega_n(jmin:jmax))
  nnurad_enubar_t      = SUM(nnurad_enubar(jmin:jmax) * d_omega_n(jmin:jmax))
  totlpn_t             = SUM(totlpn       (jmin:jmax) * d_omega_n(jmin:jmax))

END IF

!-----------------------------------------------------------------------
!
!                   \\\\\ EDIT QUANTITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Energy edit
!-----------------------------------------------------------------------

OPEN (UNIT=27, FILE=TRIM(data_path)//'/Run_Log/Energy_cons',STATUS='new', &
& IOSTAT=istat)
IF ( istat == 0 ) WRITE (27,1001)
IF ( istat /= 0 ) OPEN (UNIT=27, FILE=TRIM(data_path)//'/Run_Log/Energy_cons', &
& STATUS='old',POSITION='append')
IF ( first ) WRITE (27,1001)
WRITE (27,1003) i_flag, ncycle, time, rhobar(1), u_ge_tot_t, u_ie_tot_t, u_ke_tot_t, &
& u_ne_tot_t, e_rad_t, radtot_t, u_tot_t
CLOSE (unit=27,status='keep')

!-----------------------------------------------------------------------
!  Lepton edit
!-----------------------------------------------------------------------

OPEN (UNIT=27, FILE=TRIM(data_path)//'/Run_Log/Lepton_cons', STATUS='new', &
& IOSTAT=istat)
IF ( istat == 0 ) WRITE (27,2001)
IF ( istat /= 0 ) OPEN (UNIT=27, FILE=TRIM(data_path)//'/Run_Log/Lepton_cons', &
& STATUS='old',POSITION='append')
WRITE (27,2003) i_flag, ncycle, time, rhobar(1), elecn_t, elec_rad_t, nnucr_enu_t, &
& nnucr_enubar_t, nnurad_enu_t, nnurad_enubar_t, totlpn_t
CLOSE (unit=27,status='keep')

RETURN
END SUBROUTINE radhyd_to_testl
