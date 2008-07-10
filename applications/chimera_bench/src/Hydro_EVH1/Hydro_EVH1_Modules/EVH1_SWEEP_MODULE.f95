MODULE evh1_sweep      
!=======================================================================
!  (formerly sweep.h)
!  data structures used in sweeps, dimensioned max_12 (in
!   load_array_module)
!
! jun90 gbl
! sep90 jmb  add forces
! oct00 wrh  add composition
! jan01 wrh  module-ized
!=======================================================================

USE kind_module

INTEGER, parameter :: nmf =2   ! Number of constituents

!-----------------------------------------------------------------------
!  VH1 arrays
!-----------------------------------------------------------------------

CHARACTER (len=1) sweep              ! direction of sweep: x,y,z

LOGICAL, ALLOCATABLE, DIMENSION(:)           :: l_shock   ! true if a shock is present; otherwise false
LOGICAL, ALLOCATABLE, DIMENSION(:)           :: l_rho     ! true if the densiy is above a preset value

INTEGER                                      :: nmin      ! number of first REAL(KIND=double) zone  
INTEGER                                      :: nmax      ! number of last REAL(KIND=double) zone  
INTEGER                                      :: nshk      ! padded index above which total energy is evolved 

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: r         ! padded density (g cm^{-3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: p_mat     ! padded material pressure [ergs cm^{-3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: p_nu      ! padded neutrino pressure [ergs cm^{-3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: p         ! padded total pressure [ergs cm^{-3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: e         ! padded total energy (ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: u         ! padded velocity (cm s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: v         ! padded velocity (cm s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: w         ! padded velocity (cm s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: xa        ! padded advanced grid edges
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: xa0       ! padded original grid edges
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dx        ! padded advanced grid thickness
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dx0       ! padded original grid thickness
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dvol      ! padded advanced grid volumes
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dvol0     ! padded original grid volumes
REAL(KIND=double)                            :: radius    ! rad. (in cm) of current column of angular zones
REAL(KIND=double)                            :: theta     ! the angular coordinate
REAL(KIND=double)                            :: ctheta    ! cosine of theta
REAL(KIND=double)                            :: stheta    ! sine of theta

!-----------------------------------------------------------------------
!  EVH1 additions
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: temp      ! padded temperature (in MeV)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: entrop    ! padded entropy (in k_Boltz)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ei        ! padded internal energy [in ergs g^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: e_nu      ! padded neutrino energy [in ergs g^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: e_v       ! padded total energy minus kinetic energy (in ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ge        ! padded gamma_e = 1 + P/(ei*rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: gc        ! padded gamma_c = (d ln P /d ln rho)= (Svel)**2 /P*rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ye        ! padded the electron fraction
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: xn        ! padded neutron mass fraction (= neutron abundance)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: xp        ! padded proton mass fraction  (= proton abundance)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: flat_s    ! padded flat variable, indicating presence of shock
REAL(KIND=double)                            :: xmin      ! minimum value of x-coordinate
REAL(KIND=double)                            :: xmax      ! maximum value of x-coordinate
REAL(KIND=double)                            :: ymin      ! minimum value of y-coordinate
REAL(KIND=double)                            :: ymax      ! maximum value of y-coordinate
REAL(KIND=double)                            :: zmin      ! minimum value of z-coordinate
REAL(KIND=double)                            :: zmax      ! maximum value of z-coordinate

!-----------------------------------------------------------------------
! MGFLD additions
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rhobar    ! padded mean density (g cm^{-3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: egrav     ! gravitational potential energy (in ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: degrav    ! gravitational potential energy difference (in ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: egrav0    ! gravitational p.e. before Lagrangian update
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ekin      ! kinetic energy (in ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: nu_strc   ! zone-centered neutrino stress (dynes g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: nu_stre   ! zone-edged neutrino stress (dynes g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: g_force_c ! zone-centered gravitational force (dynes g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: g_pot_c   ! zone-centered gravitational potential (ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: g_force_e ! zone-edged gravitational force (dynes g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: g_pot_e   ! zone-edged gravitational potential (ergs g^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: e_nu_c    ! padded angular averaged neutrino energy density (ergs cm^{-3})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f_nu_e    ! padded angular averaged neutrino energy flux (ergs cm^{-2} s^{-1})
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: lapse_c   ! padded zone-centered lapse function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: lapse_e   ! padded zone-edged lapse function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: u_edge    ! padded zone-edged velocity (used for moving grid option)

END module evh1_sweep
