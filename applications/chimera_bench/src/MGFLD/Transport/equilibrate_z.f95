SUBROUTINE equilibrate_z( kmin, kmax, nx, nz, nez, nnu, ki_ray, kj_ray, &
& ij_ray_dim, k_ray_dim, i_radial, rho_z, t_z, ye_z, psi0_z, rhobar,    &
& dtnph, rho_equilibrate, agr_z )
!-----------------------------------------------------------------------
!
!    File:         equilibrate_z
!    Module:       equilibrate_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/19/06
!
!    Purpose:
!      To equilibrate neutrinos with matter.
!
!    Subprograms called:
!  abemset_z       : sets the neutrino absorption and emission rates at cube corners
!  abemrate_z      : computes neutrino absorption and emission rates given iz,k,n
!  dgesv           : Solve for the recursion coefficients relating the unknowns to the knowns
!  eqstt_z         : computes equation of state variables
!  tgvndeye_z      : computes the temperature given rho, ei, and ye
!  tgvndsye_z      : computes the temperature given rho, s, and ye
!
!    Input arguments:
!  kmin            : minimum angular zone to equilibrate matter and neutrinos
!  kmax            : maximum angular zone to equilibrate matter and neutrinos
!  nx              : x-array extent
!  nz              : z-array extent
!  ki_ray          : x (radial) index of a specific azimuthal ray
!  kj_ray          : y (angular) index of a specific azimuthal ray
!  ij_ray_dim      : the number of radial zones on a processor before swapping with y
!  k_ray_dim       : the number of z-zones on a processor after swapping with z
!  i_radial        : the unshifted radial zone (angular ray) corresponding to ki_ray
!  nez             : neutrino energy array extent
!  nnu             : neutrino flavor array extent
!  rho_z           : unshifted density [g cm^{-3}]
!  t_z             : unshifted temperature (K)
!  ye_z            : unshifted electron fraction
!  psi0_z          : unshifted zero angular moments of the neutrino occupation number
!  rhobar          : mean density of an angular ray [g cm^{-3}]
!  dtnph           : time step
!  rho_equilibrate : density above whicn neutrinos equilibrated with matter
!  agr_z           : lapse function
!
!    Output arguments:
!  psi0_z          : unshifted zero angular moments of the neutrino occupation number
!  t_z             : unshifted temperature (K)
!  ye_z            : unshifted electron fraction
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  abem_y_module, edit_module, evh1_global, it_tol_module,
!  nu_dist_module, nu_energy_grid_module, parallel_module, prb_cntl_module
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc_y
USE numerical_module, ONLY : zero, one, ncoef, ecoef
USE physcnst_module, ONLY : rmu, cvel

USE abem_z_module, ONLY : emis, emisy, absor, absory
USE edit_module, ONLY : nprint, nlog
USE evh1_global, ONLY : degen
USE it_tol_module, ONLY : iternu, tolpsimin, tolnupsi, tolnuye
USE nu_dist_module, ONLY : stwt
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx, unui, dunui
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY : iaefnp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: kmin            ! minimum angular zone index
INTEGER, INTENT(in)              :: kmax            ! maximum angular zone index
INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: nz              ! z-array extent
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: kj_ray          ! z (azimuthal) index of a specific angular ray
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: k_ray_dim       ! number of z (azimuthal) zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial        ! the unshifted radial zone corresponding to ki_ray, kj_ray

REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim) :: rho_z            ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                      :: rhobar           ! mean density of an angular ray (g cm^{-3})
REAL(KIND=double), INTENT(in)    :: dtnph           ! time step
REAL(KIND=double), INTENT(in)    :: rho_equilibrate ! mean density of an angular ray (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim) :: agr_z            ! lapse function

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim) :: t_z            ! temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim) :: ye_z           ! electron fraction
REAL(KIND=double), INTENT(inout), DIMENSION(nz,nez,nnu,ij_ray_dim,k_ray_dim) :: psi0_z ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name
CHARACTER (len=1), PARAMETER     :: consrve = 's'

LOGICAL                          :: endit

INTEGER                          :: iz              ! radial zone index
INTEGER                          :: k               ! neutrino energy index
INTEGER                          :: kp              ! neutrino energy index
INTEGER                          :: it              ! iteration index
INTEGER                          :: jiterp0         ! radial zone for which psi0 iteration fails
INTEGER                          :: kiterp0         ! neutrino energy zone for which psi0 iteration fails
INTEGER                          :: niterp0         ! neutrino flavor for which psi0 iteration fails
INTEGER                          :: jiterye         ! radial zone number for which ye iteration fails

INTEGER                          :: istat           ! allocation status
INTEGER                          :: j_radial        ! shifted radial coordinate index

INTEGER                          :: i_ray           ! radial index of the angular ray
INTEGER                          :: k_ray           ! azimuthal index of the angular ray

INTEGER                          :: ne_1            ! nnugp(1)

INTEGER                          :: n_linear        ! number of linear equations to be solved by dgesv
INTEGER                          :: n_rhs           ! number of right-hand sides of the linear equations
INTEGER                          :: l_leading_A     ! leading dimension of the coefficient matrix A_coef
INTEGER                          :: l_leading_B     ! leading dimension of the coefficient matrix C_coef
INTEGER, DIMENSION(nez+1)        :: i_piv           ! the oivit indices that define the permutation matrix P
INTEGER                          :: info            ! 0: successful exit, < 0: |info| argument illegal, > 0: u(info,info) = 0

REAL(KIND=double)                :: dtnph_equil     ! time step used for equilibration
REAL(KIND=double)                :: cdt             ! c * dtnph_equil
REAL(KIND=double)                :: cdtinv          ! 1/c * dtnph_equil

REAL(KIND=double)                :: e               ! inernal energy per unit mass
REAL(KIND=double)                :: dedd            ! d(e)/d(rho)
REAL(KIND=double)                :: dedt            ! d(e)/d(t)
REAL(KIND=double)                :: dedy            ! d(e)/d(ye)

REAL(KIND=double)                :: s_i             ! shifted array of entropy
REAL(KIND=double)                :: dsdd            ! d(s_i)/d(rho)
REAL(KIND=double)                :: dsdt            ! d(s_i)/d(T)
REAL(KIND=double)                :: dsdy            ! d(s_i)/d(Ye)

REAL(KIND=double)                :: e_eta           ! (electron Fermi energy)/KT

REAL(KIND=double), PARAMETER     :: psi0amax = 1.d0 ! maximum value of psi0
REAL(KIND=double), PARAMETER     :: psi0amin = 0.d0 ! minimum value of psi0

REAL(KIND=double), DIMENSION(2), PARAMETER   :: n_a = (/ 1.d0, -1.d0 /)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho          ! density (g/cm^3)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t            ! 1-D shifted array of temperature [K]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t_new        ! 1-D shifted array of updated temperature [K]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye           ! 1-D shifted array of electron fraction
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye_new       ! 1-D shifted array of updated electron fraction
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: psi0         ! 3-D shifted array of psi0

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0         ! coefficients of zero-order terms on left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_p0      ! coefficients of dpsi0(iz) terms left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_y       ! coefficients of dYe(iz) terms left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: RHS0         ! coefficients of rigit-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: RHS0_y       ! coefficients of d/dYe rigit-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: RHS0_p0      ! coefficients of d/dpsi0 rigit-hand side of p0 equation

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye_0         ! initial value Ye
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dye_iph      ! increment of Ye for the ith iteration
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: psi0_i       ! initial value of psi0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: dpsi0_iph    ! increment of psi0 for the ith iteration
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: cf_e         ! energy transfer coefficients
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: cf_ye        ! lepton transfer coefficients
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: de_i         ! energy change due to n-neutrinos and ye
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: de           ! total energy change due to neutrinos
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: d_ye         ! tchange in the electron fraction due to emission and absorption

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: A            ! coefficient matrix for the iz terms in iz-equation, all iz
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: C            ! constant terms in iz-equation, all iz
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: u            ! solution vector of the linear equations

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: R            ! recursion matrix
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: A_coef       ! Coefficient matrix of the variables to be solved
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: C_coef       ! Coefficient matrix of quantities on the right-hand side

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: ncoefa       ! neutrino number density coefficient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: ecoefa       ! neutrino energy density coefficient

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in equilibrate_z')
 1501 FORMAT (' Iteration attempt',i4,' myid=',i4,' ki_ray=',i4,' kj_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=',es11.3,' ye(',i4,')=',es11.3)
 1503 FORMAT (' Iteration attempt',i4,' myid=',i4,' ki_ray=',i4,' kj_ray=',i4, &
&' dpsi0_iph(',i4,',',i4,',',i4,')=', es11.3,' D_psi=',es11.3,                 &
&' psi0(',i4,',',i4,',',i4,')=', es11.3)
 1505 FORMAT (' Iteration attempt',i4,' myid=',i4,' ki_ray=',i4,' kj_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=', es11.3,' ye(',i4,')=',es11.3,           &
& ' dpsi0_iph(',i4,',',i4,',',i4,')=',es11.3,' D_psi=',es11.3,                 &
& ' psi0(',i4,',',i4,',',i4,')=', es11.3)
 2001 FORMAT (' Deallocation problem for array ',a10,' in equilibrate_z')
 7001 FORMAT (' Convergence failure in subroutine equilibrate_z, i_ray=',i4,   &
& ' k_ray=',i4,' myid=',i4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if iaefnp = 0 or nnugp(1) = 0
!-----------------------------------------------------------------------

IF ( iaefnp == 0  .or.  nnugp(1) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Return if rhobar(i_radial) < rho_equilibrate
!-----------------------------------------------------------------------

IF ( rhobar(i_radial) < rho_equilibrate ) RETURN

!-----------------------------------------------------------------------
!  Shifted radial coordinate index
!-----------------------------------------------------------------------

j_radial                 = i_radial + 1

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rho(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_new(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_new     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_new(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_new    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (LHS0(nz,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_p0(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_y(nz,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_y    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (RHS0(nz,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_y(nz,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_y    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_p0(nz,nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_p0   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ye_0(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dye_iph(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye_iph   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0_i(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi0_iph(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi0_iph '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_e(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_e      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_ye(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_ye     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de_i(nz,nnu+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de_i      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_ye(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_ye      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (A(nz,nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C(nez+1,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(nez+1,nz+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (R(nez+1,nez+1,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'R         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (A_coef(nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A_coef    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C_coef(nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C_coef    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ncoefa(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ncoefa    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ecoefa(nz,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefa    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

i_ray                      = i_radial
k_ray                      = k_ray_dim * myid + kj_ray

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rho                       = zero
t                         = zero
ye                        = zero

!-----------------------------------------------------------------------
!  Load arrays
!-----------------------------------------------------------------------

rho (:)                   = rho_z (:,kj_ray,ki_ray)
t   (:)                   = t_z   (:,kj_ray,ki_ray)
ye  (:)                   = ye_z  (:,kj_ray,ki_ray)
psi0(:,:)                 = psi0_z(:,:,1,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Functions of the neutrino energy
!-----------------------------------------------------------------------

DO iz = kmin, kmax
  ncoefa (iz,1:nez)       = ncoef * unui(1:nez)**2 * dunui(1:nez)/( agr_z(iz,kj_ray,ki_ray)**3 )
  ecoefa (iz,1:nez)       = ecoef * unui(1:nez)**3 * dunui(1:nez)/( agr_z(iz,kj_ray,ki_ray)**4 )
END DO ! DO iz = kmin, kmax

!-----------------------------------------------------------------------
! Set absorption and emission mean free paths
!-----------------------------------------------------------------------

DO iz = kmin, kmax
  CALL abemset_z( iz, ki_ray, kj_ray, rho(iz), t(iz), ye(iz), agr_z(iz,kj_ray,ki_ray) )
END DO

!-----------------------------------------------------------------------
!  Absorption and emission inverse mean free paths
!-----------------------------------------------------------------------

CALL abemrate_z( kmin, kmax, ki_ray, kj_ray, rho, t, ye, nz )

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

LHS0                     = zero
LHS0_p0                  = zero
LHS0_y                   = zero

RHS0                     = zero
RHS0_y                   = zero
RHS0_p0                  = zero

ye_0                     = zero
dye_iph                  = zero
psi0_i                   = zero
dpsi0_iph                = zero
cf_e                     = zero
cf_ye                    = zero
de_i                     = zero
de                       = zero

A                        = zero
C                        = zero
u                        = zero

R                        = zero
A_coef                   = zero
C_coef                   = zero

ne_1                     = nnugp(1)

n_linear                 = nez+1
n_rhs                    = nez+1
l_leading_A              = nez+1
l_leading_B              = nez+1

!-----------------------------------------------------------------------
!  Store initial values of psi0 and ye
!-----------------------------------------------------------------------

psi0_i                   = psi0
ye_0                     = ye
t_new                    = t
ye_new                   = ye

!-----------------------------------------------------------------------
!  Use multiple of time stwp to ensure equilibration
!-----------------------------------------------------------------------

dtnph_equil              = 10.d0 * dtnph
cdt                      = cvel * dtnph_equil
cdtinv                   = one/cdt

!-----------------------------------------------------------------------
!  Compute Ye transfer constants
!-----------------------------------------------------------------------

DO iz = kmin, kmax
  cf_e (iz,:)            = ecoefa(iz,:)/rho(iz)
  cf_ye(iz,:)            = ncoefa(iz,:) * rmu/rho(iz)
END DO ! iz = kmin,j_equil

!-----------------------------------------------------------------------
!
!                  ||||| BEGIN INTERATION  |||||
!                  |||||      HERE         |||||
!
!-----------------------------------------------------------------------

DO it = 1,iternu

!-----------------------------------------------------------------------
!  Zero coefficient arrays
!-----------------------------------------------------------------------

  LHS0                   = zero          ! transport functions
  LHS0_y                 = zero          ! d(transport functions)/dYe(iz)
  LHS0_p0                = zero          ! d(transport functions)/dpsi0(iz,k)

  RHS0                   = zero          ! source functions
  RHS0_y                 = zero          ! d( source functions )/dYe
  RHS0_p0                = zero          ! d( source functions )/dpsi0

  dye_iph                = zero          ! increments of Ye
  dpsi0_iph              = zero          ! increments of psi0

  R                      = zero          ! recursion matrix

!-----------------------------------------------------------------------
!  Absorption and emission inverse mean free paths
!-----------------------------------------------------------------------

  CALL abemrate_z( kmin, kmax, ki_ray, kj_ray, rho, t, ye_new, nz )

!-----------------------------------------------------------------------
!
!           \\\\\ GET SOURCE FUNCTIONS FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute and store combinations of e-neutrino interaction moments
!   appearing in the p2 equation for psi0.
!
!        tot mfp^-1 : emis + abs - b1  (total inverse mean free path)
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    DO k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  RHS0 e-neutrinos
!  (Right-hand side of the zero-moment transport equation for e-neutrinos)
!-----------------------------------------------------------------------

      RHS0(iz,k)         = emis(iz,k,1) + ( - emis(iz,k,1) - absor(iz,k,1) ) * psi0(iz,k)

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to Ye
!-----------------------------------------------------------------------

      RHS0_y(iz,k)       = emisy(iz,k,1) + ( - emisy(iz,k,1) - absory(iz,k,1) ) * psi0(iz,k)

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nue
!-----------------------------------------------------------------------

      RHS0_p0(iz,k,k)  = RHS0_p0(iz,k,k) - emis(iz,k,1) - absor(iz,k,1)

    END DO ! k = 1,nnugp(1)
  END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!
!             \\\\\ EQUILIBRATION LHS FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    DO k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  LHS0 e-neutrino
!  (Left-hand side of the zero-moment transport equation for e-neutrinos)
!-----------------------------------------------------------------------

      LHS0(iz,k)         = cdtinv * ( psi0(iz,k) - psi0_i(iz,k) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_mh
!-----------------------------------------------------------------------

      LHS0_p0(iz,k)      = cdtinv

    END DO ! k = 1,nnugp(1)
  END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION LHS /////
!
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    LHS0(iz,ne_1+1)      = cdtinv * ( ye_new(iz) - ye(iz) )
    LHS0_y(iz,ne_1+1)    = cdtinv
  END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION RHS /////
!
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    DO k = 1,nnugp(1)
      RHS0(iz,ne_1+1)    = RHS0(iz,ne_1+1) &
&                        - cf_ye(iz,k) * ( emis(iz,k,1)       &
&                        - ( emis(iz,k,1) + absor(iz,k,1) ) * psi0(iz,k) )
    END DO ! k = 1,nnugp(1)
  END DO ! iz = kmin, kmax

  DO iz = kmin, kmax
    DO k = 1,nnugp(1)
      RHS0_p0(iz,ne_1+1,k) = cf_ye(iz,k) * ( emis(iz,k,1) + absor(iz,k,1) )
      RHS0_y(iz,ne_1+1) = RHS0_y(iz,ne_1+1) &
&                        - cf_ye(iz,k) * ( emisy(iz,k,1) - ( emisy(iz,k,1) + absory(iz,k,1) ) * psi0(iz,k) )
    END DO ! k = 1,nnugp(1)
  END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!
!                \\\\\ LOAD COEFFICIENT MATRICES /////
!
!       A      u     + D_m   u      + D_p   u      + C      = 0
!        k,k',iz k',iz      k.iz k,iz-1      k.iz k,iz+1    k,iz
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  A: k-equations and border involving RHS dpsi0
!-----------------------------------------------------------------------

  A                      = - RHS0_p0

!-----------------------------------------------------------------------
!  A: k-equations diagonal elements involving LHS dpsi0(iz-1/2)
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    DO k = 1,ne_1
      A(iz,k,k)          = A(iz,k,k) + LHS0_p0(iz,k)
    END DO ! k = 1,ne_1
  END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!  A: k-equations and border involving LHS and RHS dYw(iz)
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    DO k = 1,ne_1+1
      A(iz,k,ne_1+1)     = A(iz,k,ne_1+1) + LHS0_y(iz,k) - RHS0_y(iz,k)
    END DO ! k = 1,ne_1+1
  END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!  C: k-equations involving LHS and RHS constant terms
!-----------------------------------------------------------------------

  DO k = 1,ne_1+1
    DO iz = kmin, kmax
      C(k,iz)            = LHS0(iz,k) - RHS0(iz,k)
    END DO ! iz = kmin, kmax
  END DO ! k = 1,ne_1+1

!-----------------------------------------------------------------------
!
!            \\\\\ SOLVE FOR THE UP RECURRSION MATRICES /////
!
!                     N
!           u     =  Sum  R       u       + R
!            k,iz     k'=1  k,k',iz  k',iz+1    k,N+1,iz
!
!-----------------------------------------------------------------------

  DO iz = kmin, kmax

!-----------------------------------------------------------------------
!  Load coefficient matrix A_coef with coefficients of the unknowns
!   u(iz,k)
!-----------------------------------------------------------------------

    A_coef(:,:)          = A(iz,:,:)

!-----------------------------------------------------------------------
!  Initialize the right-hand side matrix C_coef
!-----------------------------------------------------------------------

    C_coef               = zero

!-----------------------------------------------------------------------
!  Load right-hand side matrix C_coef with the constants terms in the
!   equations for the unknowns u(iz,k)
!-----------------------------------------------------------------------

    C_coef(:,ne_1+1)     = RHS0(iz,:) - LHS0(iz,:)

!-----------------------------------------------------------------------
!  Solve for the recursion coefficients relating the unknowns u(k,iz)
!   to the unknowns u(k,iz+1)
!-----------------------------------------------------------------------

    CALL dgesv( n_linear, n_rhs, A_coef, l_leading_A, i_piv, C_coef, &
&    l_leading_B, info )

!-----------------------------------------------------------------------
!  Load recursion matrix
!-----------------------------------------------------------------------

    R(:,:,iz)            = C_coef(:,:)

  END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!
!             \\\\\ COMPUTE INCREMENTS OF VARIABLES /////
!
!            i            i-1    i          i-1
!           u      = SUM  (R      u      ) + S
!            k,iz      kp  k,kp,iz kp,iz+1     k,iz
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize the solution vector
!-----------------------------------------------------------------------

  u                      = zero

!-----------------------------------------------------------------------
!  Recursively compute the solution vector
!-----------------------------------------------------------------------

  DO iz = kmax,kmin,-1
    DO k = 1,ne_1+1
      DO kp = 1,ne_1
        u(k,iz)          = u(k,iz) + R(k,kp,iz) * u(kp,iz+1)
      END DO ! kp = 1,ne_1
      u(k,iz)            = u(k,iz) + R(k,ne_1+1,iz)
    END DO ! k = 1,ne_1+1
  END DO ! iz = kmax,kmin,-1

!-----------------------------------------------------------------------
!
!      \\\\\ TRANSFER INVCREMENTS AND TEST FOR CONVERGENCE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer increments
!-----------------------------------------------------------------------

  DO iz = kmin, kmax

    IF ( nnugp(1) /= 0 ) THEN
      DO k = 1,nnugp(1)
        dpsi0_iph(iz,k)  = u(k,iz)
      END DO ! k = 1,nnugp(1)
    END IF ! nnugp(1) /= 0

    dye_iph(iz)          = u(ne_1+1,iz)

  END DO ! iz = kmin, kmax
 
!-----------------------------------------------------------------------
!        Increment variables
!
!  psi0(iz,k) is allowed to exceed unity during iterations in order for
!   the iterations to converge.
!  psi0(iz,k) is not allowed to fall below zero to avoid unphysical solutions
!
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    DO k = 1,nnugp(1)
      psi0(iz,k)     = DMAX1( psi0(iz,k) + dpsi0_iph(iz,k), psi0amin )
    END DO ! k = 1,nnugp(n)
  END DO ! iz = kmin, kmax

  DO iz = kmin, kmax
    ye_new(iz)           = DMAX1( ye_new(iz) + dye_iph(iz), zero )
  END DO

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

  endit                  = .true.
  jiterp0                = 0
  kiterp0                = 0
  niterp0                = 0
  jiterye                = 0

loop_checkpsi: DO iz = kmin, kmax
    DO k = 1,nnugp(1)
      IF ( DABS( dpsi0_iph(iz,k)/DMAX1( psi0(iz,k), tolpsimin ) ) >= tolnupsi ) THEN
        jiterp0          = iz
        kiterp0          = k
        niterp0          = 1
        endit            = .false.
        EXIT loop_checkpsi
      END IF ! dpsi0_iph(iz,k)/DMAX1( psi0(iz,k) > tolnupsi
    END DO ! k = 1,nnugp(1)
  END DO loop_checkpsi

loop_checkye:  DO iz = kmin, kmax
    IF ( DABS( dye_iph(iz)/DMAX1( ye_new(iz), 0.1d0 ) ) >= tolnuye ) THEN
      jiterye            = iz
      endit              = .false.
      EXIT loop_checkye
    END IF ! DABS( dye_iph(iz)/DMAX1( ye(iz), 0.1 ) ) >= tolnuye
  END DO loop_checkye

!-----------------------------------------------------------------------
!  WRITE diagnostics if it > 5
!-----------------------------------------------------------------------

  IF ( .not. endit  .and.  it >= 5 ) THEN
    IF ( jiterp0 == 0 ) THEN
      WRITE (nlog,1501) it, myid, ki_ray, kj_ray, jiterye, dye_iph(jiterye),                 &
&      dye_iph(jiterye)/DMAX1( ye_new(jiterye), 0.1d0 ), jiterye, ye_new(jiterye)
    ELSE IF ( jiterye == 0 ) THEN
      WRITE (nlog,1503) it, myid, ki_ray, kj_ray, jiterp0, kiterp0, niterp0,                 &
&      dpsi0_iph(jiterp0,kiterp0),                                                           &
&      dpsi0_iph(jiterp0,kiterp0)/DMAX1( psi0(jiterp0,kiterp0), tolpsimin ),                 &
&      jiterp0, kiterp0, niterp0, psi0(jiterp0,kiterp0)
    ELSE
      WRITE (nlog,1505) it, myid, ki_ray, kj_ray, jiterye, dye_iph(jiterye),dye_iph(jiterye) &
&      /DMAX1( ye_new(jiterye), 0.1d0 ), jiterye, ye_new(jiterye), jiterp0, kiterp0,         &
&      niterp0, dpsi0_iph(jiterp0,kiterp0),                                                  &
&      dpsi0_iph(jiterp0,kiterp0)/DMAX1( psi0(jiterp0,kiterp0), tolpsimin ),                 &
&      jiterp0, kiterp0, niterp0, psi0(jiterp0,kiterp0)
    END IF ! jiterp0 == 0
  END IF ! .not. endit

!-----------------------------------------------------------------------
!
!                    \\\\\ END ITERATION /////
!
!-----------------------------------------------------------------------

  IF ( endit  .or.  iternu == 1 ) EXIT
  IF ( it == iternu ) THEN
    WRITE (nlog,7001) i_ray, k_ray, myid
    WRITE (nprint,7001) i_ray, k_ray, myid
    STOP
    EXIT
  END IF

END DO ! iteration loop

!-----------------------------------------------------------------------
!
!                 \\\\\ UPDATE THE TEMPERATURE /////
!
!-----------------------------------------------------------------------

IF ( consrve == 'e' ) THEN

!-----------------------------------------------------------------------
!  Compute energy transferred to neutrinos
!-----------------------------------------------------------------------

  de                     = zero

  DO k = 1,nnugp(1)
    DO iz = kmin, kmax
      de(iz)             = de(iz) + cf_e(iz,k) * stwt(1) * cdt            &
&                        * ( RHS0(iz,k) + RHS0_y(iz,k) * dye_iph(iz) )
      DO kp = 1,nnugp(1)
        de(iz)           = de(iz) + cf_e(iz,k) * stwt(1) * cdt * RHS0_p0(iz,k,     kp) * dpsi0_iph(iz,kp)
      END DO ! kp = 1,nnugp(1)
    END DO ! iz = kmin, kmax
  END DO ! k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  Compute the new temperature given by energy conservation
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    CALL eqstt_z( 2, iz, ki_ray, kj_ray, rho(iz), t(iz), ye_0(iz), e, dedd, dedt, dedy )
    CALL e_degen(  rho(iz), t_new(iz), ye_0(iz), e_eta )
    IF ( e_eta > degen ) THEN
      t_new(iz)          = t_new(iz) + ( - de(iz) - dedy * ( ye_new(iz) - ye_0(iz) ) )/dedt
    ELSE
      CALL tgvndeye_z( iz, ki_ray, kj_ray, rho(iz), e - de(iz), ye_new(iz), t(iz), t_new(iz) )
    END IF ! e_eta > degen
  END DO ! iz = kmin, kmax

ELSE ! consrve = 's'

!-----------------------------------------------------------------------
!  Preserve entropy
!-----------------------------------------------------------------------

  DO iz = kmin, kmax
    CALL eqstt_z( 3, iz, ki_ray, kj_ray, rho(iz), t(iz), ye(iz), s_i, dsdd, dsdt, dsdy )
    CALL tgvndsye_z( iz, ki_ray, kj_ray, rho(iz), s_i, ye_new(iz), t, t_new(iz) )
  END DO ! iz = kmin, kmax

END IF ! consrve = 'e'

!-----------------------------------------------------------------------
!  Conserve leptons
!-----------------------------------------------------------------------

DO iz = kmin, kmax
  d_ye(iz)                = - SUM( cf_ye(iz,:) * ( psi0(iz,:) - psi0_i(iz,:) ) )
  ye_new(iz)              = ye_0(iz) + d_ye(iz)
END DO ! iz = kmin, kmax

!-----------------------------------------------------------------------
!  Restore variable for transfer back to AZIMUTHALRAY_MODULE
!-----------------------------------------------------------------------

t_z   (kmin:kmax,kj_ray,ki_ray)     = t_new (kmin:kmax)
ye_z  (kmin:kmax,kj_ray,ki_ray)     = ye_new(kmin:kmax)
psi0_z(kmin:kmax,:,1,kj_ray,ki_ray) = psi0  (kmin:kmax,:)

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (rho, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (t_new, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_new     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye_new, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_new    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (LHS0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (LHS0_p0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (LHS0_y, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_y    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (RHS0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (RHS0_y, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_y    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (RHS0_p0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_p0   '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ye_0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dye_iph, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye_iph   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi0_i, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_i    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dpsi0_iph, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi0_iph '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (cf_e, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_e      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (cf_ye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_ye     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (de_i, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de_i      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (de, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (d_ye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_ye      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (A, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (C, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (R, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'R         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (A_coef, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A_coef    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (C_coef, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C_coef    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ncoefa, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ncoefa    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ecoefa, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefa    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE equilibrate_z
