SUBROUTINE equilibrate_y( jmin, jmax, nx, ny, nez, nnu, ji_ray, jk_ray, &
& j_ray_dim, ik_ray_dim, i_radial, rho_y, t_y, ye_y, psi0_y, rhobar,    &
& dtnph, rho_equilibrate, agr_y )
!-----------------------------------------------------------------------
!
!    File:         equilibrate_y
!    Module:       equilibrate_y
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
!  abemset_y       : sets the neutrino absorption and emission rates at cube corners
!  abemrate_y      : computes neutrino absorption and emission rates
!  dgesv           : Solve for the recursion coefficients relating the unknowns to the knowns
!  eqstt_y         : computes equation of state variables
!  tgvndeye_y      : computes the temperature given rho, ei, and ye
!  tgvndsye_y      : computes the temperature given rho, s, and ye
!
!    Input arguments:
!  jmin            : minimum angular zone to equilibrate matter and neutrinos
!  jmax            : maximum angular zone to equilibrate matter and neutrinos
!  nx              : x-array extent
!  ny              : y-array extent
!  ji_ray          : x (radial) index of a specific angular ray
!  jk_ray          : z (azimuthal) index of a specific angular ray
!  j_ray_dim       : the number of radial zones on a processor after swapping with y
!  ik_ray_dim      : the number of z-zones on a processor before swapping with z
!  i_radial        : the unshifted radial zone (angular ray) corresponding to ji_ray
!  nez             : neutrino energy array extent
!  nnu             : neutrino flavor array extent
!  rho_y           : unshifted density [g cm^{-3}]
!  t_y             : unshifted temperature (K)
!  ye_y            : unshifted electron fraction
!  psi0_y          : unshifted zero angular moments of the neutrino occupation number
!  rhobar          : mean density of an angular ray [g cm^{-3}]
!  dtnph           : time step
!  rho_equilibrate : density above whicn neutrinos equilibrated with matter
!  agr_y           : lapse function
!
!    Output arguments:
!  psi0_y          : unshifted zero angular moments of the neutrino occupation number
!  t_y             : unshifted temperature (K)
!  ye_y            : unshifted electron fraction
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

USE abem_y_module, ONLY : emis, emisy, absor, absory
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

INTEGER, INTENT(in)              :: jmin            ! minimum angular zone index
INTEGER, INTENT(in)              :: jmax            ! maximum angular zone index
INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: ny              ! y-array extent
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ji_ray          ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray          ! z (azimuthal) index of a specific angular ray
INTEGER, INTENT(in)              :: j_ray_dim       ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z (azimuthal) zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial        ! the unshifted radial zone corresponding to ji_ray, jk_ray

REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim) :: rho_y            ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                      :: rhobar           ! mean density of an angular ray (g cm^{-3})
REAL(KIND=double), INTENT(in)    :: dtnph           ! time step
REAL(KIND=double), INTENT(in)    :: rho_equilibrate ! mean density of an angular ray (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim) :: agr_y            ! lapse function

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim) :: t_y            ! temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim) :: ye_y           ! electron fraction
REAL(KIND=double), INTENT(inout), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim) :: psi0_y ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name
CHARACTER (len=1), PARAMETER     :: consrve = 's'

LOGICAL                          :: endit

INTEGER                          :: j               ! radial zone index
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
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_p0      ! coefficients of dpsi0(j) terms left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_y       ! coefficients of dYe(j) terms left-hand side of p0 equation
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

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: A            ! coefficient matrix for the j terms in j-equation, all j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: C            ! constant terms in j-equation, all j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: u            ! solution vector of the linear equations

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: R            ! recursion matrix
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: A_coef       ! Coefficient matrix of the variables to be solved
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: C_coef       ! Coefficient matrix of quantities on the right-hand side

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: ncoefa       ! neutrino number density coefficient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: ecoefa       ! neutrino energy density coefficient

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in equilibrate_y')
 1501 FORMAT (' Iteration attempt',i4,' myid=',i4,' ji_ray=',i4,' jk_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=',es11.3,' ye(',i4,')=',es11.3)
 1503 FORMAT (' Iteration attempt',i4,' myid=',i4,' ji_ray=',i4,' jk_ray=',i4, &
&' dpsi0_iph(',i4,',',i4,',',i4,')=', es11.3,' D_psi=',es11.3,                 &
&' psi0(',i4,',',i4,',',i4,')=', es11.3)
 1505 FORMAT (' Iteration attempt',i4,' myid=',i4,' ji_ray=',i4,' jk_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=', es11.3,' ye(',i4,')=',es11.3,           &
& ' dpsi0_iph(',i4,',',i4,',',i4,')=',es11.3,' D_psi=',es11.3,                 &
& ' psi0(',i4,',',i4,',',i4,')=', es11.3)
 2001 FORMAT (' Deallocation problem for array ',a10,' in equilibrate_y')
 7001 FORMAT (' Convergence failure in subroutine equilibrate_y, i_ray=',i4,   &
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

ALLOCATE (rho(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_new(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_new     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_new(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_new    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (LHS0(ny,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_p0(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_y(ny,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_y    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (RHS0(ny,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_y(ny,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_y    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_p0(ny,nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_p0   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ye_0(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dye_iph(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye_iph   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0_i(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi0_iph(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi0_iph '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_e(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_e      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_ye(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_ye     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de_i(ny,nnu+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de_i      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_ye(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_ye      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (A(ny,nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C(nez+1,ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(nez+1,ny+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (R(nez+1,nez+1,ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'R         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (A_coef(nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A_coef    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C_coef(nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C_coef    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ncoefa(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ncoefa    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ecoefa(ny,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefa    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

i_ray                    = i_radial
k_ray                    = ik_ray_dim * myid + jk_ray

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rho                       = zero
t                         = zero
ye                        = zero

!-----------------------------------------------------------------------
!  Load arrays
!-----------------------------------------------------------------------

rho (:)                   = rho_y (:,ji_ray,jk_ray)
t   (:)                   = t_y   (:,ji_ray,jk_ray)
ye  (:)                   = ye_y  (:,ji_ray,jk_ray)
psi0(:,:)                 = psi0_y(:,:,1,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Functions of the neutrino energy
!-----------------------------------------------------------------------

DO j = jmin,jmax
  ncoefa (j,1:nez)        = ncoef * unui(1:nez)**2 * dunui(1:nez)/( agr_y(j,ji_ray,jk_ray)**3 )
  ecoefa (j,1:nez)        = ecoef * unui(1:nez)**3 * dunui(1:nez)/( agr_y(j,ji_ray,jk_ray)**4 )
END DO ! DO j = jmin,jmax

!-----------------------------------------------------------------------
! Set absorption and emission mean free paths
!-----------------------------------------------------------------------

DO j = jmin,jmax
  CALL abemset_y( j, ji_ray, jk_ray, rho(j), t(j), ye(j), agr_y(j,ji_ray,jk_ray) )
END DO

!-----------------------------------------------------------------------
!  Absorption and emission inverse mean free paths
!-----------------------------------------------------------------------

CALL abemrate_y( jmin, jmax, ji_ray, jk_ray, rho, t, ye, ny )

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

DO j = jmin,jmax
  cf_e (j,:)             = ecoefa(j,:)/rho(j)
  cf_ye(j,:)             = ncoefa(j,:) * rmu/rho(j)
END DO ! j = jmin,j_equil

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
  LHS0_y                 = zero          ! d(transport functions)/dYe(j)
  LHS0_p0                = zero          ! d(transport functions)/dpsi0(j,k)

  RHS0                   = zero          ! source functions
  RHS0_y                 = zero          ! d( source functions )/dYe
  RHS0_p0                = zero          ! d( source functions )/dpsi0

  dye_iph                = zero          ! increments of Ye
  dpsi0_iph              = zero          ! increments of psi0

  R                      = zero          ! recursion matrix

!-----------------------------------------------------------------------
!  Absorption and emission inverse mean free paths
!-----------------------------------------------------------------------

  CALL abemrate_y( jmin, jmax, ji_ray, jk_ray, rho, t, ye_new, ny )

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

  DO j = jmin,jmax
    DO k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  RHS0 e-neutrinos
!  (Right-hand side of the zero-moment transport equation for e-neutrinos)
!-----------------------------------------------------------------------

      RHS0(j,k)          = emis(j,k,1) + ( - emis(j,k,1) - absor(j,k,1) ) * psi0(j,k)

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to Ye
!-----------------------------------------------------------------------

      RHS0_y(j,k)        = emisy(j,k,1) + ( - emisy(j,k,1) - absory(j,k,1) ) * psi0(j,k)

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nue
!-----------------------------------------------------------------------

      RHS0_p0(j,k,k)     = RHS0_p0(j,k,k) - emis(j,k,1) - absor(j,k,1)

    END DO ! k = 1,nnugp(1)
  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!
!             \\\\\ EQUILIBRATION LHS FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    DO k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  LHS0 e-neutrino
!  (Left-hand side of the zero-moment transport equation for e-neutrinos)
!-----------------------------------------------------------------------

      LHS0(j,k)          = cdtinv * ( psi0(j,k) - psi0_i(j,k) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_mh
!-----------------------------------------------------------------------

      LHS0_p0(j,k)       = cdtinv

    END DO ! k = 1,nnugp(1)
  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION LHS /////
!
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    LHS0(j,ne_1+1)       = cdtinv * ( ye_new(j) - ye(j) )
    LHS0_y(j,ne_1+1)     = cdtinv
  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION RHS /////
!
!-----------------------------------------------------------------------

  DO j = jmin,jmax
     DO k = 1,nnugp(1)
       RHS0(j,ne_1+1)    = RHS0(j,ne_1+1) &
&                        - cf_ye(j,k) * ( emis(j,k,1)       &
&                        - ( emis(j,k,1) + absor(j,k,1) ) * psi0(j,k) )
    END DO ! k = 1,nnugp(1)
  END DO ! j = jmin,jmax

  DO j = jmin,jmax
    DO k = 1,nnugp(1)
      RHS0_p0(j,ne_1+1,k) = cf_ye(j,k) * ( emis(j,k,1) + absor(j,k,1) )
      RHS0_y(j,ne_1+1)    = RHS0_y(j,ne_1+1) &
&                         - cf_ye(j,k) * ( emisy(j,k,1) - ( emisy(j,k,1) + absory(j,k,1) ) * psi0(j,k) )
    END DO ! k = 1,nnugp(1)
  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!
!                \\\\\ LOAD COEFFICIENT MATRICES /////
!
!       A      u     + D_m   u      + D_p   u      + C      = 0
!        k,k',j k',j      k.j k,j-1      k.j k,j+1    k,j
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  A: k-equations and border involving RHS dpsi0
!-----------------------------------------------------------------------

  A                      = - RHS0_p0

!-----------------------------------------------------------------------
!  A: k-equations diagonal elements involving LHS dpsi0(j-1/2)
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    DO k = 1,ne_1
      A(j,k,k)           = A(j,k,k) + LHS0_p0(j,k)
    END DO ! k = 1,ne_1
  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!  A: k-equations and border involving LHS and RHS dYw(j)
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    DO k = 1,ne_1+1
      A(j,k,ne_1+1)      = A(j,k,ne_1+1) + LHS0_y(j,k) - RHS0_y(j,k)
    END DO ! k = 1,ne_1+1
  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!  C: k-equations involving LHS and RHS constant terms
!-----------------------------------------------------------------------

  DO k = 1,ne_1+1
    DO j = jmin,jmax
      C(k,j)             = LHS0(j,k) - RHS0(j,k)
    END DO ! j = jmin,jmax
  END DO ! k = 1,ne_1+1

!-----------------------------------------------------------------------
!
!            \\\\\ SOLVE FOR THE UP RECURRSION MATRICES /////
!
!                     N
!           u     =  Sum  R       u       + R
!            k,j     k'=1  k,k',j  k',j+1    k,N+1,j
!
!-----------------------------------------------------------------------

  DO j = jmin,jmax

!-----------------------------------------------------------------------
!  Load coefficient matrix A_coef with coefficients of the unknowns
!   u(j,k)
!-----------------------------------------------------------------------

    A_coef(:,:)          = A(j,:,:)

!-----------------------------------------------------------------------
!  Initialize the right-hand side matrix C_coef
!-----------------------------------------------------------------------

    C_coef               = zero

!-----------------------------------------------------------------------
!  Load right-hand side matrix C_coef with the constants terms in the
!   equations for the unknowns u(j,k)
!-----------------------------------------------------------------------

    C_coef(:,ne_1+1)     = RHS0(j,:) - LHS0(j,:)

!-----------------------------------------------------------------------
!  Solve for the recursion coefficients relating the unknowns u(k,j)
!   to the unknowns u(k,j+1)
!-----------------------------------------------------------------------

    CALL dgesv( n_linear, n_rhs, A_coef, l_leading_A, i_piv, C_coef, &
&    l_leading_B, info )

!-----------------------------------------------------------------------
!  Load recursion matrix
!-----------------------------------------------------------------------

    R(:,:,j)             = C_coef(:,:)

  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!
!             \\\\\ COMPUTE INCREMENTS OF VARIABLES /////
!
!            i            i-1    i          i-1
!           u      = SUM  (R      u      ) + S
!            k,j      kp  k,kp,j kp,j+1     k,j
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize the solution vector
!-----------------------------------------------------------------------

  u                      = zero

!-----------------------------------------------------------------------
!  Recursively compute the solution vector
!-----------------------------------------------------------------------

  DO j = jmax,jmin,-1
    DO k = 1,ne_1+1
      DO kp = 1,ne_1
        u(k,j)           = u(k,j) + R(k,kp,j) * u(kp,j+1)
      END DO ! kp = 1,ne_1
      u(k,j)             = u(k,j) + R(k,ne_1+1,j)
    END DO ! k = 1,ne_1+1
  END DO ! j = jmax,jmin,-1

!-----------------------------------------------------------------------
!
!      \\\\\ TRANSFER INVCREMENTS AND TEST FOR CONVERGENCE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer increments
!-----------------------------------------------------------------------

  DO j = jmin,jmax

    IF ( nnugp(1) /= 0 ) THEN
      DO k = 1,nnugp(1)
        dpsi0_iph(j,k)   = u(k,j)
      END DO ! k = 1,nnugp(1)
    END IF ! nnugp(1) /= 0

    dye_iph(j)           = u(ne_1+1,j)

  END DO ! j = jmin,jmax
 
!-----------------------------------------------------------------------
!        Increment variables
!
!  psi0(j,k) is allowed to exceed unity during iterations in order for
!   the iterations to converge.
!  psi0(j,k) is not allowed to fall below zero to avoid unphysical solutions
!
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    DO k = 1,nnugp(1)
      psi0(j,k)          = DMAX1( psi0(j,k) + dpsi0_iph(j,k), psi0amin )
    END DO ! k = 1,nnugp(n)
  END DO ! j = jmin,jmax

  DO j = jmin,jmax
    ye_new(j)            = DMAX1( ye_new(j) + dye_iph(j), zero )
  END DO

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

  endit                  = .true.
  jiterp0                = 0
  kiterp0                = 0
  niterp0                = 0
  jiterye                = 0

loop_checkpsi: DO j = jmin,jmax
    DO k = 1,nnugp(1)
      IF ( DABS( dpsi0_iph(j,k)/DMAX1( psi0(j,k), tolpsimin ) ) >= tolnupsi ) THEN
        jiterp0          = j
        kiterp0          = k
        niterp0          = 1
        endit            = .false.
        EXIT loop_checkpsi
      END IF ! dpsi0_iph(j,k)/DMAX1( psi0(j,k) > tolnupsi
    END DO ! k = 1,nnugp(1)
  END DO loop_checkpsi

loop_checkye:  DO j = jmin,jmax
    IF ( DABS( dye_iph(j)/DMAX1( ye_new(j), 0.1d0 ) ) >= tolnuye ) THEN
      jiterye            = j
      endit              = .false.
      EXIT loop_checkye
    END IF ! DABS( dye_iph(j)/DMAX1( ye(j), 0.1 ) ) >= tolnuye
  END DO loop_checkye

!-----------------------------------------------------------------------
!  WRITE diagnostics if it > 5
!-----------------------------------------------------------------------

  IF ( .not. endit  .and.  it >= 5 ) THEN
    IF ( jiterp0 == 0 ) THEN
      WRITE (nlog,1501) it, myid, ji_ray, jk_ray, jiterye, dye_iph(jiterye),                 &
&      dye_iph(jiterye)/DMAX1( ye_new(jiterye), 0.1d0 ), jiterye, ye_new(jiterye)
    ELSE IF ( jiterye == 0 ) THEN
      WRITE (nlog,1503) it, myid, ji_ray, jk_ray, jiterp0, kiterp0, niterp0,                 &
&      dpsi0_iph(jiterp0,kiterp0),                                                           &
&      dpsi0_iph(jiterp0,kiterp0)/DMAX1( psi0(jiterp0,kiterp0), tolpsimin ),                 &
&      jiterp0, kiterp0, niterp0, psi0(jiterp0,kiterp0)
    ELSE
      WRITE (nlog,1505) it, myid, ji_ray, jk_ray, jiterye, dye_iph(jiterye),dye_iph(jiterye) &
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
    DO j = jmin,jmax
      de(j)              = de(j) + cf_e(j,k) * stwt(1) * cdt            &
&                        * ( RHS0(j,k) + RHS0_y(j,k) * dye_iph(j) )
      DO kp = 1,nnugp(1)
        de(j)            = de(j) + cf_e(j,k) * stwt(1) * cdt * RHS0_p0(j,k,     kp) * dpsi0_iph(j,kp)
      END DO ! kp = 1,nnugp(1)
    END DO ! j = jmin,jmax
  END DO ! k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  Compute the new temperature given by energy conservation
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    CALL eqstt_y( 2, j, ji_ray, jk_ray, rho(j), t(j), ye_0(j), e, dedd, dedt, dedy )
    CALL e_degen(  rho(j), t(j), ye_0(j), e_eta )
    IF ( e_eta > degen ) THEN
      t_new(j)           = t_new(j) + ( - de(j) - dedy * ( ye_new(j) - ye_0(j) ) )/dedt
    ELSE
      CALL tgvndeye_y( j, ji_ray, jk_ray, rho(j), e - de(j), ye_new(j), t(j), t_new(j) )
    END IF ! e_eta > degen
  END DO ! j = jmin,jmax

ELSE ! consrve = 's'

!-----------------------------------------------------------------------
!  Preserve entropy
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    CALL eqstt_y( 3, j, ji_ray, jk_ray, rho(j), t(j), ye(j), s_i, dsdd, dsdt, dsdy )
    CALL tgvndsye_y( j, ji_ray, jk_ray, rho(j), s_i, ye_new(j), t(j), t_new(j) )
  END DO ! j = jmin,jmax

END IF ! consrve = 'e'

!-----------------------------------------------------------------------
!  Conserve leptons
!-----------------------------------------------------------------------

DO j = jmin, jmax
  d_ye(j)                = - SUM( cf_ye(j,:) * ( psi0(j,:) - psi0_i(j,:) ) )
  ye_new(j)              = ye_0(j) + d_ye(j)
END DO ! j = jmin, jmax

!-----------------------------------------------------------------------
!  Load arrays for transferring back to ANGULAR_RAY_MODULE
!-----------------------------------------------------------------------

t_y   (jmin:jmax,ji_ray,jk_ray)     = t_new (jmin:jmax)
ye_y  (jmin:jmax,ji_ray,jk_ray)     = ye_new(jmin:jmax)
psi0_y(jmin:jmax,:,1,ji_ray,jk_ray) = psi0  (jmin:jmax,:)

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
END SUBROUTINE equilibrate_y
