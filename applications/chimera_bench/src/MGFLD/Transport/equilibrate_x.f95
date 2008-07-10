SUBROUTINE equilibrate_x( imin, imax, nx, nez, nnu, ij_ray, ik_ray, &
& ij_ray_dim, ik_ray_dim, rho_c, t_c, ye_c, psi0_c, dtnph, rho_equilibrate )
!-----------------------------------------------------------------------
!
!    File:         equilibrate_x
!    Module:       equilibrate_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/5/05
!
!    Purpose:
!      To compute iterate the change in Ye and psi0 due to absorption,
!       emission, inelastic scattering, pair production, and transport.g
!
!    Subprograms called:
!  eqstt_x   : interpolates quantities in the local EOS table
!  dtau_trns : computes transport time step
!  ddc_dpsi  : computes derivatives of the diffusion coefficient
!  flux      : computes the neutrino flux
!  ludcmp    : performs an LU decomposition of a matrix
!  lubkcsb   : inverts a matrix given an LU decomposition
!  mprove    : improves the solution of the matrix inversion
!  psi_bd    : computes outer boundary value of psi0
!  psia_cal  : computes psi1
!
!    Input arguments:
!
!  imin             : unshifted inner radial zone number
!  imax             : unshifted outer radial zone number
!  nx               : x-array extent
!  nez              : neutrino energy array extent
!  nnu              : neutrino flavor array extent
!  ij_ray           : j-index of a radial ray
!  ik_ray           : k-index of a radial ray
!  ij_ray_dim       : number of y-zones on a processor before swapping
!  ik_ray_dim       : number of z-zones on a processor before swapping
!  rho_c            : density [g cm^{-3}]
!  t_c              : temperature [K]
!  ye_c             : electron fraction
!  psi0_c           : zeroth angular moment of the neutrino distribution
!                      function
!  dtnoh            : time step
!  rho_equilibrate  : density above which to equilibrate neutrinos with matter
!
!    Output arguments:
!  t_c              : temperature [K]
!  ye_c             : electron fraction
!  psi0_c           : zeroth angular moment of the neutrino distribution
!
!    Input arguments (common):
!
!  nnugp(n)         : number of energy zones for neutrinos of type n
!  unu(j,k)         : midpoint energy of energy zone k at radial zone j
!  dunu(j,k)        : energy width of energy zone k at radial zone j
!  tolnupsi         : abs(dpsi_iph)/( psi0 + tolpsimin ) < tolnupsi for
!  tolpsimin        :  for convergence
!                   : proper time step at j-1/2
!  ncoefa(j,k)      : 4.*pi/((h*!)**3)*w**2*dw
!  ecoefa(j,k)      : 4.*pi/((h*!)**3)*w**3*dw
!  vol(j)           : volume enclosed by radial zone j and j-1
!  cf_T_psi(j,k)    : coefficient relating dT to dpsi
!  aesvt(nx,12)     : derivatives with respect to t of equation of state dependent variables
!
!    Output arguments (common):
!
!  dt_scat(j,n)
!                   : change in t due to inelastic scattering and pair production
!  dpsi_scat(j,k,n) : change in psi0 due to inelastic scattering and pair production
!  psi0(j,k,n)      : updated zeroth angular moment of the neutrino distribution function
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  abem_module, edit_module, evh1_global, it_tol_module,
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module, parallel_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc_y
USE numerical_module, ONLY : zero, one, ncoef, ecoef
USE physcnst_module, ONLY : rmu, cvel

USE abem_module, ONLY : emis, emisy, absor, absory
USE edit_module, ONLY : nprint, nlog
USE evh1_global, ONLY : degen
USE it_tol_module, ONLY : iternu, tolpsimin, tolnupsi, tolnuye
USE mdl_cnfg_module, ONLY : dmrst, u_vel=>u
USE nu_dist_module, ONLY : ncoefa, ecoefa, unu, dunu, stwt
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY : iaefnp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin            ! unshifted minimum radial zone index
INTEGER, INTENT(in)              :: imax            ! unshifted maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nx              ! radial array extent
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rho_c  ! density (g cm^{-3})
REAL(KIND=double), INTENT(in)    :: dtnph           ! time step
REAL(KIND=double), INTENT(in)    :: rho_equilibrate ! mean density of an angular ray (g cm^{-3})

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: t_c    ! initial temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: ye_c   ! electron fraction
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_c ! electron fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name
CHARACTER (len=1), PARAMETER     :: consrve = 's'

LOGICAL                          :: endit

INTEGER                          :: jr_min          ! minimum shifted radial zone index
INTEGER                          :: jr_max          ! maximum shifted radial zone index
INTEGER                          :: j_equil         ! maximum shifted radial zone index out which to equilibrate

INTEGER                          :: j               ! radial zone index
INTEGER                          :: k               ! neutrino energy index
INTEGER                          :: kp              ! neutrino energy index
INTEGER                          :: it              ! iteration index
INTEGER                          :: jiterp0         ! radial zone for which psi0 iteration fails
INTEGER                          :: kiterp0         ! neutrino energy zone for which psi0 iteration fails
INTEGER                          :: niterp0         ! neutrino flavor for which psi0 iteration fails
INTEGER                          :: jiterye         ! radial zone number for which ye iteration fails

INTEGER                          :: j_ray           ! polar index of the radial ray
INTEGER                          :: k_ray           ! azimuthal index of the radial ray

INTEGER                          :: ne_1            ! nnugp(1)

INTEGER                          :: istat           ! allocation status

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

REAL(KIND=double)                :: e_eta           ! (electron Fermi energy)/KT

REAL(KIND=double), PARAMETER     :: psi0amax = 1.d0 ! maximum value of psi0
REAL(KIND=double), PARAMETER     :: psi0amin = 0.d0 ! minimum value of psi0

REAL(KIND=double), DIMENSION(nx)                 :: rho          ! 1-D shifted array of density [g cm^{-3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t            ! 1-D shifted array of temperature [K]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t_new        ! 1-D shifted array of updated temperature [K]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye           ! 1-D shifted array of electron fraction
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye_new       ! 1-D shifted array of updated electron fraction
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: psi0         ! 2-D shifted array of psi0

REAL(KIND=double)                :: s_i             ! shifted array of entropy
REAL(KIND=double)                :: dsdd            ! d(s_i)/d(rho)
REAL(KIND=double)                :: dsdt            ! d(s_i)/d(T)
REAL(KIND=double)                :: dsdy            ! d(s_i)/d(Ye)

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

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in equilibrate_x')
 1501 FORMAT (' Iteration attempt',i4,' myid=',i4,' ij_ray=',i4,' ik_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=',es11.3,' ye(',i4,')=',es11.3)
 1503 FORMAT (' Iteration attempt',i4,' myid=',i4,' ij_ray=',i4,' ik_ray=',i4, &
&' dpsi0_iph(',i4,',',i4,',',i4,')=', es11.3,' D_psi=',es11.3,                 &
&' psi0(',i4,',',i4,',',i4,')=', es11.3)
 1505 FORMAT (' Iteration attempt',i4,' myid=',i4,' ij_ray=',i4,' ik_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=', es11.3,' ye(',i4,')=',es11.3,           &
& ' dpsi0_iph(',i4,',',i4,',',i4,')=',es11.3,' D_psi=',es11.3,                 &
& ' psi0(',i4,',',i4,',',i4,')=', es11.3)
 2001 FORMAT (' Deallocation problem for array ',a10,' in equilibrate_x')
 7001 FORMAT (' Convergence failure in subroutine equilibrate_x, j_ray=',i4,   &
& ' k_ray=',i4,' myid=',i4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if iaefnp = 0 or nnugp(1) = 0
!-----------------------------------------------------------------------

IF ( iaefnp == 0  .or.  nnugp(1) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Determine j_equil, the maximum shifted radial zone index out which
!   to equilibrate
!-----------------------------------------------------------------------

jr_min                   = imin + 1
jr_max                   = imax + 1

rho(jr_min:jr_max)       = rho_c(imin:imax,ij_ray,ik_ray)
IF ( rho(2) < rho_equilibrate ) THEN
  j_equil                = 2
ELSE
  j_equil                = MAXLOC( rho, DIM = 1, MASK = rho < rho_equilibrate )
END IF ! rho(2) < rho_equilibrate

!-----------------------------------------------------------------------
!  Return if j_equil < 3
!-----------------------------------------------------------------------

IF ( j_equil < 3 ) RETURN

!-----------------------------------------------------------------------
!
!                   \\\\\ ALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (t(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_new(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_new     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_new(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_new    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (LHS0(nx,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_p0(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_y(nx,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_y    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (RHS0(nx,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_y(nx,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_y    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_p0(nx,nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_p0   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ye_0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dye_iph(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye_iph   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0_i(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi0_iph(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi0_iph '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_e(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_e      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_ye(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_ye     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de_i(nx,nnu+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de_i      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_ye(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_ye      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (A(nx,nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C(nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (R(nez+1,nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'R         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (A_coef(nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A_coef    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C_coef(nez+1,nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C_coef    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_ray                    = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
k_ray                    = ( myid/n_proc_y ) * ik_ray_dim + ik_ray

!-----------------------------------------------------------------------
!  Load t, ye, and psi0
!-----------------------------------------------------------------------

t   (jr_min:jr_max)      = t_c   (imin:imax,ij_ray,ik_ray)
ye  (jr_min:jr_max)      = ye_c  (imin:imax,ij_ray,ik_ray)
psi0(jr_min:jr_max,:)    = psi0_c(imin:imax,:,1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Functions of the neutrino energy
!-----------------------------------------------------------------------

ncoefa (jr_min:jr_max,:) = ncoef * unu (jr_min:jr_max,:)**2 * dunu (jr_min:jr_max,:)
ecoefa (jr_min:jr_max,:) = ecoef * unu (jr_min:jr_max,:)**3 * dunu (jr_min:jr_max,:)

!-----------------------------------------------------------------------
!  Absorption and emission inverse mean free paths
!-----------------------------------------------------------------------

CALL abemrate( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )

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

DO j = jr_min,j_equil
  cf_e (j,:)             = ecoefa(j,:)/rho(j)
  cf_ye(j,:)             = ncoefa(j,:) * rmu/rho(j)
END DO ! j = jr_min,j_equil

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

  CALL abemrate( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye_new, nx )

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

  DO j = jr_min,j_equil
    DO k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  RHS0 e-neutrinos
!  (Right-hand side of the zero-moment transport equation for e-neutrinos)
!-----------------------------------------------------------------------

      RHS0(j,k)          = emis(j,k,1,ij_ray,ik_ray) + ( - emis(j,k,1,ij_ray,ik_ray) &
&                        - absor(j,k,1,ij_ray,ik_ray) ) * psi0(j,k)

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to Ye
!-----------------------------------------------------------------------

      RHS0_y(j,k)        = emisy(j,k,1) + ( - emisy(j,k,1) - absory(j,k,1) ) * psi0(j,k)

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nue
!-----------------------------------------------------------------------

      RHS0_p0(j,k,k)     = RHS0_p0(j,k,k) - emis(j,k,1,ij_ray,ik_ray) &
&                        - absor(j,k,1,ij_ray,ik_ray)

    END DO ! k = 1,nnugp(1)
  END DO ! j = jr_min,j_equil

!-----------------------------------------------------------------------
!
!             \\\\\ EQUILIBRATION LHS FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
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
  END DO ! j = jr_min,j_equil

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION LHS /////
!
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
    LHS0(j,ne_1+1)       = cdtinv * ( ye_new(j) - ye(j) )
    LHS0_y(j,ne_1+1)     = cdtinv
  END DO ! j = jr_min,j_equil

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION RHS /////
!
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
    DO k = 1,nnugp(1)
      RHS0(j,ne_1+1)     = RHS0(j,ne_1+1) &
&                        - cf_ye(j,k) * ( emis(j,k,1,ij_ray,ik_ray)       &
&                        - ( emis(j,k,1,ij_ray,ik_ray) + absor(j,k,1,ij_ray,ik_ray) ) * psi0(j,k) )
    END DO ! k = 1,nnugp(1)
  END DO ! j = jr_min,j_equil

  DO j = jr_min,j_equil
    DO k = 1,nnugp(1)
      RHS0_p0(j,ne_1+1,k) = cf_ye(j,k) * ( emis(j,k,1,ij_ray,ik_ray) + absor(j,k,1,ij_ray,ik_ray) )
      RHS0_y(j,ne_1+1)   = RHS0_y(j,ne_1+1) &
&                        - cf_ye(j,k) * ( emisy(j,k,1) - ( emisy(j,k,1) + absory(j,k,1) ) * psi0(j,k) )
    END DO ! k = 1,nnugp(1)
  END DO ! j = jr_min,j_equil

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

  A                  = - RHS0_p0

!-----------------------------------------------------------------------
!  A: k-equations diagonal elements involving LHS dpsi0(j-1/2)
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
    DO k = 1,ne_1
      A(j,k,k)       = A(j,k,k) + LHS0_p0(j,k)
    END DO ! k = 1,ne_1
  END DO ! j = jr_min,j_equil

!-----------------------------------------------------------------------
!  A: k-equations and border involving LHS and RHS dYw(j)
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
    DO k = 1,ne_1+1
      A(j,k,ne_1+1)  = A(j,k,ne_1+1) + LHS0_y(j,k) - RHS0_y(j,k)
    END DO ! k = 1,ne_1+1
  END DO ! j = jr_min,j_equil

!-----------------------------------------------------------------------
!  C: k-equations involving LHS and RHS constant terms
!-----------------------------------------------------------------------

  DO k = 1,ne_1+1
    DO j = jr_min,j_equil
      C(k,j)         = LHS0(j,k) - RHS0(j,k)
    END DO ! j = jr_min,j_equil
  END DO ! k = 1,ne_211

!-----------------------------------------------------------------------
!
!            \\\\\ SOLVE FOR THE UP RECURRSION MATRICES /////
!
!                     N
!           u     =  Sum  R       u       + R
!            k,j     k'=1  k,k',j  k',j+1    k,N+1,j
!
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil

!-----------------------------------------------------------------------
!  Load coefficient matrix A_coef with coefficients of the unknowns
!   u(j,k)
!-----------------------------------------------------------------------

    A_coef(:,:)      = A(j,:,:)

!-----------------------------------------------------------------------
!  Initialize the right-hand side matrix C_coef
!-----------------------------------------------------------------------

    C_coef           = zero

!-----------------------------------------------------------------------
!  Load right-hand side matrix C_coef with the constants terms in the
!   equations for the unknowns u(j,k)
!-----------------------------------------------------------------------

    C_coef(:,ne_1+1) = RHS0(j,:) - LHS0(j,:)

!-----------------------------------------------------------------------
!  Solve for the recursion coefficients relating the unknowns u(k,j)
!   to the unknowns u(k,j+1)
!-----------------------------------------------------------------------

    CALL dgesv( n_linear, n_rhs, A_coef, l_leading_A, i_piv, C_coef, &
&    l_leading_B, info )

!-----------------------------------------------------------------------
!  Load recursion matrix
!-----------------------------------------------------------------------

    R(:,:,j)         = C_coef(:,:)

  END DO ! j = jr_min,j_equil

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

  u                  = zero

!-----------------------------------------------------------------------
!  Recursively compute the solution vector
!-----------------------------------------------------------------------

  DO j = j_equil,jr_min,-1
    DO k = 1,ne_1+1
      DO kp = 1,ne_1
        u(k,j)       = u(k,j) + R(k,kp,j) * u(kp,j+1)
      END DO ! kp = 1,ne_1
      u(k,j)         = u(k,j) + R(k,ne_1+1,j)
    END DO ! k = 1,ne_1+1
  END DO ! j = j_equil,jr_min,-1

!-----------------------------------------------------------------------
!
!      \\\\\ TRANSFER INVCREMENTS AND TEST FOR CONVERGENCE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer increments
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil

    DO k = 1,nnugp(1)
      dpsi0_iph(j,k) = u(k,j)
    END DO ! k = 1,nnugp(1)

    dye_iph(j)       = u(ne_1+1,j)

  END DO ! j = jr_min,j_equil
 
!-----------------------------------------------------------------------
!        Increment variables
!
!  psi0(j,k,n) is allowed to exceed unity during iterations in order for
!   the iterations to converge.
!  psi0(j,k) is not allowed to fall below zero to avoid unphysical solutions
!
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
    DO k = 1,nnugp(1)
      psi0(j,k)      = DMAX1( psi0(j,k) + dpsi0_iph(j,k), psi0amin )
    END DO ! k = 1,nnugp(1)
  END DO ! j = jr_min,j_equil

  DO j = jr_min,j_equil
    ye_new(j)        = DMAX1( ye_new(j) + dye_iph(j), zero )
  END DO ! j = jr_min,j_equil

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

  endit              = .true.
  jiterp0            = 0
  kiterp0            = 0
  niterp0            = 0
  jiterye            = 0

loop_checkpsi: DO j = jr_min,j_equil
    DO k = 1,nnugp(1)
      IF ( DABS( dpsi0_iph(j,k)/DMAX1( psi0(j,k), tolpsimin ) ) >= tolnupsi ) THEN
        jiterp0      = j
        kiterp0      = k
        niterp0      = 1
        endit        = .false.
      EXIT loop_checkpsi
      END IF ! dpsi0_iph(j,k)/DMAX1( psi0(j,k) > tolnupsi
    END DO ! k = 1,nnugp(1)
  END DO loop_checkpsi

loop_checkye: DO j = jr_min,j_equil
    IF ( DABS( dye_iph(j)/DMAX1( ye_new(j), 0.1d0 ) ) >= tolnuye ) THEN
      jiterye        = j
      endit          = .false.
      EXIT loop_checkye
    END IF ! DABS( dye_iph(j)/DMAX1( ye(j), 0.1 ) ) >= tolnuye
  END DO loop_checkye

!-----------------------------------------------------------------------
!  WRITE diagnostics if it > 5
!-----------------------------------------------------------------------

  IF ( .not. endit  .and.  it >= 5 ) THEN
    IF ( jiterp0 == 0 ) THEN
      WRITE (nlog,1501) it, myid, ij_ray, ik_ray, jiterye, dye_iph(jiterye),                 &
&      dye_iph(jiterye)/DMAX1( ye_new(jiterye), 0.1d0 ), jiterye, ye_new(jiterye)
    ELSE IF ( jiterye == 0 ) THEN
      WRITE (nlog,1503) it, myid, ij_ray, ik_ray, jiterp0, kiterp0, niterp0,                 &
&      dpsi0_iph(jiterp0,kiterp0),                                                           &
&      dpsi0_iph(jiterp0,kiterp0)/DMAX1( psi0(jiterp0,kiterp0), tolpsimin ),                 &
&      jiterp0, kiterp0, niterp0, psi0(jiterp0,kiterp0)
    ELSE
      WRITE (nlog,1505) it, myid, ij_ray, ik_ray, jiterye, dye_iph(jiterye),dye_iph(jiterye) &
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
    WRITE (nlog,7001) j_ray, k_ray, myid
    WRITE (nprint,7001) j_ray, k_ray, myid
    STOP
    EXIT
  END IF

END DO ! iteration loop

!-----------------------------------------------------------------------
!
!     \\\\\ BOUNDARY VALUES AND LARGEST RELATIVE INCREMENTS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Boundary value of psi0 for jr_min-1
!-----------------------------------------------------------------------

DO k = 1,nnugp(1)
  psi0(jr_min-1,k)     = psi0(jr_min,k)
END DO ! k = 1,nnugp(n)

!-----------------------------------------------------------------------
!
!                 \\\\\ UPDATE THE TEMPERATURE /////
!
!-----------------------------------------------------------------------

IF ( consrve == 'e' ) THEN

!-----------------------------------------------------------------------
!  Compute energy transferred to neutrinos
!-----------------------------------------------------------------------

  de                   = zero

  DO k = 1,nnugp(1)
    DO j = jr_min,j_equil
      de(j)            = de(j) + cf_e(j,k) * stwt(1) * cdt                &
&                      * ( RHS0(j,k) + RHS0_y(j,k) * dye_iph(j) )
      DO kp = 1,nnugp(1)
        de(j)          = de(j) + cf_e(j,k) * stwt(1) * cdt * RHS0_p0(j,k,     kp) * dpsi0_iph(j,kp)
      END DO ! kp = 1,nnugp(1)
    END DO ! j = jr_min,j_equil
  END DO ! k = 1,nnugp(1)

!-----------------------------------------------------------------------
!  Compute the new temperature given by energy conservation
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
    CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye_0(j), e, dedd, dedt, dedy )
    CALL e_degen(  rho(j), t(j), ye_0(j), e_eta )
    IF ( e_eta > degen ) THEN
      t_new(j)         = t_new(j) + ( - de(j) - dedy * ( ye_new(j) - ye_0(j) ) )/dedt
    ELSE
      CALL tgvndeye_x( j, ij_ray, ik_ray, rho(j), e - de(j), ye_new(j), t_new(j), t_new(j) )
    END IF ! e_eta > degen
  END DO ! j = jr_min,j_equil

ELSE ! consrve = 's'

!-----------------------------------------------------------------------
!  Preserve entropy
!-----------------------------------------------------------------------

  DO j = jr_min,j_equil
    CALL eqstt_x( 3, j, ij_ray, ik_ray, rho(j), t(j), ye(j), s_i, dsdd, dsdt, dsdy )
    CALL tgvndsye_x( j, ij_ray, ik_ray, rho(j), s_i, ye_new(j), t(j), t_new(j) )
  END DO ! j = jr_min,j_equil

END IF ! consrve = 'e'

!-----------------------------------------------------------------------
!  Conserve leptons
!-----------------------------------------------------------------------

DO j = jr_min,j_equil
  d_ye(j)              = - SUM( cf_ye(j,:) * ( psi0(j,:) - psi0_i(j,:) ) )
  ye_new(j)            = ye_0(j) + d_ye(j)
END DO ! j = jr_min,j_equil

!-----------------------------------------------------------------------
!  Restore variable for transfer back to RADIAL_RAY_MODULE
!-----------------------------------------------------------------------

t_c   (imin:imax,ij_ray,ik_ray)     = t_new (jr_min:jr_max)
ye_c  (imin:imax,ij_ray,ik_ray)     = ye_new(jr_min:jr_max)
psi0_c(imin:imax,:,1,ij_ray,ik_ray) = psi0  (jr_min:jr_max,:)

!-----------------------------------------------------------------------
!
!                  \\\\\ DEALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

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

RETURN
END SUBROUTINE equilibrate_x
