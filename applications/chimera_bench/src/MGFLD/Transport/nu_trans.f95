SUBROUTINE nu_trans( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& rho, t, t_new, ye, ye_new, radius, rstmss, nx, nez, nnu, jdt, dtnph_trans )
!-----------------------------------------------------------------------
!
!    File:         nu_trans
!    Module:       nu_trans
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/5/05
!
!    Purpose:
!      To compute iterate the change in Ye and psi0 due to absorption,
!       emission, inelastic scattering, pair production, and transport.
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
!  jr_min           : inner radial zone number
!  jr_max           : outer radial zone number
!  ij_ray           : j-index of a radial ray
!  ik_ray           : k-index of a radial ray
!  ij_ray_dim       : number of y-zones on a processor before swapping
!  ik_ray_dim       : number of z-zones on a processor before swapping
!  rho              : density (g cm^{-3})
!  t                : temperature (K)
!  ye               : electron fraction
!  radius           : radius of zone j (cm)
!  rstmss           : enclosed rest mass (g)
!  nx               : x-array extent
!  nez              : neutrino energy array extent
!  nnu              : neutrino flavor array extent
!  dtnph_trans      : source and transport time step
!
!    Output arguments:
!  t_new            : updated temperature (K)
!  ye_new           : updated electron fraction
!  jdt              : radial zone causing minimum time step for criteria i, radial ray ij_ray, ik_ray
!
!    Input arguments (common):
!
!  nnugp(n)         : number of energy zones for neutrinos of type n
!  unu(j,k)         : midpoint energy of energy zone k at radial zone j
!  dunu(j,k)        : energy width of energy zone k at radial zone j
!  tolnupsi         : abs(dpsi_iph)/( psi0 + tolpsimin ) < tolnupsi for
!  tolpsimin        :  for convergence
!  psi0(j,k,n)      : zeroth angular moment of the neutrino distribution
!                function
!  dtau_nutrns(j)
!                   : proper time step at j-1/2
!  cdt(j)           : ! * dtau_nutrns(j)
!  cdtinv(j)        : 1/cdt(j)
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
!  kind_module, array_module, numerical_module, physcnst_module
!  abem_module, brem_module, edit_module incrmnt_module,
!  it_tol_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module, pair_module, prb_cntl_module, scat_a_module,
!  scat_e_module, scat_n_module, scat_nn_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc_y
USE numerical_module, ONLY : zero, one, ncoef, epsilon, ergmev
USE physcnst_module, ONLY : rmu, cvel

USE abem_module, ONLY : emis, emisy, absor, absory
USE brem_module, ONLY : baf, bafy, bafp0
USE edit_module, ONLY : nprint, nlog
USE evh1_global, ONLY : degen
USE incrmnt_module, ONLY : dtmpnn, dye
USE it_tol_module, ONLY : iternu, tolpsimin, tolnupsi, tolnuye
USE mdl_cnfg_module, ONLY : dmrst, u_vel=>u
USE nu_dist_module, ONLY : psi0, psi1, ncoefa, ecoefa, vol, drjmh_inv, &
& agrjmh_nu, area, agr_nu, gamgr_nu, dc, dcr, stwt, unujrad, nnujrad, &
& dunue, dnurad, nnurad, unukrad, nnukrad, unulum, unurad, fluxnuk, &
& fluxnu, unue, unujea, unucrea, dunujeadt, dudt_nu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE pair_module, ONLY : paf, pafy, pafp0
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY : idiff
USE scat_e_module, ONLY : scef, scefy, scefp0
USE scat_n_module, ONLY : scnf, scnfy, scnfp0
USE scat_nn_module, ONLY : scnnf, scnnfy, scnnfp0
USE scat_nA_module, ONLY : scnAf, scnAfy, scnAfp0
USE t_cntrl_module, ONLY : dpsisp, psimin, cdt, cdtinv, dtauj_nutrns, &
& dtau_nutrns

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nx            ! radial array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rho         ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: t           ! initial temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: ye          ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: radius      ! radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rstmss      ! radius (cm)
REAL(KIND=double), INTENT(in)                   :: dtnph_trans ! source and transport time step

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim) :: jdt    ! zone causing dt

REAL(KIND=double), INTENT(out), DIMENSION(nx)   :: t_new  ! updated temperature (K)
REAL(KIND=double), INTENT(out), DIMENSION(nx)   :: ye_new ! updated electron fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

LOGICAL                          :: endit

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: kp            ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: it            ! iteration index
INTEGER                          :: jiterp0       ! radial zone for which psi0 iteration fails
INTEGER                          :: kiterp0       ! neutrino energy zone for which psi0 iteration fails
INTEGER                          :: niterp0       ! neutrino flavor for which psi0 iteration fails
INTEGER                          :: jiterye       ! radial zone number for which ye iteration fails

INTEGER                          :: j_ray         ! polar index of the radial ray
INTEGER                          :: k_ray         ! azimuthal index of the radial ray

INTEGER                          :: ne_1          ! nnugp(1)
INTEGER                          :: ne_2          ! nnugp(2)
INTEGER                          :: ne_3          ! nnugp(3)
INTEGER                          :: ne_4          ! nnugp(4)

INTEGER                          :: istat         ! allocation status
INTEGER                          :: lwork         ! work area flag

INTEGER                          :: n_linear      ! number of linear equations to be solved by dgesv
INTEGER                          :: n_rhs         ! number of right-hand sides of the linear equations
INTEGER                          :: l_leading_A   ! leading dimension of the coefficient matrix A_coef
INTEGER                          :: l_leading_B   ! leading dimension of the coefficient matrix C_coef
INTEGER, DIMENSION(4*nez+1)      :: i_piv         ! the oivit indices that define the permutation matrix P
INTEGER                          :: info          ! 0: successful exit, < 0: |info| argument illegal, > 0: u(info,info) = 0

REAL(KIND=double)                :: a0w           ! coefficient for computing change in psi0
REAL(KIND=double)                :: b0w           ! coefficient for computing change in psi0
REAL(KIND=double)                :: c0w           ! coefficient for computing change in psi0
REAL(KIND=double)                :: a1w           ! coefficient for computing change in psi1
REAL(KIND=double)                :: b1w           ! coefficient for computing change in psi1
REAL(KIND=double)                :: c1w           ! coefficient for computing change in psi1

REAL(KIND=double)                :: da0wdy        ! d/dy coefficient for computing change in psi0
REAL(KIND=double)                :: db0wdy        ! d/dy coefficient for computing change in psi0
REAL(KIND=double)                :: dc0wdy        ! d/dy coefficient for computing change in psi0
REAL(KIND=double)                :: da1wdy        ! d/dy coefficient for computing change in psi1
REAL(KIND=double)                :: db1wdy        ! d/dy coefficient for computing change in psi1
REAL(KIND=double)                :: dc1wdy        ! d/dy coefficient for computing change in psi1

REAL(KIND=double)                :: da0wdp0       ! d/dpsi0 coefficient for computing change in psi0
REAL(KIND=double)                :: db0wdp0       ! d/dpsi0 coefficient for computing change in psi0
REAL(KIND=double)                :: dc0wdp0       ! d/dpsi0 coefficient for computing change in psi0
REAL(KIND=double)                :: da1wdp0       ! d/dpsi0 coefficient for computing change in psi1
REAL(KIND=double)                :: db1wdp0       ! d/dpsi0 coefficient for computing change in psi1
REAL(KIND=double)                :: dc1wdp0       ! d/dpsi0 coefficient for computing change in psi1

REAL(KIND=double), DIMENSION(nnu):: psi_ratio     ! ratio of psi0 in outer ghost to outer psi0

REAL(KIND=double)                :: denu          ! maximum relative change of psi0 unutc

REAL(KIND=double)                :: e             ! inernal energy per unit mass
REAL(KIND=double)                :: dedd          ! d(e)/d(rho)
REAL(KIND=double)                :: dedt          ! d(e)/d(t)
REAL(KIND=double)                :: dedy          ! d(e)/d(ye)

REAL(KIND=double)                :: dunudt        ! rate of energy transfer to neutrinos by emission and absorption

REAL(KIND=double)                :: e_eta         ! (electron Fermi energy)/KT

REAL(KIND=double), PARAMETER     :: psi0amax = 1.d0 ! maximum value of psi0
REAL(KIND=double), PARAMETER     :: psi0amin = 0.d0 ! minimum value of psi0

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0         ! coefficients of zero-order terms on left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_p0p     ! coefficients of dpsi0(j+1) terms on left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_p0m     ! coefficients of dpsi0(j-1) terms left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_p0      ! coefficients of dpsi0(j) terms left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_yp      ! coefficients of dYe(j+1) terms on left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_ym      ! coefficients of dYe(j-1) terms left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: LHS0_y       ! coefficients of dYe(j) terms left-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: RHS0         ! coefficients of rigit-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: RHS0_y       ! coefficients of d/dYe rigit-hand side of p0 equation
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: RHS0_p0      ! coefficients of d/dpsi0 rigit-hand side of p0 equation

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye_0         ! initial value Ye
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dye_iph      ! increment of Ye for the ith iteration
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi0_i       ! initial value of psi0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dpsi0_iph    ! increment of psi0 for the ith iteration
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: cf_e         ! energy transfer coefficients
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: cf_ye        ! lepton transfer coefficients
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: de_i         ! energy change due to n-neutrinos and ye
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: de           ! total energy change due to neutrinos

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: A            ! coefficient matrix for the j terms in j-equation, all j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: D_m          ! coefficient matrix for the j-1 terms in j-equation, all j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: D_p          ! coefficient matrix for the j+1 terms in j-equation, all j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: C            ! constant terms in j-equation, all j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: u            ! solution vector of the linear equations

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: R            ! recursion matrix
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: A_coef       ! Coefficient matrix of the variables to be solved
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: C_coef       ! Coefficient matrix of quantities on the right-hand side


REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: a2vol_inv    ! agrh**2/vol
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: agr2_inv     ! 1/agr**2

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: ddcpsjph     !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: ddcpsjmh     !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: ddcdyjph     !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: ddcdyjmh     !

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in nu_trans')
 1501 FORMAT (' Iteration attempt',i4,' myid=',i4,' ij_ray=',i4,' ik_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=',es11.3,' ye(',i4,')=',es11.3)
 1503 FORMAT (' Iteration attempt',i4,' myid=',i4,' ij_ray=',i4,' ik_ray=',i4, &
&' dpsi0_iph(',i4,',',i4,',',i4,')=', es11.3,' D_psi=',es11.3,                 &
&' psi0(',i4,',',i4,',',i4,')=', es11.3)
 1505 FORMAT (' Iteration attempt',i4,' myid=',i4,' ij_ray=',i4,' ik_ray=',i4, &
& ' dye_iph(',i4,')=',es11.3,' D_ye=', es11.3,' ye(',i4,')=',es11.3,           &
& ' dpsi0_iph(',i4,',',i4,',',i4,')=',es11.3,' D_psi=',es11.3,                 &
& ' psi0(',i4,',',i4,',',i4,')=', es11.3)
 2001 FORMAT (' Deallocation problem for array ',a10,' in nu_trans')
 7001 FORMAT (' Convergence failure in subroutine nu_trans, j_ray=',i4,        &
& ' k_ray=',i4,' myid=',i4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Initialize independent variables and increments................

dtmpnn                   = zero
dye                      = zero

!........Return if nnugpmx = 0..........................................

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!                   \\\\\ ALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (LHS0(nx,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_p0p(nx,4*nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0p  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_p0m(nx,4*nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0m  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_p0(nx,4*nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_yp(nx,4*nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_yp   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_ym(nx,4*nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_ym   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (LHS0_y(nx,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_y    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (RHS0(nx,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_y(nx,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_y    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (RHS0_p0(nx,4*nez+1,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'RHS0_p0   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ye_0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dye_iph(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye_iph   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0_i(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi0_iph(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi0_iph '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_e(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_e      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cf_ye(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cf_ye     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de_i(nx,nnu+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de_i      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (A(nx,4*nez+1,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (D_m(4*nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'D_m       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (D_p(4*nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'D_p       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C(4*nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(4*nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (R(4*nez+1,4*nez+1,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'R         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (A_coef(4*nez+1,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A_coef    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (C_coef(4*nez+1,4*nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'C_coef    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (a2vol_inv(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a2vol_inv '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (agr2_inv(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr2_inv  '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ddcpsjph(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcpsjph  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ddcpsjmh(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcpsjmh  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ddcdyjph(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcdyjph  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ddcdyjmh(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcdyjmh  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_ray              = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
k_ray              = ( myid/n_proc_y ) * ik_ray_dim + ik_ray

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

LHS0                     = zero
LHS0_p0p                 = zero
LHS0_p0m                 = zero
LHS0_p0                  = zero
LHS0_yp                  = zero
LHS0_ym                  = zero
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
D_m                      = zero
D_p                      = zero
C                        = zero
u                        = zero

R                        = zero
A_coef                   = zero
C_coef                   = zero

a2vol_inv                = zero
agr2_inv                 = zero

ddcpsjph                 = zero
ddcpsjmh                 = zero
ddcdyjph                 = zero
ddcdyjmh                 = zero

lwork                    = 0

ne_1                     = nnugp(1)
ne_2                     = nnugp(1) + nnugp(2)
ne_3                     = nnugp(1) + nnugp(2) + nnugp(3)
ne_4                     = nnugp(1) + nnugp(2) + nnugp(3) + nnugp(4)

n_linear                 = 4*nez+1
n_rhs                    = 4*nez+1
l_leading_A              = 4*nez+1
l_leading_B              = 4*nez+1

!-----------------------------------------------------------------------
!  Store initial values of psi0 and ye
!-----------------------------------------------------------------------

psi0_i                   = psi0
ye_0                     = ye

!-----------------------------------------------------------------------
!  Get time step
!-----------------------------------------------------------------------

CALL dtau_trns( jr_min, jr_max, dtnph_trans )

!-----------------------------------------------------------------------
!  Compute Ye transfer constants
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  DO j = jr_min,jr_max
    cf_e(j,k)            = ecoefa(j,k)/rho(j)
    cf_ye(j,k)           = ncoefa(j,k) * rmu/rho(j)
  END DO ! j
END DO ! k

!-----------------------------------------------------------------------
!  Compute a2vol_inv
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  a2vol_inv(j)           = agrjmh_nu(j) * agrjmh_nu(j)/vol(j)
  agr2_inv(j)            = one/( agr_nu(j) * agr_nu(j) )
END DO ! j

!-----------------------------------------------------------------------
!  Initialize maximum relative change of psi0
!-----------------------------------------------------------------------

dpsisp                   = zero

ddcpsjph                 = zero
ddcpsjmh                 = zero
ddcdyjph                 = zero
ddcdyjmh                 = zero

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
  LHS0_p0p               = zero          ! d(transport functions)/dpsi0(j+1,k)
  LHS0_p0m               = zero          ! d(transport functions)/dpsi0(j-1,k)
  LHS0_p0                = zero          ! d(transport functions)/dpsi0(j,k)
  LHS0_yp                = zero          ! d(transport functions)/dYe(j+1)
  LHS0_ym                = zero          ! d(transport functions)/dYe(j-1)
  LHS0_y                 = zero          ! d(transport functions)/dYe(j)

  RHS0                   = zero          ! source functions
  RHS0_y                 = zero          ! d( source functions )/dYe
  RHS0_p0                = zero          ! d( source functions )/dpsi0

  dye_iph                = zero          ! increments of Ye
  dpsi0_iph              = zero          ! increments of psi0
  
  D_m                    = zero          ! equations involving LHS dpsi0(j-3/2)
  D_p                    = zero          ! equations involving LHS dpsi0(j+1/2)

  R                      = zero          ! recursion matrix

!-----------------------------------------------------------------------
!
!         \\\\\ APPLY INITIAL SURFACE BOUNDARY CONDITION /////
!
!               i                             i
!           psi0           = psi_ratio  * psi0
!               k,jr_max+1            k        k,jr_max
!
!-----------------------------------------------------------------------

  IF ( nnugp(1) /= 0 ) THEN
    DO k = 1,nnugp(1)
      CALL psi_bd( jr_max, k, 1, radius, nx, psi_ratio(1) )
      psi0(jr_max+1,k,1) = psi_ratio(1) * psi0(jr_max,k,1)
    END DO ! k
  END IF ! nnugp(1) /= 0
    
  IF ( nnugp(2) /= 0 ) THEN
    DO k = 1,nnugp(2)
      CALL psi_bd( jr_max, k, 2, radius, nx, psi_ratio(2) )
      psi0(jr_max+1,k,2) = psi_ratio(2) * psi0(jr_max,k,2)
    END DO ! k
  END IF ! nnugp(2) /= 0
    
  IF ( nnugp(3) /= 0 ) THEN
    DO k = 1,nnugp(3)
      CALL psi_bd( jr_max, k, 3, radius, nx, psi_ratio(3) )
      psi0(jr_max+1,k,3) = psi_ratio(3) * psi0(jr_max,k,3)
    END DO ! k
  END IF ! nnugp(3) /= 0
    
  IF ( nnugp(4) /= 0 ) THEN
    DO k = 1,nnugp(4)
      CALL psi_bd( jr_max, k, 4, radius, nx, psi_ratio(4) )
      psi0(jr_max+1,k,4) = psi_ratio(4) * psi0(jr_max,k,4)
    END DO ! k
  END IF ! nnugp(4) /= 0

!-----------------------------------------------------------------------
!  Compute psi1 and refresh all the absoprtion and scattering rates
!-----------------------------------------------------------------------

  CALL psi1_cal( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&   rho, t_new, ye_new, radius, rstmss, u_vel, psi0, psi1, nx, nez, nnu, it )

!-----------------------------------------------------------------------
!  Compute dc_pjph and ddcpsjmh
!-----------------------------------------------------------------------

  CALL ddc_dpsi( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&  radius, ddcpsjph, ddcpsjmh, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!           \\\\\ GET SOURCE FUNCTIONS FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(1) /= 0 ) THEN

!-----------------------------------------------------------------------
!  Compute and store combinations of e-neutrino interaction moments
!   appearing in the p2 equation for psi0.
!
!        tot mfp^-1 : emis + abs - b1  (total inverse mean free path)
!-----------------------------------------------------------------------

    DO k = 1,nnugp(1)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  RHS0 e-neutrinos
!  (Right-hand side of the zero-moment transport equation for e-neutrinos)
!-----------------------------------------------------------------------

        a0w              = scef (1,j,k,1) + paf(1,j,k,1) + scnnf(1,j,k,1) &
&                        + scnAf(1,j,k,1) + baf(1,j,k,1) + scnf(1,j,k,1)
        b0w              = scef (2,j,k,1) + paf(2,j,k,1) + scnnf(2,j,k,1) &
&                        + scnAf(2,j,k,1) + baf(2,j,k,1) + scnf(2,j,k,1)
        c0w              = scef (3,j,k,1) + paf(3,j,k,1) + scnnf(3,j,k,1) &
&                        + scnAf(3,j,k,1) + baf(3,j,k,1) + scnf(3,j,k,1)
        a1w              = scef (4,j,k,1) + paf(4,j,k,1) + scnnf(4,j,k,1) &
&                        + scnAf(4,j,k,1) + baf(4,j,k,1) + scnf(4,j,k,1)
        b1w              = scef (5,j,k,1) + paf(5,j,k,1) + scnnf(5,j,k,1) &
&                        + scnAf(5,j,k,1) + baf(5,j,k,1) + scnf(5,j,k,1)
        c1w              = scef (6,j,k,1) + paf(6,j,k,1) + scnnf(6,j,k,1) &
&                        + scnAf(6,j,k,1) + baf(6,j,k,1) + scnf(6,j,k,1)

        RHS0(j,k)        = emis(j,k,1,ij_ray,ik_ray) + ( a0w - emis(j,k,1,ij_ray,ik_ray) &
&                        - absor(j,k,1,ij_ray,ik_ray) ) * psi0(j,k,1) &
&                        + b0w * psi1(j,k,1) + c0w

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to Ye
!-----------------------------------------------------------------------

        da0wdy           = scefy (1,j,k,1) + pafy(1,j,k,1) + scnnfy(1,j,k,1) &
&                        + scnAfy(1,j,k,1) + bafy(1,j,k,1) + scnfy(1,j,k,1)
        db0wdy           = scefy (2,j,k,1) + pafy(2,j,k,1) + scnnfy(2,j,k,1) &
&                        + scnAfy(2,j,k,1) + bafy(2,j,k,1) + scnfy(2,j,k,1)
        dc0wdy           = scefy (3,j,k,1) + pafy(3,j,k,1) + scnnfy(3,j,k,1) &
&                        + scnAfy(3,j,k,1) + bafy(3,j,k,1) + scnfy(3,j,k,1)
        da1wdy           = scefy (4,j,k,1) + pafy(4,j,k,1) + scnnfy(4,j,k,1) &
&                        + scnAfy(4,j,k,1) + bafy(4,j,k,1) + scnfy(4,j,k,1)
        db1wdy           = scefy (5,j,k,1) + pafy(5,j,k,1) + scnnfy(5,j,k,1) &
&                        + scnAfy(5,j,k,1) + bafy(5,j,k,1) + scnfy(5,j,k,1)
        dc1wdy           = scefy (6,j,k,1) + pafy(6,j,k,1) + scnnfy(6,j,k,1) &
&                        + scnAfy(6,j,k,1) + bafy(6,j,k,1) + scnfy(6,j,k,1)

        RHS0_y(j,k)      = emisy(j,k,1) + ( da0wdy - emisy(j,k,1) - absory(j,k,1) ) * psi0(j,k,1) &
&                        + db0wdy * psi1(j,k,1) + dc0wdy

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nue
!-----------------------------------------------------------------------

        DO kp = 1,nnugp(1)

          da0wdp0        = scefp0 (1,kp,j,k,1) + scnnfp0(1,kp,j,k,1) &
&                        + scnAfp0(1,kp,j,k,1) + scnfp0(1,kp,j,k,1)
          db0wdp0        = scefp0 (2,kp,j,k,1) + scnnfp0(2,kp,j,k,1) &
&                        + scnAfp0(2,kp,j,k,1) + scnfp0(2,kp,j,k,1)
          dc0wdp0        = scefp0 (3,kp,j,k,1) + scnnfp0(3,kp,j,k,1) &
&                        + scnAfp0(3,kp,j,k,1) + scnfp0(3,kp,j,k,1)
          da1wdp0        = scefp0 (4,kp,j,k,1) + scnnfp0(4,kp,j,k,1) &
&                        + scnAfp0(4,kp,j,k,1) + scnfp0(4,kp,j,k,1)
          db1wdp0        = scefp0 (5,kp,j,k,1) + scnnfp0(5,kp,j,k,1) &
&                        + scnAfp0(5,kp,j,k,1) + scnfp0(5,kp,j,k,1)
          dc1wdp0        = scefp0 (6,kp,j,k,1) + scnnfp0(6,kp,j,k,1) &
&                        + scnAfp0(6,kp,j,k,1) + scnfp0(6,kp,j,k,1)

          RHS0_p0(j,k,kp) = da0wdp0 * psi0(j,k,1) + db0wdp0 * psi1(j,k,1) + dc0wdp0

        END DO ! kp

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nuebar
!-----------------------------------------------------------------------

        IF ( nnugp(2) /= 0 ) THEN
          DO kp = 1,nnugp(2)

            da0wdp0      = pafp0(1,kp,j,k,1) + bafp0(1,kp,j,k,1)
            db0wdp0      = pafp0(2,kp,j,k,1) + bafp0(2,kp,j,k,1)
            dc0wdp0      = pafp0(3,kp,j,k,1) + bafp0(3,kp,j,k,1)
            da1wdp0      = pafp0(4,kp,j,k,1) + bafp0(4,kp,j,k,1)
            db1wdp0      = pafp0(5,kp,j,k,1) + bafp0(5,kp,j,k,1)
            dc1wdp0      = pafp0(6,kp,j,k,1) + bafp0(6,kp,j,k,1)

            RHS0_p0(j,k,ne_1+kp) = da0wdp0 * psi0(j,k,1) + db0wdp0 * psi1(j,k,1) + dc0wdp0

          END DO ! kp
        END IF ! nnugp(2) /= 0

!-----------------------------------------------------------------------
!  Diagonal contribution to RHS0 (k' = k)
!-----------------------------------------------------------------------

        RHS0_p0(j,k,k) = RHS0_p0(j,k,k) + a0w - emis(j,k,1,ij_ray,ik_ray) &
&                      - absor(j,k,1,ij_ray,ik_ray)

      END DO ! j
    END DO ! k
  END IF ! nnugp(1) /= 0

!-----------------------------------------------------------------------
!
!         \\\\\ GET SOURCE FUNCTIONS FOR E-ANTINEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(2) /= 0 ) THEN

!-----------------------------------------------------------------------
!  Compute and store combinations of e-antineutrino interaction moments
!   appearing in the p2 equation for psi0.
!
!        tot mfp^-1 : emis + abs - b1  (total inverse mean free path)
!-----------------------------------------------------------------------

    DO k = 1,nnugp(2)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  RHS0 e-antineutrinos
!  (Right-hand side of the zero-moment transport equation for e-antineutrinos)
!-----------------------------------------------------------------------

        a0w              = scef (1,j,k,2) + paf(1,j,k,2) + scnnf(1,j,k,2) &
&                        + scnAf(1,j,k,2) + baf(1,j,k,2) + scnf(1,j,k,2)
        b0w              = scef (2,j,k,2) + paf(2,j,k,2) + scnnf(2,j,k,2) &
&                        + scnAf(2,j,k,2) + baf(2,j,k,2) + scnf(2,j,k,2)
        c0w              = scef (3,j,k,2) + paf(3,j,k,2) + scnnf(3,j,k,2) &
&                        + scnAf(3,j,k,2) + baf(3,j,k,2) + scnf(3,j,k,2)
        a1w              = scef (4,j,k,2) + paf(4,j,k,2) + scnnf(4,j,k,2) &
&                        + scnAf(4,j,k,2) + baf(4,j,k,2) + scnf(4,j,k,2)
        b1w              = scef (5,j,k,2) + paf(5,j,k,2) + scnnf(5,j,k,2) &
&                        + scnAf(5,j,k,2) + baf(5,j,k,2) + scnf(5,j,k,2)
        c1w              = scef (6,j,k,2) + paf(6,j,k,2) + scnnf(6,j,k,2) &
&                        + scnAf(6,j,k,2) + baf(6,j,k,2) + scnf(6,j,k,2)

        RHS0(j,ne_1+k) = emis(j,k,2,ij_ray,ik_ray) + ( a0w - emis(j,k,2,ij_ray,ik_ray) &
&                        - absor(j,k,2,ij_ray,ik_ray) ) * psi0(j,k,2)     &
&                        + b0w * psi1(j,k,2) + c0w

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to Ye
!-----------------------------------------------------------------------

        da0wdy           = scefy (1,j,k,2) + pafy(1,j,k,2) + scnnfy(1,j,k,2) &
&                        + scnAfy(1,j,k,2) + bafy(1,j,k,2) + scnfy(1,j,k,2)
        db0wdy           = scefy (2,j,k,2) + pafy(2,j,k,2) + scnnfy(2,j,k,2) &
&                        + scnAfy(2,j,k,2) + bafy(2,j,k,2) + scnfy(2,j,k,2)
        dc0wdy           = scefy (3,j,k,2) + pafy(3,j,k,2) + scnnfy(3,j,k,2) &
&                        + scnAfy(3,j,k,2) + bafy(3,j,k,2) + scnfy(3,j,k,2)
        da1wdy           = scefy (4,j,k,2) + pafy(4,j,k,2) + scnnfy(4,j,k,2) &
&                        + scnAfy(4,j,k,2) + bafy(4,j,k,2) + scnfy(4,j,k,2)
        db1wdy           = scefy (5,j,k,2) + pafy(5,j,k,2) + scnnfy(5,j,k,2) &
&                        + scnAfy(5,j,k,2) + bafy(5,j,k,2) + scnfy(5,j,k,2)
        dc1wdy           = scefy (6,j,k,2) + pafy(6,j,k,2) + scnnfy(6,j,k,2) &
&                        + scnAfy(6,j,k,2) + bafy(6,j,k,2) + scnfy(6,j,k,2)

        RHS0_y(j,ne_1+k) = emisy(j,k,2) + ( da0wdy - emisy(j,k,2) - absory(j,k,2) ) * psi0(j,k,2) &
&                        + db0wdy * psi1(j,k,2) + dc0wdy

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nue
!-----------------------------------------------------------------------

        DO kp = 1,nnugp(2)

          da0wdp0        = scefp0 (1,kp,j,k,2) + scnnfp0(1,kp,j,k,2) &
&                        + scnAfp0(1,kp,j,k,2) + scnfp0 (1,kp,j,k,2)
          db0wdp0        = scefp0 (2,kp,j,k,2) + scnnfp0(2,kp,j,k,2) &
&                        + scnAfp0(2,kp,j,k,2) + scnfp0 (2,kp,j,k,2)
          dc0wdp0        = scefp0 (3,kp,j,k,2) + scnnfp0(3,kp,j,k,2) &
&                        + scnAfp0(3,kp,j,k,2) + scnfp0 (3,kp,j,k,2)
          da1wdp0        = scefp0 (4,kp,j,k,2) + scnnfp0(4,kp,j,k,2) &
&                        + scnnfp0(4,kp,j,k,2) + scnfp0 (4,kp,j,k,2)
          db1wdp0        = scefp0 (5,kp,j,k,2) + scnnfp0(5,kp,j,k,2) &
&                        + scnAfp0(5,kp,j,k,2) + scnfp0 (5,kp,j,k,2)
          dc1wdp0        = scefp0 (6,kp,j,k,2) + scnnfp0(6,kp,j,k,2) &
&                        + scnAfp0(6,kp,j,k,2) + scnfp0 (6,kp,j,k,2)

          RHS0_p0(j,ne_1+k,ne_1+kp) = da0wdp0 * psi0(j,k,2) + db0wdp0 * psi1(j,k,2) + dc0wdp0

        END DO ! kp

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nuebar
!-----------------------------------------------------------------------

        IF ( nnugp(1) /= 0 ) THEN
          DO kp = 1,nnugp(1)

            da0wdp0      = pafp0(1,kp,j,k,2) + bafp0(1,kp,j,k,2)
            db0wdp0      = pafp0(2,kp,j,k,2) + bafp0(2,kp,j,k,2)
            dc0wdp0      = pafp0(3,kp,j,k,2) + bafp0(3,kp,j,k,2)
            da1wdp0      = pafp0(4,kp,j,k,2) + bafp0(4,kp,j,k,2)
            db1wdp0      = pafp0(5,kp,j,k,2) + bafp0(5,kp,j,k,2)
            dc1wdp0      = pafp0(6,kp,j,k,2) + bafp0(6,kp,j,k,2)

            RHS0_p0(j,ne_1+k,kp) = da0wdp0 * psi0(j,k,2) + db0wdp0 * psi1(j,k,2) + dc0wdp0

          END DO ! kp
        END IF ! nnugp(2) /= 0

!-----------------------------------------------------------------------
!  Diagonal contribution to RHS0 (k' = k)
!-----------------------------------------------------------------------

        RHS0_p0(j,ne_1+k,ne_1+k) = RHS0_p0(j,ne_1+k,ne_1+k) + a0w       &
&                                - emis(j,k,2,ij_ray,ik_ray) - absor(j,k,2,ij_ray,ik_ray)

      END DO ! j
    END DO ! k
  END IF ! nnugp(2) /= 0

!-----------------------------------------------------------------------
!
!           \\\\\ GET SOURCE FUNCTIONS FOR X-NEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(3) /= 0 ) THEN

!-----------------------------------------------------------------------
!  Compute and store combinations of x-neutrino interaction moments
!   appearing in the p2 equation for psi0.
!
!        tot mfp^-1 : emis + abs - b1  (total inverse mean free path)
!-----------------------------------------------------------------------

    DO k = 1,nnugp(3)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  RHS0 x-neutrinos
!  (Right-hand side of the zero-moment transport equation for x-neutrinos)
!-----------------------------------------------------------------------

        a0w              = scef (1,j,k,3) + paf(1,j,k,3) + scnnf(1,j,k,3) &
&                        + scnAf(1,j,k,3) + baf(1,j,k,3) + scnf(1,j,k,3)
        b0w              = scef (2,j,k,3) + paf(2,j,k,3) + scnnf(2,j,k,3) &
&                        + scnAf(2,j,k,3) + baf(2,j,k,3) + scnf(2,j,k,3)
        c0w              = scef (3,j,k,3) + paf(3,j,k,3) + scnnf(3,j,k,3) &
&                        + scnAf(3,j,k,3) + baf(3,j,k,3) + scnf(3,j,k,3)
        a1w              = scef (4,j,k,3) + paf(4,j,k,3) + scnnf(4,j,k,3) &
&                        + scnAf(4,j,k,3) + baf(4,j,k,3) + scnf(4,j,k,3)
        b1w              = scef (5,j,k,3) + paf(5,j,k,3) + scnnf(5,j,k,3) &
&                        + scnAf(5,j,k,3) + baf(5,j,k,3) + scnf(5,j,k,3)
        c1w              = scef (6,j,k,3) + paf(6,j,k,3) + scnnf(6,j,k,3) &
&                        + scnAf(6,j,k,3) + baf(6,j,k,3) + scnf(6,j,k,3)

        RHS0(j,ne_2+k) = a0w * psi0(j,k,3) + b0w * psi1(j,k,3) + c0w

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to Ye
!-----------------------------------------------------------------------

        da0wdy           = scefy (1,j,k,3) + pafy(1,j,k,3) + scnnfy(1,j,k,3) &
&                        + scnAfy(1,j,k,3) + bafy(1,j,k,3) + scnfy(1,j,k,3)
        db0wdy           = scefy (2,j,k,3) + pafy(2,j,k,3) + scnnfy(2,j,k,3) &
&                        + scnAfy(2,j,k,3) + bafy(2,j,k,3) + scnfy(2,j,k,3)
        dc0wdy           = scefy (3,j,k,3) + pafy(3,j,k,3) + scnnfy(3,j,k,3) &
&                        + scnAfy(3,j,k,3) + bafy(3,j,k,3) + scnfy(3,j,k,3)
        da1wdy           = scefy (4,j,k,3) + pafy(4,j,k,3) + scnnfy(4,j,k,3) &
&                        + scnAfy(4,j,k,3) + bafy(4,j,k,3) + scnfy(4,j,k,3)
        db1wdy           = scefy (5,j,k,3) + pafy(5,j,k,3) + scnnfy(5,j,k,3) &
&                        + scnAfy(5,j,k,3) + bafy(5,j,k,3) + scnfy(5,j,k,3)
        dc1wdy           = scefy (6,j,k,3) + pafy(6,j,k,3) + scnnfy(6,j,k,3) &
&                        + scnAfy(6,j,k,3) + bafy(6,j,k,3) + scnfy(6,j,k,3)

        RHS0_y(j,ne_2+k) = da0wdy * psi0(j,k,3) + db0wdy * psi1(j,k,3) + dc0wdy

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nue
!-----------------------------------------------------------------------

        DO kp = 1,nnugp(3)

          da0wdp0        = scefp0 (1,kp,j,k,3) + scnnfp0(1,kp,j,k,3) &
&                        + scnAfp0(1,kp,j,k,3) + scnfp0(1,kp,j,k,3)
          db0wdp0        = scefp0 (2,kp,j,k,3) + scnnfp0(2,kp,j,k,3) &
&                        + scnAfp0(2,kp,j,k,3) + scnfp0(2,kp,j,k,3)
          dc0wdp0        = scefp0 (3,kp,j,k,3) + scnnfp0(3,kp,j,k,3) &
&                        + scnAfp0(3,kp,j,k,3) + scnfp0(3,kp,j,k,3)
          da1wdp0        = scefp0 (4,kp,j,k,3) + scnnfp0(4,kp,j,k,3) &
&                        + scnAfp0(4,kp,j,k,3) + scnfp0(4,kp,j,k,3)
          db1wdp0        = scefp0 (5,kp,j,k,3) + scnnfp0(5,kp,j,k,3) &
&                        + scnAfp0(5,kp,j,k,3) + scnfp0(5,kp,j,k,3)
          dc1wdp0        = scefp0 (6,kp,j,k,3) + scnnfp0(6,kp,j,k,3) &
&                        + scnAfp0(6,kp,j,k,3) + scnfp0(6,kp,j,k,3)

          RHS0_p0(j,ne_2+k,ne_2+kp) = da0wdp0 * psi0(j,k,3) + db0wdp0 * psi1(j,k,3) + dc0wdp0

        END DO ! kp

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nuebar
!-----------------------------------------------------------------------

        IF ( nnugp(4) /= 0 ) THEN
          DO kp = 1,nnugp(4)

            da0wdp0      = pafp0(1,kp,j,k,3) + bafp0(1,kp,j,k,3)
            db0wdp0      = pafp0(2,kp,j,k,3) + bafp0(2,kp,j,k,3)
            dc0wdp0      = pafp0(3,kp,j,k,3) + bafp0(3,kp,j,k,3)
            da1wdp0      = pafp0(4,kp,j,k,3) + bafp0(4,kp,j,k,3)
            db1wdp0      = pafp0(5,kp,j,k,3) + bafp0(5,kp,j,k,3)
            dc1wdp0      = pafp0(6,kp,j,k,3) + bafp0(6,kp,j,k,3)

            RHS0_p0(j,ne_2+k,ne_3+kp) = da0wdp0 * psi0(j,k,3) + db0wdp0 * psi1(j,k,3) + dc0wdp0

          END DO ! kp
        END IF ! nnugp(4) /= 0

!-----------------------------------------------------------------------
!  Diagonal contribution to RHS0 (k' = k)
!-----------------------------------------------------------------------

        RHS0_p0(j,ne_2+k,ne_2+k) = RHS0_p0(j,ne_2+k,ne_2+k) + a0w

      END DO ! j
    END DO ! k
  END IF ! nnugp(3) /= 0

!-----------------------------------------------------------------------
!
!         \\\\\ GET SOURCE FUNCTIONS FOR X-ANTINEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(4) /= 0 ) THEN

!-----------------------------------------------------------------------
!  Compute and store combinations of x-antineutrino interaction moments
!   appearing in the p2 equation for psi0.
!
!        tot mfp^-1 : emis + abs - b1  (total inverse mean free path)
!-----------------------------------------------------------------------

    DO k = 1,nnugp(4)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  RHS0 x-antineutrino
!  (Right-hand side of the zero-moment transport equation for x-antineutrinos)
!-----------------------------------------------------------------------

        a0w              = scef (1,j,k,4) + paf(1,j,k,4) + scnnf(1,j,k,4) &
&                        + scnAf(1,j,k,4) + baf(1,j,k,4) + scnf(1,j,k,4)
        b0w              = scef (2,j,k,4) + paf(2,j,k,4) + scnnf(2,j,k,4) &
&                        + scnAf(2,j,k,4) + baf(2,j,k,4) + scnf(2,j,k,4)
        c0w              = scef (3,j,k,4) + paf(3,j,k,4) + scnnf(3,j,k,4) &
&                        + scnAf(3,j,k,4) + baf(3,j,k,4) + scnf(3,j,k,4)
        a1w              = scef (4,j,k,4) + paf(4,j,k,4) + scnnf(4,j,k,4) &
&                        + scnAf(4,j,k,4) + baf(4,j,k,4) + scnf(4,j,k,4)
        b1w              = scef (5,j,k,4) + paf(5,j,k,4) + scnnf(5,j,k,4) &
&                        + scnAf(5,j,k,4) + baf(5,j,k,4) + scnf(5,j,k,4)
        c1w              = scef (6,j,k,4) + paf(6,j,k,4) + scnnf(6,j,k,4) &
&                        + scnAf(6,j,k,4) + baf(6,j,k,4) + scnf(6,j,k,4)

        RHS0(j,ne_3+k) = a0w * psi0(j,k,4) + b0w * psi1(j,k,4) + c0w

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to Ye
!-----------------------------------------------------------------------

        da0wdy           = scefy (1,j,k,4) + pafy(1,j,k,4) + scnnfy(1,j,k,4) &
&                        + scnAfy(1,j,k,4) + bafy(1,j,k,4) + scnfy(1,j,k,4)
        db0wdy           = scefy (2,j,k,4) + pafy(2,j,k,4) + scnnfy(2,j,k,4) &
&                        + scnAfy(2,j,k,4) + bafy(2,j,k,4) + scnfy(2,j,k,4)
        dc0wdy           = scefy (3,j,k,4) + pafy(3,j,k,4) + scnnfy(3,j,k,4) &
&                        + scnAfy(3,j,k,4) + bafy(3,j,k,4) + scnfy(3,j,k,4)
        da1wdy           = scefy (4,j,k,4) + pafy(4,j,k,4) + scnnfy(4,j,k,4) &
&                        + scnAfy(4,j,k,4) + bafy(4,j,k,4) + scnfy(4,j,k,4)
        db1wdy           = scefy (5,j,k,4) + pafy(5,j,k,4) + scnnfy(5,j,k,4) &
&                        + scnAfy(5,j,k,4) + bafy(5,j,k,4) + scnfy(5,j,k,4)
        dc1wdy           = scefy (6,j,k,4) + pafy(6,j,k,4) + scnnfy(6,j,k,4) &
&                        + scnAfy(6,j,k,4) + bafy(6,j,k,4) + scnfy(6,j,k,4)

        RHS0_y(j,ne_3+k) = da0wdy * psi0(j,k,4) + db0wdy * psi1(j,k,4) + dc0wdy

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nue
!-----------------------------------------------------------------------

        DO kp = 1,nnugp(4)

          da0wdp0        = scefp0 (1,kp,j,k,4) + scnnfp0(1,kp,j,k,4) &
&                        + scnAfp0(1,kp,j,k,4) + scnfp0(1,kp,j,k,4)
          db0wdp0        = scefp0 (2,kp,j,k,4) + scnnfp0(2,kp,j,k,4) &
&                        + scnAfp0(2,kp,j,k,4) + scnfp0(2,kp,j,k,4)
          dc0wdp0        = scefp0 (3,kp,j,k,4) + scnnfp0(3,kp,j,k,4) &
&                        + scnAfp0(3,kp,j,k,4) + scnfp0(3,kp,j,k,4)
          da1wdp0        = scefp0 (4,kp,j,k,4) + scnnfp0(4,kp,j,k,4) &
&                        + scnAfp0(4,kp,j,k,4) + scnfp0(4,kp,j,k,4)
          db1wdp0        = scefp0 (5,kp,j,k,4) + scnnfp0(5,kp,j,k,4) &
&                        + scnAfp0(5,kp,j,k,4) + scnfp0(5,kp,j,k,4)
          dc1wdp0        = scefp0 (6,kp,j,k,4) + scnnfp0(6,kp,j,k,4) &
&                        + scnAfp0(6,kp,j,k,4) + scnfp0(6,kp,j,k,4)

          RHS0_p0(j,ne_3+k,ne_3+kp) = da0wdp0 * psi0(j,k,4) + db0wdp0 * psi1(j,k,4) + dc0wdp0

        END DO ! kp

!-----------------------------------------------------------------------
!  Derivatives of RHS0 with respect to psi0_nuebar
!-----------------------------------------------------------------------

        IF ( nnugp(3) /= 0 ) THEN
          DO kp = 1,nnugp(3)

            da0wdp0      = pafp0(1,kp,j,k,4) + bafp0(1,kp,j,k,4)
            db0wdp0      = pafp0(2,kp,j,k,4) + bafp0(2,kp,j,k,4)
            dc0wdp0      = pafp0(3,kp,j,k,4) + bafp0(3,kp,j,k,4)
            da1wdp0      = pafp0(4,kp,j,k,4) + bafp0(4,kp,j,k,4)
            db1wdp0      = pafp0(5,kp,j,k,4) + bafp0(5,kp,j,k,4)
            dc1wdp0      = pafp0(6,kp,j,k,4) + bafp0(6,kp,j,k,4)

            RHS0_p0(j,ne_3+k,ne_2+kp) = da0wdp0 * psi0(j,k,4) + db0wdp0 * psi1(j,k,4) + dc0wdp0

          END DO ! kp
        END IF ! nnugp(3) /= 0

!-----------------------------------------------------------------------
!  Diagonal contribution to RHS0 (k' = k)
!-----------------------------------------------------------------------

        RHS0_p0(j,ne_3+k,ne_3+k) = RHS0_p0(j,ne_3+k,ne_3+k) + a0w

      END DO ! j
    END DO ! k
  END IF ! nnugp(4) /= 0

!-----------------------------------------------------------------------
!
!             \\\\\ TRANSPORT LHS FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(1) /= 0 ) THEN
    DO k = 1,nnugp(1)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  LHS0 e-neutrino
!  (Left-hand side of the zero-moment transport equation for e-neutrinos)
!-----------------------------------------------------------------------

        LHS0(j,k)      = cdtinv(j) * ( psi0(j,k,1) - psi0_i(j,k,1) )                        &
&                      + a2vol_inv(j) * area(j  ) * agr2_inv(j  ) * psi1(j  ,k,1)           &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1) * psi1(j-1,k,1)

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_ph
!-----------------------------------------------------------------------

        IF ( j < jr_max ) THEN
        
          LHS0_p0p(j,k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                        &
&                      * ( psi1(j  ,k,1) * ddcpsjph(j  ,k,1)/ ( dc(j  ,k,1) + epsilon )     &
&                      - dc(j  ,k,1) * drjmh_inv(j  ) )

        ELSE ! j == jr_max

          LHS0_p0p(j,k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                         &
&                      * ( psi1(j  ,k,1) * ddcpsjph(j  ,k,1)/ ( dc(j  ,k,1) + epsilon ) )

        END IF ! j < jr_max

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_m3h
!-----------------------------------------------------------------------

        LHS0_p0m(j,k)  = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                         &
&                      * ( psi1(j-1,k,1) * ddcpsjmh(j-1,k,1)/ ( dc(j-1,k,1) + epsilon )     &
&                      + dc(j-1,k,1) * drjmh_inv(j-1) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_mh
!-----------------------------------------------------------------------

        IF ( j < jr_max ) THEN

          LHS0_p0(j,k) = cdtinv(j)                                                          &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( + dc(j  ,k,1) * drjmh_inv(j  )                                &
&                             + psi1(j  ,k,1) * ddcpsjmh(j  ,k,1)/( dc(j  ,k,1) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,1) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,1) * ddcpsjph(j-1,k,1)/( dc(j-1,k,1) + epsilon ) &
&                           )                                                               &
&                        )

        ELSE ! j == jr_max

          LHS0_p0(j,k) = cdtinv(j)                                                          &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( + psi_ratio(1)                                                &
&                             + psi1(j  ,k,1) * ddcpsjmh(j  ,k,1)/( dc(j  ,k,1) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,1) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,1) * ddcpsjph(j-1,k,1)/( dc(j-1,k,1) + epsilon ) &
&                           )                                                               &
&                        )

        END IF ! j < jr_max

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_ph
!-----------------------------------------------------------------------

        LHS0_yp(j,k)   = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                         &
&                      * psi1(j  ,k,1) * ddcdyjph(j  ,k,1)/ ( dc(j  ,k,1) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_m3h
!-----------------------------------------------------------------------

        LHS0_ym(j,k)   = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                         &
&                      * psi1(j-1,k,1) * ddcdyjmh(j-1,k,1)/ ( dc(j-1,k,1) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_mh
!-----------------------------------------------------------------------

        LHS0_y(j,k)    = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                         &
&                      * psi1(j  ,k,1) * ddcdyjmh (j  ,k,1)/ ( dc(j  ,k,1) + epsilon )      &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                           &
&                      * psi1(j-1,k,1) * ddcdyjph (j-1,k,1)/ ( dc(j-1,k,1) + epsilon )

      END DO ! j
    END DO ! k
  END IF ! nnugp(1) /= 0

!-----------------------------------------------------------------------
!
!            \\\\\ TRANSPORT LHS FOR E-ANTINEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(2) /= 0 ) THEN
    DO k = 1,nnugp(2)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  LHS0 e-antineutrino
!  (Left-hand side of the zero-moment transport equation for e-antineutrinos)
!-----------------------------------------------------------------------

        LHS0(j,ne_1+k) = cdtinv(j) * ( psi0(j,k,2) - psi0_i(j,k,2) )                        &
&                      + a2vol_inv(j) * area(j  ) * agr2_inv(j  ) * psi1(j  ,k,2)           &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1) * psi1(j-1,k,2)

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_ph
!-----------------------------------------------------------------------

        IF ( j < jr_max ) THEN

          LHS0_p0p(j,ne_1+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                   &
&                      * ( psi1(j  ,k,2) * ddcpsjph(j  ,k,2)/ ( dc(j  ,k,2) + epsilon )     &
&                      - dc(j  ,k,2) * drjmh_inv(j  ) )

        ELSE ! j == jr_max

          LHS0_p0p(j,ne_1+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                   &
&                      * ( psi1(j  ,k,2) * ddcpsjph(j  ,k,2)/ ( dc(j  ,k,2) + epsilon ) )

        END IF ! j < jr_max

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_m3h
!-----------------------------------------------------------------------

        LHS0_p0m(j,ne_1+k) = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                     &
&                      * ( psi1(j-1,k,2) * ddcpsjmh(j-1,k,2)/ ( dc(j-1,k,2) + epsilon )     &
&                      + dc(j-1,k,2) * drjmh_inv(j-1) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_mh
!-----------------------------------------------------------------------

        IF ( j < jr_max ) THEN

          LHS0_p0(j,ne_1+k) = cdtinv(j)                                                     &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( + dc(j  ,k,2) * drjmh_inv(j  )                                &
&                             + psi1(j  ,k,2) * ddcpsjmh(j  ,k,2)/( dc(j  ,k,2) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,2) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,2) * ddcpsjph(j-1,k,2)/( dc(j-1,k,2) + epsilon ) &
&                           )                                                               &
&                        )

        ELSE ! j = jr_max

          LHS0_p0(j,ne_1+k) = cdtinv(j)                                                     &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( psi_ratio(2)                                                  &
&                             + psi1(j  ,k,2) * ddcpsjmh(j  ,k,2)/( dc(j  ,k,2) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,2) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,2) * ddcpsjph(j-1,k,2)/( dc(j-1,k,2) + epsilon ) &
&                           )                                                               &
&                        )

        END IF ! j < jr_max

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_ph
!-----------------------------------------------------------------------

        LHS0_yp(j,ne_1+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                      &
&                      * psi1(j  ,k,2) * ddcdyjph(j  ,k,2)/ ( dc(j  ,k,2) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_m3h
!-----------------------------------------------------------------------

        LHS0_ym(j,ne_1+k) = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                      &
&                      * psi1(j-1,k,2) * ddcdyjmh(j-1,k,2)/ ( dc(j-1,k,2) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_mh
!-----------------------------------------------------------------------

        LHS0_y(j,ne_1+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                       &
&                      * psi1(j  ,k,2) * ddcdyjmh (j  ,k,2)/ ( dc(j  ,k,2) + epsilon )      &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                           &
&                      * psi1(j-1,k,2) * ddcdyjph (j-1,k,2)/ ( dc(j-1,k,2) + epsilon )

      END DO ! j
    END DO ! k
  END IF ! nnugp(2) /= 0

!-----------------------------------------------------------------------
!
!             \\\\\ TRANSPORT LHS FOR X-NEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(3) /= 0 ) THEN
    DO k = 1,nnugp(3)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  LHS0 x-neutrino
!  (Left-hand side of the zero-moment transport equation for x-neutrinos)
!-----------------------------------------------------------------------

        LHS0(j,ne_2+k) = cdtinv(j) * ( psi0(j,k,3) - psi0_i(j,k,3) )                        &
&                      + a2vol_inv(j) * area(j  ) * agr2_inv(j  ) * psi1(j  ,k,3)           &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1) * psi1(j-1,k,3)

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_ph
!-----------------------------------------------------------------------

        LHS0_p0p(j,ne_2+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                     &
&                      * ( psi1(j  ,k,3) * ddcpsjph(j  ,k,3)/ ( dc(j  ,k,3) + epsilon )     &
&                      - dc(j  ,k,3) * drjmh_inv(j  ) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_m3h
!-----------------------------------------------------------------------

        LHS0_p0m(j,ne_2+k) = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                     &
&                      * ( psi1(j-1,k,3) * ddcpsjmh(j-1,k,3)/ ( dc(j-1,k,3) + epsilon )     &
&                      + dc(j-1,k,3) * drjmh_inv(j-1) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_mh
!-----------------------------------------------------------------------

        IF ( j < jr_max ) THEN

          LHS0_p0(j,ne_2+k) = cdtinv(j)                                                     &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( + dc(j  ,k,3) * drjmh_inv(j  )                                &
&                             + psi1(j  ,k,3) * ddcpsjmh(j  ,k,3)/( dc(j  ,k,3) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,3) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,3) * ddcpsjph(j-1,k,3)/( dc(j-1,k,3) + epsilon ) &
&                           )                                                               &
&                        )

        ELSE ! j = jr_max

          LHS0_p0(j,ne_2+k) = cdtinv(j)                                                     &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( psi_ratio(3)                                                  &
&                             + psi1(j  ,k,3) * ddcpsjmh(j  ,k,3)/( dc(j  ,k,3) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,3) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,3) * ddcpsjph(j-1,k,3)/( dc(j-1,k,3) + epsilon ) &
&                           )                                                               &
&                        )

        END IF ! j < jr_max

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_ph
!-----------------------------------------------------------------------

        LHS0_yp(j,ne_2+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                      &
&                      * psi1(j  ,k,3) * ddcdyjph(j  ,k,3)/ ( dc(j  ,k,3) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_m3h
!-----------------------------------------------------------------------

        LHS0_ym(j,ne_2+k) = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                      &
&                      * psi1(j-1,k,3) * ddcdyjmh(j-1,k,3)/ ( dc(j-1,k,3) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_mh
!-----------------------------------------------------------------------

        LHS0_y(j,ne_2+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                       &
&                      * psi1(j  ,k,3) * ddcdyjmh (j  ,k,3)/ ( dc(j  ,k,3) + epsilon )      &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                           &
&                      * psi1(j-1,k,3) * ddcdyjph (j-1,k,3)/ ( dc(j-1,k,3) + epsilon )

      END DO ! j
    END DO ! k
  END IF ! nnugp(3) /= 0

!-----------------------------------------------------------------------
!
!            \\\\\ TRANSPORT LHS FOR X-ANTINEUTRINOS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(4) /= 0 ) THEN
    DO k = 1,nnugp(4)
      DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  LHS0 x-antineutrino
!  (Left-hand side of the zero-moment transport equation for x-antineutrinos)
!-----------------------------------------------------------------------

        LHS0(j,ne_3+k) = cdtinv(j) * ( psi0(j,k,4) - psi0_i(j,k,4) )                        &
&                      + a2vol_inv(j) * area(j  ) * agr2_inv(j  ) * psi1(j  ,k,4)           &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1) * psi1(j-1,k,4)

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_ph
!-----------------------------------------------------------------------

        LHS0_p0p(j,ne_3+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                     &
&                      * ( psi1(j  ,k,4) * ddcpsjph(j  ,k,4)/ ( dc(j  ,k,4) + epsilon )     &
&                      - dc(j  ,k,4) * drjmh_inv(j  ) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_m3h
!-----------------------------------------------------------------------

        LHS0_p0m(j,ne_3+k) = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                     &
&                      * ( psi1(j-1,k,4) * ddcpsjmh(j-1,k,4)/ ( dc(j-1,k,4) + epsilon )     &
&                      + dc(j-1,k,4) * drjmh_inv(j-1) )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to psi0_nue_mh
!-----------------------------------------------------------------------

        IF ( j < jr_max ) THEN

          LHS0_p0(j,ne_3+k) = cdtinv(j)                                                     &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( + dc(j  ,k,4) * drjmh_inv(j  )                                &
&                             + psi1(j  ,k,4) * ddcpsjmh(j  ,k,4)/( dc(j  ,k,4) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,4) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,4) * ddcpsjph(j-1,k,4)/( dc(j-1,k,4) + epsilon ) &
&                           )                                                               &
&                        )

        ELSE ! j = jr_max

          LHS0_p0(j,ne_3+k) = cdtinv(j)                                                     &
&                      + a2vol_inv(j)                                                       &
&                      * ( area(j  ) * agr2_inv(j  )                                        &
&                         * ( psi_ratio(4)                                                  &
&                             + psi1(j  ,k,4) * ddcpsjmh(j  ,k,4)/( dc(j  ,k,4) + epsilon ) &
&                            )                                                              &
&                       +  area(j-1) * agr2_inv(j-1)                                        &
&                         * ( + dc(j-1,k,4) * drjmh_inv(j-1)                                &
&                             - psi1(j-1,k,4) * ddcpsjph(j-1,k,4)/( dc(j-1,k,4) + epsilon ) &
&                           )                                                               &
&                        )

        END IF ! j < jr_max

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_ph
!-----------------------------------------------------------------------

        LHS0_yp(j,ne_3+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                      &
&                      * psi1(j  ,k,4) * ddcdyjph(j  ,k,4)/ ( dc(j  ,k,4) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_m3h
!-----------------------------------------------------------------------

        LHS0_ym(j,ne_3+k) = - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                      &
&                      * psi1(j-1,k,4) * ddcdyjmh(j-1,k,4)/ ( dc(j-1,k,4) + epsilon )

!-----------------------------------------------------------------------
!  Derivatives of LHS0 with respect to Y_mh
!-----------------------------------------------------------------------

        LHS0_y(j,ne_3+k) = + a2vol_inv(j) * area(j  ) * agr2_inv(j  )                       &
&                      * psi1(j  ,k,4) * ddcdyjmh (j  ,k,4)/ ( dc(j  ,k,4) + epsilon )      &
&                      - a2vol_inv(j) * area(j-1) * agr2_inv(j-1)                           &
&                      * psi1(j-1,k,4) * ddcdyjph (j-1,k,4)/ ( dc(j-1,k,4) + epsilon )

      END DO ! j
    END DO ! k
  END IF ! nnugp(4) /= 0

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION LHS /////
!
!-----------------------------------------------------------------------

  DO j = jr_min,jr_max
    LHS0(j,ne_4+1)   = cdtinv(j) * ( ye_new(j) - ye(j) )
    LHS0_y(j,ne_4+1) = cdtinv(j)
  END DO ! j

!-----------------------------------------------------------------------
!
!                \\\\\ LEPTON CONSERVATION RHS /////
!
!-----------------------------------------------------------------------

  IF ( nnugp(1) /= 0 ) THEN
    DO k = 1,nnugp(1)
      DO j = jr_min,jr_max
        RHS0(j,ne_4+1) = RHS0(j,ne_4+1) &
&                      - cf_ye(j,k) * ( emis(j,k,1,ij_ray,ik_ray)       &
&                      - ( emis(j,k,1,ij_ray,ik_ray) + absor(j,k,1,ij_ray,ik_ray) ) * psi0(j,k,1) )
      END DO ! j
    END DO ! k
  END IF ! nnugp(1) /= 0

  IF ( nnugp(2) /= 0 ) THEN
    DO k = 1,nnugp(2)
      DO j = jr_min,jr_max
        RHS0(j,ne_4+1) = RHS0(j,ne_4+1) &
&                      + cf_ye(j,k) * ( emis(j,k,2,ij_ray,ik_ray)       &
&                      - ( emis(j,k,2,ij_ray,ik_ray) + absor(j,k,2,ij_ray,ik_ray) ) * psi0(j,k,2) )
      END DO ! j
    END DO ! k
  END IF ! nnugp(2) /= 0

  IF ( nnugp(1) /= 0 ) THEN
    DO k = 1,nnugp(1)
      DO j = jr_min,jr_max
        RHS0_p0(j,ne_4+1,k) = cf_ye(j,k) * ( emis(j,k,1,ij_ray,ik_ray) + absor(j,k,1,ij_ray,ik_ray) )
        RHS0_y(j,ne_4+1) = RHS0_y(j,ne_4+1) &
&                        - cf_ye(j,k) * ( emisy(j,k,1) - ( emisy(j,k,1) + absory(j,k,1) ) * psi0(j,k,1) )
      END DO ! j
    END DO ! k
  END IF ! nnugp(1) /= 0

  IF ( nnugp(2) /= 0 ) THEN
    DO k = 1,nnugp(2)
      DO j = jr_min,jr_max
        RHS0_p0(j,ne_4+1,ne_1+k) = - cf_ye(j,k) * ( emis(j,k,2,ij_ray,ik_ray) + absor(j,k,2,ij_ray,ik_ray) )
        RHS0_y(j,ne_4+1) = RHS0_y(j,ne_4+1) &
&                        + cf_ye(j,k) * ( emisy(j,k,2) - ( emisy(j,k,2) + absory(j,k,2) ) * psi0(j,k,2) )
      END DO ! j
    END DO ! k
  END IF ! nnugp(2) /= 0

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

  DO k = 1,ne_4
    DO j = jr_min,jr_max
      A(j,k,k)       = A(j,k,k) + LHS0_p0(j,k)
    END DO ! j
  END DO ! k

!-----------------------------------------------------------------------
!  A: k-equations and border involving LHS and RHS dYw(j)
!-----------------------------------------------------------------------

  DO k = 1,ne_4+1
    DO j = jr_min,jr_max
      A(j,k,ne_4+1)  = A(j,k,ne_4+1) + LHS0_y(j,k) - RHS0_y(j,k)
    END DO ! j
  END DO ! k

!-----------------------------------------------------------------------
!  D_m: k-equations involving LHS dpsi0(j-3/2)
!-----------------------------------------------------------------------

  DO k = 1,ne_4
    DO j = jr_min,jr_max
      D_m(k,j)       = LHS0_p0m(j,k)
    END DO ! j
  END DO ! k

!-----------------------------------------------------------------------
!  D_p: k-equations involving LHS dpsi0(j+1/2)
!-----------------------------------------------------------------------

  DO k = 1,ne_4
    DO j = jr_min,jr_max
      D_p(k,j)       = LHS0_p0p(j,k)
    END DO ! j
  END DO ! k

!-----------------------------------------------------------------------
!  C: k-equations involving LHS and RHS constant terms
!-----------------------------------------------------------------------

  DO k = 1,ne_4+1
    DO j = jr_min,jr_max
      C(k,j)         = LHS0(j,k) - RHS0(j,k)
    END DO ! j
  END DO ! k

!-----------------------------------------------------------------------
!
!               \\\\\ CENTRAL BOUNDARY CONDITION /////
!
!            i-1
!           R         = delta
!            k,k',1           k,k'
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Impose symmetrical boundary conditions at the center
!-----------------------------------------------------------------------

  R(:,:,jr_min-1)     = zero

  DO k = 1,ne_4+1
    R(k,k,jr_min-1)   = one
  END DO ! k

!-----------------------------------------------------------------------
!
!            \\\\\ SOLVE FOR THE UP RECURRSION MATRICES /////
!
!                     N
!           u     =  Sum  R       u       + R
!            k,j     k'=1  k,k',j  k',j+1    k,N+1,j
!
!-----------------------------------------------------------------------

  DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Load coefficient matrix A_coef with coefficients of the unknowns
!   u(j,k)
!-----------------------------------------------------------------------

    A_coef(:,:)      = A(j,:,:)

!-----------------------------------------------------------------------
!  Add to A_coef matrix the previously determined recursion coefficients
!   relating the unknowns u(j-1,k) to the unknowns u(j,k)
!-----------------------------------------------------------------------

    DO kp = 1,ne_4
      DO k = 1,ne_4
        A_coef(k,kp) = A_coef(k,kp) + R(k,kp,j-1) * d_m(k,j)
      END DO ! k
    END DO ! kp

!-----------------------------------------------------------------------
!  Initialize the right-hand side matrix C_coef
!-----------------------------------------------------------------------

    C_coef           = zero

!-----------------------------------------------------------------------
!  Load right-hand side matrix C_coef with the constants terms in the
!   equations for the unknowns u(j,k)
!-----------------------------------------------------------------------

    C_coef(:,ne_4+1) = RHS0(j,:) - LHS0(j,:)

!-----------------------------------------------------------------------
!  Add to C_coef the recursion coefficients relating the unknowns
!   u(k,j-1) to the constants in the equations for u(k,j)
!-----------------------------------------------------------------------

    C_coef(:,ne_4+1) = C_coef(:,ne_4+1) - D_m(:,j) * R(:,ne_4+1,j-1)

!-----------------------------------------------------------------------
!  Load right-hand side matrix C_coef with the coefficients of the
!   unknowns u(k,j+1)
!-----------------------------------------------------------------------

    DO k = 1,ne_4
      C_coef(k,k)    = - D_p(k,j)
    END DO

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

  END DO ! j

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

  DO j = jr_max,jr_min,-1
    DO k = 1,ne_4+1
      DO kp = 1,ne_4
        u(k,j)       = u(k,j) + R(k,kp,j) * u(kp,j+1)
      END DO ! kp
      u(k,j)         = u(k,j) + R(k,ne_4+1,j)
    END DO ! k
  END DO ! j

!-----------------------------------------------------------------------
!
!      \\\\\ TRANSFER INVCREMENTS AND TEST FOR CONVERGENCE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer increments
!-----------------------------------------------------------------------

  DO j = jr_min,jr_max

    IF ( nnugp(1) /= 0 ) THEN
      DO k = 1,nnugp(1)
        dpsi0_iph(j,k,1) = u(k,j)
      END DO ! k
    END IF ! nnugp(1) /= 0
    
    IF ( nnugp(2) /= 0 ) THEN
      DO k = 1,nnugp(2)
        dpsi0_iph(j,k,2) = u(ne_1+k,j)
      END DO ! k
    END IF ! nnugp(2) /= 0
    
    IF ( nnugp(3) /= 0 ) THEN
      DO k = 1,nnugp(3)
        dpsi0_iph(j,k,3) = u(ne_2+k,j)
      END DO ! k
    END IF ! nnugp(3) /= 0
    
    IF ( nnugp(4) /= 0 ) THEN
      DO k = 1,nnugp(4)
        dpsi0_iph(j,k,4) = u(ne_3+k,j)
      END DO ! k
    END IF ! nnugp(4) /= 0

    dye_iph(j)       = u(ne_4+1,j)

  END DO ! j
 
!-----------------------------------------------------------------------
!        Increment variables
!
!  psi0(j,k,n) is allowed to exceed unity during iterations in order for
!   the iterations to converge.
!  psi0(j,k,n) is not allowed to fall below zero to avoid unphysical solutions
!
!-----------------------------------------------------------------------

  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN
      DO k = 1,nnugp(n)
        DO j = jr_min,jr_max
          psi0(j,k,n) = DMAX1( psi0(j,k,n) + dpsi0_iph(j,k,n), psi0amin )
        END DO ! j
      END DO ! k
    END IF ! nnugp(n) /= 0
  END DO ! n

  DO j = jr_min,jr_max
    ye_new(j)        = DMAX1( ye_new(j) + dye_iph(j), zero )
  END DO

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

  endit              = .true.
  jiterp0            = 0
  kiterp0            = 0
  niterp0            = 0
  jiterye            = 0

loop_checkpsi:  DO n = 1,nnu
    IF ( nnugp(n) /= 0 ) THEN
      DO k = 1,nnugp(n)
        DO j = jr_min,jr_max
          IF ( DABS( dpsi0_iph(j,k,n)/DMAX1( psi0(j,k,n), tolpsimin ) ) >= tolnupsi ) THEN
            jiterp0  = j
            kiterp0  = k
            niterp0  = n
            endit    = .false.
            EXIT loop_checkpsi
          END IF ! dpsi0_iph(j,k,n)/DMAX1( psi0(j,k,n) > tolnupsi
        END DO ! j
      END DO ! k
    END IF ! nnugp(n) /= 0
  END DO loop_checkpsi

loop_checkye:  DO j = jr_min,jr_max
    IF ( DABS( dye_iph(j)/DMAX1( ye_new(j), 0.1d0 ) ) >= tolnuye ) THEN
      jiterye        = j
      endit          = .false.
      EXIT loop_checkye
    END IF ! DABS( dye_iph(j)/DMAX1( ye(j), 0.1 ) ) >= tolnuye
  END DO loop_checkye

  IF ( .not. endit ) THEN
    IF ( jiterp0 == 0 ) THEN
      WRITE (nlog,1501) it, myid, ij_ray, ik_ray, jiterye, dye_iph(jiterye),                 &
&      dye_iph(jiterye)/DMAX1( ye_new(jiterye), 0.1d0 ), jiterye, ye_new(jiterye)
    ELSE IF ( jiterye == 0 ) THEN
      WRITE (nlog,1503) it, myid, ij_ray, ik_ray, jiterp0, kiterp0, niterp0,                 &
&      dpsi0_iph(jiterp0,kiterp0,niterp0),                                                   &
&      dpsi0_iph(jiterp0,kiterp0,niterp0)/DMAX1( psi0(jiterp0,kiterp0,niterp0), tolpsimin ), &
&      jiterp0, kiterp0, niterp0, psi0(jiterp0,kiterp0,niterp0)
    ELSE
      WRITE (nlog,1505) it, myid, ij_ray, ik_ray, jiterye, dye_iph(jiterye),dye_iph(jiterye)      &
&      /DMAX1( ye_new(jiterye), 0.1d0 ), jiterye, ye_new(jiterye), jiterp0, kiterp0,         &
&      niterp0, dpsi0_iph(jiterp0,kiterp0,niterp0),                                          &
&      dpsi0_iph(jiterp0,kiterp0,niterp0)/DMAX1( psi0(jiterp0,kiterp0,niterp0), tolpsimin ), &
&      jiterp0, kiterp0, niterp0, psi0(jiterp0,kiterp0,niterp0)
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

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      psi0(jr_min-1,k,n) = psi0(jr_min,k,n)
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!  Determine the largest relative change of psi0
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      DO j = jr_min,jr_max
        denu         = DABS( ( psi0(j,k,n) - psi0_i(j,k,n) ) )/( psi0(j,k,n) + psimin(n) )
        denu         = denu/DMAX1( 2.d-14 * rho(j), one )
        IF ( denu > dpsisp(n) ) THEN
          dpsisp(n)  = denu
          jdt(20+n,ij_ray,ik_ray) = j
        END IF ! denu > dpsisp(n)
      END DO ! j
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!
!                 \\\\\ UPDATE THE TEMPERATURE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute energy transferred to neutrinos in ith iteration
!-----------------------------------------------------------------------

de                   = zero

IF ( nnugp(1) /= 0 ) THEN
  DO k = 1,nnugp(1)
    DO j = jr_min,jr_max
      de(j)          = de(j) + cf_e(j,k) * stwt(1) * cdt(j) &
&                    * ( RHS0(j,k) + RHS0_y(j,k) * dye_iph(j) )
      DO kp = 1,nnugp(1)
        de(j)        = de(j) + cf_e(j,k) * stwt(1) * cdt(j) * RHS0_p0(j,k,     kp) * dpsi0_iph(j,kp,1)
      END DO
      IF ( nnugp(2) /= 0 ) THEN
        DO kp = 1,nnugp(2)
          de(j)      = de(j) + cf_e(j,k) * stwt(1) * cdt(j) * RHS0_p0(j,k,ne_1+kp) * dpsi0_iph(j,kp,2)
        END DO ! kp
      END IF ! nnugp(2) /= 0
      IF ( nnugp(3) /= 0 ) THEN
        DO kp = 1,nnugp(3)
          de(j)      = de(j) + cf_e(j,k) * stwt(1) * cdt(j) * RHS0_p0(j,k,ne_2+kp) * dpsi0_iph(j,kp,3)
        END DO ! kp
      END IF ! nnugp(3) /= 0
      IF ( nnugp(4) /= 0 ) THEN
        DO kp = 1,nnugp(4)
          de(j)      = de(j) + cf_e(j,k) * stwt(1) * cdt(j) * RHS0_p0(j,k,ne_3+kp) * dpsi0_iph(j,kp,4)
        END DO ! kp
      END IF ! nnugp(4) /= 0
    END DO ! j
  END DO ! k
END IF ! nnugp(1) /= 0

IF ( nnugp(2) /= 0 ) THEN
  DO k = 1,nnugp(2)
    DO j = jr_min,jr_max
      de(j)          = de(j) + cf_e(j,k) * stwt(2) * cdt(j) &
&                    * ( RHS0(j,ne_1+k) + RHS0_y(j,ne_1+k) * dye_iph(j) )
      IF ( nnugp(1) /= 0 ) THEN
        DO kp = 1,nnugp(1)
          de(j)      = de(j) + cf_e(j,k) * stwt(2) * cdt(j) * RHS0_p0(j,ne_1+k,     kp) * dpsi0_iph(j,kp,1)
        END DO ! 1,nnugp(1)
      END IF ! nnugp(1) /= 0
      DO kp = 1,nnugp(2)
        de(j)        = de(j) + cf_e(j,k) * stwt(2) * cdt(j) * RHS0_p0(j,ne_1+k,ne_1+kp) * dpsi0_iph(j,kp,2)
      END DO ! kp
     IF ( nnugp(3) /= 0 ) THEN
        DO kp = 1,nnugp(3)
          de(j)      = de(j) + cf_e(j,k) * stwt(2) * cdt(j) * RHS0_p0(j,ne_1+k,ne_2+kp) * dpsi0_iph(j,kp,3)
        END DO ! 1,nnugp(3)
      END IF ! nnugp(3) /= 0
      IF ( nnugp(4) /= 0 ) THEN
        DO kp = 1,nnugp(4)
          de(j)      = de(j) + cf_e(j,k) * stwt(2) * cdt(j) * RHS0_p0(j,ne_1+k,ne_3+kp) * dpsi0_iph(j,kp,4)
        END DO ! 1,nnugp(4)
      END IF ! nnugp(4) /= 0
   END DO ! j
  END DO ! k
END IF ! nnugp(2) /= 0

IF ( nnugp(3) /= 0 ) THEN
  DO k = 1,nnugp(3)
    DO j = jr_min,jr_max
      de(j)          = de(j) + cf_e(j,k) * stwt(3) * cdt(j) &
&                    * ( RHS0(j,ne_2+k) + RHS0_y(j,ne_2+k) * dye_iph(j) )
      IF ( nnugp(1) /= 0 ) THEN
        DO kp = 1,nnugp(1)
          de(j)      = de(j) + cf_e(j,k) * stwt(3) * cdt(j) * RHS0_p0(j,ne_2+k,     kp) * dpsi0_iph(j,kp,1)
        END DO ! kp
      END IF ! nnugp(1) /= 0
      IF ( nnugp(2) /= 0 ) THEN
        DO kp = 1,nnugp(2)
          de(j)      = de(j) + cf_e(j,k) * stwt(3) * cdt(j) * RHS0_p0(j,ne_2+k,ne_1+kp) * dpsi0_iph(j,kp,2)
        END DO ! kp
      END IF ! nnugp(2) /= 0
      DO kp = 1,nnugp(3)
        de(j)        = de(j) + cf_e(j,k) * stwt(3) * cdt(j) * RHS0_p0(j,ne_2+k,ne_2+kp) * dpsi0_iph(j,kp,3)
      END DO
      IF ( nnugp(4) /= 0 ) THEN
        DO kp = 1,nnugp(4)
          de(j)      = de(j) + cf_e(j,k) * stwt(3) * cdt(j) * RHS0_p0(j,ne_2+k,ne_3+kp) * dpsi0_iph(j,kp,4)
        END DO ! kp
      END IF ! nnugp(4) /= 0
    END DO ! j
  END DO ! k
END IF ! nnugp(3) /= 0

IF ( nnugp(4) /= 0 ) THEN
  DO k = 1,nnugp(4)
    DO j = jr_min,jr_max
      de(j)          = de(j) + cf_e(j,k) * stwt(4) * cdt(j) &
&                    * ( RHS0(j,ne_3+k) + RHS0_y(j,ne_3+k) * dye_iph(j) )
      IF ( nnugp(1) /= 0 ) THEN
        DO kp = 1,nnugp(1)
          de(j)      = de(j) + cf_e(j,k) * stwt(4) * cdt(j) * RHS0_p0(j,ne_3+k,     kp) * dpsi0_iph(j,kp,1)
        END DO ! kp
      END IF ! nnugp(1) /= 0
      IF ( nnugp(2) /= 0 ) THEN
        DO kp = 1,nnugp(2)
          de(j)      = de(j) + cf_e(j,k) * stwt(4) * cdt(j) * RHS0_p0(j,ne_3+k,ne_1+kp) * dpsi0_iph(j,kp,2)
        END DO ! kp
      END IF ! nnugp(2) /= 0
      IF ( nnugp(3) /= 0 ) THEN
        DO kp = 1,nnugp(3)
          de(j)      = de(j) + cf_e(j,k) * stwt(4) * cdt(j) * RHS0_p0(j,ne_3+k,ne_2+kp) * dpsi0_iph(j,kp,3)
        END DO
      END IF
      DO kp = 1,nnugp(4)
        de(j)        = de(j) + cf_e(j,k) * stwt(4) * cdt(j) * RHS0_p0(j,ne_3+k,ne_3+kp) * dpsi0_iph(j,kp,4)
      END DO ! kp
    END DO ! j
  END DO ! k
END IF ! nnugp(4) /= 0

!-----------------------------------------------------------------------
!  Compute the new temperature given by energy conservation
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t_new(j), ye_0(j), e, dedd, dedt, dedy )
  CALL e_degen(  rho(j), t_new(j), ye_0(j), e_eta )
  IF ( e_eta > degen ) THEN
    t_new(j)         = t_new(j) + ( - de(j) - dedy * ( ye_new(j) - ye_0(j) ) )/dedt
  ELSE
    CALL tgvndeye_x( j, ij_ray, ik_ray, rho(j), e - de(j), ye_new(j), t_new(j), t_new(j) )
  END IF ! e_eta > degen
END DO ! j

!-----------------------------------------------------------------------
!  Compute the rate of neutrino energy deposition
!-----------------------------------------------------------------------

dudt_nu(:,ij_ray,ik_ray)             = zero
dudt_nu(jr_min:jr_max,ij_ray,ik_ray) = - de(jr_min:jr_max)/( dtau_nutrns(jr_min:jr_max) + epsilon )

!-----------------------------------------------------------------------
!
!            \\\\\ UPDATE FUNCTIONS OF THE SOLUTION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Update dcr(j,k,n)
!-----------------------------------------------------------------------

IF ( idiff == 1  .or. idiff == 3 ) THEN
  dcr(:,:,:)         = dc(:,:,:)
END IF ! idiff == 1  .or. idiff == 3

!-----------------------------------------------------------------------
!  Update psi1..(Done last for consistent bookeeping)
!-----------------------------------------------------------------------

CALL psi1_cal( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, rho, t_new, ye_new, radius, &
& rstmss, u_vel, psi0, psi1, nx, nez, nnu, 2 )

!-----------------------------------------------------------------------
!  Compute the flux
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    CALL flux( jr_min, jr_max, n )
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!        Compute unujrad(j,n,ij_ray,ik_ray) and nnujrad(j,n,ij_ray,ik_ray).
!
!  unujrad(j,n,ij_ray,ik_ray) : total cumulative energy transported across
!                  radial zone j by n-type neutrinos (ergs)
!  nnujrad(j,n,ij_ray,ik_ray) : net number of n-type neutrinos transported
!                  across radial zone j
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      DO j = jr_min,jr_max
        unujrad(j,n,ij_ray,ik_ray) = unujrad(j,n,ij_ray,ik_ray) + fluxnuk(j,k,n) * area(j) * dtauj_nutrns(j)
        nnujrad(j,n,ij_ray,ik_ray) = nnujrad(j,n,ij_ray,ik_ray) + fluxnuk(j,k,n) * area(j) * dtauj_nutrns(j) &
&                          / ( ergmev * unue(j,k) )
      END DO ! j
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!        Compute dnurad(j,k,n,ij_ray,ik_ray).
!
!  dnurad(j,k,n,ij_ray,ik_ray) : cumulative number of n-type neutrinos
!   per unit energy that have crossed radial zone j in radial ray ij_ray, 
!  ik_ray (Mev^{-1})
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      DO j = jr_min,jr_max
        dnurad(j,k,n,ij_ray,ik_ray) = dnurad(j,k,n,ij_ray,ik_ray) + fluxnuk(j,k,n) * area(j) * dtauj_nutrns(j) &
&                           / ( ergmev * unue(j,k) * dunue(j,k) )
      END DO ! j
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!        Compute unukrad(k,n,ij_ray,ik_ray), nnukrad(k,n,ij_ray,ik_ray), 
!         and nnurad(n,ij_ray,ik_ray).
!
!  unukrad(k,n,ij_ray,ik_ray) : cumulative energy emitted from the core
!                  by n-type neutrinos of group k
!  nnukrad(k,n,ij_ray,ik_ray) : cumulative number of n-type neutrinos
!                  of energy k emitted by the core
!  nnurad(n,ij_ray,ik_ray)    : cumulative number of n-type neutrinos
!                  emitted by the core
!-----------------------------------------------------------------------

DO n = 1,nnu
  nnurad(n,ij_ray,ik_ray)        = zero
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      unukrad(k,n,ij_ray,ik_ray) = unukrad(k,n,ij_ray,ik_ray) + fluxnuk(jr_max,k,n) * area(jr_max) * dtauj_nutrns(jr_max)
      nnukrad(k,n,ij_ray,ik_ray) = unukrad(k,n,ij_ray,ik_ray)/( ergmev * unue(jr_max,k) )
      nnurad(n,ij_ray,ik_ray)    = nnurad(n,ij_ray,ik_ray) + nnukrad(k,n,ij_ray,ik_ray)
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!        Compute unulum(n) and unurad(n,ij_ray,ik_ray).
!
!  unulum(n) : core n-type neutrino luminosity (ergs/sec)
!  unurad(n) : cumulative energy emitted from the core in n-type
!               neutrinos (ergs)
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    unulum(n)                = fluxnu(jr_max,n) * area(jr_max)
    unurad(n,ij_ray,ik_ray)  = unurad(n,ij_ray,ik_ray) + unulum(n) * dtauj_nutrns(jr_max)
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!        Compute unujea(j,n), dunujeadt(j,n), and unucrea(n)
!
!  unujea(j,n)   : energy transferred to n-type neutrinos by emission
!                   and absorption in zone j (ergs)
! dunujeadt(j,n) : rate of energy transferred to n-type neutrinos by
!                   emission and absorption in zone j (ergs/s zone)
!  unucrea(n)    : total energy transferred to n-type neutrinos by
!                   emission and absorption (ergs)
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO j = jr_min,jr_max
    dunudt               = SUM( ecoefa(j,:) * ( emis(j,:,n,ij_ray,ik_ray) * ( one - psi0(j,:,n) ) &
&                        - absor(j,:,n,ij_ray,ik_ray) * psi0(j,:,n) ) * vol(j) * cvel )
    unujea(j,n,ij_ray,ik_ray)    = unujea(j,n,ij_ray,ik_ray) + dunudt *  dtau_nutrns(j)
    unucrea(n,ij_ray,ik_ray)     = unucrea(n,ij_ray,ik_ray) + dunudt * dtau_nutrns(j)
    dunujeadt(j,n,ij_ray,ik_ray) = dunudt
  END DO ! j
END DO ! n

!-----------------------------------------------------------------------
!
!                  \\\\\ DEALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

DEALLOCATE (LHS0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (LHS0_p0p, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0p  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (LHS0_p0m, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0m  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (LHS0_p0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_p0   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (LHS0_yp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_yp   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (LHS0_ym, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'LHS0_ym   '; WRITE (nlog,2001) var_name; END IF
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

DEALLOCATE (A, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'A         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (D_m, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'D_m       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (D_p, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'D_p       '; WRITE (nlog,2001) var_name; END IF
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

DEALLOCATE (a2vol_inv, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a2vol_inv '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (agr2_inv, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'agr2_inv  '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ddcpsjph, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcpsjph  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ddcpsjmh, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcpsjmh  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ddcdyjph, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcdyjph  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ddcdyjmh, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ddcdyjmh  '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE nu_trans
