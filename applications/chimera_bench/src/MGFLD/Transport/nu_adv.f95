SUBROUTINE nu_adv( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& rho, t, ye, r, rstmss, nx, nez, nnu, jdt, dtnph_trans )
!-----------------------------------------------------------------------
!
!    File:         nu_adv
!    Module:       nu_adv
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/29/01
!
!    Purpose:
!      To advance advance the psi0(j,k,n). This subroutine links
!       together the following subroutines:
!
!    Subprograms called:
!  pre_trans    : computes some of the quantities needed in nu_abemtr
!  nu_abemtr    : computes the change in n-neutrinos occupation numbers
!                 due to absorption, emission, and transport
!  nu_sphere    : computes the n-neutrino spheres, needed for pair and brem
!  nu_scat      : computes the change in n-neutrinos occupation numbers
!                 due to inelastic scattering and pair production
!
!    Input arguments:
!
!  jr_min       : minimum radial zone index
!  jr_max       : maximum radial zone index
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  rho          : density (g cm^{-3})
!  t            : temperature (K)
!  ye           : electron fraction
!  r            : radius (cm)
!  rstmss       : enclosed rest mass (g)
!  nx           : x-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  dtnph_trans  : source and transport time step
!
!    Output arguments:
!  jdt          : radial zone causing minimum time step for criteria i, 
!                  radial ray ij_ray, ik_ray
!
!    Input arguments (common):
!
!  jr_min       : inner radial zone for neutrino transport
!  jr_max       : outer radial zone for neutrino transport
!  agr_nu(j)    : value of the lapse function at time m
!  rho(j)       : value of the density at the end of the hydro cycle
!  gamgra(j)    : value of gamgr(j) at the end of the hydro cycle
!
!    Output arguments (common):
!
!  psi0(j,k,n)  : updated zero moment neutrino occupation numbers
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, edit_module, eos_snc_x_module, incrmnt_module,
!  nucbrn_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one

USE cycle_module, ONLY : nutrans_trns
USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY : nse, nuc_number
USE incrmnt_module, ONLY : dye, dtmpnn
USE nucbrn_module, ONLY: a_name
USE nu_dist_module, ONLY : j_sphere, r_sphere, d_sphere, t_sphere, &
& m_sphere, xn
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : inutrn

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nx            ! radial array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rho         ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: r           ! radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rstmss      ! radius (cm)
REAL(KIND=double), INTENT(in)                   :: dtnph_trans ! source and transport time step

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------
 
INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim) :: jdt    ! zone causing dt

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx) :: t            ! temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(nx) :: ye           ! electron fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

LOGICAL                          :: first = .true.

INTEGER                          :: i             ! composition index
INTEGER                          :: j             ! shifted radial index
INTEGER                          :: istat         ! allocation status
INTEGER                          :: i_n           ! neutron abundance index
INTEGER                          :: i_p           ! proton abundance index

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: t_new    ! final temperature (K)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ye_new   ! final electron fraction (K)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in nu_adv')
 2001 FORMAT (' Deallocation problem for array ',a10,' in nu_adv')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if the problem does not involve neutrinos of any type
!-----------------------------------------------------------------------

dye   (jr_min:jr_max,1,ij_ray,ik_ray)   = zero
dtmpnn(jr_min:jr_max,1,ij_ray,ik_ray)   = zero

IF ( nnugpmx == 0  .or.  inutrn == 0 ) RETURN         

!-----------------------------------------------------------------------
!
!                   \\\\\ ALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (t_new(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_new     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_new(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_new    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!               \\\\\ PRELIMINARIES FOR TRANSPORT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Composition indexing
!-----------------------------------------------------------------------

IF ( first ) THEN
  first                         = .false.
  i_n                           = nuc_number + 1
  i_p                           = nuc_number + 1
  DO i = 1,nuc_number
    IF (      a_name(i) == '  n  ' ) THEN
      i_n                 = i
    END IF ! a_name(i) == '  n  '
    IF (      a_name(i) == '  p  ' ) THEN
      i_p                 = i
    END IF ! a_name(i) == '  p  '
  END DO ! i = 1,nuc_number
END IF ! first

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

t_new                           = t
ye_new                          = ye

!-----------------------------------------------------------------------
!
!         ||||| EXECUTE TRANSPORT IF NUTRANS_TRNS = TRUE |||||
!
!-----------------------------------------------------------------------

IF ( nutrans_trns ) THEN

!-----------------------------------------------------------------------
!  Calculate quantities needed for the neutrino transport
!-----------------------------------------------------------------------

  CALL pre_trans( jr_min, jr_max, rho, r, nx, nnu )

!-----------------------------------------------------------------------
!  Calculate neutrino-spheres
!-----------------------------------------------------------------------

  CALL nu_sphere( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&  r, rho, t, rstmss,  nx, nez, nnu, j_sphere, r_sphere, d_sphere, t_sphere, &
&  m_sphere )

!-----------------------------------------------------------------------
!
!                      \\\\\ TRANSPORT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Advance neutrino occupation probabilities due to source and
!   transport
!  Input temperatures and electron fractions are in t and ye.
!  Updated temperatures and electron fractions are in t_new and ye_new.
!-----------------------------------------------------------------------

  CALL nu_trans( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&  rho, t, t_new, ye, ye_new, r, rstmss, nx, nez, nnu, jdt, dtnph_trans )

!-----------------------------------------------------------------------
!
!                   \\\\\ INCREMENT T AND YE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Increment dye(j,1,ij_ray,ik_ray) and dtmpnn(j,1,ij_ray,ik_ray)
!-----------------------------------------------------------------------

  dye   (jr_min:jr_max,1,ij_ray,ik_ray) = ye_new(jr_min:jr_max) - ye(jr_min:jr_max)
  dtmpnn(jr_min:jr_max,1,ij_ray,ik_ray) = t_new(jr_min:jr_max) - t(jr_min:jr_max)
  dye   (jr_min:jr_max,1,ij_ray,ik_ray) = dye   (jr_min:jr_max,1,ij_ray,ik_ray)         &
&                               / DMAX1( 2.d-14 * rho(jr_min:jr_max), one )
  dtmpnn(jr_min:jr_max,1,ij_ray,ik_ray) = dtmpnn(jr_min:jr_max,1,ij_ray,ik_ray)         &
&                               / DMAX1( 2.d-14 * rho(jr_min:jr_max), one )
  ye(jr_min:jr_max)             = ye_new(jr_min:jr_max)
  t(jr_min:jr_max)              = t_new(jr_min:jr_max)

!-----------------------------------------------------------------------
!  Increment neutron and proton fractions in nonNSE material
!-----------------------------------------------------------------------

  DO j = jr_min, jr_max
    IF ( nse(j,ij_ray,ik_ray) == 0 ) THEN
      xn(j,i_n)                 = DMAX1( xn(j,i_n) - dye(j,1,ij_ray,ik_ray), zero )
      xn(j,i_p)                 = DMAX1( xn(j,i_p) + dye(j,1,ij_ray,ik_ray), zero )
    END IF ! nse(j,i_ray) == 0
  END DO ! j = jr_min, jr_max

END IF ! nutrans_trns

!-----------------------------------------------------------------------
!
!                 \\\\\ DEALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

DEALLOCATE (t_new, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_new     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye_new, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_new    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE nu_adv
