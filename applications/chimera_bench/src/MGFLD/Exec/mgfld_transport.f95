SUBROUTINE mgfld_transport( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, nx, &
& nez, nnu, jdt, dtnph_trans )
!-----------------------------------------------------------------------
!
!    File:         mgfld_transport
!    Module:       mgfld_transport
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/28/00
!
!    Purpose:
!      To execute the transort code after program initialization
!
!    Subprograms called:
!  extrap       : extrapolates variables
!  w_cal        : computes the relativistic enthalpy
!  gamgr_nu_cal : computes the relativistic gammas for transport
!  eqstz_x      : updates the thermodynamic variables
!  nu_adv       : performs the neutrino transport step
!  gammaz_x     : computes the adiabatic gammas
!  nu_number    : computes the neutrino number and energy
!
!    Input arguments:
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  nx           : x_array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  dtnph_trans  : source and transport time step
!
!    Output arguments:
!  jdt          : radial zone causing minimum time step for criteria i, radial ray
!                 ij_ray, ik_ray
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!  nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE cycle_module, ONLY: nutrans_trns
USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY: aesv
USE mdl_cnfg_module, ONLY: jr_min, jr_max, r, u, rho, rstmss, t, gamgrr, &
& gamgr, ye, wgr, wgrr, dmrst
USE nu_dist_module, ONLY: unujcr, rhonur, rhonun, apnur, apnun, apmnur, &
& apmnun, psi0, psi1

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nx            ! radial array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in)    :: dtnph_trans   ! source and transport time step

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)  :: jdt   ! zone causing dt

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: jmaxp           ! jr_max+1
INTEGER                          :: n               ! neutrino flavor index

INTEGER                          :: istat           ! allocation status

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: unun   ! neutrino energy/volume

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in mgfld_transport')
 2001 FORMAT (' Deallocation problem for array ',a10,' in mgfld_transport')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ ALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (unun(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unun      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                \\\\\ INITIALIZE FOR TRANSPORT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

jmaxp                   = jr_max + 1
unun                    = zero

!-----------------------------------------------------------------------
!  Extrapolate independent thermodynamic variables
!-----------------------------------------------------------------------

CALL extrap

!-----------------------------------------------------------------------
!  Add one ghost zone to independent thermodnamic variables
!-----------------------------------------------------------------------

rho (jr_max+1)          = rho (jr_max)
t   (jr_max+1)          = t   (jr_max)
ye  (jr_max+1)          = ye  (jr_max)

!-----------------------------------------------------------------------
!  Update the enthalpies
!-----------------------------------------------------------------------

CALL w_cal( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, wgr, nx )

!-----------------------------------------------------------------------
!  Relativistic gammas
!-----------------------------------------------------------------------

CALL gamgr_nu_cal( jr_min, jr_max )

!-----------------------------------------------------------------------
!
!                 \\\\\ SOURCE AND TRANSPORT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Execute source and transport
!-----------------------------------------------------------------------

CALL nu_adv( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, rho, &
& t, ye, r, rstmss, nx, nez, nnu, jdt, dtnph_trans )

!-----------------------------------------------------------------------
!
!            \\\\\ SWITCH ARRAYS AND UPDATE QUANTITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Switch arrays
!-----------------------------------------------------------------------

wgrr  (jr_min:jr_max)   = wgr  (jr_min:jr_max)
gamgrr(jr_min:jr_max)   = gamgr(jr_min:jr_max)

IF ( nutrans_trns) THEN

  unun(jr_min:jr_max)   = zero

  DO n = 1,nnu
    CALL nu_number( jr_min, jr_max, n, ij_ray, ik_ray, nx, nez, nnu, r, &
&    u, psi0, psi1 )
    unun(jr_min:jr_max) = unun(jr_min:jr_max) + unujcr(jr_min:jr_max,n)/dmrst(jr_min:jr_max)
  END DO

  rhonur(jr_min:jr_max) = rhonun(jr_min:jr_max)
  rhonun(jr_min:jr_max) = rho   (jr_min:jr_max)
  apnur (jr_min:jr_max) = apnun (jr_min:jr_max)
  apnun (jr_min:jr_max) = rho   (jr_min:jr_max) * unun(jr_min:jr_max)/3.d0
  apmnur(jr_min:jr_max) = apmnun(jr_min:jr_max)
  apmnun(jr_min:jr_max) = aesv  (jr_min:jr_max,1,ij_ray,ik_ray)

END IF ! nutrans = true

rho(jmaxp)              = rho(jr_max)
t  (jmaxp)              = t  (jr_max)
ye (jmaxp)              = ye (jr_max)

!-----------------------------------------------------------------------
!  Recompute thermodynamic quantities
!-----------------------------------------------------------------------

CALL eqstz_x( jr_min, jmaxp, rho, t, ye, ij_ray, ik_ray )
CALL gammaz_x( jr_min, jr_max, rho, t, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!
!                  \\\\\ REALLOCATE ARRAYS /////
!
!-----------------------------------------------------------------------

DEALLOCATE (unun, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unun      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE mgfld_transport
