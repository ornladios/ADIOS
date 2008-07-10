SUBROUTINE stress( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         stress
!    Module:       stress
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/23/04
!
!    Purpose:
!      To calculate the various stresses.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min    : inner zone for which calculation agr is to be made
!  jr_max    : outer zone for which calculation agr is to be made
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!        none
!
!    Output arguments (common):
!  gstrss    : force per unit mass due to gravity
!  nustrss   : orce per unit mass due to neutrinos of all types
!  pstrss    : force per unit mass due to material pressure
!  rstrss    : force per unit mass due to relativistic terms
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module, nu_dist_module
!  prb_cntl_module, shock_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nx
USE numerical_module, ONLY : zero, half, sxtnpi
USE physcnst_module, ONLY : pi, cvel, g

USE edit_module, ONLY : nprint, gstrss, nustrss, pstrss, rstrss
USE eos_snc_x_module, ONLY : aesv
USE mdl_cnfg_module, ONLY : r, gamgr, wgr, dmrst, grvmss, agr
USE nu_dist_module, ONLY : stress_nu=>stress_x, apnu
USE prb_cntl_module, ONLY : irelhy
USE shock_module, ONLY : pq_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: j             ! radial zone index
INTEGER                          :: jp1           ! j+1

REAL(KIND=double)                :: rj2           ! r^{2}
REAL(KIND=double)                :: tgdc2         ! 2.*g/c^{2}
REAL(KIND=double)                :: tpigc2        ! 2.*pi*g/c^{2}

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------


IF ( first ) THEN
  first            = .false.
  tgdc2            = 2.d0 * g/cvel**2
  IF ( irelhy == 0  .or.  irelhy == 1 ) THEN
    tpigc2         = zero
  ELSE
    tpigc2         = pi * tgdc2
  END IF ! irelhy = 0
END IF ! first

!-----------------------------------------------------------------------
!  Forces at zone j and time n
!-----------------------------------------------------------------------

DO j = jr_min,jr_max-1
  jp1              = j + 1
  rj2              = r(j) * r(j)
  nustrss(j)       = ( stress_nu(j,1,ij_ray,ik_ray) + stress_nu(j,2,ij_ray,ik_ray) &
&                  +   stress_nu(j,3,ij_ray,ik_ray) + stress_nu(j,4,ij_ray,ik_ray) ) * gamgr(j) &
&                  / ( half * ( wgr(jp1) + wgr(j) ) )
  pstrss(j)        = - sxtnpi * rj2 * gamgr(j) &
&                  * ( ( aesv(jp1,1,ij_ray,ik_ray) + pq_x(jp1,ij_ray,ik_ray) &
&                  - aesv(j,1,ij_ray,ik_ray) - pq_x(j,ij_ray,ik_ray) ) &
&                  / ( ( wgr(jp1) + wgr(j) ) * ( dmrst(jp1) + dmrst(j) ) ) )
  IF ( irelhy > 1 ) THEN
    gstrss(j)      = - g * grvmss(j)/rj2
  END IF
  rstrss(j)        = - tpigc2 * ( aesv(jp1,1,ij_ray,ik_ray) + pq_x(jp1,ij_ray,ik_ray) + apnu(jp1) &
&                  + aesv(j,1,ij_ray,ik_ray) + pq_x(j,ij_ray,ik_ray) + apnu(j) ) * r(j)
END DO

nustrss(jr_max)    = nustrss(jr_max-1)
pstrss(jr_max)     = pstrss(jr_max-1)
gstrss(jr_max)     = gstrss(jr_max-1)
rstrss(jr_max)     = rstrss(jr_max-1)

RETURN
END SUBROUTINE stress
