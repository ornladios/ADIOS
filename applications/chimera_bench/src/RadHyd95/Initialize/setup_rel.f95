SUBROUTINE setup_rel( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         setup_rel
!    Module:       setup_rel
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!        To compute the GR variables dmgrv, gamgr, and dr
!
!    Subprograms called:
!      w_cal
!
!    Input arguments:
!  jr_min    : minimum radial index
!  jr_max    : maximum radial index
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module,
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nx
USE numerical_module, ONLY : half, one, frpith, third
USE physcnst_module, ONLY : g, cvel

USE edit_module, ONLY : nprint
USE eos_snc_x_module, ONLY : aesv
USE mdl_cnfg_module, ONLY : r, dr, u, rho, t, ye, dmrst, rstmss, gamgr, &
& dmgrv, grvmss, gamgrr, wgr

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: it            ! iteration index
INTEGER                          :: j             ! radial zone index

REAL(KIND=double)                :: csqinv        ! 1/c^{2}
REAL(KIND=double)                :: tgdc2         ! 2.*G/c^{2}
REAL(KIND=double)                :: rt            ! r(j) (temporary)
REAL(KIND=double)                :: drt           ! dr(j) (temporary)
REAL(KIND=double)                :: dmjmht        ! dmgrv(j) (temporary)
REAL(KIND=double)                :: gamgrt        ! gamgr(j) (temporary)
REAL(KIND=double)                :: rmasst        ! grvmss(j) (temporary)
REAL(KIND=double)                :: rmass         ! grvmss(j) (temporary)

REAL(KIND=double), PARAMETER     :: tol = 1.d-5

 4101 format (' The quantities dmgrv(j), dmrst(j), and gamgr(j) for zone',i3,' converged after',i3,' iterations')
 4103 format (' General relativistic quantities will not converge in subroutine genst')
 4105 format (' dmjmht=',1pe12.5,' dmjmh=',1pe12.5,' gamgrt=',1pe12.5,' gamgr(j)=',1pe12.5, &
& ' dmugrt=',1pe12.5,' dmrst(j)=',1pe12.5)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

csqinv             = 1.d0/cvel**2
tgdc2              = 2.d0 * g * csqinv

!-----------------------------------------------------------------------
!  Relativistic enthalpy
!-----------------------------------------------------------------------

CALL w_cal( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, wgr, nx )

!-----------------------------------------------------------------------
!  First approximation for gamgr(j)
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  r(j)             = r(j-1) + dr(j)
  gamgr(j)         = DSQRT( one + csqinv * ( u(j) * u(j) - 2.d0 * g * rstmss(j)/r(j) ) )
END DO

!-----------------------------------------------------------------------
!  Iterate to determine relativistic variables
!-----------------------------------------------------------------------

Radial: DO j = jr_min,jr_max

  DO it = 1,30

!........r(j)

    rt             = ( r(j-1)**3 + half * ( gamgr(j) + gamgr(j-1) ) * dmrst(j)/( frpith * rho(j) ) )**third

!........dr(j)

    drt            = rt - r(j-1)

!........dmgrv(j)

    dmjmht         = half * ( gamgr(j) + gamgr(j-1) ) * dmrst(j) * ( one + csqinv * aesv(j,2,ij_ray,ik_ray) )
    rmasst         = rmass + dmjmht

!........gamgr(j)

    gamgrt         = dsqrt( one + csqinv * ( u(j) * u(j) - 2.d0 * g * rmasst/rt ) )

!........Convergence tests..............................................

    IF (      it                        > 1               .and.    &
&             DABS( drt    - dr(j)    ) < tol * dr(j)     .and.    &
&             DABS( dmjmht - dmgrv(j) ) < tol * dmgrv(j)  .and.    &
&             DABS( gamgrt - gamgr(j) ) < tol * gamgr(j)        ) THEN

!........Successful convergence

      WRITE (nprint,4101) j,it

      dr(j)        = drt
      r(j)         = r(j-1) + drt
      rmass        = rmass + dmjmht
      dmgrv(j)     = dmjmht
      grvmss(j)    = rmass
      gamgr(j)     = gamgrt
      gamgrr(j)    = gamgrt

      CYCLE Radial

    END IF

    dmgrv(j)       = dmjmht
    gamgr(j)       = gamgrt
    dr(j)           = drt

  END DO ! it

!........Failure to converge

  WRITE (nprint,4103)
  WRITE (nprint,4105)
  STOP

END DO Radial

RETURN
END SUBROUTINE setup_rel
