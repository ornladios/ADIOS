SUBROUTINE editps( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editps
!    Module:       editps
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/25/96
!
!    Purpose:
!      To edit pressure and stress data.
!
!    Subprograms called:
!      date_and_time_print
!
!    Input arguments:
!  jr_min     : inner radial zone of region for which configuration edit is to be made.
!  jr_max     : outer radial zone of region for which configuration edit is to be made.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!  prnttest   : true  - test to see IF printing criteria is satisfied.
!               false - bypass printing criteria test.
!  iprint     : 0    - do not print to print file.
!               ne 0 - print to print file.
!  nprint     : unit number of print file.
!  iplot      : 0    - do not print to plot file.
!               ne 0 - print to plot file.
!  nplot      : unit number of plot file.
!  nedps(i)   : editps counter for data set i.
!  intdps(i)  : number of cycles between edits of data set i.
!  idxdps(i)  : edit jr_min, jr_max, and every idxdps(i) radial zone between them for data set i.
!
!    Include files:
!      kind_module, numerical_module
!      boundary_module, edit_module, eos_snc_x_module,
!      mdl_cnfg_module, nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : epsilon

USE boundary_module, ONLY : pbound
USE edit_module, ONLY : nprint, prnttest, nedps, idxeps, intdps, head, &
& gstrss, gstrss_cx, gstrss_cy, pstrss
USE eos_snc_x_module, ONLY : aesv
USE mdl_cnfg_module, ONLY : rho, r
USE nu_dist_module, ONLY : stress_nu=>stress_x, nu_str_ex, nu_str_cx

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

INTEGER                          :: j             ! radial zone index
INTEGER                          :: jv            ! do index
INTEGER                          :: n_time        ! iteration index

REAL(KIND=double)                :: strser        ! Ratio of e-neutrino stress to gravitational stress
REAL(KIND=double)                :: strsebr       ! Ratio of e-antineutrino stress to gravitational stress
REAL(KIND=double)                :: strsxr        ! Ratio of x-neutrino stress to gravitational gravitational stress
REAL(KIND=double)                :: strsxbr       ! Ratio of x-antineutrino stress to gravitational gravitational stress
REAL(KIND=double)                :: strsner       ! Ratio of zone edgeed total neutrino stress to gravitational gravitational stress
REAL(KIND=double)                :: strsncr       ! Ratio of zone center total neutrino stress to gravitational gravitational stress
REAL(KIND=double)                :: strspr        ! Ratio of material pressure stress to gravitational stress
REAL(KIND=double)                :: r2rho         ! (radius)**2*(density) (g/cm)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (/)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
  101 FORMAT (48x,'Pressure and stress data')
  103 FORMAT (47x,26('-')/)
  105 FORMAT ('   j      p        g-stress    g-str_xc    g-str_yc   e-stress/g &
& eb-strss/g  x-stress/g  xb-strss/g  ne-strss/g  nc-strss/g  p-stress/g    r2rho'/)
  107 FORMAT (1x,i4,12(es12.4))
  109 FORMAT (1x,' pbound=',es11.4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!**********************************************************************!        
!                                                                      !
!                  Pressure and stress data                            !
!                                                                      !
!----------------------------------------------------------------------!
! p           'p'          Pressure (dynes/cm**2)                      !
! gstrss(j)   'g-stress'   Gravitational stress (force/mass) in zone   ! 
!                           j (dynes/g)                                !
! gstrss_cx(j)'g-str_xc'   Gravitational x-stress (force/mass) in      ! 
!                           center of zone j (dynes/g)                 !
! gstrss_cy(j)'g-str_yc'   Gravitational y-stress (force/mass) in      ! 
!                           center of zone j (dynes/g)                 !
! strser      'e-stress/g' Ratio of e-neutrino stress to gravitational !
!                           stress                                     !
! strsebr     'eb-strss/g' Ratio of e-antineutrino stress to gravita-  !
!                           tional stress                              !
! strsxr      'x-stress/g' Ratio of x-neutrino stress to gravitational !
!                           stress                                     !
! strsxbr     'xb-strss/g' Ratio of x-antineutrino stress to gravita-  !
!                           tional stress                              !
! strsner    'ne-stress/g' Ratio of total zone-edgedn eutrino stress   !
!                           to gravitational stress                    !
! strsncr    'nc-stress/g' Ratio of total zone-center neutrino stress  !
!                           to gravitational stress                    !
! strspr      'p-stress/g' Ratio of material pressure stress to        !
!                           gravitational stress                       !
! r2rho       'r2rho'      (radius)**2*(density) (g/cm)                !
!----------------------------------------------------------------------!
!        Print header                                                  !
!----------------------------------------------------------------------!

n_time             = nprint
IF ( prnttest ) THEN
  nedps(1)         = nedps(1) + 1
  IF ( nedps(1) < intdps(1) ) RETURN
  nedps(1)         = 0
END IF !  prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print(n_time)
WRITE (nprint,101)
WRITE (nprint,103)
WRITE (nprint,105)

CALL stress( jr_min, jr_max, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Stress data
!-----------------------------------------------------------------------

DO jv = jr_min,jr_max,idxeps(1)
  j                = jr_max - jv + jr_min

  r2rho            = r(j) * r(j) * rho(j)
  strser           = stress_nu(j,1,ij_ray,ik_ray)/( gstrss(j) + epsilon )
  strsebr          = stress_nu(j,2,ij_ray,ik_ray)/( gstrss(j) + epsilon )
  strsxr           = stress_nu(j,3,ij_ray,ik_ray)/( gstrss(j) + epsilon )
  strsxbr          = stress_nu(j,4,ij_ray,ik_ray)/( gstrss(j) + epsilon )
  strsner          = nu_str_ex(j)                /( gstrss(j) + epsilon )
  strsncr          = nu_str_cx(j)                /( gstrss(j) + epsilon )
  strspr           = pstrss(j)                   /( gstrss(j) + epsilon )

!-----------------------------------------------------------------------
!  Print stress data
!-----------------------------------------------------------------------

  WRITE (nprint,107) j, aesv(j,1,ij_ray,ik_ray), gstrss(j), gstrss_cx(j), &
&  gstrss_cy(j), strser, strsebr, strsxr, strsxbr, strsner, strsncr, &
&  strspr, r2rho

END DO

WRITE (nprint,109) pbound

RETURN
END SUBROUTINE editps
