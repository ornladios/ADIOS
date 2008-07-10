SUBROUTINE agr_cal( jr_min, jr_max, ij_ray, ik_ray, rho, p, r, agr, agrh, &
& nx )
!-----------------------------------------------------------------------
!
!    File:         agr_cal
!    Module:       agr_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/15/03
!
!    Purpose:
!      To compute the GR lapse function.
!
!    Variables that must be passed through common:
!  rho(j)      : matter density of zone j (g/cm**3)
!  r(j)        : radius of zone j (cm)
!  grvmss(j)   : gravitational mass enclosed by zone j (g)
!  gamgr(j)    : inverse radial metric
!  dmrst(j)    : rest mass betwee zone j and zone j-1
!  wgr(j)      : relativistic enthalpy
!
!    Subprograms called:
!  eqstt_x : interpolates quantities in the local EOS table
!
!    Input arguments:
!  jr_min      : inner zone for which calculation agr is to be made
!  jr_max      : outer zone for which calculation agr is to be made
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!  rho         : density (g cm^{-3})
!  p           : pressure (ergs cm^{-3})
!  r           : radial distance from center (cm)
!  nx          : x-array extent
!
!    Output arguments:
!  agr(j)      : shifted zone-edged lapse function
!  agrh(j)     : shifted zone-centered lapse function
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, mdl_cnfg.cmn, nu_dist.cmn, prb_cntl.cmn
!  tov_potential_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: half, one, frpi
USE physcnst_module, ONLY: g, cvel

USE edit_module, ONLY : nprint, nlog
USE mdl_cnfg_module, ONLY: grvmss, gamgr, wgr, dmrst
USE nu_dist_module, ONLY: stress_x
USE prb_cntl_module, ONLY: irelhy
USE tov_potential_module, ONLY : effpot_c, effpot_e

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: nx            ! radial array dimension

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rho      ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: p        ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: r        ! radius (cm)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx)   :: agr      ! shifted zone-edged lapse function
REAL(KIND=double), INTENT(out), DIMENSION(nx)   :: agrh     ! shifted zone-centered lapse function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: j             ! radial zone index
INTEGER                          :: jp1           ! jp1
INTEGER                          :: jmaxm         ! jr_max-1
INTEGER                          :: jminm         ! jr_min-1

REAL(KIND=double)                :: tgdc2         ! 2.*g/!**2
REAL(KIND=double)                :: coef          ! 2./(4*pi*c**2)
REAL(KIND=double)                :: csq           ! c**2
REAL(KIND=double)                :: wgrj          ! zone-edged relativistic enthalpy
REAL(KIND=double)                :: rhoj          ! zone-edged density
REAL(KIND=double)                :: dmrstj        ! zone-edged zone mass
REAL(KIND=double)                :: scoef         ! coefficient
REAL(KIND=double)                :: arg           ! argument of an exponential
REAL(KIND=double), EXTERNAL      :: fexp          ! exponential function

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  201 FORMAT (' irelhy =',i5,' not supported by agr_cal')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  tgdc2              = 2.d0 * g/cvel**2
  csq                = cvel**2
  coef               = 2.d0/csq
  first              = .false.
END IF

jmaxm                = jr_max - 1
jminm                = jr_min - 1

GR: SELECT CASE (irelhy)

!-----------------------------------------------------------------------
!
!                 \\\\\ NONRELATIVISTIC HYDRO /////
!
!-----------------------------------------------------------------------

 CASE(0) GR

  agr (jminm:jr_max)  = one
  agrh(jminm:jr_max)  = one
  
  RETURN

!-----------------------------------------------------------------------
!
!           \\\\\ 1ST POST-NEWTONIAN  RELATIVISTIC HYDRO /////
!
!-----------------------------------------------------------------------

 CASE(1) GR

  agr (jminm:jr_max)  = SQRT( one + coef * effpot_e(jminm:jr_max) )
  agrh(jr_min:jr_max) = SQRT( one + coef * effpot_c(jminm:jmaxm ) ) 
  
  RETURN

!-----------------------------------------------------------------------
!
!                  \\\\\ RELATIVISTIC HYDRO /////
!
!-----------------------------------------------------------------------

 CASE(2) GR

!-----------------------------------------------------------------------
!  Calculate agr at the outer boundary
!-----------------------------------------------------------------------

  agr(jr_max)         = ( one - tgdc2 * grvmss(jr_max)/r(jr_max) )/gamgr(jr_max)
  agrh(jr_max)        = agr(jr_max)

  IF ( jr_max > jr_min + 1 ) THEN

!-----------------------------------------------------------------------
!  Calculate agr at intermediate zone centers
!-----------------------------------------------------------------------

    DO j = jmaxm,jr_min,-1

      jp1             = j + 1
      wgrj            = half * ( wgr(j) + wgr(jp1) )
      rhoj            = half * ( rho(j) + rho(jp1) )
      dmrstj          = half * ( dmrst(j) + dmrst(jp1) )
      scoef           = dmrstj/( frpi * r(j) * r(j) * rhoj )
      arg             = ( p(jp1) - p(j) &
&                     - scoef * ( stress_x(j,1,ij_ray,ik_ray) + stress_x(j,2,ij_ray,ik_ray) &
&                     + stress_x(j,3,ij_ray,ik_ray) + stress_x(j,4,ij_ray,ik_ray) ) ) &
&                     / ( rhoj * wgrj * csq )
      agrh(j)         = agrh(jp1) * fexp(arg)

    END DO

  END IF

!-----------------------------------------------------------------------
!  Calculate agr at the inner-zone edged boundary
!-----------------------------------------------------------------------

  agr(jminm)          = agrh(jr_min) * fexp( half * arg )

!-----------------------------------------------------------------------
!  Calculate agr at intermediate zone edges
!-----------------------------------------------------------------------

  agr(jr_min:jmaxm)   = half * ( agrh(jr_min:jmaxm) + agrh(jr_min+1:jmaxm+1) )

  RETURN

!-----------------------------------------------------------------------
!  Stop if irelhy /= 0 and irelhy /= 1
!-----------------------------------------------------------------------

 CASE DEFAULT GR

  WRITE (nprint, 201) irelhy
  WRITE (nlog  , 201) irelhy
  STOP

END SELECT GR

END SUBROUTINE agr_cal
