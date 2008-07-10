SUBROUTINE editsc( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editsc
!    Module:       editsc
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/30/96
!
!    Purpose:
!      To edit entropy and chemical potential data.
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
!  iprint     : 0    - DO not print to print file.
!               ne 0 - print to print file.
!  nprint     : unit number of print file.
!  iplot      : 0    - DO not print to plot file.
!               ne 0 - print to plot file.
!  nplot      : unit number of plot file.
!  nedsc(i)   : editsc counter for data set i.
!  intdsc(i)  : number of cycles between edits of data set i.
!  idxdsc(i)  : edit jr_min, jr_max, and every idxdsc(i) radial zone between them for data set i.
!
!    Include files:
!      kind_module, physcnst_module
!      edit_module, eos_snc_x_module, mdl_cnfg_module, nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero
USE physcnst_module, ONLY : kmev, hbar, pi, rmu, cvel, dmnp, me

USE edit_module, ONLY : prnttest, nedsc, intdsc, idxesc, nprint, head
USE eos_bck_module, ONLY : se, sneu, sd, sh
USE eos_snc_x_module, ONLY : aesv, eos_i, nse
USE mdl_cnfg_module, ONLY : rho, t, ye
USE nu_dist_module, ONLY : stwt

USE el_eos_module, ONLY : ES, PS
USE eos_m4c_module, ONLY : BSOUT, BSNUC

IMPLICIT NONE
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
INTEGER                          :: jv            ! do index
INTEGER                          :: n_time        ! used for data-and-time

INTEGER, PARAMETER               :: n1 = 1
INTEGER, PARAMETER               :: n2 = 2
INTEGER, PARAMETER               :: n3 = 3

REAL(KIND=double)                :: scoef1        ! coefficient for calculating entropy of equl neutrinos
REAL(KIND=double)                :: scoef2        ! coefficient for calculating entropy of equl neutrinos

REAL(KIND=double)                :: s_mat         ! entropy/k per baryon of the matter
REAL(KIND=double)                :: s_tot         ! entropy/k per baryon sum of matter and all neutrino types
REAL(KIND=double)                :: s_e           ! entropy/k per baryon of electrons plus positrons
REAL(KIND=double)                :: s_d           ! entropy/k per baryon of drip nucleons plus heavy nucleus translation
REAL(KIND=double)                :: s_h           ! entropy/k per baryon of excited states of heavy nuclei
REAL(KIND=double)                :: s_rad         ! entropy/k per baryon of photons
REAL(KIND=double)                :: senu          ! entropy/k per baryon of e-neutrinos
REAL(KIND=double)                :: sanu          ! entropy/k per baryon of e-antineutrinos
REAL(KIND=double)                :: sxnu          ! entropy/k per baryon of x-neuitrinos
REAL(KIND=double)                :: seaeq         ! entropy/k per baryon sum of equilib e-neutrinos and e-antineutrinos
REAL(KIND=double)                :: sxeq          ! entropy/k per baryon sum of equilib m-neutrinos assuming

REAL(KIND=double)                :: cmpe          ! electron chemical potential (MeV)
REAL(KIND=double)                :: cmpn          ! neutron chemical potential (MeV)
REAL(KIND=double)                :: cmpp          ! proton chemical potential (MeV)
REAL(KIND=double)                :: cmp_np        ! neutron - proton chemical potential (MeV)
REAL(KIND=double)                :: cmpenu        ! e-neutrino chemical potential (MeV)
REAL(KIND=double)                :: etan          ! neutron chemical potential/kt
REAL(KIND=double)                :: etap          ! proton chemical potential/kt
REAL(KIND=double)                :: eta_np        ! neutron - proton chemical potential/kt
REAL(KIND=double)                :: etae          ! electron chemical potential/kt
REAL(KIND=double)                :: etaenu        ! e-neutrino chemical potential/kt

REAL(KIND=double)                :: v             ! dummy variable
REAL(KIND=double)                :: vd            ! dummy variable
REAL(KIND=double)                :: vt            ! dummy variable
REAL(KIND=double)                :: vy            ! dummy variable
REAL(KIND=double)                :: coefs         ! coefficient for computing equilibrium entropy

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (/)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
  101 FORMAT (53x,'Entropy data')
  103 FORMAT (52x,14('-')/)
  105 FORMAT ('   j   smat/b     tot s/b      se/b    sb-trns/b   sa-ex/b    srad/b     senu/b     sanu/b &
&    sxnu/b    (se+sa)eq  (sxnu)eq'/)
  107 FORMAT (1x,i4,11(1pe11.3))
  201 FORMAT (46x,'Chemical potential data')
  203 FORMAT (45x,25('-')/)
  205 FORMAT ('   j  n cm pot   p cm pot  n-p cm pot  e cm pot  enu cm pot    eta n      eta p     eta n-p &
&     eta e     eta enu'/)
  207 FORMAT (1x,i4,10(1pe11.3))
  209 FORMAT (1x,i4,11(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  scoef1           = ( 1.d+00/6.d+00 ) * ( kmev/ hbar * cvel )**3
  scoef2           = ( 7.d+00/15.d+00 ) * pi**2
END IF ! first
n_time             = nprint

!**********************************************************************!
!                                                                      !
!                  Entropy data                                        !
!                                                                      !
!----------------------------------------------------------------------!
! s_mat       'smat/b'     (Entropy/k per baryon of the matter)        !
! s_tot       'tot s/b'    (Entropy/k per baryon sum of matter and all !
!                           neutrino types)                            !
! s_e         'se/b'       (Entropy/k per baryon of electrons plus     !
!                           positrons)                                 !
! s_d         'sb-trns/b'  (Entropy/k per baryon of drip nucleons plus !
!                           heavy nucleus translation)                 !
! s_h         'sa-ex/b'    (Entropy/k per baryon of excited states of  !
!			                heavy nuclei)                              !
! s_rad       'srad/b'     (Entropy/k per baryon of photons)           !
! senu        'senu/b'     (Entropy/k per baryon of e-neutrinos)       !
! sanu        'sanu/b'     (Entropy/k per baryon of e-antineutrinos)   !
! sxnu        'sxnu/b'     (Entropy/k per baryon of m-neuitrinos)      !
! seaeq       '(se+sa)eq'  (Entropy/k per baryon sum of e-neutrinos    !
!                           and e-antineutrinos assuming thermal and   !
!                           beta equilibrium)                          !
! sxeq        '(sxnu)eq'   (Entropy/k per baryon sum of m-neutrinos    !
!                           assuming thermal and beta equilibrium      !
!----------------------------------------------------------------------!
!        Print header.                                                 !
!----------------------------------------------------------------------!

IF ( prnttest ) THEN
  nedsc(1)          = nedsc(1) + 1
  IF ( nedsc(1) < intdsc(1) ) GO TO 1100
  nedsc(1)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print(n_time)
WRITE (nprint,101)
WRITE (nprint,103)
WRITE (nprint,105)

!-----------------------------------------------------------------------
!  Entropies
!-----------------------------------------------------------------------

DO jv = jr_min,jr_max,idxesc(1)
  j                = jr_max - jv + jr_min

  se               = zero
  sneu             = zero
  sd               = zero
  sh               = zero

  ES               = zero
  PS               = zero
  BSOUT            = zero
  BSNUC            = zero

  CALL eqstta_x(3, j, ij_ray, ik_ray, rho(j), t(j), ye(j), v, vd, vt, vy )

  IF ( nse(j,ij_ray,ik_ray) == 0 ) THEN
    s_e            = se
    s_rad          = sneu
    s_d            = sd
    s_h            = sh
  ELSE
    s_e            = ES
    s_rad          = PS
    s_d            = BSOUT
    s_h            = BSNUC
  END IF

  s_mat            = aesv(j,3,ij_ray,ik_ray)
  CALL snu(j,n1,senu)
  CALL snu(j,n2,sanu)
  CALL snu(j,n3,sxnu)
  s_tot            = s_mat + senu + sanu + sxnu
  
  etaenu           = ( aesv(j,6,ij_ray,ik_ray) + aesv(j,5,ij_ray,ik_ray) &
&                  - aesv(j,4,ij_ray,ik_ray) )/( kmev * t(j) )

  coefs            = scoef1 * ( t(j)**3 ) * rmu/rho(j)
  seaeq            = coefs * ( scoef2 + etaenu * etaenu )
  sxeq             = coefs * scoef2 * ( stwt(3)/2.d+00 )

!-----------------------------------------------------------------------
!  Print entropy data
!-----------------------------------------------------------------------

  WRITE (nprint,107) j, s_mat, s_tot, s_e, s_d, s_h, s_rad, senu, sanu, &
&  sxnu, seaeq, sxeq

END DO

!**********************************************************************!
!                                                                      !
!                  Chemical potential data                             !
!                                                                      !
!----------------------------------------------------------------------!
! cmpn        'n cm pot'   (Neutron chemical potential (mev))          !
! cmpp        'p cm pot'   (Proton chemical potential (mev))           !
! cmp_np      'n-p cm pot' (Neutron chemical potential - proton        !
!                           chemical potential (mev))                  !
! cmpe        'e cm pot'   (Electron chemical potential (mev))         !
! cmpenu      'enu cm pot' (e-neutrino chemical potential assuming     !
!                           thermal and beta equilibrium (mev))        !
! etan        'eta n'      (Neutron chemical potential/kt)             !
! etap        'eta p'      (Proton chemical potential/kt)              !
! eta_np      'eta n-p'    ((Neutron chemical potential - proton       !
!                           chemical potential)/kt)                    !
! etae        'eta e'      (Electron chemical potential/kt)            !
! etaenu      'eta enu'    (e-neutrino chemical potential/kt)          !
!----------------------------------------------------------------------!
!        Print header.                                                 !
!----------------------------------------------------------------------!

 1100 CONTINUE

IF ( prnttest ) THEN
  nedsc(2)         = nedsc(2) + 1
  IF ( nedsc(2) < intdsc(2) ) RETURN
  nedsc(2)         = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print(n_time)
WRITE (nprint,201)
WRITE (nprint,203)
WRITE (nprint,205)
 
!-----------------------------------------------------------------------
!  Chemical potentials
!-----------------------------------------------------------------------

DO jv = jr_min,jr_max,idxesc(1)
  j                = jr_max - jv + jr_min

  cmpn             = aesv(j,4,ij_ray,ik_ray)
  cmpp             = aesv(j,5,ij_ray,ik_ray)
  cmp_np           = cmpn - cmpp
  cmpe             = aesv(j,6,ij_ray,ik_ray)
  cmpenu           = cmpe + ( cmpp - cmpn ) - dmnp

  etan             = cmpn/( kmev * t(j) )
  etap             = cmpp/( kmev * t(j) )
  eta_np           = etan - etap
  etae             = ( cmpe - me )/( kmev * t(j) )
  etaenu           = cmpenu/( kmev * t(j) )

!-----------------------------------------------------------------------
!  Print chemical potential data
!-----------------------------------------------------------------------

  WRITE (nprint,207) j, cmpn, cmpp, cmp_np, cmpe, cmpenu, etan, etap, &
&  eta_np, etae, etaenu

END DO

RETURN
END SUBROUTINE editsc
