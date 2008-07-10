SUBROUTINE edity( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         edity
!    Module:       edity
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/2/00
!
!    Purpose:
!      To edit composition data.
!
!    Subprograms called:
!  date_and_time_print : prints date and time
!  w_cal               : computes relativistic gammas
!  flux                : computes neutrino fluxes
!  lectron             : computes the electron EOS quantities
!  eosnuc_x            : computes the non-NES EOS quantities
!  nuc_energy          : computes functions of the nuclear abundance distribution
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
!  nprint     : unit number of print file.
!  nplot      : unit number of plot file.
!  nedy(i)    : edity counter for data set i.
!  intedy(i)  : number of cycles between edits of data set i.
!  idxedy(i)  : edit jr_min, jr_max, and every idxedy(i) radial zone
!                between them for data set i.
!
!    Include files:
!  array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, eos_bck_module, eos_snc_x_module,
!  incrmnt_module, mdl_cnfg_module, nucbrn_module, nu_dist_module,
!  nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nnu
USE numerical_module, ONLY : zero, one, ncoef, epsilon
USE physcnst_module, ONLY : cm3fm3, rmu, kmev, ergfoe

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nlog, prnttest, nedy, intedy, head, idxedy
USE eos_bck_module, ONLY : dbck, tbck, yebck, yeplus, jshel
USE eos_snc_x_module, ONLY : aesv, nse
USE incrmnt_module, ONLY : dtmpmn
USE mdl_cnfg_module, ONLY : r, rho, t, ye, dmrst, ye0, wgr
USE nucbrn_module, ONLY : xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, dudt_nuc, &
& uburn, nuc_number, a_name, a_nuc, z_nuc
USE nu_dist_module, ONLY : unu, dunu, psi0
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY :  dtnmh

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

CHARACTER (len=10)               :: var_name

LOGICAL                          :: first = .true.
LOGICAL                          :: print_sub = .false.

INTEGER                          :: i             ! do index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: jv            ! do index
INTEGER                          :: l             ! nuclear specie index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: n_time        ! iteration index
INTEGER                          :: n_nucp1       ! nuc_number + 1
INTEGER                          :: i_mlt         ! number of times to repeat composition write
INTEGER                          :: i_rmn         ! remainder
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: kfm           ! ( # nucleons/gram )( cm3/fm3 )
REAL(KIND=double)                :: yl0           ! initial lepton number

REAL(KIND=double)                :: yn            ! neutron number per baryon
REAL(KIND=double)                :: yp            ! proton number per baryon
REAL(KIND=double)                :: yhe           ! helium number per baryon
REAL(KIND=double)                :: ynuc          ! heavy nucleus number per baryon
REAL(KIND=double)                :: yenu          ! e-neutrino number per baryon
REAL(KIND=double)                :: yanu          ! e-antineutrino number per baryon
REAL(KIND=double)                :: yxnu          ! x-neutrino number per baryon
REAL(KIND=double)                :: yl            ! lepton number per baryon
REAL(KIND=double)                :: ratlep        ! Ratio of current lepton number to original lepton number

REAL(KIND=double)                :: xnt           ! neutron mass fraction
REAL(KIND=double)                :: xp            ! proton mass fraction
REAL(KIND=double)                :: xhe           ! helium mass fraction 
REAL(KIND=double)                :: xnuc          ! heavy nucleus mass fraction
REAL(KIND=double)                :: a             ! heavy nucleus mass number
REAL(KIND=double)                :: rn            ! heavy nucleus neutron number
REAL(KIND=double)                :: z             ! heavy nucleus proton number
REAL(KIND=double)                :: ye_m          ! total number of electrons/nucleon
REAL(KIND=double)                :: ye_p          ! total number of positrons/nucleon 
REAL(KIND=double)                :: ye_pm         ! ratio of positrons to electrons 
REAL(KIND=double)                :: ypairye       ! ratio of pairs to net electrons
REAL(KIND=double)                :: xn_tot        ! sum of themass fractions in a zone
REAL(KIND=double)                :: enb_mev       ! mean binding energy per particle
REAL(KIND=double)                :: enm           ! rest mass energy

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rnnu   ! neutrino number density
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: utburn ! Total energy released by nuclear burning interior to zone j (foes)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (/)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    9 FORMAT (51x,'Composition data')
   11 FORMAT (50x,18('-')/)
  101 FORMAT ('   j     yn         yp        yhe        ynuc        ye      yenu       yanu &
&      yxnu        yl       yl/yl0'/)
  103 FORMAT (1x,i4,11(es11.3))
  105 FORMAT (1x,i4,10(es11.3))
  201 FORMAT ('   j     xn         xp        xhe        xnuc         a          n          z &
&       ye_p       ye_m    ye_p/ye_me  ypair/ye'/)
  203 FORMAT (1x,i4,11(es11.3))
  301 FORMAT ('   j',11(3x,a5,3x))
  303 FORMAT (1x,i4,11(es11.3))                                                  
  305 FORMAT ()
  401 FORMAT (35x,'Composition data - Auxiliary heavy nucleus and burn data')
  403 FORMAT (34x,58('-')/)
  405 FORMAT ('   j   A_rep      Z_rep      be_rep     xn_rep     U_burn   U_burn_encl dUdt nuc&
&   dT nuc       be     1 - xn_tot'/)
  407 FORMAT (1x,i4,11(es11.3))
 1001 FORMAT (' Allocation problem for array ',a10,' in edity')
 2001 FORMAT (' Deallocation problem for array ',a10,' in edity')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                     \\\\\ TEST FOR EDIT /////
!
!  If prnttest == false, or if nedy(i) + 1 >= intedy(i) for any i,
!   proceed through edit. 
!  Otherwise, increment the nedy(i)'s.
!
!-----------------------------------------------------------------------

print_sub           = .false.

IF ( prnttest ) THEN
  DO i = 1,6
    IF ( nedy(i) + 1 >= intedy(i) ) print_sub = .true.
  END DO
ELSE
  print_sub         = .true.
END IF

!-----------------------------------------------------------------------
!  If none of the nedy(i) + 1 >= intedy(i), increment the nedy(i)
!-----------------------------------------------------------------------

IF ( .not. print_sub ) THEN
  DO i = 1,6
    nedy(i)         = nedy(i) + 1
  END DO
END IF ! .not. print_sub

!-----------------------------------------------------------------------
!
!                  \\\\\ PROCEED THROUGH EDIT /////
!
!-----------------------------------------------------------------------

IF ( print_sub ) THEN

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

  ALLOCATE (rnnu(nnu), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rnnu      '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (utburn(nx), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'utburn    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize arrays
!-----------------------------------------------------------------------

  rnnu             = zero
  utburn           = zero

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

  IF ( first ) THEN
    kfm            = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
    n_nucp1        = nuc_number + 1
    first          = .false.
  END IF
  n_time           = nprint

  CALL w_cal( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, wgr, nx )
  wgr(jr_max+1)      = wgr(jr_max)
  DO n = 1,nnu
    CALL flux( jr_min, jr_max, n )
  END DO

ELSE
  RETURN
END IF

!**********************************************************************!
!                                                                      !
!                  Composition data                         part 1     !
!                                                                      !
!----------------------------------------------------------------------!
! yn          'yn'         Neutron number per baryon                   !
! yp          'yp'         Proton number per baryon                    !
! yhe         'yhe'        Helium number per baryon                    !
! ynuc        'ynuc'       Heavy nucleus number per baryon             !
! ye          'ye'         Net electron number per baryon              !
! yenu        'yenu'       e-neutrino number per baryon                !
! yanu        'yanu'       e-antineutrino number per baryon            !
! yxnu        'yxnu'       t-neutrino number per baryon                !
! yl          'yl'         Lepton number per baryon                    !
! ratlep      'yl/yl0'     Ratio of current lepton number to original  !
!                          lepton number                               !
!----------------------------------------------------------------------!
!        Print header.                                                 !
!----------------------------------------------------------------------!        

IF ( prnttest ) THEN
  nedy(1)          = nedy(1) + 1
  IF ( nedy(1) < intedy(1) ) GO TO 1100
  nedy(1)        = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,9)
WRITE (nprint,11)
WRITE (nprint,101)

!-----------------------------------------------------------------------
!  Composition number fractions
!-----------------------------------------------------------------------

DO jv = jr_min,jr_max,idxedy(1)
  j                = jr_max - jv + jr_min

  yn               = aesv(j,7,ij_ray,ik_ray)
  yp               = aesv(j,8,ij_ray,ik_ray)
  ynuc             = aesv(j,9,ij_ray,ik_ray)/( aesv(j,10,ij_ray,ik_ray) + 1.d-2 )
  yhe              = DMAX1( one - yn - yp - aesv(j,9,ij_ray,ik_ray), zero )/4.d0

  DO n = 1,nnu
    rnnu(n)        = zero
    IF ( nnugp(n) /= 0 ) THEN
      rnnu(n)      = SUM( unu(j,:) * unu(j,:) * dunu(j,:) * psi0(j,:,n) ) * ( ncoef/rho(j) )
    END IF !  nnugp(n) ne 0
  END DO

  yenu             = rnnu(1) * rmu
  yanu             = rnnu(2) * rmu
  yxnu             = rnnu(3) * rmu

  yl               = ye(j) + yenu - yanu
  yl0              = ye0(j)
  ratlep           = yl/( yl0 + epsilon )

!-----------------------------------------------------------------------
!  Print composition number fractions
!-----------------------------------------------------------------------

  WRITE (nprint,105) j,yn,yp,yhe,ynuc,ye(j),yenu,yanu,yxnu,yl,ratlep

END DO

!**********************************************************************!
!                                                                      !
!                  Composition data                         part 2     !
!                                                                      !
!----------------------------------------------------------------------!
! xnt         'xn'         Neutron mass fraction                       !
! xp          'xp'         Proton mass fraction                        !
! xhe         'xhe'        Helium mass fraction                        !
! xnuc        'xnuc'       Heavy nucleus mass fraction                 !
! a           'a'          Heavy nucleus mass number                   !
! rn          'n'          Heavy nucleus neutron number                !
! z           'z'          Heavy nucleus proton number                 !
! ye_m        'ye_m'       Total number of electrons/nucleon           !
! ye_p        'ye_p'       Total number of positrons/nucleon           !
! ye_pm       'ye_p/ye_m'  Ratio of positrons to electrons             !
! ypairye     'ypair/ye'   Ratio of pairs to net electrons             !
!----------------------------------------------------------------------!
!        Print header                                                  !
!----------------------------------------------------------------------!

 1100 CONTINUE

IF ( prnttest ) THEN
  nedy(2)          = nedy(2) + 1
  IF ( nedy(2) < intedy(2) ) GO TO 2100
  nedy(2)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,9)
WRITE (nprint,11)
WRITE (nprint,201)

!-----------------------------------------------------------------------
!  Composition mass fractions
!-----------------------------------------------------------------------

DO jv = jr_min,jr_max,idxedy(2)
  j                = jr_max - jv + jr_min

  xnt              = aesv(j,7,ij_ray,ik_ray)
  xp               = aesv(j,8,ij_ray,ik_ray)
  xnuc             = aesv(j,9,ij_ray,ik_ray)
  xhe              = DMAX1( one - xnt - xp - xnuc, zero )
  a                = aesv(j,10,ij_ray,ik_ray)
  z                = aesv(j,11,ij_ray,ik_ray)
  rn               = a - z

!-----------------------------------------------------------------------
!  Call lectron for electron-positron quantities
!-----------------------------------------------------------------------

  dbck             = rho(j) * kfm
  tbck             = t(j) * kmev
  yebck            = ye(j)
  CALL lectron
  ye_p             = yeplus
  ye_m             = yeplus + ye(j)
  ye_pm            = ye_p/ye_m
  ypairye          = yeplus/ye(j)

!-----------------------------------------------------------------------
!  Print composition mass fractions
!-----------------------------------------------------------------------

  WRITE (nprint,203) j,xnt,xp,xhe,xnuc,a,rn,z,ye_p,ye_m,ye_pm,ypairye

END DO

!**********************************************************************!
!                                                                      !
!                  Composition data                         part 3     !
!                                                                      !
!----------------------------------------------------------------------!
! xn(j,i)     'a_name(i)'  (mass fraction of nucleus i)                !
!----------------------------------------------------------------------!
!        Print header.                                                 !
!----------------------------------------------------------------------!
 2100 CONTINUE

IF ( prnttest ) THEN
  nedy(3)          = nedy(3) + 1
  IF ( nedy(3) < intedy(3) ) go to 3100
  nedy(3)          = 0
END IF ! prnttest

i_mlt              = n_nucp1/11
i_rmn              = MOD( n_nucp1, 11 )

IF ( i_mlt /= 0 ) THEN
  DO l = 1,i_mlt

    WRITE (nprint,1)
    WRITE (nprint,3) head
    WRITE (nprint,5)
    CALL date_and_time_print( n_time )
    WRITE (nprint,9)
    WRITE (nprint,11)
    WRITE (nprint,301) (a_name((l-1)*11+n),n=1,11)
    WRITE (nprint,305)

    DO jv = jr_min,jr_max,idxedy(3)
      j            = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Only nse = 0 stuff considered here
!-----------------------------------------------------------------------

      IF ( nse(j,ij_ray,ik_ray) == 0 ) THEN

!-----------------------------------------------------------------------
!  Call eosnuc to evaluate non-nse mass fractions
!-----------------------------------------------------------------------

        dbck       = rho(j) * kfm
        tbck       = t(j) * kmev
        yebck      = ye(j)
        jshel      = j
        CALL eosnuc_x

!-----------------------------------------------------------------------
!  Print non-nse nuclear composition mass fractions
!-----------------------------------------------------------------------

        WRITE (nprint,303) j,(xn(j,(l-1)*11+n),n=1,11)

      END IF ! nse(j,ij_ray,ik_ray) = 0
    END DO ! jv
  END DO ! l
END IF

IF ( i_rmn /= 0 ) THEN

  WRITE (nprint,1)
  WRITE (nprint,3) head
  WRITE (nprint,5)
  CALL date_and_time_print( n_time )
  WRITE (nprint,9)
  WRITE (nprint,11)
  WRITE (nprint,301) (a_name(i_mlt*11+n),n=1,i_rmn)
  WRITE (nprint,305)

  DO jv = jr_min,jr_max,idxedy(3)
    j              = jr_max - jv + jr_min

!-----------------------------------------------------------------------
!  Only nse = 0 stuff considered here
!-----------------------------------------------------------------------

    IF ( nse(j,ij_ray,ik_ray) == 0 ) THEN

!-----------------------------------------------------------------------
!  Call eosnuc to evaluate non-nse mass fractions
!-----------------------------------------------------------------------

      dbck       = rho(j) * kfm
      tbck       = t(j) * kmev
      yebck      = ye(j)
      jshel      = j
      CALL eosnuc_x

!-----------------------------------------------------------------------
!  Print non-nse nuclear composition mass fractions
!-----------------------------------------------------------------------

      WRITE (nprint,303) j,(xn(j,i_mlt*11+n),n=1,i_rmn)

    END IF ! nse(j,ij_ray,ik_ray) = 0
  END DO ! jv
END IF ! i_rmn /= 0

!**********************************************************************!
!                                                                      !
!                  Composition data                         part 4     !
!                                                                      !
!----------------------------------------------------------------------!
! a_nuc_rep(j)  'A_rep"      Mass number of representative heavy       !
!                             nucleus                                  !
! z_nuc_rep(j)  'Z_rep"      Charge number of representative heavy     !
!                             nucleus                                  !
! be_nuc_rep(j) 'be_rep'     Binding energy of representative heavy    !
!                             nucleus                                  !
! xn(j,n_nucp1) 'xn_rep'     Mass fraction of representative heavy     !
!                             nucleus                                  !
! uburn(j)      'uburn'      Cumulative energy released by nuclear     !
!                             burning in zone j (ergs/gm)              !
! utburn(j)     'uburn_encl  Total energy released by nuclear burning  !
!                             interior to zone j (foes)                !
! dtmpmn(j,4,ij_ray,ik_ray)                                            !
!               'dT nuc'     Temperature change due to nuclear         !
!                             reactions                                !
! enb_mev       'be'         Mean nuclear binding energy (MeV)         !
! 1.d0-xn_tot   '1 - xn_tot' 1 - total mass fraction                   !
!----------------------------------------------------------------------!
!        Print header.                                                 !
!----------------------------------------------------------------------!

 3100 CONTINUE
 
IF ( prnttest ) THEN
  nedy(4)          = nedy(4) + 1
  IF ( nedy(4) < intedy(4) ) GO TO 4100
  nedy(4)          = 0
END IF ! prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( n_time )
WRITE (nprint,9)
WRITE (nprint,11)
WRITE (nprint,401)
WRITE (nprint,403)
WRITE (nprint,405)

utburn(1)           = zero
DO j = 2,jr_max
  utburn(j)         = utburn(j-1) + uburn(j,ij_ray,ik_ray) * dmrst(j) * ergfoe
END DO

DO jv = jr_min,jr_max,idxedy(4)
  j                 = jr_max - jv + jr_min
  xn_tot            = 0
  DO n = 1,n_nucp1
    xn_tot          = xn_tot + xn(j,n)
  END DO
  CALL nuc_energy( j, enb_mev, enm )
  WRITE (nprint,407) j, a_nuc_rep(j), z_nuc_rep(j), be_nuc_rep(j), xn(j,n_nucp1), &
& uburn(j,ij_ray,ik_ray), utburn(j), dudt_nuc(j,ij_ray,ik_ray), dtmpmn(j,4,ij_ray,ik_ray), &
& enb_mev, 1.d0-xn_tot
END DO

 4100 CONTINUE

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (rnnu, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnnu      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (utburn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'utburn    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE edity
