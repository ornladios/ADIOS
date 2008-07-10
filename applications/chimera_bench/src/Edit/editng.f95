SUBROUTINE editng( n, jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editng
!    Module:       editng
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/05
!
!    Purpose:
!      To edit integral n-type neutrino transport data
!
!    Subprograms call:
!      date_and_time_print
!
!    Input arguments:
!  n          : neutrino flavor
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
!  nedng(i)   : editng counter for data set i.
!  intdng(i)  : number of cycles between edits of data set i.
!  idxdng(i)  : edit jr_min, jr_max, and every idxdng(i) radial zone between them for data set i.
!
!    Include files:
!      kind_module, array_module
!      edit_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu
USE numerical_module, ONLY : zero

USE edit_module, ONLY : prnttest, nedng, intdng, nprint, head, dudt_ABEM, &
& dudt_Brem, dudt_NAS, dudt_NES, dudt_NNS, dudt_PR, dudt_NET
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : isctn, nes, isctnn, ipair, ibrem, iaefnp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: n             ! neutrino flavor
INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: exit_edit
LOGICAL                          :: l_key

INTEGER                          :: i             ! edit index

REAL(KIND=double), PARAMETER     :: fthird = 4.d0/3.d0

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' n=',i10,' does not have an allowed value')

!-----------------------------------------------------------------------
!
!                  \\\\\ RETURN CRITERIA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Return if
!   prnttest = true
!  and
!   nedng(i,n) < intdng(i,n) for all i
!-----------------------------------------------------------------------

IF ( prnttest ) THEN
  DO i = 1,18
    nedng(i,n)     = nedng(i,n) + 1
  END DO
  exit_edit        = .TRUE.
  DO i = 1,18
    IF ( nedng(i,n) >= intdng(i,n) ) exit_edit = .FALSE.
  END DO
  IF ( exit_edit ) RETURN
END IF ! prnttest

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( n == 1 ) THEN
  dudt_ABEM          = zero
  dudt_Brem          = zero
  dudt_NAS           = zero
  dudt_NES           = zero
  dudt_NNS           = zero
  dudt_PR            = zero
  dudt_NET           = zero
END IF ! n == 1

l_key                = .false.

!-----------------------------------------------------------------------
!
!           \\\\\ CALL NG_EDIT SUBROUTINES AS REQUIRED /////
!
!-----------------------------------------------------------------------

IF ( .not. prnttest ) THEN
  CALL editng_set( n )
  CALL editng_01( n, jr_min, jr_max, ij_ray, ik_ray )
END IF
 
IF ( prnttest ) THEN

  IF ( nedng(1,n) >= intdng(1,n) ) THEN
    CALL editng_set( n )
    CALL editng_01( n, jr_min, jr_max, ij_ray, ik_ray )
    nedng(1,n)      = 0
  END IF ! nedng(1,n) >= intdng(1,n)

END IF ! prnttest

!-----------------------------------------------------------------------
!  Emission and absorption on nucleons
!-----------------------------------------------------------------------

IF ( .not. prnttest  .and.  iaefnp /= 0  .and.  n <= 3 ) THEN
  CALL editng_set( n )
  CALL editng_N_ABEM( n, jr_min, jr_max, ij_ray, ik_ray )
END IF
 
IF ( prnttest  .and.  iaefnp /= 0 ) THEN

  IF ( nedng(5,n) >= intdng(5,n) ) THEN
    CALL editng_set( n )
    CALL editng_N_ABEM( n, jr_min, jr_max, ij_ray, ik_ray )
    nedng(5,n)      = 0
    l_key           = .true.
  END IF ! nedng(5,n) >= intdng(5,n)

END IF ! prnttest .and. iaefnp /= 0

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering
!-----------------------------------------------------------------------

IF ( .not. prnttest  .and.  isctn /= 0 ) THEN
  CALL editng_set( n )
  CALL editng_NNS( n, jr_min, jr_max, ij_ray, ik_ray )
 END IF
 
IF ( prnttest  .and.  isctn /= 0 ) THEN

  IF ( nedng(6,n) >= intdng(6,n) ) THEN
    CALL editng_set( n )
    CALL editng_NNS( n, jr_min, jr_max, ij_ray, ik_ray )
    nedng(6,n)      = 0
    l_key           = .true.
  END IF ! nedng(6,n) >= intdng(6,n)

END IF ! prnttest .and. nes /= 0

!-----------------------------------------------------------------------
!  Neutrino-electron scattering
!-----------------------------------------------------------------------

IF ( .not. prnttest  .and.  nes /= 0 ) THEN
  CALL editng_set( n )
  CALL editng_NES( n, jr_min, jr_max, ij_ray, ik_ray )
END IF
 
IF ( prnttest   .and.  nes /= 0 ) THEN

  IF ( nedng(7,n) >= intdng(7,n) ) THEN
    CALL editng_set( n )
    CALL editng_NES( n, jr_min, jr_max, ij_ray, ik_ray )
    nedng(7,n)      = 0
    l_key           = .true.
  END IF ! nedng(7,n) >= intdng(7,n)

END IF ! prnttest .and. nes /= 0

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering
!-----------------------------------------------------------------------

IF ( .not. prnttest  .and.  isctnn /= 0 ) THEN
  CALL editng_set( n )
  CALL editng_NAS( n, jr_min, jr_max, ij_ray, ik_ray )
END IF

IF ( prnttest  .and.  isctnn /= 0 ) THEN

  IF ( nedng(8,n) >= intdng(8,n) ) THEN
    CALL editng_set( n )
    CALL editng_NAS( n, jr_min, jr_max, ij_ray, ik_ray )
    nedng(8,n)      = 0
    l_key           = .true.
  END IF ! nedng(8,n) >= intdng(8,n)

END IF ! prnttest .and. isctnn /= 0

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation
!-----------------------------------------------------------------------

IF ( .not. prnttest  .and.  ipair /= 0 ) THEN
  CALL editng_set( n )
  CALL editng_PR( n, jr_min, jr_max, ij_ray, ik_ray )
END IF

IF ( prnttest   .and.  ipair /= 0 ) THEN

  IF ( nedng(9,n) >= intdng(9,n) ) THEN
    CALL editng_set( n )
    CALL editng_PR( n, jr_min, jr_max, ij_ray, ik_ray )
    nedng(9,n)      = 0
    l_key           = .true.
  END IF ! nedng(8,n) >= intdng(8,n)

END IF ! prnttest .and. isctnn /= 0

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation
!-----------------------------------------------------------------------

IF ( .not. prnttest  .and.  ibrem /= 0 ) THEN
  CALL editng_set( n )
  CALL editng_Brem( n, jr_min, jr_max, ij_ray, ik_ray )
END IF

IF ( prnttest  .and.  ibrem /= 0 ) THEN

  IF ( nedng(10,n) >= intdng(10,n) ) THEN
    CALL editng_set( n )
    CALL editng_Brem( n, jr_min, jr_max, ij_ray, ik_ray )
    nedng(10,n)     = 0
    l_key           = .true.
  END IF ! nedng(8,n) >= intdng(8,n)

END IF ! prnttest .and. isctnn /= 0

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation
!-----------------------------------------------------------------------

IF ( .not. prnttest  .and.                                             &
&    n == nnu        .and.                                             &
&    ( iaefnp /= 0  .or.                                               &
&      isctn  /= 0  .or.                                               &
&      nes    /= 0  .or.                                               &
&      isctnn /= 0  .or.                                               &
&      ipair  /= 0  .or.                                               &
&      ibrem  /= 0 ) ) THEN
  CALL editng_DUDT( n, jr_min, jr_max )
END IF

IF ( prnttest  .and.  n == nnu  .and.  l_key ) THEN
  CALL editng_DUDT( n, jr_min, jr_max )
END IF ! prnttest

RETURN
END SUBROUTINE editng
