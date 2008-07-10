SUBROUTINE deflash_x( j, rho, t, ye, ij_ray, ik_ray, eos_reset )
!-----------------------------------------------------------------------
!
!    File:         deflash_x
!    Module:       deflash_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/07/03
!
!    Purpose:
!      To deflash a mass zone to non-nse. The new temperature
!       of the mass zone is computed such that the internal
!       energy of the zone after being switched from nse is the
!       same as before. 
!
!    Subprograms called:
!  esrgn_x  : reganerates local EOS table if necessary
!  eqstt_x  : evaluates EOS quantities by interpolation
!  eqstta_x : evaluates EOS quantities directly rather than by interpolation
!
!    Input arguments:
!  j        : radial zone index of zone to be flashed
!  rho      : shifted matter density array (g/cm**3).
!  t        : shifted matter matter temperature array (K).
!  ye       : shifted matter matter electron fraction array.
!  ij_ray   : j-index of a radial ray
!  ik_ray   : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_bck_module, eos_snc_x_module, nucbrn_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nnc
USE numerical_module, ONLY: zero, one, epsilon

USE edit_module, ONLY : nlog
USE eos_bck_module, ONLY: b
USE eos_snc_x_module, ONLY: aesv, eos_i, idr, itr, iyr, duesrc, nse, xn_e=>xn, &
& a_nuc_rep_e=>a_nuc_rep, z_nuc_rep_e=>z_nuc_rep, be_nuc_rep_e=>be_nuc_rep, &
& nuc_number, be_nuc
USE nucbrn_module, ONLY: xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, a_name

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER (len=1), INTENT(in)    :: eos_reset       ! EOS reset flag

INTEGER, INTENT(in)              :: j               ! radial zone index
INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t   ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ye  ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1)                :: eos_is          ! temporary eos identifier

LOGICAL                          :: first = .true.

INTEGER, PARAMETER               :: ncd = 300       ! composition dimension
INTEGER                          :: nse_save        ! initial value of nse(j,ij_ray,ik_ray)
INTEGER                          :: i               ! composition index
INTEGER                          :: ii              ! composition index
INTEGER                          :: n_nucp1         ! nuc_number + 1
INTEGER                          :: n_hvy           ! number of heavy nuclei (not counting representative heavy nucleus)
INTEGER                          :: i_He            ! helium index
INTEGER                          :: i_neut          ! neutron index
INTEGER                          :: i_prot          ! proton index
INTEGER                          :: i_Ni            ! 56Ni index
INTEGER, DIMENSION(ncd)          :: i_hvy           ! heavy nucleus index

REAL(KIND=double)                :: e_NSE           ! internal energy for matter in NSE
REAL(KIND=double)                :: e_nonNSE        ! internal energy for matter not in NSE
REAL(KIND=double)                :: dedd            ! d(energy)/d(density)
REAL(KIND=double)                :: dedt            ! d(energy)/d(temperature)
REAL(KIND=double)                :: dedy            ! d(energy)/d(electron fraction)
REAL(KIND=double), PARAMETER     :: x_Fe = 0.25d0   ! maxximum allowed mass fraction of 56Fe on deflashing
!REAL(KIND=double), PARAMETER     :: x_Fe = 0.85d0  ! maxximum allowed mass fraction of 56Fe on deflashing
REAL(KIND=double)                :: xn_A            ! mass fraction heavy nucleus before deflashing
REAL(KIND=double)                :: ZA_min          ! minimum ratio of Z to A for deflashing
REAL(KIND=double)                :: ZA_Fe           ! charge to mass ratio of 56Fe
REAL(KIND=double)                :: A_rep           ! mass number of representative heavy nucleus
REAL(KIND=double)                :: ZA_rep          ! charge to mass ratio of representative heavy nucleus
REAL(KIND=double)                :: be_rep          ! binding energy of representative heavy nucleus

!-----------------------------------------------------------------------
!        Formats.
!-----------------------------------------------------------------------

  101 FORMAT (' nnc =',i4,' > ncd=',i4,' in subroutine deflash_x.')
  201 FORMAT (' i_Ni cannot be foound in subroutine deflash_x.')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nse_save                  = nse(j,ij_ray,ik_ray)

IF ( nnc > ncd ) THEN
  WRITE (nlog,101) nnc,ncd
  STOP
END IF

IF ( first ) THEN
  first                   = .false.

  ZA_Fe                   = 26.d0/56.d0
  ZA_min                  = ZA_Fe * x_Fe + 0.5d0 * ( 1.d0 - x_Fe )

!........indexing

  n_hvy                   = 0
  i_He                    = 0
  i_neut                  = 0
  i_prot                  = 0
  i_hvy                   = 0
  i_Ni                    = 0
  ii                      = 0
  n_nucp1                 = nuc_number + 1

  DO i = 1,nuc_number
    IF (      a_name(i) == '  n  ' ) THEN
      i_neut              = i
    ELSE IF ( a_name(i) == '  p  ' ) THEN
      i_prot              = i
    ELSE IF ( a_name(i) == '  4He' ) THEN
      i_He                = i
    ELSE
      ii                  = ii + 1
      n_hvy               = n_hvy + 1
      i_hvy(ii)           = i
    END IF
  END DO ! i
  
  DO i = 1,nuc_number
    IF ( a_name(i) == ' 56Ni' ) THEN
      i_Ni                = i
      EXIT
    END IF ! a_name(i) == ' 56Ni'
    IF ( i == nuc_number ) THEN
      WRITE (nlog,201)
      STOP
    END IF ! i == nuc_number
  END DO ! i

END IF ! first

!-----------------------------------------------------------------------
!
!             \\\\\ DETERMINE NON-NSE COMPOSITION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Update nse composition to ensure that abundances are consistent with
!   current ye
!-----------------------------------------------------------------------

CALL eqstz_x( j, j, rho, t, ye, ij_ray, ik_ray )

!........initialize

xn(j,:)                   = zero
xn_e(j,:)                 = zero

!........helium

xn(j,i_He)                = DMAX1( one - aesv(j,7,ij_ray,ik_ray) - aesv(j,8,ij_ray,ik_ray) &
&                         - aesv(j,9,ij_ray,ik_ray), zero )
xn_e(j,i_He)              = xn(j,i_He)

!........neutrons

xn(j,i_neut)              = aesv(j,7,ij_ray,ik_ray)
xn_e(j,i_neut)            = xn(j,i_neut)

!........protons

xn(j,i_prot)              = aesv(j,8,ij_ray,ik_ray)
xn_e(j,i_prot)            = xn(j,i_prot)

!-----------------------------------------------------------------------
!  Split heavy nucleus into 56Ni and "56Fe-like" conserving mass and
!   charge if a 56Fe mass fraction less than x_Fe is added.
!
!     xn(j,i_Ni) + xn_Fe = xn_A
!     0.5*xn(j,i_Ni) + ZA_Fe*xn_Fe = ZA_rep*xn_A
!
! ==> 0.5*xn(j,i_Ni) + ZA_Fe*( ZA_rep - n(j,i_Ni) = ZA_rep*xn_A
!-----------------------------------------------------------------------

A_rep                     = aesv(j,10,ij_ray,ik_ray)
ZA_rep                    = aesv(j,11,ij_ray,ik_ray)/( A_rep + epsilon )
IF ( ZA_rep > ZA_min  .and.  ZA_rep <= 0.5d0 ) THEN

  xn_A                    = aesv(j,9,ij_ray,ik_ray)
  xn(j,i_Ni)              = xn_A * ( ZA_rep - ZA_Fe )/( 0.5d0 - ZA_Fe )
  xn_e(j,i_Ni)            = xn(j,i_Ni)
  xn(j,n_nucp1)           = xn_A - xn(j,i_Ni)
  xn_e(j,n_nucp1)         = xn(j,n_nucp1)

!-----------------------------------------------------------------------
!  Determine heavy nucleus binding energy if
!   ZA_rep > ZA_min  .and.  ZA_rep <= 0.5d0
!
!  Calling eqstta_x with nse = 1 and eos_i = 'S' causes b, the binding
!   energy per baryon of the representative heavy nucleus to be evaluated
!   in eos via hvbub.
!
!     be_nuc_rep(j)*xn(j,n_nucp1)/A_Fe + be_nuc(i_Ni)*xn(j,i_Ni)/A_ni
!           = be_rep*xn_A/xn_A
!-----------------------------------------------------------------------

  eos_is                  = eos_i
  eos_i                   = 'S'
  CALL eqstta_x( 6, j, ij_ray, ik_ray, rho(j), t(j), ye(j), e_nonNSE, &
&  dedd, dedt, dedy )
  eos_i                   = eos_is
  be_rep                  =  - b * aesv(j,10,ij_ray,ik_ray)
  be_nuc_rep(j)           = ( xn_A * be_rep/( A_rep + epsilon )       &
&                         - xn(j,i_Ni) * be_nuc(i_Ni)/56.d0 )         &
&                         / ( xn(j,n_nucp1)/56.d0 + epsilon )
  be_nuc_rep_e(j)         = be_nuc_rep(j)
  a_nuc_rep(j)            = 56.d0
  a_nuc_rep_e(j)          = 56.d0
  z_nuc_rep(j)            = 26.d0
  z_nuc_rep_e(j)          = 26.d0

  IF ( a_nuc_rep(j) < 1.d0 ) THEN
    a_nuc_rep(j)          = 56.d0
    a_nuc_rep_e(j)        = 56.d0
    z_nuc_rep(j)          = 28.d0
    z_nuc_rep_e(j)        = 28.d0
    be_nuc_rep(j)         = 492.3d0
    be_nuc_rep_e(j)       = 492.3d0
  END IF

ELSE

  xn(j,n_nucp1)           = aesv(j,9,ij_ray,ik_ray)
  xn_e(j,n_nucp1)         = xn(j,n_nucp1)

  a_nuc_rep(j)            = aesv(j,10,ij_ray,ik_ray)
  z_nuc_rep(j)            = aesv(j,11,ij_ray,ik_ray)
  a_nuc_rep_e(j)          = a_nuc_rep(j)
  z_nuc_rep_e(j)          = z_nuc_rep(j)

!-----------------------------------------------------------------------
!  Determine heavy nucleus binding energy.
!
!  Calling eqstta_x with nse = 1 and eos_i = 'S' causes b, the binding
!   energy per baryon of the representative heavy nucleus to be evaluated
!   in eos via hvbub.
!-----------------------------------------------------------------------

  eos_is                  = eos_i
  eos_i                   = 'S'
  CALL eqstta_x( 6, j, ij_ray, ik_ray, rho(j), t(j), ye(j), e_nonNSE, &
&  dedd, dedt, dedy )
  eos_i                   = eos_is
  be_nuc_rep(j)           =  - b * a_nuc_rep(j)
  be_nuc_rep_e(j)         =  - b * a_nuc_rep(j)

END IF

!-----------------------------------------------------------------------
!  Reset EOS table with non-nse values if eos_reset = 'y'.
!-----------------------------------------------------------------------

IF ( eos_reset == 'y' ) THEN

  nse(j,ij_ray,ik_ray)    = 1
  CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), e_NSE, dedd, &
&  dedt, dedy )

  nse(j,ij_ray,ik_ray)    = 0
  idr(j,ij_ray,ik_ray)    = 0
  itr(j,ij_ray,ik_ray)    = 0
  iyr(j,ij_ray,ik_ray)    = 0

  CALL esrgn_x( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), e_nonNSE, dedd, &
&  dedt, dedy )
  duesrc(j,ij_ray,ik_ray) = duesrc(j,ij_ray,ik_ray) + e_nonNSE - e_NSE

!-----------------------------------------------------------------------
!  Restore original value of nse
!  (Change this value in nse_test if a zone is deflashed)
!-----------------------------------------------------------------------

  nse(j,ij_ray,ik_ray)    = nse_save

END IF ! eos_reset == 'y'

RETURN
END SUBROUTINE deflash_x
