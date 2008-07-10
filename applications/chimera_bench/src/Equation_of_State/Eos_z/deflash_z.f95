SUBROUTINE deflash_z( k, rho, t, ye, ki_ray, kj_ray )
!-----------------------------------------------------------------------
!
!    File:         deflash_z
!    Module:       deflash_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To deflash a mass zone to non-nse. The new temperature
!       of the mass zone is computed such that the internal
!       energy of the zone after being switched from nse is the
!       same as before. 
!
!    Subprograms called:
!  eqstta_z : called to compute binding energy of auxiliary nucleus
!
!    Input arguments:
!  k        : z (azimuthal) zone index of zone to be flashed
!  rho      : azimuthal array of matter density (g/cm**3)
!  t        : azimuthal array of matter temperature (K)
!  ye       : azimuthal array of matter electron fraction
!  ki_ray   : x (radial) index of a specific z (azimuthaal) ray
!  kj_ray   : y (angular) index of a specific z (azimuthaal) ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_bck_module, eos_snc_x_module, eos_snc_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nz, nnc
USE numerical_module, ONLY: zero, one

USE edit_module, ONLY : nlog
USE eos_bck_module, ONLY: b
USE eos_snc_x_module, ONLY: a_name
USE eos_snc_z_module, ONLY: aesv, eos_i, xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, &
& a_nuc, z_nuc, be_nuc, nuc_number

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: k               ! z (azimuthal) zone index of zone to be flashed
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthaal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthaal) ray

REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rho ! azimuthal array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: t   ! azimuthal array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: ye  ! azimuthal array of matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1)                :: eos_is          ! temporary eos identifier

LOGICAL                          :: first = .true.

INTEGER, PARAMETER               :: ncd = 300       ! composition dimension
INTEGER                          :: i               ! composition index
INTEGER                          :: ii              ! composition index
INTEGER                          :: n_nucp1         ! nuc_number + 1
INTEGER                          :: n_hvy           ! number of heavy nuclei (not counting representative heavy nucleus)
INTEGER                          :: i_He            ! helium index
INTEGER                          :: i_neut          ! neutron index
INTEGER                          :: i_prot          ! proton index
INTEGER, DIMENSION(ncd)          :: i_hvy           ! heavy nucleus index

REAL(KIND=double)                :: u1              ! dummy variable
REAL(KIND=double)                :: u2              ! dummy variable
REAL(KIND=double)                :: u3              ! dummy variable
REAL(KIND=double)                :: u4              ! dummy variable

!-----------------------------------------------------------------------
!        Formats.
!-----------------------------------------------------------------------

  101 FORMAT (' nnc =',i4,' > ncd=',i4,' in subroutine deflash_z.')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.

!........dimension consistency

  IF ( nnc > ncd ) THEN
    WRITE (nlog,101) nnc,ncd
    STOP
  END IF

!........indexing

  n_hvy            = 0
  i_He             = 0
  i_neut           = 0
  i_prot           = 0
  i_hvy            = 0
  ii               = 0
  n_nucp1          = nuc_number + 1

  DO i = 1,nuc_number
    IF (      a_name(i) == '  n  ' ) THEN
      i_neut       = i
    ELSE IF ( a_name(i) == '  p  ' ) THEN
      i_prot       = i
    ELSE IF ( a_name(i) == '  4He' ) THEN
      i_He         = i
    ELSE
      ii           = ii + 1
      n_hvy        = n_hvy + 1
      i_hvy(ii)    = i
    END IF
  END DO ! i

END IF ! first

!-----------------------------------------------------------------------
!  Determine non-nse composition.
!-----------------------------------------------------------------------

DO i = 1,n_nucp1
  xn(k,i)          = zero
END DO

xn(k,n_nucp1)      = aesv(k,9,kj_ray,ki_ray)
xn(k,i_He)         = DMAX1( one - aesv(k,7,kj_ray,ki_ray)               &
&                  - aesv(k,8,kj_ray,ki_ray) - aesv(k,9,kj_ray,ki_ray), zero )
xn(k,i_neut)       = aesv(k,7,kj_ray,ki_ray)
xn(k,i_prot)       = aesv(k,8,kj_ray,ki_ray)

a_nuc_rep(k)       = aesv(k,10,kj_ray,ki_ray)
z_nuc_rep(k)       = aesv(k,11,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Determine heavy nucleus binding energy.
!
!  Calling eqstta_z with nse = 1 and eos_i = 'S' causes b, the binding
!   energy per baryon of the representative heavy nucleus to be evaluated
!   in eos via hvbub.
!-----------------------------------------------------------------------

eos_is             = eos_i
eos_i              = 'S'
CALL eqstta_z( 6, k, ki_ray, kj_ray, rho(k), t(k), ye(k), u1, u2, u3, u4 )
eos_i              = eos_is
be_nuc_rep(k)      =  - b * a_nuc_rep(k)

RETURN
END SUBROUTINE deflash_z
