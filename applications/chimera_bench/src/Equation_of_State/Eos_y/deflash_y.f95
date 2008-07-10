SUBROUTINE deflash_y( j, rho, t, ye, ji_ray, jk_ray )
!-----------------------------------------------------------------------
!
!    File:         deflash_y
!    Module:       deflash_y
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
!  eqstta_y : called to compute binding energy of auxiliary nucleus
!
!    Input arguments:
!  j        : y-array zone index of zone to be flashed
!  rho      : angular array of matter density (g/cm**3)
!  t        : angular array of matter temperature (K)
!  ye       : angular array of matter electron fraction
!  ji_ray   : i (radial) index of a specific angular ray
!  jk_ray   : k (azimuthal) index of a specific angular ray
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module, array_module, numerical_module
!  edit_module, eos_bck_module, eos_snc_x_module, eos_snc_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : ny, nnc
USE numerical_module, ONLY: zero, one

USE edit_module, ONLY : nlog
USE eos_bck_module, ONLY: b
USE eos_snc_x_module, ONLY: a_name
USE eos_snc_y_module, ONLY: aesv, eos_i, xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, &
& a_nuc, z_nuc, be_nuc, nuc_number

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j               ! y-array zone index of zone to be flashed
INTEGER, INTENT(in)              :: ji_ray          ! i (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray          ! k (azimuthal) index of a specific angular ray

REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rho ! angular array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: t   ! angular array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: ye  ! angular array of matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1)                :: eos_is         ! temporary eos identifier

LOGICAL                          :: first = .true.

INTEGER, PARAMETER               :: ncd = 300      ! composition dimension
INTEGER                          :: i              ! composition index
INTEGER                          :: ii             ! composition index
INTEGER                          :: n_nucp1        ! nuc_number + 1
INTEGER                          :: n_hvy          ! number of heavy nuclei (not counting representative heavy nucleus)
INTEGER                          :: i_He           ! helium index
INTEGER                          :: i_neut         ! neutron index
INTEGER                          :: i_prot         ! proton index
INTEGER, DIMENSION(ncd)          :: i_hvy          ! heavy nucleus index

REAL(KIND=double)                :: u1             ! dummy variable
REAL(KIND=double)                :: u2             ! dummy variable
REAL(KIND=double)                :: u3             ! dummy variable
REAL(KIND=double)                :: u4             ! dummy variable

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' nnc =',i4,' > ncd=',i4,' in subroutine deflash_y.')

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
  xn(j,i)          = zero
END DO

xn(j,n_nucp1)      = aesv(j,9,ji_ray,jk_ray)
xn(j,i_He)         = DMAX1( one - aesv(j,7,ji_ray,jk_ray)               &
&                  - aesv(j,8,ji_ray,jk_ray) - aesv(j,9,ji_ray,jk_ray), zero )
xn(j,i_neut)       = aesv(j,7,ji_ray,jk_ray)
xn(j,i_prot)       = aesv(j,8,ji_ray,jk_ray)

a_nuc_rep(j)       = aesv(j,10,ji_ray,jk_ray)
z_nuc_rep(j)       = aesv(j,11,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Determine heavy nucleus binding energy.
!
!  Calling eqstta_y with nse = 1 and eos_i = 'S' causes b, the binding
!   energy per baryon of the representative heavy nucleus to be evaluated
!   in eos via hvbub.
!-----------------------------------------------------------------------

eos_is             = eos_i
eos_i              = 'S'
CALL eqstta_y( 6, j, ji_ray, jk_ray, rho(j), t(j), ye(j), u1, u2, u3, u4 )
eos_i              = eos_is
be_nuc_rep(j)      =  - b * a_nuc_rep(j)

RETURN
END SUBROUTINE deflash_y
