SUBROUTINE enuvplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, &
& kmin, kmax, ny, nz, nnu )
!-----------------------------------------------------------------------
!
!    File:         enuvplot
!    Module:       enuvplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/11/00
!
!    Purpose:
!      To create files for neutrino mean energy vs time plots.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  jmin         : minimum y-array index for the edit
!  jmax         : maximum y-array index for the edit
!  kmin         : minimum z-array index for the edit
!  kmax         : maximum z-array index for the edit
!  ny           : y_array extent
!  nz           : z_array extent
!  nnu          : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  cycle_module, edit_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero, epsilon, frpi

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, ienuplt, ncyenu, intenu, nenu, dtimeplot, &
& r_e_rms, d_e_rms, nenuplt1, nenuplt2, d_omega, nlog, data_path
USE mdl_cnfg_module, ONLY : jr_min, jr_max, r, rho
USE nu_dist_module, ONLY : unu, dunu, psi0
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE t_cntrl_module, ONLY : time, dtnph

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)               :: ik_ray        ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)               :: jmin          ! minimum y-array index for the edit
INTEGER, INTENT(in)               :: jmax          ! maximum y-array index for the edit
INTEGER, INTENT(in)               :: kmin          ! minimum z-array index for the edit
INTEGER, INTENT(in)               :: kmax          ! maximum z-array index for the edit
INTEGER, INTENT(in)               :: ny            ! y-array extent
INTEGER, INTENT(in)               :: nz            ! z-array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                           :: lprint
LOGICAL                           :: l_r_e_rms
LOGICAL                           :: l_d_e_rms

INTEGER                           :: j             ! radial zone index
INTEGER                           :: jr            ! particular radial zone index
INTEGER                           :: jd            ! particular radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index
INTEGER                           :: itime         ! used to determine when to WRITE a file
INTEGER                           :: itimeprev     ! used to determine when to WRITE a file
INTEGER                           :: istat         ! open and close file flag
INTEGER                           :: nenu_test     ! test for editing data
INTEGER, PARAMETER                :: n_zero = 0

REAL(KIND=double)                 :: tmult         ! used to determine when to WRITE a file

REAL(KIND=double)                 :: w2dw          ! unu(jd)**2 * dunu(jd)
REAL(KIND=double)                 :: w3dw          ! unu(jd)**3 * dunu(jd)
REAL(KIND=double)                 :: w4dw          ! unu(jd)**4 * dunu(jd)
REAL(KIND=double)                 :: w5dw          ! unu(jd)**5 * dunu(jd)
REAL(KIND=double)                 :: w2dwm         ! unu(jd-1)**2 * dunu(jd-1)
REAL(KIND=double)                 :: w3dwm         ! unu(jd-1)**3 * dunu(jd-1)
REAL(KIND=double)                 :: w4dwm         ! unu(jd-1)**4 * dunu(jd-1)
REAL(KIND=double)                 :: w5dwm         ! unu(jd-1)**5 * dunu(jd-1)

REAL(KIND=double)                 :: psi0jp2       ! psi0(jd,k,n) * unu(jd)**2 * dunu(jd)
REAL(KIND=double)                 :: psi0jp3       ! psi0(jd,k,n) * unu(jd)**3 * dunu(jd)
REAL(KIND=double)                 :: psi0jp4       ! psi0(jd,k,n) * unu(jd)**4 * dunu(jd)
REAL(KIND=double)                 :: psi0jp5       ! psi0(jd,k,n) * unu(jd)**5 * dunu(jd)
REAL(KIND=double)                 :: psi0jm2       ! psi0(jd-1,k,n) * unu(jd-1)**2 * dunu(jd-1)
REAL(KIND=double)                 :: psi0jm3       ! psi0(jd-1,k,n) * unu(jd-1)**3 * dunu(jd-1)
REAL(KIND=double)                 :: psi0jm4       ! psi0(jd-1,k,n) * unu(jd-1)**4 * dunu(jd-1)
REAL(KIND=double)                 :: psi0jm5       ! psi0(jd-1,k,n) * unu(jd-1)**5 * dunu(jd-1)

REAL(KIND=double)                 :: enuvjp1       ! psi0jp3/psi0jp2
REAL(KIND=double)                 :: enuvjm1       ! psi0jm3/psi0jp2
REAL(KIND=double)                 :: enuvjp2       ! psi0jp4/psi0jp2
REAL(KIND=double)                 :: enuvjm2       ! psi0jm4/psi0jm2
REAL(KIND=double)                 :: enuvjp3       ! psi0jp5/psi0jp3
REAL(KIND=double)                 :: enuvjm3       ! psi0jm5/psi0jm3

REAL(KIND=double)                 :: enuv22        ! psi0j4/psi0j2 at jd
REAL(KIND=double)                 :: enuv32        ! psi0j5/psi0j3 at jd

REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: enuv1r   ! psi0j3/psi0j2 at jr
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: enuv2r   ! square root of psi0j4/psi0j2 at jr
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: enuv3r   ! square root of psi0j5/psi0j3 at jr
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: enuv1rho ! psi0j3/psi0j2 at jd
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: enuv2rho ! square root of psi0j4/psi0j2 at jd
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: enuv3rho ! square root of psi0j5/psi0j3 at jd

REAL(KIND=double), DIMENSION(nnu,ny,nz) :: enuv1r_a    ! psi0j3/psi0j2 at jr
REAL(KIND=double), DIMENSION(nnu,ny,nz) :: enuv2r_a    ! square root of psi0j4/psi0j2 at jr
REAL(KIND=double), DIMENSION(nnu,ny,nz) :: enuv3r_a    ! square root of psi0j5/psi0j3 at jr
REAL(KIND=double), DIMENSION(nnu,ny,nz) :: enuv1rho_a  ! psi0j3/psi0j2 at jd
REAL(KIND=double), DIMENSION(nnu,ny,nz) :: enuv2rho_a  ! square root of psi0j4/psi0j2 at jd
REAL(KIND=double), DIMENSION(nnu,ny,nz) :: enuv3rho_a  ! square root of psi0j5/psi0j3 at jd

REAL(KIND=double), DIMENSION(nnu) :: enuv1r_bar    ! angular mean of psi0j3/psi0j2 at jr
REAL(KIND=double), DIMENSION(nnu) :: enuv2r_bar    ! angular mean of square root of psi0j4/psi0j2 at jr
REAL(KIND=double), DIMENSION(nnu) :: enuv3r_bar    ! angular mean of square root of psi0j5/psi0j3 at jr
REAL(KIND=double), DIMENSION(nnu) :: enuv1rho_bar  ! angular mean of psi0j3/psi0j2 at jd
REAL(KIND=double), DIMENSION(nnu) :: enuv2rho_bar  ! angular mean of square root of psi0j4/psi0j2 at jd
REAL(KIND=double), DIMENSION(nnu) :: enuv3rho_bar  ! angular mean of square root of psi0j5/psi0j3 at jd

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Need MPI version of enuvplot since n_proc=',i4,' > 1')
  201 FORMAT ('     time         e_mn (3:2)   eb_mn (3:2)   x_mn (3:2)  xb_mn (3:2) &
& e_rms (4:2) eb_rms (4:2)  x_rms (4:2) xb_rms (4:2)  e_rms (5:3) eb_rms (5:3)  x_rms (5:3) xb_rms ((5:3)')
  203 FORMAT (es15.8,12es13.4)
 1001 FORMAT (' jd cannot be found in subroutine enuvplot for r')
 4001 FORMAT (' jd cannot be found in subroutine enuvplotfor rho')
 8001 FORMAT (' File nenuplt1 cannot be opened in subroutime enuvplot')
 8501 FORMAT (' File nenuplt2 cannot be opened in subroutime enuvplot')

!-----------------------------------------------------------------------
!
!                \\\\\ CRITERIA FOR PRINTING DATA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Stop if n_proc > 1.
!-----------------------------------------------------------------------

IF ( n_proc > 1 ) THEN
  WRITE (nprint,101) n_proc
  WRITE (nlog,101) n_proc
  STOP
END IF ! n_proc > 1

!-----------------------------------------------------------------------
!  Return if ienuplt = 0.
!-----------------------------------------------------------------------

IF ( ienuplt == 0 ) RETURN

!-----------------------------------------------------------------------
!  Return if no neutrinos.
!-----------------------------------------------------------------------

IF (nnugpmx == 0) RETURN

!-----------------------------------------------------------------------
!  Evaluate criterion for writing entries to luminosity plot files.
!
!  Cycle criterion.
!-----------------------------------------------------------------------

lprint             = .false.

IF ( ncycle > ncyenu(2) ) THEN
  nenu_test        = intenu(3)
ELSE IF ( ncycle > ncyenu(1) ) THEN
  nenu_test        = intenu(2)
ELSE
  nenu_test        = intenu(1)
END IF ! ncycle > ncyenu(2)

nenu               = nenu + 1

IF ( nenu >= nenu_test ) THEN
  lprint           = .true.
  IF ( ij_ray == ij_ray_dim  .and.  ij_ray == ij_ray_dim ) nenu = n_zero
END IF

!-----------------------------------------------------------------------
!  Time criterion.
!-----------------------------------------------------------------------

tmult              = 1.d+3/dtimeplot
itime              = INT( time * tmult )
itimeprev          = INT( ( time - dtnph ) * tmult )
IF ( itime > itimeprev ) lprint = .true.

!-----------------------------------------------------------------------
!  Return if lprint = false.
!-----------------------------------------------------------------------

IF ( .not. lprint ) RETURN

!-----------------------------------------------------------------------
!
!         \\\\\ MEAN NEUTRINO ENERGIES AT RADIUS R_E_RMS /////
!
!-----------------------------------------------------------------------

enuv1r                = zero
enuv2r                = zero
enuv3r                = zero
enuv1rho              = zero
enuv2rho              = zero
enuv3rho              = zero

!-----------------------------------------------------------------------
!  Find jr such that r(jr) > r_e_rms and r(jr-1) < r_e_rms.
!-----------------------------------------------------------------------

l_r_e_rms             = .false.
DO j = jr_min,jr_max
  IF ( r(j) >= r_e_rms ) THEN
    jr             = j
    l_r_e_rms         = .true.
    EXIT
  END IF ! r(j) > r_e_rms
END DO

IF ( .not. l_r_e_rms ) WRITE (nprint,1001)

!-----------------------------------------------------------------------
!  Compute mean energies if jr found.
!-----------------------------------------------------------------------

IF ( l_r_e_rms ) THEN

  DO n = 1,nnu

    psi0jp2         = zero
    psi0jp3         = zero
    psi0jp4         = zero
    psi0jp5         = zero
    psi0jm2         = zero
    psi0jm3         = zero
    psi0jm4         = zero
    psi0jm5         = zero

    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)

      w2dw          = unu(jr  ,k)**2 * dunu(jr  ,k)
      w3dw          = unu(jr  ,k)**3 * dunu(jr  ,k)
      w4dw          = unu(jr  ,k)**4 * dunu(jr  ,k)
      w5dw          = unu(jr  ,k)**5 * dunu(jr  ,k)
      w2dwm         = unu(jr-1,k)**2 * dunu(jr-1,k)
      w3dwm         = unu(jr-1,k)**3 * dunu(jr-1,k)
      w4dwm         = unu(jr-1,k)**4 * dunu(jr-1,k)
      w5dwm         = unu(jr-1,k)**5 * dunu(jr-1,k)

      psi0jp2       = psi0jp2 + w2dw  * psi0(jr  ,k,n)
      psi0jp3       = psi0jp3 + w3dw  * psi0(jr  ,k,n)
      psi0jp4       = psi0jp4 + w4dw  * psi0(jr  ,k,n)
      psi0jp5       = psi0jp5 + w5dw  * psi0(jr  ,k,n)
      psi0jm2       = psi0jm2 + w2dwm * psi0(jr-1,k,n)
      psi0jm3       = psi0jm3 + w3dwm * psi0(jr-1,k,n)
      psi0jm4       = psi0jm4 + w4dwm * psi0(jr-1,k,n)
      psi0jm5       = psi0jm5 + w5dwm * psi0(jr-1,k,n)

    END DO ! k

    enuvjp1         = psi0jp3 /( psi0jp2 + epsilon )
    enuvjm1         = psi0jm3 /( psi0jm2 + epsilon )
    enuvjp2         = psi0jp4 /( psi0jp2 + epsilon )
    enuvjm2         = psi0jm4 /( psi0jm2 + epsilon )
    enuvjp3         = psi0jp5 /( psi0jp3 + epsilon )
    enuvjm3         = psi0jm5 /( psi0jm3 + epsilon )

    enuv1r(n,ij_ray,ik_ray) = rinterp( enuvjp1, enuvjm1, r(jr), r_e_rms, r(jr-1) )
    enuv22          = rinterp( enuvjp2, enuvjm2, r(jr), r_e_rms, r(jr-1) )
    enuv2r(n,ij_ray,ik_ray) = DSQRT( DABS(enuv22) + epsilon )
    enuv32          = rinterp( enuvjp3, enuvjm3, r(jr), r_e_rms, r(jr-1) )
    enuv3r(n,ij_ray,ik_ray) = DSQRT( DABS(enuv32) + epsilon )

  END DO ! n

!-----------------------------------------------------------------------
!  End computation of mean energies at r_e_rms.
!-----------------------------------------------------------------------

END IF ! l_r_e_rms

!-----------------------------------------------------------------------
!
!         \\\\\ MEAN NEUTRINO ENERGIES AT DENSITY D_E_RMS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find jd such that rho(jd) < d_e_rms and rho(jd-1) > d_e_rms.
!-----------------------------------------------------------------------

l_d_e_rms           = .false.
DO j = jr_min,jr_max
  IF ( rho(j) <= d_e_rms ) THEN
    jd             = j
    l_d_e_rms       = .true.
    EXIT
  END IF ! rrh(j) <= r_e_rms
END DO

IF ( .not. l_d_e_rms ) WRITE (nprint,4001)

!-----------------------------------------------------------------------
!  Compute mean energies if jd found.
!-----------------------------------------------------------------------

IF ( l_d_e_rms ) THEN

  DO n = 1,nnu

    psi0jp2         = zero
    psi0jp3         = zero
    psi0jp4         = zero
    psi0jp5         = zero
    psi0jm2         = zero
    psi0jm3         = zero
    psi0jm4         = zero
    psi0jm5         = zero

    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)

      w2dw          = unu(jd  ,k)**2 * dunu(jd  ,k)
      w3dw          = unu(jd  ,k)**3 * dunu(jd  ,k)
      w4dw          = unu(jd  ,k)**4 * dunu(jd  ,k)
      w5dw          = unu(jd  ,k)**5 * dunu(jd  ,k)
      w2dwm         = unu(jd-1,k)**2 * dunu(jd-1,k)
      w3dwm         = unu(jd-1,k)**3 * dunu(jd-1,k)
      w4dwm         = unu(jd-1,k)**4 * dunu(jd-1,k)
      w5dwm         = unu(jd-1,k)**5 * dunu(jd-1,k)

      psi0jp2       = psi0jp2 + w2dw  * psi0(jd  ,k,n)
      psi0jp3       = psi0jp3 + w3dw  * psi0(jd  ,k,n)
      psi0jp4       = psi0jp4 + w4dw  * psi0(jd  ,k,n)
      psi0jp5       = psi0jp5 + w5dw  * psi0(jd  ,k,n)
      psi0jm2       = psi0jm2 + w2dwm * psi0(jd-1,k,n)
      psi0jm3       = psi0jm3 + w3dwm * psi0(jd-1,k,n)
      psi0jm4       = psi0jm4 + w4dwm * psi0(jd-1,k,n)
      psi0jm5       = psi0jm5 + w5dwm * psi0(jd-1,k,n)

    END DO ! k = 1,nnugp(n)

    enuvjp1         = psi0jp3 /( psi0jp2 + epsilon )
    enuvjm1         = psi0jm3 /( psi0jm2 + epsilon )
    enuvjp2         = psi0jp4 /( psi0jp2 + epsilon )
    enuvjm2         = psi0jm4 /( psi0jm2 + epsilon )
    enuvjp3         = psi0jp5 /( psi0jp3 + epsilon )
    enuvjm3         = psi0jm5 /( psi0jm3 + epsilon )

    enuv1rho(n,ij_ray,ik_ray) = rinterp( enuvjp1, enuvjm1, rho(jd), d_e_rms, rho(jd-1) )
    enuv22          = rinterp( enuvjp2, enuvjm2, rho(jd), d_e_rms, rho(jd-1) )
    enuv2rho(n,ij_ray,ik_ray) = DSQRT( DABS(enuv22) + epsilon )
    enuv32          = rinterp( enuvjp3, enuvjm3, rho(jd), d_e_rms, rho(jd-1) )
    enuv3rho(n,ij_ray,ik_ray) = DSQRT( DABS(enuv32) + epsilon )

  END DO ! n = 1,nnu

!-----------------------------------------------------------------------
!  End computation of mean energies at d_e_rms.
!-----------------------------------------------------------------------

END IF ! l_d_e_rms

!-----------------------------------------------------------------------
!
!     \\\\\ ANGULARLY AVERAGE NEUTRINO LUMINOSITIES AND EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if ij_ray /= ij_ray_dim or ij_ray /= ij_ray_dim
!-----------------------------------------------------------------------

IF ( ij_ray /= ij_ray_dim  .or.  ij_ray /= ij_ray_dim ) RETURN

!-----------------------------------------------------------------------
!  Angularly average enuv1r, enuv2r, and enuv3r.
!-----------------------------------------------------------------------

IF ( l_r_e_rms ) THEN

!-----------------------------------------------------------------------
!  Load angular arrays
!-----------------------------------------------------------------------

  enuv1r_a(:,:,:)   = enuv1r(:,:,:)
  enuv2r_a(:,:,:)   = enuv2r(:,:,:)
  enuv3r_a(:,:,:)   = enuv3r(:,:,:)

!-----------------------------------------------------------------------
!  Calculate angular means
!-----------------------------------------------------------------------

  DO n = 1, nnu
    enuv1r_bar(n)   = SUM( enuv1r_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    enuv2r_bar(n)   = SUM( enuv2r_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    enuv3r_bar(n)   = SUM( enuv3r_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! n = 1, nnu

!-----------------------------------------------------------------------
!  Open enuv_r file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nenuplt1,FILE=TRIM(data_path)//'/Plot_Files/enuv_r.d', &
&  STATUS='new',IOSTAT=istat)
  IF ( istat == 0 ) WRITE (nenuplt1,201)
  IF ( istat /= 0 ) OPEN (UNIT=nenuplt1,FILE=TRIM(data_path)//'/Plot_Files/enuv_r.d', &
&  STATUS='old', POSITION='append',IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8001)
    STOP
  END IF ! istat /= 0

!-----------------------------------------------------------------------
!  Write to enuv_r.d file.
!-----------------------------------------------------------------------

  WRITE (nenuplt1,203) time, (enuv1r_bar(n),n=1,nnu), (enuv2r_bar(n),n=1,nnu), &
&  (enuv3r_bar(n),n=1,nnu)

!-----------------------------------------------------------------------
!  CLose enuv.r file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nenuplt1, STATUS='keep')

!-----------------------------------------------------------------------
!  End editing of mean energies at r_e_rms.
!-----------------------------------------------------------------------

END IF ! l_r_e_rms

!-----------------------------------------------------------------------
!  Angularly average enuv1rho, enuv2rho, and enuv3rho.
!-----------------------------------------------------------------------

IF ( l_d_e_rms ) THEN

!-----------------------------------------------------------------------
!  Load angular arrays
!-----------------------------------------------------------------------

  enuv1rho_a(:,:,:) = enuv1rho(:,:,:)
  enuv2rho_a(:,:,:) = enuv2rho(:,:,:)
  enuv3rho_a(:,:,:) = enuv3rho(:,:,:)

!-----------------------------------------------------------------------
!  Calculate angular means
!-----------------------------------------------------------------------

  DO n = 1, nnu
    enuv1rho_bar(n) = SUM( enuv1rho_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    enuv2rho_bar(n) = SUM( enuv2rho_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    enuv3rho_bar(n) = SUM( enuv3rho_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! n = 1, nnu

!-----------------------------------------------------------------------
!  Open enuv_rho file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nenuplt2,FILE=TRIM(data_path)//'/Plot_Files/enuv_rho.d', &
&  STATUS='new',IOSTAT=istat)
  IF ( istat == 0 ) WRITE (nenuplt2,201)
  IF ( istat /= 0 ) OPEN (UNIT=nenuplt2,FILE=TRIM(data_path)//'/Plot_Files/enuv_rho.d', &
&  STATUS='old', POSITION='append',IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8501)
    STOP
  END IF

!-----------------------------------------------------------------------
!  Write to enuv.rho file.
!-----------------------------------------------------------------------

  WRITE (nenuplt2,203) time,(enuv1rho_bar(n),n=1,nnu),(enuv2rho_bar(n),n=1,nnu), &
&  (enuv3rho_bar(n),n=1,nnu)

!-----------------------------------------------------------------------
!  CLose enuv.rho file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nenuplt2, STATUS='keep')

!-----------------------------------------------------------------------
!  End computation of mean energies at d_e_rms.
!-----------------------------------------------------------------------

END IF ! l_d_e_rms

!-----------------------------------------------------------------------
!  Return.
!-----------------------------------------------------------------------

RETURN

CONTAINS

  REAL (KIND=double) FUNCTION rinterp(a,b,x,y,z)

  REAL (KIND=double) :: a
  REAL (KIND=double) :: b
  REAL (KIND=double) :: x
  REAL (KIND=double) :: y
  REAL (KIND=double) :: z

  rinterp          = b + ( a - b ) * ( y - z )/( x - z )

  END FUNCTION rinterp

END SUBROUTINE enuvplot
