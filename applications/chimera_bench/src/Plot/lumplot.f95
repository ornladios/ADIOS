SUBROUTINE lumplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, &
& kmin, kmax, ny, nz, nnu )
!-----------------------------------------------------------------------
!
!    File:         lumplot
!    Module:       lumplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To create files for angularly summed luminosity vs time plots.
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
!  kind_module, array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero, frpi
USE physcnst_module, ONLY : ergfoe

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, ilumplt, ncylum, intlum, nlum, dtimeplot, &
& r_lum, nlumplt1, d_lum, nlumplt2, d_omega, nlog, data_path
USE mdl_cnfg_module, ONLY : jr_min, jr_max, r, rho
USE nu_dist_module, ONLY : fluxnu
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
LOGICAL                           :: l_r_lum
LOGICAL                           :: l_d_lum

INTEGER                           :: j             ! radial zone index
INTEGER                           :: jd            ! particular radial zone index
INTEGER                           :: n             ! neutrino flavor index
INTEGER                           :: itime         ! used to determine when to write a file
INTEGER                           :: itimeprev     ! used to determine when to write a file
INTEGER                           :: istat         ! open and close file flag
INTEGER                           :: nlum_test     ! test for editing data
INTEGER, PARAMETER                :: n_zero = 0

REAL(KIND=double)                 :: tmult         ! used to determine when to write a file

REAL(KIND=double)                 :: area_jd       ! area at jd
REAL(KIND=double)                 :: area_jdm      ! area at zone jd-1
REAL(KIND=double)                 :: r_lumjd       ! luminosity at zone jd (determined by r_lum)
REAL(KIND=double)                 :: r_lumjdm      ! luminosity at zone jd-1 (determined by r_lum)
REAL(KIND=double)                 :: d_lumjd       ! luminosity at zone jd (determined by d_lum)
REAL(KIND=double)                 :: d_lumjdm      ! luminosity at zone jd-1 (determined by d_lum)

REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: lum_r   ! luminosity at radius r_lum
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: lum_rho ! luminosity at density d_lum

REAL(KIND=double), DIMENSION(nnu,ny,nz) :: lum_r_a   ! luminosity at radius r_lum
REAL(KIND=double), DIMENSION(nnu,ny,nz) :: lum_rho_a ! luminosity at density d_lum

REAL(KIND=double), DIMENSION(nnu) :: lum_r_bar     ! angular sum of luminosity at radius r_lum
REAL(KIND=double), DIMENSION(nnu) :: lum_rho_bar   ! angular sum of luminosity at density d_lum

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Need MPI version of lumplot since n_proc=',i4,' > 1')
  201 FORMAT ('     time           lum_e        lum_eb       lum_x        lumxb')
  203 FORMAT (es15.8,4es13.4)
 1001 FORMAT (' jd cannot be found in subroutine lumplot for r')
 4001 FORMAT (' jd cannot be found in subroutine lumplot for rho')
 8001 FORMAT (' File nlumplt1 cannot be opened in subroutime lumplot')
 8501 FORMAT (' File nlumplt2 cannot be opened in subroutime lumplot')

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

!----------------------------------------------------------------------!
!        Return if ilumplt = 0.
!----------------------------------------------------------------------!

IF ( ilumplt == 0 ) RETURN

!-----------------------------------------------------------------------
!        Return if no neutrinos.
!-----------------------------------------------------------------------

IF (nnugpmx == 0) RETURN

!-----------------------------------------------------------------------
!  Evaluate criterion for writing entries to luminosity plot files.
!
!  Cycle criterion.
!-----------------------------------------------------------------------

lprint             = .false.
IF ( ncycle > ncylum(2) ) THEN
  nlum_test        = intlum(3)
ELSE IF ( ncycle > ncylum(1) ) THEN
  nlum_test        = intlum(2)
ELSE
  nlum_test        = intlum(1)
END IF
nlum               = nlum + 1
IF ( nlum >= nlum_test ) THEN
  lprint           = .true.
  IF ( ij_ray == ij_ray_dim  .and.  ij_ray == ij_ray_dim ) nlum = n_zero
END IF

!-----------------------------------------------------------------------
!  Time criterion.
!-----------------------------------------------------------------------

tmult              = 1.d+3/dtimeplot
itime              = int( time * tmult )
itimeprev          = int( ( time - dtnph )*tmult )
IF ( itime > itimeprev ) lprint = .true.

IF ( .not. lprint ) RETURN

!-----------------------------------------------------------------------
!
!        \\\\\ MEAN NEUTRINO LUMINOSITIES AT RADIUS RLUM /////
!
!-----------------------------------------------------------------------

lum_r              = zero
lum_rho            = zero

!-----------------------------------------------------------------------
!  Find jd such that r(jd) > r_lum and r(jd-1) < r_lum.
!-----------------------------------------------------------------------

l_r_lum            = .false.
DO j = jr_min,jr_max
  IF ( r(j) >= r_lum ) THEN
    jd             = j
    l_r_lum        = .true.
    EXIT
  END IF ! r(j) > r_lum
END DO

IF ( .not. l_r_lum ) WRITE (nprint,1001)

!-----------------------------------------------------------------------
!  Compute the luminosities at radius r_lum.
!-----------------------------------------------------------------------

IF ( l_r_lum ) THEN

  area_jd          = frpi * r(jd  ) * r(jd  )
  area_jdm         = frpi * r(jd-1) * r(jd-1)
  DO n = 1,nnu
    r_lumjd        = fluxnu(jd  ,n) * area_jd  * ergfoe
    r_lumjdm       = fluxnu(jd-1,n) * area_jdm * ergfoe
    lum_r(n,ij_ray,ik_ray) = rinterp( r_lumjd, r_lumjdm, r(jd), r_lum, r(jd-1) )
  END DO

!-----------------------------------------------------------------------
!  End computation of luminosities at r_lum.
!-----------------------------------------------------------------------

END IF ! l_r_lum

!-----------------------------------------------------------------------
!
!      \\\\\ MEAN NEUTRINO LUMINOSITIES AT DENSITY D_LUM /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find jd such that rho(jd) < d_lum and rho(jd-1) > d_lum.
!-----------------------------------------------------------------------

l_d_lum           = .false.
DO j = jr_min,jr_max
  IF ( rho(j) <= d_lum ) THEN
    jd            = j
    l_d_lum       = .true.
    EXIT
  END IF ! rho(j) <= d_lum
END DO

IF ( .not. l_d_lum ) WRITE (nprint,4001)

!-----------------------------------------------------------------------
!  Compute the luminosities at density d_lum.
!-----------------------------------------------------------------------

IF ( l_d_lum ) THEN

  area_jd          = frpi * r(jd  ) * r(jd  )
  area_jdm         = frpi * r(jd-1) * r(jd-1)
  DO n = 1,nnu
    d_lumjd        = fluxnu(jd  ,n) * area_jd  * ergfoe
    d_lumjdm       = fluxnu(jd-1,n) * area_jdm * ergfoe
    lum_rho(n,ij_ray,ik_ray) = rinterp( d_lumjd, d_lumjdm, rho(jd), d_lum, rho(jd-1) )
  END DO

!-----------------------------------------------------------------------
!  End computation of luminosities at d_lum.
!-----------------------------------------------------------------------

END IF ! l_d_lum

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
!  Angularly average and edit lum_r.
!-----------------------------------------------------------------------

IF ( l_r_lum ) THEN

!-----------------------------------------------------------------------
!  Load angular arrays
!-----------------------------------------------------------------------

  lum_r_a(:,:,:)   = lum_r(:,:,:)

!-----------------------------------------------------------------------
!  Calculate angular means
!-----------------------------------------------------------------------

  DO n = 1, nnu
    lum_r_bar(n)   = SUM( lum_r_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! n = 1, nnu

!-----------------------------------------------------------------------
!  Open lum_r.d file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nlumplt1,FILE=TRIM(data_path)//'/Plot_Files/lum_r.d', &
&  STATUS='new',IOSTAT=istat)
  IF ( istat == 0 ) WRITE (nlumplt1,201)
  IF ( istat /= 0 ) OPEN (UNIT=nlumplt1,FILE=TRIM(data_path)//'/Plot_Files/lum_r.d', &
&  STATUS='old', POSITION='append',IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8001)
    STOP
  END IF
  
!-----------------------------------------------------------------------
!  Write to lum_r.d file.
!-----------------------------------------------------------------------

  WRITE (nlumplt1,203) time,(lum_r_bar(n),n=1,nnu)
  
!-----------------------------------------------------------------------
!  Close lum_r.d file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nlumplt1, STATUS='keep')

!-----------------------------------------------------------------------
!  End editing of luminosities at r_lum.
!-----------------------------------------------------------------------

END IF ! l_r_lum

!-----------------------------------------------------------------------
!  Angularly average and edit lum_rho.
!-----------------------------------------------------------------------

IF ( l_d_lum ) THEN

!-----------------------------------------------------------------------
!  Load angular arrays
!-----------------------------------------------------------------------

  lum_rho_a(:,:,:) = lum_rho(:,:,:)

!-----------------------------------------------------------------------
!  Calculate angular means
!-----------------------------------------------------------------------

  DO n = 1, nnu
    lum_rho_bar(n) = SUM( lum_rho_a(n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! n = 1, nnu

!-----------------------------------------------------------------------
!  Open lum_rho.d file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nlumplt2,FILE=TRIM(data_path)//'/Plot_Files/lum_rho.d', &
&  STATUS='new',IOSTAT=istat)
  IF ( istat == 0 ) WRITE (nlumplt2,201)
  IF ( istat /= 0 ) OPEN (UNIT=nlumplt2,FILE=TRIM(data_path)//'/Plot_Files/lum_rho.d', &
&  STATUS='old', POSITION='append',IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,8501)
    STOP
  END IF
  
!-----------------------------------------------------------------------
!  Write to lum_rho.d file.
!-----------------------------------------------------------------------

  WRITE (nlumplt2,203) time,(lum_rho_bar(n),n=1,nnu)

!-----------------------------------------------------------------------
!  Close lum_rho.d file.
!-----------------------------------------------------------------------

  CLOSE (unit=nlumplt2, status='keep')

!-----------------------------------------------------------------------
!  End editing of luminosities at d_lum.
!-----------------------------------------------------------------------

END IF ! l_d_lum

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

  rinterp      = b + ( a - b ) * ( y - z )/( x - z )

END FUNCTION rinterp

END
