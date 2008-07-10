SUBROUTINE nuradplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, &
& kmin, kmax, ny, nz, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         nuradplot
!    Module:       nuradplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To create files the angularly averaged cumulative neutrino number
!       radiated by each energy group.
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
!  nez          : neutrino energy array extent
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
USE numerical_module, ONLY : zero, half, frpi, epsilon, ncoef, ecoef

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nu_r, nu_rt, nu_rho, nu_rhot, unu_r, unu_rt, &
& unu_rho, unu_rhot, nprint, iplotnurad, nplotnurad, dtnuradplot,    &
& rho_nurad, r_nurad, d_omega, nlog, data_path
USE mdl_cnfg_module, ONLY : jr_max, r, rho
USE nu_dist_module, ONLY : psi1
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, nnugpmx
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
INTEGER, INTENT(in)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                           :: first = .true.

INTEGER                           :: j             ! radial zone index
INTEGER                           :: jd            ! particular radial zone index
INTEGER                           :: jdm2          ! jd - 2
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index
INTEGER                           :: istat         ! open and close file flag

INTEGER                           :: itime         ! parameter for criterion for writing to nuradplot files
INTEGER                           :: itimeprev     ! parameter for criterion for writing to nuradplot files

REAL(KIND=double)                 :: tmult         ! multiplier for criterion for writing to nuradplot files

REAL(KIND=double)                 :: rjdmh         ! radius midway between zone jd-1 and jd
REAL(KIND=double)                 :: rjdm3h        ! radius midway between zone jd-2 and jd-1
REAL(KIND=double)                 :: psi1mh        ! psi1 midway between zone jd-1 and jd
REAL(KIND=double)                 :: psi1m3h       ! psi1 midway between zone jd-2 and jd-1
REAL(KIND=double)                 :: areajd        ! area at radius r(jd)
REAL(KIND=double)                 :: areajdm       ! area at radius r(jd-1)
REAL(KIND=double)                 :: nuradjd       ! neutrino across areajd in time dtnph
REAL(KIND=double)                 :: nuradjdm      ! neutrino across areajdm in time dtnph
REAL(KIND=double)                 :: dnu_r         ! updated neutrino number radiated across r_nurad
REAL(KIND=double)                 :: dnu_rho       ! updated neutrino number radiated across rho_nurad

REAL(KIND=double)                 :: dtime         ! time between plot file dumps
REAL(KIND=double), DIMENSION(500) :: ncoefap       ! 4.*pi*c/((h*c)**3)*w**2*dw  
REAL(KIND=double), DIMENSION(500) :: ecoefap       ! 4.*pi*c*ergmev/((h*c)**3)*w**3   

REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: nu_r_a     ! full angular array of number of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: nu_rt_a    ! full angular array of cumulative number of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: nu_rho_a   ! full angular array of number of neutrinos radiated across rho_nurad
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: nu_rhot_a  ! full angular array of cumulative number of neutrinos radiated across rho_nurad
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: unu_r_a    ! full angular array of energy of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: unu_rt_a   ! full angular array of cumulative energy of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: unu_rho_a  ! full angular array of energy of neutrinos radiated across rho_nurad
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: unu_rhot_a ! full angular array of cumulative energy of neutrinos radiated across rho_nurad

REAL(KIND=double), DIMENSION(nez,nnu) :: nu_r_bar         ! angular sum of number of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu) :: nu_rt_bar        ! angular sum of cumulative number of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu) :: nu_rho_bar       ! angular sum of number of neutrinos radiated across rho_nurad
REAL(KIND=double), DIMENSION(nez,nnu) :: nu_rhot_bar      ! angular sum of cumulative number of neutrinos radiated across rho_nurad
REAL(KIND=double), DIMENSION(nez,nnu) :: unu_r_bar        ! angular sum of energy of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu) :: unu_rt_bar       ! angular sum of cumulative energy of neutrinos radiated across r_nurad
REAL(KIND=double), DIMENSION(nez,nnu) :: unu_rho_bar      ! angular sum of energy of neutrinos radiated across rho_nurad
REAL(KIND=double), DIMENSION(nez,nnu) :: unu_rhot_bar     ! angular sum of cumulative energy of neutrinos radiated across rho_nurad

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 format (' Need MPI version of nuradplot since n_proc=',i4,' > 1')
  201 FORMAT (1x/20x,' Differential Number of Neutrinos Radiated'/)
  203 FORMAT ('   k   w (MeV)   dw (MeV)    time (s)  dtime (s)   dN_nue &
&    dN_nuebar   dN_nux   dN_nuxbar   dNd_nue  dNd_nuebar  dNc_nux   dNc_nuxbar'/)
  205 FORMAT (1x,i3,12(es11.3))
  211 FORMAT (1x/20x,' Cumulative Number of Neutrinos Radiated'/)
  213 FORMAT ('   k   w (MeV)   dw (MeV)    time (s)  dtime (s)    N_nue &
&    N_nuebar     N_nux    N_nuxbar    Nd_nue   Nd_nuebar    Nd_nux   Nd_nuxbar'/)
  215 FORMAT (1x,i3,12(es11.3))
  221 FORMAT (1x/20x,' Differential Neutrino Energy Spectrum'/)
  223 FORMAT ('   k   w (MeV)   dw (MeV)    time (s)  dtime (s)   dE_nue &
&    dE_nuebar   dE_nux    dE_nuxbar  dEd_nue  dEd_nuebar   dEd_nux  dEd_nuxbar'/)
  225 FORMAT (1x,i3,12(es11.3))
  231 FORMAT (1x/20x,' Cumulative Neutrino Energy Spectrum'/)
  233 FORMAT ('   k   w (MeV)   dw (MeV)    time (s)  dtime (s)    E_nue &
&    E_nuebar     E_nux   E_nuxbar   Ed_nue   Ed_nuebar    Ed_nux   Ed_nuxbar'/)
  235 FORMAT (1x,i3,12(es11.3))
 1001 FORMAT (' jd cannot be found in subroutine nuradplot for r')
 3001 FORMAT (' jd cannot be found in subroutine nuradplot for rho')
 8001 FORMAT (' File nurad.d cannot be opened in subroutime nuradplot')
 8501 FORMAT (' File unurad.d cannot be opened in subroutime nuradplot')

!-----------------------------------------------------------------------
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
!  Return iplotnurad = 0
!-----------------------------------------------------------------------

IF ( iplotnurad == 0 ) RETURN

!-----------------------------------------------------------------------
!  Return if no neutrinos
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!                 \\\\\ UPDATE NEUTRINO NUMBERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  ncoefap             = zero
  ecoefap             = zero
  ncoefap(1:nnugpmx)  = ncoef * unui(1:nnugpmx)**2 * dunui(1:nnugpmx)
  ecoefap(1:nnugpmx)  = ecoef * unui(1:nnugpmx)**3
  first               = .false.
END IF ! first

!-----------------------------------------------------------------------
!  Zero arrays if ncycle = 0
!-----------------------------------------------------------------------

IF ( ncycle == 0 ) THEN
  nu_r                = zero
  nu_rt               = zero
  nu_rho              = zero
  nu_rhot             = zero
END IF ! ncycle = 0

!-----------------------------------------------------------------------
!  Find jd such that r(jd) > r_nurad and r(jd-1) < r_nurad
!-----------------------------------------------------------------------

DO j = 2,jr_max
  IF ( r(j) >= r_nurad ) THEN
    jd                = j
    EXIT
  END IF ! r(j) > r_nurad
END DO ! j

IF ( jd == jr_max ) THEN
  WRITE (nprint,1001)
  GO TO 2500
END IF ! jd = jr_max!

!-----------------------------------------------------------------------
!  Update neutrino numbers radiated across r_nurad
!-----------------------------------------------------------------------

areajd                = frpi * r(jd  ) * r(jd  )
areajdm               = frpi * r(jd-1) * r(jd-1)
DO k = 1,nnugpmx
  DO n = 1,nnu
    nuradjd           = ncoefap(k) * psi1(jd  ,k,n) * areajd  * dtnph
    nuradjdm          = ncoefap(k) * psi1(jd-1,k,n) * areajdm * dtnph
    dnu_r             = rinterp( nuradjd, nuradjdm, r(jd), r_nurad, r(jd-1) )
    nu_r (k,n,ij_ray,ik_ray) = nu_r (k,n,ij_ray,ik_ray) + dnu_r
    nu_rt(k,n,ij_ray,ik_ray) = nu_rt(k,n,ij_ray,ik_ray) + dnu_r
  END DO ! n = 1,nnu
END DO ! k = 1,nnugpmx

!-----------------------------------------------------------------------
!  Find jd such that rho(jd  ) < rho_nurad and
!                    rho(jd-1) > rho_nurad
!-----------------------------------------------------------------------
 2500 CONTINUE

DO j = 2,jr_max
  IF ( rho(j) <= rho_nurad ) THEN
    jd                = j
    EXIT
  END IF ! rrho(j) <= rho_nurad
END DO ! j = 2,jr_max

IF ( jd == jr_max ) THEN
  WRITE (nprint,3001)
  RETURN
END IF ! jd = jr_max

!-----------------------------------------------------------------------
!  Update neutrino numbers radiated across rho_nurad
!-----------------------------------------------------------------------

jdm2                  = MAX( jd - 2, 1 )
rjdmh                 = half * ( r(jd  ) + r(jd-1) )
rjdm3h                = half * ( r(jd-1) + r(jdm2) )
areajd                = frpi * rjdmh  * rjdmh
areajdm               = frpi * rjdm3h * rjdm3h

DO k = 1,nnugpmx
  DO n = 1,nnu
    psi1mh            = half * ( psi1(jd  ,k,n) + psi1(jd-1,k,n) )
    psi1m3h           = half * ( psi1(jd-1,k,n) + psi1(jdm2,k,n) )
    nuradjd           = ncoefap(k) * psi1mh  * areajd  * dtnph
    nuradjdm          = ncoefap(k) * psi1m3h * areajdm * dtnph
    dnu_rho           = rinterp( nuradjd, nuradjdm, rho(jd), rho_nurad, rho(jd-1) )
    nu_rho (k,n,ij_ray,ik_ray) = nu_rho (k,n,ij_ray,ik_ray) + dnu_rho
    nu_rhot(k,n,ij_ray,ik_ray) = nu_rhot(k,n,ij_ray,ik_ray) + dnu_rho
  END DO ! n = 1,nnu
END DO ! k = 1,nnugpmx

!-----------------------------------------------------------------------
!
!         \\\\\ CRITERIA FOR WRITING TO NURADPLOT FILES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Evaluate criterion for writing to nuradplot files, return if
!   criterion not satisfied
!-----------------------------------------------------------------------

tmult                 = 1.d+3/dtnuradplot
itime                 = int( time * tmult )
itimeprev             = int( ( time - dtnph ) * tmult )
IF ( itime == itimeprev ) RETURN

!-----------------------------------------------------------------------
!
!         \\\\\ ANGULARLY AVERAGE NEUTRINO DATA AND EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if ij_ray /= ij_ray_dim or ij_ray /= ij_ray_dim
!-----------------------------------------------------------------------

IF ( ij_ray /= ij_ray_dim  .or.  ij_ray /= ij_ray_dim ) RETURN

!-----------------------------------------------------------------------
!  Compute neutrino energy data
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  DO n = 1,nnu
    unu_r(k,n,ij_ray,ik_ray)    = ecoefap(k) * nu_r(k,n,ij_ray,ik_ray)   /( ncoefap(k) + epsilon )
    unu_rt(k,n,ij_ray,ik_ray)   = ecoefap(k) * nu_rt(k,n,ij_ray,ik_ray)  /( ncoefap(k) + epsilon )
    unu_rho(k,n,ij_ray,ik_ray)  = ecoefap(k) * nu_rho(k,n,ij_ray,ik_ray) /( ncoefap(k) + epsilon )
    unu_rhot(k,n,ij_ray,ik_ray) = ecoefap(k) * nu_rhot(k,n,ij_ray,ik_ray)/( ncoefap(k) + epsilon )
  END DO ! n = 1,nnu
END DO ! k = 1,nnugpmx

!-----------------------------------------------------------------------
!  Load angular arrays
!-----------------------------------------------------------------------

nu_r_a    (:,:,:,:)   = nu_r    (:,:,:,:)
nu_rt_a   (:,:,:,:)   = nu_rt   (:,:,:,:)
nu_rho_a  (:,:,:,:)   = nu_rho  (:,:,:,:)
nu_rhot_a (:,:,:,:)   = nu_rhot (:,:,:,:)
unu_r_a   (:,:,:,:)   = unu_r   (:,:,:,:)
unu_rt_a  (:,:,:,:)   = unu_rt  (:,:,:,:)
unu_rho_a (:,:,:,:)   = unu_rho (:,:,:,:)
unu_rhot_a(:,:,:,:)   = unu_rhot(:,:,:,:)

!-----------------------------------------------------------------------
!  Calculate angular means
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  DO n = 1,nnu
    nu_r_bar    (k,n) = SUM( nu_r_a    (k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    nu_rt_bar   (k,n) = SUM( nu_rt_a   (k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    nu_rho_bar  (k,n) = SUM( nu_rho_a  (k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    nu_rhot_bar (k,n) = SUM( nu_rhot_a (k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    unu_r_bar   (k,n) = SUM( unu_r_a   (k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    unu_rt_bar  (k,n) = SUM( unu_rt_a  (k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    unu_rho_bar (k,n) = SUM( unu_rho_a (k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    unu_rhot_bar(k,n) = SUM( unu_rhot_a(k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! n = 1,nnu
END DO ! k = 1,nnugpmx

!-----------------------------------------------------------------------
!  Open nurad file
!-----------------------------------------------------------------------

OPEN ( UNIT=nplotnurad, FILE=TRIM(data_path)//'/Plot_Files/nurad.d', &
& STATUS='new', IOSTAT=istat )
IF ( istat /= 0 ) OPEN ( UNIT=nplotnurad, FILE=TRIM(data_path)//'/Plot_Files/nurad.d', &
& STATUS='old', POSITION='append', IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,8001)
  STOP
END IF ! istat /= 0

!-----------------------------------------------------------------------
!  Write to nurad.d file.
!-----------------------------------------------------------------------

WRITE (nplotnurad,201)
WRITE (nplotnurad,203)

dtime                 = dtnuradplot/1.d+3
DO k = 1,nnugpmx
  WRITE (nplotnurad,205) k, unui(k), dunui(k), time, dtime, nu_r_bar(k,1), &
&  nu_r_bar(k,2), nu_r_bar(k,3), nu_r_bar(k,4), nu_rho_bar(k,1),           &
&  nu_rho_bar(k,2), nu_rho_bar(k,3), nu_rho_bar(k,4)
END DO ! k = 1,nnugpmx

nu_r                  = zero
nu_rho                = zero

WRITE (nplotnurad,211)
WRITE (nplotnurad,213)

DO k = 1,nnugpmx
  WRITE (nplotnurad,215) k, unui(k), dunui(k), time, dtime ,nu_rt_bar(k,1), &
&  nu_rt_bar(k,2), nu_rt_bar(k,3), nu_rt_bar(k,4), nu_rhot_bar(k,1),        &
&  nu_rhot_bar(k,2), nu_rhot_bar(k,3), nu_rhot_bar(k,4)
END DO ! k = 1,nnugpmx

!-----------------------------------------------------------------------
!  Close nurad file
!-----------------------------------------------------------------------

CLOSE (UNIT=nplotnurad, STATUS='keep')

!-----------------------------------------------------------------------
!  Open unurad file
!-----------------------------------------------------------------------

OPEN ( UNIT=nplotnurad,FILE=TRIM(data_path)//'/Plot_Files/unurad.d', &
& STATUS='new', IOSTAT=istat)
IF ( istat /= 0 ) OPEN ( UNIT=nplotnurad, FILE=TRIM(data_path)//'/Plot_Files/unurad.d', &
& STATUS='old', POSITION='append', IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,8501)
  STOP
END IF ! stat /= 0

!-----------------------------------------------------------------------
!  Write to nurad.d file
!-----------------------------------------------------------------------

WRITE (nplotnurad,221)
WRITE (nplotnurad,223)
dtime                 = dtnuradplot/1.d+3

DO k = 1,nnugpmx
  WRITE (nplotnurad,225) k, unui(k), dunui(k), time, dtime, unu_r_bar(k,1), &
&  unu_r_bar(k,2), unu_r_bar(k,3), unu_r_bar(k,4), unu_rho_bar(k,1),        &
&  unu_rho_bar(k,2), unu_rho_bar(k,3), unu_rho_bar(k,4)
END DO ! k = 1,nnugpmx

WRITE (nplotnurad,231)
WRITE (nplotnurad,233)

DO k = 1,nnugpmx
  WRITE (nplotnurad,235) k, unui(k), dunui(k), time, dtime, unu_rt_bar(k,1), &
&  unu_rt_bar(k,2), unu_rt_bar(k,3), unu_rt_bar(k,4), unu_rhot_bar(k,1),     &
&  unu_rhot_bar(k,2), unu_rhot_bar(k,3), unu_rhot_bar(k,4)
END DO ! k = 1,nnugpmx

!-----------------------------------------------------------------------
!  Close unurad file
!-----------------------------------------------------------------------

CLOSE (UNIT=nplotnurad, STATUS='keep')

!-----------------------------------------------------------------------
!  Return
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
