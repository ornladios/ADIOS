SUBROUTINE shock_plot( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& jmin, jmax, kmin, kmax, nx, ny, nz, nnu )
!-----------------------------------------------------------------------
!
!    File:         shock_plot_MPI
!    Module:       shock_plot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/11/00
!
!    Purpose:
!      To edit the mean shock and gain radii to a file at selected
!       intervals.
!
!    Subprograms called:
!      eqstta_x
!
!    Input arguments:
!  jr_min       : shifted minimum radial index
!  jr_max       : shifted maximum radial index
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  jmin         : minimum y-array index for the edit
!  jmax         : maximum y-array index for the edit
!  kmin         : minimum z-array index for the edit
!  kmax         : maximum z-array index for the edit
!  nx           : x_array extent
!  ny           : y_array extent
!  nz           : z_array extent
!  nnu          : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!  nu_dist_module, parallel_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc, n_proc_y, n_proc_z
USE numerical_module, ONLY : zero, half, frpi

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : dtimeplot, rinnerb, nprint, iplotshk, &
& nplotshk, d_omega, data_path
USE eos_snc_x_module, ONLY : nse
USE mdl_cnfg_module, ONLY : r, rstmss
USE nu_dist_module, ONLY : dunujeadt
USE parallel_module, ONLY : myid, ierr
USE t_cntrl_module, ONLY : time, dtnph

USE mpi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                   :: jr_min         ! shifted minimum radial index
INTEGER, INTENT(in)                   :: jr_max         ! shifted maximum radial index
INTEGER, INTENT(in)                   :: ij_ray         ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                   :: ik_ray         ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                   :: ij_ray_dim     ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                   :: ik_ray_dim     ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                   :: jmin           ! minimum y-array index for the edit
INTEGER, INTENT(in)                   :: jmax           ! maximum y-array index for the edit
INTEGER, INTENT(in)                   :: kmin           ! minimum z-array index for the edit
INTEGER, INTENT(in)                   :: kmax           ! maximum z-array index for the edit
INTEGER, INTENT(in)                   :: nx             ! x-array extent
INTEGER, INTENT(in)                   :: ny             ! y-array extent
INTEGER, INTENT(in)                   :: nz             ! z-array extent
INTEGER, INTENT(in)                   :: nnu            ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                               :: j              ! angular zone index
INTEGER                               :: k              ! azimuthal zone index
INTEGER                               :: itime          ! integer of ( time * tmult )
INTEGER                               :: itimeprev      ! previous value of int( time * tmult )
INTEGER                               :: istat          ! open - close flag

INTEGER                               :: m              ! processor index
INTEGER                               :: mj             ! y-block index
INTEGER                               :: mk             ! z-block index
INTEGER                               :: jsk            ! y-array index of gathered array
INTEGER                               :: ksk            ! z-array index of gathered array

INTEGER                               :: j_shock        ! shifted radial zone index of shock maximum
INTEGER                               :: j_shock_mn     ! shifted radial zone index of minimum estimated shock radius
INTEGER                               :: j_shock_mx     ! shifted radial zone index of maximum estimated shock radius
INTEGER                               :: jr_nse         ! shifted radial zone index at the nse - nonnse boundary

INTEGER                               :: c_gath_recv    ! gather recv buffer count
INTEGER                               :: c_gath_send    ! gather send buffer count

REAL(KIND=double)                     :: tmult          ! used to determine when to write a file

REAL(KIND=double), DIMENSION(nx)      :: rjmh           ! zone centered radii
REAL(KIND=double), DIMENSION(nx)      :: mjmh           ! zone centered enclosed rest mass

REAL(KIND=double), PARAMETER          :: pqcrit = 0.2d0 ! minimum pq_x/p for presence of shock

REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim) :: r_shk    ! quadraticly interpolated shock radius
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim) :: m_shk    ! quadraticly interpolated mass enclosed by shock
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim) :: rgainemx ! maximum gain radius for nue's
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim) :: rgainemn ! minimum gain radius for nue's
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim) :: rgainamx ! maximum gain radius for nuebar's
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim) :: rgainamn ! minimum gain radius for nuebar's
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim) :: r_nse    ! radius of the nse-nonnse boundary

REAL(KIND=double), DIMENSION(ny,nz)   :: r_shk_a        ! quadraticly interpolated shock radius
REAL(KIND=double), DIMENSION(ny,nz)   :: m_shk_a        ! quadraticly interpolated mass enclosed by shock
REAL(KIND=double), DIMENSION(ny,nz)   :: rgainemx_a     ! maximum gain radius for nue's
REAL(KIND=double), DIMENSION(ny,nz)   :: rgainemn_a     ! minimum gain radius for nue's
REAL(KIND=double), DIMENSION(ny,nz)   :: rgainamx_a     ! maximum gain radius for nuebar's
REAL(KIND=double), DIMENSION(ny,nz)   :: rgainamn_a     ! minimum gain radius for nuebar's
REAL(KIND=double), DIMENSION(ny,nz)   :: r_nse_a        ! radius of the nse-nonnse boundary

REAL(KIND=double)                     :: r_shk_bar      ! angular mean of quadraticly interpolated shock radius
REAL(KIND=double)                     :: m_shk_bar      ! angular mean of quadraticly interpolated mass enclosed by shock
REAL(KIND=double)                     :: rgainemx_bar   ! angular mean of maximum gain radius for nue's
REAL(KIND=double)                     :: rgainemn_bar   ! angular mean of minimum gain radius for nue's
REAL(KIND=double)                     :: rgainamx_bar   ! angular mean of maximum gain radius for nuebar's
REAL(KIND=double)                     :: rgainamn_bar   ! angular mean of minimum gain radius for nuebar's
REAL(KIND=double)                     :: r_nse_bar      ! angular mean of radius of the nse-nonnse boundary

REAL(KIND=double), DIMENSION(7,ij_ray_dim,ik_ray_dim)                   :: send_buf   ! send buffer for gathering data from processors
REAL(KIND=double), DIMENSION(7,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z) :: recv_buf   ! receive buffer for gathering data from processors

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT ('     time         r_shock     m_shock      r_nse     r_gainemn  &
& r_gainemx   r_gainamn   r_gainamx')
  103 FORMAT (es15.8,7(es12.4))
 1001 FORMAT (' File nplotshk cannot be opened in subroutime shock_plot')

!-----------------------------------------------------------------------
!
!           \\\\\ CRITERIA FOR WRITING TO SHOCK FILES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Evaluate criterion for writing entries to shock file.
!-----------------------------------------------------------------------

IF ( iplotshk == 0 ) RETURN

!-----------------------------------------------------------------------
!  Evaluate time criterion for writing to plot files.
!-----------------------------------------------------------------------

tmult              = 1.d+3/dtimeplot
itime              = INT( time * tmult )
itimeprev          = INT( ( time - dtnph ) * tmult )
if ( itime == itimeprev ) RETURN

!-----------------------------------------------------------------------
!
!                   \\\\\ SHOCK AND GAIN RADII /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Compute the zone-centered radii and masses.
!-----------------------------------------------------------------------

rjmh(2:jr_max)     = half * ( r     (2:jr_max) + r     (1:jr_max-1) )
mjmh(2:jr_max)     = half * ( rstmss(2:jr_max) + rstmss(1:jr_max-1) )

!-----------------------------------------------------------------------
!  Find shock by quadratic interpolation of pq_x(:,i_ray)/aesv(:,1,i_ray)
!-----------------------------------------------------------------------

CALL findshock( jr_min, jr_max, ij_ray, ik_ray, pqcrit, j_shock, j_shock_mx, &
& j_shock_mn, nx, rjmh, mjmh, r_shk(ij_ray,ik_ray), m_shk(ij_ray,ik_ray) )

!-----------------------------------------------------------------------
!  Find nse - nonnse boundary
!-----------------------------------------------------------------------

jr_nse             = 1
DO j = jr_max, jr_min, -1
  IF ( nse(j,ij_ray,ik_ray) == 1 ) THEN
    jr_nse         = j
    EXIT
  END IF ! nse(j,i_ray) == 1
END DO ! j = jr_max, jr_min, -1
r_nse(ij_ray,ik_ray) = r(jr_nse)
    
!-----------------------------------------------------------------------
!        e-neutrino and e-antineutrino gain radii.
!-----------------------------------------------------------------------

rgainemx(ij_ray,ik_ray) = zero
DO j = j_shock_mn,2,-1
  IF ( dunujeadt (j,1,ij_ray,ik_ray) > zero  .and.  dunujeadt (j+1,1,ij_ray,ik_ray) < zero ) THEN
    rgainemx(ij_ray,ik_ray) = r(j)
    EXIT
  END IF ! dunujeadt (j,1,ij_ray,ik_ray) > 0, dunujeadt (j+1,1,ij_ray,ik_ray) < 0
END DO ! jshock(2),2,-1

rgainemn(ij_ray,ik_ray) = zero
DO j = j_shock_mn,4,-1
  IF ( dunujeadt (j-2,1,ij_ray,ik_ray) > zero     .and. &
&      dunujeadt (j-1,1,ij_ray,ik_ray) > zero     .and. &
&      dunujeadt (j  ,1,ij_ray,ik_ray) > zero     .and. &
&      dunujeadt (j+1,1,ij_ray,ik_ray) < zero           ) THEN
    rgainemn(ij_ray,ik_ray) = r(j)
    EXIT
  END IF ! dunujeadt
END DO ! j = jshock(2),4,-1

rgainamx(ij_ray,ik_ray) = zero
DO j = j_shock_mn,2,-1
  IF ( dunujeadt (j,2,ij_ray,ik_ray) > zero  .and.  dunujeadt (j+1,2,ij_ray,ik_ray) < zero ) THEN
  rgainamx(ij_ray,ik_ray) = r(j)
    EXIT
  END IF ! dunujeadt (j,2),ij_ray,ik_ray > 0, dunujeadt (j+1,2,ij_ray,ik_ray) < 0 
END DO ! j = jshock(2),2,-1

rgainamn(ij_ray,ik_ray) = zero
DO j = j_shock_mn,4,-1
  IF ( dunujeadt (j-2,2,ij_ray,ik_ray) > zero     .and. &
&      dunujeadt (j-1,2,ij_ray,ik_ray) > zero     .and. &
&      dunujeadt (j  ,2,ij_ray,ik_ray) > zero     .and. &
&      dunujeadt (j+1,2,ij_ray,ik_ray) < zero           ) THEN
    rgainamn(ij_ray,ik_ray) = r(j)
    EXIT
  END IF ! dunujeadt
END DO ! jshock(2),4,-1

!-----------------------------------------------------------------------
!
!   \\\\\ ANGULARLY AVERAGE SHOCK AND GAIN RADII DATA AND EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if ij_ray /= ij_ray_dim or ij_ray /= ij_ray_dim
!-----------------------------------------------------------------------

IF ( ij_ray /= ij_ray_dim  .or.  ij_ray /= ij_ray_dim ) RETURN

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

send_buf(1,:,:)    = r_shk   (:,:)
send_buf(2,:,:)    = m_shk   (:,:)
send_buf(3,:,:)    = rgainemx(:,:)
send_buf(4,:,:)    = rgainemn(:,:)
send_buf(5,:,:)    = rgainamx(:,:)
send_buf(6,:,:)    = rgainamn(:,:)
send_buf(7,:,:)    = r_nse   (:,:)

!-----------------------------------------------------------------------
!  Gather quantities
!-----------------------------------------------------------------------

c_gath_recv        = 7 * ij_ray_dim * ik_ray_dim
c_gath_send        = 7 * ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( send_buf, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Unpack receive buffer
!-----------------------------------------------------------------------

  DO m = 0,n_proc-1
    mj                    = MOD( m, n_proc_y )
    mk                    = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                 = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk               = mj * ij_ray_dim + j
        r_shk_a   (jsk,ksk) = recv_buf(1,j,k,m+1)
        m_shk_a   (jsk,ksk) = recv_buf(2,j,k,m+1)
        rgainemx_a(jsk,ksk) = recv_buf(3,j,k,m+1)
        rgainemn_a(jsk,ksk) = recv_buf(4,j,k,m+1)
        rgainamx_a(jsk,ksk) = recv_buf(5,j,k,m+1)
        rgainamn_a(jsk,ksk) = recv_buf(6,j,k,m+1)
        r_nse_a   (jsk,ksk) = recv_buf(7,j,k,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!  Calculate angular means
!-----------------------------------------------------------------------

  r_shk_bar        = SUM( r_shk_a   (jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  m_shk_bar        = SUM( m_shk_a   (jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  rgainemx_bar     = SUM( rgainemx_a(jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  rgainemn_bar     = SUM( rgainemn_a(jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  rgainamx_bar     = SUM( rgainamx_a(jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  rgainamn_bar     = SUM( rgainamn_a(jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  r_nse_bar        = SUM( r_nse_a   (jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi

!-----------------------------------------------------------------------
!  Open shock file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nplotshk, FILE=TRIM(data_path)//'/Plot_Files/shock.d', STATUS='new', IOSTAT=istat)
  IF ( istat == 0 ) WRITE (nplotshk,101)
  IF ( istat /= 0 ) OPEN (UNIT=nplotshk, FILE=TRIM(data_path)//'/Plot_Files/shock.d', &
&  STATUS='old', POSITION='append', IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (nprint,1001)
    STOP
  END IF ! istat /= 0

!-----------------------------------------------------------------------
!  Write to shock file.
!-----------------------------------------------------------------------

  WRITE (nplotshk,103) time, r_shk_bar, m_shk_bar, r_nse_bar, rgainemn_bar, &
&  rgainemx_bar, rgainamn_bar, rgainamx_bar

!-----------------------------------------------------------------------
!  Close shock file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nplotshk, STATUS='keep')

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

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

rinterp            = b + ( a - b ) * ( y - z )/( x - z )
END FUNCTION rinterp

END SUBROUTINE shock_plot
