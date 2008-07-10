SUBROUTINE nuplot( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& jmin, jmax, kmin, kmax, ny, nz, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         nuplot_MPI
!    Module:       nuplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/15/98
!
!    Purpose:
!      To create files of the angularly averaged psi0 and psi1 for each
!       energy group at selected times for post processing.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  jr_min         : minimum radial zone index
!  jr_max         : maximum radial zone index
!  ij_ray         : index denoting the j-index of a specific radial ray
!  ik_ray         : index denoting the k-index of a specific radial ray
!  ij_ray_dim     : number of y-zones on a processor before swapping
!  ik_ray_dim     : number of z-zones on a processor before swapping
!  jmin           : minimum y-array index for the edit
!  jmax           : maximum y-array index for the edit
!  kmin           : minimum z-array index for the edit
!  kmax           : maximum z-array index for the edit
!  ny             : y_array extent
!  nz             : z_array extent
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!  dtnudata       : time (milliseconds) between successive data dumps.
!  r_nudata       : radius at which psi0 and psi1 are evaluated and
!  t_nudata       : time when the last dump was made.
!                    printed.
!  inudata        : 0, omit dumping neutrino data,
!                   1, dump neutrino data.
!  nnudata        : unit number for the dump files.
!  psi0dat(j,k,ij_ray,ik_ray) : time integrated psi0, to be time averaged
!   on dumping.
!  psi1dat(j,k,ij_ray,ik_ray) : time integrated psi1, to be time averaged
!   on dumping.
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, mdl_cnfg_module, nu_dist_module, nu_energy_grid_module,
!  parallel_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc, n_proc_y, n_proc_z
USE numerical_module, ONLY : zero, frpi

USE edit_module, ONLY : nprint, head, inudata, nnudata, dtnudata, t_nudata, &
& r_nudata, psi0dat, psi1dat, d_omega, data_path
USE mdl_cnfg_module, ONLY : r
USE nu_dist_module, ONLY : psi0, psi1, runu, rjmh
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, unubi, nnugpmx
USE parallel_module, ONLY : myid, ierr
USE t_cntrl_module, ONLY : dtnph, time

USE mpi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                   :: jr_min      ! minimum radial zone index
INTEGER, INTENT(in)                   :: jr_max      ! maximum radial zone index
INTEGER, INTENT(in)                   :: ij_ray      ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                   :: ik_ray      ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                   :: ij_ray_dim  ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                   :: ik_ray_dim  ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                   :: jmin        ! minimum y-array index for the edit
INTEGER, INTENT(in)                   :: jmax        ! maximum y-array index for the edit
INTEGER, INTENT(in)                   :: kmin        ! minimum z-array index for the edit
INTEGER, INTENT(in)                   :: kmax        ! maximum z-array index for the edit
INTEGER, INTENT(in)                   :: ny          ! y-array extent
INTEGER, INTENT(in)                   :: nz          ! z-array extent
INTEGER, INTENT(in)                   :: nez         ! neutrino energy array extent
INTEGER, INTENT(in)                   :: nnu         ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                               :: first = .true.

INTEGER                               :: k           ! neutrino energy zone index
INTEGER                               :: itime       ! integer of ( time * tmult )
INTEGER                               :: itimeprev   ! previous value of int( time * tmult )
INTEGER                               :: j           ! radial zone index
INTEGER                               :: jd          ! r(jd) > r_nudata and r(jd-1) < r_nudata
INTEGER                               :: jdmh        ! rjmh(jdmh) > r_nudata and rjmh(jdmh-1) < r_nudata
INTEGER                               :: n           ! neutrino flavor index
INTEGER                               :: istat       ! open - close flag

INTEGER                               :: m           ! processor index
INTEGER                               :: mj          ! y-block index
INTEGER                               :: mk          ! z-block index
INTEGER                               :: jsk         ! y-array index of gathered array
INTEGER                               :: ksk         ! z-array index of gathered array

INTEGER                               :: c_gath_recv ! gather recv buffer count
INTEGER                               :: c_gath_send ! gather send buffer count

REAL(KIND=double)                     :: tmult       ! used to determine when to write a file
REAL(KIND=double)                     :: r0_scale    ! scale factor for extrapolating psi0
REAL(KIND=double)                     :: r1_scale    ! scale factor for extrapolating psi1
REAL(KIND=double)                     :: psi0_rnu    ! psi0 interpolated or extrapolated to r_nudata
REAL(KIND=double)                     :: psi1_rnu    ! psi1 interpolated or extrapolated to r_nudata

REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: psi0dat_a ! full angular array of time integrated psi0
REAL(KIND=double), DIMENSION(nez,nnu,ny,nz) :: psi1dat_a ! full angular array of time integrated psi1

REAL(KIND=double), DIMENSION(nez,nnu) :: psi0dat_bar     ! angular sum of time integrated psi0
REAL(KIND=double), DIMENSION(nez,nnu) :: psi1dat_bar     ! angular sum of time integrated psi1

REAL(KIND=double), DIMENSION(nez,nnu) :: psi0ave         ! mean psi0 interpolated or extrapolated to r_nudata
REAL(KIND=double), DIMENSION(nez,nnu) :: psi1ave         ! mean psi1 interpolated or extrapolated to r_nudata

REAL(KIND=double), DIMENSION(2,nez,nnu,ij_ray_dim,ik_ray_dim)                   :: send_buf   ! send buffer for gathering data from processors
REAL(KIND=double), DIMENSION(2,nez,nnu,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z) :: recv_buf   ! receive buffer for gathering data from processors

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  201 format (es15.8)
  203 format (10(es11.3))
  801 format (1x,a128)
  803 format (1x,128('*')/)
  805 format (10x,' Neutrino group energies')
  807 format (10x,32('-')/)
  809 format (' unui(',i2,')=',es11.4,10x,' unubi(',i2,')=',es11.4,10x, &
& 'dunui(',i2,')=',es11.4)
  811 format (31x,' unubi(',i2,')=',es11.4)
  813 format (' runu=',es11.4//)
 2001 format (' jd cannot be found in subroutine nuplot for r_nudata')
 2003 format (' jdmh cannot be found in subroutine nuplot for r_nudata')
 8001 format (' File psi0.d cannot be opened in subroutime nuplot')
 8501 format (' File psi1.d cannot be opened in subroutime nuplot')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return inudata = 0.
!-----------------------------------------------------------------------

IF ( inudata == 0 ) RETURN

!-----------------------------------------------------------------------
!  Return if no neutrinos.
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!  Zero arrays if first = true.
!-----------------------------------------------------------------------

IF ( first ) THEN
  psi0dat             = zero
  psi1dat             = zero
  first               = .false.
END IF ! first

!-----------------------------------------------------------------------
!
!                    \\\\\ INITIALIZE FILES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Print header and energy group data if opening files for first time.
!
!  Open psi0.d file.
!-----------------------------------------------------------------------

OPEN (UNIT=nnudata,FILE=TRIM(data_path)//'/Plot_Files/psi0.d', &
& STATUS='new',IOSTAT=istat)

IF ( istat == 0 ) THEN
  WRITE (nnudata,801) head
  WRITE (nnudata,803)
  WRITE (nnudata,805)
  WRITE (nnudata,807)
  DO k = 1,nnugpmx 
    WRITE (nnudata,809) k,unui(k),k,unubi(k),k,dunui(k)
  END DO
  k                   = nnugpmx + 1
  WRITE (nnudata,811) k,unubi(k)
  WRITE (nnudata,813) runu

!-----------------------------------------------------------------------
!  Close psi0.d file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nnudata, STATUS='keep')
END IF

!-----------------------------------------------------------------------
!  Open psi1.d file.
!-----------------------------------------------------------------------

OPEN (UNIT=nnudata,FILE=TRIM(data_path)//'/Plot_Files/psi1.d', &
& STATUS='new',IOSTAT=istat)

IF ( istat .eq. 0 ) THEN
  WRITE (nnudata,801) head
  WRITE (nnudata,803)
  WRITE (nnudata,805)
  WRITE (nnudata,807)
  DO k = 1,nnugpmx
    WRITE (nnudata,809) k,unui(k),k,unubi(k),k,dunui(k)
  END DO ! k = 1,nnugpmx
  k                   = nnugpmx + 1
  WRITE (nnudata,811) k,unubi(k)
  WRITE (nnudata,813) runu

!-----------------------------------------------------------------------
!  Close psi1.d file.
!-----------------------------------------------------------------------

  CLOSE (unit=nnudata, status='keep')
END IF

!-----------------------------------------------------------------------
!
!             \\\\\ PSI0 AND PSI1 AT RADIUS R_NUDATA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine psi0 and psi1 at r_nudata.
!
!  If rjmh(jr_max-1) < r_nudata, extrapolate psi0 and psi1 to r_nudata.
!-----------------------------------------------------------------------

IF ( rjmh(jr_max-1) < r_nudata ) THEN

  r0_scale            = rjmh(jr_max-1)**2/r_nudata**2
  r1_scale            = r   (jr_max-1)**2/r_nudata**2
  DO k = 1,nnugpmx
    DO n = 1,nnu
      psi0_rnu        = r0_scale * psi0(jr_max-1,k,n)
      psi0dat(k,n,ij_ray,ik_ray) = psi0dat(k,n,ij_ray,ik_ray) + psi0_rnu * dtnph
      psi1_rnu        = r1_scale * psi1(jr_max-1,k,n)
      psi1dat(k,n,ij_ray,ik_ray) = psi1dat(k,n,ij_ray,ik_ray) + psi1_rnu * dtnph
    END DO ! n = 1,nnu
  END DO ! k = 1,nnugpmx

ELSE ! rjmh(jr_max-1) > r_nudata

!-----------------------------------------------------------------------
!  Find jd such that r(jd) > r_nudata and r(jd-1) < r_nudata.
!-----------------------------------------------------------------------

  jd                  = 0
  DO j = jr_min,jr_max
    IF ( r(j) >= r_nudata ) THEN
      jd              = j
      EXIT
    END IF ! r(j) > r_nudata
  END DO
  
  IF ( jd == 0 ) THEN
    WRITE (nprint,2001)
    RETURN
  END IF

!-----------------------------------------------------------------------
! Find jdmh such that rjmh(jdmh) > r_nudata and rjmh(jdmh-1) < r_nudata.
!-----------------------------------------------------------------------

  jdmh                = 0
  DO j = jr_min,jr_max
    IF ( rjmh(j) >= r_nudata ) THEN
      jdmh            = j
      EXIT
    END IF ! rjmh(j) > r_nudata
  END DO ! j = 2,jr_max

  IF ( jdmh == 0 ) THEN
    WRITE (nprint,2003)
    RETURN
  END IF

!-----------------------------------------------------------------------
!  Update neutrino numbers radiated across r_nudata.
!-----------------------------------------------------------------------

  DO k = 1,nnugpmx
    DO n = 1,nnu
      psi0_rnu        = rinterp( psi0(jdmh,k,n), psi0(jdmh-1,k,n),  &
&                      rjmh(jdmh) * rjmh(jdmh), r_nudata * r_nudata, &
&                      rjmh(jdmh-1) * rjmh(jdmh-1) )
      psi0dat(k,n,ij_ray,ik_ray) = psi0dat(k,n,ij_ray,ik_ray) + psi0_rnu * dtnph
      psi1_rnu        = rinterp( psi1(jd,k,n), psi1(jd-1,k,n), r(jd) * r(jd), &
&                      r_nudata*r_nudata, r(jd-1) * r(jd-1) )
      psi1dat(k,n,ij_ray,ik_ray) = psi1dat(k,n,ij_ray,ik_ray) + psi1_rnu * dtnph
    END DO ! n = 1,nnu
  END DO ! k = 1,nnugpmx

END IF ! rjmh(jr_max-1) < r_nudata

!-----------------------------------------------------------------------
!
!           \\\\\ CRITERIA FOR WRITING TO NUPLOT FILES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Evaluate criterion for writing to nuplot files.
!-----------------------------------------------------------------------

tmult                 = 1.d+3/dtnudata
itime                 = int( time * tmult )
itimeprev             = int( ( time - dtnph ) * tmult )
IF ( itime == itimeprev ) RETURN

!-----------------------------------------------------------------------
!
!         \\\\\ ANGULARLY AVERAGE NEUTRINO DATA AND EDIT /////
!
!-----------------------------------------------------------------------

IF ( ij_ray /= ij_ray_dim  .or.  ij_ray /= ij_ray_dim ) RETURN

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

send_buf(1,:,:,:,:)   = psi0dat(:,:,:,:)
send_buf(2,:,:,:,:)   = psi1dat(:,:,:,:)

!-----------------------------------------------------------------------
!  Gather quantities
!-----------------------------------------------------------------------

c_gath_recv           = 2 * nez * nnu * ij_ray_dim * ik_ray_dim
c_gath_send           = 2 * nez * nnu * ij_ray_dim * ik_ray_dim

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
    mj                = MOD( m, n_proc_y )
    mk                = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk             = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk           = mj * ij_ray_dim + j
        psi0dat_a(:,:,jsk,ksk) = recv_buf(1,:,:,j,k,m+1)
        psi1dat_a(:,:,jsk,ksk) = recv_buf(2,:,:,j,k,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!  Calculate angular means
!-----------------------------------------------------------------------

  DO k = 1,nnugpmx
    DO n = 1,nnu
      psi0dat_bar(k,n)  = SUM( psi0dat_a(k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
      psi1dat_bar(k,n)  = SUM( psi1dat_a(k,n,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
    END DO ! n = 1,nnu
  END DO ! k = 1,nnugpmx

!-----------------------------------------------------------------------
!  Average psi0 and psi1 over time dtnudata.
!-----------------------------------------------------------------------

  DO n = 1,nnu
    DO k = 1,nnugpmx
      psi0ave(k,n)    = psi0dat_bar(k,n)/( time - t_nudata )
      psi1ave(k,n)    = psi1dat_bar(k,n)/( time - t_nudata )
    END DO ! k = 1,nnugpmx
  END DO ! n = 1,nnu

  t_nudata            = time

!-----------------------------------------------------------------------
!  Open psi0.d file.
!-----------------------------------------------------------------------

  OPEN (UNIT=nnudata,FILE=TRIM(data_path)//'/Plot_Files/psi0.d', &
&  STATUS='old',POSITION='append', IOSTAT=istat)

  IF ( istat /= 0 ) THEN
    WRITE (nprint,8001)
    STOP
  END IF

!-----------------------------------------------------------------------
!  Write to psi0.d file.
!-----------------------------------------------------------------------

  WRITE (nnudata,201) time
  DO n = 1,nnu
    WRITE (nnudata,203) (psi0ave(k,n),k=1,10)
    WRITE (nnudata,203) (psi0ave(k,n),k=11,20)
  END DO ! n = 1,nnu

!-----------------------------------------------------------------------
!  Close psi0.d file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nnudata, STATUS='keep')

!-----------------------------------------------------------------------
!  Open psi1.d file.
!-----------------------------------------------------------------------

  OPEN (unit=nnudata,file=TRIM(data_path)//'/Plot_Files/psi1.d', &
&  STATUS='old',POSITION='append', IOSTAT=istat)

  IF ( istat /= 0 ) THEN
    WRITE (nprint,8501)
    STOP
  END IF

!-----------------------------------------------------------------------
!  Write to psi1.d file.
!-----------------------------------------------------------------------

  WRITE (nnudata,201) time
  DO n = 1,nnu
    WRITE (nnudata,203) (psi0ave(k,n),k=1,10)
    WRITE (nnudata,203) (psi0ave(k,n),k=11,20)
  END DO ! n = 1,nnu

!-----------------------------------------------------------------------
!  Close psi0.d file.
!-----------------------------------------------------------------------

  CLOSE (UNIT=nnudata, STATUS='keep')

!-----------------------------------------------------------------------
!                ||||| All processors used now |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

psi0dat               = zero
psi1dat               = zero

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

END SUBROUTINE nuplot
