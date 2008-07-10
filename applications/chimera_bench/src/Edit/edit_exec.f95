SUBROUTINE edit_exec( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim,     &
& ik_ray_dim, jmin, jmax, kmin, kmax, i_edit, first, nx, ny, nz, nez, &
& nnu )
!-----------------------------------------------------------------------
!
!    File:         edit_exec
!    Module:       edit_exec
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/28/00
!
!    Purpose:
!      To perform edits of the problem at preselected intervals
!
!    Subprograms called:
!  bnuplot      : writes to the the boundary and shock files
!  complot      : writes to the comparison plot files
!  editc        : performs a configuration edit
!  edit_e       : edits detailed energy conservation data
!  edith        : performs an edit of hydrodynamic variables
!  editma       : performs an edit of mass averaged variables
!  editn        : performs an edit of neutrino energy group quantities
!  editng       : performs an edit of energy integrated neutrino quantities
!  editpmtr     : optionally changes selected edit parameters
!  editps       : edits pressure and stress data
!  editsc       : performs an edit of entropy and chemical potential data
!  editu        : performs an edit of energy data
!  edity        : performs an edit of compositionvariables
!  enuvplot     : writes to neutrino energy moments file
!  lagrangeplot : writes to lagrangian evolution files
!  lumplot      : writes to neutrino luminosities file
!  nuplot       : write to a file of psi0 and psi1 data
!  rnuplot      : writes to a file of radius evolution
!  nuradplot    : writes to a file of neutrino energy loss data
!  rlagplot     : writes to a file of lagrangian mass shell evolution
!  roextr       : finds the maximum and minimum central densities between editc edits
!  varplot      : writes to a file of configuration time slices
!
!    Input arguments:
!  jr_min       : minimum zone index
!  jr_max       : maximum zone index
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  jmin         : minimum y-array index for the edit
!  jmax         : maximum y-array index for the edit
!  kmin         : minimum z-array index for the edit
!  kmax         : maximum z-array index for the edit
!  i_edit       : edit parameter
!  first        : initialization flag
!  nx           : x_array extent
!  ny           : y_array extent
!  nz           : z_array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  cycle_module, edit_module, mdl_cnfg.cmn, nu_dist_module,
!  parallel_module, prb_cntl.cmn, t_cntrl.cmn
!      
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc_y
USE numerical_module, ONLY: zero

USE cycle_module, ONLY: ncycle, ncymax
USE edit_module, ONLY: nprint, iprtbgn, nmodel, iprint, prnttest, nprt, npltf, &
& rhoprint, intprt, iflprt, it_edit, dt_edit, data_path, nlog
USE mdl_cnfg_module, ONLY: rhor, rho, t, ye, r, rstmss
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY: jnumin, jnumax
USE t_cntrl_module, ONLY: time, t_bounce, dtnph

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

LOGICAL, INTENT(in)                      :: first         ! initial model print flag

INTEGER, INTENT(in)                      :: jr_min        ! minimum zone index
INTEGER, INTENT(in)                      :: jr_max        ! maximum zone index
INTEGER, INTENT(in)                      :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)                      :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)                      :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                      :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                      :: jmin          ! minimum y-array index for the edit
INTEGER, INTENT(in)                      :: jmax          ! maximum y-array index for the edit
INTEGER, INTENT(in)                      :: kmin          ! minimum z-array index for the edit
INTEGER, INTENT(in)                      :: kmax          ! maximum z-array index for the edit
INTEGER, INTENT(in)                      :: i_edit        ! print on demand flag
INTEGER, INTENT(in)                      :: nx            ! x-array extent
INTEGER, INTENT(in)                      :: ny            ! y-array extent
INTEGER, INTENT(in)                      :: nz            ! z-array extent
INTEGER, INTENT(in)                      :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)                      :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)                      :: intfile       ! character string containing name of model files
CHARACTER (len=128)                      :: outfile       ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_d     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_t     ! character string containing name of model files
CHARACTER (len=128)                      :: outfile_dt    ! character string containing name of model files

INTEGER                                  :: i             ! do index
INTEGER                                  :: itprint       ! used to at specfied times
INTEGER                                  :: itprintmx     ! used to at specfied times
INTEGER, PARAMETER                       :: iunit=20      ! unit number to print models
INTEGER, PARAMETER                       :: iunitt=24     ! unit number to print models at specified temperatures
INTEGER, PARAMETER                       :: iunitd=24     ! unit number to print models at specified densities
INTEGER                                  :: jmaxp         ! jr_max+1
INTEGER, PARAMETER                       :: n1 = 1        ! neutrino flavor index
INTEGER, PARAMETER                       :: n2 = 2        ! neutrino flavor index
INTEGER, PARAMETER                       :: n3 = 3        ! neutrino flavor index
INTEGER, PARAMETER                       :: n4 = 4        ! neutrino flavor index
INTEGER                                  :: nprint_save   ! temporarily stored value of nprint
INTEGER                                  :: istat         ! open file flag
INTEGER                                  :: j_ray         ! polar index of the radial ray
INTEGER                                  :: k_ray         ! azimuthal index of the radial ray

INTEGER                                  :: itime         ! used to determine when to edit
INTEGER                                  :: itimeprev     ! used to determine when to edit

REAL(KIND=double), DIMENSION(10)         :: tprint        ! times for printing

REAL(KIND=double)                        :: t_tb          ! time from bounce
REAL(KIND=double)                        :: tmult         ! used to determine when to write a file

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 2001 FORMAT (' i_edit =',i3,' does not correspond to any anticipated value')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_ray              = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
k_ray              = ( myid/n_proc_y ) * ik_ray_dim + ik_ray

!-----------------------------------------------------------------------
!
!           \\\\\ INITIALIZE AND PERFORM INITIAL EDIT /////!           
!
!-----------------------------------------------------------------------

IF ( first ) THEN

  nprint_save      = nprint
  jmaxp            = jr_max + 1

  WRITE (intfile,'(a14,i4.4,a1,i4.4,a1,i7.7,a2)') '/Log/intlmodel', &
& j_ray,'_',k_ray,'_',ncycle,'.d'
  intfile          = TRIM(data_path)//TRIM(intfile)
  OPEN (UNIT=iunit,FILE=TRIM(intfile),STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=iunit,FILE=TRIM(intfile), STATUS='old')

  itprint          = 1
  itprintmx        = 2

  DO i = 1,10
    tprint(i)      = 1.d+20
  END DO

!-----------------------------------------------------------------------
!        iprtbgn = 0, ncycle = 0 : abbreviated edit
!        iprtbgn = 0, ncycle > 0 : configuration edit
!        iprtbgn = 1, ncycle = 0 : abbreviated edit
!        iprtbgn = 1, ncycle > 0 : abbreviated edit
!        iprtbgn = 2, ncycle = 0 : full edit
!        iprtbgn = 2, ncycle > 0 : configuration edit
!        iprtbgn = 3, ncycle = 0 : full edit
!        iprtbgn = 3, ncycle > 0 : abreviated edit
!        iprtbgn = 4, ncycle = 0 : full edit
!        iprtbgn = 4, ncycle > 0 : full edit
!-----------------------------------------------------------------------

  IF ( iprint /= 0 ) THEN

    nprint         = iunit

    IF         ( iprtbgn == 4                           .or. &
&                iprtbgn == 3  .and.  ncycle == 0       .or. &
&                iprtbgn == 2  .and.  ncycle == 0              ) THEN

      prnttest     = .false.
      CALL editp ( jr_min, jr_max)
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray)
      CALL edit_e( jr_min, jr_max, ij_ray, ik_ray)
      CALL editma( jr_min, jr_max, ij_ray, ik_ray)
      CALL edith ( jr_min, jr_max, ij_ray, ik_ray)
      CALL editps( jr_min, jr_max, ij_ray, ik_ray)
      CALL editu ( jr_min, jr_max, ij_ray, ik_ray)
      CALL edity ( jr_min, jr_max, ij_ray, ik_ray)
      CALL editsc( jr_min, jr_max, ij_ray, ik_ray)
      CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray)
      CALL editng( n1, jnumin, jnumax, ij_ray, ik_ray)
      CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray)
      CALL editng( n2, jnumin, jnumax, ij_ray, ik_ray)
      CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray)
      CALL editng( n3, jnumin, jnumax, ij_ray, ik_ray)
      CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray)
      CALL editng( n4, jnumin, jnumax, ij_ray, ik_ray)

    ELSE IF    ( iprtbgn == 3  .and.  ncycle > 0        .or. &
&                iprtbgn == 2                           .or. &
&                iprtbgn == 1                                  ) THEN

      prnttest     = .false.
      CALL editp ( jr_min, jr_max )
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edit_e( jr_min, jr_max, ij_ray, ik_ray )
      CALL editma( jr_min, jr_max, ij_ray, ik_ray )
      CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editps( jr_min, jr_max, ij_ray, ik_ray )
      CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
      CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray )

    ELSE

      prnttest     = .false.
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray )

    END IF ! iprtbgn

    nprint         = nprint_save
    CLOSE (UNIT=iunit,STATUS='keep')

    IF ( ncycle == 0 ) THEN
      nmodel(ij_ray, ik_ray)  = nmodel(ij_ray, ik_ray) + 1
    END IF

  END IF ! iprint /= 0

!-----------------------------------------------------------------------
!
!        \\\\\ INITIALIZE AND PERFORM INITIAL PLOT EDITS /////
!
!-----------------------------------------------------------------------

  CALL lumplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, kmin, &
&  kmax, ny, nz, nnu )
  CALL enuvplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, kmin, &
&  kmax, ny, nz, nnu )
  CALL nuradplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, kmin, &
&  kmax, ny, nz, nez, nnu )
  CALL nuplot( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, &
&  jmax, kmin, kmax, ny, nz, nez, nnu )
  CALL shock_plot( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&  jmin, jmax, kmin, kmax, nx, ny, nz, nnu )

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

  IF ( myid == 0  .and.  ij_ray == 1  .and.  ik_ray == 1 ) THEN

    CALL rnuplot
    CALL bnuplot( ij_ray, ik_ray )
    CALL varplot( ij_ray, ik_ray )
    CALL complot( ij_ray, ik_ray )
    CALL lagrangeplot(ij_ray, ik_ray)
    CALL rlagplot

!-----------------------------------------------------------------------
!               ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

  END IF ! myid == 0  .and.  ij_ray == 1  .and.  ik_ray == 1

  RETURN

END IF ! first

!-----------------------------------------------------------------------
!
!             \\\\\ PRINT SELECTED EDITS ON DEMAND /////
!
!-----------------------------------------------------------------------

IF ( i_edit /= 0  .and.  iprint /= 0 ) THEN

  nprint_save      = nprint
  nprint           = iunit
  prnttest         = .false.
  WRITE (outfile,'(a15,i4.4,a1,i4.4,a1,i5.5,a2)') '/Models_n/model', &
&  j_ray,'_',k_ray,'_',nmodel(ij_ray, ik_ray),'.d'
  outfile          = TRIM(data_path)//TRIM(outfile)
  OPEN (UNIT=iunit,FILE=TRIM(outfile),STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=iunit,FILE=TRIM(outfile),STATUS='old',POSITION='append')

  SELECT CASE (i_edit)
 
    CASE (1) 

      CALL editc (jr_min,jr_max,ij_ray, ik_ray)
      nprint       = nprint_save
      RETURN

    CASE (2)

      CALL editp ( jr_min, jr_max )
      CALL edit_e( jr_min, jr_max, ij_ray, ik_ray)
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray)
      CALL editma( jr_min, jr_max, ij_ray, ik_ray)
      CALL edith ( jr_min, jr_max, ij_ray, ik_ray)
      CALL editps( jr_min, jr_max, ij_ray, ik_ray)
      CALL editsc( jr_min, jr_max, ij_ray, ik_ray)
      CALL editu ( jr_min, jr_max, ij_ray, ik_ray)
      CALL edity ( jr_min, jr_max, ij_ray, ik_ray)
      nprint       = nprint_save
      RETURN

    CASE (3)

      CALL editp ( jr_min, jr_max )
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edit_e( jr_min, jr_max, ij_ray, ik_ray )
      CALL editma( jr_min, jr_max, ij_ray, ik_ray )
      CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editps( jr_min, jr_max, ij_ray, ik_ray )
      CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
      CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray )
      nprint       = nprint_save
      RETURN

    CASE (4)

      CALL editp ( jr_min, jr_max )
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edit_e( jr_min, jr_max, ij_ray, ik_ray )
      CALL editma( jr_min, jr_max, ij_ray, ik_ray )
      CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editps( jr_min, jr_max, ij_ray, ik_ray )
      CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
      CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray )
      CALL editng( n1, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray )
      CALL editng( n2, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray )
      CALL editng( n3, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray )
      CALL editng( n4, jnumin, jnumax, ij_ray, ik_ray )
      nprint       = nprint_save
      RETURN

    CASE (5:)

      WRITE (nprint,2001)
      WRITE (nlog,2001)
      RETURN

  END SELECT

  CLOSE (unit=iunit,status='keep')

END IF ! i_edit /= 0

!-----------------------------------------------------------------------
!
!                     \\\\\ EDIT BY CRITERIA /////
!                     /////                  \\\\\  
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Possibly modIfy edit parameters
!-----------------------------------------------------------------------

CALL editpmtr

!-----------------------------------------------------------------------
!  Increment edit counters
!-----------------------------------------------------------------------

nprt(ij_ray, ik_ray) = nprt(ij_ray, ik_ray)  + 1
npltf              = npltf + 1

!-----------------------------------------------------------------------
!  Determine central density extrema for subroutine editc
!-----------------------------------------------------------------------

CALL roextr

!-----------------------------------------------------------------------
!
!             \\\\\ PRINT AT SPECIFIED CYCLE NUMBERS /////
!
!-----------------------------------------------------------------------

IF ( iprint /= 0 ) THEN

  WRITE (outfile,'(a15,i4.4,a1,i4.4,a1,i5.5,a2)') '/Models_n/model', &
&  j_ray,'_',k_ray,'_',nmodel(ij_ray, ik_ray),'.d'
  outfile          = TRIM(data_path)//TRIM(outfile)
  OPEN (UNIT=iunit,FILE=TRIM(outfile),STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=iunit,FILE=TRIM(outfile),STATUS='old',POSITION='append')

  prnttest         = .true.
  nprint_save      = nprint
  nprint           = iunit
  CALL editp ( jr_min, jr_max )
  CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
  CALL edit_e( jr_min, jr_max, ij_ray, ik_ray )
  CALL editmi( jr_min, jr_max )
  CALL editma( jr_min, jr_max, ij_ray, ik_ray )
  CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
  CALL editps( jr_min, jr_max, ij_ray, ik_ray )
  CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
  CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
  CALL editsc( jr_min, jr_max, ij_ray, ik_ray) 
  CALL editn ( n1,jnumin, jnumax, ij_ray, ik_ray )
  CALL editng( n1,jnumin, jnumax, ij_ray, ik_ray )
  CALL editn ( n2,jnumin, jnumax, ij_ray, ik_ray )
  CALL editng( n2,jnumin, jnumax, ij_ray, ik_ray)
  CALL editn ( n3,jnumin, jnumax, ij_ray, ik_ray )
  CALL editng( n3,jnumin, jnumax, ij_ray, ik_ray )
  CALL editn ( n4,jnumin, jnumax, ij_ray, ik_ray )
  CALL editng( n4,jnumin, jnumax, ij_ray, ik_ray )
  nprint       = nprint_save

  CLOSE (unit=iunit,status='keep')

END IF ! iprint /= 0

!-----------------------------------------------------------------------
!  Increment Models_n file counter
!-----------------------------------------------------------------------

IF ( iprint /= 0 ) THEN
  IF ( nprt(ij_ray, ik_ray) >= intprt ) THEN
    nprt(ij_ray, ik_ray)    = 0
    nmodel(ij_ray, ik_ray)  = nmodel(ij_ray, ik_ray) + 1
  END IF ! nprt ge intprt
END IF ! iprint /= 0

!-----------------------------------------------------------------------
!
!             \\\\\ PRINT AT SPECIFIED TIMES /////
!
!-----------------------------------------------------------------------

IF ( iprint /= 0 ) THEN
  IF ( time >= tprint(itprint)  .and.  itprint < itprintmx ) THEN

    WRITE (outfile_t,'(a16,i4.4,a1,i4.4,a1,i5.5,a2)') '/Models_t/modelt', &
& j_ray,'_',k_ray,'_',i,'.d'
    outfile_t      = TRIM(data_path)//TRIM(outfile_t)
    OPEN (UNIT=iunit,FILE=TRIM(outfile_t),STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=iunit,FILE=TRIM(outfile_t),STATUS='old',POSITION='append')

    itprint        = itprint + 1
    prnttest       = .false.
    nprint_save    = nprint
    nprint         = iunit
    CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
    CALL edit_e( jr_min, jr_max, ij_ray, ik_ray )
    CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
    CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
    nprint         = nprint_save

    CLOSE (UNIT=iunit,STATUS='keep')

  END IF ! time > tprint(itprint)
END IF ! iprint /= 0

!-----------------------------------------------------------------------
!
!           \\\\\ PRINT AT SPECIFIED CENTRAL DENSITIES /////
!
!-----------------------------------------------------------------------

IF ( iprint /= 0 ) THEN

  DO i = 1,10
    IF ( rhor(2) < rhoprint(i)  .and.  rho (2) >= rhoprint(i) ) THEN
      prnttest     = .false.
      nprint_save  = nprint
      nprint       = iunitd
      WRITE (outfile_d,'(a16,i4.4,a1,i4.4,a1,i5.5,a2)') '/Models_d/modeld', &
&      j_ray,'_',k_ray,'_',i,'.d'
      outfile_d    = TRIM(data_path)//TRIM(outfile_d)
      OPEN (UNIT=iunitd,FILE=TRIM(outfile_d),STATUS='new',IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=iunitd,FILE=TRIM(outfile_d),STATUS='old')

      CALL editp ( jr_min, jr_max )
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edit_e( jr_min, jr_max, ij_ray, ik_ray )
      CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
      CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray )

      nprint       = nprint_save
      CLOSE (unit=iunitd,status='keep')
    END IF ! rhor(2) < rhoprint(i)  and  rho (2) > rhoprint(i)
  END DO

END IF ! iprint /= 0

!-----------------------------------------------------------------------
!
!       \\\\\ PRINT AT SPECIFIED POST-BOUNCE TINE INTERVALS /////
!
!-----------------------------------------------------------------------

IF ( iprint /= 0 ) THEN

  IF ( it_edit /= 0  .AND.  t_bounce /= zero ) THEN

    t_tb           = time - t_bounce
    tmult          = 1.d+3/dt_edit
    itime          = int( t_tb * tmult )
    itimeprev      = int( ( t_tb - dtnph ) * tmult )
    IF ( itime /= itimeprev ) THEN 

      prnttest     = .false.
      nprint_save  = nprint
      nprint       = iunitt
      WRITE (outfile_dt,'(a17,i4.4,a1,i4.4,a1,i5.5,a2)') '/Models_t/modeldt', &
&      j_ray,'_',k_ray,'_',itime,'.d'
      outfile_dt   = TRIM(data_path)//TRIM(outfile_dt)
      OPEN (UNIT=iunitt,FILE=TRIM(outfile_dt),STATUS='new',IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=iunitt,FILE=TRIM(outfile_dt),STATUS='old')

      CALL editp ( jr_min, jr_max )
      CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edit_e( jr_min, jr_max, ij_ray, ik_ray )
      CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
      CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
      CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
      CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray )
      CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray )

      nprint       = nprint_save
      CLOSE (unit=iunitt,status='keep')
    END IF ! itime /= itimeprev

  END IF ! it_edit /= 0  .AND  t_bounce /= zero

END IF ! iprint /= 0

!-----------------------------------------------------------------------
!
!             \\\\\ UPDATE PLOT FILES /////
!
!-----------------------------------------------------------------------

CALL lumplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, kmin, &
& kmax, ny, nz, nnu )
CALL enuvplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, kmin, &
& kmax, ny, nz, nnu )
CALL nuradplot( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, jmax, kmin, &
& kmax, ny, nz, nez, nnu )
CALL nuplot( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, jmin, &
& jmax, kmin, kmax, ny, nz, nez, nnu )
CALL shock_plot( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& jmin, jmax, kmin, kmax, nx, ny, nz, nnu )

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0  .and.  ij_ray == 1  .and.  ik_ray == 1 ) THEN

  CALL rnuplot
  CALL bnuplot( ij_ray, ik_ray )
  CALL varplot( ij_ray, ik_ray )
  CALL complot( ij_ray, ik_ray )
  CALL lagrangeplot(ij_ray, ik_ray)
  CALL rlagplot

!-----------------------------------------------------------------------
!               ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0  .and.  ij_ray == 1  .and.  ik_ray == 1

RETURN
END SUBROUTINE edit_exec
