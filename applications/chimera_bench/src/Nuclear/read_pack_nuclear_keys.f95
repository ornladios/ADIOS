SUBROUTINE read_pack_nuclear_keys( nreadp, nprint, iskip, nx, ij_ray_dim, &
& ik_ray_dim, nnc, c_nuc_data, i_nuc_data, d_nuc_data, nrst, nuc_number )
!-----------------------------------------------------------------------
!
!    File:         read_pack_nuclear_keys
!    Module:       read_pack_nuclear_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the nuclear keys defining the mode and tolerances of nuclear
!       burning, and to pack them into an integer and a real*8 array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp     : unit number from which to read
!  nprint     : unit number from which to print
!  iskip      : echo data read flag
!  nx         : radial array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  nnc        : nuclear specie array dimension
!  nrst       : cycle number at start or restart
!
!    Output arguments:
!  nuc_number : number of nuclear species (not counting representative heavy nucleus)
!  c_nuc_data : character array of nuclei
!  i_nuc_data : integer array of nuclear keys
!  d_nuc_data : real*8 array of nuclear keys
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, nucbrn_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon
USE physcnst_module, ONLY : mn_ex, mp_ex

USE edit_module, ONLY : nlog
USE nucbrn_module, ONLY : nse

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nnc           ! nuclear specie array dimension
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER (len=5), INTENT(out), DIMENSION(nnc)                 :: c_nuc_data  ! character array of nuclei

INTEGER, INTENT(out)                                           :: nuc_number  ! number of nuclear species (not counting representative heavy nucleus)
INTEGER, INTENT(out), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim)                         :: i_nuc_data  ! integer array of edit keys

REAL(KIND=double), INTENT(out), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 6 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line
CHARACTER (len=10)               :: var_name
CHARACTER (len=5)                :: c_nuc         ! name of nucleus

INTEGER                          :: n             ! integer variable read in a integer datum
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl1           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl2           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl3           ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!  Nuclear network controls
!-----------------------------------------------------------------------
!  inuc : nuclear reaction network switch.
!
!     inuc = 0 : nuclear reactions (for nse(j) = 0 zones) bypassed.
!     inuc = 1 : nuclear reactions (for nse(j) = 0 zones) included.
!     inuc = 2 : no network, no flashing
!-----------------------------------------------------------------------

INTEGER                                         :: inuc

!-----------------------------------------------------------------------
!  Nuclear network adjust switch
!-----------------------------------------------------------------------
!  iadjst : nuclear reaction network switch.
!
!     iadjst = 0 : leave nuclear mass fractions unnormalized after read in
!     iadjst = 1 : normalize nuclear mass fractions after read in
!-----------------------------------------------------------------------

INTEGER                                         :: iadjst

!-----------------------------------------------------------------------
!  Nuclear reaction rate tolerances
!-----------------------------------------------------------------------
!  itnuc   : maximum number of iterations to attempt in order to obtain a convergent solution of
!   the reaction rate equations for the ion abundances (i.e., the xn(j,i)'s) and the temperature
!   in zones not assumed to be in nse. If itnuc = 1, the variables are assumed to have converged
!   after the first iteration attempt.
!
!  ttolnuc : temperature convergence parameter for nuclear reaction rate equations. The criterion
!   for temperature convergence is that
!
!     abs(dt/t) < ttolnuc .
!
!  ytolnuc : ion abundance convergence parameter for nuclear reaction rate equations. The
!   criteria for ion abundance convergence is that
!
!     abs(dy)/(y(i) + ynmin) < ytolnuc .
!-----------------------------------------------------------------------

INTEGER                                         :: itnuc

REAL(KIND=double)                               :: ttolnuc
REAL(KIND=double)                               :: ytolnuc
REAL(KIND=double)                               :: ynmin

!-----------------------------------------------------------------------
!  Time step controls
!-----------------------------------------------------------------------
!  t_cntl_burn(1) : nuclear burn temperature change time step criterion, i.e., the maximum
!   permitted abs( dT_burn(j)/t(j) ), where dT_burn(j) is the nuclear burn temperature
!   change in radial zone j.
!
!  t_cntl_burn(2) : nuclear burn composition change time step criterion, i.e., the maximum
!   permitted abs( dyn(j,i)/yn(j,i) ), where dyn(j,i) is the abundance change of specie i in
!   radial zone j.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(2)                 :: t_cntl_burn

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------
!  a_name      : Mass and symbol of nucleus.
!
!  a_nuc(n)    : mass number of the nth nuclear species.
!
!  z_nuc(n)    : charge number of the nth nuclear species.
!
!  m_ex_nuc(n) : mass excess of the nth nuclear species (MeV).
!
!  be_nuc(n)   : binding energy of the nth nuclear species (MeV).
!-----------------------------------------------------------------------

CHARACTER (len=5), ALLOCATABLE, DIMENSION(:)    :: a_name

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: a_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: z_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: m_ex_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: be_nuc

  101 FORMAT (a128)
  105 FORMAT (1x,a128)
  121 FORMAT (10x,2i10)
  123 FORMAT (1x,a6,4x,2i10,42x,a16)
  131 FORMAT (20x,i10,5x,e15.8)
  133 FORMAT (1x,a6,14x,i10,5x,es15.8,22x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  151 FORMAT (10x,i10)
  153 FORMAT (1x,a6,4x,i10,52x,a16)
  221 FORMAT (6x,a5,i9,3(5x,e15.8))
  223 FORMAT (1x,a6,a5,i9,3(5x,es15.8))
  231 FORMAT (10x,i10,3(5x,e15.8))
  233 FORMAT (1x,a6,4x,i10,3(5x,es15.8))
  341 FORMAT (10x,i10,5(i3,1x,e15.8))
  343 FORMAT (1x,a6,4x,i10,5(i3,1x,es15.8))
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)
 1001 FORMAT (' Allocation problem for array ',a10,' in read_pack_nuclear_keys')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_pack_nuclear_keys')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (a_name(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_name    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (a_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (m_ex_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm_ex_nuc  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                   \\\\\ READ NUCLEAR KEYS /////
!
!-----------------------------------------------------------------------

a_nuc                = zero
z_nuc                = zero
m_ex_nuc             = zero
be_nuc               = zero

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nuc_number           = 0
REWIND (nreadp)

READ: DO

!-----------------------------------------------------------------------
!  Read line
!-----------------------------------------------------------------------

  READ (nreadp,101,END=5000) line

  type               = line(1:6)
  name               = line(73:88)

!-----------------------------------------------------------------------
!  Ignore comment cards in the data stream
!-----------------------------------------------------------------------

  IF ( type == 'cccccc'  .or.  type(1:1) == '!'  .or.  type == '      ') THEN
    IF ( iskip == 0 ) WRITE (nprint,105) line
    CYCLE
  END IF ! type = 'cccccc' or type = '!'

!-----------------------------------------------------------------------
!
!           \\\\\ NUCLEAR NETWORK CONTROL PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  inuc
!-----------------------------------------------------------------------

  IF ( type == 'inuc  ' ) THEN
     READ (line ,121) inuc,itnuc  
     IF ( iskip .eq. 0 ) write (nprint,123) type,inuc,itnuc,name
     CYCLE
   END IF ! type = 'inuc  '

!-----------------------------------------------------------------------
!  tolnuc
!-----------------------------------------------------------------------

  IF ( type == 'tolnuc' ) THEN
    READ (line ,141) rl

    IF ( name == 'ttolnuc ' ) THEN
      ttolnuc        = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,ttolnuc,name
      CYCLE
    END IF ! name = 'ttolnuc '

    IF ( name == 'ytolnuc ' ) THEN
      ytolnuc        = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,ytolnuc,name
      CYCLE
    END IF ! name = 'ytolnuc '

    IF ( name == 'ynmin  ' ) THEN
      ynmin          = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,ynmin,name
      CYCLE
    END IF ! name = 'ynmin   '

  END IF ! type = 'tolnuc'

!-----------------------------------------------------------------------
!
!          \\\\\ NUCLEAR TIME STEP CONTROL PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  tcntrl
!-----------------------------------------------------------------------

  IF ( type == 'tcntrl' ) THEN
    READ (line ,131) n,rl

    IF ( name == 't_cntl_burn') THEN
      t_cntl_burn(n) = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,t_cntl_burn(n),name
      CYCLE
    END IF ! name = 't_cntl_burn'

  END IF ! name = 'tcntrl'

!-----------------------------------------------------------------------
!
!                   \\\\\ COMPOSITION DATA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  iadjst
!-----------------------------------------------------------------------

  IF ( type == 'iadjst' ) THEN
    READ (line ,151) iadjst
    IF ( iskip == 0 ) WRITE (nprint,153) type,iadjst
     CYCLE
   END IF ! type = 'iadjst '

!-----------------------------------------------------------------------
!  a_nuc
!-----------------------------------------------------------------------

  IF ( type == 'a_nuc ' ) THEN
    READ (line ,221) c_nuc,n,rl1,rl2,rl3
    a_name(n)        = c_nuc
    a_nuc(n)         = rl1
    z_nuc(n)         = rl2
    m_ex_nuc(n)      = rl3
    nuc_number       = MAX( nuc_number, n )
    IF ( iskip == 0 ) WRITE (nprint,223) type,a_name(n),n,a_nuc(n),z_nuc(n),m_ex_nuc(n)
     CYCLE
   END IF ! type = 'a_nuc '

!-----------------------------------------------------------------------
!  Unrecognized card
!-----------------------------------------------------------------------

  IF ( nrst == 0 ) THEN
    WRITE (nprint,401)
    WRITE (nprint,403) line
  END IF

END DO READ

!-----------------------------------------------------------------------
!
!                \\\\\ CALCULATE BINDING ENERGIES /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE

DO n = 1,nuc_number
  be_nuc(n)                  = z_nuc(n) * mp_ex + ( a_nuc(n) - z_nuc(n) ) * mn_ex - m_ex_nuc(n)
END DO

!-----------------------------------------------------------------------
!
!                 \\\\\ NAME AUXILIARY NUCLEUS /////
!
!-----------------------------------------------------------------------
 
a_name(nuc_number+1)         = ' Aux '

!-----------------------------------------------------------------------
!
!                   \\\\\ PACK NUCLEAR KEYS /////
!
!-----------------------------------------------------------------------

DO n = 1,nnc
  c_nuc_data(n)              = a_name(n)
END DO

i_nuc_data(1,:,:)            = inuc
i_nuc_data(2,:,:)            = itnuc
i_nuc_data(3,:,:)            = nuc_number
i_nuc_data(4,:,:)            = iadjst

d_nuc_data(1,:,:)            = ttolnuc
d_nuc_data(2,:,:)            = ytolnuc
d_nuc_data(3,:,:)            = ynmin
d_nuc_data(4,:,:)            = t_cntl_burn(1)
d_nuc_data(5,:,:)            = t_cntl_burn(2)

DO n = 1,nnc
  d_nuc_data(10+n,:,:)       = a_nuc(n)
  d_nuc_data(10+nnc+n,:,:)   = z_nuc(n)
  d_nuc_data(10+2*nnc+n,:,:) = m_ex_nuc(n)
  d_nuc_data(10+3*nnc+n,:,:) = be_nuc(n)
END DO ! n

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (a_name, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_name    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (a_nuc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (z_nuc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (m_ex_nuc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm_ex_nuc  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (be_nuc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc    '; WRITE (nlog,2001) var_name; END IF


RETURN
END SUBROUTINE read_pack_nuclear_keys
