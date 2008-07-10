SUBROUTINE read_nuclear_abn_to_node( nreadp, nprint, iskip, nx, nnc, nrst, &
& nuc_number, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         read_nuclear_abn_to_node
!    Module:       read_nuclear_abn_to_node
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the nuclear abundance data, and to pack them into an integer
!       and a real*8 array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp     : unit number from which to read
!  nprint     : unit number from which to print
!  iskip      : echo data read flag
!  nx         : radial array dimension
!  nnc        : nuclear specie array dimension
!  nrst       : cycle number at start or restart
!  nuc_number : number of nuclear species (not counting representative heavy nucleus)
!  ij_ray     : index denoting the j-index of a specific radial ray
!  ik_ray     : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_x_module, nucbrn_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon

USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY : nse_e=>nse
USE nucbrn_module, ONLY : nse_n=>nse
USE radial_ray_module, ONLY : xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, &
& uburn_c, nse_c

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: nnc           ! nuclear specie array dimension
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart
INTEGER, INTENT(in)              :: nuc_number    ! number of nuclear species (not counting representative heavy nucleus)
INTEGER, INTENT(in)              :: ij_ray      ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray      ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 6 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line
CHARACTER (len=10)               :: var_name

INTEGER                          :: i             ! radial array index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: n             ! integer variable read in a integer datum
INTEGER                          :: istat         ! allocation status
INTEGER                          :: iadjst        ! 0: leave nuclear abundaces as is
                                                  ! 1: renormalize nuclear abundaces so mass fractions sum to unity

INTEGER                          :: i1            ! integer data variable to read in an interger datum
INTEGER                          :: i2            ! integer data variable to read in an interger datum
INTEGER                          :: i3            ! integer data variable to read in an interger datum
INTEGER                          :: i4            ! integer data variable to read in an interger datum
INTEGER                          :: i5            ! integer data variable to read in an interger datum

REAL(KIND=double)                :: xn_tot        ! sum of themass fractions in a zone

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl1           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl2           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl3           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl4           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl5           ! 64 bit real variable read in a real datum

!.........Abundance parameters..........................................

!      nse(j) : nuclear statistical equilibrium flag.

!      xn(j,i) : mass fraction of the ith nucleus.

!      uburn(j) : cumulative energy generated in zone j by nuclear reactions (ergs/gm).

!      be_nuc_rep(j) : binding energy of the representative heavy nucleus (MeV).

!      a_nuc_rep(j) : mass number of the representative heavy nucleus.

!      z_nuc_rep(j) : charge number of the representative heavy nucleus.

INTEGER, ALLOCATABLE, DIMENSION(:)              :: nse

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)  :: xn

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: uburn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: be_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: a_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)    :: z_nuc_rep

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
 1001 FORMAT (' Allocation problem for array ',a10,' in read_nuclear_abn_to_node')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_nuclear_abn_to_node')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (nse(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (xn(nx,0:nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uburn(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc_rep(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_rep(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_rep(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE NUCLEAR ABUNDANCE DATA /////
!
!-----------------------------------------------------------------------

xn                   = zero

uburn                = zero
be_nuc_rep           = 492.3d0
a_nuc_rep            = 56.d0
z_nuc_rep            = 28.d0

nse                  = 1

!-----------------------------------------------------------------------
!
!                \\\\\ READ NUCLEAR ABUNDANCE DATA /////
!
!-----------------------------------------------------------------------

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
!  Composition data
!-----------------------------------------------------------------------

!........iadjst

  IF ( type == 'iadjst' ) THEN
    READ (line ,151) iadjst
    IF ( iskip == 0 ) WRITE (nprint,153) type,iadjst
     CYCLE
   END IF ! type = 'iadjst '

!........nccomp

  IF ( type == 'nccomp' ) THEN
    READ (line ,341) j,i1,rl1,i2,rl2,i3,rl3,i4,rl4,i5,rl5
    xn(j,i1)         = rl1
    xn(j,i2)         = rl2
    xn(j,i3)         = rl3
    xn(j,i4)         = rl4
    xn(j,i5)         = rl5
    nse(j)           = 0
    IF ( iskip == 0 ) WRITE (nprint,343) type,j,i1,xn(j,i1),i2,xn(j,i2),i3,xn(j,i3), &
&    i4,xn(j,i4),i5,xn(j,i5)
    CYCLE
  END IF ! type = 'nccomp'

!........a_nuc_rep

  IF ( type == 'a_nucr' ) THEN
    READ (line ,231) j,rl1,rl2,rl3
    a_nuc_rep(j)     = rl1
    z_nuc_rep(j)     = rl2
    be_nuc_rep(j)    = rl3
    IF ( iskip == 0 ) WRITE (nprint,233) type,j,a_nuc_rep(j),z_nuc_rep(j),be_nuc_rep(j)
     CYCLE
   END IF ! type = 'a_nuc_rep '

!........burn

   IF ( type == 'burn  ' ) THEN
     READ (line,131) j,rl
     uburn(j)        = rl
     IF ( iskip == 0 ) WRITE (nprint,133) type,j,uburn(j) 
     CYCLE
   END IF ! type = 'burn  '

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
!                \\\\\ ADJUST NUCLEAR ABUNDANCES /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE

IF ( iadjst == 1 ) THEN
  DO j = 1,nx
    IF ( nse(j) == 0 ) THEN

      xn_tot                        = zero
      DO n = 1,nuc_number+1
        xn_tot                      = xn_tot + xn(j,n)
      END DO
  
      DO n = 1,nuc_number+1
        xn(j,n)                     = xn(j,n)/( xn_tot + epsilon )
      END DO

    END IF ! nse(j) == 0
  END DO ! j
END IF ! iadjst == 1

!-----------------------------------------------------------------------
!
!         \\\\\ TRANSFER NUCLEAR ABUNDANCE DATA TO MODULES /////
!
!-----------------------------------------------------------------------

DO i = 1,nx
  nse_n(i)                          = nse(i)
  nse_e(i,ij_ray,ik_ray)            = nse(i)
END DO  

DO i = 2,nx
  nse_c(i-1,ij_ray,ik_ray)          = nse(i)
END DO  

DO n = 1,nnc
  DO i = 2,nx
    xn_c(i-1,n,ij_ray,ik_ray)       = xn(i,n)
  END DO
END DO

DO i = 2,nx
  a_nuc_rep_c (i-1,ij_ray,ik_ray)   = a_nuc_rep(i)
  z_nuc_rep_c (i-1,ij_ray,ik_ray)   = z_nuc_rep(i)
  be_nuc_rep_c(i-1,ij_ray,ik_ray)   = be_nuc_rep(i)
  uburn_c     (i-1,ij_ray,ik_ray)   = uburn(i)
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (nse, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (xn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (uburn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (be_nuc_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (a_nuc_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (z_nuc_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE read_nuclear_abn_to_node
