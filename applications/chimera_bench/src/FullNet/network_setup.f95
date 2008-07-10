SUBROUTINE network_setup
!----------------------------------------------------------------------
!  This routine loads the control flags and the network data files.
!----------------------------------------------------------------------

USE abundances, ONLY : yo, y, yt, ydot, dy
USE controls, ONLY : kstmx, knrmx, idiag, itso, iweak, iscrn, iconvc, &
& changemx, ytime, tolm, tolc, ymin, tdelmm
USE edit_module, ONLY : nlog
USE nuc_number, ONLY : ny

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (LEN=120)              :: data_dir
CHARACTER (LEN=80)               :: data_desc

INTEGER                          :: nmf

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Beginning to read nuclear data')
  103 FORMAT (' Nuclear data has been read')
  105 FORMAT (' Beginning to read reaction rate data')
  107 FORMAT (' Reaction rate data has been read')

!      use evh1_zone
!      integer :: i,j,k,n,nz,izone ! Loop indices
!      character (LEN=80) :: descript(3),data_desc,data_dir,diag_file

!-----------------------------------------------------------------------
!  Read the parameters which control the action of the network.  
!-----------------------------------------------------------------------
!      call input_search('NETWORK input')
!      Read(1,"(72x)")
!      Read(1,"(a80)") (descript(i),i=1,3) ! text description of the problem.
!      Read(1,"(72x)")
!      Read(1,*) szone    ! number of the zone with which to begin
!      Read(1,*) nzone    ! total # of zones
!      Read(1,*) kstmx    ! max # of timesteps for each zone
kstmx                = 2000        ! max # of timesteps for each zone
!      Read(1,*) knrmx             ! max # of Newton-Raphson iterations before retry
knrmx                = 5           ! max # of Newton-Raphson iterations before retry
!      Read(1,*) idiag             ! sets diagnostic output level
idiag                = 0           ! sets diagnostic output level
!      Read(1,*) itso              ! sets per timestep output level
itso                 = 0           ! sets per timestep output level
!      Read(1,*) iweak             ! controls the treatment of weak reactions 
iweak                = 1           ! controls the treatment of weak reactions 
!      Read(1,*) iscrn             ! controls the treatment of nuclear screening
iscrn                = 0           ! controls the treatment of nuclear screening
!      Read(1,*) iconvc            ! determines which convergence condition is used
iconvc               = 1           ! determines which convergence condition is used
!      Read(1,*) changemx          ! allowed abundance change used to set the timestep.
changemx             = 1.00d-01    ! allowed abundance change used to set the timestep.
!      Read(1,*) ytime             ! abundances > ytime used for timestep calculation
ytime                = 1.00d-07    ! abundances > ytime used for timestep calculation
!      Read(1,*) tolm              ! mass conservation convergence criterion
tolm                 = 1.00d-06    ! mass conservation convergence criterion
!      Read(1,*) tolc              ! convergence limit on the iterative abundance change
tolc                 = 1.00d-04    ! convergence limit on the iterative abundance change
!      Read(1,*) ymin              ! abundance < ymin is set to 0.0
ymin                 = 1.00d-25    ! abundance < ymin is set to 0.0
!      Read(1,*) tdelmm            ! max factor by which the timestep is changed
tdelmm               = 2.00d+00    ! max factor by which the timestep is changed
!      Read(1,"(72x)")
!      Read(1,"(a80)") data_dir
!      Read(1,"(72x)")
!      Read(1,"(a80)") diag_file 
!      Close(1)

ny                   = 299
nmf                  = 14

!-----------------------------------------------------------------------
! Read nuclear and reaction data
!-----------------------------------------------------------------------

data_dir             = '../../FullNet/Data_alpha'

WRITE (nlog,101)
CALL read_nuclear_data( data_dir, data_desc )
WRITE (nlog,103)
WRITE (nlog,105)
CALL read_reaction_data( data_dir )
WRITE (nlog,107)
!     call read_match_data(data_dir)
!     call flux_init

!-----------------------------------------------------------------------
!  Open diagnositic output file
!-----------------------------------------------------------------------

Open(50,file='diag_file')
!      Write(50,"(a)") (descript(i),i=1,3),data_desc

!-----------------------------------------------------------------------
!  Test network size
!-----------------------------------------------------------------------

IF ( ny /= nmf ) THEN
  WRITE(6,*) 'network size mismatch, ',ny,'/=',nmf
  STOP
END IF

!-----------------------------------------------------------------------
!  Set sizes of abundance arrays
!-----------------------------------------------------------------------

ALLOCATE ( yo(ny), y(ny), yt(ny), ydot(ny), dy(ny) )

RETURN
END
