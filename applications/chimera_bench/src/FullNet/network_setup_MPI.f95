SUBROUTINE network_setup
!=======================================================================
! This routine handles nuclear data I/O for MPI versions, broadcasting 
! the necessary nuclear and reaction data from PE0 to the production PEs
!=======================================================================


USE array_module, ONLY : npe=>n_proc
USE parallel_module, ONLY : myid, ierr

USE abundances, ONLY : yo, y, yt, ydot, dy
USE cross_sect_data, ONLY : nreac, csect1, rc1, q1, rpf1, n1i, iwk1, ires1, &
& irev1, iffn, csect2, rc2, q2, rpf2, h2, n2i, iwk2, ires2, irev2, csect3, &
& rc3, q3, n3i, iwk3, ires3, irev3, h3
USE controls, ONLY : kstmx, knrmx, idiag, itso, iweak, iscrn, iconvc, &
& changemx, ytime, tolm, tolc, ymin, tdelmm
USE edit_module, ONLY : nlog
USE ffn_data, ONLY : nffn, rf, r1, r2, ffnsum, ffnenu
USE nuc_number, ONLY : ny
USE nuclear_data, ONLY : nname, zz, nn,  aa, be
USE part_funct_data, ONLY : g, gg, angm, t9i
USE reac_rate_data, ONLY : la, le, nan, mu1, a1, b1, n11, mu2, a2, b2, &
& n21, n22, mu3, a3, b3, n31, n32, n33

USE mpi

IMPLICIT none
SAVE

COMMON /ceign/ nds, ndi,idel, jdel, ijd

!IMPLICIT none
!SAVE

CHARACTER (LEN=120)              :: data_dir
CHARACTER (LEN=80)               :: data_desc

INTEGER                          :: i
INTEGER                          :: j
INTEGER                          :: nbc
INTEGER                          :: nr, idel, jdel, ijd, ndi, nds

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Number species=',i5)
  103 FORMAT (' i=',i4,' nname=',a5,' zz=',i4,' aa=',i4,' be=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

data_dir             = '../../FullNet/Data_alpha'

!-----------------------------------------------------------------------
!  Parameters which control the action of the network.  
!-----------------------------------------------------------------------

kstmx             = 2000        ! max # of timesteps for each zone
knrmx             = 5           ! max # of Newton-Raphson iterations before retry
idiag             = 0           ! sets diagnostic output level
itso              = 0           ! sets per timestep output level
iweak             = 1           ! controls the treatment of weak reactions 
iscrn             = 0           ! controls the treatment of nuclear screening
iconvc            = 1           ! determines which convergence condition is used
changemx          = 1.00d-01    ! allowed abundance change used to set the timestep.
ytime             = 1.00d-07    ! abundances > ytime used for timestep calculation
tolm              = 1.00d-06    ! mass conservation convergence criterion
tolc              = 1.00d-04    ! convergence limit on the iterative abundance change
ymin              = 1.00d-25    ! abundance < ymin is set to 0.0
tdelmm            = 2.00d+00    ! max factor by which the timestep is changed

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Read the nuclear and reaction data
!-----------------------------------------------------------------------

  CALL read_nuclear_data( data_dir, data_desc )
  CALL read_reaction_data( data_dir )

!-----------------------------------------------------------------------
!  Write descriptions of nuclei to MC header
!-----------------------------------------------------------------------

  WRITE (nlog,"(a)") trim(data_desc)
  WRITE (nlog,101) ny
  WRITE (nlog,103) (i,nname(i),int(zz(i)),int(aa(i)),8.07144*nn(i)+7.28899*zz(i)-be(i),i=1,ny)

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid

!-----------------------------------------------------------------------
!  Share data for nuc_number module
!-----------------------------------------------------------------------

CALL mpi_bcast(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


!-----------------------------------------------------------------------
!  Share data for nuclear_data module
!-----------------------------------------------------------------------

If ( myid /= 0 ) ALLOCATE(nname(0:ny),aa(ny),zz(ny),nn(ny),be(ny))

CALL MPI_BCAST(aa,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(zz,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nn,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(be,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = 5 * ( ny + 1) 

!    nbc=ny

CALL MPI_BCAST(nname,nbc,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------------------
!  Share data for the part_funct_data module
!-----------------------------------------------------------------------

IF (myid/=0) ALLOCATE(g(24,ny),gg(0:ny),angm(0:ny))

CALL MPI_BCAST(t9i,24,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = 24 * ny
CALL MPI_BCAST(g,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = ny + 1
CALL MPI_BCAST(angm,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------------------
!  Share data for the cross_sect_data module
!-----------------------------------------------------------------------

CALL MPI_BCAST(nreac,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

IF ( myid /= 0 ) THEN
  nr              = nreac(1)
  ALLOCATE (csect1(nr),rc1(7,nr),q1(nr),rpf1(nr))
  ALLOCATE (n1i(3,nr),iwk1(nr),ires1(nr),irev1(nr),iffn(nr))
END IF ! myid /= 0

nbc               = 7 * nreac(1)
CALL MPI_BCAST(rc1,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = 3 * nreac(1)
CALL MPI_BCAST(n1i,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
nbc               = nreac(1)
CALL MPI_BCAST(iwk1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
CALL MPI_BCAST(ires1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
CALL MPI_BCAST(irev1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 

IF ( myid /= 0 ) THEN
  nr              = nreac(2)
  ALLOCATE (csect2(nr),rc2(7,nr),q2(nr),rpf2(nr),h2(nr))
  ALLOCATE (n2i(4,nr),iwk2(nr),ires2(nr),irev2(nr))
END IF ! myid /= 0

nbc               = 7 * nreac(2)
CALL MPI_BCAST(rc2,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = 4 * nreac(2)
CALL MPI_BCAST(n2i,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
nbc               = nreac(2)
CALL MPI_BCAST(iwk2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
CALL MPI_BCAST(ires2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
CALL MPI_BCAST(irev2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 

IF ( myid /= 0 ) THEN
  nr              = nreac(3)
  ALLOCATE (csect3(nr),rc3(7,nr),q3(nr),h3(nr))
  ALLOCATE (n3i(3,nr),iwk3(nr),ires3(nr),irev3(nr))
END IF ! myid /= 0

nbc               = 7 * nreac(3)
CALL MPI_BCAST(rc3,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = 3 * nreac(3)
CALL MPI_BCAST(n3i,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
nbc               = nreac(3)
CALL MPI_BCAST(iwk3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
CALL MPI_BCAST(ires3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
CALL MPI_BCAST(irev3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 

!-----------------------------------------------------------------------
!  Share the data for the reac_rate_data module
!-----------------------------------------------------------------------

CALL MPI_BCAST(nan,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

IF (myid/=0) THEN
  ALLOCATE (la(3,ny),le(3,ny))
  ALLOCATE (mu1(nan(1)),a1(nan(1)),b1(nan(1)))
  ALLOCATE (mu2(nan(2)),a2(nan(2)),b2(nan(2)))
  ALLOCATE (mu3(nan(3)),a3(nan(3)),b3(nan(3)))
END IF ! myid /= 0

nbc               = 3 * ny
CALL MPI_BCAST(la,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(le,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
nbc               = nan(1)
CALL MPI_BCAST(mu1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(a1,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = nan(2)
CALL MPI_BCAST(mu2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(a2,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
nbc               = nan(3)
CALL MPI_BCAST(mu3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(a3,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

IF (myid/=0) THEN
  ALLOCATE (n11(nan(1)),n21(nan(2)),n22(nan(2)))
  ALLOCATE (n31(nan(3)),n32(nan(3)),n33(nan(3)))
  DO j = 1,nan(1)
    n11(j)        = n1i(1,mu1(j))
  END DO
  DO j = 1,nan(2)
    n21(j)        = n2i(1,mu2(j))
    n22(j)        = n2i(2,mu2(j))
  END DO
  DO j = 1,nan(3)
    n31(j)        = n3i(1,mu3(j))
    n32(j)        = n3i(2,mu3(j))
    n33(j)        = n3i(3,mu3(j))
  END DO
END IF ! myid /= 0

!-----------------------------------------------------------------------
!  Share the matrix shape parameters
!-----------------------------------------------------------------------

nbc               = 5
CALL MPI_BCAST(nds,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------------------
!  Share the data for the ffn_data
!-----------------------------------------------------------------------

CALL MPI_BCAST(nffn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

IF ( nffn > 0 ) THEN

  IF ( myid /= 0 ) THEN
    ALLOCATE (rf(nffn),r1(nffn),r2(nffn))
    ALLOCATE (ffnsum(nffn,143),ffnenu(nffn,143))
  END IF ! myid /= 0

  CALL MPI_BCAST(rf,nffn,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(r1,nffn,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(r2,nffn,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc             = nffn * 143
  CALL MPI_BCAST(ffnsum,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ffnenu,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

END IF ! nffn > 0

!-----------------------------------------------------------------------
!  Set sizes of abundance arrays
!-----------------------------------------------------------------------

ALLOCATE ( yo(ny), y(ny), yt(ny), ydot(ny), dy(ny) )

RETURN
END SUBROUTINE network_setup
