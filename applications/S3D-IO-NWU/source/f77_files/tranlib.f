C     CVS $Revision: 1.4 $ created $Date: 2005/09/28 20:58:09 $
C
!  Changes by Mark Fahey, ORNL
!  Update by Evatt Hawkes, SNL
!
      SUBROUTINE MCABS
C
C     Revision 4.6, 1999/06/03 (E. Meeks)
C     1) Fixed problem with KERR not declared LOGICAL in MCORDF
C     Revision 4.5, 1999/04/10 (E. Meeks)
C     1)  Added 4 new subroutines to allow a user to modify
C         data stored in a Transport linking file through
C         library calls.  The new routines are:  MCREWR, which
C         writes a linking file based on data stored in the
C         passed-in work arrays;  MCCDEX, MCCCEX, and
C         MCCVEX which allow extraction of the fitting
C         coefficients for diffusion coefficients, conductivity,
C         and viscosity, respectively, as well as
C         replacement of those values (similar to CKRAEX)
C     Revision 4.4, 1998/12/04 (E. Meeks)
C     1)  Per Action#195:  Added error messages for later parsing
C         into error-handling system, for all fatal errors.
C     2)  Replaced two occurences of 'STOP' statement in MCLMDT,
C         with RETURN; added logical variable in call list KERR
C         that is returned as .TRUE. where the 'STOP's were.
C     3)  Replaced two occurences of 'STOP' statement in MCORDF,
C         with RETURN; added logical variable in call list KERR
C         that is returned as .TRUE. where the 'STOP's were.
C     4)  Added logical KERR to calls to MCORDF and MCLMDT in
C         user-callable subroutines MCMDIF and MCMCDT,respectively.
C     Revision 4.3, 1998/03/03 (E. Meeks)
C     1)  Action #116: Put a change block for LINPACK around
C         the declaration and initialization of DET and JOB
C         variables in routine MCORDF to avoid unused variable warning.
C     Revision: 4.2, Date: 1997/06/10, Author: J. Grcar
C     1)  Fixed bug #035.  Corrected spelling of "corporation".
C     2)  Edited "CVS $Revision" line to fit in 72 columns.
C     3)  Updated output version number and creation date.
C     V. 4.1, 96/05/24
C     1. initial sccs version
C     CHANGES FOR VERSION 4.0 (Initial CHEMKIN-III version,
C                              January 1996 F. Rupley)
C     1.  Add/modify PROLOGUE sections
C     2.  Reverse change comment order in MCABS.
C     3.  binary/ascii linkfile option.
C     4.  separate linkfile, code version numbers.
C
C     CHANGES FOR VERSION 3.9 (2/27/95 F. Rupley)
C     1.  Change character index "(:" to "(1:"
C     CHANGES FOR VERSION 3.8 (1/20/95 F. Rupley)
C     CHANGES FOR VERSION 3.7 (10/3/94) H.K. Moffat
C     1.  Made the small parameter even small for 64 bit systems
C         (change block).
C     CHANGES FOR VERSION 3.6 (8/31/94) H.K. Moffat
C           - Passed check against v3.5 with tranlibtest
C     1.  Performance improvements for MCEVAL and MCEDIF
C     2.  MCEDIF now special cases KK = 1
C     3.  Added a hard-coded Horner's polynomial evaluation routine
C         for NO=4, MCEVAL4
C     4.  Changed MCEVAL so that the MAX function is used to massage
C         the mole fractions instead of the addition function.
C     5.  Changed a statement in MCAVIS to use more exact arithmetic,
C         causing change in 8th digit to occur
C     6.  Performance improvement to MCAVIS - SQRT calls used instead
C         of powers of floating point variables.
C     CHANGES FOR VERSION 3.5 (8/10/94 H.K. Moffat
C     1.  Accepts to 3.2 version number
C     CHANGES FOR VERSION 3.4 (6/12/94) H.K. Moffat
C     1.  Added change blocks for LAPACK linear algebra
C     CHANGES FOR VERSION 3.3 (6/9/94) H.K. Moffat
C     1.  Changed the Gas constant to agree with 1986 CODATA value.
C         (two extra significant digits of accuracy)
C     2.  Changed parameter values for pi and its powers to
C         so that the values are accurate to machine precision.
C     CHANGES FOR VERSION 3.2 (6/8/94 H.K. Moffat
C     1.  Made MCLMDT roughly 30 times faster on a sample problem.
C         Opcount in the function now scales like KK**2, instead
C         of the previous KK**3.
C     2.  No longer need a delta function in the library -
C         CHANGES to MCLMDT and MCORDF
C     CHANGES FOR VERSION 3.1 (4/5/94 M. Coltrin)
C     1.  Fix error in MCLMDT around line 1674 ("ZERO" replaced by
C         "ONE").
C     2.  Indices in formula for BINDIF(I,I) in loop 150 should have
C         been (K,K).
C     3.  A factor of PFAC was omitted from the statement before
C         statement number 1600 in MCLMDT.
C     4.  The factor PIFAC and the statement in which it was used in
C         loop 1600 in MCLMDT were both incorrect.
C     CHANGES FOR VERSION 3.0 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C     CHANGES FOR VERSION 1.7 (10/1/92 F. Rupley per M. Coltrin)
C     1. Created MCABS to hold version and change information
C     2. COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR eliminates
C        the need for argument LINKMC in the MCSAVE call list
C     CHANGES FOR VERSION 1.6
C     1. Versions 1.8 and 1.9 added (TRANLIB V.1.5
C        and TRANFIT V.1.8 were intermediate versions which may
C        not be legitimate; TRANFIT V.1.9 is actually a
C        correction to V.1.7, and TRANLIB 1.6 is an update of
C        V.1.4)
C     CHANGES FOR VERSION 1.4
C     1. Additional record to binary file indicates
C        version, machine precision, and error status
C     2. Additional record to binary file has required lengths for
C        integer, real work arrays
C     3. New Subroutines MCPNT, MCSAVE read, write binary
C        file information, work arrays
C     CHANGES FOR VERSION 1.3
C     1. SUBROUTINE MCLEN
C     CHANGES FROM VERSION 1.1:
C     1. Eliminated many GOTO's
C     CHANGES FROM VERSION 1.0:
C     1. Changed REAL*8 to real
C
C     end of SUBROUTINE MCABS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK,
     1           IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK,
C                     IFLAG)
C  This subroutine reads the transport linkfile from the fitting code
C  and creates the internal storage and work arrays, IMCWRK(*) and
C  RMCWRK(*).  MCINIT must be called before any other transport
C  subroutine is called.  It must be called after the CHEMKIN package
C  is initialized.
C
C  INPUT
C  LINKMC    - Integer scalar, transport linkfile input unit number.
C  LOUT      - Integer scalar, formatted output file unit number.
C  LENIMC    - Integer scalar, minimum dimension of the integer
C              storage and workspace array IMCWRK(*);
C              LENIMC must be at least:
C              LENIMC = 4*KK + NLITE,
C              where KK is the total species count, and
C                    NLITE is the number of species with molecular
C                          weight less than 5.
C  LENRMC    - Integer scalar, minimum dimension of the real storage
C              and workspace array RMCWRK(*);
C              LENRMC must be at least:
C              LENRMC = KK*(19 + 2*NO + NO*NLITE) + (NO+15)*KK**2,
C              where KK is the total species count,
C                    NO is the order of the polynomial fits (NO=4),
C                    NLITE is the number of species with molecular
C                          weight less than 5.
C
C  OUTPUT
C  IMCWRK(*) - Integer workspace array; dimension at least LENIMC.
C  RMCWRK(*) - Real    workspace array; dimension at least LENRMC.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST, NCST, NXL,
     4                NR, NWRK, K3
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
C
      DIMENSION IMCWRK(*), RMCWRK(*)
      CHARACTER*16 PRVERS, VERS, PREC, PRDATE, IFMT, RFMT, CFMT, LFMT
      CHARACTER*80 MSGSTR
      PARAMETER
     1(CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)', RFMT='(1P,5E24.16)')
C
      LOGICAL IOK, ROK, KERR, LBIN
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C     The following number SMALL is used in the mixture diffusion
C     coefficient calculation; its use allows a smooth and well-
C     defined diffusion coefficient as the mixture approaches a pure
C     species, even though stictlyh speaking there does not exist a
C     diffusion coefficient in this case.  The value of SMALL should
C     be small relative to any species mole fraction of importance,
C     but large enough to be represented on the computer.
C
C*****SMALL 1) 64 bit floats
C      SMALL = 1.0E-50
C*****END SMALL 1) 64 bit floats
C*****SMALL 2) 32 bit floats
      SMALL = 1.0E-20
C*****END SMALL 2) 32 bit floats
C
C     Gas constant as reported in 1993 CRC, (J. Research of
C     National Bureau of Standards, 92, 85, 1987).
C     ( 8.314510(70)E+07 Joules mol-1 K-1)
C
      RU    = 8.314510E+07
C
C     Standard atmosphere (defined as an exact quantity)
C
      PATMOS= 1.01325E+06
C
C     Write version number
C
C*****precision > double
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      PREC = 'SINGLE'
C*****END precision > single
C
      PRVERS ='4.6'
      PRDATE ='1999/06/03'
C
      WRITE (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
     1 ' TRANLIB:  CHEMKIN-III MULTICOMPONENT TRANSPORT LIBRARY,',
     2 PREC(1:CKLSCH(PREC)), ' PRECISION Vers. ',
     3 PRVERS(1:CKLSCH(PRVERS)+1), PRDATE,
     4 ' Copyright 1995, Sandia Corporation.',
     5' The U.S. Government retains a limited license in this software.'
C
C     Read the problem size
      CALL MCLEN (LINKMC, LOUT, LI, LR, IFLAG)
      IOK = (LENIMC .GE. LI)
      ROK = (LENRMC .GE. LR)
C
      IF (.NOT.IOK .OR. .NOT.ROK) THEN
         IF (.NOT. IOK) WRITE (LOUT, 300) LI
C <error module="tranlib" severity="error">
C <id>1</id>
C <message>The size of the TRANSPORT integer work array
C (variable IMCWRK) is too small.  You must increase the value of
C LENIMC from %1 to at least %2 in the transport driver and rebuild this
C program.
C </message>
C </error>
         IDERR = 1
         MSGSTR = ' '
         WRITE (MSGSTR, '(I12,A1)') LENIMC,'$', LI
         CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
         IF (.NOT. ROK) WRITE (LOUT, 350) LR
C <error module="tranlib" severity="error">
C <id>2</id>
C <message>The size of the TRANSPORT real work array
C (variable RMCWRK) is too small.  You must increase the value of
C LENRMC from %1 to at least %2 in the transport driver and rebuild this
C program.
C </mesage>
C </error>
         IDERR = 2
         MSGSTR = ' '
         WRITE (MSGSTR, '(I12,A1)') LENRMC,'$', LR
         CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
         REWIND (LINKMC)
         IFLAG = 1
         RETURN
      ENDIF
C
      REWIND LINKMC
C*****linkfile (transport) > binary
C      LBIN = .TRUE.
C*****END linkfile (transport) > binary
C*****linkfile (transport) > ascii
      LBIN = .FALSE.
C*****END linkfile (transport) > ascii
C
      NREC = 1
      IF (LBIN) THEN
         READ (LINKMC, ERR=999) VERS
         NREC = 2
         READ (LINKMC, ERR=999) PRVERS
         NREC = 3
         READ (LINKMC, ERR=999) PREC
         NREC = 4
         READ (LINKMC, ERR=999) KERR
         NREC = 5
         READ (LINKMC, ERR=999) LI, LR, NO, NKK, NLITE
         NREC = 6
         READ (LINKMC, ERR=999) PATMOS
      ELSE
         READ (LINKMC, CFMT, ERR=999) VERS
         NREC = 2
         READ (LINKMC, CFMT, ERR=999) PRVERS
         NREC = 3
         READ (LINKMC, CFMT, ERR=999) PREC
         NREC = 4
         READ (LINKMC, LFMT, ERR=999) KERR
         NREC = 5
         READ (LINKMC, IFMT, ERR=999) LI, LR, NO, NKK, NLITE
         NREC = 6
         READ (LINKMC, RFMT, ERR=999) PATMOS
      ENDIF
C
      NK  = NO*NKK
      NK2 = NO*NKK*NKK
      K2  = NKK*NKK
      K3  = 3*NKK
      K32 = K3*K3
      NKT = NO*NKK*NLITE
C
C     APPORTION THE REAL WORKSPACE:
C
C     molecular weights for the species
      NWT  = 1
C     the epsilon/k well depth for the species
      NEPS = NWT + NKK
C     the collision diameter for the species
      NSIG = NEPS + NKK
C     the dipole moments for the species
      NDIP = NSIG + NKK
C     the polarizabilities for the species
      NPOL = NDIP + NKK
C     the rotational relaxation collision numbers
      NZROT= NPOL + NKK
C     the fit coefficients for conductivity
      NLAM = NZROT + NKK
C     the fit coefficients for viscosity
      NETA = NLAM + NK
      NDIF = NETA + NK
C     the fit coefficients for thermal diffusion ratio
      NTDIF= NDIF + NK2
C     mole fractions of the mixture
      NXX  = NTDIF + NO*NKK*NLITE
C     species viscosities
      NVIS = NXX + NKK
C     rotational relaxation collision numbers before Parker coffection
      NXI  = NVIS + NKK
C     species specific heats
      NCP  = NXI + NKK
C     rotational parts of the specific heats
      NCROT= NCP + NKK
C     internal parts of the specific heats
      NCINT= NCROT + NKK
C     the binary diffusion coefficients
      NBIND= NCINT + NKK
C     the matrix of reduced well depths
      NEOK = NBIND + K2
C     the matrix of reduced collision diameters
      NSGM = NEOK + K2
C     the matrix of A* collision integrals for each species pair
      NAST = NSGM + K2
C     the matrix of B* collision integrals for each species pair
      NBST = NAST + K2
C     the matrix of C* collision integrals for each species pair
      NCST = NBST + K2
C     the "L" matrix
      NXL  = NCST + K2
C     the right-hand sides of the linear system involving the
C     "L" matrix
      NR   = NXL + K32
C     the workspace needed by LINPACK to solve the "L" matrix linear
C     system
      NWRK = NR + K3
C      NTOT = NWRK + K3 - 1
C
C     APPORTION THE INTEGER WORKSPACE:
C
C     the indicators for the molecule linearity
      INLIN = 1
C     the species indices for the "light" species
      IKTDIF= INLIN + NKK
C     the pivot indices for LINPACK calls
      IPVT  = IKTDIF + NLITE
C      ITOT  = IPVT + K3 - 1
C
C     Read the data from the linkfile
C
      IF (LBIN) THEN
         NREC = 7
         READ (LINKMC, ERR=999) (RMCWRK(NWT+N-1), N = 1, NKK)
         NREC = 8
         READ (LINKMC, ERR=999) (RMCWRK(NEPS+N-1), N = 1, NKK)
         NREC = 9
         READ (LINKMC, ERR=999) (RMCWRK(NSIG+N-1), N = 1, NKK)
         NREC = 10
         READ (LINKMC, ERR=999) (RMCWRK(NDIP+N-1), N = 1, NKK)
         NREC = 11
         READ (LINKMC, ERR=999) (RMCWRK(NPOL+N-1), N = 1, NKK)
         NREC = 12
         READ (LINKMC, ERR=999) (RMCWRK(NZROT+N-1), N = 1, NKK)
         NREC = 13
         READ (LINKMC, ERR=999) (IMCWRK(INLIN+N-1), N = 1, NKK)
         NREC = 14
         READ (LINKMC, ERR=999) (RMCWRK(NLAM+N-1), N = 1, NK)
         NREC = 15
         READ (LINKMC, ERR=999) (RMCWRK(NETA+N-1), N = 1, NK)
         NREC = 16
         READ (LINKMC, ERR=999) (RMCWRK(NDIF+N-1), N = 1, NK2)
         NREC = 17
         READ (LINKMC, ERR=999) (IMCWRK(IKTDIF+N-1), N = 1, NLITE)
         NREC = 18
         READ (LINKMC, ERR=999) (RMCWRK(NTDIF+N-1), N = 1, NKT)
      ELSE
         NREC = 7
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NWT+N-1), N = 1, NKK)
         NREC = 8
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NEPS+N-1), N = 1, NKK)
         NREC = 9
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NSIG+N-1), N = 1, NKK)
         NREC = 10
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NDIP+N-1), N = 1, NKK)
         NREC = 11
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NPOL+N-1), N = 1, NKK)
         NREC = 12
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NZROT+N-1), N = 1, NKK)
         NREC = 13
         READ (LINKMC, IFMT, ERR=999) (IMCWRK(INLIN+N-1), N = 1, NKK)
         NREC = 14
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NLAM+N-1), N = 1, NK)
         NREC = 15
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NETA+N-1), N = 1, NK)
         NREC = 16
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NDIF+N-1), N = 1, NK2)
         NREC = 17
         READ (LINKMC, IFMT, ERR=999) (IMCWRK(IKTDIF+N-1), N = 1, NLITE)
         NREC = 18
         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NTDIF+N-1), N = 1, NKT)
         NREC = 19
      ENDIF
C
C     Set EPS/K and SIG for all I,J pairs
C
      CALL MCEPSG (NKK, RMCWRK(NEPS), RMCWRK(NSIG), RMCWRK(NDIP),
     1            RMCWRK(NPOL), RMCWRK(NEOK), RMCWRK(NSGM) )
C
  300 FORMAT (10X,'IMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
  350 FORMAT (10X,'RMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
      RETURN
  999 CONTINUE
C <error module="tranlib" severity="error">
C <id>3</id>
C <message>An error was encountered while trying to read
C the TRANSPORT linking file.  Check the output from the TRAN
C fitting routine for errors or rerun TRAN.
C </message>
C </error>
      IDERR = 3
      MSGSTR = ' '
      CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
      WRITE (LOUT, *) ' Error reading Transport linkfile...'
      REWIND (LINKMC)
      IFLAG = NREC
C
C     end of SUBROUTINE MCINIT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCLEN (LINKMC, LOUT, LI, LR, IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE MCLEN (LINKMC, LOUT, LI, LR, IFLAG)
C  Returns the lengths required for work arrays.
C
C  INPUT
C  LINKMC   - Integer scalar, input file unit for the linkfile.
C  LOUT     - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  LI       - Integer scalar, minimum length required for the
C             integer work array.
C  LR       - Integer scalar, minimum length required for the
C             real work array.
C  IFLAG    - Integer scalar, indicates successful reading of
C             linkfile; IFLAG>0 indicates error type.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
      PARAMETER (NLIST = 1)
      LOGICAL KERR, VOK, POK, LBIN
      CHARACTER*16 LIST(NLIST), VERS, PREC, PRVERS, IFMT, RFMT, LFMT,
     1             CFMT
      CHARACTER*80 MSGSTR
      PARAMETER
     1(CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)', RFMT='(1P,5E24.16)')
      DATA LIST/'1.0'/
C
      IFLAG = 0
      VERS = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LI   = LENI
      LR   = LENR
C
C*****linkfile (transport) > binary
C      LBIN = .TRUE.
C*****END linkfile (transport) > binary
C*****linkfile (transport) > ascii
      LBIN = .FALSE.
C*****END linkfile (transport) > ascii
C
      REWIND (LINKMC)
      NREC = 1
      IF (LBIN) THEN
         READ (LINKMC, ERR=999) VERS
         NREC = 2
         READ (LINKMC, ERR=999) PRVERS
         NREC = 3
         READ (LINKMC, ERR=999) PREC
         NREC = 4
         READ (LINKMC, ERR=999) KERR
      ELSE
         READ (LINKMC, CFMT, ERR=999) VERS
         NREC = 2
         READ (LINKMC, CFMT, ERR=999) PRVERS
         NREC = 3
         READ (LINKMC, CFMT, ERR=999) PREC
         NREC = 4
         READ (LINKMC, LFMT, ERR=999) KERR
      ENDIF
C
      VOK = .FALSE.
      DO 5 N = 1, NLIST
         IF (VERS .EQ. LIST(N)) VOK = .TRUE.
    5 CONTINUE
C
      POK = .FALSE.
C*****precision > double
      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
C*****END precision > double
C*****precision > single
C      IF (INDEX(PREC, 'SING') .GT. 0) POK = .TRUE.
C*****END precision > single
C
      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
C <error module="tranlib" severity="error">
C <id>4</id>
C <message>There is an error in the transport linkfile...',
C Check TRAN output for error conditions.'
C </message>
C </error>
            IDERR = 4
            MSGSTR = ' '
            CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
            WRITE (LOUT,'(/A,/A)')
     1      ' There is an error in the transport linkfile...',
     2      ' Check TRANFIT output for error conditions.'
         ENDIF
         IF (.NOT. VOK) THEN
C <error module="tranlib" severity="error">
C <id>5</id>
C <message>The transport linkfile is incompatible with ',
C this version of the Transport Library.  Rerun the TRAN program.
C </message>
C </error>
            IDERR = 5
            MSGSTR = ' '
            CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
            WRITE (LOUT,'(/A,A)')
     1      ' Transport linkfile is incompatible with Transport',
     2      ' Library Version 3.6'
         ENDIF
         IF (.NOT. POK) THEN
C <error module="tranlib" severity="error">
C <id>6</id>
C <message>The floating-point precision in the transport
C linkfile is incompatible with that of the Transport Library.
C Rerun the TRAN program.
C </message>
C </error>
            IDERR = 6
            MSGSTR = ' '
            CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
            WRITE (LOUT, '(/A,A)')
     1      ' Precision of Transport linkfile does not agree with',
     2      ' precision of Transport Library'
         ENDIF
C <error module="tranlib" severity="error">
C <id>7</id>
C <message>Returning from MCLEN with fatal error condition.
C </message>
C </error>
         IDERR = 7
         MSGSTR = ' '
         CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
         IFLAG = 21
         RETURN
      ENDIF
C
      NREC = 5
      IF (LBIN) THEN
         READ (LINKMC, ERR=1111) LENIMC, LENRMC, NO, NKK, NLITE
      ELSE
         READ (LINKMC, IFMT, ERR=1111) LENIMC, LENRMC, NO, NKK, NLITE
      ENDIF
      REWIND (LINKMC)
      LENI = LENIMC
      LENR = LENRMC
      LI   = LENI
      LR   = LENR
      RETURN
C
  999 CONTINUE
C <error module="tranlib" severity="error">
C <id>8</id>
C <message>An error was encountered during READ of the
C records 1-4 of the transport linkfile in subroutine MCLEN.
C Rerun the TRAN program and check for errors in the output.
C </message>
C </error>
      IDERR = 8
      MSGSTR = ' '
      CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
      WRITE (LOUT, 50)
   50 FORMAT
     1 (' Error reading Transport linkfile.')
      REWIND (LINKMC)
      IFLAG = NREC
      RETURN
 1111 CONTINUE
C <error module="tranlib" severity="error">
C <id>9</id>
C <message>An error was encountered during READ of
C record #5 of the transport linkfile in subroutine MCLEN.
C Rerun the TRAN program and check for errors in the output.
C </message>
C </error>
      IDERR = 9
      MSGSTR = ' '
      CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
      WRITE (LOUT, 50)
      REWIND (LINKMC)
      IFLAG = 22
C
C     end of SUBROUTINE MCLEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCPNT (LSAVE, LOUT, NPOINT, V, P, LI, LR, IERR)
C
C  START PROLOGUE
C
C  SUBROUTINE MCPNT (LSAVE, LOUT, NPOINT, V, P, LI, LR, IERR)
C  Reads from a binary file information about a Transport linkfile,
C  pointers for the Transport Library, and returns lengths of work
C  arrays.
C
C  INPUT
C  LSAVE  - Integer scalar, input unit for binary data file.
C  LOUT   - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  NPOINT - Integer scalar, total number of pointers.
C  V      - Real scalar, version number of the Transport linkfile.
C  P      - Character string, machine precision of the linkfile.
C  LI     - Integer scalar, minimum dimension required for integer
C           workspace array.
C  LR     - Integer scalar, minimumm dimension required for real
C           workspace array.
C  IERR   - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT,  NEPS,  NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF,  NTDIF, NXX, NVIS, NXI,  NCP,  NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST, NCST,
     4                NXL, NR, NWRK, K3
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
C
      LOGICAL KERR, IERR
      CHARACTER PREC*16, VERS*16, P*16, V*16
      CHARACTER*80 MSGSTR
C
      KERR = .FALSE.
      IERR = .FALSE.
C
      READ (LSAVE, ERR=100) NPOINT, VERS, PREC, LENI, LENR,
     *                RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP,
     3                NCROT, NCINT, NBIND, NEOK, NSGM,
     4                NAST, NBST, NCST, NXL, NR, NWRK, K3
      LI = LENI
      LR = LENR
      IERR = KERR
      V = VERS
      P = PREC
      RETURN
C
  100 CONTINUE
C <error module="tranlib" severity="error">
C <id>10</id>
C <message>An error was encountered during READ of the
C the transport linkfile in subroutine MCPNT.
C Rerun the TRAN program and check for errors in the output.
C </message>
C </error>
      IDERR = 10
      MSGSTR = ' '
      CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
      WRITE (LOUT, *) ' Error reading Transport linkfile data...'
      KERR   = .TRUE.
      IERR   = KERR
      NPOINT = 0
      VERS   = ' '
      V      = VERS
      PREC   = ' '
      P      = PREC
C
C     end of SUBROUTINE MCPNT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCPRAM (IMCWRK, RMCWRK, EPS, SIG, DIP, POL, ZROT, NLIN)
C
C  START PROLOGUE
C
C  SUBROUTINE MCPRAM (IMCWRK, RMCWRK, EPS, SIG, DIP, POL, ZROT, NLIN)
C  Returns the arrays of molecular parameters as read from the
C  transport database.
C
C  INPUT
C  IMCWRK(*) - Integer workspace array; dimension at least LENIMC.
C  RMCWRK(*) - Real    workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  EPS(*)    - Real array, Lennard-Jones Potential well depths for
C              the species;
C              dimension at least KK, the total species count.
C                 cgs units, K
C  SIG(*)    - Real array, Lennary-Jones collision diameters for
C              the species;
C              dimension at least KK, the total species count.
C                 cgs units, Angstrom
C  DIP(*)    - Real array, dipole moments for the species;
C              dimension at least KK, the total species count.
C                 cgs units, Debye
C  POL(*)    - Real array, polarizabilities for the species;
C              dimension at least KK, the total species count.
C                 cgs units, Angstrom**3
C  ZROT(*)   - Real array, rotational collision numbers evaluated at
C              298K for the species;
C              dimension at least KK, the total species count.
C  NLIN(*)   - Integer array, flags for species linearity;
C              dimension at least KK, the total species count.
C              NLIN=0, single atom,
C              NLIN=1, linear molecule,
C              NLIN=2, linear molecule.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST, NCST,
     4                NXL, NR, NWRK, K3
      DIMENSION IMCWRK(*), RMCWRK(*), EPS(*), SIG(*), DIP(*), POL(*),
     1          ZROT(*), NLIN(*)
C
      DO 100 K = 1, NKK
         KIND = K - 1
         EPS(K) = RMCWRK(NEPS + KIND)
         SIG(K) = RMCWRK(NSIG + KIND)
         DIP(K) = RMCWRK(NDIP + KIND)
         POL(K) = RMCWRK(NPOL + KIND)
         ZROT(K)= RMCWRK(NZROT+ KIND)
         NLIN(K)= IMCWRK(INLIN+ KIND)
  100 CONTINUE
C
C     end of SUBROUTINE MCPRAM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCREWR (LINKMC, LOUT, IMCWRK, RMCWRK, IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE MCREWR (LINKMC, LOUT, IMCWRK, RMCWRK, IFLAG)
C  This subroutine writes a new the transport linkfile from
C  the data stored in the integer and real work arrays,
C  IMCWRK(*) and RMCWRK(*).
C
C  INPUT
C  LINKMC    - Integer scalar, transport linkfile output unit number.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  IMCWRK(*) - Integer workspace array; dimension at least LENIMC.
C  RMCWRK(*) - Real    workspace array; dimension at least LENRMC.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST, NCST, NXL,
     4                NR, NWRK, K3
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
C
      DIMENSION IMCWRK(*), RMCWRK(*)
      CHARACTER*16 PRVERS, VERS, PREC, IFMT, RFMT, CFMT, LFMT
      CHARACTER*80 MSGSTR
      PARAMETER
     1(CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)', RFMT='(1P,5E24.16)')
C
      LOGICAL KERR, LBIN
C
C     Write version number
C
C*****precision > double
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      PREC = 'SINGLE'
C*****END precision > single
C
      PRVERS ='4.6'
C
      REWIND LINKMC
C*****linkfile (transport) > binary
C      LBIN = .TRUE.
C*****END linkfile (transport) > binary
C*****linkfile (transport) > ascii
      LBIN = .FALSE.
C*****END linkfile (transport) > ascii
C
      NREC = 1
      IF (LBIN) THEN
         WRITE (LINKMC, ERR=999) VERS
         NREC = 2
         WRITE (LINKMC, ERR=999) PRVERS
         NREC = 3
         WRITE (LINKMC, ERR=999) PREC
         NREC = 4
         WRITE (LINKMC, ERR=999) KERR
         NREC = 5
         WRITE (LINKMC, ERR=999) LENI, LENR, NO, NKK, NLITE
         NREC = 6
         WRITE (LINKMC, ERR=999) PATMOS
      ELSE
         WRITE (LINKMC, CFMT, ERR=999) VERS
         NREC = 2
         WRITE (LINKMC, CFMT, ERR=999) PRVERS
         NREC = 3
         WRITE (LINKMC, CFMT, ERR=999) PREC
         NREC = 4
         WRITE (LINKMC, LFMT, ERR=999) KERR
         NREC = 5
         WRITE (LINKMC, IFMT, ERR=999) LENI, LENR, NO, NKK, NLITE
         NREC = 6
         WRITE (LINKMC, RFMT, ERR=999) PATMOS
      ENDIF
C
      NK  = NO*NKK
      NK2 = NO*NKK*NKK
      K2  = NKK*NKK
      K3  = 3*NKK
      K32 = K3*K3
      NKT = NO*NKK*NLITE
C
C     Write the data to the linkfile
C
      IF (LBIN) THEN
         NREC = 7
         WRITE (LINKMC, ERR=999) (RMCWRK(NWT+N-1), N = 1, NKK)
         NREC = 8
         WRITE (LINKMC, ERR=999) (RMCWRK(NEPS+N-1), N = 1, NKK)
         NREC = 9
         WRITE (LINKMC, ERR=999) (RMCWRK(NSIG+N-1), N = 1, NKK)
         NREC = 10
         WRITE (LINKMC, ERR=999) (RMCWRK(NDIP+N-1), N = 1, NKK)
         NREC = 11
         WRITE (LINKMC, ERR=999) (RMCWRK(NPOL+N-1), N = 1, NKK)
         NREC = 12
         WRITE (LINKMC, ERR=999) (RMCWRK(NZROT+N-1), N = 1, NKK)
         NREC = 13
         WRITE (LINKMC, ERR=999) (IMCWRK(INLIN+N-1), N = 1, NKK)
         NREC = 14
         WRITE (LINKMC, ERR=999) (RMCWRK(NLAM+N-1), N = 1, NK)
         NREC = 15
         WRITE (LINKMC, ERR=999) (RMCWRK(NETA+N-1), N = 1, NK)
         NREC = 16
         WRITE (LINKMC, ERR=999) (RMCWRK(NDIF+N-1), N = 1, NK2)
         NREC = 17
         WRITE (LINKMC, ERR=999) (IMCWRK(IKTDIF+N-1), N = 1, NLITE)
         NREC = 18
         WRITE (LINKMC, ERR=999) (RMCWRK(NTDIF+N-1), N = 1, NKT)
      ELSE
         NREC = 7
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NWT+N-1), N = 1, NKK)
         NREC = 8
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NEPS+N-1), N = 1, NKK)
         NREC = 9
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NSIG+N-1), N = 1, NKK)
         NREC = 10
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NDIP+N-1), N = 1, NKK)
         NREC = 11
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NPOL+N-1), N = 1, NKK)
         NREC = 12
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NZROT+N-1), N = 1, NKK)
         NREC = 13
         WRITE (LINKMC, IFMT, ERR=999) (IMCWRK(INLIN+N-1), N = 1, NKK)
         NREC = 14
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NLAM+N-1), N = 1, NK)
         NREC = 15
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NETA+N-1), N = 1, NK)
         NREC = 16
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NDIF+N-1), N = 1, NK2)
         NREC = 17
         WRITE (LINKMC, IFMT, ERR=999)
     &         (IMCWRK(IKTDIF+N-1), N = 1, NLITE)
         NREC = 18
         WRITE (LINKMC, RFMT, ERR=999) (RMCWRK(NTDIF+N-1), N = 1, NKT)
         NREC = 19
      ENDIF
C
      RETURN
  999 CONTINUE
C <error module="tranlib" severity="error">
C <id>3</id>
C <message>An error was encountered while trying to write
C a new TRANSPORT linking file, in subroutine MCREWR.
C </message>
C </error>
      IDERR = 11
      MSGSTR = ' '
      CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
      WRITE (LOUT, *) ' Error writing new Transport linkfile...'
      REWIND (LINKMC)
      IFLAG = NREC
C
C     end of SUBROUTINE MCREWR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCSAVE (LOUT, LSAVE, IMCWRK, RMCWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE MCSAVE (LOUT, LSAVE, IMCWRK, RMCWRK)
C  Writes to a binary file information about a Transport linkfile,
C  pointers for the Transport library, and Transport work arrays.
C
C  INPUT
C  LOUT      - Integer scalar, formatted output file unit number.
C  LSAVE     - Integer scalar, unformatted output file unit number.
C  IMCWRK(*) - Integer workspace array; dimension at least LENIMC.
C  RMCWRK(*) - Real    workspace array; dimension at least LENRMC.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP,
     3                NCROT, NCINT, NBIND, NEOK, NSGM,
     4                NAST, NBST, NCST, NXL, NR, NWRK, K3
      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR
      DIMENSION IMCWRK(*), RMCWRK(*)
      CHARACTER VERS*16, PREC*16
      CHARACTER*80 MSGSTR
      LOGICAL KERR
C
      NPOINT = 41
      WRITE (LSAVE,ERR=999)   NPOINT, VERS,   PREC,   LENI,   LENR,
     *                RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP,
     3                NCROT, NCINT, NBIND, NEOK, NSGM,
     4                NAST, NBST, NCST, NXL, NR, NWRK, K3
      WRITE (LSAVE,ERR=999) (IMCWRK(L), L = 1, LENI)
      WRITE (LSAVE,ERR=999) (RMCWRK(L), L = 1, LENR)
C
C     end of SUBROUTINE MCSAVE
      RETURN
  999 CONTINUE
C
C <error module="tranlib" severity="error">
C <id>11</id>
C <message>An error was encountered during WRITE of the
C the transport linkfile information in subroutine MCSAVE.
C Rerun the TRAN program and check for errors in the output.
C </message>
C </error>
      IDERR = 12
      MSGSTR = ' '
      CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
      WRITE (LOUT, *)
     1    ' Error writing Transport linkfile information...'
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCSVIS (T, RMCWRK, VIS)
C
C  START PROLOGUE
C
C  SUBROUTINE MCSVIS (T, RMCWRK, VIS)
C  Returns the array of pure species viscosities, given temperature.
C
C  INPUT
C  T         - Real scalar, temperature.
C                 cgs units, K
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  VIS(*)    - Real array, species viscosities;
C              dimension at least KK, the total species count.
C                 cgs units, gm/cm*s
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION VIS(*), RMCWRK(*)
C
      ALOGT = LOG(T)
      IF (NO .EQ. 4) THEN
        CALL MCEVAL4 (ALOGT, NKK, RMCWRK(NETA), VIS)
      ELSE
        CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NETA), VIS)
      ENDIF
C
      DO 25 K = 1, NKK
         VIS(K) = EXP(VIS(K))
   25 CONTINUE
C
C     end of SUBROUTINE MCSVIS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCAVIS (T, X, RMCWRK, VISMIX)
C
C  START PROLOGUE
C
C  SUBROUTINE MCAVIS (T, X, RMCWRK, VISMIX)
C  Returns mixture viscosity, given temperature and species mole
C  fractions.  It uses modification of the Wilke semi-empirical
C  formulas.
C
C  INPUT
C  T         - Real scalar, temperature.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  VISMIX    - Real scalar, mixture viscosity.
C                 cgs units, gm/cm*s
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, EIGHT=8.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0, ONE=1.0, EIGHT=8.0)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION X(*), RMCWRK(*)
C
C     In the following call, the species molecular weights are stored
C     in RMCWRK(NWT) and the pure species viscosities are in
C     RMCWRK(NVIS)
C
      ALOGT = LOG(T)
      IF (NO .EQ. 4) THEN
        CALL MCEVAL4 (ALOGT, NKK, RMCWRK(NETA), RMCWRK(NVIS))
      ELSE
        CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NETA), RMCWRK(NVIS))
      ENDIF
C
      DO 25 K = 1, NKK
         RMCWRK(NVIS+K-1) = EXP(RMCWRK(NVIS+K-1))
   25 CONTINUE
C
      SUMO = ZERO
      DO 200 K = 1, NKK
C
         SUMI = ZERO
         DO 100 J = 1, NKK
            SUMI = SUMI + X(J) *
     $             (ONE +
     $              SQRT( RMCWRK(NVIS+K-1) / RMCWRK(NVIS+J-1) *
     $                    SQRT( RMCWRK(NWT+J-1)/RMCWRK(NWT+K-1))
     $                  )
     $             )**2   /
     $             SQRT(ONE + RMCWRK(NWT+K-1) / RMCWRK(NWT+J-1))
  100    CONTINUE
C
         SUMO = SUMO + X(K) * RMCWRK(NVIS+K-1) / SUMI
  200 CONTINUE
C
      VISMIX = SUMO * SQRT(EIGHT)
C
C     end of SUBROUTINE MCAVIS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCSCON (T, RMCWRK, CON)
C
C  START PROLOGUE
C
C  SUBROUTINE MCSCON (T, RMCWRK, CON)
C  Returns the array of pur species conductivities given temperature.
C
C  INPUT
C  T         - Real scalar, temperature.
C                cgs units, K
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  CON(*)    - Real array, species thermal conductivities;
C              dimension at least KK, the total species count.
C                cgs units, erg/cm*K*s
C  END PROLOGUE
C
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION CON(*), RMCWRK(*)
C
      ALOGT = LOG(T)
      IF (NO .EQ. 4) THEN
        CALL MCEVAL4 (ALOGT, NKK, RMCWRK(NLAM), CON)
      ELSE
        CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NLAM), CON)
      ENDIF
C
      DO 25 K = 1, NKK
         CON(K) = EXP(CON(K))
   25 CONTINUE
C
C     end of SUBROUTINE MCSCON
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCACON (T, X, RMCWRK, CONMIX)
C
C  START PROLOGUE
C
C  SUBROUTINE MCACON (T, X, RMCWRK, CONMIX)
C  Returns the mixture thermal conductivity given temperature and
C  species mole fractions.
C
C  INPUT
C  T         - Real scalar, temperature.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  CONMIX    - Real scalar, mixture thermal conductivity.
C                 cgs units, erg/cm*K*s
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0, ONE=1.0)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION X(*), RMCWRK(*)
C
C     In the following call, the pure species conductivities are in
C     RMCWRK(NXI)
C
      ALOGT = LOG(T)
      IF (NO .EQ. 4) THEN
        CALL MCEVAL4 (ALOGT, NKK, RMCWRK(NLAM), RMCWRK(NXI))
      ELSE
        CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(NLAM), RMCWRK(NXI))
      ENDIF
C
      SUM = ZERO
      SUMR = ZERO
      DO 100 K = 1, NKK
         RMCWRK(NXI+K-1) = EXP(RMCWRK(NXI+K-1))
         SUM =  SUM  + X(K)*RMCWRK(NXI+K-1)
         SUMR = SUMR + X(K)/RMCWRK(NXI+K-1)
  100 CONTINUE
C
      CONMIX = 0.5 * (SUM + ONE/SUMR)
C
C     end of SUBROUTINE MCACON
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCSDIF (P, T, KDIM, RMCWRK, DJK)
C
C  START PROLOGUE
C
C  SUBROUTINE MCSDIF (P, T, KDIM, RMCWRK, DJK)
C  Returns the binary diffusion coefficients given pressure and
C  temperature.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T         - Real scalar, temperature.
C                 cgs units, K
C  KDIM      - Integer scalar, actual first dimension of DJK(KDIM,KK).
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  DJK(*,*)  - Real matrix, binary diffusion coefficients;
C              dimension at least KK, the total species count, for
C              both the first and second dimensions.
C                 cgs units, cm**2/s
C              CJK(J,K) is the diffusion coefficient of species J
C              in species K.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION RMCWRK(*), DJK(KDIM,*)
C
      PFAC = PATMOS/P
      ALOGT = LOG(T)
C
C     Use the fact that DJK(J,K) = DJK(K,J) here.
C     Also, call a hard-wired horner's routine for NO=4.
C
      NOKK = NO * NKK
      ISTART = NDIF
      IF (NO .EQ. 4) THEN
        DO 100 K = 1, NKK
           CALL MCEVAL4 (ALOGT, K, RMCWRK(ISTART), DJK(1,K))
! note from evatt: mark's changes here are ok because the
! entire djk matrix is evaluated in the previous line 
#ifdef VECTORVERSION
CMRF
C Cray optimizations by Mark Fahey(ORNL)
           DO 90 J = 1, NKK
              DJK(J,K) = EXP(DJK(J,K)) * PFAC
   90      CONTINUE
#else
           DO 90 J = 1, K
              DJK(J,K) = EXP(DJK(J,K)) * PFAC
              DJK(K,J) = DJK(J,K)
   90      CONTINUE
#endif
           ISTART = ISTART + NOKK
  100   CONTINUE
      ELSE
C
         DO 200 K = 1, NKK
            CALL MCEVAL (ALOGT, NKK, NO, RMCWRK(ISTART), DJK(1,K))
            DO 190 J = 1, NKK            
               DJK(J,K) = EXP(DJK(J,K)) * PFAC
  190       CONTINUE
            ISTART = ISTART + NOKK
  200    CONTINUE
      ENDIF
C
C     end of SUBROUTINE MCSDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCADIF (P, T, X, RMCWRK, D)
C
C  START PROLOGUE
C
C  SUBROUTINE MCADIF (P, T, X, RMCWRK, D)
C  Returns mixture-averaged diffusion coefficients given pressure,
C  temperature, and species mole fractions.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T         - Real scalar, temperature.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  D(*)      - Real array, mixture diffusion coefficients;
C              dimension at least KK, the total species count.
C                 cgs units, cm**2/s
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION X(*), D(*), RMCWRK(*)
C
      CALL MCEDIF (T, NO, NKK, X, RMCWRK(NDIF), RMCWRK(NWT), SMALL,
     1             RMCWRK(NXX), RMCWRK(NBIND), D)
C
      PFAC = PATMOS / P
      DO 100 K = 1, NKK
         D(K) = D(K) * PFAC
  100 CONTINUE
C
C     end of SUBROUTINE MCADIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCEDIF(T, NO, KK, X, COFD, WT, SMALL, XX, DJK, D)
C
C  START PROLOGUE
C
C  SUBROUTINE MCEDIF(T, NO, KK, X, COFD, WT, SMALL, XX, DJK, D)
C  This subroutine is used internally to compute the mixture
C  diffusion coefficients; normally not called by the package user.
C
C  INPUT
C  T          - Real scalar, temperature.
C                  cgs units, K
C  NO         - Integer scalar, order of fit.
C  KK         - Integer scalar, total species count.
C  X(*)       - Real array, mole fractions of the mixture;
C               dimension at least KK, the total species count.
C  COFD(*,*,*)- Real three-dimensional array, coefficients of the
C               fits for binary diffusion coefficients;
C               dimension at least NO for the first dimension,
C               the fit order, and at least KK, the total species
C               count, for both the second and last dimensions.
C  WT(*)      - Real array, species molecular weights;
C               dimension at least KK, the total species count.
C  SMALL      - Real scalar, a small number added to all mole fractions
C               before computing the mixture diffusion coefficients;
C               this process avoids an undefined situation when a pure
C               species condition is approached.
C  XX(*)      - Real array, mole fractions plus SMALL to avoid the
C               problem of a pure species;
C               dimension at least KK, the total species count.
C  RMCWRK(*)  - Real workspace array; dimension at LENRMC.
C
C  OUTPUT
C  D(*)       - Real array, mixture diffusion coefficients.
C                  cgs units, cm**2/s.
C  DJK(*,*)   - Real matrix, binary diffusion coefficients;
C               dimension at least KK, the total species count, for
C               both dimensions.
C                  cgs units, cm**2/s
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO = 0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (ZERO = 0.0E0)
C*****END precision > single
C
      DIMENSION X(KK), COFD(NO,KK,KK), WT(KK), XX(KK),
     $          DJK(KK,KK), D(KK)
C
      ALOGT = LOG(T)
C
C Special Case for K = 1 - return the self-diffusion coefficient
C
      IF (KK .EQ. 1) THEN
        CALL MCEVAL (ALOGT, 1, NO, COFD(1,1,1), DJK(1,1))
        D(1) = EXP(DJK(1,1))
        RETURN
      ENDIF
C
C Use the fact that DJK is symmetric to cut the work down by 1/2
C  - also we don't need the self-diffusion coefficient evaluated
C
      IF (NO .EQ. 4) THEN
#ifdef VECTORVERSION
C  Added by Mark Fahey, ORNL.
C  Manually inline for Cray.
C  BUG FIX Evatt Hawkes - the whole DJK must be evaluated
C        DO 90 K = 2, KK
C           DO J = 1, K-1
        DO 90 K = 1, KK
           DO J = 1, KK
             DJK(J,K) = (((COFD(4,J,K) * ALOGT) + COFD(3,J,K)) * ALOGT +
     $                        COFD(2,J,K)) * ALOGT + COFD(1,J,K)
           ENDDO
#else
        DO 90 K = 2, KK
         CALL MCEVAL4 (ALOGT, K-1, COFD(1,1,K), DJK(1,K) )
#endif
 90     CONTINUE
      ELSE
#ifdef VECTORVERSION
C  Added by Mark Fahey, ORNL.
C  Manually inline for Cray.
C  BUG FIX Evatt Hawkes - the whole DJK must be evaluated
C        DO 100 K = 2, KK 
C          DO J = 1, K-1
        DO 100 K = 1, KK 
          DO J = 1, KK
            DJK(J,K) = COFD(NO,J,K)
          ENDDO
          DO I = 1, NO-1
C            DO J = 1, K-1
            DO J = 1, KK
              DJK(J,K) = DJK(J,K) * ALOGT + COFD(NO-I,J,K)
            ENDDO
          ENDDO
#else
        DO 100 K = 2, KK
          CALL MCEVAL (ALOGT, K-1, NO, COFD(1,1,K), DJK(1,K) )
#endif
100     CONTINUE
      ENDIF
C
C Fill in the entire DJK, only after the exponential !
C - actually, evaluate and store the inverse of the
C   binary diffusion coefficient - this is what's needed.
C
#ifdef VECTORVERSION
C Added by Mark Fahey, ORNL.
C If DJK is symmetric, can we do something more clever
      DO K = 1, KK
         DO J = 1, KK
            DJK(J,K) = EXP(-DJK(J,K))
         ENDDO
      ENDDO
      DO K = 1, KK
         DJK(K,K) = ZERO
      ENDDO
#else
      DO 150 K = 1, KK
         DO 140 J = 1, K-1
            DJK(J,K) = EXP(-DJK(J,K))
            DJK(K,J) = DJK(J,K)
  140    CONTINUE
         DJK(K,K) = ZERO
  150 CONTINUE
#endif
C
      WTM = ZERO
      DO 175 K = 1, KK
         WTM = WTM + WT(K)*X(K)
         XX(K) = MAX (X(K), SMALL)
  175 CONTINUE
C
      DO 300 K = 1, KK
C
         SUMXW  = - XX(K) * WT(K)
         SUMXOD = ZERO
C
         DO 200 J = 1, KK
             SUMXW  = SUMXW  + XX(J)*WT(J)
             SUMXOD = SUMXOD + XX(J)*DJK(J,K)
  200    CONTINUE
C
         D(K) = SUMXW/(WTM*SUMXOD)
  300 CONTINUE
C
C     end of SUBROUTINE MCEDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCMDIF (P, T, X, KDIM, IMCWRK, RMCWRK, D)
C
C  START PROLOGUE
C
C  SUBROUTINE MCMDIF (P, T, X, KDIM, IMCWRK, RMCWRK, D)
C  Returns the ordinary multicomponent diffusion coefficients,
C  given pressure, temperature, and mole fractions.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T         - Real scalar, temperature.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  KDIM      - Integer scalar, actual first dimension of D(KDIM,KK);
C              KDIM must be at least KK, the total species count.
C  IMCWRK(*) - Integer workspace array; dimension at least LENIMC.
C  RMCWRK(*) - Real    workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  D(*,*)    - Real matrix, ordinary multicomponent diffusion
C              coefficients;
C              dimension at least KK, the total species count, for
C              both the first and second dimensions.
C                 cgs units, cm**2/s
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION X(*), IMCWRK(*), RMCWRK(*), D(KDIM,*)
      LOGICAL   KERR
C
      CALL MCORDF (P, T, X, NKK, KDIM, SMALL, RMCWRK(NWT), RMCWRK,
     1             RMCWRK(NXX), RMCWRK(NBIND), RMCWRK(NXL),
     2             RMCWRK(NWRK), IMCWRK(IPVT), D, KERR)
C
C     end of SUBROUTINE MCMDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCMCDT (P, T, X, IMCWRK, RMCWRK, ICKWRK, CKWRK,
     1                   DT, COND)
C
C  START PROLOGUE
C
C  SUBROUTINE MCMCDT (P, T, X, IMCWRK, RMCWRK, ICKWRK, CKWRK,
C                     DT, COND)
C  Returns thermal diffusion coefficients, and mixture thermal
C  conductivities, given pressure, temperature, and mole fraction.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T         - Real scalar, temperature.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C
C  IMCWRK(*) - Integer TRANSPORT workspace array;
C              dimension at least LENIMC.
C  RMCWRK(*) - Real    TRANSPORT workspace array;
C              dimension at least LENRMC.
C  ICKWRK(*) - Integer CHEMKIN workspace array;
C              dimension at least LENICK.
C  RCKWRK(*) - Real    CHEMKIN workspace array;
C              dimension at least LENRCK.
C
C  OUTPUT
C  DT(*)     - Real array, thermal multicomponent diffusion
C              coefficients;
C              dimension at least KK, the total species count.
C                 cgs units, gm/(cm*sec)
C                  CGS UNITS - GM/(CM*SEC)
C  COND      - Real scalar, mixture thermal conductivity.
C                 cgs units, erg/(cm*K*s)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION X(*), IMCWRK(*), RMCWRK(*), ICKWRK(*), CKWRK(*), DT(*)
      LOGICAL   KERR
C
      CALL MCLMDT (P, T, X, NKK, K3, SMALL, RMCWRK(NWT), RMCWRK(NEOK),
     1             RMCWRK(NZROT), IMCWRK(INLIN), RMCWRK(NEPS),
     2             ICKWRK, CKWRK, RMCWRK, RMCWRK(NXX), RMCWRK(NVIS),
     3             RMCWRK(NAST), RMCWRK(NBST), RMCWRK(NCST),
     4             RMCWRK(NXI),  RMCWRK(NCP), RMCWRK(NCROT),
     5             RMCWRK(NCINT), RMCWRK(NXL), RMCWRK(NR),
     6             RMCWRK(NBIND), IMCWRK(IPVT), DT, COND, KERR)
C
C     end of SUBROUTINE MCMCDT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCEVAL (TF, KK, NO, COF, VAL)
C
C  START PROLOGUE
C
C  SUBROUTINE MCEVAL (TF, KK, NO, COF, VAL)
C  This subroutine uses Horner's algorithm to evaluate a polynomial
C  fit.  This routine is not normally called by the package user.
C
C  INPUT
C  TF        - Real scalar, independent variable of fit;
C              either temperature or log(temperature).
C  KK        - Integer scalar, total species count.
C  NO        - Integer scalar, order of fit.
C  COF(*,*)  - Real matrix, fit coefficients;
C              dimension exactly NO for the first dimension and at
C              least KK for the second.
C              COF(N,K) is the Nth coefficient of a property fit for
C              species K.
C
C  OUTPUT
C  VAL(*)    - Real array, evaluations of the fit at TF;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION COF(NO,KK), VAL(KK)
C
      NOM1 = NO-1
C
      DO 10 K = 1, KK
        VAL(K) = COF(NO,K)
   10 CONTINUE
      DO 200 I = 1, NOM1
        DO 150 K = 1, KK
          VAL(K) = VAL(K) * TF + COF(NO-I,K)
  150   CONTINUE
200   CONTINUE
C
C     end of SUBROUTINE MCEVAL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCEVAL4 (TF, KK, COF, VAL)
C
C  START PROLOGUE
C
C  SUBROUTINE MCEVAL4 (TF, KK, COF, VAL)
C  This subroutine uses Horner's algorithm to evaluate a polynomial;
C  the polynomial fit is hard-coded for order four polynomials, NO=4.
C  This routine is not normally called by the package user.
C
C  INPUT
C  TF       - Real scalar, independent variable of fit;
C             either temperature or log(temperature).
C  KK       - Integer scalar, total species count.
C  NO       - Integer scalar, order of fit.
C  COF(*,*) - Real matrix, fit coefficients;
C             dimension exactly NO for the first dimension, and at
C             least KK, the total species count, for the second;
C             COF(N,K) is the Nth coefficient of a property fit for
C             species K.
C
C  OUTPUT
C  VAL(*)   - Real array, evaluations of the fit at TF;
C             dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NO=4)
      DIMENSION COF(NO, KK), VAL(KK)
C
      DO 10 K = 1, KK
        VAL(K) = (((COF(4,K) * TF) + COF(3,K)) * TF + COF(2,K))
     $                                         * TF + COF(1,K)
10    CONTINUE
C
C     end of SUBROUTINE MCEVAL4
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCEPSG (KK, EPS, SIG, DIP, POL, EOK, SGM)
C
C  START PROLOGUE
C
C  SUBROUTINE MCEPSG (KK, EPS, SIG, DIP, POL, EOK, SGM)
C  This subroutine computes the reduced wel depth EOK(I,J) and
C  collision diameter SGM(I,J) for each I,J species pair.  The routine
C  is called only once, by the initialization subroutine MCINIT.
C  This routine is normally not called by the user.
C
C  INPUT
C  KK      - Integer scalar, total species count.
C  EPS(*)  - Real array, Lennard-Jones potential well depths;
C            dimension at least KK, the total species count.
C               cgs units, K
C  SIG(*)  - Real array, Lennardl-Jones collision diameters;
C            dimension at least KK, the total species count.
C               cgs units, Angstrom
C  DIP(*)  - Real array, dipole moments;
C            dimension at least KK, the total species count.
C               cgs units, Debye
C  POL(*)  - Real array, polarizabilities;
C            dimension at least KK, the total species count.
C               cgs units, Angstrom**3
C
C  OUTPUT
C  EOK(*,*)- Real matrix, reduced well depths for species pairs;
C            dimension at least KK, the total species count, for
C            both the first and second dimensions.
C               cgs units, K
C  SGM(*,*)- Real matrix, reduced collision diameters for species
C            species pairs;
C            dimension at least KK, the total species count, for
C            both the first and second dimensions.
C               cgs units, Angstrom
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ONE=1.0D0, FDTCGS=1.0D-18, FATCM=1.0D8,
     1           DIPMIN=1.0D-20, BOLTZ=1.38056D-16)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ONE=1.0, FDTCGS=1.0E-18, FATCM=1.0E8,
C     1           DIPMIN=1.0E-20, BOLTZ-1.38056E-16)
C*****END precision > single
C
      DIMENSION EPS(*), SIG(*), DIP(*), POL(*), EOK(KK,*), SGM(KK,*)
C
C     Compute and store EPS/K and SIGMA for all pairs
C
      DO 1000 J = 1, KK
C
         DO 999 K = 1, J
C
           IF(DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN) THEN
C             K is polar, J is nonpolar
C
              XI = ONE + 0.25*(POL(J)/SIG(J)**3) *
     1                     (FDTCGS**2*FATCM**3/BOLTZ) *
     2                     (DIP(K)**2/(EPS(K)*SIG(K)**3)) *
     3                      SQRT(EPS(K)/EPS(J))
              SGM(K,J) = 0.5 * (SIG(J)+SIG(K)) * XI**(-ONE/6.0)
              EOK(K,J) = SQRT(EPS(J)*EPS(K)) * XI**2
C
          ELSE IF(DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN) THEN
C             J is polar, K is nonpolar
C
              XI = ONE + 0.25*(POL(K)/SIG(K)**3) *
     1                     (FDTCGS**2*FATCM**3/BOLTZ) *
     2                     (DIP(J)**2/(EPS(J)*SIG(J)**3)) *
     3                      SQRT(EPS(J)/EPS(K))
              SGM(K,J) = 0.5 * (SIG(J)+SIG(K)) * XI**(-ONE / 6.0)
              EOK(K,J) = SQRT(EPS(J)*EPS(K)) * XI**2
C
          ELSE
C             normal case, either both polar or both nonpolar
C
              SGM(K,J) = 0.5 * (SIG(J) + SIG(K))
              EOK(K,J) = SQRT(EPS(J)*EPS(K))
C
          ENDIF
          SGM(J,K) = SGM(K,J)
          EOK(J,K) = EOK(K,J)
 999    CONTINUE
1000  CONTINUE
C
C     end of SUBROUTINE MCEPSG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCORDF (P, T, X, KK, KDIM, SMALL, WT, RMCWRK, XX,
     1                   BINDIF, XL0000, WORK, IPVT, D, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE MCORDF (P, T, X, KK, KDIM, SMALL, WT, RMCWRK, XX,
C 1                   BINDIF, XL0000, WORK, IPVT, D, KERR)
C  This subroutine computes ordinary multicomponent diffusion coeffi-
C  cient matrix; it does so by computing the inverse of the L00,00
C  matrix.
C  This routine is not normally called directly by the user;
C  the user calls MCMDIF, which in turn calls MCORDF.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T         - Real scalar, temperature.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  KK        - Integer scalar, total species count.
C  KDIM      - Integer scalar, actual first dimension of D(KDIM,KK);
C              KDIM must be at least KK, the total species count.
C  SMALL     - Real scalar;  the mole fractions used in the transport
C              computation are given by XX(K) = X(K) + SMALL.
C  WT(*    - Real array, species molecular weights;
C              dimension at least KK, the total species count.
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C  XX(*)     - Real array, mole fractions used in transport computa-
C              tion.
C              XX(K) = X(K) + SMALL
C  BINDIF(*,*)-Real matrix, binary diffusion coefficients;
C              dimension at least KK, the total species count,
C              for both the first and second dimensions.
C                 cgs units, cm**2/s
C  XL0000(*,*)-Real matrix L00,00;
C              dimension at least KK, the total species count,
C              for both the first and second dimensions.
C  WORK(*)   - Real workspace array for the inversion of the L00,00
C              matrix by LINPACK routines SGEFA and SGEDI;
C              dimension at least KK, the total species count.
C  IPVT(*)   - Integer array, pivot indices for inversion of the
C              L00,00 matrix by LINPACK routines SGEFA and SGEDI;
C              dimension at least KK, the total species count.
C
C  OUTPUT
C  D(*,*)    - Real matrix, ordinary multicomponent diffusion
C              coefficients;
C              dimension at least KK, the total species count, for
C              both the first and second dimensions.
C                 cgs units, cm**2/s
C  KERR      - Logical flag indicating whether an error was found
C              during matrix manipulations
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0, ONE=1.0)
C*****END precision > single
C
      LOGICAL KERR
      CHARACTER*80 MSGSTR
      DIMENSION X(*), WT(*), BINDIF(KK,*), XL0000(KK,*), IPVT(*),
     1          WORK(*), XX(*), RMCWRK(*), D(KDIM,*)
C*****precision > double - linpack
C      DIMENSION DET(2)
C      SAVE      JOB
C      DATA      JOB /1/
C*****END precision > double - linpack
C*****precision > single - linpack
C      DIMENSION DET(2)
C      SAVE      JOB
C      DATA      JOB /1/
C*****END precision > single - linpack
C
C     Set minimum mole fraction to SMALL
C
      DO 50 I = 1, KK
         XX(I) = MAX( X(I) , SMALL)
   50 CONTINUE
C
C     Evaluate the binary diffusion coefficients
C
      CALL MCSDIF (P, T, KK, RMCWRK, BINDIF)
C
C     Assemble L00,00
C
      PFAC = 16.0 * T / (25.0 * P)
      DO 200 I = 1, KK
        SUM = -XX(I) / BINDIF(I,I)
        DO 90 J = 1, KK
          SUM = SUM + XX(J) / BINDIF(I,J)
  90    CONTINUE
        SUM = SUM / WT(I)
        DO 100 J = 1, KK
          XL0000(I,J) = PFAC * XX(J) *
     $                  (WT(J) * SUM + XX(I) / BINDIF(I,J))
 100    CONTINUE
        XL0000(I,I) = ZERO
 200  CONTINUE
C
C     Invert L00,00 using LAPACK or LINPACK
C
C*****precision > double - lapack
      CALL DGETRF (KK, KK, XL0000, KK, IPVT, INFO)
      IF (INFO .NE. 0) THEN
C <error module="tranlib" severity="error">
C <id>12</id>
C <message>An error occurred in subroutine MCORDF
C during LU factorization of the L-matrix, used to calculate
C ordinary diffusion coefficients.  Check that the transport
C properties for the species are correct.
C </message>
C <message level="2">The math Lapack subroutine
C DGETRF returned the value of INFO = %1.
C </message>
C </error>
          KERR = .TRUE.
          IDERR = 13
          MSGSTR = ' '
          WRITE (MSGSTR,'(I5)') INFO
          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
          WRITE (6,*) ' ERROR IN DGETRF, INFO = ', INFO
          RETURN
CEM         STOP
      ENDIF
      CALL DGETRI(KK, XL0000, KK, IPVT, WORK, KK, INFO)
      IF (INFO .NE. 0) THEN
C <error module="tranlib" severity="error">
C <id>13</id>
C <message>An error occurred in subroutine MCORDF
C during inversion of the L-matrix, used to calculate
C ordinary diffusion coefficients.  Check that the transport
C properties for the species are correct.
C </message>
C <message level="2">The math Lapack
C subroutine DGETRI returned the value of INFO = %1.
C </message>
C </error>
          KERR = .TRUE.
          IDERR = 14
          MSGSTR = ' '
          WRITE (MSGSTR,'(I5)') INFO
          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
          WRITE (6,*) ' ERROR IN DGETRI, INFO = ', INFO
          RETURN
CEM         STOP
      ENDIF
C*****END precision > double - lapack
C
C*****precision > double - linpack
C      CALL DGEFA (XL0000, KK, KK, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
CC <error module="tranlib" severity="error">
CC <id>14</id>
CC <message>An error occurred in subroutine MCORDF
CC during LU factorization of the L-matrix, used to calculate
CC ordinary diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Linpack subroutine
CC DGEFA returned the value of INFO = %1.
CC </message>
CC </error>
C        KERR = .TRUE.
C        IDERR = 15
C        MSGSTR = ' '
C        WRITE(MSGSTR,'(I5)') INFO
C        CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C        WRITE (6, *) ' ERROR IN DGEFA, INFO = ', INFO
C        STOP
C      ENDIF
C      CALL DGEDI (XL0000, KK, KK, IPVT, DET, WORK, JOB)
C*****END precision > double - linpack
C
C*****precision > single - lapack
C      CALL SGETRF (KK, KK, XL0000, KK, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
CC <error module="tranlib" severity="error">
CC <id>15</id>
CC <message>An error occurred in subroutine MCORDF
CC during LU factorization of the L-matrix, used to calculate
CC ordinary diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Lapack subroutine
CC SGETRF returned the value of INFO = %1.
CC </message>
CC </error>
C          KERR = .TRUE.
C          IDERR = 16
C          MSGSTR = ' '
C          WRITE (MSGSTR,'(I5)') INFO
C          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C          WRITE (6,*) ' ERROR IN SGETRF, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL SGETRI(KK, XL0000, KK, IPVT, WORK, KK, INFO)
CC <error module="tranlib" severity="error">
CC <id>16</id>
CC <message>An error occurred in subroutine MCORDF
CC during inversion of the L-matrix, used to calculate
CC ordinary diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Lapack
CC subroutine SGETRI returned the value of INFO = %1.
CC </message>
CC </error>
C      KERR = .TRUE.
C      IDERR = 17
C      MSGSTR = ' '
C      WRITE (MSGSTR,'(I5)') INFO
C      CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C      IF (INFO .NE. 0) THEN
C          WRITE (6,*) ' ERROR IN SGETRI, INFO = ', INFO
C          STOP
C      ENDIF
C*****END precision > single - lapack
C
C*****precision > single - linpack
C      CALL SGEFA (XL0000, KK, KK, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
CC <error module="tranlib" severity="error">
CC <id>17</id>
CC <message>An error occurred in subroutine MCORDF
CC during LU factorization of the L-matrix, used to calculate
CC ordinary diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Linpack subroutine
CC SGEFA returned the value of INFO = %1.
CC </message>
CC </error>
C         KERR = .TRUE.
C         IDERR = 18
C         MSGSTR = ' '
C         WRITE (MSGSTR,'(I5)') INFO
C         CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C        WRITE (6, *) ' ERROR IN SGEFA, INFO = ', INFO
C        STOP
C      ENDIF
C      CALL SGEDI (XL0000, KK, KK, IPVT, DET, WORK, JOB)
C*****END precision > single - linpack
C
C     Compute the ordinary multicomponent diffusion coefficients
C
      SUM = ZERO
      DO 400 I = 1, KK
        SUM = SUM + WT(I) * X(I)
  400 CONTINUE
      PFAC = PFAC * SUM
C
      DO 500 J = 1, KK
         PFAC_J = PFAC / WT(J)
         DO 450 I = 1, KK
            D(I,J) = PFAC_J * XX(I) * (XL0000(I,J)-XL0000(I,I))
  450    CONTINUE
  500 CONTINUE
C
C     end of SUBROUTINE MCORDF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCLMDT (P, T, X, KK, KK3, SMALL, WT, EOK, ZROT, LIN,
     1                   EPS, ICKWRK, CKWRK, RMCWRK, XX, VIS, ASTAR,
     2                   BSTAR, CSTAR, XI, CPOR, CROTOR, CINTOR, XL,
     3                   R, BINDIF, IPVT, DT, COND, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE MCLMDT (P, T, X, KK, KK3, SMALL, WT, EOK, ZROT, LIN,
C 1                   EPS, ICKWRK, CKWRK, RMCWRK, XX, VIS, ASTAR,
C 2                   BSTAR, CSTAR, XI, CPOR, CROTOR, CINTOR, XL,
C 3                   R, BINDIF, IPVT, DT, COND, KERR)
C  This subroutine computes the thermal conductivity, and the thermal
C  diffusion coefficient array; it does so by first forming the L
C  matrix, and then solving Eq. 24A.
C  This routine is not normally called directly by the user;
C  the user calls MCMCDT, which in turn calls MCLMDT.
C
C  INPUT
C  P        - Real scalar, pressure.
C                cgs units, dynes/cm**2
C  T        - Real scalar, temperature.
C                cgs units, K
C  X(*)     - Real array, mole fractions of the mixture;
C             dimension at least KK, the total species count.
C  KK       - Integer scalar, total species count.
C  KK3      - Integer scalar, three times the species count;
C             the size of the L matrix is KK3 * KK3.
C  SMALL    - Real scalar; the mole fractions used in the transport
C             computation are given by XX(K) = X(K) + SMALL.
C  WT(*)    - Real array, species molecular weights;
C             dimension at least KK, the total species count.
C  EOK(*,*) - Real matrix, reduced well depths for species pairs;
C             dimension at least KK, the total species count, for
C             both the first and second dimensions.
C                cgs units, K
C  ZROT(*)  - Real array, rotational collision numbers evaluated at
C             298K;
C             dimension at least KK, the total species count.
C  LIN(*)   - Integer array; flags indicating linearity of species;
C             dimension at least KK, the total species count.
C             NLIN=0, single atom,
C             NLIN=1, linear molecule,
C             NLIN=2, nonlinear molecule.
C  EPS(*)   - Real array, Lennard-Jones potential well depths;
C             dimension at least KK, the total species count.
C                cgs units, K
C  ICKWRK(*)- Integer CHEMKIN workspace array;
C             dimension at least LENICK.
C  RCKWRK(*)- Real   CHEMKIN workspace array;
C             dimension at least LENRCK.
C  RMCWRK(*)- Real   TRANSPORT workspace array;
C             dimension at least LENRMC.
C  XX(*)    - Real array, species mole fractions used in the
C             transport computation;
C             dimension at least KK, the total species count.
C             XX(K) = X(K) + SMALL.
C  VIS(*)   - Real array, species viscosities evaluated from MCSVIS;
C             dimension at least KK, the total species count.
C                cgs units, cm/cm-s
C  ASTAR(*,*) Real matrix, collision integrals A* for species pairs;
C             dimension at least KK, the total species count, for
C             both the first and second dimensions.
C  BSTAR(*,*) Real matrix, collision integrals B* for species pairs;
C             dimension at least KK, the total species count, for
C             both the first and second dimensions.
C  CSTAR(*,*) Real matrix, collision integrals C* for species pairs;
c             dimension at least KK, the total species count, for
C             both the first and second dimensions.
C  XI(*)    - Real array, collision numbers for the transfer of
C             rotational energy of species I into translational
C             energy of species J (Eq. 42), assuming that all
C             XI(I,J) = XI(I,I) (see p.132 for discussion);
C             dimension at least KK, the total species count.
C  CPOR(*)  - Real array, dimensionless specific heats CP/R for the
C             species;
C             dimension at least KK, the total species count.
C  CROT(*)  - Real array, dimensionless rotational contributions to
C             the specific heats of the species;
C             cimension at least KK, the total species count.
C  CINT(*)  - Real array, dimensionless internal contributions to
C             the specific heats of the species;
C             dimension at least KK, the total species count.
C  XL(*,*)  - Real matrix, the L matrix, Eq. 43 and 49;
C             dimension at least 3*KK, where KK is the total species
C             count, for both the first and second dimensions.
C  R(*)     - Real array, the right-hand sides of Eq. 24A;
C             dimension at least 3*KK, where KK is the total species
C             count.
C  BINDIF(*,*)-Real matrix, binary diffusion coefficients;
C             dimension at least KK, the total species count, for
C             both the first and second dimensions.
C                cgs units, cm**2/s
C  IPVT(*)  - Integer array, pivot indices for inversion of the XL
C             matrix by LINPACK routines SGEFA and SGEDI;
C             dimension at least 3*KK, where KK is the total species
C             count.
C
C  OUTPUT
C  DT(*)    - Real array, thermal multicomponent diffusion
C             coefficients;
C             dimension at least KK, the total species count.
C                cgs units, gm/(cm*sec)
C  COND     - Real scalar, mixture thermal conductivity.
C                cgs units, erg/(cm*K*s)
C  KERR     - Logical flag indicating whether an error was found
C             during matrix manipulations
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
      PARAMETER       (ONE = 1.0D0, ZERO = 0.0D0)
      PARAMETER (RU=8.314510D+07, PI= 3.1415926535897932D+00,
     $                        PI32O2= 2.7841639984158539D+00,
     $                        P2O4P2= 4.4674011002723397D+00,
     $                          PI32= 5.5683279968317078D+00)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ONE = 1.0E0, ZERO = 0.0E0)
C      PARAMETER (RU=8.314510E+07, PI= 3.1415926535897932E+00,
C     $                        PI32O2= 2.7841639984158539E+00,
C     $                        P2O4P2= 4.4674011002723397E+00,
C     $                          PI32= 5.5683279968317078E+00)
C*****END precision > single
C
      DIMENSION X(*), ICKWRK(*), CKWRK(*), RMCWRK(*), WT(*), XX(*),
     1          VIS(*), EOK(KK,*), ZROT(*), LIN(*), EPS(*),
     2          ASTAR(KK,*), BSTAR(KK,*), CSTAR(KK,*), XI(*), CPOR(*),
     3          CROTOR(*), CINTOR(*), XL(KK3,*), R(*), BINDIF(KK,*),
     4          IPVT(*), DT(*), FITAST(7), FITBST(7), FITCST(7)
      LOGICAL   KERR
      CHARACTER*80 MSGSTR
      SAVE      FITAST, FITBST, FITCST
C
C     Fits of A*, B*, and C* as functions of LN(T*)
C
      DATA FITAST / .1106910525E+01, -.7065517161E-02,
     1             -.1671975393E-01,  .1188708609E-01,
     2              .7569367323E-03, -.1313998345E-02,
     3              .1720853282E-03/
C
      DATA FITBST / .1199673577E+01, -.1140928763E+00,
     1             -.2147636665E-02,  .2512965407E-01,
     2             -.3030372973E-02, -.1445009039E-02,
     3              .2492954809E-03/
C
      DATA FITCST / .8386993788E+00,  .4748325276E-01,
     1              .3250097527E-01, -.1625859588E-01,
     2             -.2260153363E-02,  .1844922811E-02,
     3             -.2115417788E-03/
C
C     Set minimum mole fraction to SMALL
C     (Note, the possibility of negative mole fractions
C     necessitates the use of the MAX function ).
C
      DO 50 K = 1, KK
        XX(K) = MAX( X(K) ,  SMALL)
   50 CONTINUE
C
C     Determine A*, B*, and C* for each species
C     Note, these are symmetric matrices because EOK(I,J)
C     is symmetric
C
C      TLOG = LOG(T)
      DO 100 J = 1, KK
         DO 90 I = 1, J
C
         TSLOG = LOG ( T/EOK(I,J) )
C         TSLOG = TLOG - LOG(EOK(I,J))

         T1 = TSLOG
         T2 = TSLOG*T1
         T3 = TSLOG*T2
         T4 = TSLOG*T3
         T5 = TSLOG*T4
         T6 = TSLOG*T5
         ASTAR(I,J) = FITAST(1)    + FITAST(2)*T1 + FITAST(3)*T2 +
     1                FITAST(4)*T3 + FITAST(5)*T4 + FITAST(6)*T5 +
     2                FITAST(7)*T6
         ASTAR(J,I) = ASTAR(I,J)
         BSTAR(I,J) = FITBST(1)    + FITBST(2)*T1 + FITBST(3)*T2 +
     1                FITBST(4)*T3 + FITBST(5)*T4 + FITBST(6)*T5 +
     2                FITBST(7)*T6
         BSTAR(J,I) = BSTAR(I,J)
         CSTAR(I,J) = FITCST(1)    + FITCST(2)*T1 + FITCST(3)*T2 +
     1                FITCST(4)*T3 + FITCST(5)*T4 + FITCST(6)*T5 +
     2                FITCST(7)*T6
         CSTAR(J,I) = CSTAR(I,J)
   90    CONTINUE
  100 CONTINUE
C
C     Evaluate the binary diffusion coefficients and viscosity
C
      CALL MCSDIF (P, T, KK, RMCWRK, BINDIF)
      CALL MCSVIS (T, RMCWRK, VIS)
C
      PFAC = 1.2 * RU * T / P
      DO 150 K = 1, KK
C
C        Evaluate binary self-diffusion coefficients from viscosity
C
         BINDIF(K,K) = PFAC * ASTAR(K,K) * VIS(K) / WT(K)
C
C        Compute Parker correction for ZROT
C
         DD = EPS(K) / T
         DR = EPS(K) / 298.0
         SQRTDD = SQRT(DD)
         SQRTDR = SQRT(DR)
         DD32 = SQRTDD*DD
         DR32 = SQRTDR*DR
         XI(K) = ( (ONE + PI32O2*SQRTDR + P2O4P2*DR + PI32*DR32) /
     1             (ONE + PI32O2*SQRTDD + P2O4P2*DD + PI32*DD32)  )
     2            * MAX(ONE, ZROT(K))
  150 CONTINUE
C
C     Rotational and internal parts of specific heat
C
      CALL CKCPOR (T, ICKWRK, CKWRK, CPOR)
      DO 400 K = 1, KK
         IF (LIN(K) .EQ. 0) THEN
            CROTOR(K) = ZERO
            CINTOR(K) = ZERO
         ELSEIF (LIN(K) .EQ. 1) THEN
            CROTOR(K) = ONE
            CINTOR(K) = CPOR(K) - 2.5
         ELSEIF (LIN(K) .EQ. 2) THEN
            CROTOR(K) = 1.5
            CINTOR(K) = CPOR(K) - 2.5
         ENDIF
  400 CONTINUE
C
C     Assemble L00,00
C
      PFAC = 16.0 * T / (25.0 * P)
      DO 600 I = 1, KK
         SUM = - XX(I) / BINDIF(I,I)
         DO 450 K = 1, KK
           SUM = SUM + XX(K) / BINDIF(I,K)
  450    CONTINUE
         SUM = SUM / WT(I)
         DO 500 J = 1, KK
           XL(I,J) =   PFAC * XX(J) *
     $                        (WT(J) * SUM + XX(I) / BINDIF(J,I))
500      CONTINUE
         XL(I,I) = ZERO
  600 CONTINUE
C
C     Assemble L00,10 and L10,00
C
      PFAC = 8.0 * T / (5.0 * P)
      DO 1200 J = 1, KK
         WTJ_TMP = WT(J)
         XJ_TMP  = X(J)
         SUM     = ZERO
         DO 1150 I = 1, KK
            XL(I, J+KK) = -  PFAC * XX(I) * XJ_TMP * WT(I)
     1                       * (1.2*CSTAR(J,I) - ONE) /
     2                         ((WTJ_TMP + WT(I)) * BINDIF(J,I))
            XL(J+KK, I) = XL(I, J+KK)
            SUM = SUM   + XL(I, J+KK)
 1150    CONTINUE
         XL(J, J+KK) = XL(J, J+KK) - SUM
         XL(J+KK, J) = XL(J, J+KK)
 1200 CONTINUE
C
C     Assemble L01,00 and L00,01
C
      DO 1400 J = 1, KK
         DO 1300 I = 1, KK
            XL(2*KK+I, J) = ZERO
            XL(I, 2*KK+J) = ZERO
 1300    CONTINUE
 1400 CONTINUE
C
C     Assemble diagonal and off-diagonal elements of L10,10
C
      PFAC = 16.0D0 * T / (25.0 * P)
      PIFAC = 5.0 / (3.0*PI)
      DO 1600 J = 1, KK
        WTJ_TMP = WT(J)
        CROT_J  = CROTOR(J) / XI(J)
        PFAC_J  = PFAC * XX(J) * WTJ_TMP
        SUM     = ZERO
        DO 1550 I = 1, KK
          FAC_1 = XX(I) / ((WT(I) + WTJ_TMP)**2 * BINDIF(I,J))
          FAC_2 = 4.0*ASTAR(I,J)*
     $              (ONE + PIFAC*(CROTOR(I)/XI(I) + CROT_J))
          XL(I+KK, J+KK) = PFAC_J * WT(I) * FAC_1
     $                    * ( 13.75 - 3.0*BSTAR(I,J) - FAC_2 )
          SUM = SUM + FAC_1
     $              * (   7.5*WTJ_TMP**2
     $                  + WT(I)*(  WT(I)*(6.25 - 3.0*BSTAR(J,I))
     $                           + WTJ_TMP * FAC_2 )
     $                )
 1550   CONTINUE
        XL(J+KK, J+KK) = XL(J+KK, J+KK) - PFAC*XX(J)*SUM
 1600 CONTINUE
C
C     Assemble L10,01 and L01,10, both the off-diagonal entries
C     and the on-diagonal entries.
C
      NN = 2*KK
      PFAC = 32.0 * T / (5.0 * PI * P)
      DO 1850 J = 1, KK
         IF (LIN(J) .NE. 0) THEN
            NN = NN + 1
            SUM = ZERO
            WTJ_TMP = WT(J)
            PFAC_J =   ( PFAC * WTJ_TMP * XX(J) * CROTOR(J) )
     $               / ( CINTOR(J) * XI(J) )
            DO 1800 I = 1, KK
C                             The L10,01 term:
              XL(I+KK, NN) = ( PFAC_J * ASTAR(J,I) * XX(I)     )
     $                     / ( (WTJ_TMP + WT(I)) * BINDIF(J,I) )
C                             The L01,10 term:
              XL(NN, I+KK) =  XL(I+KK, NN)
C                             The extra term that get's stuck
C                             on the diagonal:
              SUM    = SUM +  XL(I+KK, NN)
 1800       CONTINUE
C
C           Extra diagonal entries:
C               (These use the viscosity, eq. 49, in their formulation,
C                because the self-diffusion coefficient has been
C                reevaluated to be consistent with the viscosity.)
C
            XL(J+KK, NN) = XL(J+KK, NN) + SUM
            XL(NN, J+KK) = XL(NN, J+KK) + SUM
         ENDIF
 1850 CONTINUE
C
C     Assemble L01,01 using viscosity Eq. (49).
C
      DO 2000 J = 1, KK
         DO 1900 I = 1, KK
            XL(2*KK+I, 2*KK+J) = ZERO
 1900    CONTINUE
 2000 CONTINUE
C
      NN = 2*KK
      PFAC  = 4.0 * T / P
      PIFAC = 12.0 / (5.0 * PI)
      PIRU  = - 8.0 / (PI * RU)
      DO 2200 I = 1, KK
         IF (LIN(I) .NE. 0) THEN
            NN = NN + 1
            SUM = ZERO
            FAC_1 = ( PIFAC * WT(I) * CROTOR(I) )
     $            / ( CINTOR(I) * XI(I)         )
            DO 2100 K = 1, KK
               FAC_2 = XX(K) / BINDIF(I,K)
               SUM = SUM + FAC_2
               IF (I .NE. K) THEN
                  SUM = SUM + (FAC_1 * FAC_2 * ASTAR(I,K)) / WT(K)
               ENDIF
 2100       CONTINUE
            FAC_2 = XX(I) / CINTOR(I)
            XL(NN, NN) =     ( PIRU * WT(I) * FAC_2 * CROTOR(I) )
     1                     / ( VIS(I) * XI(I)                   )
     2                   - PFAC * SUM
            XL(NN, NN) = FAC_2 * XL(NN, NN)
         ENDIF
 2200 CONTINUE
C
C     Assemble the right-hand side for solving Eq. (24)
C
      NN = 2*KK
      DO 3300 I = 1, KK
         R(I)    = ZERO
         R(I+KK) = XX(I)
         IF (LIN(I) .NE. 0) THEN
            NN  = NN + 1
            R(NN) = XX(I)
         ENDIF
 3300 CONTINUE
C
C     Factor and solve Eq. (24).
C
C*****precision > double - lapack
      CALL DGETRF (NN, NN, XL, KK3, IPVT, INFO)
      IF (INFO .NE. 0) THEN
C <error module="tranlib" severity="error">
C <id>18</id>
C <message>An error occurred in subroutine MCLMDT
C during LU factorization of the L-matrix, used to calculate
C thermal diffusion coefficients.  Check that the transport
C properties for the species are correct.
C </message>
C <message level="2">The math Lapack subroutine
C DGETRF returned the value of INFO = %1.
C </message>
C </error>
          IDERR = 19
          MSGSTR = ' '
          WRITE (MSGSTR,'(I5)') INFO
          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
          KERR = .TRUE.
          WRITE (6,*) ' ERROR IN DGETRF, INFO = ', INFO
          RETURN
CEM          STOP
      ENDIF
      CALL DGETRS('N', NN, 1, XL, KK3, IPVT, R, NN, INFO)
      IF (INFO .NE. 0) THEN
C <error module="tranlib" severity="error">
C <id>19</id>
C <message>An error occurred in subroutine MCLMDT
C during solution of the L-matrix equation for determination of
C thermal diffusion coefficients.  Check that the transport
C properties for the species are correct.
C </message>
C <message level="2">The math Lapack subroutine
C DGETRS returned the value of INFO = %1.
C </message>
C </error>
          IDERR = 20
          MSGSTR = ' '
          WRITE (MSGSTR,'(I5)') INFO
          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
          KERR = .TRUE.
          WRITE (6,*) ' ERROR IN DGETRS, INFO = ', INFO
          RETURN
CEM          STOP
      ENDIF
C*****END precision > double - lapack
C
C*****precision > double - linpack
C      CALL DGEFA (XL, KK3, NN, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
CC <error module="tranlib" severity="error">
CC <id>20</id>
CC <message>An error occurred in subroutine MCLMDT
CC during LU factorization of the L-matrix, used to calculate
CC thermal diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Linpack subroutine
CC DGEFA returned the value of INFO = %1.
CC </message>
CC </error>
C          IDERR = 21
C          MSGSTR = ' '
C          WRITE (MSGSTR,'(I5)') INFO
C          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C          WRITE (6,*) ' ERROR IN DGEFA, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL DGESL (XL, KK3, NN, IPVT, R, 0)
C*****END precision > double - linpack
C
C*****precision > single - lapack
C      CALL SGETRF (NN, NN, XL, KK3, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
CC <error module="tranlib" severity="error">
CC <id>21</id>
CC <message>An error occurred in subroutine MCLMDT
CC during LU factorization of the L-matrix, used to calculate
CC thermal diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Lapack subroutine
CC SGETRF returned the value of INFO = %1.
CC </message>
CC </error>
C          IDERR = 22
C          MSGSTR = ' '
C          WRITE (MSGSTR,'(I5)') INFO
C          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C          WRITE (6,*) ' ERROR IN SGETRF, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL SGETRS('N', NN, 1, XL, KK3, IPVT, R, NN, INFO)
C      IF (INFO .NE. 0) THEN
CC <error module="tranlib" severity="error">
CC <id>22</id>
CC <message>An error occurred in subroutine MCLMDT
CC during solution of the L-matrix equation for determination of
CC thermal diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Lapack subroutine
CC SGETRS returned the value of INFO = %1.
CC </message>
CC </error>
C          IDERR = 23
C          MSGSTR = ' '
C          WRITE (MSGSTR,'(I5)') INFO
C          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C          WRITE (6,*) ' ERROR IN SGETRS, INFO = ', INFO
C          STOP
C      ENDIF
C*****END precision > single - lapack
C
C*****precision > single - linpack
C      CALL SGEFA (XL, KK3, NN, IPVT, INFO)
C      IF (INFO .NE. 0) THEN
CC <error module="tranlib" severity="error">
CC <id>23</id>
CC <message>An error occurred in subroutine MCLMDT
CC during LU factorization of the L-matrix, used to calculate
CC thermal diffusion coefficients.  Check that the transport
CC properties for the species are correct.
CC </message>
CC <message level="2">The math Linpack subroutine
CC SGEFA returned the value of INFO = %1.
CC </message>
CC </error>
C          IDERR = 24
C          MSGSTR = ' '
C          WRITE (MSGSTR,'(I5)') INFO
C          CALL ERSET( 'tranlib', IDERR, 3, MSGSTR )
C          WRITE (6,*) ' ERROR IN SGEFA, INFO = ', INFO
C          STOP
C      ENDIF
C      CALL SGESL (XL, KK3, NN, IPVT, R, 0)
C*****END precision > single - linpack
C
C     Form thermal diffusion coefficients
C
      PFAC = 1.6 / RU
      DO 4000 K = 1, KK
         DT(K)  = PFAC * WT(K) * XX(K) * R(K)
 4000 CONTINUE
C
C     Form the thermal conductivity
C
      CONDTR = ZERO
      DO 4100 K = 1, KK
         CONDTR = CONDTR + X(K) * R(KK+K)
 4100 CONTINUE
C
      NN = 2*KK
      CONDIN = ZERO
      DO 4200 K = 1, KK
         IF (LIN(K) .NE. 0) THEN
            NN = NN + 1
            CONDIN = CONDIN + X(K) * R(NN)
         ENDIF
 4200 CONTINUE
C
      COND = -4.0 * (CONDTR + CONDIN)
C
C     end of SUBROUTINE MCLMDT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCATDR (T, X, IMCWRK, RMCWRK, TDR)
C
C  START PROLOGUE
C
C  SUBROUTINE MCATDR (T, X, IMCWRK, RMCWRK, TDR)
C  This subroutine computes the thermal diffusion ratios for the light
C  species into the mixture.
C
C  INPUT
C  T         - Real scalar, temperature.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  IMCWRK(*) - Integer workspace array; dimension at least LENIMC.
C  RMCWRK(*) - Real    workspace array; dimension at least LENRMC.
C
C  OUTPUT
C  TDR(*)    - Real array, thermal diffusion ratios for the species;
C              dimension at least KK, the total species count.
C              TDR(K) = 0 for any species with molecular weight less
C              than 5.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO = 0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO = 0.0)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION X(*), IMCWRK(*), RMCWRK(*), TDR(*)
C
C     In this subroutine, temporary storage is assigned as:
C     a vector of the "fitted" parts of the thermal diffusion ratios
C     are stored in RMCWRK(NXI), specifically, the vector represents
C     the J components of TDR(J,K), where K is the light species
C
      DO 100 K = 1, NKK
         TDR(K) = ZERO
  100 CONTINUE
C
      IF (NO .EQ. 4) THEN
        DO 400 L = 1, NLITE
           K = IMCWRK(IKTDIF+L-1)
           ISTRT = NTDIF + (L-1)*NO*NKK
           CALL MCEVAL4 (T, NKK, RMCWRK(ISTRT), RMCWRK(NXI))
           DO 350 J = 1, NKK
              TDR(K) = TDR(K) + RMCWRK(NXI+J-1)*X(K)*X(J)
  350      CONTINUE
  400   CONTINUE
        RETURN
      ENDIF
C
      DO 500 L = 1, NLITE
         K = IMCWRK(IKTDIF+L-1)
         ISTRT = NTDIF + (L-1)*NO*NKK
         CALL MCEVAL (T, NKK, NO, RMCWRK(ISTRT), RMCWRK(NXI))
         DO 450 J = 1, NKK
            TDR(K) = TDR(K) + RMCWRK(NXI+J-1)*X(K)*X(J)
  450    CONTINUE
  500 CONTINUE
C
C     end of SUBROUTINE MCATDR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCCCEX (K, KDIM, RMCWRK, COFCON)
C
C  START PROLOGUE
C
C  SUBROUTINE MCCCEX (K, KDIM, RCKWRK, COFCON)
C  Gets or puts values of the fitting coefficients for the
C  polynomial fits to species conductivity.
C
C  INPUT
C  K         - Integer scalar, species index.
C              K > 0 gets coefficients from RMCWRK
C              K < 0 puts coefficients into RMCWRK
C  KDIM      - Dimension for COFCON - the total number of species
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  If K < 1:
C  COFCON    - Real vector of polynomial coefficients for
C              the species' conductivity; dimension at least NO,
C              usually 4.
C
C  OUTPUT
C  If K > 1:
C  COFCON    - Real vector of polynomial coefficients for
C              the species' conductivity; dimension at least NO,
C              usually 4.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION RMCWRK(*), COFCON(*)
C
      NK = IABS(K)
      IF (K .GT. 0) THEN
C
C GET the data
C
         DO 200 N = 1, NO
            COFCON(N) = RMCWRK(NLAM + (NK-1)*NO + N - 1 )
200      CONTINUE
      ELSE IF (K .LT. 0) THEN
C
C PUT the data
C
         DO 400 N = 1, NO
            RMCWRK(NLAM + (NK-1)*NO + N - 1) = COFCON(N)
400      CONTINUE
      ENDIF
C
C     end of SUBROUTINE MCCCEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCCDEX (K, KDIM, RMCWRK, COFDIF)
C
C  START PROLOGUE
C
C  SUBROUTINE MCCDEX (K, KDIM, RCKWRK, COFDIF)
C  Gets or puts values of the fitting coefficients for the
C  polynomial fits to species binary diffusion coefficients.
C
C  INPUT
C  K         - Integer scalar, species index.
C              K > 0 gets coefficients from RMCWRK
C              K < 0 puts coefficients into RMCWRK
C  KDIM      - Dimension for COFDIF - the total number of species
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  If K < 1:
C  COFDIF    - Real matrix of polynomial coefficients for
C              the species' binary diffusion coefficient with all
C              other species; The first dimension should be KK;
C              the second dimension should be NO, usually 4;
C
C  OUTPUT
C  If K > 1:
C  COFDIF    - Real matrix of polynomial coefficients for
C              the species' binary diffusion coefficient with all
C              other species; first dimension should be NO, usually 4;
C              The second dimension should be NKK
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION RMCWRK(*), COFDIF(KDIM,*)
C
      NK = IABS(K)
      IF (K .GT. 0) THEN
C
C GET the data
C
         DO 200 N = 1, NO
C
C    get diffusion coefficients
C
            DO 150 J = 1, NKK
               COFDIF(J,N) =
     &           RMCWRK(NDIF + (NK-1)*NO*NKK + (J-1)*NO + N-1)
150         CONTINUE
200      CONTINUE
      ELSE IF (K .LT. 0) THEN
C
C PUT the data
C
         DO 400 N = 1, NO
C
C    put diffusion coeffs, keep binary diffusion matrix symmetric:
C
            DO 350 J = 1, NKK
               RMCWRK(NDIF + (NK-1)*NO*NKK + (J-1)*NO + N-1)
     &           = COFDIF(J,N)
               RMCWRK(NDIF + (J-1)*NO*NKK + (NK-1)*NO + N-1)
     &             = COFDIF(J,N)
350         CONTINUE
400      CONTINUE
      ENDIF
C
C     end of SUBROUTINE MCCDEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE MCCVEX (K, KDIM, RMCWRK, COFVIS)
C
C  START PROLOGUE
C
C  SUBROUTINE MCCVEX (K, KDIM, RCKWRK, COFVIS)
C  Gets or puts values of the fitting coefficients for the
C  polynomial fits to species viscosity.
C
C  INPUT
C  K         - Integer scalar, species index.
C              K > 0 gets coefficients from RMCWRK
C              K < 0 puts coefficients into RMCWRK
C  KDIM      - Dimension for COFVIS - the total number of species
C  RMCWRK(*) - Real workspace array; dimension at least LENRMC.
C
C  If K < 1:
C  COFVIS    - Real vector of polynomial coefficients for
C              the species' viscosity; dimension at least NO, usually 4
C
C  OUTPUT
C  If K > 1:
C  COFVIS    - Real vector of polynomial coefficients; dimension
C              at least NO, usually = 4
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT real (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST,
     4                NCST, NXL, NR, NWRK, K3
      DIMENSION RMCWRK(*), COFVIS(*)
C
      NK = IABS(K)
      IF (K .GT. 0) THEN
C
C GET the data
C
         DO 200 N = 1, NO
            COFVIS(N) = RMCWRK(NETA + (NK-1)*NO + N - 1 )
200      CONTINUE
      ELSE IF (K .LT. 0) THEN
C
C PUT the data
C
         DO 400 N = 1, NO
            RMCWRK(NETA + (NK-1)*NO + N - 1) = COFVIS(N)
400      CONTINUE
         ENDIF
C
C     end of SUBROUTINE MCCVEX
      RETURN
      END

