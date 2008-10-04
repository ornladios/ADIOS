C ckerr.f
C ---------------------------------------------------------------------------
C $Header: /spin/home/rsankar/S3D-Repository/source/f77_files/ckerr.f,v 1.1.1.1 2005/01/11 00:27:29 rsankar Exp $
C ---------------------------------------------------------------------------
C
C FORTRAN callable implementation of Chemkin error handling.
C ---------------------------------------------------------------------------

C
C SUBROUTINE ERINIT( MAXERR )
C
	SUBROUTINE ERINIT( MAXERS )
C
C Dimension and initialize ERSAVE variables for use in Chemkin
C Error Handling Routines.
C
C  Input: 
C  MAXERS - User specification of maximum number of errors recorded
C           per program run
C
C  Parameters:
C  MAXERR - Maximum dimension for number of recordable errors
C  MAXPRM - Maximum number of parameters per error statement
C  MAXPCH - Maximum number of characters per parameter string
C           + 1 delimiter
C  MAXPST - Maximum length of the composit parameter string PARMST
C           (NOTE:  This must be less than 32767)
C  MAXNPR - Maximum number of parameter substrings stored in PARMST
C
C  Common block /ERSAV1/:
C  MAXSAV - User specified maximum error number for keeping (=MAXERS)
C  IDMOD  - Module ID number, dimensioned MAXERR
C  IDERR  - Error ID number, dimensioned MAXERR
C  ISEVER - Integer indicating severity of error;
C           3=ERROR, 2=WARN, 1=INFO.
C  NERPAR - Number of parameters recorded in PARMST for each error.
C  INDEXP - Indices in the PARMST of the $ delimitors for parameters
C  LASTCH - Last character used in PARMST.
C  NUMPAR - Total number of parameter strings stored in PARMST.
C  NUMERR - Number of errors found.
C
C  Common block /ERSAV2/:
C  PARMST - Character string containing Error-message parameters
C           delimited by '$' but not '$$'.
C
      INTEGER  MAXERS
	INTEGER  MAXSAV, MAXERR, MAXPRM, MAXPCH, MAXPST, MAXNPR
      PARAMETER (MAXPRM = 10, MAXPCH = 21, MAXERR=100)
      PARAMETER (MAXNPR = MAXPRM * MAXERR, MAXPST = MAXPCH * MAXNPR )
C
      INTEGER  IDMOD(MAXERR), IDERR(MAXERR), ISEVER(MAXERR), 
     &         NERPAR(MAXERR), INDEXP(MAXERR*MAXPRM)
      INTEGER  NUMERR, LASTCH, NUMPAR
      CHARACTER PARMST*(MAXPST)
      COMMON /ERSAV1/ MAXSAV,
     &                 IDMOD, IDERR, ISEVER, 
     &                 NERPAR, INDEXP,
     &                 LASTCH, NUMPAR, NUMERR 
      COMMON /ERSAV2/ PARMST
      DATA IDMOD/MAXERR*0/, IDERR/MAXERR*0/, ISEVER/MAXERR*0/, 
     &     INDEXP/MAXNPR*0/,
     &     NERPAR/MAXERR*0/, LASTCH/0/, NUMPAR/0/, NUMERR/0/, 
     &     PARMST/' '/ 
      MAXSAV = MAXERS
      RETURN
	END
C
C END SUBROUTINE ERINIT
C
C----------------------------------------------------------------------C
C
C SUBROUTINE ERSET( LOCMOD, LOCERR, LOCSEV, MSGPRM )
C
	SUBROUTINE ERSET( MODULE, LOCERR, LOCSEV, MSGPRM )
C
C Set the error information for the next error recorded in
C the ERSAV1 common block for later extraction by Chemkin
C Error Handling Routines.
C
C  Parameters:
C  MAXERR - Maximum dimension for number of recordable errors
C  MAXPRM - Maximum number of parameters per error statement
C  MAXPCH - Maximum number of characters per parameter string
C           + 1 delimiter
C  MAXPST - Maximum length of the composit parameter string PARMST
C           (NOTE:  This must be less than 32767)
C  MAXNPR - Maximum number of parameter substrings stored in PARMST
C
C  Common block /ERSAV1/:
C  MAXSAV - User specified maximum error number for keeping
C  IDMOD  - Module ID number, dimensioned MAXERR
C  IDERR  - Error ID number, dimensioned MAXERR
C  ISEVER - Integer indicating severity of error;
C           3=ERROR, 2=WARN, 1=INFO.
C  NERPAR - Number of parameters recorded in PARMST for each error.
C  INDEXP - Indices in the PARMST of the $ delimitors for parameters
C  LASTCH - Last character used in PARMST.
C  NUMPAR - Total number of parameter strings stored in PARMST.
C  NUMERR - Number of errors found.
C
C  Common block /ERSAV2/:
C  PARMST - Character string containing Error-message parameters
C           delimited by '$' but not '$$'.
C
	INTEGER  MAXSAV, MAXERR, MAXPRM, MAXPCH, MAXPST, MAXNPR
      PARAMETER (MAXPRM = 10, MAXPCH = 21, MAXERR=100)
      PARAMETER (MAXNPR = MAXPRM * MAXERR, MAXPST = MAXPCH * MAXNPR )
C
      INTEGER  IDMOD(MAXERR), IDERR(MAXERR), ISEVER(MAXERR), 
     &         NERPAR(MAXERR), INDEXP(MAXERR*MAXPRM)
      INTEGER  NUMERR, LASTCH, NUMPAR
      CHARACTER PARMST*(MAXPST)
      COMMON /ERSAV1/ MAXSAV,
     &                 IDMOD, IDERR, ISEVER, 
     &                 NERPAR, INDEXP,
     &                 LASTCH, NUMPAR, NUMERR 
      COMMON /ERSAV2/ PARMST
C
      INTEGER LOCERR, LOCSEV
      CHARACTER*(*) MSGPRM
      CHARACTER MODULE*(*)
C  external cklib routines
      INTEGER CKLSCH, CKFRCH
      EXTERNAL CKLSCH, CKFRCH
C  local variables
      INTEGER ISTBEG, ISTEND, LENSTR, NITEM, NLAST, N, NCKMOD, LOCMOD
      PARAMETER (NCKMOD=33)
      CHARACTER*20 CKMODS(NCKMOD)
      DATA CKMODS /'aurora',
     &             'block',
     &             'cdassl',
     &             'ckinterp',
     &             'cklib',
     &             'creslaf',
     &             'dasac',
     &             'dassl',
     &             'eqlib',
     &             'mach',
     &             'math',
     &             'oppdif',
     &             'ovend',
     &             'plug',
     &             'premix',
     &             'refine',
     &             'senkin',
     &             'shock',
     &             'skinterp',
     &             'sklib',
     &             'spin',
     &             'spl',
     &             'stanlib',
     &             'surftherm',
     &             'tools',
     &             'tranfit',
     &             'tranlib',
     &             'twafer',
     &             'twopnt',
     &             'vode',
     &             'xerror',
     &             'xlinpk' ,
     &             'UNRECOGNIZED'/
C
C  increment error count and store local module & error numbers
C
      NUMERR = NUMERR + 1
C
C    find the module id number
C
      CALL CKCOMP(MODULE,CKMODS,NCKMOD,LOCMOD)
      IF (LOCMOD.LE.0) LOCMOD = NCKMOD
C      
      IDMOD(NUMERR) = LOCMOD
      IDERR(NUMERR) = LOCERR
      ISEVER(NUMERR) = LOCSEV
C 
C  parse parameter string to find any parameters
C    MSGPRM should look like 'parm1$parm2$parm3...'
C    PARMST will look like '$er1parm1$er1parm2$er1parm3$er2parm1...'
C  
      ISTBEG = CKFRCH(MSGPRM)
      ISTEND = CKLSCH(MSGPRM)
      LENSTR = ISTEND - ISTBEG + 1
C
      IF (LENSTR .LE. 0) THEN
         NERPAR(NUMERR) = 0
      ELSE
         NUMPAR = NUMPAR+1
         NITEM = 1
         NLAST = 0
         LASTCH = LASTCH + 1
         PARMST(LASTCH:LASTCH)='$'
         INDEXP(NUMPAR) = LASTCH
         DO 10 N = ISTBEG, ISTEND
            IF (MSGPRM(N:N).EQ.'$'.AND.MSGPRM(N:N+1).NE.'$$'
     &                            .AND.MSGPRM(N-1:N).NE.'$$') THEN
               PARMST(LASTCH+1:LASTCH+N)=MSGPRM(NLAST+1:NLAST+N)
               NLAST = N
               LASTCH = LASTCH + N
               NUMPAR = NUMPAR+1
               INDEXP(NUMPAR) = LASTCH
               NITEM = NITEM + 1
            ENDIF
10       CONTINUE
         NERPAR(NUMERR) = NITEM
      ENDIF
C
      RETURN
	END
C
C END SUBROUTINE ERRSET
C
C----------------------------------------------------------------------C
C
C SUBROUTINE ERCHEK()
C
	INTEGER FUNCTION ERCHEK()
C
C Check errors recorded by calls to ERSET, and return 
C highest integer severity value.
C
C  Parameters:
C  MAXERR - Maximum dimension for number of recordable errors
C  MAXPRM - Maximum number of parameters per error statement
C  MAXPCH - Maximum number of characters per parameter string
C           + 1 delimiter
C  MAXPST - Maximum length of the composit parameter string PARMST
C           (NOTE:  This must be less than 32767)
C  MAXNPR - Maximum number of parameter substrings stored in PARMST
C
C  Common block /ERSAV1/:
C  MAXSAV - User specified maximum error number for keeping
C  IDMOD  - Module ID number, dimensioned MAXERR
C  IDERR  - Error ID number, dimensioned MAXERR
C  ISEVER - Integer indicating severity of error;
C           3=ERROR, 2=WARN, 1=INFO.
C  NERPAR - Number of parameters recorded in PARMST for each error.
C  INDEXP - Indices in the PARMST of the $ delimitors for parameters
C  LASTCH - Last character used in PARMST.
C  NUMPAR - Total number of parameter strings stored in PARMST.
C  NUMERR - Number of errors found.
C
C  Common block /ERSAV2/:
C  PARMST - Character string containing Error-message parameters
C           delimited by '$' but not '$$'.
C
	INTEGER  MAXSAV, MAXERR, MAXPRM, MAXPCH, MAXPST, MAXNPR
      PARAMETER (MAXPRM = 10, MAXPCH = 21, MAXERR=100)
      PARAMETER (MAXNPR = MAXPRM * MAXERR, MAXPST = MAXPCH * MAXNPR )
C
      INTEGER  IDMOD(MAXERR), IDERR(MAXERR), ISEVER(MAXERR), 
     &         NERPAR(MAXERR), INDEXP(MAXERR*MAXPRM)
      INTEGER  NUMERR, LASTCH, NUMPAR
      CHARACTER PARMST*(MAXPST)
      COMMON /ERSAV1/ MAXSAV,
     &                 IDMOD, IDERR, ISEVER, 
     &                 NERPAR, INDEXP,
     &                 LASTCH, NUMPAR, NUMERR 
      COMMON /ERSAV2/ PARMST
C  local variables
      INTEGER N, IERROR
C      
      IERROR = 0
      DO 10 N = 1, NUMERR
         IERROR = MAX(IERROR,ISEVER(N))
10    CONTINUE
      ERCHEK = IERROR
      RETURN
	END
C
C END FUNCTION ERCHEK
C
C----------------------------------------------------------------------C
C
C SUBROUTINE ERGET(STRING)
C
	SUBROUTINE ERGET( NUMFND, ERRMSG )
C
C Get the error messages recorded by calls to ERSET in one
C string, delimited by '$' characters.
C
C  Common block /ERSAV1/:
C  MAXERR - Maximum number of recordable errors per program run
C  MAXPRM - Maximum number of parameters per error statement
C  MAXPCH - Maximum number of characters per parameter string
C           + 1 delimiter
C  IDMOD  - Module ID number, dimensioned MAXERR
C  IDERR  - Error ID number, dimensioned MAXERR
C  ISEVER - Integer indicating severity of error;
C           3=ERROR, 2=WARN, 1=INFO.
C  NERPAR - Number of parameters recorded in PARMST for each error.
C  INDEXP - Indices in the PARMST of the $ delimitors for parameters
C  LASTCH - Last character used in PARMST.
C  NUMPAR - Total number of parameter strings stored in PARMST.
C  NUMERR - Number of errors found.
C
C  Common block /ERSAV2/:
C  PARMST - Character string containing Error-message parameters
C           delimited by '$' but not '$$'.
C
C  Return value:
C    0 - no errors were recorded
C    1 - only information messages were recorded
C    2 - only warnings or information messages were recorded
C    3 - fatal errors were found.
C
      INTEGER  MAXSAV
	INTEGER  MAXERR, MAXPRM, MAXPCH, MAXPST, MAXNPR
      PARAMETER (MAXPRM = 10, MAXPCH = 21, MAXERR=100)
      PARAMETER (MAXNPR = MAXPRM * MAXERR, MAXPST = MAXPCH * MAXNPR )
      INTEGER  IDMOD(MAXERR), IDERR(MAXERR), ISEVER(MAXERR), 
     &         NERPAR(MAXERR), INDEXP(MAXERR*MAXPRM)
      INTEGER  NUMERR, LASTCH, NUMPAR
      CHARACTER PARMST*(MAXPST)
      COMMON /ERSAV1/ MAXSAV,
     &                 IDMOD, IDERR, ISEVER, 
     &                 NERPAR, INDEXP,
     &                 LASTCH, NUMPAR, NUMERR 
      COMMON /ERSAV2/ PARMST
      CHARACTER*(*) ERRMSG
      ERRMSG = ' '
      NUMFND = NUMERR
C  local variables
C TODO:  call an xml? parser to get error messages that
C        correspond to error numbers and fill in parameters, 
C        then return a string.
      RETURN
	END
C
C END SUBROUTINE ERGET
C
C----------------------------------------------------------------------C
C
C SUBROUTINE ERTERM()
C
	SUBROUTINE ERTERM()
      RETURN
	END
C
C END SUBROUTINE ERTERM
C
C ---------------------------------------------------------------------------
C end ckerr.f

