C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       GETFNM.FOR
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         8/5/89
C
C    PURPOSE:      OBTAINS A FILE NAME FROM USER
C
C    CALL LINE:    CALL GETFNM(FNAME,FNML,PROMPT,PROMTL)
C
C    INPUTS:       PROMPT = STRING TO POMPT USER WITH (C*60)
C                  PROMTL = LENGTH OF PROMPT STRING (I)
C
C    OUTPUTS:      FNAME = FILE NAME (C*60)
C                  FNML = FILE NAME LENGTH (I)
C*************************************************************************
C
      SUBROUTINE GETFNM(FNAME,FNML,PROMPT,PROMTL)
C
      IMPLICIT NONE
C
      INTEGER FNML, PROMTL
      CHARACTER*60 FNAME, PROMPT
C
C                       Local variables
C
      INTEGER STDOUT, STDIN, I
      DATA STDOUT/6/, STDIN/5/
C
C                       Prompt user for file name
C
      WRITE(STDOUT,'(T2,A,$)') PROMPT(1:PROMTL)
      READ(STDIN,'(A)') FNAME
C
C                        Figure out input file name length
      DO 10 I=1,20,1
        IF(FNAME(I:I).EQ.' ') GOTO 20
 10   CONTINUE
C
 20   CONTINUE
      FNML = I-1
C
C
 999  RETURN
      END
