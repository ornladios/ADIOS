!C-----------------------------------------------------------------------
!C  IMSL Name:  B22DR/DB22DR (Single/Double precision version)
!C
!C  Computer:   CONVEX/DOUBLE
!C
!C  Revised:    August 7, 1986
!C
!C  Purpose:    Evaluate the derivatives of a two-dimensional tensor
!C              product spline, given its tensor product B-spline
!C              representation.
!C
!C  Usage:      B22DR(IXDER, IYDER, X, Y, KXORD, KYORD, XKNOT, YKNOT,
!C                    NXCOEF, NYCOEF, BSCOEF, WK)
!C
!C  Arguments:  (See BS2DR)
!C
!C  Chapter:    MATH/LIBRARY Interpolation and Approximation
!C
!C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
!C
!C  Warranty:   IMSL warrants only that IMSL testing has been applied
!C              to this code.  No other warranty, expressed or implied,
!C              is applicable.
!C
!C-----------------------------------------------------------------------
!C
real (kind=8) function my_DB22DR (IXDER, IYDER, X, Y, KXORD, &
     KYORD, XKNOT, YKNOT, NXCOEF, NYCOEF, BSCOEF)
  implicit none
  !                                  SPECIFICATIONS FOR ARGUMENTS
  integer  ::  IXDER, IYDER, KXORD, KYORD, NXCOEF, NYCOEF
  real (kind=8)  :: X, Y, XKNOT(*), YKNOT(*), BSCOEF(NXCOEF,*)
  !                                 SPECIFICATIONS FOR LOCAL VARIABLES
  integer    IAJ, IDL, IDR, IYBSCF, J, LEFTY, MFLAG, MXKORD
  real (kind=8) ::  VALUE
  !                                  SPECIFICATIONS FOR INTRINSICS
  INTRINSIC  MAX0
  !                                  SPECIFICATIONS FOR FUNCTIONS
  real (kind=8) , external ::   my_DB3DER
  real (kind=8) :: WK1(3+10), WK2(3+10), WK3(3+10), WK4(3+10)
  ! the size of wk1, wk2, wk3, wk4 should be larger then kxord and kyord.

  save WK1, WK2, WK3,WK4

  VALUE = 0.0D0

  MXKORD = MAX0(KXORD,KYORD)
!  IAJ = 1
!  IDL = IAJ + MXKORD
!  IDR = IDL + MXKORD
!  IYBSCF = IDR + MXKORD

  call my_DB4DER (YKNOT, NYCOEF+KYORD, Y, LEFTY, MFLAG)
  if (MFLAG .ne. 0) GO TO 9000
!C                                  Interpolate to compute B-spline
!C                                  coefs.
  do  J=1, KYORD
     WK4(J) = my_DB3DER(IXDER,X,KXORD,XKNOT,NXCOEF,BSCOEF(1, &
          LEFTY-KYORD+J),WK1,WK2,WK3)
!     if (N1RTY(0) .ne. 0) GO TO 9000
  enddo
!C                                  Interpolate to find the value
!C                                  at (X,Y).
  VALUE = my_DB3DER(IYDER,Y,KYORD,YKNOT(LEFTY-KYORD+1),KYORD, &
       WK4,WK1,WK2,WK3)
!C

9000 continue
  my_DB22DR = VALUE

end function MY_DB22DR



!C-----------------------------------------------------------------------
!C  IMSL Name:  B3DER/DB3DER (Single/Double precision version)
!C
!C  Computer:   CONVEX/DOUBLE
!C
!C  Revised:    August 7, 1986
!C
!C  Purpose:    Calculate the value of the I-th derivative of a spline
!C              given its B-spline representation.
!C
!C  Usage:      B3DER(IDERIV, X, KORDER, XKNOT, NCOEF, BSCOEF, AJ,
!C                    DL, DR)
!C
!C  Arguments:  (See BSDER)
!C
!C  Chapter:    MATH/LIBRARY Interpolation and Approximation
!C
!C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
!C
!C  Warranty:   IMSL warrants only that IMSL testing has been applied
!C              to this code.  No other warranty, expressed or implied,
!C              is applicable.
!C
!C-----------------------------------------------------------------------
!C
real (kind=8)  function my_DB3DER (IDERIV, X, KORDER, XKNOT, &
     NCOEF, BSCOEF, AJ, DL, DR)
  implicit none
  !C                                  SPECIFICATIONS FOR ARGUMENTS
  integer  ::  IDERIV, KORDER, NCOEF
  real (kind=8) ::  X, XKNOT(*), BSCOEF(*), AJ(*), DL(*), DR(*)
  !C                                  SPECIFICATIONS FOR LOCAL VARIABLES
  integer    I, ILO, J, JCMAX, JCMIN, JJ, MFLAG, NKNOT
  real (kind=8) ::  VALUE
  !C                                  SPECIFICATIONS FOR INTRINSICS
  !C     INTRINSIC  DBLE
  !C                                  SPECIFICATIONS FOR SUBROUTINES
  !external   DCOPY, DSET, DB4DER
  
  VALUE = 0.0D0
  if (IDERIV .ge. KORDER) GO TO 9000
  !C                                  Find I such that 1 .LE. I .LT. NKNOT
  !C                                  and XKNOT(I) .LT. XKNOT(I+1) and
  !C                                  XKNOT(I) .LE. X .LT. XKNOT(I+1) .
  !C                                  If no such I can be found, X lies
  !C                                  outside the support of the spline
  !C                                  F and VALUE = 0. (The asymmetry in
  !C                                  this choice of I makes F
  !C                                  right continuous)
  NKNOT = NCOEF + KORDER
  call my_DB4DER (XKNOT, NKNOT, X, I, MFLAG)
  if (MFLAG .ne. 0) GO TO 9000
  if (KORDER .le. 1) then
     VALUE = BSCOEF(I)
     GO TO 9000
  end if
  !C                                  Store the KORDER B-spline
  !C                                  coefficients relevant for the
  !C                                  XKNOT interval
  !C                                  (XKNOT(I),XKNOT(I+1)) in
  !C                                  AJ(1),...,AJ(KORDER) and compute
  !C                                  DL, DR. Set any of the AJ not
  !C                                  obtainable. From input to zero.
  !C                                  Set any XKNOTS not obtainable
  !C                                  equal to XKNOT(1) to XKNOT(NKNOT)
  !C                                  appropriately.
  JCMIN = 1
  if (I .lt. KORDER) then
     JCMIN = KORDER - I + 1
     do   J=1, I
        DL(J) = X - XKNOT(I+1-J)
     enddo
     AJ(1:KORDER-I) = 0D0        !     call DSET (KORDER-I, 0.0D0, AJ, 1)
     DL(I:korder)=DL(I)          !     call DSET (KORDER-I, DL(I), DL(I), 1)


  else
     do  J=1, KORDER - 1
        DL(J) = X - XKNOT(I+1-J)
     enddo
  end if
  !C
  JCMAX = KORDER
  if (NCOEF .lt. I) then
     JCMAX = NKNOT - I
     do  J=1, JCMAX
        DR(J) = XKNOT(I+J) - X
     enddo
     AJ(JCMAX+1:KORDER)=0D0        !call DSET (KORDER-JCMAX, 0.0D0, AJ(JCMAX+1), 1)
     DR(JCMAX: KORDER-1)=DR(JCMAX) !call DSET (KORDER-JCMAX, DR(JCMAX), DR(JCMAX), 1)
  else
     do  J=1, KORDER - 1
        DR(J) = XKNOT(I+J) - X
     enddo
  end if
  !call DCOPY (JCMAX-JCMIN+1, BSCOEF(I-KORDER+JCMIN), 1, AJ(JCMIN),1)
  AJ(JCMIN: JCMAX) = BSCOEF(I-KORDER+JCMIN: I-KORDER + JCMAX)


  !C                                  Difference the coefficients IDERIV
  !C                                  times.
  do J=1, IDERIV
     ILO = KORDER - J
     do   JJ=1, KORDER - J
        AJ(JJ) = ((AJ(JJ+1)-AJ(JJ))/(DL(ILO)+DR(JJ)))* dble(KORDER-J)
        ILO = ILO - 1
     enddo
  enddo
  !C                                  Compute value at X in
  !C                                  (XKNOT(I),XKNOT(I+1)) of IDERIV-th
  !C                                  derivative, given its relevant
  !C                                  B-spline coeffs in
  !C                                  AJ(1),...,AJ(KORDER-IDERIV).
  do  J=IDERIV + 1, KORDER - 1
     ILO = KORDER - J
     do   JJ=1, KORDER - J
        AJ(JJ) = (AJ(JJ+1)*DL(ILO)+AJ(JJ)*DR(JJ))/(DL(ILO)+DR(JJ))
        ILO = ILO - 1
     enddo
  enddo
  
  VALUE = AJ(1)
  
9000 continue
  my_DB3DER = VALUE
  return
end function MY_DB3DER



!C-----------------------------------------------------------------------
!C  IMSL Name:  B4DER/DB4DER (Single/Double precision version)
!C
!C  Computer:   CONVEX/DOUBLE
!C
!C  Revised:    August 7, 1986
!C
!C  Purpose:    Compute LEFT = MAX (I, 1 .LE. NKNOT. .AND. XKNOT(I) .LE.
!C              X
!C
!C  Usage:      CALL B4DER (XKNOT, NKNOT, X, LEFT, MFLAG)
!C
!C  Arguments:
!C     XKNOT  - The knot sequence.  (Input)
!C     NKNOT  - Number of knots.  (Input)
!C     X      - The point whose location in XKNOT is to be found.
!C              (Input)
!C     LEFT   - Integer whose value is given below.  (Output)
!C     MFLAG  - Flag defined below.  (Output)
!C               LEFT   MFLAG
!C                 1     -1      IF                    X .LT.  XKNOT(1)
!C                 I      0      IF   XKNOT(I)    .LE. X .LT. XKNOT(I+1)
!C              NKNOT     1      IF  XKNOT(NKNOT) .LE. X
!C              In particular, MFLAG = 0 is the usual case.  MFLAG .NE. 0
!C              indicates that X lies outside the half open interval
!C              XKNOT(1) .LE. X .LT. XKNOT(NKNOT). The asymmetric
!C              treatment of the interval is due to the decision to make
!C              all PP functions continuous from the right.
!C
!C  Chapter:    MATH/LIBRARY Interpolation and Approximation
!C
!C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
!C
!C  Warranty:   IMSL warrants only that IMSL testing has been applied
!C              to this code.  No other warranty, expressed or implied,
!C              is applicable.
!C
!C-----------------------------------------------------------------------
!C
SUBROUTINE my_DB4DER (XKNOT, NKNOT, X, LEFT, MFLAG)
  !C                                  SPECIFICATIONS FOR ARGUMENTS
  INTEGER ::   NKNOT, LEFT, MFLAG
  real (kind=8) ::  X, XKNOT(*)
  !C                                  SPECIFICATIONS FOR LOCAL VARIABLES
  INTEGER ::   IHI, ISTEP, MIDDLE
  !C                                  SPECIFICATIONS FOR SAVE VARIABLES
  INTEGER ::   ILO
  SAVE       ILO
  !C
  DATA ILO/1/
  !C
  IHI = ILO + 1
  IF (IHI .GE. NKNOT) THEN
     IF (X .GE. XKNOT(NKNOT)) THEN
        MFLAG = 1
        LEFT = NKNOT
        GO TO 9000
     ELSE IF (NKNOT .LE. 1) THEN
        MFLAG = -1
        LEFT = 1
        GO TO 9000
     END IF
     ILO = NKNOT - 1
     IHI = NKNOT
  END IF
  !C
  IF (X .LT. XKNOT(IHI)) THEN
     IF (X .GE. XKNOT(ILO)) THEN
        MFLAG = 0
        LEFT = ILO
        GO TO 9000
     END IF
     !C                                  Now X .LT. XKNOT(ILO) . Decrease ILO
     !C                                  to capture X .
     ISTEP = 1
10   CONTINUE
     IHI = ILO
     ILO = IHI - ISTEP
     IF (ILO .GT. 1) THEN
        IF (X .GE. XKNOT(ILO)) GO TO 30
        ISTEP = ISTEP*2
        GO TO 10
     END IF
     ILO = 1
     IF (X .LT. XKNOT(1)) THEN
        MFLAG = -1
        LEFT = 1
        GO TO 9000
     END IF
     GO TO 30
  END IF
!  C                                  Now X .GE. XKNOT(IHI) . Increase IHI
!  C                                  to capture X .
  ISTEP = 1
20 CONTINUE
  ILO = IHI
  IHI = ILO + ISTEP
  IF (IHI .LT. NKNOT) THEN
     IF (X .LT. XKNOT(IHI)) GO TO 30
     ISTEP = ISTEP*2
     GO TO 20
  END IF
  IF (X .GE. XKNOT(NKNOT)) THEN
     MFLAG = 1
     LEFT = NKNOT
     GO TO 9000
  END IF
  IHI = NKNOT
!  C                                  Now XKNOT(ILO) .LE. X .LT.
!  C                                  XKNOT(IHI) . Narrow the inteval.
30 CONTINUE
  MIDDLE = (ILO+IHI)/2
  IF (MIDDLE .EQ. ILO) THEN
     MFLAG = 0
     LEFT = ILO
     GO TO 9000
  END IF
!  C                                  It is assumed that MIDDLE = ILO in
!  C                                  case IHI = ILO+1 .
  IF (X .LT. XKNOT(MIDDLE)) THEN
     IHI = MIDDLE
  ELSE
     ILO = MIDDLE
  END IF
  GO TO 30
9000 CONTINUE
  RETURN
END SUBROUTINE MY_DB4DER


!---------------

!!$C-----------------------------------------------------------------------
!!$C  IMSL Name:  PPDER/DPPDER (Single/Double precision version)
!!$C
!!$C  Computer:   CONVEX/DOUBLE
!!$C
!!$C  Revised:    August 11, 1986
!!$C
!!$C  Purpose:    Evaluate the derivative of a piecewise polynomial.
!!$C
!!$C  Usage:      PPDER(IDERIV, X, KORDER, NINTV, BREAK, PPCOEF)
!!$C
!!$C  Arguments:
!!$C     IDERIV - Order of the derivative to be evaluated.  (Input)
!!$C              In particular, IDERIV = 0 returns the value of the
!!$C              polynomial.
!!$C     X      - Point at which the polynomial is to be evaluated.
!!$C              (Input)
!!$C     KORDER - Order of the polynomial.  (Input)
!!$C     NINTV  - Number of polynomial pieces.  (Input)
!!$C     BREAK  - Array of length NINTV+1 containing the breakpoints of
!!$C              the piecewise polynomial representation.  (Input)
!!$C              BREAK must be strictly increasing.
!!$C     PPCOEF - Array of size KORDER*NINTV containing the
!!$C              local coefficients of the piecewise polynomial pieces.
!!$C              (Input)
!!$C              PPCOEF is treated internally as a matrix of size
!!$C              KORDER by NINTV.
!!$C     PPDER  - Value of the IDERIV-th derivative of the piecewise
!!$C              polynomial at X.  (Output)
!!$C
!!$C  GAMS:       E3; K6
!!$C
!!$C  Chapter:    MATH/LIBRARY Interpolation and Approximation
!!$C
!!$C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
!!$C
!!$C  Warranty:   IMSL warrants only that IMSL testing has been applied
!!$C              to this code.  No other warranty, expressed or implied,
!!$C              is applicable.
!!$C
!!$C-----------------------------------------------------------------------
!!$C
real (kind=8)  FUNCTION my_DPPDER (IDERIV, X, KORDER, NINTV,BREAK, PPCOEF)
!C                                  SPECIFICATIONS FOR ARGUMENTS
   INTEGER    IDERIV, KORDER, NINTV
   real (kind=8) :: X, BREAK(*), PPCOEF(KORDER,*)
!C                                  SPECIFICATIONS FOR LOCAL VARIABLES
   INTEGER    J, LEFT
   real (kind=8) :: FMM, H, VALUE
   !C                                  SPECIFICATIONS FOR SUBROUTINES
   !   EXTERNAL   E1MES, E1POP, E1PSH, E1STI, DP3DER
   !C                                  SPECIFICATIONS FOR FUNCTIONS
   !      EXTERNAL   N1RTY
   !      INTEGER    N1RTY
   !C
   !      CALL E1PSH ('DPPDER ')
   !C                                  IN CASE OF ERRORS
   VALUE = 0.0D0
   !C                                  Check argument NINTV
   !      IF (NINTV .LT. 1) THEN
   !         CALL E1STI (1, NINTV)
   !         CALL E1MES (5, 1, 'The number of intervals must be '//
   !     &               'at least 1 while NINTV = %(I1) is given. ')
   !      END IF
   !C                                  Check argument IDERIV
   !      IF (IDERIV .LT. 0) THEN
   !         CALL E1STI (1, IDERIV)
   !         CALL E1MES (5, 2, 'The order of the derivative must '//
   !     &               'be positive while IDERIV = %(I1) is given.')
   !      END IF
   !C                                  Check argument KORDER
   !      IF (KORDER .LE. 0) THEN
   !         CALL E1STI (1, KORDER)
   !         CALL E1MES (5, 3, 'The order of the interpolating '//
   !     &               'polynomial must be positive while KORDER = '//
   !     &               '%(I1) is given.')
   !      END IF
   !C                                  Check for errors
   !      IF (N1RTY(0) .NE. 0) GO TO 9000
   !C                                  Derivatives of order KORDER or
   !C                                  higher are identically zero
   IF (IDERIV .GE. KORDER) GO TO 9000
   !C                                  Find index I of largest breakpoint
   !C                                  to the left of X
   CALL my_DP3DER (KORDER, NINTV, BREAK, X, LEFT)
      !C                                  Evaluate jderiv-th derivative of
   !C                                  I-th polynomial piece at X
   FMM = KORDER - IDERIV
   H = X - BREAK(LEFT)
   DO  J=KORDER, IDERIV + 1, -1
      VALUE = (VALUE/FMM)*H + PPCOEF(J,LEFT)
      FMM = FMM - 1.0D0
   enddo
   !C
9000 CONTINUE
   my_DPPDER = VALUE
!   CALL E1POP ('DPPDER ')
   RETURN
 END
    
!!$C-----------------------------------------------------------------------
!!$C  IMSL Name:  P3DER/DP3DER (Single/Double precision version)
!!$C
!!$C  Computer:   CONVEX/DOUBLE
!!$C
!!$C  Revised:    August 11, 1986
!!$C
!!$C  Purpose:    Compute MAX (I, 1 .LE. NINTV. .AND. BREAK(I) .LE. X)
!!$C
!!$C  Usage:      CALL P3DER (KORD, NINTV, BREAK, X, LEFT)
!!$C
!!$C  Arguments:
!!$C     KORD   - Order of the polynomial.  (Input)
!!$C     NINTV  - Number of polynomial pieces.  (Input)
!!$C     BREAK  - Vector of length NINTV+1 containing the breakpoints
!!$C              of the piecewise polynomial representation.  (Input)
!!$C     X      - The point whose location in BREAK is to be found.
!!$C     LEFT   - Integer whose value is
!!$C                LEFT
!!$C                  1      IF                       X .LT.  BREAK(1)
!!$C                  I      IF         BREAK(I).LE. X .LT. BREAK(I+1)
!!$C                 NINTV   IF                    BREAK(NINTV) .LE. X
!!$C              The asymmetric treatment of the interval is due to the
!!$C              decision to make all PP functions continuous from the
!!$C              right.  (Output)
!!$C
!!$C  Remark:
!!$C     This routine is based in INTERV in de Boor, p92-93.
!!$C
!!$C  Chapter:    MATH/LIBRARY Interpolation and Approximation
!!$C
!!$C  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
!!$C
!!$C  Warranty:   IMSL warrants only that IMSL testing has been applied
!!$C              to this code.  No other warranty, expressed or implied,
!!$C              is applicable.
!!$C
!!$C-----------------------------------------------------------------------
!C
      SUBROUTINE my_DP3DER (KORD, NINTV, BREAK, X, LEFT)
!C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    KORD, NINTV, LEFT
      real (kind=8) ::  X, BREAK(*)
!C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IHI, ISTEP, MIDDLE
!C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    ILO
      SAVE       ILO
!C
      DATA ILO/1/
!C
      IHI = ILO + 1
      IF (IHI .GE. NINTV) THEN
         IF (X .GE. BREAK(NINTV)) THEN
            LEFT = NINTV
            GO TO 9000
         ELSE IF (NINTV .LE. 1) THEN
            LEFT = 1
            GO TO 9000
         END IF
         ILO = NINTV - 1
         IHI = NINTV
      END IF
!C
      IF (X .LT. BREAK(IHI)) THEN
         IF (X .GE. BREAK(ILO)) THEN
            LEFT = ILO
            GO TO 9000
         END IF
!C                                  Now X .LT. BREAK(ILO) - decrease ILO
!C                                  to capture X
         ISTEP = 1
   10    CONTINUE
         IHI = ILO
         ILO = IHI - ISTEP
         IF (ILO .GT. 1) THEN
            IF (X .GE. BREAK(ILO)) GO TO 30
            ISTEP = ISTEP*2
            GO TO 10
         END IF
         ILO = 1
         IF (X .LT. BREAK(1)) THEN
            LEFT = 1
            GO TO 9000
         END IF
         GO TO 30
      END IF
!C                                  Now X .GE. BREAK(IHI) - increase IHI
!C                                  to capture X
      ISTEP = 1
   20 CONTINUE
      ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .LT. NINTV) THEN
         IF (X .LT. BREAK(IHI)) GO TO 30
         ISTEP = ISTEP*2
         GO TO 20
      END IF
      IF (X .GE. BREAK(NINTV)) THEN
         LEFT = NINTV
         GO TO 9000
      END IF
      IHI = NINTV
!C                                  Now BREAK(ILO) .LE. X .LT.
!C                                  BREAK(IHI) - narrow the inteval
   30 CONTINUE
      MIDDLE = (ILO+IHI)/2
      IF (MIDDLE .EQ. ILO) THEN
         LEFT = ILO
         GO TO 9000
      END IF
!C                                  It is assumed that MIDDLE = ILO in
!C                                  case IHI = ILO+1
      IF (X .LT. BREAK(MIDDLE)) THEN
         IHI = MIDDLE
      ELSE
         ILO = MIDDLE
      END IF
      GO TO 30
 9000 CONTINUE
      RETURN
      END


!C-----------------------------------------------------------------------
!C  IMSL Name:  B2NAK/DB2NAK (Single/Double precision version)
!C
!C  Computer:   CONVEX/DOUBLE
!C
!C  Revised:    August 1, 1986
!C
!C  Purpose:    Compute the 'not-a-knot' spline knot sequence.
!C
!C  Usage:      CALL B2NAK (NDATA, XDATA, KORDER, XKNOT, XSRT, IWK)
!C
!C  Arguments:  (See BSNAK)
!C
!C  Chapter:    MATH/LIBRARY Interpolation and Approximation
!C
!C  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
!C
!C  Warranty:   IMSL warrants only that IMSL testing has been applied
!C              to this code.  No other warranty, expressed or implied,
!C              is applicable.
!C
!C-----------------------------------------------------------------------
!C
subroutine my_DB2NAK (NDATA, XDATA, KORDER, XKNOT)
  implicit none
!C                                  SPECIFICATIONS FOR ARGUMENTS
  integer ::  NDATA, KORDER
  real (kind=8) :: XDATA(*), XKNOT(*)
!C                                  SPECIFICATIONS FOR LOCAL VARIABLES
  integer    I, J
  real (kind=8) ::  EPS
!                                  SPECIFICATIONS FOR INTRINSICS
!     intrinsic  MOD

  double precision DMACH

!                                  Check argument KORDER
  if (KORDER .le. 1) then
     print *, 'The order of the spline must be at least 2 while KORDER =' &
          , korder, 'is given.'
  end if
  !                                  Check NDATA
  if (NDATA .lt. KORDER) then
     print *, 'The number of data points must be at least as large as the'
     print *, 'order of the spline while NDATA = ',NDATA, 'and KORDER = '&
          ,KORDER,' are given.'
  end if

  !   Check argument XDATA
  do 10  I=2, NDATA
     if (XDATA(I-1) .ge. XDATA(I)) then
        !   Check that XDATA values are distinct
        if (XDATA(I-1) .eq. XDATA(I)) then
           J = I - 1
           print *, 'Points in the data point abscissas array, XDATA, must be'
           print *, 'distinct, but XDATA(',J,') = XDATA(',I,') = ',XDATA(I)
           GO TO 9000
        else
           GO TO 20
        end if
     end if
10   continue
!C                                  data is already sorted.  Move
!C                                  XDATA to XSRT.
!     XSRT=XDATA
     GO TO 50
!C                                  Set initial permutation
20   print *, 'Sort XDATA. This is not implimented my imsl routine'
     goto 9000
!   20 do 30  I=1, NDATA
!         IWK(I) = I
!   30 continue
!C                                  Find sorting permutation
!      call DSVRGP (NDATA, XDATA, XSRT, IWK)
!C                                  Check the XDATA values are distinct
!      do 40  I=2, NDATA
!         if (XSRT(I-1) .eq. XSRT(I)) then
!            call E1STI (1, IWK(I-1))
!            call E1STI (2, IWK(I))
!            call E1STD (1, XSRT(I))
!            call E1MES (5, 5, 'Points in the data point abscissas '//
!     &                  'array, XDATA, must be distinct, '//
!     &                  'but XDATA(%(I1)) = XDATA(%(I2)) = %(D1).')
!            GO TO 9000
!         end if
!   40 continue
50   continue
!C                                  Move the last endpoint slightly to
!C                                  the right of the last data point.
     EPS = 100.0D0*DMACH(4)
60   if (XDATA(NDATA)+EPS .le. XDATA(NDATA)) then
        EPS = EPS*10.0D0
        GO TO 60
     end if
!C                                  The first and last knots have
!C                                  multiplicity KORDER.
     do 70  I=1, KORDER
        XKNOT(I) = XDATA(1)
        XKNOT(NDATA+I) = XDATA(NDATA) + EPS
70      continue
!C                                  The middle knots have multiplicty one
      if (mod(KORDER,2) .eq. 0) then
!C                                  KORDER even - knots at XDATA points.
!         call DCOPY (NDATA-KORDER, XDATA(KORDER/2+1), 1, XKNOT(KORDER+1), 1)
         XKNOT(KORDER+1:NDATA)=XDATA(KORDER/2+1:NDATA-KORDER/2)
      else
!C                                  KORDER odd - knots between XDATA
!C                                  points.
         do 80  I=KORDER + 1, NDATA
            XKNOT(I) = 0.5D0*(XDATA(I-KORDER/2-1)+XDATA(I-KORDER/2))
80          continue
      end if

9000     continue
      return
end subroutine 
