      SUBROUTINE GJR(A,NC,NR,N,MC,*,JC,V)                               MATH3501
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NR,NC),JC(1),V(2)                                     MATH3502
C              ---------------------------------------------------------MATH3503
C               JC IS THE PERMUTATION VECTOR                            MATH3504
C              KD IS THE OPTION KEY FOR DETERMINANT EVALUATION          MATH3505
C              KI IS THE OPTION KEY FOR MATRIX INVERSION                MATH3506
C              L IS THE COLUMN CONTROL FOR AX=B                         MATH3507
C              M IS THE COLUMN CONTOL FOR MATRIX INVERSION              MATH3508
C              ---------------------------------------------------------MATH3509
C              INITIALIZATION                                           MATH3510
C              ---------------------------------------------------------MATH3511
C                  -----------------------------------------------------MATH3512
C                  DUMMY CALL TO TURN OFF OVERFLOW SWITCH               MATH3513
C                  -----------------------------------------------------MATH3514
      CALL OVERFL(IFL)                                                  MATH3515
      IW=V(1)                                                           MATH3516
      M=1                                                               MATH3517
      S=1.                                                              MATH3518
      L=N+(MC-N)*(IW/4)                                                 MATH3519
      KD=2-MOD(IW/2,2)                                                  MATH3520
      IF(KD.EQ.1) V(2)=0.                                               MATH3521
      KI=2-MOD(IW,2)                                                    MATH3522
      GO TO (5,20),KI                                                   MATH3523
C              ---------------------------------------------------------MATH3524
C                  INITIALIZE JC FOR INVERSION                          MATH3525
C              ---------------------------------------------------------MATH3526
5     DO 10  I=1,N                                                      MATH3527
10    JC(I)=I                                                           MATH3528
C              ---------------------------------------------------------MATH3529
C              SEARCH FOR PIVOT ROW                                     MATH3530
C              ---------------------------------------------------------MATH3531
20    DO 91  I=1,N                                                      MATH3532
      GO TO (22,21),KI                                                  MATH3533
21    M=I                                                               MATH3534
22    IF (I.EQ.N)  GO TO 60                                             MATH3535
      X=-1.                                                             MATH3536
      DO 30  J=I,N                                                      MATH3537
      IF (X.GT.DABS(A(J,I)))  GO TO 30                                   MATH3538
      X=DABS(A(J,I))                                                     MATH3539
      K=J                                                               MATH3540
30    CONTINUE                                                          MATH3541
      IF(K.EQ.I) GO TO 60                                               MATH3542
      S=-S                                                              MATH3543
      V(1)=-V(1)                                                        MATH3544
      GO TO (35,40),KI                                                  MATH3545
 35   MU=JC(I)                                                          MATH3546
      JC(I)=JC(K)                                                       MATH3547
      JC(K)=MU                                                          MATH3548
C              ---------------------------------------------------------MATH3549
C              INTERCHANGE ROW I AND ROW K                              MATH3550
C              ---------------------------------------------------------MATH3551
40    DO 50  J=M,L                                                      MATH3552
      X=A(I,J)                                                          MATH3553
      A(I,J)=A(K,J)                                                     MATH3554
50    A(K,J)=X                                                          MATH3555
C              ---------------------------------------------------------MATH3556
C              TEST FOR SINGULARITY                                     MATH3557
C              ---------------------------------------------------------MATH3558
60    IF (DABS(A(I,I)).GT.0.)  GO TO 70                                  MATH3559
C              ---------------------------------------------------------MATH3560
C              MATRIX IS SINGULAR                                       MATH3561
C              ---------------------------------------------------------MATH3562
      IF(KD.EQ.1) V(1)=0.                                               MATH3563
      JC(1)=I-1                                                         MATH3564
      RETURN 1
70    GO TO (71,72),KD                                                  MATH3566
C              ---------------------------------------------------------MATH3567
C              COMPUTE THE DETERMINANT                                  MATH3568
C              ---------------------------------------------------------MATH3569
71    IF(A(I,I).LT.0.) S=-S                                             MATH3570
      V(2)=V(2)+DLOG(DABS(A(I,I)))                                       MATH3571
72    X=A(I,I)                                                          MATH3572
      A(I,I)=1.                                                         MATH3573
C              ---------------------------------------------------------MATH3574
C              REDUCTION OF THE I-TH ROW                                MATH3575
C              ---------------------------------------------------------MATH3576
      DO 80 J=M,L                                                       MATH3577
      A(I,J)=A(I,J)/X                                                   MATH3578
C              ---------------------------------------------------------MATH3579
C              TEST OVERFLOW SWITCH. IF ON                              MATH3580
C              RETURN NEGATIVE VALUE OF I IN JC(1)                      MATH3581
C              ---------------------------------------------------------MATH3582
      CALL OVERFL (IFL)                                                 MATH3583
      IF(IFL.EQ.1) GO TO 150                                            MATH3584
80    CONTINUE                                                          MATH3585
C              ---------------------------------------------------------MATH3586
C              REDUCTION OF ALL REMAINING ROWS                          MATH3587
C              ---------------------------------------------------------MATH3588
      DO 91  K=1,N                                                      MATH3589
      IF (K.EQ.I)  GO TO 91                                             MATH3590
      X=A(K,I)                                                          MATH3591
      A(K,I)=0.                                                         MATH3592
      DO 90  J =M,L                                                     MATH3593
      A(K,J)=A(K,J)-X*A(I,J)                                            MATH3594
C              ---------------------------------------------------------MATH3595
C              TEST OVERFLOW SWITCH. IF ON                              MATH3596
C              RETURN NEGATIVE VALUE OF I IN JC(1)                      MATH3597
C              ---------------------------------------------------------MATH3598
      CALL OVERFL (IFL)                                                 MATH3599
      IF(IFL.EQ.1) GO TO 150                                            MATH3100
90    CONTINUE                                                          MATH3101
91    CONTINUE                                                          MATH3102
C              ---------------------------------------------------------MATH3103
C              AX=B AND DET.(A) ARE NOW COMPUTED                        MATH3104
C              ---------------------------------------------------------MATH3105
      GO TO (95,140),KI                                                 MATH3106
C              ---------------------------------------------------------MATH3107
C              PERMUTATION OF THE COLUMNS FOR MATRIX INVERSION          MATH3108
C              ---------------------------------------------------------MATH3109
95    DO 130  J=1,N                                                     MATH3110
      CALL OVERFL(IFL)                                                  MATH3515
      IF (JC(J).EQ.J)  GO TO 130                                        MATH3111
      JJ=J+1                                                            MATH3112
      DO 100  I=JJ,N                                                    MATH3113
      IF (JC(I).EQ.J)  GO TO 110                                        MATH3114
100   CONTINUE                                                          MATH3115
110   JC(I)=JC(J)                                                       MATH3116
      DO 120  K=1,N                                                     MATH3117
      X=A(K,I)                                                          MATH3118
      A(K,I)=A(K,J)                                                     MATH3119
120   A(K,J)=X                                                          MATH3120
130   CONTINUE                                                          MATH3121
140    JC(1)=N                                                          MATH3122
      IF(KD.EQ.1) V(1)=S                                                MATH3123
      RETURN                                                            MATH3124
150   JC(1)=1-I                                                         MATH3125
      IF(KD.EQ.1) V(1)=S                                                MATH3126
      RETURN 1
      END                                                               MATH3128
