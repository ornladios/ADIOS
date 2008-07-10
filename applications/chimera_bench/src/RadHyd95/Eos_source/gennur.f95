SUBROUTINE gennur
!-----------------------------------------------------------------------
!
!    File:         gennur
!    Module:       gennur
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991!
!
!    Date:         5/28/93
!
!    Purpose:
!      To set up the inelastic neutrino-nucleus absorption and
!       scattering function arrays from the data stored in files
!       wicknuem, wicknuep, wicknunu, wicknbnb
!
!    Variables that must be passed through common:
!  nprint      : unit number for print
!  nnugp(n)    : number of energy zones for neutrinos of type n
!  unu(k)      : energy of energy zone k (MeV)
!  dunu(k)     : energy width of energy zone k (MeV)
!  rncem(k,kp) : e-neutrino - nucleus absorption rate (per iron-like
!          nucleus) and per incident neutrino of energy k and final
!          electrons of energy kp as given by Haxton. Absorption and
!          emission rates are computed from this array in subroutine
!          abemhnc, and the results contribute to the arrays ab and
!          em above.
!
!  rncep(k,kp) : e-antineutrino - nucleus absorption rate (per
!          iron-like nucleus) and per incident neutrino of energy k
!          and final electrons of energy kp as given by Haxton.
!          Absorption and emission rates are computed from this array
!          in subroutine abemhnc, and the results contribute to the
!          arrays ab and em above.
!
!  rncnu0(k,kp) - zero moment of the n-neutrino-nucleus non-
!          isoenergetic scattering functions (per iron-like nucleus)
!          and per incident neutrino of energy k and final neutrino
!          energy kp as given by Haxton. n-neutrino-nucleus non-
!          isoenergetic scattering functions are computed from this
!          array in subroutine scatacal, and the results are stored
!          in the array scta0 below.
!
!  rncnb0(k,kp) - zero moment of the n-antineutrino-nucleus non-
!          isoenergetic scattering functions (per iron-like nucleus)
!          and per incident neutrino of energy k and final neutrino
!          energy kp as given by Haxton. n-antineutrino-nucleus non-
!          isoenergetic scattering functions are computed from this
!          array in subroutine scatacal, and the results are stored
!          in the array scta0 below.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!       none
!
!    Include files:
!        numerical_module
!        abem_module, edit_module, nu_dist_module, scat_a_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero

USE abem_module, ONLY : rncem, rncep
USE edit_module, ONLY : nprint
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, unubi
USE scat_a_module, ONLY : rncnu0, rncnb0

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                            :: i             ! do index
INTEGER                            :: j             ! do index
INTEGER                            :: ie            ! neutrino energy do index
INTEGER                            :: je            ! neutrino energy do index
INTEGER                            :: imaxm         ! array indices
INTEGER                            :: iminp         ! array indices
INTEGER                            :: nnugpmx       ! maximum number of neutrino energy groups
INTEGER, PARAMETER                 :: nrdata = 31   ! unit number for accessing the data files
INTEGER                            :: istat         ! open-close file flag

REAL(KIND=double)                  :: ri            ! ratio of energy overlap
REAL(KIND=double)                  :: rj            ! ratio of energy overlap

REAL(KIND=double)                  :: deinu = 1.d0  ! half-width of neutrino energy bin in file WICKNUNU
REAL(KIND=double)                  :: defnu = 1.d0  ! neutrino energy width
REAL(KIND=double)                  :: deinb = 1.d0  ! neutrino energy width
REAL(KIND=double)                  :: defnb = 1.d0  ! neutrino energy width
REAL(KIND=double)                  :: deiem = 0.5d0 ! half-width of neutrino energy bin in file WICKNCEM
REAL(KIND=double)                  :: defem = 0.5d0 ! neutrino energy width
REAL(KIND=double)                  :: deiep = 0.5d0 ! neutrino energy width
REAL(KIND=double)                  :: defep = 0.5d0 ! neutrino energy width
REAL(KIND=double)                  :: dwnunu = 2.d0 ! neutrino energy width
REAL(KIND=double)                  :: dwnbnb = 2.d0 ! neutrino energy width
REAL(KIND=double)                  :: dwnuem = 1.d0 ! neutrino energy width
REAL(KIND=double)                  :: dwnuep = 1.d0 ! neutrino energy width

REAL(KIND=double), DIMENSION(1322) :: eiem          ! centered neutrino initial energy in file WICKNCEM
REAL(KIND=double), DIMENSION(1322) :: efem          ! centered neutrino final energy in file WICKNCEM
REAL(KIND=double), DIMENSION(1322) :: fenuem0       ! zero moment of the neutrino absorption function
REAL(KIND=double), DIMENSION(1322) :: fenuem1       ! first moment of the neutrino absorption function
REAL(KIND=double), DIMENSION(1322) :: fenuem2       ! second moment of the neutrino absorption function

REAL(KIND=double), DIMENSION(580)  :: eiep          ! centered antineutrino initial energy in file WICKNCEP
REAL(KIND=double), DIMENSION(580)  :: efep          ! centered antineutrino final energy in file WICKNCEP
REAL(KIND=double), DIMENSION(580)  :: fenuep0       ! zero moment of the antineutrino absorption function
REAL(KIND=double), DIMENSION(580)  :: fenuep1       ! first moment of the antineutrino absorption function
REAL(KIND=double), DIMENSION(580)  :: fenuep2       ! second moment of the antineutrino absorption function

REAL(KIND=double), DIMENSION(682)  :: einu          ! centered neutrino initial energy in file WICKNUNU
REAL(KIND=double), DIMENSION(682)  :: efnu          ! centered neutrino final energy in file WICKNUNU
REAL(KIND=double), DIMENSION(682)  :: fenunu0       ! zero moment of the neutrino scattering function
REAL(KIND=double), DIMENSION(682)  :: fenunu1       ! first moment of the neutrino scattering function
REAL(KIND=double), DIMENSION(682)  :: fenunu2       ! second moment of the neutrino scattering function

REAL(KIND=double), DIMENSION(682)  :: einb          ! centered antineutrino initial energy in file WICKNUNU
REAL(KIND=double), DIMENSION(682)  :: efnb          ! centered antineutrino final energy in file WICKNUNU
REAL(KIND=double), DIMENSION(682)  :: fenbnb0       ! zero moment of the antineutrino scattering function
REAL(KIND=double), DIMENSION(682)  :: fenbnb1       ! first moment of the antineutrino scattering function
REAL(KIND=double), DIMENSION(682)  :: fenbnb2       ! second moment of the antineutrino scattering function

 1001 format (3X,F5.1,5X,F5.1,2X,3E12.5)
 1003 format (4x,F5.1,5X,F5.1,2X,3(1PE12.5))
 5001 format (' The quantity iminm for the quantity rncnu0 in subroutine gennur cannot be found')
 5002 format (' The quantity imaxm for the quantity rncnu0 in subroutine gennur cannot be found')
 6001 format (' The quantity iminm for the quantity rncnb0 in subroutine gennur cannot be found')
 6002 format (' The quantity imaxm for the quantity rncnb0 in subroutine gennur cannot be found')
 7001 format (' The quantity iminm for the quantity rncem in subroutine gennur cannot be found')
 7002 format (' The quantity imaxm for the quantity rncem in subroutine gennur cannot be found')
 8001 format (' The quantity iminm for the quantity rncep in subroutine gennur cannot be found')
 8002 format (' The quantity imaxm for the quantity rncep in subroutine gennur cannot be found')
 9001 format (' Error in opening file WICKNUEM')
 9101 format (' Error in closing file WICKNUEM')
 9201 format (' Error in opening file WICKNUEP')
 9301 format (' Error in closing file WICKNUEP')
 9401 format (' Error in opening file WICKNUNU')
 9501 format (' Error in closing file WICKNUNU')
 9601 format (' Error in opening file WICKNBNB')
 9701 format (' Error in closing file WICKNBNB')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Initialize.....................................................

nnugpmx            = MAX( nnugp(1) , nnugp(2) , nnugp(3) )

DO i = 1,nnugpmx
  DO j = 1,nnugpmx
    rncnu0(i,j)    = zero
  END DO
END DO

!........Read in data from file wicknuem.d..............................

!-----------------------------------------------------------------------
!        Read in data from files
!         wicknuem.d
!         wicknuep.d
!         wicknunu.d
!         wicknbnb.d
!-----------------------------------------------------------------------

OPEN (UNIT=nrdata,FILE='../../MGFLD/Scat_as/wicknuem.d',STATUS='old',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9001)
  STOP
END IF

DO i = 1,1322
  READ (nrdata,1001) eiem(i),efem(i),fenuem0(i),fenuem1(i),fenuem2(i)
END DO

CLOSE (unit=nrdata,status='KEEP',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9101)
  STOP
END IF

!........Read in data from file wicknuep.d..............................

OPEN (UNIT=nrdata,FILE='../../MGFLD/Scat_as/wicknuep.d',STATUS='old',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9201)
  STOP
END IF

DO i = 1,580
  READ (nrdata,1001) eiep(i),efep(i),fenuep0(i),fenuep1(i),fenuep2(i)
END DO

CLOSE (unit=nrdata,status='KEEP',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9301)
  STOP
END IF

!........Read in data from file wicknuep.d..............................

OPEN (UNIT=nrdata,FILE='../../MGFLD/Scat_as/wicknunu.d',STATUS='old',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9401)
  STOP
END IF

DO i = 1,682
  READ (nrdata,1001) einu(i),efnu(i),fenunu0(i),fenunu1(i),fenunu2(i)
END DO

CLOSE (unit=nrdata,status='KEEP',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9501)
  STOP
END IF

!........Read in data from file wicknuep.d..............................

OPEN (UNIT=nrdata,FILE='../../MGFLD/Scat_as/wicknbnb.d',STATUS='old',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9601)
  STOP
END IF

DO i = 1,682
  READ (nrdata,1001) einb(i),efnb(i),fenbnb0(i),fenbnb1(i),fenbnb2(i)
END DO

CLOSE (unit=nrdata,status='KEEP',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9701)
  STOP
END IF

!........Compute the array rncnu0(nez,nez)..............................
!.......................................................................
!.......................................................................

DO ie = 1,nnugpmx
   DO je = 1,nnugpmx

!-----------------------------------------------------------------------
!        Compute iminp and imaxm, where
!
!         einu(iminp-1) + deinu < unub(ie)   < einu(iminp)   + deinu
!         einu(imaxm)   - deinu < unub(ie+1) < einu(imaxm+1) - deinu
!
!         einu    : centered neutrino initial energy in file WICKNUNU
!         deinu   : half-width of neutrino energy bin in file WICKNUNU
!-----------------------------------------------------------------------

    IF ( ( einu(1  ) - deinu ) < unubi(ie+1) ) THEN
      IF ( ( einu(682) + deinu ) > unubi(ie  ) ) THEN

        DO i = 1,682
          IF ( ( einu(i) + deinu ) > unubi(ie) ) THEN
            iminp  = i
            GO TO 51
          END IF
        END DO
        WRITE (nprint,5001)
        STOP

   51   CONTINUE

        DO i = 682,1,-1
          IF ( ( einu(i) - deinu ) < unubi(ie+1) ) THEN
            imaxm  = i
            GO TO 53
          END IF
        END DO
        WRITE (nprint,5002)
        STOP

   53   CONTINUE

!-----------------------------------------------------------------------
!        Compute ri, the ratio of the overlap of the incident neu-
!         trino bin width and the incident neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------

        ri         = zero
        DO i = iminp,imaxm

          IF ( ( einu(i) - deinu ) < unubi(ie)  .and.  ( einu(i) + deinu ) <= unubi(ie+1) ) THEN 
            ri     = ( einu(i) + deinu - unubi(ie) )/( 2.d0 * deinu )
            GO TO 54
          END IF

          IF ( ( einu(i) - deinu ) < unubi(ie)  .and.  ( einu(i) + deinu ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - unubi(ie) )/( 2.d0 * deinu )
            GO TO 54
          END IF

          IF ( ( einu(i) - deinu ) >= unubi(ie)  .and.  ( einu(i) + deinu ) <= unubi(ie+1) ) THEN 
            ri     = 1.d0
            GO TO 54
          END IF

          IF ( ( einu(i) - deinu ) >= unubi(ie)  .and.  ( einu(i) + deinu ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - ( einu(i) - deinu ) )/( 2.d0 * deinu )
            GO TO 54
          END IF

   54     CONTINUE

!-----------------------------------------------------------------------
!        Compute rj, the ratio of the overlap of the final neu-
!         trino bin width and the final neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------

          rj       = zero

          IF ( ( efnu(i) + defnu ) <= unubi(je) ) GO TO 55

          IF ( ( efnu(i) - defnu ) >= unubi(je+1) ) GO TO 55

          IF ( ( efnu(i) - defnu ) < unubi(je)  .and.  ( efnu(i) + defnu ) <= unubi(je+1) ) THEN 
            rj     = ( efnu(i) + defnu - unubi(je) )/( 2.d0 * defnu )
            GO TO 55
          END IF

          IF ( ( efnu(i) - defnu ) < unubi(je)  .and.  ( efnu(i) + defnu ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - unubi(je) )/( 2.d0 * defnu )
            GO TO 55
          END IF

          IF ( ( efnu(i) - defnu ) >= unubi(je)  .and.  ( efnu(i) + defnu ) <= unubi(je+1) ) THEN 
            rj     = 1.d0
            GO TO 55
          END IF

          IF ( ( efnu(i) - defnu ) >= unubi(je  )  .and.  ( efnu(i) + defnu ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - ( efnu(i) - defnu ) )/( 2.d0 * defnu )
            GO TO 55
          END IF

   55     CONTINUE

!-----------------------------------------------------------------------
!        Sum over final energy bins; divide by unu(je)**2*dunu(je) to
!         convert cross section into a quantity proportional to a
!         scattering function.
!
!        Average over incident energy, i.e., multipy fenunui by the
!         incident energy data spacing, dwnunu, and divide
!         rncnu0(ie,je) by dunu(ie) to obtain the average.
!-----------------------------------------------------------------------
          rncnu0(ie,je) = rncnu0(ie,je) + ri * rj * fenunu0(i) * dwnunu

        END DO

      END IF ! (einu(682) + deinu ) > unubi(ie  )
    END IF ! (einu(1  ) - deinu ) < unubi(ie+1)
  END DO ! end loop over je
END DO ! end loop over ie

DO ie = 1,nnugpmx
  DO je = 1,nnugpmx
    IF ( je > ie ) rncnu0(ie,je) = zero
    rncnu0(ie,je)= rncnu0(ie,je)/( ( dunui(ie) ) * ( dunui(je) * unui(je) * unui(je) ) )
  END DO
END DO

!........Compute the array rncnb0(nez,nez)..............................
!.......................................................................
!.......................................................................

DO ie = 1,nnugpmx
  DO je = 1,nnugpmx

!-----------------------------------------------------------------------
!        Compute iminp and imaxm, where
!
!         einb(iminp-1) + deinb < unub(ie)   < einb(iminp)   + deinb
!         einb(imaxm)   - deinb < unub(ie+1) < einb(imaxm+1) - deinb
!
!         einb    : centered neutrino initial energy in file WICKNBNB
!         deinb   : half-width of neutrino energy bin in file WICKNBNB
!-----------------------------------------------------------------------

    IF ( ( einb(1) - deinb) < unubi(ie+1) ) THEN
      IF ( ( einb(682) + deinb) > unubi(ie) ) THEN

        DO 60 i = 1,682
          IF ( (einb(i) + deinb) > unubi(ie) ) THEN
            iminp  = i
            GO TO 61
          END IF
   60   CONTINUE
        WRITE (nprint,6001)
        STOP

   61   CONTINUE

        DO 62 i = 682,1,-1
          IF ( ( einb(i) - deinb ) < unubi(ie+1) ) THEN
            imaxm  = i
            GO TO 63
          END IF
   62   CONTINUE
        WRITE (nprint,6002)
        STOP

   63   CONTINUE

!-----------------------------------------------------------------------
!        Compute ri, the ratio of the overlap of the incident neu-
!         trino bin width and the incident neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------

        ri         = zero
        DO i = iminp,imaxm

          IF ( ( einb(i) - deinb ) < unubi(ie)  .and.  ( einb(i) + deinb ) <= unubi(ie+1) ) THEN 
            ri     = ( einb(i) + deinb - unubi(ie) )/( 2.d0 * deinb )
            GO TO 64
          END IF

          IF ( ( einb(i) - deinb ) < unubi(ie)  .and.  ( einb(i) + deinb ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - unubi(ie) )/( 2.d0 * deinb )
            GO TO 64
          END IF

          IF ( ( einb(i) - deinb ) >= unubi(ie)  .and.  ( einb(i) + deinb ) <= unubi(ie+1) ) THEN 
            ri     = 1.d0
            GO TO 64
          END IF

          IF ( ( einb(i) - deinb ) >= unubi(ie)  .and.  ( einb(i) + deinb ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - ( einb(i) - deinb ) )/( 2.d0 * deinb )
            GO TO 64
          END IF

   64     CONTINUE

!-----------------------------------------------------------------------
!        Compute rj, the ratio of the overlap of the final neu-
!         trino bin width and the final neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------

          rj       = zero

          IF ( ( efnb(i) + defnb ) <= unubi(je) ) GO TO 65

          IF ( ( efnb(i) - defnb ) >= unubi(je+1) ) GO TO 65

          IF ( ( efnb(i) - defnb ) < unubi(je)  .and.  ( efnb(i) + defnb ) <= unubi(je+1) ) THEN 
            rj     = ( efnb(i) + defnb - unubi(je) )/( 2.d0 * defnb )
            GO TO 65
          END IF

          IF ( ( efnb(i) - defnb ) < unubi(je)  .and.  ( efnb(i) + defnb ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - unubi(je) )/( 2.d0 * defnb )
            GO TO 65
          END IF

          IF ( ( efnb(i) - defnb ) >= unubi(je)  .and.  ( efnb(i) + defnb ) <= unubi(je+1) ) THEN 
            rj     = 1.d0
            GO TO 65
          END IF

          IF ( ( efnb(i) - defnb ) >= unubi(je)  .and.  ( efnb(i) + defnb ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - ( efnb(i) - defnb ) )/( 2.d0 * defnb )
            GO TO 65
          END IF

   65     CONTINUE
!-----------------------------------------------------------------------
!        Sum over final energy bins; divide by unu(je)**2*dunu(je) to
!         convert cross section into a quantity proportional to a
!         scattering function.
!
!        Average over incident energy, i.e., multipy fenbnbi by the
!         incident energy data spacing, dwnbnb, and divide
!         rncnb0(ie,je) by dunu(ie) to obtain the average.
!-----------------------------------------------------------------------

          rncnb0(ie,je) = rncnb0(ie,je) + ri * rj * fenbnb0(i) * dwnbnb

        END DO

      END IF ! (einb(682) + deinb) > unubi(ie)
    END IF ! (einb(1) - deinb)   < unubi(ie+1)
  END DO ! END loop over je
END DO ! END loop over ie

DO  ie = 1,nnugpmx
  DO je = 1,nnugpmx
    IF ( je > ie ) rncnb0(ie,je) = zero
    rncnb0(ie,je) = rncnb0(ie,je)/( ( dunui(ie) ) * ( dunui(je) * unui(je) * unui(je) ) )
  END DO
END DO

!........Compute the array rncem(nez,nez)...............................
!.......................................................................
!.......................................................................

DO ie = 1,nnugpmx
   DO je = 1,nnugpmx

!-----------------------------------------------------------------------
!        Compute iminp and imaxm, where
!
!         eiem(iminp-1) + deiem < unub(ie)   < eiem(iminp)   + deiem
!         eiem(imaxm)   - deiem < unub(ie+1) < eiem(imaxm+1) - deiem
!
!         eiem    : centered neutrino initial energy in file WICKNCEM
!         deiem   : half-width of neutrino energy bin in file WICKNCEM
!-----------------------------------------------------------------------

    IF ( ( eiem(1) - deiem ) < unubi(ie+1) ) THEN
      IF ( ( eiem(1322) + deiem) > unubi(ie) ) THEN

        DO i = 1,1322
          IF ( ( eiem(i) + deiem ) > unubi(ie) ) THEN
            iminp  = i
            GO TO 71
          END IF
        END DO
        WRITE (nprint,7001)
        STOP

   71   CONTINUE

        DO i = 1322,1,-1
          IF ( (eiem(i) - deiem) < unubi(ie+1) ) THEN
            imaxm  = i
            GO TO 73
          END IF
        END DO
        WRITE (nprint,7002)
        STOP

   73   CONTINUE

!-----------------------------------------------------------------------
!        Compute ri, the ratio of the overlap of the incident neu-
!         trino bin width and the incident neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------

        ri         = zero
        DO i = iminp,imaxm

          IF ( ( eiem(i) - deiem ) < unubi(ie)  .and.  ( eiem(i) + deiem ) <= unubi(ie+1) ) THEN 
            ri     = ( eiem(i) + deiem - unubi(ie) )/( 2.d0 * deiem )
            GO TO 74
          END IF

          IF ( ( eiem(i) - deiem ) < unubi(ie)  .and.  ( eiem(i) + deiem ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - unubi(ie) )/( 2.d0 * deiem )
            GO TO 74
          END IF

          IF ( ( eiem(i) - deiem ) >= unubi(ie)  .and.  ( eiem(i) + deiem ) <= unubi(ie+1) ) THEN 
            ri     = 1.d0
            GO TO 74
          END IF

          IF ( ( eiem(i) - deiem ) >= unubi(ie)  .and.  ( eiem(i) + deiem ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - ( eiem(i) - deiem ) )/( 2.d0 * deiem )
            GO TO 74
          END IF

   74     CONTINUE

!-----------------------------------------------------------------------
!        Compute rj, the ratio of the overlap of the final neu-
!         trino bin width and the final neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------

          rj       = zero

          IF ( ( efem(i) + defem ) <= unubi(je) ) GO TO 75

          IF ( ( efem(i) - defem ) >= unubi(je+1) ) GO TO 75

          IF ( ( efem(i) - defem ) < unubi(je)  .and.  ( efem(i) + defem ) <= unubi(je+1) ) THEN 
            rj     = ( efem(i) + defem - unubi(je) )/( 2.d0 * defem )
            GO TO 75
          END IF

          IF ( ( efem(i) - defem ) < unubi(je)  .and.  ( efem(i) + defem ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - unubi(je) )/( 2.d0 * defem )
            GO TO 75
          END IF

          IF ( ( efem(i) - defem ) >= unubi(je)  .and.  ( efem(i) + defem ) <= unubi(je+1) ) THEN 
            rj     = 1.d0
            GO TO 75
          END IF

          IF ( ( efem(i) - defem ) >= unubi(je)  .and.  ( efem(i) + defem ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - ( efem(i) - defem ) )/( 2.d0 * defem )
            GO TO 75
          END IF

   75     CONTINUE

!-----------------------------------------------------------------------
!        Sum over final energy bins
!
!        Average over incident energy, i.e., multipy fenuemi by the
!         incident energy data spacing, dwnuem, and divide
!         rncem(ie,je) by dunu(ie) to obtain the average.
!-----------------------------------------------------------------------

          rncem(ie,je) = rncem(ie,je) + ri*rj*fenuem0(i)*dwnuem

        END DO

      END IF ! (eiem(1322) + deiem) > unubi(ie  )
    END IF ! (eiem(1) - deiem   ) < unubi(ie+1)
  END DO ! END loop over je
END DO ! END loop over ie

DO ie = 1,nnugpmx
  DO je = 1,nnugpmx
    IF ( je > ie ) rncem(ie,je) = zero
    rncem(ie,je) = rncem(ie,je)/dunui(ie)
  END DO
END DO

!........Compute the array rncep(nez,nez)...............................
!.......................................................................
!.......................................................................

DO ie = 1,nnugpmx
  DO je = 1,nnugpmx

!-----------------------------------------------------------------------
!        Compute iminp and imaxm, where
!
!         eiep(iminp-1) + deiep < unub(ie)   < eiep(iminp)   + deiep
!         eiep(imaxm)   - deiep < unub(ie+1) < eiep(imaxm+1) - deiep
!
!         eiep    : centered neutrino initial energy in file WICKNCEP
!         deiep   : half-width of neutrino energy bin in file WICKNCEP
!-----------------------------------------------------------------------

    IF ( ( eiep(1) - deiep ) < unubi(ie+1) ) THEN
      IF ( ( eiep(580) + deiep ) > unubi(ie) ) THEN

        DO i = 1,580
          IF ( ( eiep(i) + deiep ) > unubi(ie) ) THEN
            iminp  = i
            GO TO 81
          END IF
        END DO
        WRITE (nprint,8001)
        STOP

   81   CONTINUE

        DO i = 580,1,-1
          IF ( ( eiep(i) - deiep ) < unubi(ie+1) ) THEN
            imaxm  = i
            GO TO 83
          END IF
        END DO
        WRITE (nprint,8002)
        STOP

   83   CONTINUE

!-----------------------------------------------------------------------
!        Compute ri, the ratio of the overlap of the incident neu-
!         trino bin width and the incident neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------

        ri         = zero
        DO i = iminp,imaxm

          IF ( ( eiep(i) - deiep ) < unubi(ie  )  .and.  ( eiep(i) + deiep ) <= unubi(ie+1) ) THEN 
            ri     = ( eiep(i) + deiep - unubi(ie) )/( 2.d0 * deiep )
            GO TO 84
          END IF

          IF ( ( eiep(i) - deiep ) < unubi(ie  )  .and.  ( eiep(i) + deiep ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - unubi(ie) )/( 2.d0 * deiep )
            GO TO 84
          END IF

          IF ( ( eiep(i) - deiep ) >= unubi(ie  )  .and.  ( eiep(i) + deiep ) <= unubi(ie+1) ) THEN 
            ri     = 1.d0
            GO TO 84
          END IF

          IF ( ( eiep(i) - deiep ) >= unubi(ie)  .and.  ( eiep(i) + deiep ) > unubi(ie+1) ) THEN 
            ri     = ( unubi(ie+1) - ( eiep(i) - deiep ) )/( 2.d0 * deiep )
            GO TO 84
          END IF

   84     CONTINUE

!-----------------------------------------------------------------------
!        Compute rj, the ratio of the overlap of the final neu-
!         trino bin width and the final neutrino energy group
!         width to the incident neutrino bin width.
!-----------------------------------------------------------------------
          rj       = zero

          IF ( ( efep(i) + defep ) <= unubi(je) ) GO TO 85

          IF ( ( efep(i) - defep ) >= unubi(je+1) ) GO TO 85

          IF ( ( efep(i) - defep ) < unubi(je)  .and.  ( efep(i) + defep ) <= unubi(je+1) ) THEN 
            rj     = ( efep(i) + defep - unubi(je) )/( 2.d0 * defep )
            GO TO 85
          END IF

          IF ( ( efep(i) - defep ) < unubi(je)  .and.  ( efep(i) + defep ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - unubi(je) )/( 2.d0 * defep )
            GO TO 85
          END IF

          IF ( ( efep(i) - defep ) >= unubi(je)  .and.  ( efep(i) + defep ) <= unubi(je+1) ) THEN 
            rj     = 1.d0
            GO TO 85
          END IF

          IF ( ( efep(i) - defep ) >= unubi(je)  .and.  ( efep(i) + defep ) > unubi(je+1) ) THEN 
            rj     = ( unubi(je+1) - ( efep(i) - defep ) )/( 2.d0 * defep )
            GO TO 85
          END IF

   85     CONTINUE

!-----------------------------------------------------------------------
!        Sum over final energy bins
!
!        Average over incident energy, i.e., multipy fenuepi by the
!         incident energy data spacing, dwnuep, and divide
!         rncep(ie,je) by dunu(ie) to obtain the average
!-----------------------------------------------------------------------

          rncep(ie,je) = rncep(ie,je) + ri * rj * fenuep0(i) * dwnuep

        END DO

      END IF ! (eiep(580) + deiep) > unubi(ie)
    END IF ! (eiep(1) - deiep) < unubi(ie+1)
  END DO ! END loop over je
END DO ! END loop over ie

DO ie = 1,nnugpmx
  DO je = 1,nnugpmx
    IF ( je > ie ) rncep(ie,je) = zero
    rncep(ie,je) = rncep(ie,je)/dunui(ie)
  END DO
END DO

RETURN
END SUBROUTINE gennur
