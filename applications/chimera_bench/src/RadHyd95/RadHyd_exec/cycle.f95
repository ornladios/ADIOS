SUBROUTINE cycle( time, dtnph )
!-----------------------------------------------------------------------
!
!    File:         cycle
!    Module:       cycle
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To update the cycle number
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  time  : elapsed time
!  dtnph : current time step
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  cycle_module, edit_module, parallel_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, data_path, nlog
USE parallel_module, ONLY : myid
USE radial_ray_module, ONLY : ncycle_r=>ncycle
USE t_cntrl_module, ONLY : jrdt, jadt, jzdt, dt

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: time            ! elapsed time
REAL(KIND=double), INTENT(in)    :: dtnph           ! current time step

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: iunit2 = 22     ! unit number to open cycle file
INTEGER                          :: istat           ! open file flag
INTEGER, DIMENSION(50)           :: itcntrl         ! time step criteria index
INTEGER                          :: i               ! do index
INTEGER                          :: i_min           ! index of criterion giving minimum time step

REAL(KIND=double)                :: dtmin           ! minimum of i time step criteria

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (9x,'dt (x-hydro Courant criterion)                   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  103 FORMAT (9x,'dt (x-hydro density change criterion)            =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  105 FORMAT (9x,'dt (x-hydro temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  107 FORMAT (9x,'dt (y-hydro Courant criterion)                   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  109 FORMAT (9x,'dt (y-hydro density change criterion)            =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  111 FORMAT (9x,'dt (y-hydro temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  113 FORMAT (9x,'dt (z-hydro Courant criterion)                   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  115 FORMAT (9x,'dt (z-hydro density change criterion)            =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  117 FORMAT (9x,'dt (z-hydro temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  119 FORMAT (9x,'dt (convection crossing criterion)               =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  121 FORMAT (9x,'dt (nuclear temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  123 FORMAT (9x,'dt (nuclear xn change criterion)                 =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  125 FORMAT (9x,'dt (max increase criterion)                      =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  127 FORMAT (9x,'dt (all nu temperature change criterion)         =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  129 FORMAT (9x,'dt (all nu electron fraction change criterion)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  131 FORMAT (9x,'dt (e_nu    psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  133 FORMAT (9x,'dt (e_nubar psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  135 FORMAT (9x,'dt (x_nu    psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  137 FORMAT (9x,'dt (x_nubar psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  139 FORMAT (9x,'dt (e_nu    psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  141 FORMAT (9x,'dt (e_nubar psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  143 FORMAT (9x,'dt (x_nu    psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  145 FORMAT (9x,'dt (x_nubar psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
 1001 FORMAT (' Current cycle number is',i7)
 1003 FORMAT (' Current cycle number is',i7,' elapsed time is', es15.8, &
 & ' current time step is',es15.8)
 8501 FORMAT (' Error in closing iunit2 ("cycle") in subroutine cycle')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Index the time step control criteria
!-----------------------------------------------------------------------

IF ( first ) THEN
  first           = .false.
  itcntrl( 1)     = 1
  itcntrl( 2)     = 2
  itcntrl( 3)     = 3
  itcntrl( 4)     = 4
  itcntrl( 5)     = 5
  itcntrl( 6)     = 6
  itcntrl( 7)     = 7
  itcntrl( 8)     = 8
  itcntrl( 9)     = 9
  itcntrl(10)     = 31
  itcntrl(11)     = 33
  itcntrl(12)     = 34
  itcntrl(13)     = 10
  itcntrl(14)     = 11
  itcntrl(15)     = 16
  itcntrl(16)     = 21
  itcntrl(17)     = 22
  itcntrl(18)     = 23
  itcntrl(19)     = 24
  itcntrl(20)     = 41
  itcntrl(21)     = 42
  itcntrl(22)     = 43
  itcntrl(23)     = 44
END IF ! first

!-----------------------------------------------------------------------
!  Determine controling criterion
!-----------------------------------------------------------------------

dtmin             = 1.d+20
i_min             = 10
DO i = 1, 23
  IF ( dt(itcntrl(i)) < dtmin ) THEN
    dtmin         = dt(itcntrl(i))
    i_min         = itcntrl(i)
  END IF ! dt(itcntrl(i)) < dtmin
END DO ! i = 1, 19

!-----------------------------------------------------------------------
!  Update the cycle number
!-----------------------------------------------------------------------

ncycle            = ncycle + 1
ncycle_r          = ncycle

IF ( myid == 0 ) THEN
  OPEN (UNIT=iunit2, FILE=TRIM(data_path)//'/Run_Log/cycle.d', STATUS='new', IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=iunit2, FILE=TRIM(data_path)//'/Run_Log/cycle.d', STATUS='old')
  WRITE (iunit2,1001) ncycle
  WRITE (nprint,1001) ncycle
  CLOSE (UNIT=iunit2,STATUS='keep',ERR=8500)
END IF ! myid == 0
WRITE (nlog,1003) ncycle, time, dtnph

tcrit: SELECT CASE (i_min)

  CASE(1)
    WRITE (nlog,101) dt( 1),jrdt( 1),jadt( 1)
  CASE(2)
    WRITE (nlog,103) dt( 2),jrdt( 2),jadt( 2)
  CASE(3)
    WRITE (nlog,105) dt( 3),jrdt( 3),jadt( 3)
  CASE(4)
    WRITE (nlog,107) dt( 4),jrdt( 4),jadt( 4)
  CASE(5)
    WRITE (nlog,109) dt( 5),jrdt( 5),jadt( 5)
  CASE(6)
    WRITE (nlog,111) dt( 6),jrdt( 6),jadt( 6)
  CASE(7)
    WRITE (nlog,113) dt( 7),jrdt( 7),jadt( 7)
  CASE(8)
    WRITE (nlog,115) dt( 8),jrdt( 8),jadt( 8)
  CASE(9)
    WRITE (nlog,117) dt( 9),jrdt( 9),jadt( 9)
  CASE(31)
    WRITE (nlog,119) dt(31),jrdt(31),jadt(31)
  CASE(33)
    WRITE (nlog,121) dt(33),jrdt(33),jadt(33)
  CASE(34)
    WRITE (nlog,123) dt(34),jrdt(34),jadt(34)
  CASE(10)
    WRITE (nlog,125) dt(10),jrdt(10),jadt(10)
  CASE(11)
    WRITE (nlog,127) dt(11),jrdt(11),jadt(11)
  CASE(16)
    WRITE (nlog,129) dt(16),jrdt(16),jadt(16)
  CASE(21)
    WRITE (nlog,131) dt(21),jrdt(21),jadt(21)
  CASE(22)
    WRITE (nlog,133) dt(22),jrdt(22),jadt(22)
  CASE(23)
    WRITE (nlog,135) dt(23),jrdt(23),jadt(23)
  CASE(24)
    WRITE (nlog,137) dt(24),jrdt(25),jadt(24)
  CASE(41)
    WRITE (nlog,139) dt(41),jrdt(41),jadt(41)
  CASE(42)
    WRITE (nlog,141) dt(42),jrdt(42),jadt(42)
  CASE(43)
    WRITE (nlog,143) dt(43),jrdt(43),jadt(43)
  CASE(44)
    WRITE (nlog,145) dt(44),jrdt(44),jadt(44)

END SELECT tcrit

RETURN

 8500 write (nprint,8501)
STOP

END SUBROUTINE cycle
