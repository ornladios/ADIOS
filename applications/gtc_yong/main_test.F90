!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!                 Gyrokinetic Toroidal Code (GTC)                            !
!                          Version 2                                         !
!          Zhihong Lin, Stephane Ethier, Jerome Lewandowski                  !
!              Princeton Plasma Physics Laboratory                           !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gtc
  use global_parameters
  use particle_array
  use particle_tracking
  use field_array
  use diagnosis_array
  implicit none

  integer i,ierror,m
  real(wp) :: avg_weight,avg_abs_weight
  real(doubleprec) time(8),timewc(8),t0,dt,t0wc,dtwc,loop_time
  real(doubleprec) tracktcpu,tracktwc,tr0,tr0wc
  character(len=10) ic(8)

! MPI initialize
  call mpi_init(ierror)

! input parameters, setup equilibrium, allocate memory 
  CALL SETUP

! initialize particle position and velocity
  CALL LOAD

  avg_weight=0.0_wp
  do m=1,mi
     avg_weight=avg_weight+zion(5,m)
     avg_abs_weight=avg_abs_weight+abs(zion(5,m))
  enddo
  avg_weight=avg_weight/real(mi)
  avg_abs_weight=avg_abs_weight/real(mi)

  write(0,*)mype,'  avg_weight =',avg_weight,'   avg_abs_weight =',avg_abs_weight

! MPI finalize
  call mpi_finalize(ierror)

end program gtc

!=========================================
subroutine timer(t0,dt,t0wc,dtwc)
!=========================================
  use precision
  implicit none
  real(doubleprec) t0,dt,t0wc,dtwc
  real(doubleprec) t1,t1wc

! Get cpu usage time since the beginning of the run and subtract value
! from the previous call
  call cpu_time(t1)
  dt=t1-t0
  t0=t1

! Get wall clock time and subtract value from the previous call
  t1wc=MPI_WTIME()
  dtwc=t1wc-t0wc
  t0wc=t1wc

end subroutine timer

