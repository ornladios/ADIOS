!-----------------------------------------------------------------------
!    Module:       boundary_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE boundary_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Outer boundary
!-----------------------------------------------------------------------
!  ipbnd   : outer pressure boundary codition switch.
!
!     ipbnd = 0 : outer pressure boundary condition not imposed.
!     ipbnd = 1 : outer pressure boundary condition imposed.
!
!  pbound  : the value of the outer pressure boundary condition.
!
!  iubcjmx : outer velocity boundary codition switch.
!
!     iubcjmx = 0 : outer velocity (j=jmax) not imposed.
!     iubcjmx = 1 : outer velocity (j=jmax) imposed.
!
!  ubcjmx : the value of the outer velocity.
!-----------------------------------------------------------------------

INTEGER                                        :: ipbnd
REAL(KIND=double)                              :: pbound

INTEGER                                        :: iubcjmx
REAL(KIND=double)                              :: ubcjmx

!-----------------------------------------------------------------------
!  Inner boundary
!-----------------------------------------------------------------------
!  iubcjmn : inner velocity boundary codition switch.
!
!     iubcjmn = 0 : inner velocity (j=1) not imposed.
!     iubcjmn = 1 : inner velocity (j=1) imposed.
!
!  ubcjmn  : the value of the inner velocity.
!-----------------------------------------------------------------------

INTEGER                                        :: iubcjmn
REAL(KIND=double)                              :: ubcjmn

!-----------------------------------------------------------------------
!  Intermediate zones
!-----------------------------------------------------------------------
!  iubcjmn: intermediate zone velocity boundary codition switch.
!
!  iuset = 0 : no velocities imposed for 1 < j < jmax.
!  iuset = 1 : velocities imposed for 1 < j < jmax.
!
!  uset(j) : the value of the velocity imposed for zone j.
!-----------------------------------------------------------------------

INTEGER                                        :: iuset
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: uset

!-----------------------------------------------------------------------
!  zero velocity cutoff
!-----------------------------------------------------------------------
!  r_u : zero velocity cutoff
!
!     r_u >= 0.: velocity computed normally
!     r_u <  0.: velocity set to zero for r(j) > abs(r_u)
!-----------------------------------------------------------------------

REAL(KIND=double)                              :: r_u


END module boundary_module
