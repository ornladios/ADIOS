!-----------------------------------------------------------------------
!    Module:       rezone_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE rezone_module

USE kind_module
SAVE


!-----------------------------------------------------------------------
!
!         Lagrangian rezoning controls
!
!       n_lgrgrid : 1 - generate a smooth Lagrangian grid
!                   2 - generate three separate grids from m_1 separated m_2 and m_3
!
!       m_1 : mass of first zone
!
!       m_2 : mass separating grid 1 from grid 2
!
!       m_3 : mass separating grid 2 from grid 3
!
!       n1zoom : number of zones between m_1 and m_2
!
!       n2zoom : number of zones between m_2 and m_3
!
!       n3zoom : number of zones between m_3 and the outer edge
!
!       zoome1 : mass ratio for rezoning the first n1 zones
!
!       zoome2 : mass ratio for rezoning the next n2zoom zones
!
!       zoome3 : mass ratio for rezoning the remaining n3zoom zones
!
!
!         Eulerian rezoning controls
!
!       n_eulgrid : 1 - generate a smooth Eulerian grid
!                   2 - generate three separate grids from r_1 separated r_2 and r_3
!
!       r_1 : radius of first zone
!
!       r_2 : radius separating grid 1 from grid 2
!
!       r_3 : radius separating grid 2 from grid 3
!
!       n1zoom : number of zones between r_1 and r_2
!
!       n2zoom : number of zones between r_2 and r_3
!
!       n3zoom : number of zones between r_3 and the outer edge
!
!       zoome1 : radius ratio for rezoning the first n1zoom zones
!
!       zoome2 : radius ratio for rezoning the next n2zoom zones
!
!       zoome3 : radius ratio for rezoning the remaining n3zoom zones
!
!
!-----------------------------------------------------------------------

INTEGER                                        :: n_eulgrid
INTEGER                                        :: n_lgrgrid
INTEGER                                        :: n1zoom
INTEGER                                        :: n2zoom
INTEGER                                        :: n3zoom

REAL(KIND=double)                              :: r_1
REAL(KIND=double)                              :: r_2
REAL(KIND=double)                              :: r_3
REAL(KIND=double)                              :: m_1
REAL(KIND=double)                              :: m_2
REAL(KIND=double)                              :: m_3
REAL(KIND=double)                              :: zoome1
REAL(KIND=double)                              :: zoome2
REAL(KIND=double)                              :: zoome3

END module rezone_module
