!-----------------------------------------------------------------------
!    Module:       bomb_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE bomb_module

USE kind_module

SAVE


!-----------------------------------------------------------------------
!  Explosion Simulation variables
!-----------------------------------------------------------------------
!  e_bomb: energy added to generate explosion (ergs).
!
!  bomb_time: time during which energy is added.
!
!  t_start_bomb: time at which energy addition begins.
!
!  jexpl_min: inner zone for energy addition.
!
!  jexpl_max: outer zone for energy addition.
!-----------------------------------------------------------------------

REAL(KIND=double)                              :: e_bomb
REAL(KIND=double)                              :: bomb_time
REAL(KIND=double)                              :: t_start_bomb
INTEGER                                        :: jexpl_min
INTEGER                                        :: jexpl_max

END module bomb_module
