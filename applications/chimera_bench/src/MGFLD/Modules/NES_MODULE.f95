!-----------------------------------------------------------------------
!    Module:       nes_module
!    Author:       S. W. Bruenn
!    Date:         6/26/03
!-----------------------------------------------------------------------

MODULE nes_module

USE kind_module
SAVE

!-----------------------------------------------------------------------
!  Neutrino energies
!-----------------------------------------------------------------------

REAL(KIND=double)                :: w             ! e_in/kt
REAL(KIND=double)                :: w2            ! w**2
REAL(KIND=double)                :: w3            ! w**3
REAL(KIND=double)                :: w4            ! w**4
REAL(KIND=double)                :: w5            ! w**5
REAL(KIND=double)                :: w6            ! w**6

REAL(KIND=double)                :: wp            ! e_out/kt
REAL(KIND=double)                :: wp2           ! wp**2
REAL(KIND=double)                :: wp3           ! wp**3
REAL(KIND=double)                :: wp4           ! wp**4
REAL(KIND=double)                :: wp5           ! wp**5
REAL(KIND=double)                :: wp6           ! wp**6

REAL(KIND=double)                :: wwp           ! w*wp
REAL(KIND=double)                :: w2pwp2        ! wwp**2

REAL(KIND=double)                :: wp_w          ! wp - w
REAL(KIND=double)                :: w_wp          ! w - wp
REAL(KIND=double)                :: w_wp2         ! w_wp**2
REAL(KIND=double)                :: w_wp3         ! w_wp**3
REAL(KIND=double)                :: w_wp4         ! w_wp**4
REAL(KIND=double)                :: w_wp5         ! w_wp**5

REAL(KIND=double)                :: w_p_wp        ! w + wp
REAL(KIND=double)                :: w_p_wp2       ! w_p_wp**2
REAL(KIND=double)                :: w_p_wp3       ! w_p_wp**3

REAL(KIND=double)                :: e             ! incident electron energy/kt
REAL(KIND=double)                :: e2            ! e**2
REAL(KIND=double)                :: e3            ! e**3
REAL(KIND=double)                :: e4            ! e**4
REAL(KIND=double)                :: e5            ! e**5
REAL(KIND=double)                :: e6            ! e**6
REAL(KIND=double)                :: e7            ! e**7
REAL(KIND=double)                :: e_min         ! minimum incident electron energy/kt

!-----------------------------------------------------------------------
!  Constants
!-----------------------------------------------------------------------

REAL(KIND=double)                :: r2_3         = 2.d+00/3.d+00
REAL(KIND=double)                :: r2_15        = 2.d+00/15.d+00
REAL(KIND=double)                :: r2_105       = 2.d+00/105.d+00
REAL(KIND=double)                :: r4_3         = 4.d+00/3.d+00
REAL(KIND=double)                :: r4_5         = 4.d+00/5.d+00
REAL(KIND=double)                :: r4_15        = 4.d+00/15.d+00
REAL(KIND=double)                :: r8_3         = 8.d+00/3.d+00
REAL(KIND=double)                :: r8_5         = 8.d+00/5.d+00
REAL(KIND=double)                :: r8_7         = 8.d+00/7.d+00
REAL(KIND=double)                :: r8_15        = 8.d+00/15.d+00
REAL(KIND=double)                :: r12_5        = 12.d+00/5.d+00
REAL(KIND=double)                :: r12_35       = 12.d+00/35.d+00
REAL(KIND=double)                :: r16_3        = 16.d+00/3.d+00
REAL(KIND=double)                :: r16_5        = 16.d+00/5.d+00
REAL(KIND=double)                :: r16_30       = 16.d+00/30.d+00
REAL(KIND=double)                :: r16_35       = 16.d+00/35.d+00
REAL(KIND=double)                :: r28_5        = 28.d+00/5.d+00
REAL(KIND=double)                :: r28_15       = 28.d+00/15.d+00
REAL(KIND=double)                :: r36_105      = 36.d+00/105.d+00


END module nes_module
