SUBROUTINE edit_Global_energetics( imin, imax, nx, time, t_tb, x_in,  &
& grav_pot_a_ray, ke_a_ray, e_no_bind_a_ray, e_bind_a_ray, e_bind_fnl_a_ray, &
& dudt_nuc_min, dudt_nuc_max, dudt_nuc_mean, dudt_nuc_zone, dudt_nu_min, &
& dudt_nu_max, dudt_nu_mean, dudt_nu_zone, c_shock, nprint )
!-----------------------------------------------------------------------
!
!    File:         edit_Global_energetics
!    Module:       edit_Global_energetics
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/2/00
!
!    Purpose:
!      To edit the minimum, maximum, and mean values of quantities as a
!      function of the radius.
!
!    Subprograms called:
!      date_and_time_print
!
!    Input arguments:
!  imin             : minimum x-array index for the edit
!  imax             : maximum x-array index for the edit
!  nx               : x-array extent
!  nnu              : neutrino flavor array extent
!  time             : elapsed time
!  t_tb             : time from core bounce
!  x_in             : radial midpoint of zone (cm)
!  grav_pot_a_ray   : angular ray gravitational potential (ergs)
!  ke_a_ray         : angular ray kinetic energy (ergs)
!  e_no_bind_a_ray  : angular ray internal minus binding energy (ergs)
!  e_bind_a_ray     : angular ray current binding energy (ergs)
!  e_bind_fnl_a_ray : angular ray final binding energy (ergs)
!  dudt_nuc_min     : minimum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
!  dudt_nuc_max     : maximum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
!  dudt_nuc_mean    : mass averaged energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
!  dudt_nuc_zone    : energy generation rate by nuclear reactions (B s^{-1} radial-zone^{-1})
!  dudt_nu_min      : minimum neutrino energy deposition rate (ergs g^{-1} s^{1})
!  dudt_nu_max      : maximum neutrino energy deposition rate (ergs g^{-1} s^{1})
!  dudt_nu_mean     : mass averaged neutrino energy deposition rate (ergs g^{-1} s^{1})
!  dudt_nu_zone     : neutrino energy deposition rate  (B s^{-1} radial-zone^{-1})
!  c_shock          : shock location
!  nprint           : unit numberto print data
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module, numerical_module
!  edit_head
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE edit_module, ONLY : head


IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                              :: imin              ! minimum unshifted radial zone index
INTEGER, INTENT(in)                              :: imax              ! maximum unshifted radial zone index
INTEGER, INTENT(in)                              :: nx                ! x-array extent
INTEGER, INTENT(in)                              :: nprint            ! unit numberto print data

REAL(KIND=double), INTENT(in)                    :: time              ! elapsed time
REAL(KIND=double), INTENT(in)                    :: t_tb              ! time from core bounce
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: x_in              ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: grav_pot_a_ray    ! angular ray gravitational potential (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: ke_a_ray          ! angular ray kinetic energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: e_no_bind_a_ray   ! angular ray internal minus binding energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: e_bind_a_ray      ! angular ray current binding energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: e_bind_fnl_a_ray  ! angular ray final binding energy (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_mean     ! mass averaged energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_min      ! minimum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_max      ! maximum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nuc_zone     ! energy generation rate by nuclear reactions (B s^{-1} radial-zone^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_mean      ! mass averaged neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_min       ! minimum neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_max       ! maximum neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: dudt_nu_zone      ! neutrino energy deposition rate  (B s^{-1} radial-zone^{-1})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(in), DIMENSION(nx)     :: c_shock           ! shock location

INTEGER                                          :: i                 ! do index

REAL(KIND=double), DIMENSION(nx)                 :: de_bind_a_ray     ! current - final binding energy
REAL(KIND=double), DIMENSION(nx)                 :: e_tot_no_bind     ! total energy - binding energy
REAL(KIND=double), DIMENSION(nx)                 :: e_tot_no_bind_sum ! e_tot_no_bind summed from outer edge
REAL(KIND=double), DIMENSION(nx)                 :: de_bind_sum       ! current - final binding energy summed from outer edge
REAL(KIND=double), DIMENSION(nx)                 :: e_tot_sum         ! total energy summed from outer edge

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1a)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (9x,'Elapsed time=',es14.7,10x,' Time from bounce=',es14.7/)
    9 FORMAT (41x,'Energetics')
   11 FORMAT (40x,12('-')/)
  101 FORMAT ('   i      r        e_grav/z    e_kin/z   e_int-b/z  e_tot-b/z  e_bind/z   de_bind/z &
&  e_tot-b    de_bind     e_tot    dudt_nuc   dudt_ncmn  dudt_ncmx   dudt_nu   dudt_numn &
& dudt_numx  dudt_nczn  dudtnuzn'/)
  103 FORMAT (1x,i4,es11.3,a2,17(es11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                          \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,7) time,t_tb
WRITE (nprint,9)
WRITE (nprint,11)
WRITE (nprint,101)

e_tot_no_bind(imin:imax)    = grav_pot_a_ray(imin:imax) + ke_a_ray(imin:imax) &
&                           + e_no_bind_a_ray
de_bind_a_ray(imin:imax)    = e_bind_a_ray(imin:imax) - e_bind_fnl_a_ray(imin:imax)

e_tot_no_bind_sum           = zero
DO i = imax, imin, -1
  e_tot_no_bind_sum(i)      = e_tot_no_bind_sum(i+1) + e_tot_no_bind(i)
END DO

de_bind_sum                 = zero
DO i = imax, imin, -1
  de_bind_sum(i)            = de_bind_sum(i+1) + de_bind_a_ray(i)
END DO

e_tot_sum(imin:imax)        = e_tot_no_bind_sum(imin:imax) + de_bind_sum(imin:imax)


DO i = imax, imin, -1
  WRITE (nprint,103) i, x_in(i), c_shock(i), grav_pot_a_ray(i), ke_a_ray(i), &
&  e_no_bind_a_ray(i), e_tot_no_bind(i), e_bind_a_ray(i), de_bind_a_ray(i), &
&  e_tot_no_bind_sum(i), de_bind_sum(i), e_tot_sum(i), dudt_nuc_mean(i), &
&  dudt_nuc_min(i), dudt_nuc_max(i), dudt_nu_mean(i), dudt_nu_min(i), &
&  dudt_nu_max(i), dudt_nuc_zone(i), dudt_nu_zone(i)
END DO ! i

RETURN
END SUBROUTINE edit_Global_energetics
