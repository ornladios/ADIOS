SUBROUTINE tgvndeye_sweep_x( ni, nf, ij_ray, ik_ray, rho_f, rho_i )
!-----------------------------------------------------------------------
!
!    File:         tgvndeye_sweep_x
!    Module:       tgvndeye_sweep_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/29/03
!
!    Purpose:
!      To compute the temperature given rho, e, and ye. Iteration
!       is by means of the bisection method if an initial guess
!       of the temperature is unavailable or is not within fraction
!       of the previous value, otherwise by Newton-Rhapson. Once
!       the temperature is found, the pressure entropy, and the
!       adiabatic exponents are computed.
!
!    Variables that must be passed through common:
!  padded atate variables
!
!    Subprograms called:
!  e_degen    : determines the electron degeneracy parameter
!  tgvndsye_x : unpdates the temperature by iterating on the entropy
!  tgvndeye_x : unpdates the temperature by iterating on the energy
!  eqstt_x    : interpolates EOS quantities from cube corners to a given state
!
!    Input arguments:
!  ni         : minimum paddded array index
!  nf         : maximum padded array index
!  ij_ray     : n-index of a radial ray
!  ik_ray     : k-index of a radial ray
!  rho_f      : final density (for the degenerate case)
!  rho_i      : initial density (for the degenerate case)
!
!    Output arguments:
!  t          : temperature (K)
!
!    Include files:
!  kind_module, array_module, physcnst_module
!  edit_module, evh1_global, evh1_global, evh1_sweep,
!  parallel_module
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : max_12, n_proc_y, ij_ray_dim, ik_ray_dim
USE physcnst_module, ONLY : mbary=>rmu, ergmev, bok=>kmev

USE edit_module, ONLY : nlog
USE evh1_global, ONLY : degen
USE evh1_sweep, ONLY : re=>r, yee=>ye, we=>w, se=>entrop, temp, ei, e_nu, &
& p_mat, p_nu, p, gc, ge, l_shock
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: ni       ! minimum zone index
INTEGER, INTENT(in)                    :: nf       ! maximum zone index
INTEGER, INTENT(in)                    :: ij_ray   ! j-index of a radial ray
INTEGER, INTENT(in)                    :: ik_ray   ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(max_12) :: rho_f    ! final rho (padded zones)
REAL(KIND=double), INTENT(in), DIMENSION(max_12) :: rho_i    ! initial rho (padded zones)

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                :: n        ! padded zone index
INTEGER                                :: j_ray    ! polar index of the radial ray
INTEGER                                :: k_ray    ! azimuthal index of the radial ray

REAL(KIND=double)                      :: t_out    ! temperature given the internal energy
REAL(KIND=double)                      :: dummy1   ! dummy eos argument
REAL(KIND=double)                      :: dummy2   ! dummy eos argument
REAL(KIND=double)                      :: dummy3   ! dummy eos argument
REAL(KIND=double)                      :: gam_s    ! gamma1
REAL(KIND=double)                      :: e_eta    ! (electron Fermi energy)/KT

!-----------------------------------------------------------------------
!        Formats.
!-----------------------------------------------------------------------

  101 FORMAT (' Warning: Small gamma ',es11.3,' reset to ',es11.3, &
&  ' in zone ',i4,' radial ray',2i4,' proc', i4,' in tgvndeye_sweep_x; re(n)=', &
&  es11.3,' temp(n)=', es11.3,' yee(n)=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                \\\\\ BEGIN RADIAL ZONE LOOP /////
!
!-----------------------------------------------------------------------

DO n = ni, nf

!-----------------------------------------------------------------------
!  Get degeneracy parameter
!-----------------------------------------------------------------------

  CALL e_degen( re(n), temp(n), yee(n), e_eta )

  IF ( e_eta > degen  .and.  ( .not. l_shock(n) ) ) THEN

!-----------------------------------------------------------------------
!  Iterate on entropy to update temperature for high electron
!   degeneracy (i.e., if e_eta > degen)
!-----------------------------------------------------------------------

    CALL tgvndsye_x( n-5, ij_ray, ik_ray, re(n), se(n), yee(n), temp(n), &
&    t_out )
    temp(n)         = t_out

  ELSE ! e_eta <= degen

!-----------------------------------------------------------------------
!  Iterate on energy to update temperature for low electron
!   degeneracy (i.e., if e_eta < degen)
!-----------------------------------------------------------------------

    CALL tgvndeye_x( n-5, ij_ray, ik_ray, re(n), ei(n), yee(n), temp(n), &
&    t_out )
    temp(n)         = t_out
  
  END IF ! e_eta > degen

!-----------------------------------------------------------------------
!  Get the new pressure, entropy, and gamma
!-----------------------------------------------------------------------

  CALL eqstt_x( 1, n-5, ij_ray, ik_ray, re(n), temp(n), yee(n), p_mat(n),   &
&  dummy1, dummy2, dummy3 )

  CALL eqstt_x( 3, n-5, ij_ray, ik_ray, re(n), temp(n), yee(n), se(n),  &
&  dummy1, dummy2, dummy3 )

  CALL eqstt_x( 12, n-5, ij_ray, ik_ray, re(n), temp(n), yee(n), gam_s, &
&  dummy1, dummy2, dummy3 )

!-----------------------------------------------------------------------
!  Place a lower bound on gc(n)
!-----------------------------------------------------------------------

  IF ( gam_s < 0.01d0 ) THEN
    gc(n)           = .01d0
    j_ray          = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
    k_ray          = ( myid/n_proc_y ) * ik_ray_dim + ik_ray
    WRITE (6,101) gam_s, gc(n), n, j_ray, k_ray, myid, re(n), temp(n), yee(n)
    WRITE (nlog,101) gam_s, gc(n), n, j_ray, k_ray, myid, re(n), temp(n), yee(n)
  ELSE
    gc(n)           = gam_s
  END IF ! gam_s < 0.01d0

!-----------------------------------------------------------------------
!  Compute sum of matter pressure and neutrino pressure
!-----------------------------------------------------------------------

  p(n)              = p_mat(n) + p_nu(n)

!-----------------------------------------------------------------------
!  Compute ge(n)
!-----------------------------------------------------------------------

  ge(n)             = 1.0d0 + p(n)/( ( ei(n) + e_nu(n) ) * re(n) )

!-----------------------------------------------------------------------
!
!                 \\\\\ END RADIAL ZONE LOOP /////
!
!-----------------------------------------------------------------------

END DO ! n = ni,nf

RETURN
END SUBROUTINE tgvndeye_sweep_x
