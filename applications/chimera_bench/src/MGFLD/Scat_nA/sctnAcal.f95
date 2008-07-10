SUBROUTINE sctnAcal( e_in, e_out, e_in_l, e_in_u, e_out_l, e_out_u, tmev, &
& XH, A, Z, rho, phi_fm )
!-----------------------------------------------------------------------
!
!    File:         sctnAcal
!    Module:       sctnAcal
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/07/02
!
!    Purpose:
!      To compute moments of the neutrino-bucleus inelastic scattering
!       kernals.
!
!    Input variables passed through the calling statement
!  e_in            : zone-centered incident neutrino energy (MeV)
!  e_in_l          : inner boundary of the e_in energy zone (MeV)
!  e_in_u          : outer boundary of the e_in energy zone (MeV)
!  e_out           : zone-centered scattered neutrino energy (MeV)
!  e_out_l         : inner boundary of the e_out energy zone (MeV)
!  e_out_u         : outer boundary of the e_out energy zone (MeV)
!  tmev            : temperature (MeV)
!  XH              : mean heavy nucleus mass fraction
!  A               : mean heavy nucleus mass number
!  Z               : mean heavy nucleus charge number
!  rho             : density (g cm^{-2}
!
!    Output variables passed through the calling statement
!  phi_fm          : zero moment of the Fuller & Meyer (1991) inelastic neutrino-
!                     nucleus scattering functions
!
!    Subprograms called:
!  beta_down_mgfld : computes the exothermic neutrino nucleus scattering kernel
!  beta_up_mgfld   : computes the endothermic neutrino nucleus scattering kernel
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero, half, one
USE physcnst_module, ONLY: rmu

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: e_out         ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: e_in_l        ! inner zone edge of incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: e_in_u        ! outer zone edge of incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: e_out_l       ! inner zone edge of outgoing neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: e_out_u       ! outer zone edge of outgoing neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: tmev          ! temperature (MeV)
REAL(KIND=double), INTENT(in)    :: XH            ! mean heavy nucleus mass fraction
REAL(KIND=double), INTENT(in)    :: A             ! mean heavy nucleus mass number
REAL(KIND=double), INTENT(in)    :: Z             ! mean heavy nucleus charge number
REAL(KIND=double), INTENT(in)    :: rho           ! density (g cm^{-2}

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: phi_fm        ! zero moment of the innelastic neutrino-nucleus scattering function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER, PARAMETER               :: nquad_e =  8    ! number of points of energy Gauss-Lagendre quadrature
INTEGER                          :: i_e             ! summation index of energy Gauss-Lagendre quadrature

REAL(KIND=double), DIMENSION(nquad_e) :: x_e        ! points of energy Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nquad_e) :: wt_e       ! weights of energy Gauss-Lagendre quadrature
REAL(KIND=double)                :: xu_e            ! upper limit of energy integration
REAL(KIND=double)                :: xl_e            ! lower limit of energy integration
REAL(KIND=double)                :: mid_e           ! midpoint of energy integration
REAL(KIND=double)                :: width_e         ! half-width of energy integration
REAL(KIND=double)                :: c_e             ! scaled points of energy quadrature
REAL(KIND=double)                :: e_out_q         ! points of energy quadrature
REAL(KIND=double)                :: phi_e           ! energy integral for a given angle

INTEGER, PARAMETER               :: nquad_ein1 = 4  ! number of points of energy Gauss-Lagendre quadrature
INTEGER                          :: i_ein           ! summation index of energy Gauss-Lagendre quadrature

REAL(KIND=double), DIMENSION(nquad_ein1) :: x_ein1  ! points of energy Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nquad_ein1) :: wt_ein1 ! weights of energy Gauss-Lagendre quadrature
REAL(KIND=double)                :: xu_ein          ! upper limit of energy integration
REAL(KIND=double)                :: xl_ein          ! lower limit of energy integration
REAL(KIND=double)                :: mid_ein         ! midpoint of energy integration
REAL(KIND=double)                :: width_ein       ! half-width of energy integration
REAL(KIND=double)                :: c_ein           ! scaled points of energy quadrature
REAL(KIND=double)                :: e_in_q          ! points of energy quadrature
REAL(KIND=double)                :: phi0_ein        ! energy integral for a given angle
REAL(KIND=double)                :: phi1_ein        ! energy integral for a given angle

INTEGER, PARAMETER               :: nquad_eot1 = 4  ! number of points of energy Gauss-Lagendre quadrature
INTEGER                          :: i_eot           ! summation index of energy Gauss-Lagendre quadrature

REAL(KIND=double), DIMENSION(nquad_eot1) :: x_eot1  ! points of energy Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nquad_eot1) :: wt_eot1 ! weights of energy Gauss-Lagendre quadrature
REAL(KIND=double)                :: xu_eot          ! upper limit of energy integration
REAL(KIND=double)                :: xl_eot          ! lower limit of energy integration
REAL(KIND=double)                :: mid_eot         ! midpoint of energy integration
REAL(KIND=double)                :: width_eot       ! half-width of energy integration
REAL(KIND=double)                :: c_eot           ! scaled points of energy quadrature
REAL(KIND=double)                :: e_ot            ! points of energy quadrature
REAL(KIND=double)                :: phi0_eot        ! energy integral for a given angle
REAL(KIND=double)                :: phi1_eot        ! energy integral for a given angle

INTEGER, PARAMETER               :: nquad_ein2 = 4  ! number of points of energy Gauss-Lagendre quadrature

REAL(KIND=double), DIMENSION(nquad_ein2) :: x_ein2  ! points of energy Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nquad_ein2) :: wt_ein2 ! weights of energy Gauss-Lagendre quadrature

INTEGER, PARAMETER               :: nquad_eot2 = 4  ! number of points of energy Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nquad_eot2) :: x_eot2  ! points of energy Gauss-Lagendre quadrature
REAL(KIND=double), DIMENSION(nquad_eot2) :: wt_eot2 ! weights of energy Gauss-Lagendre quadrature

REAL(KIND=double)                :: e_representative ! mean energy of nucleus
REAL(KIND=double)                :: beta_kernel      ! neutrino-nucleus inelastic pre-scattering kernel
REAL(KIND=double), PARAMETER     :: sigma_zero = 2.584d-44

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                      \\\\\ INITIALIZE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Get quadrature points and weights
!-----------------------------------------------------------------------

IF ( first ) THEN
  first           = .false.
  CALL gquad(nquad_e,x_e,wt_e,nquad_e)
  CALL gquad(nquad_ein1,x_ein1,wt_ein1,nquad_ein1)
  CALL gquad(nquad_eot1,x_eot1,wt_eot1,nquad_eot1)
  CALL gquad(nquad_ein2,x_ein2,wt_ein2,nquad_ein2)
  CALL gquad(nquad_eot2,x_eot2,wt_eot2,nquad_eot2)
END IF

!-----------------------------------------------------------------------
!  Compute e_representative
!-----------------------------------------------------------------------

e_representative  = ( A/8.d0 ) * tmev

!-----------------------------------------------------------------------
!  Initialize scattering kernel
!-----------------------------------------------------------------------

phi_fm            = zero

!-----------------------------------------------------------------------
!
!                 \\\\\ CASE 1 : E_IN = E_OUT /////
!
!-----------------------------------------------------------------------

IF ( e_in == e_out ) THEN

!-----------------------------------------------------------------------
!  Integrate outgoing neutrino energy from e_out_l to e_out_u.
!   This is the inner integration.
!-----------------------------------------------------------------------

  xu_e            = e_out_u
  xl_e            = e_out_l

!-----------------------------------------------------------------------
!  Set integration points and weights
!-----------------------------------------------------------------------

  phi_e           = zero
  mid_e           = half * ( xu_e + xl_e )
  width_e         = half * ( xu_e - xl_e )

  energy: DO i_e = 1,nquad_e
    c_e           = x_e(i_e) * width_e
    e_out_q       = mid_e + c_e

!-----------------------------------------------------------------------
!  Calculate the differential cross section
!-----------------------------------------------------------------------

    CALL beta( e_in, e_out_q, Z, A, tmev, e_representative, beta_kernel )

!-----------------------------------------------------------------------
!  Sum outgoing neutrino energies.
!-----------------------------------------------------------------------

    phi_e         = phi_e + beta_kernel * wt_e(i_e)

  END DO energy

!-----------------------------------------------------------------------
!  Assemble final result
!-----------------------------------------------------------------------

  phi_fm          = XH * ( rho/( A * rmu) ) * sigma_zero * phi_e * width_e &
&                 /( e_out**2 * ( e_in_u - e_in_l ) )

  RETURN

END IF ! e_in == e_out

!-----------------------------------------------------------------------
!
!      \\\\\ CASE 2 : E_IN_L = E_OUT_U OR E_IN_U = E_OUT_L /////
!
!-----------------------------------------------------------------------

IF ( e_in_l == e_out_u  .or.  e_in_l == e_out_l ) THEN

!-----------------------------------------------------------------------
!  Integrate incoming neutrino energy from e_in_l to e_in_u.
!   This is the outer integration.
!-----------------------------------------------------------------------

  phi0_ein        = zero
  xu_ein          = e_in_u
  xl_ein          = e_in_l
  mid_ein         = half * ( xu_ein + xl_ein )
  width_ein       = half * ( xu_ein - xl_ein )

  outer_e: DO i_ein = 1,nquad_ein1

    c_ein         = x_ein1(i_ein) * width_ein 
    e_in_q        = mid_ein + c_ein

!-----------------------------------------------------------------------
!  Integrate outgoing neutrino energy from e_out_l to e_out_u.
!   This is the inner integration.
!-----------------------------------------------------------------------

    phi0_eot      = zero
    xu_eot        = e_out_u
    xl_eot        = e_out_l
    mid_eot       = half * ( xu_eot + xl_eot )
    width_eot     = half * ( xu_eot - xl_eot )

    inner_e: DO i_eot = 1,nquad_eot1

      c_eot       = x_eot1(i_eot) * width_eot
      e_ot        = mid_eot + c_eot

!-----------------------------------------------------------------------
!  Calculate the differential cross section
!-----------------------------------------------------------------------

      CALL beta( e_in_q, e_ot, Z, A, tmev, e_representative, beta_kernel )

!-----------------------------------------------------------------------
!  Sum outgoing neutrino energies.
!-----------------------------------------------------------------------

      phi0_eot    = phi0_eot + beta_kernel * wt_eot1(i_eot) * e_ot**2

    END DO inner_e

!-----------------------------------------------------------------------
!  Sum incoming neutrino energies
!-----------------------------------------------------------------------

    phi0_ein      = phi0_ein + phi0_eot * wt_ein1(i_ein) * width_eot * e_in_q**3

  END DO outer_e

!-----------------------------------------------------------------------
!  Assemble final result
!-----------------------------------------------------------------------

  phi0_ein        = phi0_ein * width_ein

  phi_fm          = XH * ( rho/( A * rmu) ) * sigma_zero * phi0_ein/e_out**2 &
&                 / ( e_in**3  * ( e_in_u  - e_in_l  ) * e_out**2 * ( e_out_u - e_out_l ) )

  RETURN

END IF ! e_in_l == e_out_u

!-----------------------------------------------------------------------
!
!                \\\\\ CASE 3 : E_IN_L > E_OUT_U /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Integrate incoming neutrino energy from e_in_l to e_in_u.
!   This is the outer integration.
!-----------------------------------------------------------------------

phi0_ein          = zero
phi1_ein          = zero
xu_ein            = e_in_u
xl_ein            = e_in_l
mid_ein           = half * ( xu_ein + xl_ein )
width_ein         = half * ( xu_ein - xl_ein )

outer_e2: DO i_ein = 1,nquad_ein2

  c_ein           = x_ein2(i_ein) * width_ein 
  e_in_q          = mid_ein + c_ein

!-----------------------------------------------------------------------
!  Integrate outgoing neutrino energy from e_out_l to e_out_u.
!   This is the inner integration.
!-----------------------------------------------------------------------

  phi0_eot        = zero
  phi1_eot        = zero
  xu_eot          = e_out_u
  xl_eot          = e_out_l
  mid_eot         = half * ( xu_eot + xl_eot )
  width_eot       = half * ( xu_eot - xl_eot )

  inner_e2: DO i_eot = 1,nquad_eot2

    c_eot         = x_eot2(i_eot) * width_eot
    e_ot          = mid_eot + c_eot
      
!-----------------------------------------------------------------------
!  Calculate the differential cross section
!-----------------------------------------------------------------------

    CALL beta( e_in_q, e_ot, Z, A, tmev, e_representative, beta_kernel )

!-----------------------------------------------------------------------
!  Sum outgoing neutrino energies.
!-----------------------------------------------------------------------

    phi0_eot      = phi0_eot + beta_kernel * wt_eot2(i_eot) * e_ot**2

  END DO inner_e2

!-----------------------------------------------------------------------
!  Sum incoming neutrino energies
!-----------------------------------------------------------------------

  phi0_ein        = phi0_ein + phi0_eot * wt_ein2(i_ein) * width_eot * e_in_q**3

END DO outer_e2

!-----------------------------------------------------------------------
!  Assemble final result
!-----------------------------------------------------------------------

phi0_ein          = phi0_ein * width_ein

phi_fm            = XH * ( rho/( A * rmu) ) * sigma_zero * phi0_ein/e_out**2  &
&                 / ( e_in**3  * ( e_in_u  - e_in_l  ) * e_out**2 * ( e_out_u - e_out_l ) )

RETURN
END SUBROUTINE sctnAcal
