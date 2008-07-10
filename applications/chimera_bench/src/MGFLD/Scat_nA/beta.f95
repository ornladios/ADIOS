SUBROUTINE beta( e_in, e_out, Z, A, tmev, rho, phi_fm )
!-----------------------------------------------------------------------
!
!    File:         beta
!    Module:       beta
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/23/02
!
!    Purpose:
!      To compute moments of the neutrino-nucleus inelastic scattering
!       kernals.
!
!    Input variables passed through the calling statement
!  e_in            : zone-centered incident neutrino energy (MeV)
!  e_out           : zone-centered scattered neutrino energy (MeV)
!  tmev            : temperature (MeV)
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
!      kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: e_out         ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: tmev          ! temperature (MeV)
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

REAL(KIND=double)                :: e_representative
REAL(KIND=double)                :: delta_e
REAL(KIND=double)                :: beta_up_mgfld
REAL(KIND=double)                :: beta_down_mgfld

e_representative      = ( A/8.d0 ) * tmev
delta_e               = DABS( e_out - e_in )


IF ( e_out >= e_in ) THEN
  phi_fm              = beta_down_mgfld( Z, A, tmev, e_representative, delta_e )
ELSE
  phi_fm              = beta_up_mgfld( Z, A, tmev, e_representative, delta_e )
END IF

RETURN
END SUBROUTINE beta