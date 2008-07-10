SUBROUTINE mfp_cal( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         mfp_cal
!    Module:       mfp_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/22/00
!
!    Purpose:
!      To compute the transport mean free paths.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min     : inner radial zone for computing the diffusion coefficient.
!  jr_max     : outer radial zone for computing the diffusion coefficient.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Modules:
!  kind_module, numerical_module, array_module
!  abem_module, brem_module, nu_dist_module, nu_energy_grid_module,
!  pair_module, scat_a_module, scat_e_module, scat_i_module,
!  scat_n_module, scat_nA_module, scat_nn_module
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : one, epsilon
USE array_module, ONLY : nnu

USE abem_module, ONLY : emis, absor
USE brem_module, ONLY : baf
USE nu_dist_module, ONLY : rhs1, tmfp, e_mfp_inv
USE nu_energy_grid_module, ONLY : nnugp
USE pair_module, ONLY : paf
USE scat_a_module, ONLY : scaf
USE scat_e_module, ONLY : scef
USE scat_i_module, ONLY : scti
USE scat_n_module, ONLY : scnf
USE scat_nn_module, ONLY : scnnf
USE scat_nA_module, ONLY : scnAf

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: jr_min     ! minimum radial zone index
INTEGER, INTENT(in)            :: jr_max     ! maximum radial zone index
INTEGER, INTENT(in)            :: ij_ray     ! j-index of a radial ray
INTEGER, INTENT(in)            :: ik_ray     ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                        :: j          ! radial zone index
INTEGER                        :: k          ! neutrino energy index
INTEGER                        :: n          ! neutrino flavor index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       Compute transport and inverse transport mean free paths.
!
!        1/lambda^t = 1/tmfp = - rhs1
!-----------------------------------------------------------------------
 
DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max

      rhs1(j,k,n)      = - emis(j,k,n,ij_ray,ik_ray)                    &
&                      - absor(j,k,n,ij_ray,ik_ray)                     &
&                      - scti(j,k,n)                                    &
&                      + scaf(5,j,k,n) + scef(5,j,k,n)                  &
&                      + paf(5,j,k,n)  + baf(5,j,k,n)                   &
&                      + scnf(5,j,k,n) + scnnf(5,j,k,n) + scnAf(5,j,k,n)

      tmfp(j,k,n)      = one/( - rhs1(j,k,n) + epsilon )

    END DO

  tmfp(jr_max+1,k,n)     = tmfp(jr_max,k,n)

  END DO
END DO

!-----------------------------------------------------------------------
!       Compute inverse energy mean free path.
!
!        1/lambda^t = 1/tmfp = - rhs1
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max

      e_mfp_inv(j,k,n) = emis(j,k,n,ij_ray,ik_ray)                      &
&                      + absor(j,k,n,ij_ray,ik_ray)                     &
&                      + scaf(5,j,k,n) + scef(5,j,k,n)                  &
&                      + paf(5,j,k,n)  + baf(5,j,k,n)                   &
&                      + scnf(5,j,k,n) + scnnf(5,j,k,n) + scnAf(5,j,k,n)
    END DO

    e_mfp_inv(jr_max+1,k,n) = e_mfp_inv(jr_max,k,n)

  END DO
END DO

RETURN
END
