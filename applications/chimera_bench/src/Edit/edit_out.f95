SUBROUTINE edit_out( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, nedc_ot, &
& nedmi_ot, nedma_ot, nedh_ot, nedps_ot, nedu_ot, nedy_ot, nedsc_ot, &
& nedn_ot, nedng_ot, t_bouncep )
!-----------------------------------------------------------------------
!
!    File:         edit_out
!    Module:       edit_out
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/27/04
!
!    Purpose:
!      To load variables for the mgfld edits.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!  ij_ray_dim  : number of y-zones on a processor before swapping
!  ik_ray_dim  : number of z-zones on a processor before swapping
!
!    Output arguments:
!
!  nedc_ot     : the number of cycles since the last implementation of subsection i of subroutine editc
!  nedmi_ot    : the number of cycles since the last implementation of subsection i of subroutine editmi
!  nedma_ot    : the number of cycles since the last implementation of subsection i of subroutine editma
!  nedh_ot     : the number of cycles since the last implementation of subsection i of subroutine edith
!  nedps_ot    : the number of cycles since the last implementation of subsection i of subroutine editps
!  nedu_ot     : the number of cycles since the last implementation of subsection i of subroutine editu
!  nedy_ot     : the number of cycles since the last implementation of subsection i of subroutine edity
!  nedsc_ot    : the number of cycles since the last implementation of subsection i of subroutine editsc
!  nedn_ot     : the number of cycles since the last implementation of subroutine editn for neutrinos
!                 of type n.
!  nedng_ot    : the number of cycles since the last implementation of subsection i of subroutine 
!                 editng for neutrinos of type n.
!
!    Include files:
!  kind_module, array_module
!  edit_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu

USE edit_module, ONLY : nedc, nedmi, nedma, nedh, nedps, nedu, nedy, nedsc, nedn, nedng
USE t_cntrl_module, ONLY : t_bounce

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                   :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                   :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                   :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)                   :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedc_ot    ! edit counter for subroutine editc
INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedmi_ot   ! edit counter for subroutine editmi
INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedma_ot   ! edit counter for subroutine editma
INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedh_ot    ! edit counter for subroutine edith
INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedps_ot   ! edit counter for subroutine editps
INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedu_ot    ! edit counter for subroutine editu
INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedy_ot    ! edit counter for subroutine edity
INTEGER, INTENT(out), DIMENSION(20,ij_ray_dim,ik_ray_dim)    :: nedsc_ot   ! edit counter for subroutine editsc

INTEGER, INTENT(out), DIMENSION(6,ij_ray_dim,ik_ray_dim)     :: nedn_ot    ! edit counter for subroutine editn
INTEGER, INTENT(out), DIMENSION(100,6,ij_ray_dim,ik_ray_dim) :: nedng_ot   ! edit counter for subroutine editng

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)                :: t_bouncep       ! time of core bounce

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer edit counters to call argument arrays
!-----------------------------------------------------------------------

nedc_ot (1:20,ij_ray,ik_ray)       = nedc (1:20)
nedmi_ot(1:20,ij_ray,ik_ray)       = nedmi(1:20)
nedma_ot(1:20,ij_ray,ik_ray)       = nedma(1:20)
nedh_ot (1:20,ij_ray,ik_ray)       = nedh (1:20)
nedps_ot(1:20,ij_ray,ik_ray)       = nedps(1:20)
nedu_ot (1:20,ij_ray,ik_ray)       = nedu (1:20)
nedy_ot (1:20,ij_ray,ik_ray)       = nedy (1:20)
nedsc_ot(1:20,ij_ray,ik_ray)       = nedsc(1:20)

nedn_ot(1:nnu,ij_ray,ik_ray)       = nedn(1:nnu)
nedng_ot(1:40,1:nnu,ij_ray,ik_ray) = nedng(1:40,1:nnu)

t_bouncep                  = t_bounce

RETURN
END SUBROUTINE edit_out
