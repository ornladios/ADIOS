SUBROUTINE e_zone
!-----------------------------------------------------------------------
!
!    File:         e_zone_e
!    Module:       e_zone
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/8/99
!
!    Purpose:
!      To compute the energy zoning for the neutrinos.
!
!    Variables that must be passed through common:
!  agr(j)      : GR lapse function
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  nnugp(n)    : number of neutrino energy groups for n-neutrinos
!  unumn       : minimum neutrino zone centered energy at infinity
!  unumx       : maximum neutrino zone centered energy at infinity
!
!    Output arguments (common):
!  unubi(k)    : energy of outer edge of neutrino energy zone k as measured by an observer at infinity
!  unui(k)     : zone centered energy of neutrino energy zone k as measured by an observer at infinity
!  dunui(k)    : unubi(k+1) - unubi(k)
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, nu_dist_module, nu_energy_grid_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half

USE edit_module, ONLY : nprint, nlog
USE nu_dist_module, ONLY : runu, unumn, unumx
USE nu_energy_grid_module, ONLY : nnugp, unui, unubi, dunui, nnugpmx
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: k             ! energy zone index
INTEGER                          :: np            ! nnugpmx
INTEGER                          :: npm           ! nnugpmx-1
INTEGER                          :: npp           ! nnugpmx+1

REAL(KIND=double)                :: rngpm1        ! nnugpmx-1
REAL(KIND=double)                :: runul         ! log(unumx/unumn)/rngpm1
REAL(KIND=double)                :: runuh         ! sqrt(runul)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 4001 format (1x)
 4003 format (10x,'Neutrino group energies')
 4005 format (10x,32('-')/)
 4007 format (' unui(',i2,')=',es11.4,10x,' unubi(',i2,')=',es11.4,10x, 'dunui(',i2,')=',es11.4)
 4009 format (31x,' unubi(',i2,')=',es11.4)
 4011 format (' runu=',es11.4//)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!               \\\\\ NEUTRINO GRPUP ENERGIES /////
!
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

np                   = nnugpmx
npm                  = nnugpmx - 1
npp                  = nnugpmx + 1

IF ( unui(1) == zero ) THEN

  IF ( nnugpmx == 1 ) THEN

!-----------------------------------------------------------------------
!  One energy group
!-----------------------------------------------------------------------

    unubi(1)         = unumn
    unubi(2)         = unumx
    unui(1)          = half * ( unumn + unumx )

!-----------------------------------------------------------------------
!  More than one energy group
!-----------------------------------------------------------------------

  ELSE ! nnugpmx > 1

    rngpm1           = DBLE(nnugpmx-1)                                           

    runul            = DLOG(unumx/unumn)/rngpm1
    runu             = DEXP(runul)
    runuh            = DSQRT(runu)

    unui(1)          = unumn
    unubi(1)         = unumn/runuh

    DO k = 2,nnugpmx
      unui(k)        = runu * unui(k-1)
      unubi(k)       = runu * unubi(k-1)
    END DO

    unubi(npp)       = runu * unubi(np)

!-----------------------------------------------------------------------
!  Choose zone-centered energy to reproduce density of states
!-----------------------------------------------------------------------

    unubi(1)         = zero
    unui(1:nnugpmx)  = DSQRT(unubi(1:nnugpmx)**2 + unubi(1:nnugpmx) * unubi(2:nnugpmx+1) + unubi(2:nnugpmx+1)**2 )/DSQRT(3.d0)

  END IF ! nnugpmx > 1

END IF ! unui(1) = 0

!-----------------------------------------------------------------------
!  Weights for the energy group summation
!-----------------------------------------------------------------------

 dunui(1:nnugpmx)    = unubi(2:nnugpmx+1) - unubi(1:nnugpmx)

!-----------------------------------------------------------------------
!  Print energy group data
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN
  WRITE (nprint,4001)
  WRITE (nprint,4003)
  WRITE (nprint,4005)
  WRITE (nprint,4007) (k,unui(k),k,unubi(k),k,dunui(k),k=1,nnugpmx)
  WRITE (nprint,4009) nnugpmx+1,unubi(nnugpmx+1)
  WRITE (nprint,4011) runu
END IF

WRITE (nlog,4001)
WRITE (nlog,4003)
WRITE (nlog,4005)
WRITE (nlog,4007) (k,unui(k),k,unubi(k),k,dunui(k),k=1,nnugpmx)
WRITE (nlog,4009) nnugpmx+1,unubi(nnugpmx+1)
WRITE (nlog,4011) runu

RETURN
END SUBROUTINE e_zone
