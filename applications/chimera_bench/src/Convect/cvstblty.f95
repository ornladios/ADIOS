SUBROUTINE cvstblty( j, jr_min, jr_max, ij_ray, ik_ray, sch, ldx, sfr )
!-----------------------------------------------------------------------
!
!    File:         cvstblty
!    Module:       cvstblty
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/17/95
!
!    Purpose:
!      To determine the local stability of two adjacent zones
!       against convection using Ledoux, Schwarzchild, and Thorn
!       and Wilson's Neutron Finger criteria.
!
!    Subprograms called:
!  tdgvnspye_x : computes the temperature given s, p, and ye
!
!    Input arguments:
!  j           : radial zone for which neutrino occupation adjustment
!                 is to be computed.
!  jr_min      : minimum radial zone index
!  jr_max      : maximum radial zone index
!  ij_ray      : index denoting the j-index of a specific radial ray
!  ik_ray      : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!  sch         : 'yes': unstable to Schwarzschild convection
!                'no' : stable to Schwarzschild convection
!  ldx         : 'yes': unstable to Ledoux convection
!                'no' : stable to Ledoux convection
!  sfr         : 'yes': unstable to 'salt-fingers'
!                'no' : stable to 'salt-fingers'
!
!    Include files:
!  kind_module, numerical_module
!  eos_snc_x_module, mdl_cnfg_module, shock_module
!
!-----------------------------------------------------------------------

USE kind_module

USE eos_snc_x_module, ONLY : aesv, nse
USE mdl_cnfg_module, ONLY : rho, t, ye
USE shock_module, ONLY : pqcrit

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

 CHARACTER(LEN=1), INTENT(out)   :: sch           !
 CHARACTER(LEN=1), INTENT(out)   :: ldx
 CHARACTER(LEN=1), INTENT(out)   :: sfr

INTEGER                          :: jp            ! radial zone index to be passed
INTEGER                          :: jshockmx      ! outer radial index of shock
INTEGER                          :: jshockmn      ! inner radial index of shock
INTEGER, PARAMETER               :: itmax = 40    ! radial zone index
INTEGER, PARAMETER               :: itrgn = 20    ! radial zone index

REAL(KIND=double)                :: pqmin         ! value of pq_x to test for the presence of a shock
REAL(KIND=double), PARAMETER     :: tol = 1.d-5   ! zone centered incoming neutrino energy (MeV)

REAL(KIND=double)                :: yeblob        ! electron fraction of blob
REAL(KIND=double)                :: pblob         ! pressure of blob
REAL(KIND=double)                :: sblob         ! entropy of blob
REAL(KIND=double)                :: rhoblob       ! density of blob (determined by iteration)
REAL(KIND=double)                :: tblob         ! temperature of blob (determined by iteration)
REAL(KIND=double)                :: rho_prev      ! density of blob guess
REAL(KIND=double)                :: t_prev        ! temperature of blob guess

!-----------------------------------------------------------------------
!  Return if nse(j,i_ray) = 0
!-----------------------------------------------------------------------

IF ( nse(j,ij_ray,ik_ray) .eq. 0 ) THEN
  sch              = 'u'
  ldx              = 'u'
  sfr              = 'u'
  RETURN
END IF

!-----------------------------------------------------------------------
!  Return if above shock
!-----------------------------------------------------------------------

pqmin              = pqcrit
CALL findshock_min( jr_min, jr_max, ij_ray, ik_ray, pqmin, jshockmn, jshockmx)
IF ( j >= jshockmn ) THEN
  sch              = 'u'
  ldx              = 'u'
  sfr              = 'u'
  RETURN
END IF

!-----------------------------------------------------------------------
!
!                 \\\\\ SCHWARZSCHILD STABILITY /////
!
!-----------------------------------------------------------------------

yeblob             = ye(j+1)
pblob              = aesv(j+1,1,ij_ray,ik_ray)
sblob              = aesv(j  ,3,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Determine rhoblob and tblob of a blob displaced at constant entropy
!   but equilibrated with the background in pressure and composition
!-----------------------------------------------------------------------

jp                 = 1
rho_prev           = rho(j)
t_prev             = t(j)

CALL tdgvnspye_x( jp, ij_ray, ik_ray, sblob, pblob, yeblob, rho_prev, &
& t_prev, rhoblob, tblob )

IF ( rhoblob < rho(j+1) ) THEN
  sch              = 'y'
ELSE
  sch              = 'n'
END IF ! rhoblob < rho(j+1)

!-----------------------------------------------------------------------
!
!                    \\\\\ LEDOUX STABILITY /////
!
!-----------------------------------------------------------------------

yeblob             = ye(j)
pblob              = aesv(j+1,1,ij_ray,ik_ray)
sblob              = aesv(j  ,3,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Determine rhoblob and tblob of a blob displaced at constant entropy
!   and composition but equilibrated the background in pressure and
!   composition
!-----------------------------------------------------------------------

jp                 = 1
rho_prev           = rho(j)
t_prev             = t(j)

CALL tdgvnspye_x( jp, ij_ray, ik_ray, sblob, pblob, yeblob, rho_prev, &
& t_prev, rhoblob, tblob )

IF ( rhoblob < rho(j+1) ) THEN
  ldx              = 'y'
ELSE
  ldx              = 'n'
END IF ! rhoblob < rho(j+1)

!-----------------------------------------------------------------------
!
!                 \\\\\ SALT-FINGER STABILITY /////
!
!-----------------------------------------------------------------------

tblob              = t(j+1)
yeblob             = ye(j)
pblob              = aesv(j+1,1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Determine rhoblob of a blob displaced at constant composition but
!   equilibrated with the background in pressure and temperature
!-----------------------------------------------------------------------

jp                 = 1
rho_prev           = rho(j)

CALL dgvntpye_x( jp, ij_ray, ik_ray, tblob, pblob, yeblob, rho_prev, rhoblob )

IF ( rhoblob < rho(j+1) ) THEN
  sfr         = 'y'
ELSE
  sfr         = 'n'
END IF


RETURN
END SUBROUTINE cvstblty
