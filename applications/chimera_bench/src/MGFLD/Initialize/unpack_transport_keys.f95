SUBROUTINE unpack_transport_keys( nez, nezp1, nnu, i_trans_data, d_trans_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_transport_keys
!    Module:       unpack_transport_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/23/04
!
!    Purpose:
!      To unpack the transport key arrays and restore the values
!       to the appropriate variables in the appropriate modules.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nez          : energy array dimension
!  nez+1        : nez + 1
!  nnu          : neutrino flavor dimension
!  i_trans_data : integer array of transport keys
!  d_trans_data : real*8 array of transport keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  cycle_module, eos_snc_x_module, it_tol_module, nu_dist_module,
!  nu_energy_grid_module, prb_cntl_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE cycle_module, ONLY : intnu_trns, ncynu_trns, intnur_trns
USE eos_snc_x_module, ONLY : rhopnu
USE it_tol_module, ONLY : tolnut, tolnuye, tolnupsi, tolpsimin, iternu, &
& itfail, a_prec
USE nu_dist_module, ONLY : unumn, unumx, stwt
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, unubi, nnugpmx
USE prb_cntl_module, ONLY : iyenu, itnu, jnumin, jnumax, idiff, inutrn, &
& ireltrns, iaence, edmpe, iaenca, edmpa, iaefnp, rhoaefnp, iaencnu, &
& roaencnu, i_aeps, icc_wm, iscat, in, ip, ihe, iheavy, iicor, ietann, &
& inc_wm, nes, rhonesmn, rhonesmx, nncs, rhonncsmn, rhonncsmx, ipair, &
& rhopairemn, rhopairemx, rhopairtmn, rhopairtmx, ibrem, rhobrememn, &
& rhobrememx, rhobremtmn, rhobremtmx, isctn, rhosctnemn, rhosctnemx, &
& rhosctntmn, rhosctntmx, isctnn, rhosctnnemn, rhosctnnemx, rhosctnntmn, &
& rhosctnntmx, isctnA, rhosctnAemn, rhosctnAemx, rhosctnAtmn, rhosctnAtmx, &
& iaenct, roaenct
USE t_cntrl_module, ONLY : dtnph_trans, dtnmhn_trans, psimin, psipmin, &
& tcntrl, rdtmax

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nez          ! energy array dimension
INTEGER, INTENT(in)              :: nezp1        ! nez + 1
INTEGER, INTENT(in)              :: nnu          ! neutrino flavor dimension

INTEGER, INTENT(in), DIMENSION(40+2*nnu)                      :: i_trans_data  ! integer array of transport keys

REAL(KIND=double), INTENT(in), DIMENSION((110+3*nnu+3*nez+1)) :: d_trans_data  ! 64 bit real array of transport keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i             ! do i ndex
INTEGER                          :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                \\\\\ UNPACK TRANSPORT KEYS /////
!
!-----------------------------------------------------------------------

iyenu                   = i_trans_data(1)
jnumin                  = i_trans_data(2)
jnumax                  = i_trans_data(3)
idiff                   = i_trans_data(4)
inutrn                  = i_trans_data(5)
ireltrns                = i_trans_data(6)
iternu                  = i_trans_data(7)
itfail                  = i_trans_data(8)
intnu_trns              = i_trans_data(9)
ncynu_trns              = i_trans_data(11)
intnur_trns             = i_trans_data(13)
iaence                  = i_trans_data(15)
iaenca                  = i_trans_data(16)
iaefnp                  = i_trans_data(17)
iaencnu                 = i_trans_data(18)
i_aeps                  = i_trans_data(19)
icc_wm                  = i_trans_data(20)
iscat                   = i_trans_data(21)
in                      = i_trans_data(22)
ip                      = i_trans_data(23)
ietann                  = i_trans_data(24)
ihe                     = i_trans_data(25)
iheavy                  = i_trans_data(26)
iicor                   = i_trans_data(27)
inc_wm                  = i_trans_data(28)
nes                     = i_trans_data(29)
nncs                    = i_trans_data(30)
ipair                   = i_trans_data(31)
ibrem                   = i_trans_data(32)
isctn                   = i_trans_data(33)
isctnn                  = i_trans_data(34)
isctnA                  = i_trans_data(35)
iaenct                  = i_trans_data(36)

DO i = 1,nnu
  itnu(i)               = i_trans_data(40+0*nnu+i)
END DO

DO i = 1,nnu
  nnugp(i)              = i_trans_data(40+1*nnu+i)
  IF ( inutrn == 0 ) nnugp(i) = 0
END DO

nnugpmx                 = 0
DO n = 1,nnu
  nnugpmx               = MAX( nnugpmx, nnugp(n) )
END DO

tolnut                  = d_trans_data(1)
tolnuye                 = d_trans_data(2)
tolnupsi                = d_trans_data(3)
tolpsimin               = d_trans_data(4)
a_prec                  = d_trans_data(5)
dtnph_trans             = d_trans_data(13)
dtnmhn_trans            = d_trans_data(15)
rdtmax                  = d_trans_data(19)
unumn                   = d_trans_data(24)
unumx                   = d_trans_data(25)
edmpe                   = d_trans_data(26)
edmpa                   = d_trans_data(27)
rhoaefnp                = d_trans_data(28)
roaencnu                = d_trans_data(29)
rhonesmn                = d_trans_data(30)
rhonesmx                = d_trans_data(31)
rhonncsmn               = d_trans_data(32)
rhonncsmx               = d_trans_data(33)
rhopairemn              = d_trans_data(34)
rhopairemx              = d_trans_data(35)
rhopairtmn              = d_trans_data(36)
rhopairtmx              = d_trans_data(37)
rhosctnemn              = d_trans_data(38)
rhosctnemx              = d_trans_data(39)
rhosctntmn              = d_trans_data(40)
rhosctntmx              = d_trans_data(41)
rhosctnnemn             = d_trans_data(42)
rhosctnnemx             = d_trans_data(43)
rhosctnntmn             = d_trans_data(44)
rhosctnntmx             = d_trans_data(45)
rhosctnAemn             = d_trans_data(46)
rhosctnAemx             = d_trans_data(47)
rhosctnAtmn             = d_trans_data(48)
rhosctnAtmx             = d_trans_data(49)
rhopnu                  = d_trans_data(50)
roaenct                 = d_trans_data(51)

DO i = 11,40
  tcntrl(i)             = d_trans_data(50 +0*nnu+0*nez+i)
END DO

tcntrl(50)              = d_trans_data(50 +0*nnu+0*nez+50)

DO i = 1,nnu
  psimin(i)             = d_trans_data(100+0*nnu+0*nez+i)
END DO

DO i = 1,nnu
  psipmin(i)            = d_trans_data(100+1*nnu+0*nez+i)
END DO

DO i = 1,nnu
  stwt(i)               = d_trans_data(100+2*nnu+0*nez+i)
END DO

DO i = 1,nez
  unui(i)               = d_trans_data(100+3*nnu+0*nez+i)
END DO

DO i = 1,nez
  dunui(i)              = d_trans_data(100+3*nnu+1*nez+i)
END DO

DO i = 1,nez+1
  unubi(i)              = d_trans_data(100+3*nnu+2*nez+i)
END DO

END SUBROUTINE unpack_transport_keys
