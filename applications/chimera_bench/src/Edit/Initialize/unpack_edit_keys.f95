SUBROUTINE upack_edit_keys( ij_ray_dim, ik_ray_dim, nez, nnu, i_edit_data, &
& d_edit_data )
!-----------------------------------------------------------------------
!
!    File:         upack_edit_keys
!    Module:       upack_edit_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/23/04
!
!    Purpose:
!      To unpack the edit key arrays and restore the values
!       to the appropriate variables in edit_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray_dim  : number of y-zones on a processor before swapping with y
!  ik_ray_dim  : number of z-zones on a processor before swapping with z
!  nez         : energy array dimension
!  nnu         : neutrino flavor dimension
!  i_edit_data : integer array of edit keys
!  d_edit_data : real*8 array of edit keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  edit_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE edit_module
USE radial_ray_module, ONLY : nedc_c=>nedc, nede_c=>nede, nedmi_c=>nedmi, &
& nedma_c=>nedma, nedh_c=>nedh, nedps_c=>nedps, nedu_c=>nedu, nedy_c=>nedy, &
& nedsc_c=>nedsc, nedn_c=>nedn, nedng_c=>nedng

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray_dim   ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim   ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nez          ! energy array dimension
INTEGER, INTENT(in)              :: nnu          ! neutrino flavor dimension

INTEGER, INTENT(in), DIMENSION(1200+3*40*nnu) :: i_edit_data  ! integer array of edit keys

REAL(KIND=double), INTENT(in), DIMENSION(50)  :: d_edit_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ij_ray       ! index denoting the j-index of a specific radial ray
INTEGER                          :: ik_ray       ! index denoting the k-index of a specific radial ray
INTEGER                          :: i            ! do i ndex
INTEGER                          :: k            ! neutrino eneregy index
INTEGER                          :: n            ! neutrino flavor index

!-----------------------------------------------------------------------
!
!                   \\\\\ UNPACK EDIT KEYS /////
!
!-----------------------------------------------------------------------

iflprt                             = i_edit_data(1)
noutfl                             = i_edit_data(2)
iprint                             = i_edit_data(3)
nrstd1                             = i_edit_data(4)
nrstd2                             = i_edit_data(5)
noutpmt                            = i_edit_data(6)
nouttmp                            = i_edit_data(7)
iplot                              = i_edit_data(8)
nplotc                             = i_edit_data(9)
nplote                             = i_edit_data(10)
nplota                             = i_edit_data(11)
nplott                             = i_edit_data(12)
intplf                             = i_edit_data(13)
npltf                              = i_edit_data(14)
nmodel(:,:)                        = i_edit_data(15)
nprt  (:,:)                        = i_edit_data(17)
intprt                             = i_edit_data(16)
iprtbgn                            = i_edit_data(18)
intrst                             = i_edit_data(19)
nnrst                              = i_edit_data(20)
nrstfl                             = i_edit_data(21)
irstbgn                            = i_edit_data(22)
ncychg                             = i_edit_data(23)
intprm                             = i_edit_data(24)
nprm                               = i_edit_data(25)
ilumplt                            = i_edit_data(26)
nlumplt1                           = i_edit_data(27)
nlumplt2                           = i_edit_data(28)
nlum                               = i_edit_data(29)
intlum(1)                          = i_edit_data(30)
intlum(2)                          = i_edit_data(31)
intlum(3)                          = i_edit_data(32)
ncylum(1)                          = i_edit_data(33)
ncylum(2)                          = i_edit_data(34)
ienuplt                            = i_edit_data(35)
nenuplt1                           = i_edit_data(36)
nenuplt2                           = i_edit_data(37)
nenu                               = i_edit_data(38)
intenu(1)                          = i_edit_data(39)
intenu(2)                          = i_edit_data(40)
intenu(3)                          = i_edit_data(41)
ncyenu(1)                          = i_edit_data(42)
ncyenu(2)                          = i_edit_data(43)
irnuplt                            = i_edit_data(44)
nrnuplt                            = i_edit_data(45)
nrnu                               = i_edit_data(46)
intrnu(1)                          = i_edit_data(47)
intrnu(2)                          = i_edit_data(48)
intrnu(3)                          = i_edit_data(49)
ncyrnu(1)                          = i_edit_data(50)
ncyrnu(1)                          = i_edit_data(51)
ivarplt                            = i_edit_data(52)
nvar                               = i_edit_data(53)
nvarint                            = i_edit_data(54)
nvarplt                            = i_edit_data(55)
nvardump                           = i_edit_data(56)
icomplt                            = i_edit_data(57)
ncomplt                            = i_edit_data(58)
ncomdump                           = i_edit_data(59)
nplotinnerb                        = i_edit_data(60)
iplotinnerb                        = i_edit_data(61)
nplotouterb                        = i_edit_data(62)
iplotouterb                        = i_edit_data(63)
nplotlum                           = i_edit_data(64)
iplotlum                           = i_edit_data(65)
nplotshk                           = i_edit_data(66)
iplotshk                           = i_edit_data(67)
nplotcnv                           = i_edit_data(68)
iplotcnv                           = i_edit_data(69)
nplotnurad                         = i_edit_data(70)
iplotnurad                         = i_edit_data(71)
nnudata                            = i_edit_data(72)
inudata                            = i_edit_data(73)
nlagplt                            = i_edit_data(74)
ilagplt                            = i_edit_data(75)
nlagdump                           = i_edit_data(76)
nrlagplt                           = i_edit_data(77)
irlagplt                           = i_edit_data(78)
n_eplt                             = i_edit_data(79)
i_eplt                             = i_edit_data(80)
iedMDu                             = i_edit_data(81)
nedMDu                             = i_edit_data(82)
intedMDu                           = i_edit_data(83)
n_editMDu                          = i_edit_data(84)
iedMDv                             = i_edit_data(85)
nedMDv                             = i_edit_data(86)
intedMDv                           = i_edit_data(87)
n_editMDv                          = i_edit_data(88)
iedMDw                             = i_edit_data(89)
nedMDw                             = i_edit_data(90)
intedMDw                           = i_edit_data(91)
n_editMDw                          = i_edit_data(92)
iedMDs                             = i_edit_data(93)
nedMDs                             = i_edit_data(94)
intedMDs                           = i_edit_data(95)
n_editMDs                          = i_edit_data(96)
iedMDd                             = i_edit_data(97)
nedMDd                             = i_edit_data(98)
intedMDd                           = i_edit_data(99)
n_editMDd                          = i_edit_data(100)
iedMDe                             = i_edit_data(101)
nedMDe                             = i_edit_data(102)
intedMDe                           = i_edit_data(103)
n_editMDe                          = i_edit_data(104)
iedMDp                             = i_edit_data(105)
nedMDp                             = i_edit_data(106)
intedMDp                           = i_edit_data(107)
n_editMDp                          = i_edit_data(108)
iedMDenu                           = i_edit_data(109)
nedMDenu                           = i_edit_data(110)
intedMDenu                         = i_edit_data(111)
n_editMDenu                        = i_edit_data(112)
iedMDfnu                           = i_edit_data(113)
nedMDfnu                           = i_edit_data(114)
intedMDfnu                         = i_edit_data(115)
n_editMDfnu                        = i_edit_data(116)
iedMDa                             = i_edit_data(117)
nedMDa                             = i_edit_data(118)
intedMDa                           = i_edit_data(119)
n_editMDa                          = i_edit_data(120)
iedMDx                             = i_edit_data(121)
nedMDx                             = i_edit_data(122)
intedMDx                           = i_edit_data(123)
n_editMDx                          = i_edit_data(124)
iedMDye                            = i_edit_data(125)
nedMDye                            = i_edit_data(126)
intedMDye                          = i_edit_data(127)
n_editMDye                         = i_edit_data(128)
iedMDcm                            = i_edit_data(129)
nedMDcm                            = i_edit_data(130)
intedMDcm                          = i_edit_data(131)
n_editMDcm                         = i_edit_data(132)
iedMDnu                            = i_edit_data(133)
nedMDnu                            = i_edit_data(134)
intedMDnu                          = i_edit_data(135)
n_editMDnu                         = i_edit_data(136)
iedMDnc                            = i_edit_data(137)
nedMDnc                            = i_edit_data(138)
intedMDnc                          = i_edit_data(139)
n_editMDnc                         = i_edit_data(140)
iedMDnl                            = i_edit_data(141)
nedMDnl                            = i_edit_data(142)
intedMDnl                          = i_edit_data(143)
n_editMDnl                         = i_edit_data(144)
iedMDgx                            = i_edit_data(149)
nedMDgx                            = i_edit_data(150)
intedMDgx                          = i_edit_data(151)
n_editMDgx                         = i_edit_data(152)
iedMDgy                            = i_edit_data(153)
nedMDgy                            = i_edit_data(154)
intedMDgy                          = i_edit_data(155)
n_editMDgy                         = i_edit_data(156)
iedMDBVw                           = i_edit_data(157)
nedMDBVw                           = i_edit_data(158)
intedMDBVw                         = i_edit_data(159)
n_editMDBVw                        = i_edit_data(160)
iedMDyl                            = i_edit_data(161)
nedMDyl                            = i_edit_data(162)
intedMDyl                          = i_edit_data(163)
n_editMDyl                         = i_edit_data(164)
i_editMD                           = i_edit_data(165)
n_editMD                           = i_edit_data(166)
i_HDFedit                          = i_edit_data(167)
nd_HDFedit                         = i_edit_data(168)
n_HDFedit                          = i_edit_data(169)
it_edit                            = i_edit_data(170)
ied_global_n                       = i_edit_data(171)
nned_global                        = i_edit_data(172)
inted_global                       = i_edit_data(173)
n_edit_global                      = i_edit_data(174)
ied_global_t                       = i_edit_data(175)
nted_global                        = i_edit_data(176)

DO i           = 1,20
  intedc(i)                        = i_edit_data(200+i)
  nedc(i)                          = i_edit_data(220+i)
  idxedc(i)                        = i_edit_data(240+i)
  intede(i)                        = i_edit_data(260+i)
  nede(i)                          = i_edit_data(280+i)
  idxede(i)                        = i_edit_data(300+i)
  intdmi(i)                        = i_edit_data(320+i)
  nedmi(i)                         = i_edit_data(340+i)
  idxemi(i)                        = i_edit_data(360+i)
  intdma(i)                        = i_edit_data(380+i)
  nedma(i)                         = i_edit_data(400+i)
  idxema(i)                        = i_edit_data(420+i)
  intedh(i)                        = i_edit_data(440+i)
  nedh(i)                          = i_edit_data(460+i)
  idxedh(i)                        = i_edit_data(480+i)
  intdps(i)                        = i_edit_data(500+i)
  nedps(i)                         = i_edit_data(520+i)
  idxeps(i)                        = i_edit_data(540+i)
  intedu(i)                        = i_edit_data(560+i)
  nedu(i)                          = i_edit_data(580+i)
  idxedu(i)                        = i_edit_data(600+i)
  intedy(i)                        = i_edit_data(620+i)
  nedy(i)                          = i_edit_data(640+i)
  idxedy(i)                        = i_edit_data(660+i)
  intdsc(i)                        = i_edit_data(680+i)
  nedsc(i)                         = i_edit_data(700+i)
  idxesc(i)                        = i_edit_data(720+i)
END DO

DO i           = 1,60
  niedn(i)                         = i_edit_data(800+i)
  neden(i)                         = i_edit_data(900+i)
END DO

DO i           = 1,100
  ncyrst(i)                        = i_edit_data(1000+i)
END DO

DO n           = 1,nnu
  intedn(n)                        = i_edit_data(1100+0*nnu+n)
  nedn(n)                          = i_edit_data(1100+1*nnu+n)
  idxedn(n)                        = i_edit_data(1100+2*nnu+n)
END DO

DO n = 1,nnu
  DO k = 1,40
    intdng(k,n)                    = i_edit_data(1200+0*40+(n-1)*3*40+k)
    nedng(k,n)                     = i_edit_data(1200+1*40+(n-1)*3*40+k)
    idxeng(k,n)                    = i_edit_data(1200+2*40+(n-1)*3*40+k)
  END DO
END DO

r_lum                              = d_edit_data(1)
d_lum                              = d_edit_data(2)
r_e_rms                            = d_edit_data(3)
d_e_rms                            = d_edit_data(4)
dtvarplot                          = d_edit_data(5)
dtcomplot                          = d_edit_data(6)
dtimeplot                          = d_edit_data(7)
rinnerb                            = d_edit_data(8)
routerb                            = d_edit_data(9)
r_lumerms                          = d_edit_data(10)
dtnuradplot                        = d_edit_data(11)
r_nurad                            = d_edit_data(12)
rho_nurad                          = d_edit_data(13)
dtnudata                           = d_edit_data(14)
r_nudata                           = d_edit_data(15)
t_nudata                           = d_edit_data(16)
msslag                             = d_edit_data(17)
dmlag                              = d_edit_data(18)
r_pinch                            = d_edit_data(19)
d_pinch                            = d_edit_data(20)
dt_eplot                           = d_edit_data(21)
dt_MDedit1                         = d_edit_data(22)
dt_MDedit2                         = d_edit_data(23)
dt_HDFedit1                        = d_edit_data(24)
dt_HDFedit2                        = d_edit_data(25)
dt_edit                            = d_edit_data(26)
dt_global_ed1                      = d_edit_data(27)
dt_global_ed2                      = d_edit_data(28)

DO i = 1,10
  rhoprint(i)                      = d_edit_data(40+i)
END DO

!-----------------------------------------------------------------------
!
!                  \\\\\ LOAD RADHYD MODULE /////
!
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    DO i = 1,20
      nedc_c (i,ij_ray,ik_ray)     = i_edit_data(220+i)
      nede_c (i,ij_ray,ik_ray)     = i_edit_data(280+i)
      nedmi_c(i,ij_ray,ik_ray)     = i_edit_data(340+i)
      nedma_c(i,ij_ray,ik_ray)     = i_edit_data(400+i)
      nedh_c (i,ij_ray,ik_ray)     = i_edit_data(460+i)
      nedps_c(i,ij_ray,ik_ray)     = i_edit_data(520+i)
      nedu_c (i,ij_ray,ik_ray)     = i_edit_data(580+i)
      nedy_c (i,ij_ray,ik_ray)     = i_edit_data(640+i)
      nedsc_c(i,ij_ray,ik_ray)     = i_edit_data(700+i)
    END DO ! i
  END DO ! ij_ray
END DO ! ik_ray

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    DO n = 1,nnu
      nedn_c(n,ij_ray,ik_ray)      = i_edit_data(1100+1*nnu+n)
    END DO ! n
  END DO ! ij_ray
END DO ! ik_ray

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    DO n = 1,nnu
      DO k = 1,40
        nedng_c(k,n,ij_ray,ik_ray) = i_edit_data(1200+1*40+(n-1)*3*40+k)
      END DO ! k
    END DO ! n
  END DO ! ij_ray
END DO ! ik_ray


RETURN
END SUBROUTINE upack_edit_keys
