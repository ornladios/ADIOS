      module actual_types
!
      logical*4   ::  logical4
      integer*4   ::  integer4
      integer*8   ::  integer8
      real*4      ::  real4
      real*8      ::  real8
      real*16     ::  real16
      complex*8   ::  complex4
      complex*16  ::  complex8
      complex*32  ::  complex16
!
      end module actual_types
!
      module real_prec
!
      use actual_types
!
      integer, parameter :: rl = kind(real8)
      integer, parameter :: rl4= kind(real4)
!
      end module real_prec
!
      module param
!
      use real_prec
!
      integer :: in
      integer :: jn
      integer :: kn
      integer :: ijkn
      integer, parameter :: neqm = 1
!
      integer, parameter :: nbvar = 14
!
      real(rl), parameter :: pi = 3.14159265358979324
      real(rl), parameter :: tiny = 1.0d-99
      real(rl), parameter :: huge = 1.0d+99
!
      real(rl), parameter :: zro  = 0.0d0
      real(rl), parameter :: one  = 1.0d0
      real(rl), parameter :: two  = 2.0d0
      real(rl), parameter :: haf  = 0.5d0
!
      integer, parameter :: nbuff = 40
      integer, parameter :: mreq=300
!
      end module param
!
      module config
!
      use real_prec
      use param
!
!======================================================================
!     ZEUS-MP CONFIGURATION PARAMETERS -- Replacing the CPP macros
!     contained in zeusmp def
!======================================================================
!
!----------------------------------------------------------------------
!     geometry control
!----------------------------------------------------------------------
!
      integer  ::  lgeom, ldimen
      real(rl) :: r_min
      logical  :: xmap_to_twod
!
!----------------------------------------------------------------------
!     dimensions of local process data blocks
!----------------------------------------------------------------------
!
      integer :: izones, jzones, kzones, maxijk
!
!----------------------------------------------------------------------
!     physics control
!----------------------------------------------------------------------
!
      integer :: lrad, leos, nspec, ngrp, nflv
      logical :: xhydro , xgrav, xgrvfft, xptmass, xmhd  , xsphgrv,&
                 xtotnrg, xiso,  xvgrid , xsubav , xforce, xrot
!
      logical, dimension(:), allocatable, save :: x1dflag

      logical :: rhostop
!
!----------------------------------------------------------------------
!     output control
!----------------------------------------------------------------------
!
      logical  :: xascii, xhdf4, xhdf5, xrestart, xtsl
      integer  :: iomode
      real(rl) :: t_out(nbuff-8)
!
!----------------------------------------------------------------------
!     precision
!----------------------------------------------------------------------
!
      real(rl) :: small_no, large_no
!
!----------------------------------------------------------------------
!     special O/S switches
!----------------------------------------------------------------------
!
      logical :: xdec, xt3e, xibm, xunicos

      end module config
      module root

      use real_prec

      real(rl) :: b1floor  , b2floor  , b3floor  ,ciso,                &
                  courno   , dfloor   ,                                &
                  dtal     , dtcs     , dtv1     ,dtv2    , dtv3,      &
                  dtqq     , dtnew    , dtprev   ,                     &
                  avisc_dt ,                                           &
                  dtrd     ,                                           &
                  dt       , dtdump   ,                                &
                  dthdf    , dthist   , dtmin    ,dttsl,               &
                  dtqqi2   , dtnri2   , dtrdi2   ,dtimrdi2,            &
                  dtusr    ,                                           &
                  efloor   , erfloor  , gamma    ,gamm1,               &
                  qcon     , qlin     ,                                &
                  tdump    ,                                           &
                  thdf     , thist    , time     ,tlim      ,cpulim,   &
                  trem     , tsave    , ttsl,                          &
                  tused    , tusr     ,                                &
                  v1floor  , v2floor  , v3floor  ,                     &
                  emf1floor, emf2floor, emf3floor,                     & 
                  gpfloor  , estart   , rad_loss

      integer :: ifsen(6)

      integer ::idebug ,&
       iordb1 , iordb2 , iordb3, iordd ,&
       iorde  , iorder , iords1, iords2, iords3 ,&
       istpb1 , istpb2 , istpb3, istpd , istpe     ,istper,&
       istps1 , istps2 , istps3,&
       ix1x2x3, jx1x2x3,&
       nhy    , nlim   , nred  , mbatch,&
       nwarn  , nseq   , flstat,&
       ioinp  , iotsl  , iolog , iohst , iomov     ,iores ,&
       ioshl

      integer :: nedc_in(20), nedmi_in(20), nedma_in(20),&
                 nedh_in(20), nedps_in(20), nedu_in(20),&
                 nedy_in(20), nedsc_in(20), nedn_in( 6),&
                 nedng_in(100,6)

      character*2  :: id
      character*15 :: hdffile, hstfile, resfile, usrfile
      character*8  :: tslfile

      end module root
      module field

      use real_prec

      real(rl), dimension(:,:,:), allocatable, save :: &
                d, e, p, tt, ye, v1, v2, v3, b1, b2, b3, gp, ggp, &
                oldgp, q11, q22, q33, a_nu, s

      real(rl), dimension(:,:,:,:,:), allocatable, save:: er

      real(rl), dimension(:,:,:,:), allocatable, save ::   abun
 
      real(rl), dimension(:), allocatable, save :: intm, intm_half

      real(rl) :: oldtime

      end module field

      module grid

      use real_prec

      integer :: is, ie, js, je, ks, ke, iga, jga, kga, igcon
      integer :: nx1z, nx2z, nx3z
 
      real(rl), dimension(:), allocatable, save ::&
             x1a  , x2a   ,  x3a   ,&
             x1ai , x2ai  ,  x3ai  ,&
            dx1a  , dx2a  , dx3a   ,&
            dx1ai , dx2ai , dx3ai  ,&
            vol1a , vol2a , vol3a  ,&
            dvl1a , dvl2a , dvl3a  ,&
            dvl1ai, dvl2ai, dvl3ai ,&
             g2a  , g31a  , dg2ad1 ,&
             g2ai , g31ai , dg31ad1,&
             g32a , g32ai , dg32ad2,&
             g4a 
 
      real(rl), dimension(:), allocatable, save ::&
             x1b  ,  x2b  ,  x3b   ,&
             x1bi ,  x2bi ,  x3bi  ,&
            dx1b  , dx2b  , dx3b   ,&
            dx1bi , dx2bi , dx3bi  ,&
            vol1b , vol2b , vol3b  ,&
            dvl1b , dvl2b , dvl3b  ,&
            dvl1bi, dvl2bi, dvl3bi ,&
             g2b  , g31b  , dg2bd1 ,&
             g2bi , g31bi , dg31bd1,&
             g32b , g32bi , dg32bd2,&
             g4b 
 
      real(rl), dimension(:), allocatable, save :: vg1, vg2, vg3

      real(rl) :: x1fac, x2fac, x3fac
!
!     dimension -- "in"
!
      real(rl), dimension(:), allocatable, save ::&
                     x1ah   , x1an   ,&
                     dx1ah  , dx1an  ,&
                     dvl1ah , dvl1an ,&
                     g2ah   , g2an   ,&
                     g31ah  , g31an  ,&
                     x1ahi  , x1ani  ,&
                     dx1ahi , dx1ani ,&
                     dvl1ahi, dvl1ani,&
                     g2ahi  , g2ani  ,&
                     g31ahi , g31ani ,&
                     x1bh   , x1bn   ,&
                     dx1bh  , dx1bn  ,&
                     dvl1bh , dvl1bn ,&
                     g2bh   , g2bn   ,&
                     g31bh  , g31bn  ,&
                     x1bhi  , x1bni  ,&
                     dx1bhi , dx1bni ,&
                     dvl1bhi, dvl1bni,&
                     g2bhi  , g2bni  ,&
                     g31bhi , g31bni

      real(rl), dimension(:), allocatable, save ::&
!
!     dimension -- "jn"
!
                     x2ah   , x2an   ,&
                     dx2ah  , dx2an  ,&
                     dvl2ah , dvl2an ,&
                     g32ah  , g32an  ,&
                     g4ah  , g4an  ,&
                     x2ahi  , x2ani  ,&
                     dx2ahi , dx2ani ,&
                     dvl2ahi, dvl2ani,&
                     g32ahi , g32ani ,&
                     x2bh   , x2bn   ,&
                     dx2bh  , dx2bn  ,&
                     dvl2bh , dvl2bn ,&
                     g32bh  , g32bn  ,&
                     g4bh  , g4bn  ,&
                     x2bhi  , x2bni  ,&
                     dx2bhi , dx2bni ,&
                     dvl2bhi, dvl2bni,&
                     g32bhi , g32bni 

      real(rl), dimension(:), allocatable, save ::&
!
!     dimension -- "kn"
!
                     x3ah   , x3an   ,&
                     dx3ah  , dx3an  ,&
                     dvl3ah , dvl3an ,&
                     x3ahi  , x3ani  ,&
                     dx3ahi , dx3ani ,&
                     dvl3ahi, dvl3ani,&
                     x3bh   , x3bn   ,&
                     dx3bh  , dx3bn  ,&
                     dvl3bh , dvl3bn ,&
                     x3bhi  , x3bni  ,&
                     dx3bhi , dx3bni ,&
                     dvl3bhi, dvl3bni

      end module grid
      module bndry

      use real_prec
      use param

      real(rl) ::  fiis(nbvar), fois(nbvar),&
                   fijs(nbvar), fojs(nbvar),&
                   fiks(nbvar), foks(nbvar)

      integer :: niis( 3),  nois( 3),&
                 nijs( 3),  nojs( 3),&
                 niks( 3),  noks( 3)

      integer :: bvstat(8,nbvar)

      integer, dimension(:,:), allocatable, save ::&
               niib, niib2, niib3, niib23,&
               noib, noib2, noib3, noib23,&
               nijb, nijb3, nijb1, nijb31,&
               nojb, nojb3, nojb1, nojb31,&
               nikb, nikb1, nikb2, nikb12,&
               nokb, nokb1, nokb2, nokb12
 
      integer, dimension(:,:), allocatable, save ::&
               liib,  loib,&
               lijb,  lojb,&
               likb,  lokb

!
!     dimension (:,:,3)
!
      real(rl), dimension(:,:,:), allocatable, save ::&
              diib,  doib,&
              dijb,  dojb,&
              dikb,  dokb

!
!     dimension (:,:,2)
!
      real(rl), dimension(:,:,:), allocatable, save ::&
              eiib,  eoib,&
              eijb,  eojb,&
              eikb,  eokb,&
             v1iib, v1oib,&
             v1ijb, v1ojb,&
             v1ikb, v1okb,&
             v2iib, v2oib,&
             v2ijb, v2ojb,&
             v2ikb, v2okb,&
             v3iib, v3oib,&
             v3ijb, v3ojb,&
             v3ikb, v3okb

!
!     dimension (:,;,:,2)
!
      real(rl), dimension(:,:,:,:), allocatable, save::&
             abiib, aboib, abijb, abojb, abikb, abokb

!
!     dimension (:,:,2)
!
      real(rl), dimension(:,:,:), allocatable, save ::&
             b1iib, b1oib,&
             b1ijb, b1ojb,&
             b1ikb, b1okb,&
             b2iib, b2oib,&
             b2ijb, b2ojb,&
             b2ikb, b2okb,&
             b3iib, b3oib,&
             b3ijb, b3ojb,&
             b3ikb, b3okb

!
!     dimension (:,:,3)
!
      real(rl), dimension(:,:,:), allocatable, save ::&
             emf1iib, emf1oib, emf1ijb, emf1ojb, emf1ikb, emf1okb,&
             emf2iib, emf2oib, emf2ijb, emf2ojb, emf2ikb, emf2okb,&
             emf3iib, emf3oib, emf3ijb, emf3ojb, emf3ikb, emf3okb
 
!
!     dimension (:,:,2)
!
      real(rl), dimension(:,:,:), allocatable, save ::&
             gpiib, gpoib,&
             gpijb, gpojb,&
             gpikb, gpokb
!
!     dimension (:,:,:,2)
!
      real(rl), dimension(:,:,:), allocatable, save ::&
             eriib, eroib, erijb, erojb, erikb, erokb

      end module bndry
      module gravmod

      use real_prec

      real(rl) :: tgrav, ptmass, x1ptm, x2ptm, x3ptm, gp_ext,&
                  m_tot

      integer :: gsup, ifirst

      logical :: xwedge

      end module gravmod

      module cons

      use real_prec

!   physical constants (cgs)
      real(rl), parameter :: clight  = 2.99792d10
      real(rl), parameter :: me      = 9.10953d-28
      real(rl), parameter :: mp      = 1.67265d-24
      real(rl), parameter :: mh      = 1.66053d-24
      real(rl), parameter :: boltz   = 1.38066d-16
      real(rl), parameter :: rad_con = 7.56000d-15
      real(rl), parameter :: hplanck = 6.62618d-27
!
!      real(rl), parameter :: guniv   = 6 6720d-8
!
!     G is now user specified!
!
      real(rl) :: guniv
!
!      real(rl), parameter :: mmw     = 1 0
!
!     MMW is now user specified!
!
      real(rl) :: mmw
      real(rl), parameter :: gasc    = 8.625d7
      real(rl), parameter :: sbc     = 1.8044d-5

!   astronomical constants (cgs)
      real(rl), parameter :: msol = 1.989d33
      real(rl), parameter :: lsol = 3.826d33

!   conversion factors
      real(rl), parameter :: cmau  = 1.49597d13
      real(rl), parameter :: cmpc  = 3.084d18
      real(rl), parameter :: cmkpc = 3.084d21
      real(rl), parameter :: cmkm  = 1.0d5
      real(rl), parameter :: everg = 1.60219d-12

      end module cons

!#ifdef MPI_USED
!      module mpiyes
!
!      use param
!
!      include "mpif h"
!      integer :: stat(MPI_STATUS_SIZE,mreq)
!      integer :: req(mreq)
!
!      end module mpiyes
!#endif /* MPI_USED */

      module mpino

      integer :: stat, req

      end module mpino

      module mpipar

      use param
      use real_prec

      logical :: periodic(3)
      logical :: reorder
      integer :: myid, myid_w, nprocs, nprocs_w, coords(3)
      integer :: ierr, nreq, nsub
      integer :: comm3d
      integer :: ntiles(3)
      integer :: n1m, n1p, n2m, n2p, n3m, n3p
      integer :: i_slice,j_slice,k_slice
      integer :: ier_slice,jer_slice,ker_slice
      integer :: ils_slice,jls_slice,kls_slice
      integer :: iab_slice,jab_slice,kab_slice
      integer :: ilsm_slice,jlsm_slice,klsm_slice
      integer :: ibuf_in(nbuff), ibuf_out(nbuff)

      real(rl) :: buf_in(nbuff), buf_out(nbuff)

      end module mpipar

      module clockmod

      use real_prec

      real(rl4) :: tarray(2), cputime0,&
                   wclock0  

      integer   :: iarray

      end module clockmod

      module scratch

      use param
      use real_prec

!
!     dimension = "ijkn"
!
      real(rl), dimension(:), allocatable, save ::&
            w1da, w1db, w1dc,&
            w1dd, w1de, w1df,&
            w1dg, w1dh, w1di,&
            w1dj, w1dk, w1dl,&
            w1dm, w1dn, w1do,&
            w1dp, w1dq, w1dr,&
            w1ds, w1dt, w1du

!
!     dimension = "in,jn,kn"
!
      real(rl), dimension(:,:,:), allocatable, save ::&
            w3da, w3db, w3dc,&
            w3dd, w3de, w3df,&
            w3dg,&
            w3di, w3dj,&
            w3dh

!
!     dimension = "in,jn,kn,nspec"
!
      real(rl), dimension(:,:,:,:), allocatable, save ::    w4da

!
!     dimension = "in,jn,kn,ngrp,nflv"
!
      real(rl), dimension(:,:,:,:,:), allocatable, save :: w5da

      end module scratch

      module lor_scr

      use real_prec
  
      real(rl), dimension(:,:,:), allocatable, save:: srd1, srd2, srd3

      end module lor_scr

      module restart_arrays

      use real_prec

      integer :: ngridr,ngridi,mgridr
      real(rl), dimension(:), allocatable :: rlgrdvr ! ngridr + mgridr
      integer , dimension(:), allocatable :: ntgrdvr ! ngridi
!
      integer :: nfieldr
      real(rl), dimension(:), allocatable :: rlfldvr ! nfieldr
!
      integer :: nbdryr,nbdryi
      real(rl), dimension(:), allocatable :: rlbdryvr ! nbdryr
      integer , dimension(:), allocatable :: ntbdryvr ! nbdryi
!
      integer, parameter :: nrootr = 53,&
                            nrooti = 41,&
                            nrootch= 4*15+8
      integer, parameter :: nmgfldi = 766
      real(rl) :: rlrtvr(nrootr)
      integer  :: ntrtvr(nrooti)
      integer  :: ntmgfld(nmgfldi)
!
      character*68:: chrtvr
!
      integer, parameter :: ngravr = 7, ngravi = 1
      real(rl) :: rlgrvvr(ngravr)
      integer  :: ntgrvvr(ngravi)

      end module restart_arrays

      module domain
!
      use real_prec
!
      real(rl) :: x1min, x1max, x2min, x2max, x3min, x3max
!
      end module domain
!
      module eos_par
!
      use real_prec
!
      real(rl) :: brydns, rho_edge, bdcgs
      real(rl), parameter :: a_in=28.0D0
      real(rl), parameter :: bunucin=8.447
      real(rl), parameter :: rhoscl=6.022045d-16
      real(rl), parameter :: escl=6.2415d5
      real(rl), parameter :: dscl=1.0d39
      real(rl), parameter :: mn_cgs=1.6605655d-24
!
      end module eos_par

      module impsoln

      use real_prec  

      real(rl) :: totlsit, totnrit , nrpert  , lspert, tot_trout,&
                  tot_tsr, tot_tsyn, tot_tcom, lstime

      integer  :: cgerrcrit
      integer  :: nits, maxitr, ipcflag, totcgit

      end module impsoln
!
!#ifdef GA_USED
!      module glb_arrays

!#include "mafdecls fh"
!#include "global fh"

!      integer :: nga_procs, myid_ga, ga_d
!      logical :: OK
!
!      end module glb_arrays
!#endif /* GA_USED */

      module mgfld_rays

      use real_prec

      integer :: is_g, ie_g, in_g, js_g, je_g, jn_g, ks_g, ke_g, kn_g

      real(rl), dimension(:), allocatable :: d_ray, t_ray, y_ray, &
                                             e_ray, v_ray, r_ray,&
                                             a_ray,od_ray,rn_ray

      real(rl), dimension(:,:,:), allocatable :: psi_ray

      end module mgfld_rays

      module eos_data

      use real_prec

      real(rl),dimension(:,:,:),allocatable,save :: gam_ad, xhv,  xp, xn, xal, machno

      end module eos_data

      module echeck

      use real_prec

      integer  :: i_kinetic,i_grav
      real(rl) :: u_ge_tot,u_ie_tot,u_ke_tot,u_ne_tot,radtot

      end module echeck
