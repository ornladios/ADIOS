SUBROUTINE editc( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editc
!    Module:       editc
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/25/96
!
!    Purpose:
!      To edit the model configuration.
!
!    Subprograms called:
!  date_and_time_print : prints date and time
!  roextrrd            : compute minimum and maximum central densities since last edit
!
!    Input arguments:
!  jr_min     : inner radial zone of region for which configuration
!                edit is to be made.
!  jr_max     : outer radial zone of region for which configuration
!                edit is to be made.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!  prnttest   : true  - test to see if printing criteria is satisfied.
!               false - bypass printing criteria test.
!  iprint     : 0    - do not print to print file.
!               ne 0 - print to print file.
!  nprint     : unit number of print file.
!  iplot      : 0    - do not print to plot file.
!               ne 0 - print to plot file.
!  nplot      : unit number of plot file.
!  nedc(i)    : editc counter for data set i.
!  intedc(i)  : number of cycles between edits of data set i.
!  idxedc(i)  : edit jr_min, jr_max, and every idxedc(i) radial zone
!                between them for data set i.
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module, nu_dist_module,
!  parallel_module, shock_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, n_proc_y, ij_ray_dim, ik_ray_dim
USE numerical_module, ONLY : zero, epsilon, frpi
USE physcnst_module, ONLY : ergfoe, msolar, kmev

USE edit_module, ONLY : nprint, prnttest, nedc, intedc, head, idxedc, &
& mass_ns, vel_ns
USE eos_snc_x_module, ONLY : aesv, nse
USE mdl_cnfg_module, ONLY : rhor, rho, t, ye, u, v, w, r, dr, rojmin, &
& rojmax, grvmss, rstmss
USE nu_dist_module, ONLY : fluxnu
USE parallel_module, ONLY : myid
USE shock_module, ONLY : pq_x
USE t_cntrl_module, ONLY : time, dtnph, dtnmh, dtnph_trans, dtnmhn_trans, &
& tcntrl, jrdt, jadt, jzdt, dt, t_bounce, dtime_hydro, dtime_nuadvct, &
& dtime_nutrns

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1)                :: shock_flag
CHARACTER (len=1)                :: nnse_flag

LOGICAL                          :: first = .true.

INTEGER                          :: j             ! radial zone index
INTEGER                          :: jv            ! do index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: n_time        ! nprint, used for date_and_time
INTEGER                          :: j_ray         ! polar index of the radial ray
INTEGER                          :: k_ray         ! azimuthal index of the radial ray

REAL(KIND=double)                :: flux_nu          ! neutrno energy flux_nu
REAL(KIND=double)                :: lum           ! neutrino luminosity
REAL(KIND=double)                :: pqpv          ! ratio of pseudoviscous pressure to real pressure
REAL(KIND=double)                :: rltmss        ! total rest mass
REAL(KIND=double)                :: rstmssj       ! enclosed mass (solar masses)
REAL(KIND=double)                :: slmass        ! total gravitational mass
REAL(KIND=double)                :: tmev          ! temperature (K)
REAL(KIND=double)                :: t_tb          ! time from  bounce

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1x)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (1x,'ij_ray=',i4,', ik_ray=',i4,',',4x,' ray(',i4,',',i4,')')
  101 FORMAT (1x,'Time and time step control'/)
  103 FORMAT (9x,'Elapsed time=',es14.7,10x,' Time from bounce=',es14.7)
  105 FORMAT (9x,'dt      (  n+1/2)=',es12.5,' dt',6x,'(  n-1/2)=',es12.5/)
  111 FORMAT (9x,'dt (x-hydro Courant criterion)                   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  113 FORMAT (9x,'dt (x-hydro density change criterion)            =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  115 FORMAT (9x,'dt (x-hydro temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  117 FORMAT (9x,'dt (y-hydro Courant criterion)                   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  119 FORMAT (9x,'dt (y-hydro density change criterion)            =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  121 FORMAT (9x,'dt (y-hydro temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  123 FORMAT (9x,'dt (z-hydro Courant criterion)                   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  125 FORMAT (9x,'dt (z-hydro density change criterion)            =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  127 FORMAT (9x,'dt (z-hydro temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  131 FORMAT (9x,'dt (convection crossing criterion)               =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  133 FORMAT (9x,'dt (nuclear temperature change criterion)        =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  135 FORMAT (9x,'dt (nuclear xn change criterion)                 =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  139 FORMAT (9x,'dt (max increase criterion)                      =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  141 FORMAT (9x,'dt (all nu temperature change criterion)         =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  151 FORMAT (9x,'dt (all nu electron fraction change criterion)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  161 FORMAT (9x,'dt (e_nu    psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  163 FORMAT (9x,'dt (e_nubar psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  165 FORMAT (9x,'dt (x_nu    psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  167 FORMAT (9x,'dt (x_nubar psi0_change criterion - transport)   =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  171 FORMAT (9x,'dt (e_nu    psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  173 FORMAT (9x,'dt (e_nubar psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  175 FORMAT (9x,'dt (x_nu    psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  177 FORMAT (9x,'dt (x_nubar psi0_change criterion - e_advection) =',es12.5,' jdrt=',i4,' jadt=',i4,' jzdt=',i4)
  207 FORMAT (9x,'dtime_trans(  n+1/2)=',es12.5,' dtime_trans(  n-1/2)=',es12.5)
  231 FORMAT (9x,'dt(minimum hydro time step                        ) =',es12.5)
  233 FORMAT (9x,'dt(minimum neutrino energy advection time step    ) =',es12.5)
  235 FORMAT (9x,'dt(minimum neutrino source and transport time step) =',es12.5)
  239 FORMAT (132('.')/)
  251 FORMAT (51x,'Model configuration')
  253 FORMAT (50x,21('-')/)
  255 FORMAT ('   j     u          v          w           r         dr         q/p     &
  & lum (f/s)    rstmss        rho      T (MeV)        s         ah         ye'/)
  257 FORMAT (1x,i4,6es11.3,a1,es11.3,es14.6,3(es11.3),es10.2,a1,es11.3)
  261 FORMAT ('0The gravitational mass of the star is ',es14.7,' gms, or ', es14.7,' solar masses')
  263 FORMAT (' The     rest      mass of the star is ',es14.7,' gms, or ', es14.7,' solar masses')
  265 FORMAT (' The mass and velocity of the neutron star are m=',es11.3,' v=',es11.3)
  267 FORMAT (1x,'The maximum and minimum values of rho(2) attained since the last editc are',es10.3, &
  & ' g/cm3 and',es10.3,' g/cm3, respectively')
  301 FORMAT (1x,' cycle =',i5,' time=',es14.7)
  303 FORMAT (1x,i4,11(es11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!**********************************************************************!
!                                                                      !
!                  Model configuration                                 !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!        Time and time step control                                    !
!                                                                      !
!----------------------------------------------------------------------!
! elapsed time     : Total elapsed time from problem initiation        !
! dt(n + 1/2)      : Time step of upcoming hydro cycle                 !
! dt(n - 1/2)      : Time step of just completed hydro cycle           !
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!                                                                      !
!        Model configuration data                                      !
!                                                                      !
!----------------------------------------------------------------------!
! u(j)        'u'          x-component of velocity of zone j           !
!                           (cm s^{-1})                                !
! v(j)        'v'          y-component of velocity of zone j           !
!                           (cm s^{-1})                                !
! w(j)        'w'          z-component of velocity of zone j           !
!                           (cm s^{-1})                                !
! r(j)        'r'          Coordinate radius of zone j (cm)            !
! dr(j)       'dr'         Coordinate width of zone j (cm)             !
! pqpv        'q/p'        Ratio of the pseudoviscosity pressure to    !
!                           the matter pressure in zone j              !
! lum         'lum'        Neutrino (all types) luminosity across      !
!                           zone j (foes/sec)                          !
! rstmssj     'rstmss'     Total rest mass enclosed by zone j (solar   !
!                           masses)                                    !
! rho(j)      'rho'        Proper density of zone j (g/cm**3)          !
! tmev        't (MeV)     Temperature of zone j - 1/2 (MeV)           !
! aesv(j,3,ij_ray,ik_ray)  's' Entropy per baryon (k)                  !
! aesv(j,10,ij_ray,ik_ray) 'ah' Mean atomic mass number of heavy       !
!                           nuclei in zone j - 1/2                     !
! ye(j)       'ye'         Electron number per baryon in zone j - 1/2  !
!----------------------------------------------------------------------!

IF ( prnttest ) THEN
  nedc(1)          = nedc(1) + 1
  IF ( nedc(1) < intedc(1) ) RETURN
  nedc(1)          = 0
END IF

!-----------------------------------------------------------------------
!  Time from bounce
!-----------------------------------------------------------------------

n_time             = nprint
IF ( t_bounce == zero ) THEN
  t_tb             = zero
ELSE
  t_tb             = time - t_bounce
END IF ! t_bounce = 0

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_ray              = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
k_ray              = ( myid/n_proc_y ) * ik_ray_dim + ik_ray

!-----------------------------------------------------------------------
!  Print header
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
WRITE (nprint,7) ij_ray, ik_ray, j_ray, k_ray

CALL date_and_time_print( n_time )

!-----------------------------------------------------------------------
!  Print time step controls
!-----------------------------------------------------------------------

WRITE (nprint,101)
WRITE (nprint,103) time,t_tb
WRITE (nprint,105) dtnph,dtnmh
WRITE (nprint,111) dt( 1),jrdt( 1),jadt( 1),jzdt( 1)
WRITE (nprint,113) dt( 2),jrdt( 2),jadt( 2),jzdt( 2)
WRITE (nprint,115) dt( 3),jrdt( 3),jadt( 3),jzdt( 3)
WRITE (nprint,117) dt( 4),jrdt( 4),jadt( 4),jzdt( 4)
WRITE (nprint,119) dt( 5),jrdt( 5),jadt( 5),jzdt( 5)
WRITE (nprint,121) dt( 6),jrdt( 6),jadt( 6),jzdt( 6)
WRITE (nprint,123) dt( 7),jrdt( 7),jadt( 7),jzdt( 7)
WRITE (nprint,125) dt( 8),jrdt( 8),jadt( 8),jzdt( 8)
WRITE (nprint,127) dt( 9),jrdt( 9),jadt( 9),jzdt( 9)
WRITE (nprint,131) dt(31),jrdt(31),jadt(31),jzdt(31)
WRITE (nprint,133) dt(33),jrdt(33),jadt(33),jzdt(33)
WRITE (nprint,135) dt(34),jrdt(34),jadt(34),jzdt(34)
WRITE (nprint,139) dt(10),jrdt(10),jadt(10),jzdt(10)
WRITE (nprint,141) dt(11),jrdt(11),jadt(11),jzdt(11)
WRITE (nprint,151) dt(16),jrdt(16),jadt(16),jzdt(16)
WRITE (nprint,161) dt(21),jrdt(21),jadt(21),jzdt(21)
WRITE (nprint,163) dt(22),jrdt(22),jadt(22),jzdt(22)
WRITE (nprint,165) dt(23),jrdt(23),jadt(23),jzdt(23)
WRITE (nprint,167) dt(24),jrdt(24),jadt(24),jzdt(24)
WRITE (nprint,171) dt(41),jrdt(41),jadt(41),jzdt(41)
WRITE (nprint,173) dt(42),jrdt(42),jadt(42),jzdt(42)
WRITE (nprint,175) dt(43),jrdt(43),jadt(43),jzdt(43)
WRITE (nprint,177) dt(44),jrdt(44),jadt(44),jzdt(44)
WRITE (nprint,207) dtnph_trans,dtnmhn_trans
WRITE (nprint,231) dtime_hydro
WRITE (nprint,233) dtime_nuadvct
WRITE (nprint,235) dtime_nutrns
WRITE (nprint,239)
WRITE (nprint,251)
WRITE (nprint,253)
WRITE (nprint,255)

!-----------------------------------------------------------------------
!  Model configuration quantities
!-----------------------------------------------------------------------

DO n = 1,nnu
  CALL flux( jr_min, jr_max, n )
END DO

DO jv = jr_min,jr_max,idxedc(1)
  j                = jr_max - jv + jr_min

  pqpv             = pq_x(j,ij_ray,ik_ray)/aesv(j,1,ij_ray,ik_ray)
  flux_nu          = SUM(fluxnu(j,:))
  IF ( DABS(flux_nu) < epsilon ) flux_nu = zero
  lum              = flux_nu * ( frpi * r(j) * r(j) * ergfoe )
  rstmssj          = rstmss(j)/msolar
  tmev             = kmev * t(j)
  shock_flag       = ' '
  IF ( pqpv >= 1.d-01 ) shock_flag = '*'
  nnse_flag        = ' '
  IF ( nse(j,ij_ray,ik_ray) == 0 ) nnse_flag = '*'

!-----------------------------------------------------------------------
!  Print model configuration
!-----------------------------------------------------------------------

  WRITE (nprint,257) j,u(j),v(j),w(j),r(j),dr(j),pqpv,shock_flag,lum, &
& rstmssj,rho(j),tmev,aesv(j,3,ij_ray,ik_ray),aesv(j,10,ij_ray,ik_ray),&
& nnse_flag,ye(j)

END DO

slmass             = grvmss(jr_max)/msolar
rltmss             = rstmss(jr_max)/msolar

WRITE (nprint,261) grvmss(jr_max),slmass
WRITE (nprint,263) rstmss(jr_max),rltmss
WRITE (nprint,265) mass_ns, vel_ns
CALL roextrrd( rojmin, rojmax )
WRITE (nprint,267) rojmax,rojmin

RETURN
END SUBROUTINE editc
