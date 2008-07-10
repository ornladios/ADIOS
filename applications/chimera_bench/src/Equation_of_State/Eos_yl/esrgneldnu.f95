SUBROUTINE esrgneldnu( jr_min, jr_max, nx, rho, t, ye, ij_ray, ik_ray, yl )
!-----------------------------------------------------------------------
!
!    File:         esrgneldnu
!    Module:       esrgneldnu
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/14/02
!
!    Purpose:
!      Given rho(j), t(j), and ye(j), to test each radial zone j
!       between jr_min and jr_max to determine if its thermodynamic
!       state point (rho(j),t(j),yl(j)) is still within
!       the local table of nearest neighbor entries, and derivatives
!       of these entries, created for that zone in a prior cycle.
!       If not, a new local table of nearest neighbor entries, and
!       there derivatives, is created. Subroutine eostlgen is
!       called to fill this table with eos quantities.
!
!      E-neutrino-antineutrino lepton number (leptons/fm**3) and
!       neutrino energy (MeV/fm**3) are stored in i = 10 and i = 11,
!       respectively. The eos is iterated on ye until the yl
!       computed equals the yl inputed.
!
!      The variables that must be passed through common are
!
!  yl(j)          : lepton fraction of radial zone j
!  dgrid(idty(j)) : number of table entries per decade in rho for zone j                                                   !
!  tgrid(idty(j)) : number of table entries per decade in t for zone j
!  ygrid(idty(j)) : number of table entries in ye - yl between ye - yl
!                    = 1.0 and ye - yl = 0 for zone j
!  idty(j)        : index for dgrid, tgrid, and ygrid for zone j
!  rhoes(i)       : density boundary between idty=i and idty=i+1
!  idr(j)         : rho grid index for zone j
!  itr(j)         : t grid index for zone j
!  iyr(j)         : ye - yl grid index for zone j
!  estble(i,j,ida,ita,iya)
!                 : equation of state table array
!
!    Subprograms called:
!        eostlgen, eqstta (module eqstlsa)
!
!    Input arguments:
!  jr_min         : minimum radial zone for which thermodynamic variables
!                    are to be evaluated
!  jr_max         : minimum radial zone for which thermodynamic variables
!                    are to be evaluated
!  nx             : x (radual) array extent
!  rho            : shifted matter density array (g/cm**3).
!  t              : shifted matter temperature array (K).
!  ye             : shifted matter electron fraction array.
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!
!    Output arguments:
!  yl             : shifted lepton fraction array.
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  eos_drv_module, eos_snc_x_module
!
!-----------------------------------------------------------------------	

USE kind_module
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : pi, cm3fm3, rmu, kmev, hbarc, dmnp

USE eos_drv_module, ONLY:  idr, itr, iyr, estble, escnst, estbled, estblet, &
& estbley, escnstd, escnstt, escnsty, tblt, yl_drv=>yl, rho_couple
USE eos_snc_x_module, ONLY : idty, dgrid, tgrid, ygrid, eos_i, rhoes

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min           ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max           ! maximum radial zone index
INTEGER, INTENT(in)              :: nx               ! x (radual) array extent
INTEGER, INTENT(in)              :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray           ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx)  :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx)  :: t   ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)  :: ye  ! shifted matter matter electron fraction array

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx) :: yl  ! shifted matter matter lepton fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (LEN=1)                :: eos_iprev

LOGICAL                          :: first = .true.   ! used for initialization

INTEGER                          :: j                ! radial zone index
INTEGER                          :: i                ! equation of state variable index

INTEGER                          :: id               ! density index
INTEGER                          :: it               ! temperature index
INTEGER                          :: iy               ! lepton fraction index

INTEGER                          :: idd              ! density indices bordering id
INTEGER                          :: itt              ! temperature indices bordering it
INTEGER                          :: iyy              ! lepton fraction indices bordering iy

INTEGER                          :: ida              ! density indices of table
INTEGER                          :: ita              ! temperature indices of table
INTEGER                          :: iya              ! lepton fraction indices of table

INTEGER                          :: idp              ! density cube index
INTEGER                          :: itp              ! temperature cube index
INTEGER                          :: iyp              ! electron fraction cube index

INTEGER                          :: sf               ! eos flag


REAL(KIND=double)                :: pi2              ! pi**2
REAL(KIND=double)                :: kfm              ! ( # nucleons/gram )( cm3/fm3 )
REAL(KIND=double)                :: dbar             ! baryon density ( # /fermi**3 )
REAL(KIND=double)                :: tmev             ! temperature (MeV)

REAL(KIND=double)                :: cmpn             ! neutron chemical potential (MeV)
REAL(KIND=double)                :: cmpp             ! proton chemical potential (MeV)
REAL(KIND=double)                :: cmpe             ! electron chemical potential (MeV)

REAL(KIND=double)                :: dum1,dum2,dum3   ! dummy variables
REAL(KIND=double)                :: etanu            ! ( electron neutrino chemical potential )/kT
REAL(KIND=double)                :: nnuex            ! number of nue's - nuebar's per cubic fermi
REAL(KIND=double)                :: yneu             ! number of nue's - nuebar's per baryon
REAL(KIND=double)                :: ye_yl            ! electron or lepton number deoending on the density

REAL(KIND=double)                :: fd               ! location of log(rho) in the table
REAL(KIND=double)                :: ft               ! location of log(t) in the table
REAL(KIND=double)                :: fy               ! location of ye - yl in the table

REAL(KIND=double), DIMENSION(2,2)      :: estmn      ! minimum values of energy
REAL(KIND=double), DIMENSION(2,2)      :: estmx      ! maximum values of energy
REAL(KIND=double), DIMENSION(4)        :: rhod       ! density at cube corners
REAL(KIND=double), DIMENSION(4)        :: td         ! temperature at cube corners
REAL(KIND=double), DIMENSION(4)        :: yld        ! lepton fraction at cube corners
REAL(KIND=double), DIMENSION(12,2,2,2) :: tbltd      ! density derivatives of EOS quantities
REAL(KIND=double), DIMENSION(12,2,2,2) :: tbltt      ! temperature derivatives of EOS quantities
REAL(KIND=double), DIMENSION(12,2,2,2) :: tblty      ! lepton fraction derivatives of EOS quantities

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set eos_iprev = eos_i for the initial entry into this subroutine
!   and initialize constants.
!-----------------------------------------------------------------------

IF ( first ) THEN
  eos_iprev        = eos_i
  pi2              = pi * pi       ! pi**2
  kfm              = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
  first            = .false.
END IF

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Begin radial loop
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

j_loop: DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Compute equilibrated value of yl(j) if rho > rho_couple
!-----------------------------------------------------------------------

  IF ( rho(j) > rho_couple ) THEN

    CALL eqstta_x( 4, j, ij_ray, ik_ray, rho(j), t(j), ye(j), cmpn, dum1, dum2, dum3 )
    CALL eqstta_x( 5, j, ij_ray, ik_ray, rho(j), t(j), ye(j), cmpp, dum1, dum2, dum3 )
    CALL eqstta_x( 6, j, ij_ray, ik_ray, rho(j), t(j), ye(j), cmpe, dum1, dum2, dum3 )

    dbar            = kfm * rho(j)
    tmev            = kmev * t(j)
    etanu           = ( cmpe - ( cmpn + dmnp - cmpp ) )/tmev
    nnuex           = half * tmev**3/( pi2 * hbarc**3 ) * ( etanu**3 + pi2 * etanu )/3.0d0
    yneu            = nnuex/dbar
    ye_yl           = ye(j) + yneu
    yl_drv(j,ij_ray,ik_ray) = ye_yl
    yl(j)           = ye_yl

  ELSE ! rho <= rho_couple

    ye_yl           = ye(j)
    yl(j)           = ye_yl
    yl_drv(j,ij_ray,ik_ray) = ye_yl

  END IF ! rho > rho_couple

!-----------------------------------------------------------------------
!  Determine which grid to use, then store in the grid index array
!   idty(j,ij_ray,ik_ray).
!-----------------------------------------------------------------------

  IF ( rho(j) < rhoes(1) ) THEN
    idty(j,ij_ray,ik_ray)  = 1
  ELSE IF ( rho(j) < rhoes(2) ) THEN
    idty(j,ij_ray,ik_ray)  = 2
  ELSE
    idty(j,ij_ray,ik_ray)  = 3
  END IF

!-----------------------------------------------------------------------
!  Compute independent variable grid indices
!-----------------------------------------------------------------------

  fd               = dgrid(idty(j,ij_ray,ik_ray)) * dlog10(rho(j) )
  ft               = tgrid(idty(j,ij_ray,ik_ray)) * dlog10(t(j)   )
  fy               = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye_yl )

  id               = IDINT(fd)
  it               = IDINT(ft)
  iy               = IDINT(fy)

!-----------------------------------------------------------------------
!  Test whether the thermodynamic state of radial zone j is still
!   within its local table of nearest neighbor entries, and if the
!   equation of state identifier is the same.
!-----------------------------------------------------------------------

  IF ( id    == idr(j,ij_ray,ik_ray)  .and.   &
&      it    == itr(j,ij_ray,ik_ray)  .and.   &
&      iy    == iyr(j,ij_ray,ik_ray)  .and.   &
&      eos_i == eos_iprev           ) CYCLE

!-----------------------------------------------------------------------
!  Recompute local table if the thermodynamic state of radial zone j is
!   no longer inside its local table, or if the equation of state
!  identifier is different.
!-----------------------------------------------------------------------

  DO idd = id-1,id+2
    DO  itt = it-1,it+2
      DO iyy = iy-1,iy+2

        idp        = ( idd - id + 2 )
        itp        = ( itt - it + 2 )
        iyp        = ( iyy - iy + 2 )
        CALL eostldnugen( j, ij_ray, ik_ray, idd, itt, iyy, idp, itp, iyp, sf, rho(j) )

      END DO ! iyy = iy-1,iy+2
    END DO ! itt = it-1,it+2
  END DO ! idd = id-1,id+2

!-----------------------------------------------------------------------
!  Save the thermodymaic state indices for radial zone j
!-----------------------------------------------------------------------

  idr(j,ij_ray,ik_ray) = id
  itr(j,ij_ray,ik_ray) = it
  iyr(j,ij_ray,ik_ray) = iy

!-----------------------------------------------------------------------
!  Compute independent variables at grid points for use in constructing
!   finite differences
!-----------------------------------------------------------------------

  DO idd = id-1,id+2
    ida            = ( idd - id + 2 )
    rhod(ida)      = 10.d0**( FLOAT(idd)/dgrid(idty(j,ij_ray,ik_ray)) )
  END DO ! idd = id-1,id+2

  DO itt = it-1,it+2
    ita            = ( itt - it + 2 )
    td(ita)        = 10.d0**( FLOAT(itt)/tgrid(idty(j,ij_ray,ik_ray)) )
  END DO ! itt = it-1,it+2

  DO iyy = iy-1,iy+2
    iya            = ( iyy - iy + 2 )
    yld(iya)       = one - ( FLOAT(iyy)/ygrid(idty(j,ij_ray,ik_ray)) )
  END DO ! iyy = iy-1,iy+2

!-----------------------------------------------------------------------
!  Compute derivatives by centered differences
!-----------------------------------------------------------------------

  DO i = 1,3
    DO ida = 2,3
      DO ita = 2,3
        DO iya = 2,3
           tbltd(i,ida-1,ita-1,iya-1) = ( tblt(i,ida+1,ita,iya) - tblt(i,ida-1,ita,iya) ) &
&                  / ( rhod(ida+1) - rhod(ida-1) )
           tbltt(i,ida-1,ita-1,iya-1) = ( tblt(i,ida,ita+1,iya) - tblt(i,ida,ita-1,iya) ) &
&                  / ( td(ita+1) - td(ita-1) )
           tblty(i,ida-1,ita-1,iya-1) = ( tblt(i,ida,ita,iya+1) - tblt(i,ida,ita,iya-1) ) &
&                  / ( yld(iya+1) - yld(iya-1) )
        END DO ! iya = 2,3
      END DO ! ita = 2,3
    END DO ! ida = 2,3
  END DO ! i = 1,3

!-----------------------------------------------------------------------
!  Compute escnst(i,j), which is added to the thermodynamic table
!   entries for the ith quantity of radial zone j to ensure that all 
!   entries are positive so that logarithms can be taken. This quantity
!   will be subtracted from the interpolated value of the ith
!   thermodynmaic quantity to restore its original value. Do the same
!   for escnstd(i,j) escnstt(i,j), and escnsty(i,j) for the derivatives
!   of the ith thermodynamic quantity with respect to rho, t, and ye - yl.
!-----------------------------------------------------------------------

  DO i = 1,3
    escnst (i,j,ij_ray,ik_ray)   = zero
    escnstd(i,j,ij_ray,ik_ray)   = zero
    escnstt(i,j,ij_ray,ik_ray)   = zero
    escnsty(i,j,ij_ray,ik_ray)   = zero
  END DO !  i = 1,3

  DO i = 1,3
    DO ida = 2,3
      DO ita = 2,3
        DO iya = 2,3
          escnst(i,j,ij_ray,ik_ray)  = MIN( escnst (i,j,ij_ray,ik_ray), tblt (i,ida  ,ita  ,iya  ) )
          escnstd(i,j,ij_ray,ik_ray) = MIN( escnstd(i,j,ij_ray,ik_ray), tbltd(i,ida-1,ita-1,iya-1) )
          escnstt(i,j,ij_ray,ik_ray) = MIN( escnstt(i,j,ij_ray,ik_ray), tbltt(i,ida-1,ita-1,iya-1) )
          escnsty(i,j,ij_ray,ik_ray) = MIN( escnsty(i,j,ij_ray,ik_ray), tblty(i,ida-1,ita-1,iya-1) )
        END DO ! iya = 2,3
      END DO ! ita = 2,3
    END DO ! ida = 2,3
  END DO ! i = 1,3

  DO i = 1,3
    escnst (i,j,ij_ray,ik_ray)   = -2.d0 * escnst (i,j,ij_ray,ik_ray)
    escnstd(i,j,ij_ray,ik_ray)   = -2.d0 * escnstd(i,j,ij_ray,ik_ray)
    escnstt(i,j,ij_ray,ik_ray)   = -2.d0 * escnstt(i,j,ij_ray,ik_ray)
    escnsty(i,j,ij_ray,ik_ray)   = -2.d0 * escnsty(i,j,ij_ray,ik_ray)
  END DO

!-----------------------------------------------------------------------
!  Store logarithms of table entries
!-----------------------------------------------------------------------

  DO ida = 1,2
    DO ita = 1,2
      DO iya = 1,2
         DO i = 1,3
           estble (i,j,ida,ita,iya,ij_ray,ik_ray) = LOG10( tblt (i,ida+1,ita+1,iya+1) + escnst (i,j,ij_ray,ik_ray) + epsilon )
           estbled(i,j,ida,ita,iya,ij_ray,ik_ray) = LOG10( tbltd(i,ida  ,ita  ,iya  ) + escnstd(i,j,ij_ray,ik_ray) + epsilon )
           estblet(i,j,ida,ita,iya,ij_ray,ik_ray) = LOG10( tbltt(i,ida  ,ita  ,iya  ) + escnstt(i,j,ij_ray,ik_ray) + epsilon )
           estbley(i,j,ida,ita,iya,ij_ray,ik_ray) = LOG10( tblty(i,ida  ,ita  ,iya  ) + escnsty(i,j,ij_ray,ik_ray) + epsilon )
        END DO
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!  Rearrange for monotonicity of energy
!-----------------------------------------------------------------------

  DO ida = 1,2
    DO iya = 1,2
      estmn(ida,iya) = MIN( estble(2,j,ida,1,iya,ij_ray,ik_ray), estble(2,j,ida,2,iya,ij_ray,ik_ray) )
      estmx(ida,iya) = MAX( estble(2,j,ida,1,iya,ij_ray,ik_ray), estble(2,j,ida,2,iya,ij_ray,ik_ray) )
      estble(2,j,ida,1,iya,ij_ray,ik_ray) = estmn(ida,iya)
      estble(2,j,ida,2,iya,ij_ray,ik_ray) = estmx(ida,iya)
    END DO
  END DO

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  End radial loop
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END DO j_loop

!-----------------------------------------------------------------------
!  Table generation complete. Reset eos_iprev
!-----------------------------------------------------------------------

eos_iprev          = eos_i

RETURN
END SUBROUTINE esrgneldnu

