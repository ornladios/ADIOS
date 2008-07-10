SUBROUTINE varplot( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         varplot
!    Module:       varplot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To creates files for important variable profiles at selected
!       cycle numbers.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!  nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nx, nnu
USE numerical_module, ONLY : zero, half, third, epsilon, frpi, ncoef
USE physcnst_module, ONLY : kmev, rmu, msolar, ergfoe

USE cycle_module
USE edit_module, ONLY : ivarplt, nvar, nvarint, nvarplt, dtvarplot, &
& rhoprint, nvardump, nprint, head, data_path
USE eos_snc_x_module, ONLY : aesv
USE mdl_cnfg_module, ONLY : rhor, rho, r, t, ye, rstmss, dmrst, u, jr_max
USE nu_dist_module, ONLY : unu, dunu, stwt, psi0, fluxnu, psi1
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY : t_bounce, time, dtnph

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray      ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray      ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)              :: varfile       ! file to WRITE com data to

LOGICAL                          :: l_varplt

INTEGER                          :: i             ! radial zone index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: itime         ! used to determine when to WRITE a file
INTEGER                          :: itimeprev     ! used to determine when to WRITE a file
INTEGER                          :: n_time        ! unit number to WRITE date & time
INTEGER                          :: istat         ! open and close file flag

REAL(KIND=double)                :: t_tb          ! time from bounce
REAL(KIND=double)                :: tmult         ! used to determine when to WRITE a file
REAL(KIND=double)                :: rsave         ! stored value of r(1)
REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: yenu          ! e-neutrino fraction
REAL(KIND=double)                :: yanu          ! e-antineutrino fraction
REAL(KIND=double)                :: yl            ! lepton fraction
REAL(KIND=double)                :: rstmssjmh     ! rest mass enclosed by zone center (m_solar)
REAL(KIND=double)                :: rstmssj       ! rest mass enclosed by zone edge (m_solar)

REAL(KIND=double)                :: dennjd        ! proportional to neutrino number
REAL(KIND=double)                :: dene2jd       ! proportional to square of neutrino energy
REAL(KIND=double)                :: eden          ! proportional to neutrino number at zone edge
REAL(KIND=double)                :: flux          ! proportional to neutrino flux at zone edge
REAL(KIND=double)                :: coefn         ! coefficient for computing neutrino number
REAL(KIND=double)                :: coefe2        ! coefficient for computing square of neutrino energy
REAL(KIND=double)                :: enuv2jd       ! neutrino energy sqyared

REAL(KIND=double), DIMENSION(nx)  :: rjmh         ! mass centered radius
REAL(KIND=double), DIMENSION(nnu) :: rnnu         ! neutrino number per gram
REAL(KIND=double), DIMENSION(nnu) :: lum_r         ! neutrino luminosity
REAL(KIND=double), DIMENSION(nnu) :: enuvrms      ! neutrino rms energy (MeV)
REAL(KIND=double), DIMENSION(nnu) :: flxfacinv    ! neutrino inverse flux factor

    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (1x,'time=',1pe15.8,'time from bounce=',1pe15.8)
  101 FORMAT ('   j     u           r        rjmh        rho          t &
&         s        ye         yl          m          me         dm'/)
  103 FORMAT ('   j     p           u        elum       alum       tlum &
&      erms       arms       trms      eflxfac-1  aflxfac-1  tflxfac-1'/)
  105 FORMAT (1x/)
  201 FORMAT (1x,i3,11(1pe11.3))
  203 FORMAT (1x,i3,11(1pe11.3))
 8001 FORMAT (' File nvarplt cannot be opened in subroutime varplot')
 9001 FORMAT (' File nvarplt cannot be closed in subroutime varplot')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if ivarplt = 0
!-----------------------------------------------------------------------

IF ( ivarplt == 0 ) RETURN

!-----------------------------------------------------------------------
!  Initialize for print criteria
!-----------------------------------------------------------------------

l_varplt           = .false.
t_tb               = time - t_bounce

!-----------------------------------------------------------------------
!  Create varplot file at bounce
!-----------------------------------------------------------------------

IF ( t_tb - dtnph <= zero  .and.  t_tb > zero ) l_varplt = .true.

!-----------------------------------------------------------------------
!  Create varplot file 1 ms after bounce
!-----------------------------------------------------------------------

IF ( t_tb - dtnph <= 1.d-3  .and.  t_tb > 1.d-3 ) l_varplt = .true.

!-----------------------------------------------------------------------
!  Create varplot files every multiple of dtvarplot ms after bounce
!-----------------------------------------------------------------------

tmult              = 1.d+3/dtvarplot
itime              = int( t_tb * tmult )
itimeprev          = int( ( t_tb - dtnph ) * tmult )
IF ( itime /= itimeprev ) l_varplt = .true.

!-----------------------------------------------------------------------
!  Create varplot file if nvar = nvarint
!-----------------------------------------------------------------------

nvar               = nvar + 1
IF ( nvar >= nvarint ) THEN
  nvar             = 0
  l_varplt         = .true.
END IF ! nvar >= nvarint

!-----------------------------------------------------------------------
!  Create varplot file if rho = rhoprint(i)
!-----------------------------------------------------------------------

DO i = 1,10
  IF ( rhor(2) < rhoprint(i)  .and.  rho (2) >= rhoprint(i) ) l_varplt = .true.
END DO

!-----------------------------------------------------------------------
!  Return if l_varplt = .false
!-----------------------------------------------------------------------

IF ( .not. l_varplt ) RETURN

!-----------------------------------------------------------------------
!
!                   \\\\\ PRINT VAR FILES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

n_time             = nvarplt
IF ( t_bounce == zero ) THEN
  t_tb             = zero
ELSE
  t_tb             = time - t_bounce
END IF ! t_bounce = 0

!-----------------------------------------------------------------------
!  Generate var file name
!-----------------------------------------------------------------------

nvardump           = nvardump + 1
WRITE (varfile,'(a13,i4.4,a2)') '/Var_Files/var',nvardump,'.d'
varfile            = TRIM(data_path)//TRIM(varfile)

OPEN (UNIT=nvarplt,FILE=varfile,STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nvarplt, FILE=TRIM(varfile), STATUS='old', &
& POSITION='append', IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,8001)
  STOP
END IF

!-----------------------------------------------------------------------
!  Print header
!-----------------------------------------------------------------------

WRITE (nvarplt,3) head
WRITE (nvarplt,5)
WRITE (nvarplt,7) time,t_tb
CALL date_and_time_print(n_time)
WRITE (nvarplt,101)

rsave              = r(1)
r(1)               = epsilon

DO j = jr_max,2,-1

  rjmh(j)          = ( half * ( r(j)**3.d+00 + r(j-1)**3.d+00 ) )**third
  tmev             = kmev * t(j)

  DO n = 1,nnu
    rnnu(n)        = zero
    IF ( nnugp(n) == 0 ) CYCLE
    DO k = 1,nnugp(n)
      rnnu(n)      = rnnu(n) + ( ncoef/rho(j) ) * unu(j,k)**2 * dunu(j,k) * psi0(j,k,n)
    END DO
  END DO

  yenu             = rnnu(1) * rmu
  yanu             = rnnu(2) * rmu
  yl               = ye(j) + yenu - yanu
  rstmssjmh        = ( rstmss(j-1) + half * dmrst(j) )/msolar
  rstmssj          = rstmss(j)/msolar

!-----------------------------------------------------------------------
!  Print model configuration
!-----------------------------------------------------------------------

  WRITE (nvarplt,201) j, u(j), r(j), rjmh(j), rho(j), tmev, aesv(j,3,ij_ray,ik_ray), &
&  ye(j), yl, rstmssjmh, rstmssj, dmrst(j)

END DO

r(1)               = rsave

!-----------------------------------------------------------------------
!  Print header for neutrino data
!-----------------------------------------------------------------------

WRITE (nvarplt,105)
WRITE (nvarplt,103)

!-----------------------------------------------------------------------
!  Compute neutrino luminosities, rms energies, and mean inverse flux
!   factors
!-----------------------------------------------------------------------

DO j = jr_max,2,-1

  DO n = 1,nnu
        lum_r(n)    = zero
  END DO

  DO n = 1,nnu
    lum_r(n)        = frpi * r(j) * r(j) * fluxnu(j,n) * ergfoe/stwt(n)
  END DO

  DO n = 1,nnu

    dennjd         = zero
    dene2jd        = zero
    eden           = zero
    flux           = zero
    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)
      coefn        = unu(j,k)**2 * dunu(j,k)
      coefe2       = coefn * unu(j,k)**2
      dennjd       = dennjd    + coefn  * psi0(j,k,n)
      dene2jd      = dene2jd   + coefe2 * psi0(j,k,n)
      eden         = eden  + unu(j,k)**3 * dunu(j,k) * half * ( psi0(j,k,n) + psi0(j+1,k,n) )
      flux         = flux  + unu(j,k)**3 * dunu(j,k) * psi1(j,k,n)
    END DO

    enuv2jd        = dene2jd /( dennjd   + epsilon )
    enuvrms(n)     = DSQRT( dabs(enuv2jd) + epsilon )
    flxfacinv(n)   = eden/( flux + epsilon )

  END DO

!-----------------------------------------------------------------------
!  Print neutrino data
!-----------------------------------------------------------------------

  WRITE (nvarplt,203) j, aesv(j,1,ij_ray,ik_ray), aesv(j,2,ij_ray,ik_ray), &
&  lum_r(1), lum_r(2), lum_r(3), enuvrms(1), enuvrms(2), enuvrms(3),       &
&  flxfacinv(1), flxfacinv(2), flxfacinv(3)

END DO

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------

CLOSE (UNIT=nvarplt,STATUS='keep',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (nprint,9001)
  STOP
END IF ! istat /= 0

RETURN
END SUBROUTINE varplot
