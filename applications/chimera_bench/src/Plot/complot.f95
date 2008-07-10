SUBROUTINE complot( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         complot
!    Module:       complot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To creates files for important variable profiles at selected
!       times from bounce.
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
!  kind_module, array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!  nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nnu
USE numerical_module, ONLY : zero, third, half, epsilon, ncoef, frpi
USE physcnst_module, ONLY : kmev, rmu, msolar, ergfoe, ergmev

USE cycle_module
USE edit_module, ONLY : icomplt, ncomplt, ncomdump, dtcomplot, nnrst, &
& noutpmt, head, nprint, data_path
USE eos_snc_x_module
USE mdl_cnfg_module, ONLY : jr_max, u, r, t, rho, ye, rstmss, dmrst
USE nu_dist_module, ONLY : unu, dunu, stwt, psi0, psi1, fluxnu, &
& dunujeadt, dunujnisdt, dunujpadt, dunujbadt, dunujnsdt, dunujnnsdt
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, unubi, nnugpmx
USE t_cntrl_module, ONLY : t_bounce, time, dtnph

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)              :: comfile       ! file to write comparison data to

LOGICAL                          :: first_tb = .true.
LOGICAL                          :: l_complt

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: itime         ! used to determine when to WRITE a file
INTEGER                          :: itimeprev     ! used to determine when to WRITE a file
INTEGER                          :: n_time        ! unit number to WRITE date & time
INTEGER                          :: istat         ! open and close file flag

REAL(KIND=double)                :: t_tb          ! time from bounce
REAL(KIND=double)                :: tmult         ! used to determine when to WRITE a file
REAL(KIND=double)                :: r_save        ! saved value of r
REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double), DIMENSION(nnu) :: rnnu         ! numbeer of neutrinos per unit mass
REAL(KIND=double)                :: yenu          ! number of e-neutrinos per baryon
REAL(KIND=double)                :: yanu          ! number of e-antineutrinos per baryon
REAL(KIND=double)                :: yl            ! lepton fraction
REAL(KIND=double)                :: rstmssjmh     ! encloseed mass to j-1/2
REAL(KIND=double)                :: rstmssj       ! encloseed mass to j

REAL(KIND=double), DIMENSION(nx) :: r2            ! r^{2}
REAL(KIND=double), DIMENSION(nx) :: rjmh          ! zone-centered value of r
REAL(KIND=double), DIMENSION(nx) :: rjmh2         ! rjmh^{2}
REAL(KIND=double), DIMENSION(nx) :: lum_r          ! neutrino luminosity (foes/s)
REAL(KIND=double), DIMENSION(nx) :: enuvrms       ! rms neutrino energy (MeV)
REAL(KIND=double), DIMENSION(nx) :: flxfacinv     ! neutrino inverse flux factor
REAL(KIND=double), DIMENSION(nx) :: pinch         ! neutrino pinching parameter

REAL(KIND=double)                :: nu_number     ! proportional to neutrino number/volume
REAL(KIND=double)                :: nu_energy     ! proportional to neutrino energy/volume
REAL(KIND=double)                :: nu_energy2    ! proportional to neutrino energy^2/volume
REAL(KIND=double)                :: mean_e        ! neutrino mean energy
REAL(KIND=double)                :: mean_e2       ! neutrino mean square energy
REAL(KIND=double)                :: flux          ! neutrino flux
REAL(KIND=double)                :: w2dw          ! unu^2 * dunu
REAL(KIND=double)                :: w3dw          ! unu^3 * dunu
REAL(KIND=double)                :: w4dw          ! unu^4 * dunu
REAL(KIND=double)                :: psi0j         ! psi0 interpolated to zone edge

REAL(KIND=double)                :: dudt_ea       ! energy transfered to matter by emision and absorption (MeV/N)
REAL(KIND=double)                :: dudt_nes      ! energy transfered to matter by NES (MeV/N)
REAL(KIND=double)                :: dudt_pa       ! energy transfered to matter by pair-annihilation (MeV/N)
REAL(KIND=double)                :: dudt_ba       ! energy transfered to matter by bremsstrahlung (MeV/N)
REAL(KIND=double)                :: dudt_nns      ! energy transfered to matter by neutrino-N elastic scattering (MeV/N)
REAL(KIND=double)                :: dudt_nnns     ! energy transfered to matter by neutrino-N inelastic scattering (MeV/N)
REAL(KIND=double)                :: dudt_net      ! energy transfered to matter by all processes (MeV/N)

    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (1x,'time=',1pe15.8,' time from bounce=',1pe15.8)
  101 FORMAT ('   j     u           r        rjmh        rho          t          s        ye&
&         yl          m          me         dm'/)
  103 FORMAT ('   j     p           u        elum       alum       tlum       erms       arms&
&       trms      eflxfac-1  aflxfac-1  tflxfa!-1'/)
  105 FORMAT ('   j  pinch_e     pinch_a   pinch_x'/)
  107 FORMAT (1x/)
  201 FORMAT (1x,i3,11(1pe11.3))
  203 FORMAT (1x,i3,11(1pe11.3))
  205 FORMAT (1x,i3,3(1pe11.3))
  301 FORMAT (1x/)
  303 FORMAT (10x,'Neutrino group energies')
  305 FORMAT (10x,32('-')/)
  307 FORMAT (' unui(',i2,')=',1pe11.4,10x,' unubi(',i2,')=',1pe11.4,10x,' dunui(',i2,')=',1pe11.4)
  309 FORMAT (31x,' unubi(',i2,')=',1pe11.4)
  401 FORMAT (1x/)
  403 FORMAT (10x,'Psi0 data')
  405 FORMAT (10x,9('-')/)
  407 FORMAT (' j=',i4)
  409 FORMAT (10(1pe11.3))
  501 FORMAT (1x/)
  503 FORMAT (10x,'Psi1 data')
  505 FORMAT (10x,9('-')/)
  507 FORMAT (' j=',i4)
  509 FORMAT (10(1pe11.3))
 9001 FORMAT (' Error in closing com file in subroutine complot3')
 9051 FORMAT (' Error in closing rstfile in subroutine complot3')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!        \\\\\ EXAMINE CRITERIA FOR WRITING A COMPLOT FILE /////
!
!-----------------------------------------------------------------------

!........Return if icomplt = 0..........................................

IF ( icomplt == 0 ) RETURN

!........Return if t_bounce = 0.........................................

IF ( t_bounce == zero ) RETURN

!........Initialize.....................................................

t_tb               = time - t_bounce
l_complt           = .false.

!........Create complot file at bounce..................................

IF ( t_tb - dtnph <= zero  .and.  t_tb > zero ) l_complt = .true.

!........Create complot file 1 ms after bounce..........................

IF ( t_tb - dtnph <= 1.d-3  .and.  t_tb > 1.d-3 ) l_complt = .true.

!........Create complot file 10 ms after bounce.........................

IF ( t_tb - dtnph <= 1.d-2  .and.  t_tb > 1.d-2 ) l_complt = .true.

!........Create complot filea every multiple of dtcomplot ms............
!....... after bounce...................................................

tmult              = 1.d+3/dtcomplot
itime              = int( t_tb * tmult )
itimeprev          = int( ( t_tb - dtnph ) * tmult )
IF ( itime /= itimeprev ) l_complt = .true.

!........Return if l_complt = .false...................................

IF ( .not. l_complt ) RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ GENERATE A COMPLOT FILE /////
!
!-----------------------------------------------------------------------

n_time             = ncomplt

!........Profile variables..............................................
!.......................................................................

!........Give file a sequential name....................................

ncomdump           = ncomdump + 1
WRITE (comfile,'(a14,i3.3,a2)') '/Plot_Files/com',ncomdump,'.d'
comfile            = TRIM(data_path)//TRIM(comfile)
OPEN (UNIT=ncomplt,FILE=TRIM(comfile),STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=ncomplt,FILE=TRIM(comfile),STATUS='old', &
& POSITION='append')

!........Write headers..................................................

WRITE (ncomplt,3) head
WRITE (ncomplt,5)
CALL date_and_time_print(n_time)
WRITE (ncomplt,7) time,t_tb
WRITE (ncomplt,101)

!........Compute profile variables......................................

r_save              = r(1)
r(1)                = epsilon

DO j = jr_max,2,-1

!........Mass averaged zone-centered radius

  rjmh(j)          = ( half * ( r(j)**3.d+00 + r(j-1)**3.d+00 ) )**third
  tmev             = kmev * t(j)

  DO n = 1,nnu
    rnnu(n)        = zero

!........Lepton fraction

    IF ( nnugp(n) /= 0 ) THEN
      DO k = 1,nnugp(n)
        rnnu(n)    = rnnu(n) + ( ncoef/rho(j) ) * unu(j,k)**2 * dunu(j,k) * psi0(j,k,n)
      END DO
    END IF ! nnugp(n) ne 0

  END DO

  yenu             = rnnu(1) * rmu
  yanu             = rnnu(2) * rmu
  yl               = ye(j) + yenu - yanu

!........Zone-edged and zone-centered enclosed rest masses

  rstmssjmh        = ( rstmss(j-1) + half * dmrst(j) )/msolar
  rstmssj          = rstmss(j)/msolar

!........Print profile variables........................................

WRITE (ncomplt,201) j,u(j),r(j),rjmh(j),rho(j),tmev,aesv(j,3,ij_ray,ik_ray),ye(j), &
& yl,rstmssjmh,rstmssj,dmrst(j)

END DO

r(1)               = r_save

!........Neutrino luminosities, rms energies, and mean..................
!........inverse flux factors...........................................
!.......................................................................

WRITE (ncomplt,107)
WRITE (ncomplt,103)

!........Space averaged zone-centered radius

DO j = 2,jr_max
  rjmh(j)          = half * ( r(j) + r(j-1) )
  rjmh2(j)         = rjmh(j) * rjmh(j)
  r2(j)            = r(j) * r(j)
END DO

rjmh(jr_max+1)         = half * ( r(jr_max+1) + r(jr_max) )
rjmh2(jr_max+1)        = rjmh(jr_max+1) * rjmh(jr_max+1)

DO j = jr_max,2,-1

  DO n = 1,nnu
    lum_r(n)        = zero
  END DO

!........Neutrino luminosities

  DO n = 1,nnu
    lum_r(n)        = frpi * r(j) * r(j) * fluxnu(j,n) * ergfoe/stwt(n)
  END DO

!........Neutrino rms energies and flux factors

  DO n = 1,nnu

    nu_number      = zero
    nu_energy2     = zero
    nu_energy      = zero
    flux           = zero
    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)

      w2dw         = unu(j,k)**2 * dunu(j,k)
      w4dw         = w2dw * unu(j,k)**2
      nu_number    = nu_number    + w2dw  * psi0(j,k,n)
      nu_energy2   = nu_energy2   + w4dw * psi0(j,k,n)
      psi0j        = ( ( rjmh2(j+1) - r2(j) ) * psi0(j,k,n) + ( r2(j) - rjmh2(j) ) * psi0(j+1,k,n) ) &
&                  / ( rjmh2(j+1) - rjmh2(j) )
      nu_energy    = nu_energy  + unu(j,k)**3 * dunu(j,k) * psi0j
      flux         = flux  + unu(j,k)**3 * dunu(j,k) * psi1(j,k,n)

    END DO

    mean_e2        = nu_energy2 /( nu_number   + epsilon )
    enuvrms(n)     = DSQRT( DABS(mean_e2) + epsilon )
    flxfacinv(n)   = nu_energy/( flux + epsilon )

  END DO

  WRITE (ncomplt,203) j,aesv(j,1,ij_ray,ik_ray),aesv(j,2,ij_ray,ik_ray),lum_r(1),lum_r(2),lum_r(3),enuvrms(1), &
&                     enuvrms(2),enuvrms(3),flxfacinv(1),flxfacinv(2),flxfacinv(3)

END DO

!........Neutrino pinching parameters...................................
!.......................................................................

WRITE (ncomplt,107)
WRITE (ncomplt,105)

DO j = jr_max,2,-1
  DO n = 1,nnu

    nu_number      = zero
    nu_energy      = zero
    nu_energy2     = zero

    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)
      w2dw         = unu(j,k)**2 * dunu(j,k)
      w3dw         = w2dw * unu(j,k)
      w4dw         = w2dw * unu(j,k)**2
      nu_number    = nu_number  + w2dw  * psi0(j,k,n)
      nu_energy    = nu_energy  + w3dw * psi0(j,k,n)
      nu_energy2   = nu_energy2 + w4dw * psi0(j,k,n)
    END DO

    mean_e         = nu_energy /( nu_number   + epsilon )
    mean_e2        = nu_energy2 /( nu_number   + epsilon )
    pinch(n)       = 0.75d0 * mean_e2/mean_e**2

  END DO

  WRITE (ncomplt,205) j,pinch(1),pinch(2),pinch(3)

END DO

!........Neutrino energy deposition rates...............................
!.......................................................................

WRITE (ncomplt,107)
WRITE (ncomplt,111)
  111 format ('   j   dudt_ea   dudt_nes    dudt_pa    dudt_ba   dudt_nns   dudt_nnns  dudt_net'/)
  113 FORMAT (1x,i3,7(1pe11.3))

DO j = jr_max,2,-1

  dudt_ea          = - ( dunujeadt (j,1,ij_ray,ik_ray) + dunujeadt (j,2,ij_ray,ik_ray) ) * rmu/( dmrst(j) * ergmev )
  dudt_nes         = - ( dunujnisdt(j,1,ij_ray,ik_ray) + dunujnisdt(j,2,ij_ray,ik_ray) + dunujnisdt(j,3,ij_ray,ik_ray) ) &
&                  * rmu/( dmrst(j) * ergmev )
  dudt_pa          = - ( dunujpadt(j,1,ij_ray,ik_ray)  + dunujpadt(j,2,ij_ray,ik_ray)  + dunujpadt(j,3,ij_ray,ik_ray)  ) &
&                  * rmu/( dmrst(j) * ergmev )
  dudt_ba          = - ( dunujbadt(j,1,ij_ray,ik_ray)  + dunujbadt(j,2,ij_ray,ik_ray)  + dunujbadt(j,3,ij_ray,ik_ray)  ) &
&                  * rmu/( dmrst(j) * ergmev )
  dudt_nns         = - ( dunujnsdt(j,1,ij_ray,ik_ray)  + dunujnsdt(j,2,ij_ray,ik_ray)  + dunujnsdt(j,3,ij_ray,ik_ray)  ) &
&                  * rmu/( dmrst(j) * ergmev )
  dudt_nnns        = - ( dunujnnsdt(j,1,ij_ray,ik_ray) + dunujnnsdt(j,2,ij_ray,ik_ray) + dunujnnsdt(j,3,ij_ray,ik_ray) ) &
&                  * rmu/( dmrst(j) * ergmev )
  dudt_net         = - dudt_ea - dudt_nes - dudt_pa - dudt_ba - dudt_nns - dudt_nnns

  WRITE (ncomplt,113) j,dudt_ea,dudt_nes,dudt_pa,dudt_ba,dudt_nns,dudt_nnns,dudt_net
  
END DO

!........Print neutrino group energies..................................

WRITE (ncomplt,301)
WRITE (ncomplt,303)
WRITE (ncomplt,305)

DO k = 1,nnugpmx 
  WRITE (ncomplt,307) k,unui(k),k,unubi(k),k,dunui(k)
END DO

WRITE (ncomplt,309) nnugpmx+1,unubi(nnugpmx + 1)

!........Print psi0 data................................................

WRITE (ncomplt,401)
WRITE (ncomplt,403)
WRITE (ncomplt,405)

DO j = jr_max,2,-1
  WRITE (ncomplt,407) j
  DO n = 1,nnu
    WRITE (ncomplt,409) (psi0(j,k,n),k=1,10)
    WRITE (ncomplt,409) (psi0(j,k,n),k=11,20)
  END DO
END DO

!........Print psi1 data................................................

WRITE (ncomplt,501)
WRITE (ncomplt,503)
WRITE (ncomplt,505)

DO j = jr_max,2,-1
  WRITE (ncomplt,507) j
  DO n = 1,nnu
    WRITE (ncomplt,509) (psi1(j,k,n),k=1,10)
    WRITE (ncomplt,509) (psi1(j,k,n),k=11,20)
  END DO
END DO

CLOSE (UNIT=ncomplt,STATUS='keep',IOSTAT=istat)
IF ( istat/= 0 ) WRITE (nprint,9001)

!........Done...........................................................

RETURN

END SUBROUTINE complot
