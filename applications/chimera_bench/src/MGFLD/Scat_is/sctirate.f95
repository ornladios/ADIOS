SUBROUTINE sctirate( jr_min, jr_max, ij_ray, ik_ray, n, rho, t, ye, nx )
!-----------------------------------------------------------------------
!
!    File:         sctirate
!    Module:       sctirate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/17/00
!
!    Purpose:
!      To interpolate, for all energy zones k, all radial zones j, and all
!       neutrino types n, the logs of the neutrino and antineutrino isoenergetic
!       scattering inverse mean free paths in a local table of nearest entries
!       created for each zone.
!
!      The table consists of the eight nearest-neighbor entries of the log(rho),
!       log(t), and ye grid.
!      Derivatives of the absorption and emission inverse mean free paths are
!       obtained from the interpolation expressions for these rates by direct
!       differentiation of these expressions
!
!      The scattering functions are included in the multi-group
!       diffusion equations, which have the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2 = xw + yw*psi0 + zw*dpsi0/dr
!
!       where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!      Inelastic neutrino scattering contributes to the terms a0w, a1w, b0w,
!       b1w, c0w, c1w as follows:
!
!            a0w  = -K   Int{ w2'dw'[ phi0in(w,w')psi0(w') + phi0out(w,w')( 1 - psi0(w') ) ] }
!            b0w  = -K/3 Int{ w2'dw'[ phi1in(w,w') - phi1out(w,w') ]psi1(w') }
!            c0w  =  K   Int{ w2'dw'[ phi0in(w,w')psi0(w') ] }
!            a1w  = -K   Int{ w2'dw'[ phi1in(w,w') - phi1out(w,w') ]psi1(w') }
!            b1w  = -K   Int{ w2'dw'[ phi0in(w,w')psi0(w') + phi0out(w,w')( 1 - psi0(w') ) ] }
!            c1w  =  K   Int{ w2'dw'[ phi1in(w,w')psi1(w') ] }
!
!      Isoenergetic scattering contributes to the term b1w as follows:
!
!            b1w  = K*w2*[ phi1(w,w) - phi0(w,w) ]
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  jr_min       : inner radial zone number
!  jr_max       : outer radial zone number
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  n            : neutrino flavor index
!  rho          : density (g cm^{-3})
!  t            : temperature (K)
!  ye           : electron fraction
!  nx           : x-array extent
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  iscat        : 0, all neutrino scattering processes omitted
!                 1, neutrino scattering processes not necessarily omitted
!  in           : 0, neutrino-neutron isoenergetic scattering omitted
!                 1, neutrino-neutron isoenergetic scattering included
!  ip           : 0, neutrino-proton isoenergetic scattering omitted
!                 1, neutrino-proton isoenergetic scattering included
!  ihe          : 0, neutrino-helium isoenergetic scattering omitted
!                 1, neutrino-helium isoenergetic scattering included
!  iheavy       : 0, neutrino-heavy nucleus isoenergetic scattering omitted
!                 1, neutrino-heavy nucleus isoenergetic scattering included
!  nnugp(n)     : number of energy zones for neutrinos of type n
!  rho(j)       : density of zone j
!  t(j)         : temperature of zone j
!  ye(j)        : electron fraction of zone j
!  idris(j)     : rho grid index for zone j
!  itris(j)     : t grid index for zone j
!  iyris(j)     : ye grid index for zone j
!
!    Output arguments (common):
!
!  scti(j,k,n)  : isoenergetic neutrino scattering inverse mean free path
!                (1/cm) for radial zone j, energy zone k
!  sctid(j,k,n) : d(scti(j,k))/d(rho(j))
!  sctit(j,k,n) : d(scti(j,k))/d(t  (j))
!  sctiy(j,k,n) : d(scti(j,k))/d(ye (j))
!
!    Include files:
!  kind_module, numerical_module
!  nu_energy_grid_module, prb_cntl_module, scat_i_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : iscat, in, ip, ihe, iheavy
USE scat_i_module, ONLY : scti, sctid, sctit, sctiy

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: n             ! neutrino flavor index
INTEGER, INTENT(in)              :: nx            ! radial array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rho    ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: t      ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: ye     ! electron fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: noscati

INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: j             ! radial zone index

REAL(KIND=double)                :: coh           ! isoenergetic scattering function
REAL(KIND=double)                :: cohd          ! d(coh)/d(rho)
REAL(KIND=double)                :: coht          ! d(coh)/d(t)
REAL(KIND=double)                :: cohy          ! d(coh)/d(ye)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  No conservative scattering if nnugpmx = 0
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------        
!  No conservative scattering if
!     iscat    = 0                       or
!     in       = 0                       and
!     ip       = 0                       and
!     ihe      = 0                       and
!     iheavy   = 0                       and
!-----------------------------------------------------------------------

noscati            =   ( iscat  == 0 )    .or.  &
&                    ( ( in     == 0 )    .and. &
&                      ( ip     == 0 )    .and. &
&                      ( ihe    == 0 )    .and. &
&                      ( iheavy == 0 ) )  

!-----------------------------------------------------------------------
!  Isoenergetic neutrino scattering
!-----------------------------------------------------------------------

IF ( noscati ) THEN

  scti(jr_min:jr_max,:,:)  = zero
  sctid(jr_min:jr_max,:,:) = zero
  sctit(jr_min:jr_max,:,:) = zero
  sctiy(jr_min:jr_max,:,:) = zero

ELSE

  DO k = 1,nnugpmx
    DO j = jr_min,jr_max
      CALL sctikrnl( j, ij_ray, ik_ray, k, n, rho(j), t(j), ye(j), coh, &
&      cohd, coht, cohy )
      scti(j,k,n)   = coh
      sctid(j,k,n)  = cohd
      sctit(j,k,n)  = coht
      sctiy(j,k,n)  = cohy
    END DO ! j = jr_min,jr_max
  END DO ! k = 1,nnugpmx

END IF ! noscati

RETURN
END SUBROUTINE sctirate
