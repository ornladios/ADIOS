SUBROUTINE e_advct_write( ndump, nnu )
!-----------------------------------------------------------------------
!
!    File:         e_advct_write
!    Module:       e_advct_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/00
!
!    Purpose:
!      To dump the neutrino energy advection keys and control parameters.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!   ndump       : unit number to dump the data
!   nnu         : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!      cycle_module, e_advct_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle
USE e_advct_module, ONLY : ivc_x, ivc_y, ivc_z, tau_advct, t_cntrl_e_advct, &
& psivmin, rhomin_y_eadvect, rhomin_z_eadvect
USE radial_ray_module, ONLY : nprint

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ndump           ! unit number to write restart file
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: n               ! neutrino flavor index

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT ('!                  \\\\\                /////')
    2 FORMAT ('!                       E_ADVECTION KEYS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!-----------------------------------------------------------------------
!
!                    \\\\\ E-ADVECTION ARRAYS /////!                    
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino energy advection controls
!-----------------------------------------------------------------------

   21 FORMAT ('!      ivc_x : neutrino x energy advection switch.')
   22 FORMAT ('!      ivc_y : neutrino y energy advection switch.')
   23 FORMAT ('!      ivc_z : neutrino z energy advection switch.')
   24 FORMAT ('!      tau_advct : optical depth above which neutrinos are advected like a &
&gamma = 4/3 gas.')
   25 FORMAT ('!      rhomin_y_eadvect : density below which y-neutrino energy advection &
&is turned off.')
   26 FORMAT ('!      rhomin_z_eadvect : density below which z-neutrino energy advection &
&is turned off.')
   27 FORMAT ('ivc   ',14x,i10,42x,'ivc_x')
   28 FORMAT ('ivc   ',14x,i10,42x,'ivc_y')
   29 FORMAT ('ivc   ',14x,i10,42x,'ivc_z')
   30 FORMAT ('tau   ',29x,1pe15.8,22x,'tau_advct')
   31 FORMAT ('rhomin',29x,1pe15.8,22x,'rhomin_y_eadvect')
   32 FORMAT ('rhomin',29x,1pe15.8,22x,'rhomin_y_eadvect')

!-----------------------------------------------------------------------
!  E_advection time step control criteria
!-----------------------------------------------------------------------

   41 FORMAT ('!      t_cntrl_e_advct(n): n-neutrino zero-moment change time step criterion &
&due to advection')
   42 FORMAT ('!      psivmin(n): parameters used in determining the psi0 change time step.')
   43 FORMAT ('tcntrl',14x,i10,5x,1pe15.8,22x,'t_cntrl_e_advct')
   44 FORMAT ('psivtl',14x,i10,5x,1pe15.8,22x,'psivmin')

!-----------------------------------------------------------------------
!  Document the dump
!-----------------------------------------------------------------------

 1001 FORMAT (' ***E_advection keys dump written at cycle',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!  Header
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!-----------------------------------------------------------------------
!  Neutrino energy advection controls
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,22)
WRITE (ndump,23)
WRITE (ndump,24)
WRITE (ndump,25)
WRITE (ndump,26)
WRITE (ndump,15)
WRITE (ndump,27) ivc_x
WRITE (ndump,28) ivc_y
WRITE (ndump,29) ivc_z
WRITE (ndump,30) tau_advct
WRITE (ndump,31) rhomin_y_eadvect
WRITE (ndump,32) rhomin_z_eadvect

!-----------------------------------------------------------------------
!  E_advection time step control criteria
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,41)
WRITE (ndump,42)
WRITE (ndump,15)
WRITE (ndump,43) (n,t_cntrl_e_advct(n),n = 1,nnu)
WRITE (ndump,44) (n,psivmin(n),n = 1,nnu)

!-----------------------------------------------------------------------
!  Record the dump
!-----------------------------------------------------------------------

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE e_advct_write
