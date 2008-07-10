SUBROUTINE eos_write_adios( ndump )
!-----------------------------------------------------------------------
!
!    File:         eos_write
!    Module:       eos_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/00
!
!    Purpose:
!      To dump the equation of state parameters.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!   ndump       : unit number to dump the data
!
!    Output arguments:
!        none
!
!    Include files:
!      cycle_module, eos_snc_x_module, prb_cntl_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, rhoes, eos_i, eosrho
USE prb_cntl_module, ONLY : tnse, tdnse
USE radial_ray_module, ONLY : nprint

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER*8, INTENT(in)              :: ndump           ! unit number to write restart file

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! do index
  
    1 FORMAT ('!                  \\\\\                ///// ')
    2 FORMAT ('!                           EOS KEYS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!........Equation of state parameters...................................
!.......................................................................

!........Equation of state grid

   21 FORMAT ('!      dgrid(i) : divisions per unit change in log(rho).')
   22 FORMAT ('!      tgrid(i) : divisions per unit change in log(t).')
   23 FORMAT ('!      ygrid(i) : divisions per 0.5 change in ye.')
   24 FORMAT ('!      rhoes(k) : permits different eos gridding in different density regimes.')
   25 FORMAT ('esgrid',14x,i10,5x,1pe15.8,22x,'dgrid')
   26 FORMAT ('esgrid',14x,i10,5x,1pe15.8,22x,'tgrid')
   27 FORMAT ('esgrid',14x,i10,5x,1pe15.8,22x,'ygrid')
   28 FORMAT ('esgrid',14x,i10,5x,1pe15.8,22x,'rhoes')

!........EOS identifier

   31 FORMAT ('!      eos_i : equation of state identifier.')
   32 FORMAT ('eos_i ',14x,9x,a1,42x,'eos_i')

!........Equation of state borders

   41 FORMAT ('!      eosrho : the border density between the LS EOS and the BCK EOS (/fm3).')
   42 FORMAT ('eosrho',29x,1pe15.8,22x,'eosrho')

!........NSE flashing and deflashing controls

   51 FORMAT ('!      tnse : temperature at which material is flashed to nse (K).')
   52 FORMAT ('!      tdnse : temperature at which material is deflashed from nse (K).')
   53 FORMAT ('tnse  ',29x,1pe15.8,22x,'tnse')
   54 FORMAT ('tdnse ',29x,1pe15.8,22x,'tdnse ')

!........Document the dump

 1001 FORMAT (' ***EOS keys dump written at cycle        ',i10,' on unit',i5,'***')

!........Equation of state parameters...................................
!.......................................................................

!........Header

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!........Equation of state grid

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,22)
WRITE (ndump,23)
WRITE (ndump,24)
WRITE (ndump,15)
WRITE (ndump,25) (i,dgrid(i),i = 1,3)
WRITE (ndump,26) (i,tgrid(i),i = 1,3)
WRITE (ndump,27) (i,ygrid(i),i = 1,3)
WRITE (ndump,28) (i,rhoes(i),i = 1,2)

!........EOS identifier

WRITE (ndump,13)
WRITE (ndump,31)
WRITE (ndump,15)
WRITE (ndump,32) eos_i

!........Equation of state borders

WRITE (ndump,13)
WRITE (ndump,41)
WRITE (ndump,15)
WRITE (ndump,42) eosrho

!........Equation of state borders

WRITE (ndump,13)
WRITE (ndump,51)
WRITE (ndump,52)
WRITE (ndump,15)
WRITE (ndump,53) tnse
WRITE (ndump,54) tdnse

!........Record the dump................................................
!.......................................................................

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE eos_write_adios
