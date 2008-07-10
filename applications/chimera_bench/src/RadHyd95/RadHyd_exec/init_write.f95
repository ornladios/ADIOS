SUBROUTINE init_write( ndump )
!-----------------------------------------------------------------------
!
!    File:         init_write
!    Module:       init_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/00
!
!    Purpose:
!      To dump the model header, cycle number, and restart file name.
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
!      cycle_module, nucbrn_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : head, nouttmp, nprint
USE radial_ray_module, ONLY : ncycle

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ndump           ! unit number to write restart file

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: n_restart_key   ! restart file key

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

   11 FORMAT (a128)
   21 FORMAT (10x,i10)
   31 FORMAT ('output',14x,i10,42x,'nouttmp')
   41 FORMAT ()
   43 FORMAT ('!-----------------------------------------------------------------------')
   45 FORMAT ('!')
   47 FORMAT ('! Set nrst to 0 to read initial model')
   49 FORMAT ('! If nrst /= 0, set noutput to 1 to read in rst_tmp1_ files')
   51 FORMAT ('! If nrst /= 0, set noutput to 2 to read in rst_tmp2_ files')
   53 FORMAT ('! If nrst /= 0, set noutput to 3 to read in restart_model files at cycle number nrst')
   55 FORMAT ('! If nrst /= 0, set noutput to 4 to read in restart_final_mod files')
 1001 FORMAT (' ***Restart directory written at cycle    ',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!  Record the header
!-----------------------------------------------------------------------

WRITE (ndump,11) head

!-----------------------------------------------------------------------
!  Record the cycle number
!-----------------------------------------------------------------------

WRITE (ndump,21) ncycle

!-----------------------------------------------------------------------
!  Record the unit number
!-----------------------------------------------------------------------

n_restart_key      = 1
IF ( ndump == 21 ) n_restart_key = 2

WRITE (ndump,31) n_restart_key

!-----------------------------------------------------------------------
!  Write restart instructions
!-----------------------------------------------------------------------

WRITE (ndump,41) 
WRITE (ndump,41) 
WRITE (ndump,41) 
WRITE (ndump,43) 
WRITE (ndump,45) 
WRITE (ndump,47) 
WRITE (ndump,49) 
WRITE (ndump,51) 
WRITE (ndump,53) 
WRITE (ndump,55) 
WRITE (ndump,45) 
WRITE (ndump,43) 

!-----------------------------------------------------------------------
!  Save unit number for future restart
!-----------------------------------------------------------------------

nouttmp            = ndump

!-----------------------------------------------------------------------
!  Record the dump
!-----------------------------------------------------------------------

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE init_write
