SUBROUTINE array_dimensions_write( ndump )
!-----------------------------------------------------------------------
!
!    File:         array_dimensions_write
!    Module:       array_dimensions_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/23/05
!
!    Purpose:
!      To dump the neutrino energy advection keys and control parameters.
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
!      cycle_module, array_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle
USE array_module, ONLY : nx, ny, nz, nez, nnu, nnc, n_proc, n_proc_y, &
& n_proc_z
USE radial_ray_module, ONLY : nprint

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ndump           ! unit number to write restart file
  
    1 FORMAT ('!                  \\\\\                ///// ')
    2 FORMAT ('!                       ARRAY DIMENSIONS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

   21 FORMAT ('!      nx       :  x-array extent.')
   22 FORMAT ('!      ny       :  y-array extent.')
   23 FORMAT ('!      nz       :  z-array extent.')
   24 FORMAT ('!      nez      : neutrino energy array extent.')
   25 FORMAT ('!      nnu      : neutrino flavor array extent.')
   26 FORMAT ('!      nnc      : abundance array extent.')
   27 FORMAT ('!      n_proc   : number of processors assigned to the run.')
   28 FORMAT ('!      n_proc_y : number of processors assigned to the y-zones.')
   29 FORMAT ('!      n_proc_z : number of processors assigned to the z-zones.')
   41 FORMAT ('nx    ',14x,i10,42x,'nx    ')
   42 FORMAT ('ny    ',14x,i10,42x,'ny    ')
   43 FORMAT ('nz    ',14x,i10,42x,'nz    ')
   44 FORMAT ('nez   ',14x,i10,42x,'nez   ')
   45 FORMAT ('nnu   ',14x,i10,42x,'nnu   ')
   46 FORMAT ('nnc   ',14x,i10,42x,'nnc   ')
   47 FORMAT ('n_proc',14x,i10,42x,'n_proc')
   48 FORMAT ('n_proc',14x,i10,42x,'n_proc_y')
   49 FORMAT ('n_proc',14x,i10,42x,'n_proc_z')

!........Document the dump

 1001 FORMAT (' ***Array dimensions dump written at cycle',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!
!                    \\\\\ ARRAY DIMENSIONS /////
!
!-----------------------------------------------------------------------

!........Header

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!........Array dimensions

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,22)
WRITE (ndump,23)
WRITE (ndump,24)
WRITE (ndump,25)
WRITE (ndump,26)
WRITE (ndump,27)
WRITE (ndump,28)
WRITE (ndump,29)
WRITE (ndump,15)
WRITE (ndump,41) nx
WRITE (ndump,42) ny
WRITE (ndump,43) nz
WRITE (ndump,44) nez
WRITE (ndump,45) nnu
WRITE (ndump,46) nnc
WRITE (ndump,47) n_proc
WRITE (ndump,48) n_proc_y
WRITE (ndump,49) n_proc_z

!-----------------------------------------------------------------------
!  Record the dump
!-----------------------------------------------------------------------

WRITE (nprint,1001) ncycle, ndump

RETURN
END SUBROUTINE array_dimensions_write
