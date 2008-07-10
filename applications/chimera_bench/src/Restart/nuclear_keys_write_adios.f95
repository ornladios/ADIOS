SUBROUTINE nuclear_keys_write_adios( ndump )
!-----------------------------------------------------------------------
!
!    File:         nuclear_keys_write
!    Module:       nuclear_keys_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/31/05
!
!    Purpose:
!      To dump the nuclear keys and the present values of the nuclear
!       counters.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ndump       : unit number to dump the data
!
!    Output arguments:
!        none
!
!    Include files:
!  cycle_module, nucbrn_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle
USE nucbrn_module, ONLY : inuc, itnuc, ttolnuc, ytolnuc, ynmin, t_cntl_burn, &
& nuc_number, a_name, a_nuc, z_nuc, m_ex_nuc
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
INTEGER                          :: n               ! nuclear data index
INTEGER                          :: n_nucp1         ! nuc_number + 1
INTEGER, PARAMETER               :: i_ray = 1
  
!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT ('!                  \\\\\                ///// ')
    2 FORMAT ('!                         NUCLEAR KEYS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!........Nuclear arrays.................................................
!.......................................................................

!........Nuclear network control parameters

   21 FORMAT ('!      inuc : nuclear reaction network switch.')
   22 FORMAT ('!      itnuc : maximum number of iterations to attempt.')
   23 FORMAT ('!      ttolnuc : temperature convergence parameter for nuclear reaction rate equations.')
   24 FORMAT ('!      ytolnuc : ion abundance convergence parameter for nuclear reaction rate equations.')
   25 FORMAT ('!      ynmin : used with ytolnuc.')
   26 FORMAT ('inuc  ',4x,2i10,42x,'inuc')
   27 FORMAT ('tolnuc',29x,es15.8,22x,'ttolnuc')
   28 FORMAT ('tolnuc',29x,es15.8,22x,'ytolnuc')
   29 FORMAT ('tolnuc',29x,es15.8,22x,'ynmin')

!........Nuclear time step control parameters

   31 FORMAT ('!      t_cntl_burn(1) : nuclear burn temperature change time step criterion.')
   32 FORMAT ('!      t_cntl_burn(2) : nuclear burn composition change time step criterion.')
   33 FORMAT ('tcntrl',14x,i10,5x,es15.8,22x,'t_cntl_burn')

!........Composition data

   41 FORMAT ('!      a_nuc(n) : mass number of the nth nuclear species.')
   42 FORMAT ('!      z_nuc(n) : charge number of the nth nuclear species.')
   43 FORMAT ('!      m_ex_nuc(n) : binding energy of the nth nuclear species (MeV).')
   45 FORMAT ('a_nuc ',a5,i9,3(5x,es15.8))

!........Document the dump

 1001 FORMAT (' ***Nuclear keys dump written at cycle    ',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Initialize.....................................................

n_nucp1          = nuc_number + 1

!-----------------------------------------------------------------------
!
!                   \\\\\ NUCLEAR ARRAYS /////
!
!-----------------------------------------------------------------------

!........Header

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!........Nuclear network control parameters

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,22)
WRITE (ndump,23)
WRITE (ndump,24)
WRITE (ndump,25)
WRITE (ndump,15)
WRITE (ndump,26) inuc,itnuc
WRITE (ndump,27) ttolnuc
WRITE (ndump,28) ytolnuc
WRITE (ndump,29) ynmin

!........Nuclear time step control parameters

WRITE (ndump,13)
WRITE (ndump,31)
WRITE (ndump,32)
WRITE (ndump,15)
WRITE (ndump,33) (i,t_cntl_burn(i),i = 1,2)

!........Composition data

WRITE (ndump,13)
WRITE (ndump,41)
WRITE (ndump,42)
WRITE (ndump,43)
WRITE (ndump,15)

DO n = 1,nuc_number
  WRITE (ndump,45) a_name(n),n,a_nuc(n),z_nuc(n),m_ex_nuc(n)
END DO

!-----------------------------------------------------------------------
!
!                    \\\\\ RECORD THE DUMP /////
!
!-----------------------------------------------------------------------

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE nuclear_keys_write_adios
