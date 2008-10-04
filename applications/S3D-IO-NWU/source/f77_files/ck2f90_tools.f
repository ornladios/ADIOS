!========================================================================================
      subroutine set_number_elem_spec_reac(ne,ns,nr,myid,io)
!========================================================================================
!     subroutine sets:
!     number of elements (ne)
!     number of species (ns)
!     number of reactions (nr)
!     
!     based on current values in 'ckstrt.h'
!
!     note that there is actually a chemkin utility to do this (routine CKINDX)
!----------------------------------------------------------------------------------------
      implicit none

      include 'ckstrt.h'

      integer ne,ns,nr,myid,io
!----------------------------------------------------------------------------------------
!     set ne, ns, nr

      ne=nmm
      ns=nkk
      nr=nii

      if (myid==0) then
         write(io,1) ' number of elements in reaction mechansim = ',ne
         write(io,1) ' number of species in reaction mechansim  = ',ns
         write(io,1) ' number of steps in reaction mechansims   = ',nr
         write(io,*)
 1       format(a44,i4)
      endif

      return
      end
!========================================================================================
      subroutine set_number_third_body_reactions(ntb_reac,io)
!========================================================================================
!     subroutine sets:
!     number of elements (ne)
!     number of species (ns)
!     number of reactions (nr)
!     
!     based on current values in 'ckstrt.h'
!----------------------------------------------------------------------------------------
      implicit none

      include 'ckstrt.h'

      integer ntb_reac,io
!----------------------------------------------------------------------------------------
!     set ntb

      ntb_reac=NTHB

      write(io,1) ' number of reaction third-body reactions = ',ntb_reac
 1    format(a43,i4)
!----------------------------------------------------------------------------------------
      return
      end
