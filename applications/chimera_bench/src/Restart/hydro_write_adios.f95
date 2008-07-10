SUBROUTINE hydro_write_adios( ndump )
!-----------------------------------------------------------------------
!
!    File:         hydro_write
!    Module:       hydro_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/25/05
!
!    Purpose:
!      To dump the hydro keys and hydro parameters.
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
!  boundary_module, convect_module, cycle_module, evh1_global,
!  prb_cntl_module, radial_ray_module, shock_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE boundary_module, ONLY : iubcjmn, ubcjmn, iubcjmx, ubcjmx, iuset, uset, r_u, &
& ipbnd, pbound
USE convect_module, ONLY : adrftcv, alphlmix
USE cycle_module, ONLY : ncycle
USE evh1_global, ONLY : degen
USE prb_cntl_module, ONLY : ihydro, irelhy, ilapsehy, ilcnvct, iemix, iyemix, &
& ipsimix, ipcnvct, iscnvct
USE radial_ray_module, ONLY : imin, imax, nprint
USE shock_module, ONLY : ipq, q0_x
USE t_cntrl_module, ONLY : idtj, jtv, jtvshk, ncyrho, ncyshk, rhojtv, tcntrl

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
INTEGER                          :: j               ! radial zone index
  
    1 FORMAT ('!                  \\\\\                ///// ')
    2 FORMAT ('!                          HYDRO KEYS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!........Hydro arrays...................................................
!.......................................................................

!........Hydro switch

   21 FORMAT ('!      ihydro : hydro switch.')
   22 FORMAT ('ihydro',14x,i10,42x,'ihydro')

!........Newtonian - GR switch

   31 FORMAT ('!      irelhy : Newtonian - GR switch.')
   32 FORMAT ('irel  ',14x,i10,42x,'irelhy')
   33 FORMAT ('irel  ',14x,i10,42x,'ilapsehy')

!........Degeneracy - nondegeneracy criterion

   41 FORMAT ('!      degen : Degeneracy - nondegeneracy criterion.')
   42 FORMAT ('degen ',29x,1pe15.8,22x,'degen')

!........Individual time steps

   51 FORMAT ('!      idtj : individual time step flag.')
   52 FORMAT ('!      ncyshk : idtj is set to 1 for ncycle >= ncyshk.')
   53 FORMAT ('!      jtv : iradial zones jmin to jtv have individual hydro time steps if idtj = 1.')
   54 FORMAT ('!      jtvshk : jtv = jshock - jtvshk.')
   55 FORMAT ('!      ncyrho : idtj is set to 1 for ncycle >= ncyshk.')
   56 FORMAT ('!      rhojtv : jtv is the maximum j such that rho(j) > rhojtv.')
   57 FORMAT ('vtime ',14x,i10,42x,'idtj')
   58 FORMAT ('vtime ',14x,i10,42x,'jtv')
   59 FORMAT ('vtime ',14x,i10,42x,'jtvshk')
   60 FORMAT ('vtime ',14x,i10,42x,'ncyshk')
   61 FORMAT ('vtime ',14x,i10,42x,'ncyrho')
   62 FORMAT ('vtime ',29x,1pe15.8,22x,'rhojtv')

!........Time step control criteria

   71 FORMAT ('!      tcntrl(1)  : Courant condition for the x-hydro.')
   72 FORMAT ('!      tcntrl(2)  : density time step criterion for the x-hydro.')
   73 FORMAT ('!      tcntrl(3)  : hydro temperature change time step criterion for the x-hydro.')
   74 FORMAT ('!      tcntrl(4)  : Courant condition for the y-hydro.')
   75 FORMAT ('!      tcntrl(5)  : density time step criterion for the y-hydro.')
   76 FORMAT ('!      tcntrl(6)  : hydro temperature change time step criterion for the y-hydro.')
   77 FORMAT ('!      tcntrl(9)  : convective time step control.')
   78 FORMAT ('!      tcntrl(10) : maximum increase in the hydro time step.')
   79 FORMAT ('tcntrl',14x,i10,5x,1pe15.8,22x,'tcntrl_hydro')

!........Velocity boundary conditions

   81 FORMAT ('!      iubcjmn : inner velocity boundary codition switch.')
   82 FORMAT ('!      ubcjmn : the value of the inner velocity.')
   83 FORMAT ('!      iubcjmx : outer velocity boundary codition switch.')
   84 FORMAT ('!      ubcjmx : the value of the outer velocity.')
   85 FORMAT ('!      iuset : intermediate zone velocity boundary codition switch.')
   86 FORMAT ('!      uset(j) : the value of the velocity imposed for zone j.')
   87 FORMAT ('!      r_u : zero velocity cutoff.')
   88 FORMAT ('ubcjm ',14x,i10,5x,1pe15.8,22x,'ubcjmn')
   89 FORMAT ('ubcjm ',14x,i10,5x,1pe15.8,22x,'ubcjmx')
   90 FORMAT ('iuset ',14x,i10,42x,'iuset')
   91 FORMAT ('iuset ',14x,i10,5x,1pe15.8,22x,'uset')
   92 FORMAT ('iuset ',29x,1pe15.8,22x,'r_u')

!........Pressure boundary conditions

  101 FORMAT ('!      ipbnd : outer pressure boundary codition switch.')
  102 FORMAT ('!      pbound : the value of the outer pressure boundary condition.')
  103 FORMAT ('ipbnd ',14x,i10,5x,1pe15.8,22x,'ipbnd')

!........Mixing length convection controls

  111 FORMAT ('!      ilcnvct : mixing length convection switch.')
  112 FORMAT ('!      iemix : energy mixing switch.')
  113 FORMAT ('!      iyemix : ye mixing switch.')
  114 FORMAT ('!      ipsimix : psi0 mixing switch.')
  115 FORMAT ('!      ipcnvct : turbulent pressure switch.')
  116 FORMAT ('!      iscnvct : convective energy dissipation switch.')
  117 FORMAT ('!      alphlmix :  multiple of pressure scale height used to compute lmix.')
  118 FORMAT ('!      adrftcv : neutrino advection parameter.')
  119 FORMAT ('conv  ',14x,i10,42x,'ilcnvct')
  120 FORMAT ('conv  ',14x,i10,42x,'iemix')
  121 FORMAT ('conv  ',14x,i10,42x,'iyemix')
  122 FORMAT ('conv  ',14x,i10,42x,'ipsimix')
  123 FORMAT ('conv  ',14x,i10,42x,'ipcnvct')
  124 FORMAT ('conv  ',14x,i10,42x,'iscnvct')
  125 FORMAT ('conv  ',29x,1pe15.8,22x,'adrftcv')
  126 FORMAT ('conv  ',29x,1pe15.8,22x,'alphlmix')

!........Shock controls

  131 FORMAT ('!      ipq : pseudoviscosity switch.')
  132 FORMAT ('!      q0_x(j) : pseudoviscous pressure multiplyer for radial zone j.')
  133 FORMAT ('shock ',14x,i10,5x,1pe15.8,22x,'shock')

!........Document the dump

 1001 FORMAT (' ***Hydro keys dump written at cycle      ',i10,' on unit',i5,'***')

!........Hydro arrays...................................................
!.......................................................................

!........Header

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!........Hydro switch

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,15)
WRITE (ndump,22) ihydro

!........Newtonian - GR switch

WRITE (ndump,13)
WRITE (ndump,31)
WRITE (ndump,15)
WRITE (ndump,32) irelhy
WRITE (ndump,33) ilapsehy

!........Degeneracy - Nondegeneracy criterion

WRITE (ndump,13)
WRITE (ndump,41)
WRITE (ndump,15)
WRITE (ndump,42) degen

!........Individual time steps

WRITE (ndump,13)
WRITE (ndump,51)
WRITE (ndump,52)
WRITE (ndump,53)
WRITE (ndump,54)
WRITE (ndump,55)
WRITE (ndump,56)
WRITE (ndump,15)
WRITE (ndump,57) idtj
WRITE (ndump,58) jtv
WRITE (ndump,59) jtvshk
WRITE (ndump,60) ncyshk
WRITE (ndump,61) ncyrho
WRITE (ndump,62) rhojtv

!........Time step control criteria

WRITE (ndump,13)
WRITE (ndump,71)
WRITE (ndump,72)
WRITE (ndump,73)
WRITE (ndump,74)
WRITE (ndump,75)
WRITE (ndump,76)
WRITE (ndump,77)
WRITE (ndump,78)
WRITE (ndump,15)
WRITE (ndump,79) (i,tcntrl(i),i = 1,6)
WRITE (ndump,79) (i,tcntrl(i),i = 9,10)

!........Velocity boundary conditions

WRITE (ndump,13)
WRITE (ndump,81)
WRITE (ndump,82)
WRITE (ndump,83)
WRITE (ndump,84)
WRITE (ndump,85)
WRITE (ndump,86)
WRITE (ndump,87)
WRITE (ndump,15)
WRITE (ndump,88) iubcjmn,ubcjmn
WRITE (ndump,89) iubcjmx,ubcjmx
WRITE (ndump,90) iuset
IF ( iuset .eq. 1 ) THEN
  WRITE (ndump,91) (j,uset(j),j = imin+1,imax+1)
END IF
WRITE (ndump,92) r_u

!........Pressure boundary conditions

WRITE (ndump,13)
WRITE (ndump,101)
WRITE (ndump,102)
WRITE (ndump,15)
WRITE (ndump,103) ipbnd,pbound

!........Mixing length convection controls

WRITE (ndump,13)
WRITE (ndump,111)
WRITE (ndump,112)
WRITE (ndump,113)
WRITE (ndump,114)
WRITE (ndump,115)
WRITE (ndump,116)
WRITE (ndump,117)
WRITE (ndump,118)
WRITE (ndump,15)
WRITE (ndump,119) ilcnvct
WRITE (ndump,120) iemix
WRITE (ndump,121) iyemix
WRITE (ndump,122) ipsimix
WRITE (ndump,123) ipcnvct
WRITE (ndump,124) iscnvct
WRITE (ndump,125) adrftcv
WRITE (ndump,126) alphlmix

WRITE (ndump,13)
WRITE (ndump,131)
WRITE (ndump,132)
WRITE (ndump,15)
WRITE (ndump,133) ipq,q0_x(1)

!........Record the dump................................................
!.......................................................................

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE hydro_write_adios
