SUBROUTINE edit_e( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         edit_e
!    Module:       edit_e
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/25/96
!
!    Purpose:
!      To edit kinetic, potential, gravitational, and total energy as given by EVH1.
!
!    Subprograms called:
!      date_and_time_print
!
!    Input arguments:
!  jr_min     : inner radial zone of region for which configuration edit is to be made.
!  jr_max     : outer radial zone of region for which configuration edit is to be made.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!  prnttest   : true  - test to see IF printing criteria is satisfied.
!               false - bypass printing criteria test.
!  iprint     : 0    - do not print to print file.
!               ne 0 - print to print file.
!  nprint     : unit number of print file.
!  iplot      : 0    - do not print to plot file.
!               ne 0 - print to plot file.
!  nplot      : unit number of plot file.
!  nede(i)    : edit_e counter for data set i.
!  intede(i)  : number of cycles between edits of data set i.
!  idxede(i)  : edit jr_min, jr_max, and every idxedc(i) radial zone between them for data set i.
!
!    Include files:
!      kind_module, numerical_module
!      boundary_module, edit_module, eos_snc_x_module, hydro_module,
!      mdl_cnfg_module, nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nx
USE numerical_module, ONLY : epsilon

USE edit_module, ONLY : nprint, prnttest, head, intede, nede, idxede
USE mdl_cnfg_module, ONLY : dmrst, r

USE evh1_sweep, ONLY : e, ei, ekin, egrav, xa, rho=>r

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index
INTEGER                          :: jv            ! do index
INTEGER                          :: n_time        ! iteration index

REAL(KIND=double), DIMENSION(nx) :: ei_m          ! Internal energy enclosed (ergs)
REAL(KIND=double), DIMENSION(nx) :: eg_m          ! Gravitational energy enclosed (ergs)
REAL(KIND=double), DIMENSION(nx) :: ek_m          ! Kinetic energy enclosed (ergs)
REAL(KIND=double), DIMENSION(nx) :: et_m          ! Total energy enclosed (ergs)

    1 FORMAT (/)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
  101 FORMAT (48x,'EVH1 Energy Data')
  103 FORMAT (47x,18('-')/)
  105 FORMAT ('   j e_internal   e_grav   e_kinetic   e_total   e_int_encl e_grv_encl e_kin_encl e_tot_encl'/)
  107 FORMAT (1x,i4,8es11.4)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!**********************************************************************!        
!                                                                      !
!                  EVH1 energy data                                    !
!                                                                      !
!----------------------------------------------------------------------!
! ei          'e_internal' Internal energy (ergs/gm)                   !
! egrav       'e_grv'      Gravitational energy (ergs/gm)              ! 
! ekin        'e_kinetic'  Kinetic energy (ergs/gm)                    !
! e           'e_total'    Total energy (ergs/gm)                      !
! ei_m        'e_int_encl' Internal energy enclosed (ergs)             !
! eg_m        'e_grv_encl' Gravitational energy enclosed (ergs)        !        
! ek_m        'e_kin_encl' Kinetic energy enclosed (ergs)              !
! et_m        'e_tot_encl' Total energy enclosed (ergs)                !
!----------------------------------------------------------------------!
!        Print header                                                  !
!----------------------------------------------------------------------!

n_time             = nprint
IF ( prnttest ) THEN
  nede(1)          = nede(1) + 1
  IF ( nede(1) < intede(1) ) RETURN
  nede(1)          = 0
END IF !  prnttest

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print(n_time)
WRITE (nprint,101)
WRITE (nprint,103)
WRITE (nprint,105)

!----------------------------------------------------------------------!
!  Determine egrav
!----------------------------------------------------------------------!

xa(jr_min+5:jr_max+6) = r(jr_min:jr_max+1)

CALL e_potential( xa, rho ) 

!----------------------------------------------------------------------!
!  Enclosed energy data
!----------------------------------------------------------------------!

ei_m(jr_min)       = ei   (jr_min+5) * dmrst(jr_min)
eg_m(jr_min)       = egrav(jr_min+5) * dmrst(jr_min)
ek_m(jr_min)       = ekin (jr_min+5) * dmrst(jr_min)
et_m(jr_min)       = e    (jr_min+5) * dmrst(jr_min)

DO j = jr_min+1,jr_max
  ei_m(j)          = ei_m(j-1) + ei   (j+5) * dmrst(j)
  eg_m(j)          = eg_m(j-1) + egrav(j+5) * dmrst(j)
  ek_m(j)          = ek_m(j-1) + ekin (j+5) * dmrst(j)
  et_m(j)          = et_m(j-1) + ( ei(j+5) + egrav(j+5) + ekin(j+5) ) * dmrst(j)
END DO

!----------------------------------------------------------------------!
!  Print specific and enclosed energy data
!----------------------------------------------------------------------!

DO jv = jr_min,jr_max,idxede(1)
  j                = jr_max - jv + jr_min
  WRITE (nprint,107) j, ei(j+5), egrav(j+5), ekin(j+5), e(j+5), ei_m(j), &
&  eg_m(j), ek_m(j), et_m(j)
END DO

RETURN
END SUBROUTINE edit_e
