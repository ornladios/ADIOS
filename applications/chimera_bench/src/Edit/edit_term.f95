SUBROUTINE edit_term( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         edit_term
!    Module:       edit_term
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/31/05
!
!    Purpose:
!      To perform a final edit on termination of the simulation
!
!    Subprograms called:
!      bnuplot, complot, editc, edith, editma, editn, editng, editpmtr, editps,
!      editsc, editu, edity, enuvplot, lagrangeplot, lumplot, nuplot, rnuplot, nuradplot, 
!      rlagplot, roextr, rstdmp, varplot
!
!    Input arguments:
!  jr_min       : minimum radial index
!  jr_max       : maximum radial index
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  nx           : x-array extent
!  nnu          : nrutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!      kind_module, parallel_module
!      edit_module, mdl_cnfg.cmn, nu_dist_module, prb_cntl.cmn
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE parallel_module, ONLY : myid

USE edit_module, ONLY: iflprt, ncychg, nrstd1, prnttest, nprint, data_path
USE prb_cntl_module, ONLY: jnumin, jnumax

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                      :: jr_min         ! minimum radial index
INTEGER, INTENT(in)                      :: jr_max         ! maximum radial index
INTEGER, INTENT(in)                      :: ij_ray_dim     ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                      :: ik_ray_dim     ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                      :: ij_ray         ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                      :: ik_ray         ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                      :: nx             ! x-array extent
INTEGER, INTENT(in)                      :: nnu            ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=128)                       :: outfile       ! character string containing name of model files

INTEGER, PARAMETER                        :: n1 = 1        ! neutrino flavor index
INTEGER, PARAMETER                        :: n2 = 2        ! neutrino flavor index
INTEGER, PARAMETER                        :: n3 = 3        ! neutrino flavor index
INTEGER, PARAMETER                        :: n4 = 4        ! neutrino flavor index
INTEGER                                   :: istat         ! open file flag
INTEGER                                   :: nprint_save   ! temporary storage for nprint
INTEGER, PARAMETER                        :: iunit=20      ! unit number to print models

 2011 FORMAT ('Term: ij_ray=',i4,' ik_ray=',i4,' u_ge_tot=',es14.7,' u_ie_tot=',es14.7, &
& ' u_ke_tot=',es12.4,' u_ne_tot=',es12.4,' radtot=',es11.3,' e_adv=',es11.3,' u_tot=',es15.8)
 2013 FORMAT (' elecn=',es15.8,' elec_adv=',es11.3,' lep_nu=',es11.3,' lep_nu_radadv=', &
& es11.3,' totlpn=',es15.8)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                      \\\\\ FINAL EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        iflprt = 0:  no printout at termination
!        iflprt = 1:  configuration printout at termination
!        iflprt = 2:  full printout except editn at termination
!        iflprt = 3:  full printout except editng at termination
!        iflprt = 4:  full printout at termination
!-----------------------------------------------------------------------

!........Open model_term

WRITE (outfile,'(a20,i2.2,a2)') '/Models_n/model_term', &
& ij_ray_dim * ik_ray_dim * myid + ( ik_ray - 1 ) * ij_ray + ij_ray,'.d'
outfile            = TRIM(data_path)//TRIM(outfile)
OPEN (UNIT=iunit,FILE=TRIM(outfile),STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=iunit,FILE=TRIM(outfile),STATUS='old', &
& POSITION='append')

!........Save nprint

nprint_save        = nprint
nprint             = iunit

!........Print final model

SELECT CASE (iflprt)

  CASE (0)

  CASE (1)

    prnttest       = .false.
    CALL editc ( jr_min, jr_max, ij_ray, ik_ray )

  CASE (2)

    prnttest       = .false.
    CALL editp ( jr_min, jr_max )
    CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editma( jr_min, jr_max, ij_ray, ik_ray )
    CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editps( jr_min, jr_max, ij_ray, ik_ray )
    CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
    CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
    CALL editng( n1, jnumin, jnumax, ij_ray, ik_ray)
    CALL editng( n2, jnumin, jnumax, ij_ray, ik_ray)
    CALL editng( n3, jnumin, jnumax, ij_ray, ik_ray)
    CALL editng( n4, jnumin, jnumax, ij_ray, ik_ray)

  CASE (3)

    prnttest       = .false.
    CALL editp ( jr_min, jr_max )
    CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editma( jr_min, jr_max, ij_ray, ik_ray )
    CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editps( jr_min, jr_max, ij_ray, ik_ray )
    CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
    CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
    CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray)
    CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray)
    CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray)
    CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray)

  CASE (4)

    prnttest    = .false.
    CALL editp ( jr_min, jr_max )
    CALL editc ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editma( jr_min, jr_max, ij_ray, ik_ray )
    CALL edith ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editps( jr_min, jr_max, ij_ray, ik_ray )
    CALL editu ( jr_min, jr_max, ij_ray, ik_ray )
    CALL edity ( jr_min, jr_max, ij_ray, ik_ray )
    CALL editsc( jr_min, jr_max, ij_ray, ik_ray )
    CALL editn ( n1, jnumin, jnumax, ij_ray, ik_ray)
    CALL editng( n1, jnumin, jnumax, ij_ray, ik_ray)
    CALL editn ( n2, jnumin, jnumax, ij_ray, ik_ray)
    CALL editng( n2, jnumin, jnumax, ij_ray, ik_ray)
    CALL editn ( n3, jnumin, jnumax, ij_ray, ik_ray)
    CALL editng( n3, jnumin, jnumax, ij_ray, ik_ray)
    CALL editn ( n4, jnumin, jnumax, ij_ray, ik_ray)
    CALL editng( n4, jnumin, jnumax, ij_ray, ik_ray)

END SELECT

!........Restore nprint

nprint             = nprint_save

!........Close model_term

CLOSE (UNIT=nrstd1,STATUS='keep')

RETURN
END SUBROUTINE edit_term
