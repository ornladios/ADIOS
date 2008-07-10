SUBROUTINE full_net( izone, dtime )
!===============================================================================
!  The abundances are evolved using a Newton-Raphson iteration scheme to solve 
!  the equation yt(i)-y(i)/tdel = ydot(i), where y(i) is the abundance at the 
!  beginning of the iteration, yt(i) is the trial abundance at the end of the 
!  timestep, ydot(i) is the time derivative of the trial abundance  calculated 
!  from reaction rates and tdel is the timestep.  
!===============================================================================

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one

USE abundances, ONLY : y, yt, dy, yo, ydot
USE conditions, ONLY : t, tt, tdel, t9t, rhot
USE controls, ONLY : kstmx, knrmx, iconvc, idiag, itso, tolc, tolm, tdelmm, ymin
USE edit_module, ONLY : nlog
USE nuclear_data, ONLY : nname, aa, zz
USE nuc_number, ONLY : ny
USE thermo_data, ONLY : tstart, tstop

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)           :: izone

REAL(KIND=double), INTENT(in) :: dtime

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER            :: irdymx(1), idymx(1)
INTEGER            :: k, kstep, kts, ktsmx, knr, kout, idiag0

REAL(KIND=double)  :: reldy(ny)
REAL(KIND=double)  :: toln      !  Network convergence condition
REAL(KIND=double)  :: testc, testc2, testm, testn    ! Convergence tests
REAL(KIND=double)  :: ta, tdela, xtot, xtoto, enm, enb, enold, en0, edot, ye
REAL(KIND=double)  :: ytot, ztot, atot

!-----------------------------------------------------------------------
!  Set reaction controls not read in from control
!-----------------------------------------------------------------------

idiag0             = idiag
ktsmx              = 10
kstep              = 0
kout               = 0

!-----------------------------------------------------------------------
!  Set tdel to zero so that time step will be computed from real
!   derivatives rather than abundance changes
!-----------------------------------------------------------------------

tdel               = zero

!-----------------------------------------------------------------------
!  Normalize initial abundances and change units if necessary
!-----------------------------------------------------------------------

yt                 = y
!     call norm( yt ) 

!-----------------------------------------------------------------------
!  Calculate the total energy of the nuclei
!-----------------------------------------------------------------------

CALL benuc( yt, enb, enm, ytot, ztot, atot )
xtot               = atot - one
en0                = enm
edot               = zero

!-----------------------------------------------------------------------
!  Start evolution 
!-----------------------------------------------------------------------

t                  = tstart
tt                 = tstart
CALL t9rhofind( kstep )
IF ( itso > 0 ) CALL ts_output( kstep, (enm-en0), edot, kts, knr )

!-----------------------------------------------------------------------------  
!  For each step, an initial guess for the timestep is calculated, based on 
!  the abundance variations of the previous step if available.  
!-----------------------------------------------------------------------------  

STEP: DO kstep = 1,kstmx
  CALL timestep( kstep, dtime )
  If( idiag >= 1 ) WRITE(50,*) 'TDel',tt,tdel
  xtoto            = xtot

!-----------------------------------------------------------------------
!  Determine if this is an output step
!-----------------------------------------------------------------------
  idiag            = idiag0

!-----------------------------------------------------------------------------  
!  For each trial timestep, tdel, the Newton-Raphson iteration is attempted.
!  The possible results for a timestep are a converged set of yt or a 
!  failure to converge.  If convergence fails, iteration is retried with 
!  the trial timestep reduced by tdelmm, up to ktsmx times.
!-----------------------------------------------------------------------------  

  TS: DO kts = 1,ktsmx

!-----------------------------------------------------------------------
!  Calculate the thermodynamic factors necessary for reaction rates, 
!  including screening, and the reaction rates.
!-----------------------------------------------------------------------

    CALL cross_sect 

!-----------------------------------------------------------------------
!  The Newton-Raphson iteration occurs for at most knrmx iterations.  
!-----------------------------------------------------------------------

    NR: DO knr = 1,knrmx
                   
!-----------------------------------------------------------------------
!  Calculate the changes in abundaces, dy
!-----------------------------------------------------------------------

      CALL netmatr( kstep )

!-----------------------------------------------------------------------
!  Evolve the abundances and calculate convergence tests
!-----------------------------------------------------------------------

      yt           = yt + dy
      WHERE ( yt < ymin ) 
        yt         = 0.0d0
        reldy      = 0.0d0
      ELSEWHERE
        reldy      = DABS(dy/yt)
      END WHERE

      IF ( idiag >= 3 ) THEN
        irdymx     = maxloc(reldy)
        idymx      = maxloc(dy)
        WRITE (50,"(a2,i5,2i3,2(a5,2es12.4))") &
&              'dY',kstep,kts,knr,nname(idymx(1)),&
&              dy(idymx(1)),y(idymx(1)),nname(irdymx(1)),&
&              reldy(irdymx(1)),y(irdymx(1))
        IF ( idiag >= 4 ) WRITE (50,"(a5,5es12.4)") &
&              (nname(k),yt(k),dy(k),reldy(k),(aa(k)*dy(k)),&
&              (aa(k)*yt(k)),k=1,ny)
      END IF ! idiag >= 3

!-----------------------------------------------------------------------------  
!  There are 3 included convergence tests: testc, which measures relative
!  changes, testc2 which measures total abundance changes, and testm
!  which tests mass conservation.  
!-----------------------------------------------------------------------------  

      testc        = SUM(reldy)
      testc2       = SUM(aa*dy)
      xtot         = SUM(aa*yt) - 1.0d0
      testm        = xtot - xtoto

      IF ( idiag >= 2 ) WRITE (50,"(a3,i5,i3,3es14.6)") &
&          'KNR',kstep,knr,testm,testc,testc2

!-----------------------------------------------------------------------------  
!  testc is the most stringent test, and requires the most iterations.  
!  testm is the most lax, and therefore the fastest, often requiring only one 
!  iteration.  Considering the uncertainties of the reaction rates, it is
!  doubtful that the increased precision of testc is truly increased
!  accuracy. 
!-----------------------------------------------------------------------------  

!-----------------------------------------------------------------------
!  Ordinarily, test for true convergence
!-----------------------------------------------------------------------

      IF ( iconvc /= 0  .or.  tt >= tstop ) THEN
        testn      = testc
        toln       = tolc
  
!-----------------------------------------------------------------------
!  Otherwise, use mass conservation for convergence condition 
!-----------------------------------------------------------------------

      ELSE
        testn      = testm
        toln       = tolm           
      END IF ! iconvc /= 0  .or.  tt >= tstop
        IF ( DABS(testn) <= toln)  Exit TS

    END DO NR

!-----------------------------------------------------------------------
!  If convergence is not achieved in knrmx iterations, reset abundances 
!  and try again with the timestep reduced.  
!-----------------------------------------------------------------------

    IF ( idiag >= 1 ) WRITE(50,*) 'TS Failure',knr,kts,xtot,testn
    tdel           = tdel/tdelmm
    tt             = t + tdel
    CALL t9rhofind( kstep )
    yt             = y

  END DO TS

!-----------------------------------------------------------------------
!  If convergence is successful, update time and abundances 
!-----------------------------------------------------------------------

  IF ( kts < ktsmx ) THEN
    IF ( idiag >= 1 ) WRITE (50,"(a4,i5,i3,3es12.4)") &
&        'Conv',kstep,knr,xtot,testn,toln
    ta             = t
    tdela          = tdel
    t              = tt

!-----------------------------------------------------------------------
!     call norm(y)
!-----------------------------------------------------------------------

    yo             = y
    y              = yt
    ye             = SUM( zz * y )

    IF ( idiag >= 2 ) THEN
      WRITE (50,"(a)") 'delta Y'
      WRITE (50,"(a5,4es12.4)") (nname(k),y(k),yo(k),&
     &        (yt(k)-y(k)),(tdel*ydot(k)),k=1,ny)
    END IF ! idiag >= 2

    enold          = enm
    CALL benuc( yt, enb, enm, ytot, ztot, atot )
    edot           = ( enm - enold )/tdel

    IF ( idiag >= 1 ) THEN 
      WRITE (50,"(i5,5es14.7)") kstep,t,tdel,t9t,rhot,ye
      WRITE (50,"(5(a5,es11.4))") (nname(k),y(k),k=1,ny)
    END IF ! idiag >= 1

    IF ( itso > 0 ) CALL ts_output( kstep, ( enm - en0 ), edot, kts, knr )
    IF ( t >= tstop ) THEN
      EXIT STEP
    END IF ! t >= tstop
     
!-----------------------------------------------------------------------
!  If reduced timesteps fail to yield convergence, warn and exit
!-----------------------------------------------------------------------

  ELSE
    WRITE (nlog,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,knr,kts
    WRITE (nlog,*) 'Timestep retrys fail after ',kts,' attempts'
    EXIT STEP
  END IF ! kts < ktsmx

END DO STEP

!-----------------------------------------------------------------------
!  Test that the stop time is reached
!-----------------------------------------------------------------------

IF ( t < tstop ) THEN
  WRITE (50,"(a,3es12.4)") 'Evolution incomplete!!!',t,tstop,tdel 
  WRITE (nlog,"(a,es12.4,a,es12.4,a)") 'Evolution stopped at time=',t,&
&         '.  Stop time (',tstop,') not reached!' 
  WRITE (nlog,"(a,i6)") &
&          'Approximately',int((tstop-t)/tdel),'more steps needed' 
  WRITE (nlog,"(a,i6)") &
&          'izone =',izone
END IF ! t < tstop

!-----------------------------------------------------------------------
!  End Post Processing cycle
!-----------------------------------------------------------------------

CALL final_output(kstep)

RETURN
END SUBROUTINE full_net
