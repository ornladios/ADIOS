SUBROUTINE netmatr( kstep )
!===============================================================================
!  This routine calculates the Jacobian matrix dYdot/dY, and solves for the 
!  Newton-Raphson iteration, dy.
!===============================================================================

USE kind_module, ONLY : double

USE abundances, ONLY : y, yt, dy, ydot
USE conditions, ONLY : tdel
USE controls, ONLY : idiag
USE edit_module, ONLY : nlog
USE nuclear_data, ONLY : nname
USE nuc_number, ONLY : ny
USE reac_rate_data, ONLY : la, le, n11, n21, n22, n31, n32, n33, b1, b2, b3

IMPLICIT none
SAVE

CHARACTER (len=10)               :: var_name

LOGICAL :: first = .true.

INTEGER, DIMENSION(:), allocatable :: indx
INTEGER :: i, j, kstep, i0, la1, le1, la2, le2, la3, le3, j1, l1, l2, l3, info
INTEGER                          :: istat         ! allocation status

REAL(KIND=double), DIMENSION(:),   ALLOCATABLE :: f, ydot0, um
REAL(KIND=double), DIMENSION(:,:), ALLOCATABLE :: am ! the Jacobian Matrix
REAL(KIND=double) :: alph, beta, rdt

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in netmatr')
 2001 FORMAT (' Deallocation problem for array ',a10,' in netmatr')

ALLOCATE (indx(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'indx      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ydot0(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ydot0     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (um(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'um        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (am(ny,ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'am        '; WRITE (nlog,1001) var_name; END IF

!If ( .not. allocated(f)  ) ALLOCATE(indx(ny), f(ny), ydot0(ny), um(ny))
!If ( .not. allocated(am) ) ALLOCATE(am(ny,ny))

alph            = 0.0d0
beta            = 1.0d0
ydot0           = 0.0d0

!-----------------------------------------------------------------------
!  Calculate the reaction rates and abundance time derivatives
!-----------------------------------------------------------------------

CALL yderiv

!-----------------------------------------------------------------------
!  Build the Jacobian, row by row
!-----------------------------------------------------------------------

rdt             = 1.0d0/tdel/beta 
DO i0=1,ny 
  um            = 0.0d0
  um(i0)        = rdt  
  la1           = la(1,i0) 
  le1           = le(1,i0)  
  DO j1=la1,le1  
    l1          = n11(j1)
    um(l1)      = um(l1) - b1(j1)
  END DO 
  la2           = la(2,i0) 
  le2           = le(2,i0)  
  DO j1=la2,le2  
    l1          = n21(j1)
    l2          = n22(j1) 
    um(l1)      = um(l1)-b2(j1) * yt(l2) 
    um(l2)      = um(l2)-b2(j1) * yt(l1)  
  END DO       
  la3           = la(3,i0) 
  le3           = le(3,i0)  
  DO j1=la3,le3 
    l1          = n31(j1) 
    l2          = n32(j1) 
    l3          = n33(j1)
    um(l1)      = um(l1) - b3(j1) * yt(l2) * yt(l3)    
    um(l2)      = um(l2) - b3(j1) * yt(l1) * yt(l3)     
    um(l3)      = um(l3) - b3(j1) * yt(l1) * yt(l2)      
  END DO                 

!-----------------------------------------------------------------------------  
!  Tranfer to matrix row
!-----------------------------------------------------------------------------  

  am(i0,:)=um

!-----------------------------------------------------------------------------  
!  at(:,i0)=um ! or column if the solver wants the transpose
!-----------------------------------------------------------------------------  

END DO                                                      

!-----------------------------------------------------------------------------  
!  am = transpose(at)
!-----------------------------------------------------------------------------  

!-----------------------------------------------------------------------------  
!  Calculate equation to zero
!-----------------------------------------------------------------------------  

f               = ( y - yt ) * rdt + ydot
!     f=y*rdt-yt*rdt+ydot+alph*ydot0/beta  

IF ( idiag >= 4 ) THEN
  WRITE (50,"(a2,i5,es14.7)") 'F',kstep,rdt
  DO i=1,ny
    WRITE (50,"(a5,4es17.9)") nname(i),f(i),ydot(i),yt(i),y(i)
    WRITE (50,"(5es16.8)") (am(i,j),j=1,ny)
  END DO
END IF

!  Test the eigenvalues 
!     If(idiag>=6) Then
!       call eigen_test(kstep,am,rdt)
!     Endif
!-----------------------------------------------------------------------------  
!  The bulk of the computational cost of the network (60-95%) is the solving 
!  of the matrix equation.  Careful selection of the matrix solver is therefore 
!  very important to fast computation.  Generally, hand tuned solvers such as 
!  those supplied by the hardware manufacturer or third-parties like NAG, IMSL,
!  etc. are the fastest.  However for portability, by default we USE Numerical 
!  Recipes routines.  
!-----------------------------------------------------------------------------  

!  Use Num Rec LU Decomp.
!     call ludcmp(am,ny,ny,indx,d)
!     call lubksb(am,ny,ny,indx,f)
!     dy=f

!  Use LAPACK solver
!     call sgesv(ny,1,am,ny,indx,f,ny,info) ! Single precision version
CALL dgesv(ny,1,am,ny,indx,f,ny,info) ! Double precision version
dy              = f

!  Diagnostic output
!     If(idiag>=4) Then
!       Write(50,"(a5,3es12.4)") (nname(i),dy(i),y(i),yt(i),i=1,ny)
!     Endif
DEALLOCATE (indx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'indx      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ydot0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ydot0     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (um, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'um        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (am, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'am        '; WRITE (nlog,2001) var_name; END IF

RETURN                                                                    
END SUBROUTINE netmatr                                                                       
