SUBROUTINE read_nuclear_data( data_dir, data_desc )
!-----------------------------------------------------------------------------  
!  This routine reads, from the file netwinv, the nuclei included along 
!  with the nuclear data which will be needed for later calculations.  This 
!  data includes the atomic number, the number of protons and neutrons, and 
!  the binding energy (calculated from the tabulated mass excess).  Also the 
!  tabulations of the  partition functions, g, are read in for later 
!  interpolation. Once the set of nuclear data is read in, it is assigned 
!  to the proper nuclei.
!-----------------------------------------------------------------------------  

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE nuclear_data, ONLY : nname, aa, zz, nn, be
USE nuc_number, ONLY : ny
USE part_funct_data, ONLY : g, gg, t9i, angm

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (LEN=*),  INTENT(in)  :: data_dir

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

CHARACTER (LEN=80), INTENT(out) :: data_desc

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=5)       :: nam                                 

INTEGER                 :: i,l,n,m,na,nb                                     
INTEGER                 :: it9i(24)

REAL(KIND=double)       :: mex,a,sp 

!-----------------------------------------------------------------------
!  character (LEN=80) :: data_dir,data_desc
!  data_dir='../../FullNet/Data_alpha'
!
!  Read in the data description
!-----------------------------------------------------------------------

OPEN (13,file=TRIM(data_dir)//"/net_desc",STATUS='old')
READ (13,"(A)") data_desc
OPEN (12,FILE=TRIM(data_dir)//"/netwinv",STATUS='old')
READ (12,"(i5)") ny

!-----------------------------------------------------------------------
!  Read in the partition function iteration grid, and fix endpoints
!-----------------------------------------------------------------------

READ (12,"(24i3)") (it9i(i),i=1,24)

DO i = 1,24
  t9i(i)      = it9i(i) * 0.01d0
END DO
t9i(24)       = t9i(24) * 10.d0

DO  n= 1,ny
  READ (12,"(a5)") nam
END DO

!-----------------------------------------------------------------------
!  Set size of nuclear parameter arrays and read in nuclear parameters 
!  and partition function interpolation table. nname(0), gg(0) and angm(0) 
!  are placeholders for non-nuclei.
!-----------------------------------------------------------------------

ALLOCATE (nname(0:ny))
ALLOCATE (aa(ny),zz(ny),nn(ny),be(ny))
ALLOCATE (g(24,ny),gg(0:ny),angm(0:ny))

aa            = zero
zz            = zero
nn            = zero
be            = zero
g             = zero
gg            = zero
angm          = zero

nname(0)      = '   '
angm(0)       = 0.0d0
DO l = 1,ny
  READ (12,"(a5,f12.3,2i4,f6.1,f10.3)") nam,a,na,nb,sp,mex
  READ (12,"(8f9.2)") (g(m,l),m=1,24)
  aa(l)       = a
  zz(l)       = dble(float(na))
  nn(l)       = dble(float(nb))
  be(l)       = 8.07144 * nn(l) + 7.28899 * zz(l) - mex
  angm(l)     = 2.d0 * sp + 1.d0
  nname(l)    = nam
END DO

RETURN
END
