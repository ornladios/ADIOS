MODULE evh1_bound
!=======================================================================
!  Boundary values.  Names conform to sweep module
!  _bcl is left (inner) boundary, _bcr is right (outer) boundary
!=======================================================================

USE kind_module, ONLY : double

REAL(KIND=double)                            :: u_bcl    ! left boundary x-velocity
REAL(KIND=double)                            :: v_bcl    ! left boundary y-velocity
REAL(KIND=double)                            :: w_bcl    ! left boundary z-velocity
REAL(KIND=double)                            :: r_bcl    ! left boundary density
REAL(KIND=double)                            :: p_bcl    ! left boundary pressure
REAL(KIND=double)                            :: ei_bcl   ! left boundary density
REAL(KIND=double)                            :: ye_bcl   ! left boundary electron fraction
REAL(KIND=double)                            :: temp_bcl ! left boundary temperature
REAL(KIND=double)                            :: gc_bcl   ! left boundary gamma
REAL(KIND=double)                            :: ge_bcl   ! left boundary gamma
REAL(KIND=double)                            :: psi0_bcl ! left boundary psi0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: comp_bcl ! left boundary composition

REAL(KIND=double)                            :: u_bcr    ! right boundary x-velocity
REAL(KIND=double)                            :: v_bcr    ! right boundary y-velocity
REAL(KIND=double)                            :: w_bcr    ! right boundary z-velocity
REAL(KIND=double)                            :: r_bcr    ! right boundary density
REAL(KIND=double)                            :: p_bcr    ! right boundary pressure
REAL(KIND=double)                            :: ei_bcr   ! right boundary density
REAL(KIND=double)                            :: ye_bcr   ! right boundary electron fraction
REAL(KIND=double)                            :: temp_bcr ! right boundary temperature
REAL(KIND=double)                            :: gc_bcr   ! right boundary gamma
REAL(KIND=double)                            :: ge_bcr   ! right boundary gamma
REAL(KIND=double)                            :: psi0_bcr ! right boundary psi0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: comp_bcr ! right boundary composition

END MODULE evh1_bound
