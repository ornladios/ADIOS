SUBROUTINE unpack_arrays( c_init_data, c_radhyd_data, c_eos_data,        &
& c_nuc_data, i_init_data, i_radhyd_data, i_trans_data, i_e_advct_data,  &
& i_edit_data, i_hydro_data, i_nuc_data, i_model_data, l_rezone_data,    &
& d_radhyd_data, d_eos_data, d_trans_data, d_e_advct_data, d_edit_data,  &
& d_hydro_data, d_nuc_data, d_model_data1, d_model_data2, d_model_data3, &
& d_psi_data1, d_psi_data2, d_psi_data3, d_psi_data4, d_psi_data5,       &
& d_rezone_data, l_multi_d )
!-----------------------------------------------------------------------
!
!    File:         unpack_arrays
!    Module:       unpack_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To unpack arrays and put vartiables in appropriate modules.
!
!    Subprograms called:
!  unpack_init            : unpacks init data
!  unpack_radial_ray_keys : unpacks radhyd_ray_keys
!  unpack_rezone_arrays   : unpacks rezone keys
!  unpack_eos_keys        : unpacks eos keys
!  unpack_transport_keys  : unpacks transport keys
!  unpack_e_advct_keys    : unpacks e_advection keys
!  upack_edit_keys        : unpacks edit keys
!  unpack_hydro_keys      : unpacks hydro keys
!  unpack_nuclear_data    : unpacks nuclear keys and abundance data
!  unpack_initial_model   : unpacks the initial model
!
!    Input arguments:
!        none
!
!    Output arguments:
!  c_init_data            : character array of initial data
!  c_radhyd_data          : character array of radhyd keys
!  c_eos_data             : character array of edit keys
!  c_nuc_data             : character array of nuclei
!  l_rezone_data          : lagrangian flag
!  i_init_data            : integer array of initial data
!  i_radhyd_data          : integer array of radhyd keys
!  i_trans_data           : integer array of transport keys
!  i_e_advct_data         : integer array of e_advect keys
!  i_edit_data            : integer array of transport keys
!  i_hydro_data           : integer array of transport keys
!  i_nuc_data             : integer array of edit keys
!  i_model_data           : integer array of initial model data
!  d_radhyd_data          : 64 bit real array of radhyd keys
!  d_eos_data             : 64 bit real array of edit keys
!  d_trans_data           : 64 bit real array of transport keys
!  d_e_advct_data         : 64 bit real array of e_advect keys
!  d_edit_data            : 64 bit real array of edit keys
!  d_hydro_data           : 64 bit real array of hydro keys
!  d_nuc_data             : 64 bit real array of nuclear keys
!  d_model_data1          : 64 bit real array of initial model data
!  d_model_data2          : 64 bit real array of initial model data
!  d_model_data3          : 64 bit real array of initial model data
!  d_psi_data1            : 64 bit real array of neutrino data
!  d_psi_data2            : 64 bit real array of neutrino data
!  d_psi_data3            : 64 bit real array of neutrino data
!  d_psi_data4            : 64 bit real array of neutrino data
!  d_psi_data5            : 64 bit real array of neutrino data
!  d_rezone_data          : 64 bit real array of rezoned model data
!
!    Include files:
!  array_module, kind_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE array_module, ONLY : nx, ny, nz, nez, nezp1, nnu, nnc, ij_ray_dim, &
& ik_ray_dim
USE kind_module, ONLY : double

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER(len=128), INTENT(in), DIMENSION(1)   :: c_init_data    ! character array of initial data
CHARACTER (len=2), INTENT(in), DIMENSION(20)   :: c_radhyd_data  ! character array of radhyd keys
CHARACTER (len=1), INTENT(in), DIMENSION(1)    :: c_eos_data     ! character array of edit keys
CHARACTER (len=5), INTENT(in), DIMENSION(nnc)  :: c_nuc_data     ! character array of nuclei

LOGICAL, DIMENSION(1), INTENT(in)              :: l_rezone_data  ! lagrangian flag
LOGICAL                                        :: l_multi_d      ! multiD flag

INTEGER, INTENT(in), DIMENSION(2)              :: i_init_data    ! integer array of initial data
INTEGER, INTENT(in), DIMENSION(50)             :: i_radhyd_data  ! integer array of radhyd keys
INTEGER, INTENT(in), DIMENSION(40+2*nnu)       :: i_trans_data   ! integer array of transport keys
INTEGER, INTENT(in), DIMENSION(5)              :: i_e_advct_data ! integer array of e_advect keys
INTEGER, INTENT(in), DIMENSION(1200+3*40*nnu)  :: i_edit_data    ! integer array of edit keys
INTEGER, INTENT(in), DIMENSION(30)             :: i_hydro_data   ! integer array of transport keys
INTEGER, INTENT(in), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim)   :: i_nuc_data     ! integer array of edit keys
INTEGER, INTENT(in), DIMENSION(2)              :: i_model_data   ! integer array of initial model data

REAL(KIND=double), INTENT(in), DIMENSION(50)                  :: d_radhyd_data  ! 64 bit real array of radhyd keys
REAL(KIND=double), INTENT(in), DIMENSION(14)                  :: d_eos_data     ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION((110+3*nnu+3*nez+1)) :: d_trans_data   ! 64 bit real array of transport keys
REAL(KIND=double), INTENT(in), DIMENSION(5+2*nnu)             :: d_e_advct_data ! 64 bit real array of e_advect keys
REAL(KIND=double), INTENT(in), DIMENSION(50)                  :: d_edit_data    ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION((30+nx))             :: d_hydro_data   ! 64 bit real array of hydro keys
REAL(KIND=double), INTENT(in), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data     ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)                :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(2,nx)                                      :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(2,nx+1)                                    :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(8,nez,nnu,ij_ray_dim,ik_ray_dim)           :: d_psi_data1   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim)        :: d_psi_data2   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,nx,nnu,ij_ray_dim,ik_ray_dim)            :: d_psi_data3   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,nnu,ij_ray_dim,ik_ray_dim)               :: d_psi_data4   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,ij_ray_dim,ik_ray_dim)                   :: d_psi_data5   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(20+3*ny+3*nz+2)      :: d_rezone_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                        :: nouttmp        ! reatart read flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' init variables have been unpacked')
 103 FORMAT (' radhyd_keys have been unpacked')
 105 FORMAT (' radial_ray_keys have been unpacked')
 107 FORMAT (' rezone_arrays have been unpacked')
 109 FORMAT (' eos_keys have been unpacked')
 111 FORMAT (' transport_keys have been unpacked')
 113 FORMAT (' e_advct_keys have been unpacked')
 115 FORMAT (' edit_keys have been unpacked')
 117 FORMAT (' hydro_keys have been unpacked')
 119 FORMAT (' nuclear_keys have been unpacked')
 121 FORMAT (' initial_model has been unpacked')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Unpack init data
!-----------------------------------------------------------------------

CALL unpack_init( c_init_data, i_init_data )
IF ( myid == 0 ) WRITE (nlog,101)

!-----------------------------------------------------------------------
!  Unpack radial_ray_keys
!-----------------------------------------------------------------------

CALL unpack_radial_ray_keys( c_radhyd_data, i_radhyd_data, d_radhyd_data )
IF ( myid == 0 ) WRITE (nlog,105)

!-----------------------------------------------------------------------
!  Unpack rezone_arrays
!-----------------------------------------------------------------------

CALL unpack_rezone_arrays( l_rezone_data, d_rezone_data )
IF ( myid == 0 ) WRITE (nlog,107)

!-----------------------------------------------------------------------
!  Unpack eos_keys
!-----------------------------------------------------------------------

CALL unpack_eos_keys( c_eos_data, d_eos_data )
IF ( myid == 0 ) WRITE (nlog,109)

!-----------------------------------------------------------------------
!  Unpack ransport_keys
!-----------------------------------------------------------------------

CALL unpack_transport_keys( nez, nezp1, nnu, i_trans_data, d_trans_data )
IF ( myid == 0 ) WRITE (nlog,111)

!-----------------------------------------------------------------------
!  Unpack e_advct_keys
!-----------------------------------------------------------------------

CALL unpack_e_advct_keys( nnu, i_e_advct_data, d_e_advct_data )
IF ( myid == 0 ) WRITE (nlog,113)

!-----------------------------------------------------------------------
!  Unpack edit_keys
!-----------------------------------------------------------------------

CALL upack_edit_keys( ij_ray_dim, ik_ray_dim, nez, nnu, i_edit_data, &
& d_edit_data )
IF ( myid == 0 ) WRITE (nlog,115)

!-----------------------------------------------------------------------
!  Unpack hydro_keys
!-----------------------------------------------------------------------

CALL unpack_hydro_keys( nx, i_hydro_data, d_hydro_data )
IF ( myid == 0 ) WRITE (nlog,117)

!-----------------------------------------------------------------------
!  Unpack muclear_data
!-----------------------------------------------------------------------

CALL unpack_nuclear_data( nx, nnc, ij_ray_dim, ik_ray_dim, c_nuc_data, &
& i_nuc_data, d_nuc_data )
IF ( myid == 0 ) WRITE (nlog,119)

!-----------------------------------------------------------------------
!  Unpack model
!
!  l_multi_d = false and nouttmp /= 5 : 
!   (1) simulation initiated from cycle 0              or
!   (2) 1-D simulation initiated from cycle > 0        or
!   (3) M-D simulation initiated from a 1-D precursor
!  l_multi_d = true and nouttmp /= 6 :
!   M-D simulation restarted from binary restart files
!  l_multi_d = true and nouttmp /= 6 :
!   M-D simulation restarted from hdf restart files
!-----------------------------------------------------------------------

nouttmp            = i_init_data(2)

IF ( .not. l_multi_d ) THEN
  IF ( nouttmp /= 5 ) THEN
    CALL unpack_initial_model( nx, nez, nnu, ij_ray_dim, ik_ray_dim, &
&    i_model_data, d_model_data1, d_model_data2, d_model_data3, d_psi_data2 )
    IF ( myid == 0 ) WRITE (nlog,121)
  ELSE ! nouttmp = 5
    CALL unpack_restart_model_to_node( nx, nez, nnu, nnc, ij_ray_dim,       &
&    ik_ray_dim, i_model_data, d_model_data1, d_model_data2, d_model_data3, &
&    d_psi_data1, d_psi_data2, d_psi_data3, d_psi_data4, d_psi_data5,       &
&    i_nuc_data, d_nuc_data )
  END IF ! nouttmp /= 5
END IF ! .not. l_multi_d

RETURN
END SUBROUTINE unpack_arrays
