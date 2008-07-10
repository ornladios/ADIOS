SUBROUTINE unpack_nuclear_data( nx, nnc, ij_ray_dim, ik_ray_dim, c_nuc_data, &
& i_nuc_data, d_nuc_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_nuclear_data
!    Module:       unpack_nuclear_data
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/23/04
!
!    Purpose:
!      To unpack the nuclear key arrays and restore the values to the
!       appropriate variables in nucbrn_module and radhyd_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : radial array dimension
!  nnc        : nuclear specie array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  c_nuc_data : character array of nuclei
!  i_nuc_data : integer array of nuclear keys
!  d_nuc_data : real*8 array of nuclear keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  nucbrn_module, radial_ray_module, eos_snc_x_module, eos_snc_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE nucbrn_module, ONLY : inuc, itnuc, ttolnuc, ytolnuc, ynmin, t_cntl_burn, &
& nuc_number, a_name, a_nuc, z_nuc, m_ex_nuc, be_nuc
USE radial_ray_module, ONLY : xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, &
& nuc_number_c=>nuc_number,a_nuc_c, z_nuc_c, be_nuc_c, uburn_c, nse_c
USE eos_snc_x_module, ONLY : nse_e=>nse, nuc_number_e=>nuc_number, &
& a_nuc_e=>a_nuc, z_nuc_e=>z_nuc, be_nuc_e=>be_nuc, a_name_e=>a_name
USE eos_snc_y_module, ONLY : nuc_number_y=>nuc_number, a_nuc_y=>a_nuc, &
& z_nuc_y=>z_nuc, be_nuc_y=>be_nuc
USE eos_snc_z_module, ONLY : nuc_number_z=>nuc_number, a_nuc_z=>a_nuc, &
& z_nuc_z=>z_nuc, be_nuc_z=>be_nuc

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx           ! energy array dimension
INTEGER, INTENT(in)              :: nnc          ! neutrino flavor dimension
INTEGER, INTENT(in)              :: ij_ray_dim   ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim   ! number of z-zones on a processor before swapping

CHARACTER (len=5), INTENT(in), DIMENSION(nnc)                           :: c_nuc_data  ! character array of nuclei

INTEGER, INTENT(in), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim)                         :: i_nuc_data  ! integer array of edit keys

REAL(KIND=double), INTENT(in), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i             ! x-griddo index
INTEGER                          :: n             ! nuclear abundance index

!-----------------------------------------------------------------------
!
!                 \\\\\ UNPACK NUCLEAR KEYS /////
!
!-----------------------------------------------------------------------


inuc                  = i_nuc_data(1,1,1)
itnuc                 = i_nuc_data(2,1,1)
nuc_number            = i_nuc_data(3,1,1)
nuc_number_c          = i_nuc_data(3,1,1)
nuc_number_e          = i_nuc_data(3,1,1)
nuc_number_y          = i_nuc_data(3,1,1)
nuc_number_z          = i_nuc_data(3,1,1)

ttolnuc               = d_nuc_data(1,1,1)
ytolnuc               = d_nuc_data(2,1,1)
ynmin                 = d_nuc_data(3,1,1)
t_cntl_burn(1)        = d_nuc_data(4,1,1)
t_cntl_burn(2)        = d_nuc_data(5,1,1)

FORALL ( i = 1:nx )
  nse_c(i,:,:)        = i_nuc_data(10+i,:,:)
END FORALL

FORALL ( i = 1:nx-1 )
  nse_e(i+1,:,:)      = i_nuc_data(10+i,:,:)
END FORALL

FORALL ( n = 1:nnc )
  a_nuc_c  (n)        = d_nuc_data(10+n,1,1)
  z_nuc_c  (n)        = d_nuc_data(10+nnc+n,1,1)
  be_nuc_c (n)        = d_nuc_data(10+3*nnc+n,1,1)
  a_name   (n)        = c_nuc_data(n)
  a_name_e (n)        = c_nuc_data(n)
  a_nuc    (n)        = d_nuc_data(10+n,1,1)
  z_nuc    (n)        = d_nuc_data(10+nnc+n,1,1)
  m_ex_nuc (n)        = d_nuc_data(10+2*nnc+n,1,1)
  be_nuc   (n)        = d_nuc_data(10+3*nnc+n,1,1)
  a_nuc_e  (n)        = d_nuc_data(10+n,1,1)
  z_nuc_e  (n)        = d_nuc_data(10+nnc+n,1,1)
  be_nuc_e (n)        = d_nuc_data(10+3*nnc+n,1,1)
  a_nuc_y  (n)        = d_nuc_data(10+n,1,1)
  z_nuc_y  (n)        = d_nuc_data(10+nnc+n,1,1)
  be_nuc_y (n)        = d_nuc_data(10+3*nnc+n,1,1)
  a_nuc_z  (n)        = d_nuc_data(10+n,1,1)
  z_nuc_z  (n)        = d_nuc_data(10+nnc+n,1,1)
  be_nuc_z (n)        = d_nuc_data(10+3*nnc+n,1,1)
END FORALL

FORALL ( i = 1:nx )
  a_nuc_rep_c (i,:,:) = d_nuc_data(10+4*nnc+(nnc+0)*nx+i,:,:)
  z_nuc_rep_c (i,:,:) = d_nuc_data(10+4*nnc+(nnc+1)*nx+i,:,:)
  be_nuc_rep_c(i,:,:) = d_nuc_data(10+4*nnc+(nnc+2)*nx+i,:,:)
  uburn_c     (i,:,:) = d_nuc_data(10+4*nnc+(nnc+3)*nx+i,:,:)
END FORALL

FORALL ( i = 1:nx, n = 1:nnc )
  xn_c(i,n,:,:)       = d_nuc_data(10+4*nnc+(n-1)*nx+i,:,:)
END FORALL

RETURN
END SUBROUTINE unpack_nuclear_data
