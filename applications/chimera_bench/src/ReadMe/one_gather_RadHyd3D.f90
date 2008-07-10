subroutine prin

! collect all rows of data onto jcol=0, and then collect them onto mype=0
!---------------------------------------------------------------------------

use global
use zone

include 'mpif.h'

! LOCALS
INTEGER :: gathbuffer_size, mpierr, m, mj, mk, jsk, ksk
REAL(KIND=4), DIMENSION(ii,jj/pey,kk/pez,pey*pez) :: recv_buff
REAL(KIND=4), DIMENSION(ii,jj,kk) :: var

!------------------------------------------------------------------------------
! everybody sends out their chunk of a global array, zro(imax,ij_ray_dim,ik_ray_dim)

gathbuffer_size = imax * ij_ray_dim * ik_ray_dim 
call MPI_GATHER(zro, gathbuffer_size, MPI_REAL, recv_buff, gathbuffer_size, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)

  do m = 0, n_proc-1
   mj = mod(m,n_proc_y)   
   mk = m / n_proc_y     
   do k = 1, ik_ray_dim
    ksk = mk*ik_ray_dim+k
    do j = 1, ij_ray_dim
     jsk = mj*ij_ray_dim + j
     do i = 1, imax
       var(i,jsk,ksk) = recv_buff(i,j,k,m+1)
     enddo
    enddo
   enddo
  enddo

! Full array is now in var(imax,jmax,kmax) on the processor with jcol=krow=0 (mype=0)


return
end



