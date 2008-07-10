MODULE controls
!===============================================================================
!  This module contains the values of the flags and limits which control 
!  the behavior of the network
!===============================================================================
USE kind_module, ONLY : double

INTEGER            :: iweak     ! If =0, weak reaction are ignored
                                ! If <0, only weak reactions are used
INTEGER            :: iscrn     ! If =0, screening is ignored
INTEGER            :: itso      ! Sets level of time series output
INTEGER            :: idiag     ! Sets level of diagnostic output
INTEGER            :: iconvc    ! Controls type of convergence condition (0=mass)
INTEGER            :: kstmx     ! Max # of timesteps before exit
INTEGER            :: knrmx     ! Max # of Newton-Raphson iterations
INTEGER            :: nzone     ! Number of zones

REAL(KIND=double)  :: tolc      ! The iterative convergence test limit
REAL(KIND=double)  :: tolm      ! Max network mass error
REAL(KIND=double)  :: changemx  ! Relative abundance change used to guess timestep
REAL(KIND=double)  :: ytime     ! Min abundance included in timestep estimate
REAL(KIND=double)  :: ymin      ! Min abundance, y < ymin =0
REAL(KIND=double)  :: tdelmm    ! new timestep <= tdelmm * previous timestep

      common /flag/ iweak,iscrn,itso,idiag,iconvc,kstmx,knrmx,nzone
      common /tol/ changemx,tolm,tolc,ytime,ymin,tdelmm

END MODULE controls
