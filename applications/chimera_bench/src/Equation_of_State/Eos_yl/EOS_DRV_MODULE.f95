!-----------------------------------------------------------------------
!    Module:       eos_drv_module
!    Author:       S. W. Bruenn
!    Date:         10/26/02
!-----------------------------------------------------------------------

MODULE eos_drv_module

USE kind_module, ONLY : single, double
SAVE

!-----------------------------------------------------------------------
!  Interpolated equation of state variables
!-----------------------------------------------------------------------
!  aesv(j,i,i_ray)  : equation of state dependent variable i at radial
!   zone j.
!  aesvd(j,i,i_ray) : derivative with respect to the density of equation
!   of state dependent variable i at radial zone j.
!  aesvt(j,i,i_ray) : derivative with respect to the temperature of
!   equation of state dependent variable i at radial zone j.
!  aesvy(j,i,i_ray) : derivative with respect to the electron fraction
!   of equation of state dependent variable i at radial zone j.
!
!     i = 1   : pressure
!     i = 2   : energy
!     i = 3   : entropy
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)        :: aesv
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)        :: aesvd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)        :: aesvt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)        :: aesvy

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!  idr(j), itr(j), and iyr(j) : integers defining the location of log(rho),
!   log(t), and ye for radial zone j on the grid, or n-neutrino bremsstrahlung,
!   i.e.,
!
!           idrb(j)/dgrid < log(rho(j)) < ( idrb(j) + 1 )/dgrid 
!
!           itrb(j)/tgrid <  log(t(j))  < ( itrb(j) + 1 )/tgrid
!
!          0.5 - iyrb(j)/ygrid < ye < 0.5 - ( iyrb(j) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Equation of state quantities
!   for radial zone j are stored at the corners of unit cube j. Rates for
!   the equation of state quantities zone j are interpolated from the rates
!   stored at the corners.
!
!  idty(j)  : the density regime (i.e., 1, 2, or 3) of radial zone j as
!   given by the above inequalities.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idr
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itr
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyr

!-----------------------------------------------------------------------
!  EOS table
!-----------------------------------------------------------------------
!  estble(i,j,id,it,iy,i_ray) : equation of state variable i of radial
!   zone j y-zone ja, z-xone ka, at the unit cube corners id, it, and iy
!   (id, it, iy = 1,2). This is the table of equation of state variables
!   from which interpolations are performed.
!
!  escnst(i,j,i_ray) : constant that is added to equation of state variable
!   i of radial zone j before taking logarithms for interpolation, and
!   then subtracted after interpolation in order to avoid taking the log
!   of a negative number.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: estble
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: escnst

!-----------------------------------------------------------------------
!  Density derivative EOS table
!-----------------------------------------------------------------------
!  estbled(i,j,id,it,iy,i_ray) : density derivative of equation of state
!   variable i of radial zone j y-zone ja, z-xone ka, at the unit cube
!   corners id, it, and iy (id, it, iy = 1,2). This is the table of
!   equation of state variables from which interpolations are performed.
!
!  escnstd(i,j,i_ray) : constant that is added to the density derivative
!   of  equation of state variable i of radial zone j before taking
!   logarithms for interpolation, and then subtracted after interpolation
!   in order to avoid taking the log of a negative number.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: estbled
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: escnstd

!-----------------------------------------------------------------------
!  Temperature derivative EOS table
!-----------------------------------------------------------------------
!  estblet(i,j,id,it,iy,i_ray) : temperature derivative of equation of
!   state variable i of radial zone j y-zone ja, z-xone ka, at the unit
!   cube corners id, it, and iy (id, it, iy = 1,2). This is the table of
!   equation of state variables from which interpolations are performed.
!
!  escnstt(i,j,i_ray) : constant that is added to the temperature
!   derivative of  equation of state variable i of radial zone j before
!   taking logarithms for interpolation, and then subtracted after
!   interpolation in order to avoid taking the log of a negative number.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: estblet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: escnstt

!-----------------------------------------------------------------------
!  Lepton fraction derivative EOS table
!-----------------------------------------------------------------------
!  estblet(i,j,id,it,iy,i_ray) : lepton fraction derivative of equation of
!   state variable i of radial zone j y-zone ja, z-xone ka, at the unit
!   cube corners id, it, and iy (id, it, iy = 1,2). This is the table of
!   equation of state variables from which interpolations are performed.
!
!  escnstt(i,j,i_ray) : constant that is added to the lepton fraction
!   derivative of  equation of state variable i of radial zone j before
!   taking logarithms for interpolation, and then subtracted after
!   interpolation in order to avoid taking the log of a negative number.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: estbley
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: escnsty

REAL(KIND=double), DIMENSION(3,4,4,4)                      :: tblt

!-----------------------------------------------------------------------
!  Lepton fraction array
!-----------------------------------------------------------------------
!  yl(j,i_ray) : lepton fraction of radial zone j at timestep m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: yl

!-----------------------------------------------------------------------
!  Density above which to couple neutrinos and matter
!-----------------------------------------------------------------------
!  rho_couple : rho > rho_couple, neutrinos and matter coupled
!             : rho < rho_couple, matter only
!-----------------------------------------------------------------------

REAL(KIND=double)                                          :: rho_couple


END module eos_drv_module
