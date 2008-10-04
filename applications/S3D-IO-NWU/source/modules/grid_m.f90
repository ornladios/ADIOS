!=========================================================================================
  module grid_m
!=========================================================================================
! Change Record
!-----------------------------------------------------------------------------------------
! Evatt Hawkes NOV 2004
! Added target attribute to scale_1 arrays for in plane derivative routines.
!

  implicit none
!-----------------------------------------------------------------------------------------
! integers

  integer unif_grid_x   !uniform grid switch for x-direction
  integer unif_grid_y   !uniform grid switch for y-direction
  integer unif_grid_z   !uniform grid switch for z-direction

  integer x_center      !switch for grid compression location in x-direction
  integer y_center      !switch for grid compression location in y-direction
  integer z_center      !switch for grid compression location in z-direction

! reals

  real xmin             !minimum location of grid in x-direction
  real xmax             !maximum location of grid in x-direction
  real ymin             !minimum location of grid in y-direction
  real ymax             !maximum location of grid in y-direction
  real zmin             !minimum location of grid in z-direction
  real zmax             !maximum location of grid in z-direction

  real bx               !grid compression parameter for x-direction
  real by               !grid compression parameter for y-direction
  real bz               !grid compression parameter for z-direction

  real delx             !uniform grid spacing in x-direction
  real dely             !uniform grid spacing in y-direction
  real delz             !uniform grid spacing in z-direction

! real arrays

  real, allocatable :: x(:)   !coordinates in x-direction
  real, allocatable :: y(:)   !coordinates in y-direction
  real, allocatable :: z(:)   !coordinates in z-direction

  real, allocatable :: ad2fx(:)
  real, allocatable :: ad2fy(:)
  real, allocatable :: ad2fz(:)
  real, allocatable :: adfx(:)
  real, allocatable :: adfy(:)
  real, allocatable :: adfz(:)

! 24-NOV-2004: Evatt Hawkes
! adding target attribute for in-plane derivative operations
! needed in subroutines deriv_inplane_1 and deriv_inplane_2
!  real, allocatable :: scale_1x(:)
!  real, allocatable :: scale_1y(:)
!  real, allocatable :: scale_1z(:) 
  real, target, allocatable :: scale_1x(:)
  real, target, allocatable :: scale_1y(:)
  real, target, allocatable :: scale_1z(:)

  real, allocatable :: scale_2x(:)
  real, allocatable :: scale_2y(:)
  real, allocatable :: scale_2z(:)
  real, allocatable :: scale_3x(:)
  real, allocatable :: scale_3y(:)
  real, allocatable :: scale_3z(:)
!-----------------------------------------------------------------------------------------
  contains
!=========================================================================================
  subroutine allocate_grid_arrays(flag)
!=========================================================================================
! allocate grid arrays
!-----------------------------------------------------------------------------------------
  use param_m, only : nx, ny, nz

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed

  integer flag
!-----------------------------------------------------------------------------------------
! grid arrays

  if(flag.eq.1) then

    allocate(x(nx));    x=0.0
    allocate(y(ny));    y=0.0
    allocate(z(nz));    z=0.0

    allocate(ad2fx(nx));    ad2fx=0.0
    allocate(ad2fy(ny));    ad2fy=0.0
    allocate(ad2fz(nz));    ad2fz=0.0
    allocate(adfx(nx));     adfx=0.0
    allocate(adfy(ny));     adfy=0.0
    allocate(adfz(nz));     adfz=0.0
    allocate(scale_1x(nx));     scale_1x=0.0
    allocate(scale_1y(ny));     scale_1y=0.0
    allocate(scale_1z(nz));     scale_1z=0.0
    allocate(scale_2x(nx));     scale_2x=0.0
    allocate(scale_2y(ny));     scale_2y=0.0
    allocate(scale_2z(nz));     scale_2z=0.0
    allocate(scale_3x(nx));     scale_3x=0.0
    allocate(scale_3y(ny));     scale_3y=0.0
    allocate(scale_3z(nz));     scale_3z=0.0

  elseif(flag.eq.-1) then

    deallocate(x)
    deallocate(y)
    deallocate(z)

    deallocate(ad2fx)
    deallocate(ad2fy)
    deallocate(ad2fz)
    deallocate(adfx)
    deallocate(adfy)
    deallocate(adfz)
    deallocate(scale_1x)
    deallocate(scale_1y)
    deallocate(scale_1z)
    deallocate(scale_2x)
    deallocate(scale_2y)
    deallocate(scale_2z)
    deallocate(scale_3x)
    deallocate(scale_3y)
    deallocate(scale_3z)

  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine allocate_grid_arrays
!=========================================================================================
  subroutine initialize_grid(io)
!=========================================================================================
! generates the grid locations and the grid jacobians for the derivative evaluation
!
! b(x,y,z) - Grid compression parameters in (x,y,z) directions
! fmap     - maps uniform points in s to nonuniform points in
!    physical space
! dfmap    - computes (d(x,y,z)/ds)/x_range as a function of s
! d2fmap   - computes (d^{2}(x,y,z)/ds^2)/x_range as a function of s
! nx,ny,nz - number of grid points in (x,y,z)-direction
! x/y/z_center - -1 ==> compressed points close to left edge
!     0 ==> compressed points close to the middle
!     1 ==> compressed points close to the right edge
!
! note: the programmed transformation is given by
!
! X = (X_{max} - X_{min})*(Sinh(b*s)/Sinh(b)) + X_{min}
!   = X_{range}*(Sinh(b*s)/Sinh(b)) + X_{min}
!
! Sinh(0) = 0, so we can either have maximum grid compression
! on the left or right end of the boundary or in the middle.
! This mapping has DISCONTINUOUS derivatives at the boundaries
! when the grid is compressed.
!
! For x:
!
! ad2fx  - d^2x/ds^2
! adfx   - dx/ds
! x  - Coordinates of possibly nonuniform grid
! xmax   - Largest value of x on grid
! xmin   - Smallest value of x on grid
! unif_grid_=() = 0 - Grid is non-uniform in the ()-direction
!  d()/dx = [ d(s)/dx ]*[ d()/d(s) ]
!   [ dx/d(s) ]^-1*[ d()/d(s) ]
!
!  d^2()/dx^2 = [ dx/d(s) ]^-2*[ d^2()/d(s)^2 ] -
!   [ dx/d(s) ]^-3[ d^2x/d(s)^2 ]*[ d()/d(s) ]
!
! scale_1 - Proportional to ( dx/d(eta) )^-1 where eta is fictitious
!   uniform grid and x is nonuniform grid. When the grid is
!   uniform then scale_1 is constant along x-axis.
! scale_2 - Proportional to ( d^2x/d(eta)^2 )*( dx/d(eta) )^-3.
!    When the grid is uniform then scale_2 is zero.
! scale_3 - Proportional to ( dx/d(eta) )^-2. When the grid is
!   uniform then scale_3 is constant along x-axis.
!
! s    - position along fictitious uniform grid -
!    0 <= s <= 1, -1 <= s <=0,  or -1 <= s <= 1
! smin     - either -1 or 0
! smax     - either  0 or 1

!----------------------------------------------------------------------------------------
  use param_m
  use topology_m
  use reference_m

  implicit none
!----------------------------------------------------------------------------------------
! declarations passed in

  integer io

! local declarations

  real sxmin, symin, szmin
  real sxmax, symax, szmax
  real s, scale
  real  x_range,  y_range,  z_range
  real sx_range, sy_range, sz_range
  real dsx ,dsy ,dsz
  real sxmin_l, symin_l , szmin_l

  integer i
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    write(io,*) 'initializing grid module...'
    write(io,*)
  endif
!----------------------------------------------------------------------------------------
! allocate arrays

  call allocate_grid_arrays(1)
!----------------------------------------------------------------------------------------
! non-dimensionalize domain size
! factor of 100.0 converts input in cm to m

  xmin=(xmin/100.0)/l_ref
  xmax=(xmax/100.0)/l_ref

  ymin=(ymin/100.0)/l_ref
  ymax=(ymax/100.0)/l_ref

  zmin=(zmin/100.0)/l_ref
  zmax=(zmax/100.0)/l_ref

! if direction is not active, make min=0 and max=1.0

  if(vary_in_x.ne.1) then
    xmin=0.0; xmax=1.0
  endif

  if(vary_in_y.ne.1) then
    ymin=0.0; ymax=1.0
  endif

  if(vary_in_z.ne.1) then
    zmin=0.0; zmax=1.0
  endif
!----------------------------------------------------------------------------------------
! hard code some stuff for now

  bx=0.0
  by=0.0
  bz=0.0

  unif_grid_x=1   !always uniform grid
  unif_grid_y=1   !always uniform grid
  unif_grid_z=1   !always uniform grid
!----------------------------------------------------------------------------------------
! check

  if((unif_grid_x.ne.1 .and. bx.eq. 0.) .or.            &
     (unif_grid_y.ne.1 .and. by.eq. 0.) .or.            &
     (unif_grid_z.ne.1 .and. bz.eq. 0.)    ) then
    write(io,*) 'Beta=0 is not allowed for non-uniform grids '
    call terminate_run(io,0)  !must be called by all processors
  endif
!----------------------------------------------------------------------------------------
! hard code these for now (sdm)

  x_center=0
  y_center=0
  z_center=0
!----------------------------------------------------------------------------------------
  select case (x_center)
  case (-1)            ! compress along left edge
     sxmin = 0.;    sxmax = 1.
  case (0)             ! compress in the middle
     sxmin = -1.;    sxmax = 1.
  case (1)             ! compress along right edge
     sxmin = -1.;    sxmax = 0.
  end select

  select case (y_center)
  case (-1)            ! compress along left edge
     symin = 0.;    symax = 1.
  case (0)             ! compress in the middle
     symin = -1.;    symax = 1.
  case (1)             ! compress along right edge
     symin = -1.;    symax = 0.
  end select

  select case (z_center)
  case (-1)            ! compress along left edge
     szmin = 0.;    szmax = 1.
  case (0)             ! compress in the middle
     szmin = -1.;    szmax = 1.
  case (1)             ! compress along right edge
     szmin = -1.;    szmax = 0.
  end select
!----------------------------------------------------------------------------------------
  x_range = ( xmax - xmin )
  y_range = ( ymax - ymin )
  z_range = ( zmax - zmin )

  delx=x_range  !default value if x-direction is not active
  dely=y_range  !default value if y-direction is not active
  delz=z_range  !default value if z-direction is not active

  if(unif_grid_x.eq.1) then
     sxmin = 0
     sxmax = 1
     if(vary_in_x.eq.1) delx=x_range/real(nx_g-1)
  endif

  if(unif_grid_y.eq.1) then
     symin = 0
     symax = 1
     if(vary_in_y.eq.1) dely=y_range/real(ny_g-1)
  endif

  if(unif_grid_z.eq.1) then
     szmin = 0
     szmax = 1
     if(vary_in_z.eq.1) delz=z_range/real(nz_g-1)
  endif

  sx_range= sxmax - sxmin;
  sy_range= symax - symin;
  sz_range= szmax - szmin;

  dsx     = (sx_range)/(nxm_g)
  dsy     = (sy_range)/(nym_g)
  dsz     = (sz_range)/(nzm_g)

  sxmin_l = sxmin + xid*dsx*nx
  symin_l = symin + yid*dsy*ny
  szmin_l = szmin + zid*dsz*nz
!----------------------------------------------------------------------------------------
! X-direction
!----------------------------------------------------------------------------------------
  if(nx_g .gt. 1 ) then
    do i = 1, nx
      s = sxmin_l + dsx*real(i-1)
      x(i) = (x_range/sx_range) * ( fmap (1, unif_grid_x, s, bx ) - sxmin) +  xmin
      adfx(i) = (x_range/sx_range) * dfmap (1, unif_grid_x, s, bx )
      ad2fx(i) = (x_range/sx_range) * d2fmap (1, unif_grid_x, s, bx )
      scale_1x(i) = 1.0 / adfx(i)
      scale_3x(i) = scale_1x(i) * scale_1x(i)
      scale_2x(i) = ad2fx(i) * scale_3x(i) * scale_1x(i)
    enddo
  else
    x(1) = xmin
  endif
!----------------------------------------------------------------------------------------
! Y-direction
!----------------------------------------------------------------------------------------
  if(ny_g .gt. 1 ) then
    do i = 1, ny
      s = symin_l + dsy*real(i-1)
      y(i) = (y_range/sy_range) * ( fmap (1, unif_grid_y, s, by ) - symin) +  ymin
      adfy(i) = (y_range/sy_range) * dfmap (1, unif_grid_y, s, by )
      ad2fy(i) = (y_range/sy_range) * d2fmap (1, unif_grid_y, s, by )
      scale_1y(i) = 1.0 / adfy(i)
      scale_3y(i) = scale_1y(i) * scale_1y(i)
      scale_2y(i) = ad2fy(i) * scale_3y(i) * scale_1y(i)
    enddo
  else
    y(1) = ymin
  endif
!----------------------------------------------------------------------------------------
! Z-direction
!----------------------------------------------------------------------------------------
  if(nz_g .gt. 1 ) then
    do i = 1, nz
      s = szmin_l + dsz*real(i-1)
      z(i) = (z_range/sz_range) * ( fmap (1, unif_grid_z, s, bz ) - szmin) +  zmin
      adfz(i) = (z_range/sz_range) * dfmap (1, unif_grid_z, s, bz )
      ad2fz(i) = (z_range/sz_range) * d2fmap (1, unif_grid_z, s, bz )
      scale_1z(i) = 1.0 / adfz(i)
      scale_3z(i) = scale_1z(i) * scale_1z(i)
      scale_2z(i) = ad2fz(i) * scale_3z(i) * scale_1z(i)
    enddo
  else
    z(1) = zmin
  endif

!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif
!----------------------------------------------------------------------------------------
  return
  end subroutine initialize_grid
!=========================================================================================
  real function fmap(imap,unif_grid,s,beta)
!=========================================================================================
! defines the function mapping from a uniform grid running from
! (0 to 1), (-1 to 0), or (-1 to 1) to a (non)uniform grid in the
! also in the ideal domain.
!
! This subroutine is intimately related to DEFGRI.f (DEFine GRId).
!
! Important Note:
!
! For this mapping, beta cannot be zero for non-uniform grids; so
! error trapping for (unif_grid.ne.1 .and. Beta .eq. 0.) should
! be done in the calling program to avoid overhead. So the following
! code should appear in the caller:
!
! if (unif_grid_<x/y/z>.ne.1 .and. Beta_<x/y/z> .eq. 0.) then
!   write(6,*) 'Beta=0 is not allowed for non-uniform grids '
!   call terminate_run(6,0)
! endif
!
! imap    - Possible compression schemes
!        (currently unused - only one scheme is implemented)
! unif_grid   - Is the grid uniform (0 = no)
! s       - location on fictitious uniform grid
! beta    - Grid compression parameter
!-----------------------------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------------------------
  integer imap
  integer unif_grid
  real    beta
  real    s
!-----------------------------------------------------------------------------------------
! perform map
!-----------------------------------------------------------------------------------------
  if(unif_grid.eq.1) then
    fmap = s
  else
    fmap=sinh(beta*s)/sinh(beta)
  endif
!-----------------------------------------------------------------------------------------
  return
  end function fmap
!=========================================================================================
  real function dfmap(imap,unif_grid,s,beta)
!=========================================================================================
! defines d(x)/ds where s is a unifrom grid defined from
! -1 to 1, -1 to 0, or 0 to 1 and x is a possibly nonuniform grid
! defined from xmin to xmax. X is the physical grid.
! This subroutine is intimately related to DEFGRI.f (DEFine GRId).
! notice that this function call may be made in either direction but
! the variable names would have you think that it's only good in the
! x-direction.
!
! please see comments in fmap
!-----------------------------------------------------------------------------------------
  implicit none

  integer imap
  integer unif_grid
  real    s
  real    beta
!-----------------------------------------------------------------------------------------
! perform map
!-----------------------------------------------------------------------------------------
  if(Unif_grid.eq.1) then
     dfmap = 1.0
  else
     dfmap=beta*cosh(beta*s)/sinh(beta)
  endif
!-----------------------------------------------------------------------------------------
  return
  end function dfmap
!=========================================================================================
  real function d2fmap(imap,unif_grid,s,beta)
!=========================================================================================
! defines d^2(x)/ds^2 where s is a unifrom grid defined from
! -1 to 1, -1 to 0, or 0 to 1 and x is a possibly nonuniform grid
! defined from xmin to xmax; x is the physical grid
! this subroutine is intimately related to DEFGRI.f (DEFine GRId).
!
! please see comments in fmap
!-----------------------------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------------------------
! declarations

  integer imap
  integer unif_grid
  real    beta
  real    s
!-----------------------------------------------------------------------------------------
! perform map
!-----------------------------------------------------------------------------------------
  if(unif_grid.eq.1) then
    d2fmap = 0.0
  else
    d2fmap=beta*beta*sinh(beta*s)/sinh(beta)
  endif
!-----------------------------------------------------------------------------------------
  return
  end function d2fmap
!-----------------------------------------------------------------------------------------
  end module grid_m
