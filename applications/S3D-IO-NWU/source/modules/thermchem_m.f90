!-------------------------------------------------------------------------------
! James Sutherland made a MAJOR revision to this module 5-20-03
! Changes propagated throughout several routines in the code.
!
! BUG FIX 6-19-03:  coefTbounds was being nondimensionalized, and there was a
!      bug dealing with that.  It is now stored DIMENSIONALLY (units of KELVIN)
!
! BUG FIX 6-20-03:  In specEnth and specCp, 
!           n = int((temp-temp_lobound)*invEnthInc)
!      should have been:
!           n = int((temp-temp_lobound)*invEnthInc)+1
! BUG FIX 8-28-03:  In convergence test for temperature iteration, should be
!      if( abs(deltat) < atol ) rather than if( deltat < atol )
! BUG FIX AUG-16-04 Evatt Hawkes 
!      In Newton iteration in calc_temp, cpmix was not updated 
!      after successful convergence
!-------------------------------------------------------------------------------
module thermchem_m

  implicit none
  private

  real, public, allocatable :: mixMW(:,:,:)     ! mixture molecular weight
  real, public, allocatable :: avmolwt(:,:,:)   ! avmolwt = 1.0/mixMW
  real, public, allocatable :: cpmix(:,:,:)     ! mixture specific heat at constant pressure
  real, public, allocatable :: gamma(:,:,:)     ! ratio of specific heat

  real, public, allocatable :: h_chem(:)        ! dimensionless enthalpy of formation
  real, public :: Ru                            ! nondimensional universal gas constant
  integer, public :: i_react                    ! reaction switch

!-- Public Routines
  public :: initialize_thermchem        ! initialize the module
  public :: allocate_thermchem_arrays   ! allocate/deallocate arrays
  public :: calc_temp                   ! compute temperature from internal energy
  public :: calc_press                  ! compute pressure from Ideal Gas EOS
  public :: calc_gamma                  ! compute cp/cv
  public :: calc_inv_avg_mol_wt         ! compute mixture molecular weight and its inverse
  public :: mixEnth                     ! compute mixture enthalpy
  public :: specEnth,specEnth_all       ! compute species enthalpy
  public :: calc_specEnth_allpts        ! "New"? routine by ramanan. 
  public :: gibbsEnrg_all_dimT           ! compute Gibbs energy
  public :: mixCp                       ! compute mixture heat capacity
  public :: specCp                      ! compute species heat capacity

!-- Private Variables
  ! Change by Ramanan = Ru = 8.314"51". 
  ! Extra digits in order to be consistent with cklib
  real, parameter :: univ_gascon = 8.31451 ! J/mol-K
  real, parameter :: tempmin = 250.0     ! Low  temperature limit for enthalpy & Cp tables
  real, parameter :: tempmax = 3500.0    ! High temperature limit for enthalpy & Cp tables
  integer, parameter :: npts=2001        ! number of points in tables for Cp and h

! datatype to hold linear interpolation coefficients
!   y = aa * x + bb
!   aa is slope, bb is y-intercept
  real, allocatable, dimension(:,:) :: & 
    cpCoef_aa, cpCoef_bb, enthCoef_aa, enthCoef_bb, gibbsCoef_aa, gibbsCoef_bb

  real, allocatable, dimension(:,:,:) :: coef      ! thermo coefficients [J/g*K^n] n varies.
  real, allocatable, dimension(:,:) :: coefTbounds ! valid temperature ranges for coef.
  real :: temp_lobound, temp_hibound    ! nondimensional table limits
  real :: invEnthInc, invCpInc, invGibbsInc

  logical :: thermchem_initialized = .false.

!=========================================================================================

contains

!=========================================================================================

  subroutine calc_press( pressure, rho, temp, mixMW )
    !---------------------------------------------
    ! compute the pressure from the ideal gas EOS
    !---------------------------------------------
    use param_m, only : nx, ny, nz
    implicit none
    real, intent(in), dimension(nx,ny,nz) :: rho, temp, mixMW
    real, intent(out), dimension(nx,ny,nz) :: pressure
    real :: tmp

!!$    pressure = rho*Ru*temp/mixMW
    pressure = rho*Ru*temp*avmolwt

    return
  end subroutine calc_press

!=========================================================================================

  subroutine calc_gamma( gamma, cpmix, mixMW )
    !---------------------------------------------
    ! compute gamma = cp/cv = cp/(cp-R)
    !---------------------------------------------
    use param_m, only : nx, ny, nz
    implicit none
    real, intent(in), dimension(nx,ny,nz) :: cpmix, mixMW
    real, intent(out), dimension(nx,ny,nz) :: gamma

!!$    gamma = cpmix/(cpmix-Ru/mixMW)
    gamma = cpmix/(cpmix-Ru*avmolwt)

    return
  end subroutine calc_gamma

!=========================================================================================

  subroutine calc_inv_avg_mol_wt(yspecies,avmolwt)
    !---------------------------------------------------------------------------
    ! computes the mixture (average) molecular weight
    !
    ! avmolwt  - 1/(average molecular weight of gas mixture)
    ! yspecies - Species mass fraction -> Y_i (see massfr)
    !
    ! W^(-1) = \sum_{i=1}^{n_spec} Y_{i}/W_{i}
    !---------------------------------------------------------------------------
    use param_m, only : n_spec, nx,ny,nz
    use chemkin_m, only : molwt_c
    implicit none
    real, intent(in)  :: yspecies(nx,ny,nz,n_spec)
    real, intent(out) :: avmolwt(nx,ny,nz)
    integer :: i,j,k,m

  !-- compute inverse of mixture molecular weight
    m = n_spec
    avmolwt(:,:,:) = yspecies(:,:,:,m) * molwt_c(m)
    do m = 1, n_spec-1
       avmolwt(:,:,:) = avmolwt(:,:,:) + yspecies(:,:,:,m) * molwt_c(m)
    enddo
  !-- compute mixture molecular weight
    mixMW = 1.0/avmolwt

    return
  end subroutine calc_inv_avg_mol_wt

#ifdef VECTORVERSION
!=========================================================================================
! This routine is a vectorized version for the Cray

  subroutine calc_temp( temp, e0, u, ys )
    !---------------------------------------------------------------------------
    ! calculates the temperature from the total internal energy
    ! also sets mixture heat capacity (cpmix)
    !---------------------------------------------------------------------------
    use topology_m, only : myid
    use param_m, only : nx, ny, nz, n_spec, vary_in_x, vary_in_y, vary_in_z
    use reference_m, only : t_ref
    implicit none
    real, dimension(nx,ny,nz), intent(inout) :: temp
    real, dimension(nx,ny,nz), intent(in) :: e0
    real, dimension(nx,ny,nz,3), intent(in) :: u
    real, dimension(nx,ny,nz,n_spec), intent(in) :: ys

    real, dimension(nx,ny,nz) :: tmp1, deltat, enthmix
    integer, dimension(nx,ny,nz) :: n
    integer, parameter :: icountmax = 100 ! max number of iterations for temperature solve
    integer :: i,j,k,m,icount
    real :: atol
    logical :: itersuccess = .false. 

    atol = 0.001 / t_ref
 
  !-- remove the kinetic energy from the total energy to obtain the internal energy
    if(vary_in_x==1) then
       tmp1 = u(:,:,:,1)*u(:,:,:,1)
    else
       tmp1 = 0.0
    endif
    if(vary_in_y==1) tmp1 = tmp1 + u(:,:,:,2)*u(:,:,:,2)
    if(vary_in_z==1) tmp1 = tmp1 + u(:,:,:,3)*u(:,:,:,3)

  !-- now store internal energy in tmp
    tmp1 = e0 - 0.5*tmp1

    ITERATION: do icount = 1, icountmax
      if(minval(temp) < temp_lobound) then
         write(*,*)'T is less than temp_lobound',  &
                  minloc(temp), minval(temp), temp_lobound
         stop
      end if

      if(maxval(temp) > temp_hibound) then
         write(*,*)'T is greater than temp_hibound', &
                  maxloc(temp), maxval(temp), temp_hibound
         stop
      end if

      n = int((temp-temp_lobound)*invEnthInc)+1

      cpmix = 0.0
      enthmix = 0.0

      do m=1, n_spec
      do k=1, nz
      do j=1, ny
      do i=1, nx
         cpmix(i,j,k) = cpmix(i,j,k) +  &
             ys(i,j,k,m)*(cpCoef_aa(m,n(i,j,k)) * temp(i,j,k) + cpCoef_bb(m,n(i,j,k)) )
         enthmix(i,j,k) = enthmix(i,j,k) +  &
             ys(i,j,k,m)*(enthCoef_aa(m,n(i,j,k))*temp(i,j,k) + enthCoef_bb(m,n(i,j,k)))
      end do
      end do
      end do
      end do

      deltat = ( tmp1 - (enthmix-Ru*avmolwt*temp) )/( cpmix - Ru*avmolwt)
      temp   = temp   + deltat
      if( maxval(abs(deltat)) < atol) then
        n = int((temp-temp_lobound)*invEnthInc)+1
        cpmix = 0.0
        do m=1, n_spec
        do k=1, nz
        do j=1, ny
        do i=1, nx
           cpmix(i,j,k) = cpmix(i,j,k) +  &
               ys(i,j,k,m)*(cpCoef_aa(m,n(i,j,k)) * temp(i,j,k) + cpCoef_bb(m,n(i,j,k)) )
        end do
        end do
        end do
        end do
        itersuccess = .true. 
        exit ITERATION
      end if
    end do ITERATION

    if(.not. itersuccess) then
      write(6,*)'calc_temp cannot converge after 100 iterations'
      write(6,*) 'for processor with rank =',myid
      write(6,*) maxloc(abs(deltat))
      stop  !ugly termination but that's the way it is without doing a broadcast
    end if

    return
  end subroutine calc_temp

#else
!=========================================================================================

  subroutine calc_temp( temp, e0, u, ys )
    !---------------------------------------------------------------------------
    ! calculates the temperature from the total internal energy
    ! also sets mixture heat capacity (cpmix)
    !---------------------------------------------------------------------------
    use topology_m, only : myid
    use param_m, only : nx, ny, nz, n_spec, vary_in_x, vary_in_y, vary_in_z
    use reference_m, only : t_ref
    implicit none
    real, dimension(nx,ny,nz), intent(out) :: temp
    real, dimension(nx,ny,nz), intent(in) :: e0
    real, dimension(nx,ny,nz,3), intent(in) :: u
    real, dimension(nx,ny,nz,n_spec), intent(in) :: ys

    real, dimension(nx,ny,nz) :: tmp1
    integer, parameter :: icountmax = 100 ! max number of iterations for temperature solve
    integer :: i,j,k,m,icount
    real :: atol, enthmix, r_gas, deltat
!   Added by Ramanan - 04/07/05
!   Due to the seaborg IBM executing malloc and free we are manually creating an array
!   to use for that purpose. Hoping to reduce some CPU load.
    real, dimension(n_spec) :: yspec

    atol = 0.001 / t_ref
 
  !-- remove the kinetic energy from the total energy to obtain the internal energy
    if(vary_in_x==1) then
       tmp1 = u(:,:,:,1)*u(:,:,:,1)
    else
       tmp1 = 0.0
    endif
    if(vary_in_y==1) tmp1 = tmp1 + u(:,:,:,2)*u(:,:,:,2)
    if(vary_in_z==1) tmp1 = tmp1 + u(:,:,:,3)*u(:,:,:,3)

  !-- now store internal energy in tmp
    tmp1 = e0 - 0.5*tmp1

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx

             icount = 1
             r_gas = Ru*avmolwt(i,j,k)
             yspec(:) = ys(i, j, k, :)

             ITERATION: do

              !-- compute mixture heat capacity and enthalpy for this temperature
                cpmix(i,j,k) =   mixCp( yspec, temp(i,j,k) )
                enthmix      = mixEnth( yspec, temp(i,j,k) )

              !-- calculate deltat, new temp
              !   remember tmp1 holds the internal energy
                deltat = ( tmp1(i,j,k) - (enthmix-r_gas*temp(i,j,k)) )  &
                        /( cpmix(i,j,k) - r_gas )
                temp(i,j,k) = temp(i,j,k) + deltat

              !-- check for convergence
                if( abs(deltat) < atol ) then  ! converged
!                  BUG- FIX AUG-16-04 - cpmix was not updated after successful convergence
                   cpmix(i,j,k) =   mixCp( yspec, temp(i,j,k) )
                   exit
                elseif( icount > icountmax ) then  ! maximum count violation
                   write(6,*)'calc_temp cannot converge after 100 iterations'
                   write(6,*) 'for processor with rank =',myid
                   write(6,*) 'i=',i
                   write(6,*) 'j=',j
                   write(6,*) 'k=',k
                   stop  !ugly termination but that's the way it is without doing a broadcast
                else
                   icount = icount + 1
                endif
 

             enddo ITERATION

          enddo
       enddo
    enddo

    return
  end subroutine calc_temp

#endif

!=========================================================================================

  subroutine allocate_thermchem_arrays(flag)
    !---------------------------------------------
    ! allocates memory for thermchem_m arrays
    !---------------------------------------------
    use param_m, only : n_spec, nx, ny, nz
    implicit none
    integer, intent(in) :: flag

    select case (flag)
    case(1)
       ! Changes by Ramanan - 01/27/05
       ! The 2nd dimension is changed to be 7 from 6
       ! When calculating entropy the 7th coefficient is needed
       allocate(coef(n_spec,7,2));      coef=0.0
       allocate(coefTbounds(3,n_spec)); coefTbounds=0.0
       allocate(mixMW(nx,ny,nz));       mixMW=0.0
       allocate(avmolwt(nx,ny,nz));     avmolwt=0.0
       allocate(cpmix(nx,ny,nz));       cpmix=0.0
       allocate(gamma(nx,ny,nz));       gamma=0.0
       allocate(h_chem(n_spec));        h_chem=0.0
       allocate( enthCoef_aa( n_spec,1:npts) )
       allocate( enthCoef_bb( n_spec,1:npts) )
       allocate( cpCoef_aa( n_spec,1:npts) )
       allocate( cpCoef_bb( n_spec,1:npts) )
       allocate( gibbsCoef_aa( n_spec,1:npts) )
       allocate( gibbsCoef_bb( n_spec,1:npts) )
    case(-1)
       deallocate(coef)
       deallocate(coefTbounds)
       deallocate(mixMW)
       deallocate(avmolwt)
       deallocate(cpmix)
       deallocate(gamma)
       deallocate(h_chem)
       deallocate( enthCoef_aa)
       deallocate( enthCoef_bb)
       deallocate( cpCoef_aa)
       deallocate( cpCoef_bb)
       deallocate( gibbsCoef_aa)
       deallocate( gibbsCoef_bb)
    end select
    return
  end subroutine allocate_thermchem_arrays

!=========================================================================================

  subroutine initialize_thermchem(io)
    !---------------------------------------------
    ! initializes thermchem module
    !---------------------------------------------
    use topology_m, only : myid
    use chemkin_m, only : molwt_c
    use param_m, only : n_spec
    use reference_m, only : cp_ref, t_ref, a_ref
    implicit none
    integer, intent(in) :: io
    integer :: i, j, n_poly, i_sp
    real :: t_ref_h

    if(myid.eq.0) then
       write(io,*) 'initializing thermchem module...'
       write(io,*)
    endif

    call allocate_thermchem_arrays(1)

  !-- set nondimensional universal gas constant
    Ru = univ_gascon / cp_ref

    call setPolyCoefs

  !-- multiply coefficients by MW and R so we can compute
  !   the specific enthalpy more directly.
    do i_sp = 1, n_spec
       coef(i_sp,:,:) = coef(i_sp,:,:)*univ_gascon*molwt_c(i_sp)
    enddo

  !-- set min and max temperatures in h, cp tables
    temp_lobound = tempmin/t_ref
    temp_hibound = tempmax/t_ref

  !-- Check for possible problems on out-of-range temperatures in species polynomials
    do i_sp = 1,n_spec
       if(tempmin < coefTbounds(1,i_sp)) then
          if(myid==0) then
             write(io,*)'  ***WARNING***'
             write(io,601)'Table T is lower than valid range for polynomial on species',i_sp
             write(io,602)'polynomial lower temperature limit:',coefTbounds(1,i_sp)
             write(io,602)'cp&h table lower temperature limit:',tempmin
             write(io,*)
          endif
       endif
       if(tempmax > coefTbounds(3,i_sp)) then
          if(myid==0) then
             write(io,*)'  ***WARNING***'
             write(io,601)'Table T exceeds valid range for polynomials on species',i_sp
             write(io,602)'polynomial upper temperature limit:',coefTbounds(3,i_sp)
             write(io,602)'cp&h table upper temperature limit:',tempmax
             write(io,*)
          endif
       endif
    enddo
601 format(4x,A,I3)
602 format(4x,A,2x,1f6.1)

  !-- create h, cp tables
    call createEnthTable
    call createCpTable
    call createGibbsTable

  !-- set t_ref_h and n_poly for enthalpy of formation
    t_ref_h = 298.0     ! reference temperature for enthalpy of formation

  !-- calculate enthalpy of formation in J/kg
    do i_sp = 1,n_spec
       if( t_ref_h < coefTbounds(2,i_sp) ) then
          n_poly = 1
       else
          n_poly = 2
       endif
       h_chem(i_sp) =           ( coef(i_sp,6,n_poly)                   &
            + t_ref_h * ( coef(i_sp,1,n_poly)                   &
            + t_ref_h * ( coef(i_sp,2,n_poly) / 2.              &
            + t_ref_h * ( coef(i_sp,3,n_poly) / 3.              &
            + t_ref_h * ( coef(i_sp,4,n_poly) / 4.              &
            + t_ref_h *   coef(i_sp,5,n_poly) / 5. )))))
    enddo

  !-- non-dimensionalize enthalpy of formation
    h_chem = h_chem / (a_ref*a_ref)

  !-- set flag indicating that this module has been initialized.
    thermchem_initialized = .true.

  !-- write header
    if(myid.eq.0) then
       call write_header(io,'-')
    endif

    return
  end subroutine initialize_thermchem

!=========================================================================================

  real function mixEnth( ys, temp )
    !---------------------------------------------
    ! calculates the mixture enthalpy given
    ! the mixture composition and temperature
    !---------------------------------------------
    use param_m, only : n_spec
    implicit none
    real, intent(in), dimension(n_spec) :: ys
    real, intent(in) :: temp
    integer :: m, n

    if(temp > temp_hibound .or. temp < temp_lobound ) then
       write(*,*)'T is out of range!',temp,temp_lobound,temp_hibound
       stop
    endif

    n = int((temp-temp_lobound)*invEnthInc)+1

    mixEnth = 0.0
    do m=1,n_spec
       mixEnth = mixEnth + ys(m)*(enthCoef_aa(m,n)*temp + enthCoef_bb(m,n))
    enddo

    return
  end function mixEnth

!=========================================================================================

  real function mixCp( ys, temp )
    !---------------------------------------------
    ! calculates the mixture enthalpy given
    ! the mixture composition and temperature
    !---------------------------------------------
    use param_m, only : n_spec
    implicit none
    real, intent(in), dimension(n_spec) :: ys
    real, intent(in) :: temp
    integer :: m, n

    if(temp > temp_hibound .or. temp < temp_lobound ) then
       write(*,*)'T is out of range!',temp,temp_lobound,temp_hibound
       stop
    endif

    n = int((temp-temp_lobound)*invEnthInc)+1

    mixCp = 0.0
    do m=1,n_spec
       mixCp = mixCp +  &
           ys(m)*(cpCoef_aa(m,n) * temp + cpCoef_bb(m,n) )
    enddo

    return
  end function mixCp

!=========================================================================================
! mixenth used to call this function
! the math is directly built into mixenth to cut some cpu cycles.
! any changes to this routine must also be done to mixenth and specenth_all
  real function specEnth(i_sp,temp)
    !---------------------------------------------------------------------------
    ! calculates the species enthalpies using lookup-tables
    !
    ! temp     - NON-DIMENSIONAL temperature [input]
    ! i_sp     - species index [input]
    ! specEnth - NON-DIMENSIONAL specific enthalpy per unit mass [output]
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: i_sp
    real, intent(in) :: temp
    integer :: n

    if(temp > temp_hibound .or. temp < temp_lobound ) then
       write(*,*)'T is out of range!',temp,temp_lobound,temp_hibound
       stop
    endif

    n = int((temp-temp_lobound)*invEnthInc)+1
    specEnth = enthCoef_aa(i_sp,n) * temp + enthCoef_bb(i_sp,n)

    return
  end function specEnth

!=========================================================================================
! New function by Ramanan Sankaran - 01/20/05
! This calculates enthalpy of species for all species
!=========================================================================================
  function specEnth_all(temp)
  use param_m, only : n_spec
  implicit none
  real, intent(in) :: temp
  real :: specEnth_all(n_spec)
  integer n

  if(temp > temp_hibound .or. temp < temp_lobound ) then
     write(*,*)'T is out of range!',temp,temp_lobound,temp_hibound
     stop
  endif

  n = int((temp-temp_lobound)*invEnthInc)+1
  specEnth_all(:) = enthCoef_aa(:,n) * temp + enthCoef_bb(:,n)

  return
  end function specEnth_all

!=========================================================================================
! Calculate species enthalpies for all points in the domain.
! "New"?? routine by Ramanan - 04/07/05
! Same as what is done in specenth_all
! This will help vectorize the operation on the Cray X1
! And relieve some malloc/free operations on Seaborg.
!=========================================================================================
  subroutine calc_specEnth_allpts(temp, hs)
  use param_m, only : n_spec, nx, ny, nz
  implicit none
  real, intent(in) :: temp(nx,ny,nz)
  real, intent(out) :: hs(nx,ny,nz,n_spec)
                                                                                                                                                         
  real :: specEnth_all(n_spec)
                                                                                                                                                         
  integer, dimension(nx,ny,nz) ::  n
  integer i, j, k
                                                                                                                                                         
  if(minval(temp) < temp_lobound) then
    write(*,*)'T is less than temp_lobound',  &
               minloc(temp), minval(temp), temp_lobound
    stop
  end if
                                                                                                                                                         
  if(maxval(temp) > temp_hibound) then
    write(*,*)'T is greater than temp_hibound', &
               maxloc(temp), maxval(temp), temp_hibound
    stop
  end if
                                                                                                                                                         
  n = int((temp-temp_lobound)*invEnthInc)+1
                                                                                                                                                         
  do k=1, nz; 
    do j=1, ny; 
      do i = 1, nx
        hs(i,j,k,:) = enthCoef_aa(:,n(i,j,k)) * temp(i,j,k) + enthCoef_bb(:,n(i,j,k))
      end do; 
    end do; 
  end do
              
  return
  end subroutine calc_specEnth_allpts

!=========================================================================================
! New function by Ramanan Sankaran - 01/27/05
! This calculates gibbs free energy for all species
! It takes temperature in dimensional units
!=========================================================================================
  subroutine gibbsEnrg_all_dimT(tempdim, G)
  use reference_m, only: t_ref
  use param_m, only : n_spec
  implicit none
  real, intent(in) :: tempdim
  real,intent(out) :: G(n_spec)
  integer n
  real temp

  temp = tempdim/t_ref

  if(temp > temp_hibound .or. temp < temp_lobound ) then
     write(*,*)'T is out of range!',temp,temp_lobound,temp_hibound
     stop
  endif

  n = int((temp-temp_lobound)*invGibbsInc)+1
  G(:) = gibbsCoef_aa(:,n) * temp + gibbsCoef_bb(:,n)

  return
  end subroutine gibbsEnrg_all_dimT

!----------------------------------------------------------------------
! Ramanan - 01/23/05 - Build the specCP inside the mixcp itself
! mixcp is called a lot of times since it is inside the iteration 
! loop for temperature. To cut some cpu cycles specCp is calculated
! inside mixcp. That way the `if' logic and determining `n' 
! has to be done only once.
! The same functionality has been coded again in mixcp. 
! Any changes here must be repeated there.
!=========================================================================================

  real function specCp(i_sp,temp)
    !---------------------------------------------------------------------------
    ! calculates the species heat capacities using lookup-tables
    !
    ! temp   - NON-DIMENSIONAL temperature [input]
    ! i_sp   - species index [input]
    ! specCp - NON-DIMENSIONAL heat capacity at constant pressure [output]
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: i_sp
    real, intent(in) :: temp
    integer :: n

    if(temp > temp_hibound .or. temp < temp_lobound ) then
       write(*,*)'T is out of range!',temp,temp_lobound,temp_hibound
       stop
    endif

    n = int((temp-temp_lobound)*invEnthInc)+1
    specCp = cpCoef_aa(i_sp,n) * temp + cpCoef_bb(i_sp,n)

    return
  end function specCp

!=========================================================================================

  subroutine createEnthTable
    !---------------------------------------------------------------------------
    ! creates a NON-DIMENSIONAL enthalpy look-up table for each species.
    ! for each temperature in the table, interpolation coefficients are stored
    ! so that a linear interpolation may be done to get the enthalpy.
    !---------------------------------------------------------------------------
    use param_m, only : n_spec
    use reference_m, only : a_ref, t_ref
    implicit none
    real :: inc, temp1, temp2, f1, f2
    integer :: m, n

    temp_hibound = tempmax / t_ref
    temp_lobound = tempmin / t_ref

    ! set temperature increment and its inverse (the inverse is saved)
    inc = ( temp_hibound - temp_lobound ) / real( npts-1 )
    invEnthInc = 1.0/inc

    do m = 1, n_spec
       do n = 1, npts
          temp1 = temp_lobound + real(n-1)*inc
          temp2 = temp1 + inc
          f1 = his( m, temp1*t_ref ) / (a_ref*a_ref)
          f2 = his( m, temp2*t_ref ) / (a_ref*a_ref)
          enthCoef_aa(m,n)= (f2-f1)*invEnthInc
          enthCoef_bb(m,n) = f1 - enthCoef_aa(m,n) * temp1
       enddo
    enddo
    return
  end subroutine createEnthTable

!=========================================================================================
! New routine by Ramanan Sankaran - 01/27/05
! Repeating what James did for enthalpy
!----------------------------------------------------------------------
  subroutine createGibbsTable
    !---------------------------------------------------------------------------
    ! creates a Gibbs free energy  look-up table for each species.
    ! What is stored is G/Ru, where Ru is universal gas constant
    ! and G is same as gotten from CKGML
    ! Ths unit is 1/K
    ! for each temperature in the table, interpolation coefficients are stored
    ! so that a linear interpolation may be done to get G.
    !---------------------------------------------------------------------------
    use param_m, only : n_spec
    use reference_m, only : a_ref, t_ref
    use chemkin_m, only : molwt
    implicit none
    real :: inc, temp1, temp2, f1, f2
    integer :: m, n

    temp_hibound = tempmax / t_ref
    temp_lobound = tempmin / t_ref

    ! set temperature increment and its inverse (the inverse is saved)
    inc = ( temp_hibound - temp_lobound ) / real( npts-1 )
    invGibbsInc = 1.0/inc

    do m = 1, n_spec
       do n = 1, npts
          temp1 = temp_lobound + real(n-1)*inc
          temp2 = temp1 + inc
          f1 = (his( m, temp1*t_ref ) - temp1*t_ref*entropy(m, temp1*t_ref)) &
                    *molwt(m)/univ_gascon
          f2 = (his( m, temp2*t_ref ) - temp2*t_ref*entropy(m, temp2*t_ref)) &
                    *molwt(m)/univ_gascon
          gibbsCoef_aa(m,n)= (f2-f1)*invGibbsInc
          gibbsCoef_bb(m,n) = f1 - gibbsCoef_aa(m,n) * temp1
       enddo
    enddo
    return
  end subroutine createGibbsTable

!=========================================================================================

  subroutine createCpTable
    !---------------------------------------------------------------------------
    ! creates a NON-DIMENSIONAL heat capacity look-up table for each species.
    ! for each temperature in the table, interpolation coefficients are stored
    ! so that a linear interpolation may be done to get the heat capacity.
    !---------------------------------------------------------------------------
    use param_m, only : n_spec
    use reference_m, only : cp_ref, t_ref
    implicit none
    real :: inc, temp1, temp2, f1, f2
    integer :: m,n

    temp_hibound = tempmax / t_ref
    temp_lobound = tempmin / t_ref

    ! set temperature increment and its inverse (the inverse is saved)
    inc = ( temp_hibound - temp_lobound ) / real ( npts-1 )
    invCpInc = 1.0/inc

    do m = 1, n_spec
       do n = 1, npts
          temp1 = temp_lobound + real(n-1)*inc
          temp2 = temp1 + inc
          f1 = cps( m, temp1*t_ref ) / cp_ref
          f2 = cps( m, temp2*t_ref ) / cp_ref
          cpCoef_aa(m,n) = (f2-f1)*invCpInc
          cpCoef_bb(m,n) = f1 - cpCoef_aa(m,n) * temp1
       enddo
    enddo
    return
  end subroutine createCpTable

!=========================================================================================

  real function cps(m,temp)
    !---------------------------------------------------------------------------
    ! Computes the specific heat at constant pressure.
    ! Since the specific heat has only temperature dependence for dilute gases,
    ! (we think), it is represented as a polynomial in temperature as:
    !
    ! c_{p} = c1 + c2*T + c3*T^2 + c4*T^3 + c5*T^4
    !
    ! also, for more accuracy,there are low and high temperature coefficients
    ! for this polynomial.
    !
    ! cps   - DIMENSIONAL specific heat at constant pressure per unit mass [J/(kg*K)]
    ! coef  - DIMENSIONAL coefficients of high or low temperature polynomial fit
    ! temp  - DIMENSIONAL temperature [K]
    ! m     - mth species
    ! nt    - low or high temperature curve fit coefficients
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m
    real,    intent(in) :: temp         ! in DEGREES KELVIN (dimensional)

    integer :: nt

  !-- set temperature range nt
    if( temp < coefTbounds(2,m) ) then
       nt = 1
    else
       nt = 2
    endif

  !-- calculate specific heat (note the compounding of temperature exponents!)
    cps =           coef(m,1,nt)                &
         + temp * ( coef(m,2,nt)                &
         + temp * ( coef(m,3,nt)                &
         + temp * ( coef(m,4,nt)                &
         + temp *   coef(m,5,nt) )))

    return
  end function cps

!=========================================================================================

  real function his(m,temp)
    !---------------------------------------------------------------------------
    ! calculates the specific enthalpy for species i
    !
    ! Notes by Ramanan - 01/27/05
    ! The equation below is wrong. But the coding is right.
    !   h = c1*T + c2*T^2 + c3*T^3 + c4*T^4 + c5*T^5 + c6
    !
    ! his    - DIMENSIONAL species specific enthalpy [J/kg]
    ! coef   - DIMENSIONAL polynomial (in temperature) coefficients for curve fit
    ! temp   - DIMENSIONAL Temperature [K]
    ! m     - mth species
    ! nt    - low or high temperature curve fit coefficients
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m
    real,    intent(in) :: temp         ! in DEGREES KELVIN (dimensional)
    integer :: nt

  !-- set temperature range nt
    if( temp < coefTbounds(2,m) ) then
       nt = 1
    else
       nt = 2
    endif

  !-- calculate specific enthalpy (note the compounding of temperature exponents!)
    his =           coef(m,6,nt)                     &
         + temp * ( coef(m,1,nt)                     &
         + temp * ( coef(m,2,nt) / 2.0               &
         + temp * ( coef(m,3,nt) / 3.0               &
         + temp * ( coef(m,4,nt) / 4.0               &
         + temp *   coef(m,5,nt) / 5.0 ))))

    return
  end function his

!=========================================================================================

  real function entropy(m,temp)
    !---------------------------------------------------------------------------
    ! calculates the entropy  Sk for species i
    !
    !   Sk = c1*ln T + c2*T + c3*T^2/2 + c4*T^3/3 + c5*T^4/4 + c7
    !
    ! his    - DIMENSIONAL species specific enthalpy [J/kg]
    ! coef   - DIMENSIONAL polynomial (in temperature) coefficients for curve fit
    ! temp   - DIMENSIONAL Temperature [K]
    ! m     - mth species
    ! nt    - low or high temperature curve fit coefficients
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m
    real,    intent(in) :: temp         ! in DEGREES KELVIN (dimensional)
    integer :: nt

  !-- set temperature range nt
    if( temp < coefTbounds(2,m) ) then
       nt = 1
    else
       nt = 2
    endif

  !-- calculate specific enthalpy (note the compounding of temperature exponents!)
    entropy =       coef(m,1,nt)*log(temp)           &
         +          coef(m,7,nt)                              &
         + temp * ( coef(m,2,nt)                     &
         + temp * ( coef(m,3,nt) / 2.0               &
         + temp * ( coef(m,4,nt) / 3.0               &
         + temp *   coef(m,5,nt) / 4.0 )))

    return
  end function entropy

!=========================================================================================

  subroutine setPolyCoefs
    !---------------------------------------------------------------------------
    ! Sets the polynomial coefficients for heat capacity evaulation for each
    ! species using Chemkin's information.  This requires that Chemkin is
    ! initialized (chemkin_m)
    !---------------------------------------------------------------------------
    use param_m, only : n_spec
    use chemkin_m, only : ickwrk, rckwrk
    implicit none
    real, dimension(7,2,n_spec) :: tmpCoef1
    integer :: maxtp, icoef, ispec, itrange
    integer, dimension(n_spec) :: nt

    tmpCoef1 = 0.0

    call CKMXTP( ickwrk, maxtp )
    call CKATHM( 7, 2, ickwrk, rckwrk, maxtp, nt, coefTbounds, tmpCoef1 )

    do itrange=1,2
       do icoef=1,7
          do ispec=1,n_spec
             coef(ispec,icoef,itrange) = tmpCoef1(icoef,itrange,ispec)
          enddo
       enddo
    enddo

    return
  end subroutine setPolyCoefs

!=========================================================================================

end module thermchem_m
