!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

! 
! Query Fortran 90 API for ADIOS 
!    
! Use this module in your source code to ensure that
! you are calling the adios_* query functions with
! the correct arguments
!
module adios_query_mod

    use adios_read_mod

    type ADIOS_QUERY
        private
        integer*8  q
    end type ADIOS_QUERY

    !
    ! ADIOS Query method                                 
    !
    integer, parameter :: ADIOS_QUERY_METHOD_FASTBIT  = 1 
    integer, parameter :: ADIOS_QUERY_METHOD_ALACRITY = 2 

    !
    ! Predicate
    !
    integer, parameter :: ADIOS_LT           = 0
    integer, parameter :: ADIOS_LTEQ         = 1
    integer, parameter :: ADIOS_GT           = 2
    integer, parameter :: ADIOS_GTEQ         = 3
    integer, parameter :: ADIOS_EQ           = 4
    integer, parameter :: ADIOS_NE           = 5

    ! 
    ! Clause
    !
    integer, parameter :: ADIOS_QUERY_OP_AND = 0
    integer, parameter :: ADIOS_QUERY_OP_OR  = 1


interface

    subroutine adios_query_create (f, sel, varname, pred, value, q)
        import :: ADIOS_QUERY
        implicit none
        integer*8,         intent(in)  :: f       ! ADIOS FILE (from adios_read_open())
        integer*8,         intent(in)  :: sel     ! ADIOS_SELECTION from read API
        character(*),      intent(in)  :: varname
        integer,           intent(in)  :: pred    ! PREDICATE like ADIOS_GT
        character(*),      intent(in)  :: value   ! comparison value (integer or real)
        type(ADIOS_QUERY), intent(out) :: q       ! output variable, 0 on error
    end subroutine

    subroutine adios_query_combine (q1, op, q2, q)
        import :: ADIOS_QUERY
        implicit none
        type(ADIOS_QUERY), intent(in)  :: q1   ! Query 1
        integer,           intent(in)  :: op   ! Clasue like ADIOS_QUERY_OP_AND
        type(ADIOS_QUERY), intent(in)  :: q2   ! Query 2
        type(ADIOS_QUERY), intent(out) :: q    ! Result Query 
    end subroutine

    subroutine adios_query_set_method (q, method)
        import :: ADIOS_QUERY
        implicit none
        type(ADIOS_QUERY), intent(in)  :: q        ! Query 
        integer,           intent(in)  :: method   ! Method like ADIOS_QUERY_METHOD_FASTBIT
    end subroutine

    integer*8 function adios_query_estimate (q, timestep)
        ! return the (estimated) number of points (an upper bound)
        import :: ADIOS_QUERY
        implicit none
        type(ADIOS_QUERY), intent(in)  :: q        ! Query 
        integer,           intent(in)  :: timestep ! must be 0 in case of streaming
    end function

    subroutine adios_query_evaluate (q, sel_outboundary, timestep, batchsize, sel_result, err)
        import :: ADIOS_QUERY
        implicit none
        type(ADIOS_QUERY), intent(in)  :: q           ! Query 
        integer*8,         intent(in)  :: sel_outboundary  ! apply hits on this selection
        integer,           intent(in)  :: timestep    ! must be 0 in case of streaming
        integer*8,         intent(in)  :: batchsize   ! limit result size for one call
        integer*8,         intent(out) :: sel_result  ! result selection (ADIOS point selection)
        integer,           intent(out) :: err         ! 0 on OK
    end subroutine

    subroutine adios_query_free (q)
        import :: ADIOS_QUERY
        implicit none
        type(ADIOS_QUERY), intent(in)  :: q
    end subroutine

end interface


!
! ADIOS_QUERY_CREATE generic interface
!
! Usage: call adios_query_create (fp, varname, sel, pred, value, query, err)
!        integer*8,      intent(in)  :: fp     ! ADIOS FILE pointer
!        character(*),   intent(in)  :: varname
!        integer*8,      intent(in)  :: sel    ! ADIOS_SELECTION from read API
!        integer,        intent(in)  :: pred   ! PREDICATE like ADIOS_GT
!                        intent(in)  :: value  ! comparison value (integer or real)
!        integer*8,      intent(out) :: query  ! output variable, 0 on error
!
!interface adios_query_create
!
!    subroutine adios_query_create_int1 (fp, varname, sel, pred, value, query)
!        implicit none
!        integer*8,      intent(in)  :: fp
!        character(*),   intent(in)  :: varname
!        integer*8,      intent(in)  :: sel
!        integer,        intent(in)  :: pred 
!        integer*1,      intent(in)  :: value
!        integer*8,      intent(out) :: query
!    end subroutine
!
!    subroutine adios_query_create_int2 (fp, varname, sel, pred, value, query)
!        implicit none
!        integer*8,      intent(in)  :: fp
!        character(*),   intent(in)  :: varname
!        integer*8,      intent(in)  :: sel
!        integer,        intent(in)  :: pred 
!        integer*2,      intent(in)  :: value
!        integer*8,      intent(out) :: query
!    end subroutine
!
!    subroutine adios_query_create_int4 (fp, varname, sel, pred, value, query)
!        implicit none
!        integer*8,      intent(in)  :: fp
!        character(*),   intent(in)  :: varname
!        integer*8,      intent(in)  :: sel
!        integer,        intent(in)  :: pred 
!        integer*4,      intent(in)  :: value
!        integer*8,      intent(out) :: query
!    end subroutine
!
!    subroutine adios_query_create_int8 (fp, varname, sel, pred, value, query)
!        implicit none
!        integer*8,      intent(in)  :: fp
!        character(*),   intent(in)  :: varname
!        integer*8,      intent(in)  :: sel
!        integer,        intent(in)  :: pred 
!        integer*8,      intent(in)  :: value
!        integer*8,      intent(out) :: query
!    end subroutine
!
!    subroutine adios_query_create_real4 (fp, varname, sel, pred, value, query)
!        implicit none
!        integer*8,      intent(in)  :: fp
!        character(*),   intent(in)  :: varname
!        integer*8,      intent(in)  :: sel
!        integer,        intent(in)  :: pred 
!        real*4,         intent(in)  :: value
!        integer*8,      intent(out) :: query
!    end subroutine
!
!    subroutine adios_query_create_real8 (fp, varname, sel, pred, value, query)
!        implicit none
!        integer*8,      intent(in)  :: fp
!        character(*),   intent(in)  :: varname
!        integer*8,      intent(in)  :: sel
!        integer,        intent(in)  :: pred 
!        real*8,         intent(in)  :: value
!        integer*8,      intent(out) :: query
!    end subroutine
!
!end interface adios_query_create


contains

    logical function adios_query_is_method_available (method)
        ! returns .true. if method is available
        implicit none
        integer,           intent(in)  :: method   ! Method like ADIOS_QUERY_METHOD_FASTBIT
        integer :: ret
        logical :: avail

        call adios_query_is_method_available_f2c (method, ret)
        if (ret.ne.0) then
            avail = .true.
        else 
            avail = .false.
        endif
        adios_query_is_method_available = avail
    end function

end module adios_query_mod
