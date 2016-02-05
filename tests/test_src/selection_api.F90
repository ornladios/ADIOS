!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!/**************************************************************/
!/*      Test selection creation and information API           */
!/*      it does not write or read data                       */
!/**************************************************************/
program selection_api
    use adios_read_mod
    implicit none

    integer             :: i, j, ierr
    integer             :: comm
    integer             :: nerr = 0, lerr
    integer*8           :: an_int8

    ! Selections
    integer*8               :: selbb, selp, selwb, selauto

    ! Parameters to create the selections
    integer, parameter      ::  ndim_set = 2
    integer*8, dimension(2) ::  start_set = (/ 2, 3 /)
    integer*8, dimension(2) ::  count_set = (/ 5, 7 /)
    integer*8, parameter    ::  npoints_set  = 4
    integer*8, dimension(8) ::  points_set = (/ 1,2, 3,4, 5,6, 7,8 /)
    integer, parameter      ::  index_set  = 9
    character(len=25)       ::  hints_set = "please work"

    ! Variables to retrieve
    integer                 ::  type_get
    integer                 ::  ndim_get
    integer*8               ::  npoints_get
    integer                 ::  index_get
    character(len=25)       ::  hints_get
    integer*8, dimension(:), allocatable ::  start_get
    integer*8, dimension(:), allocatable ::  count_get
    integer*8, dimension(:), allocatable ::  points_get

    comm = 0

    write (*,'("Test Fortran selection API")') 

    write (*,'("  Create selections")') 
    call adios_selection_boundingbox (selbb, ndim_set, start_set, count_set)
    call adios_selection_points      (selp,  ndim_set, npoints_set, points_set)
    call adios_selection_writeblock  (selwb, index_set)
    call adios_selection_auto        (selauto, hints_set)

    ! Get information about the created selections

    !
    ! Bounding Box
    !
    write (*,'("  Get information on bounding box selection")') 
    call adios_selection_get_type (selbb, type_get)
    write (*,'("    type  = ",$)') 
    nerr = nerr + check_intval (type_get, ADIOS_SELECTION_TYPE_BOUNDINGBOX)

    ! Get and check NDIM
    call adios_selection_get_ndim (selbb, ndim_get)
    write (*,'("    ndim  = ",$)') 
    nerr = nerr + check_intval (ndim_get, ndim_set)

    ! Get and check start/count arrays
    allocate (start_get(ndim_get))
    allocate (count_get(ndim_get))
    call adios_selection_get_boundingbox (selbb, start_get, count_get)
    write (*,'("    start = ",$)') 
    nerr = nerr + check_array (ndim_get, start_get, start_set)
    write (*,'("    count = ",$)') 
    nerr = nerr + check_array (ndim_get, count_get, count_set)
    deallocate (start_get)
    deallocate (count_get)


    ! 
    ! Point list
    !
    write (*,'("  Get information on points selection")') 
    call adios_selection_get_type (selp, type_get)
    write (*,'("    type    = ",$)') 
    nerr = nerr + check_intval (type_get, ADIOS_SELECTION_TYPE_POINTS)

    ! Get and check NDIM
    call adios_selection_get_ndim (selp, ndim_get)
    write (*,'("    ndim    = ",$)') 
    nerr = nerr + check_intval (ndim_get, ndim_set)
     
    ! Get and check NPOINTS
    call adios_selection_get_npoints (selp, npoints_get)
    write (*,'("    npoints = ",$)') 
    nerr = nerr + check_intval8 (npoints_get, npoints_set)
     
    ! Get and check points arrays
    allocate (points_get(ndim_get*npoints_get))
    an_int8 = 0
    call adios_selection_get_points (selp, points_get, an_int8, npoints_get)
    write (*,'("    start = ",$)') 
    j = ndim_get*npoints_get
    nerr = nerr + check_array (j, points_get, points_set)
    deallocate (points_get)


    ! 
    ! Writeblock
    !
    write (*,'("  Get information on writeblock selection")') 
    call adios_selection_get_type (selwb, type_get)
    write (*,'("    type  = ",$)') 
    nerr = nerr + check_intval (type_get, ADIOS_SELECTION_TYPE_WRITEBLOCK)

    ! Get and check INDEX
    call adios_selection_get_index (selwb, index_get)
    write (*,'("    index = ",$)') 
    nerr = nerr + check_intval (index_get, index_set)
     

    ! 
    ! Auto
    !
    write (*,'("  Get information on auto selection")') 
    call adios_selection_get_type (selauto, type_get)
    write (*,'("    type  = ",$)') 
    nerr = nerr + check_intval (type_get, ADIOS_SELECTION_TYPE_AUTO)

    ! Get and check hints
    call adios_selection_get_hints (selauto, hints_get)
    write (*,'("    hints = [",a,"]")') hints_get
    if (hints_get == hints_set) then
        write (*,'("  OK")') 
    else
        write (*,'("  FAIL. Should be [",a,"]")') hints_set 
        nerr = nerr + 1
    endif
     

    !
    ! Clean-up
    !
    write (*,'("  Free up space of selections")') 
    call adios_selection_delete (selbb)
    call adios_selection_delete (selp)
    call adios_selection_delete (selwb)
    call adios_selection_delete (selauto)

    write (*,'("  Done. Number of errors: ",i3)') nerr 

    if (nerr > 0) then
        STOP 1
    endif


contains


integer function check_intval (a, b)
    implicit none
    integer, intent(in) ::  a, b
    write (*,'(i0,$)') a
    if (a == b) then
        write (*,'("  OK")') 
        check_intval = 0
    else
        write (*,'("  FAIL. Should be ",i1)') b 
        check_intval = 1
    endif
end function check_intval
    
integer function check_intval8 (a, b)
    implicit none
    integer*8, intent(in) ::  a, b
    write (*,'(i0,$)') a
    if (a == b) then
        write (*,'("  OK")') 
        check_intval8 = 0
    else
        write (*,'("  FAIL. Should be ",i1)') b 
        check_intval8 = 1
    endif
end function check_intval8

integer function check_array (n, a, b)
    implicit none
    integer, intent(in) ::  n
    integer*8, dimension(*), intent(in) :: a, b
    !character, dimension(*), intent(in) :: header

    integer             :: i, lerr

    lerr=0
    !write (*,'(a,"[",$)') header
    write (*,'("[",$)')
    do i=1,n
        write (*,'(i0," ",$)') a(i) 
        if (a(i) /= b(i)) then
            lerr = lerr+1
        endif
    enddo

    if (lerr == 0) then
        write (*,'("]  OK")') 
        check_array = 0
    else
        write (*,'("]  FAIL. Should be [",$)')
        do i=1,n
            write (*,'(i0," ",$)') b(i)
        enddo
        write (*,'("]")') 
        check_array = 1
    endif
end function check_array


end program selection_api
