! subroutine prototype

integer adios_open (integer*8 fd, character*(*) group_name, character*(*) name, character*(*) mode)

integer adios_group_size (integer*8 fd_p, integer*8 data_size, integer*8 total_size, integer comm)

integer adios_write (integer*8 fd_p, character*(*) name, void * var);

