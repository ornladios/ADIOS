28c28,30
<   #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
---
>   integer*8 :: group_prof_handle
>   integer :: adios_err
>   #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b,adios_err)
379,382c381,393
<      call adios_get_group (group_handle, "snapshot"//char(0))
<      call adios_open (file_handle, group_handle, bpdum)
<      call adios_write (file_handle, "timestep"//char(0), (mstepall+istep))
<      call adios_write (file_handle, "time"//char(0), q0+0.5*q1+0.25*q2)
---
> 
> !!! modified by zf2 for timing
>      CALL open_start_for_group(group_prof_handle, "snapshot"//char(0),istep)
> 
>      !!!call adios_get_group (group_handle, "snapshot"//char(0))
>      call adios_open (file_handle, "snapshot"//char(0), bpdum, "w"//char(0), adios_err)
> 
>      CALL open_end_for_group(group_prof_handle,istep)
> 
>      CALL write_start_for_group(group_prof_handle,istep)
> 
>      call adios_write (file_handle, "timestep"//char(0), (mstepall+istep),adios_err)
>      call adios_write (file_handle, "time"//char(0), q0+0.5*q1+0.25*q2,adios_err)
384,385c395,396
<      call adios_write (file_handle, "3+mbin_psi*6"//char(0), 3+mbin_psi*6)
<      call adios_write (file_handle, "jm+1"//char(0), jm+1)
---
>      call adios_write (file_handle, "3+mbin_psi*6"//char(0), 3+mbin_psi*6,adios_err)
>      call adios_write (file_handle, "jm+1"//char(0), jm+1, adios_err)
391c402
<      call adios_write(file_handle,"num_of_quantities"//char(0),i)
---
>      call adios_write(file_handle,"num_of_quantities"//char(0),i, adios_err)
399c410
<      call adios_write (file_handle, "modeeigen"//char(0), nmode)
---
>      call adios_write (file_handle, "modeeigen"//char(0), nmode, adios_err)
402,405c413,424
<      call adios_write (file_handle, "dataout_p"//char(0), booz_out)
<      call adios_write (file_handle, "dataout_f"//char(0), flux3darray)
<      call adios_write (file_handle, "dataeigen"//char(0), eigenmode)
<      call adios_close(file_handle)
---
>      call adios_write (file_handle, "dataout_p"//char(0), booz_out, adios_err)
>      call adios_write (file_handle, "dataout_f"//char(0), flux3darray, adios_err)
>      call adios_write (file_handle, "dataeigen"//char(0), eigenmode, adios_err)
> 
>      CALL write_end_for_group(group_prof_handle,istep)
> 
>      CALL close_start_for_group(group_prof_handle,istep)
> 
>      call adios_close(file_handle, adios_err)
> 
>      CALL close_end_for_group(group_prof_handle,istep)
> 
