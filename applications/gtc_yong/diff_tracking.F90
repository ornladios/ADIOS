177a178,179
>   integer :: merror ! MPI error flag
> 
182c184,186
<     #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
---
>     integer*8 :: group_prof_handle
>     integer :: adios_err
>     #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b,adios_err)
187a192
> !!!    write(fname,'("trackp_dir/TRACKP_",i5.5,"_",i5.5,".bp")')mstepall+istep,mype  !!(1+(mstepall+istep-ndiag)/isnap)
191,195c196,214
<     write(dirstr,'("node_",i5.5)')mype!,mstepall+istep !!!(mstepall+istep)/ndiag
<     dirstr=trim(dirstr)//char(0)
<     call adios_get_group (group_id, "particles"//char(0))
<     call adios_set_path (group_id,dirstr//char(0))
<     call adios_open (buf_id, group_id, fname)!!call adios_open_append(buf_id, group_id, fname)
---
> 
>     !!!modified by zf2
>     !!!write(dirstr,'("node_",i5.5)')mype!,mstepall+istep !!!(mstepall+istep)/ndiag
>     !dirstr=trim(dirstr)//char(0)
>     dirstr=trim("/")//char(0)
> 
> !!! modified by zf2 for timing
>    CALL MPI_BARRIER(comm,merror)
> 
>    CALL open_start_for_group(group_prof_handle, "particles"//char(0),istep)
> 
> !!!    call adios_get_group (group_id, "particles"//char(0))
>     call adios_open (buf_id, "particles"//char(0), fname, "w"//char(0),adios_err)!!call adios_open_append(buf_id, group_id, fname)
>     call adios_set_path (buf_id,dirstr//char(0),adios_err)
>     
>    CALL open_end_for_group(group_prof_handle,istep)
> 
>    CALL write_start_for_group(group_prof_handle,istep)
> 
204c223,228
<         call adios_write(buf_id,'ptrackede'//char(0),ptrackede(:,1:ntrackp(2)))
---
>         !!! modified by zf2
>         ADIOS_WRITE(buf_id,numberpe)  
>         ADIOS_WRITE(buf_id,nparam*numberpe)
>         ADIOS_WRITE(buf_id,nparam*mype)
> 
>         call adios_write(buf_id,'ptrackede'//char(0),ptrackede(:,1:ntrackp(2)), adios_err)
210c234,242
<     call adios_close (buf_id)
---
> 
>     CALL write_end_for_group(group_prof_handle,istep)
> 
>     CALL close_start_for_group(group_prof_handle,istep)
> 
>     call adios_close (buf_id,adios_err)
> 
>     CALL close_end_for_group(group_prof_handle,istep)
> 
