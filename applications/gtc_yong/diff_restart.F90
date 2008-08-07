18,20c18,21
< #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
< #define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0))
< #define ADIOS_READ(a,b) call adios_read(a,'b'//char(0),b)
---
>   integer :: adios_err
> #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b,adios_err)
> #define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0),adios_err)
> #define ADIOS_READ(a,b) call adios_read(a,'b'//char(0),b,adios_err)
23a25
>   integer*8 :: group_prof_handle
35,40c37,49
< !  if(mod(irest,2)==0)then
< !        write(restart_fname,'(a,i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal 
< !  else
< !        write(restart_fname,'(a,i5.5,".bp")')"restart_dir2/restart_",myrank_toroidal 
< !  endif
<   write(restart_fname,'(a,i5.5,"_",i5.5,".bp")')"restart_dir/restart_",myrank_toroidal, mstepall+istep
---
> !!!  if(mod(irest,2)==0)then
> !!!        write(restart_fname,'(a,i5.5,"_",i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal,mype 
> !!!  else
> !!!        write(restart_fname,'(a,i5.5,"_",i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal,mype
> !!!  endif
> 
> !!! modified by zf2
>   if(mod(irest,2)==0)then
>         write(restart_fname,'(a,i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal
>   else
>         write(restart_fname,'(a,i5.5,".bp")')"restart_dir2/restart_",myrank_toroidal
>   endif
> 
53,55c62,69
<      call adios_get_group (group_id, "restart"//char(0))
<      ! set the path for all vars in the type for proper sizing
<      call adios_set_path (group_id,dirstr//char(0));
---
> 
> !!! modified by zf2 for timing
>      CALL MPI_BARRIER(partd_comm,merror)
> 
>      CALL open_start_for_group(group_prof_handle, "restart"//char(0),istep)
> 
>      !!! modified by zf2 for new ADIOS API 
>      !!!call adios_get_group (group_id, "restart"//char(0))
58c72,74
<         call adios_open_read (buf_id, group_id, restart_fname)
---
>         !!! modified by zf2 for new ADIOS API 
>         !!!call adios_open_read (buf_id, group_id, restart_fname)
>         call adios_open (buf_id, 'restart'//char(0), restart_fname, 'r'//char(0), adios_err)
60c76,78
<         call adios_open (buf_id, group_id, restart_fname)
---
>         !!! modified by zf2 for new ADIOS API 
>         !!!call adios_open (buf_id, group_id, restart_fname)
>         call adios_open (buf_id, 'restart'//char(0), restart_fname, 'w'//char(0),adios_err)
61a80,87
> 
>      ! set the path for all vars in the type for proper sizing
>      call adios_set_path (buf_id,dirstr//char(0),adios_err);
> 
>      CALL open_end_for_group(group_prof_handle,istep)
>  
>      CALL write_start_for_group(group_prof_handle,istep)
> 
87c113
<      call adios_write(buf_id,"zion0"//char(0),zion0(6,:))
---
>      call adios_write(buf_id,"zion0"//char(0),zion0(6,:), adios_err)
92c118
<         call adios_write(buf_id,"zelectron0"//char(0),zelectron0(6,:))
---
>         call adios_write(buf_id,"zelectron0"//char(0),zelectron0(6,:), adios_err)
96a123,127
> 
>     CALL write_end_for_group(group_prof_handle,istep)
> 
>     CALL close_start_for_group(group_prof_handle,istep)
> 
98c129,132
<      call adios_close (buf_id)
---
>      call adios_close (buf_id,adios_err)
> 
>     CALL close_end_for_group(group_prof_handle,istep)
> 
