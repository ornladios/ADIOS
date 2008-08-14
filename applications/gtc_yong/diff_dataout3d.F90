23a24
>   integer :: merror               ! MPI error flag
26,27c27,29
<   #define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0))
<   #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
---
>   integer :: adios_err
>   #define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0),adios_err)
>   #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b,adios_err)
29a32
>   integer*8 :: group_prof_handle 
50a54,55
> !!!     write(fdum,'("phi_dir/RUNcoord_",i5.5,".bp")')myrank_toroidal
> !!!     fdum=trim(fdum)//char(0)
52c57,59
<      write(dirstr,'("/Coordinate_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
---
> !!! modified by zf2
> !!!     write(dirstr,'("/Coordinate_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
>      write(dirstr,'("/Coordinate")') !!!(mstepall+istep)/ndiag
55,57c62,72
<      call adios_get_group(type_id,"output3d.0"//char(0))
<      call adios_set_path (type_id,dirstr)
<      call adios_open(group_handle,type_id,fdum)
---
> !!! modified by zf2 for timing
>      CALL MPI_BARRIER(toroidal_comm,merror)
> 
>      CALL open_start_for_group(group_prof_handle, "output3d.0"//char(0),istep)
> 
>      !!!call adios_get_group(type_id,"output3d.0"//char(0))
>      call adios_open(group_handle,"output3d.0"//char(0),fdum,"w"//char(0),adios_err)
>      call adios_set_path (group_handle,dirstr,adios_err)
> 
>     CALL open_end_for_group(group_prof_handle,istep)
> 
103a119,121
> 
>     CALL write_start_for_group(group_prof_handle,istep)
> 
113,115c131,137
<         ADIOS_WRITE(group_handle,radial)
<         call adios_write(group_handle,"mtheta"//char(0),(mtheta+1))
<         ADIOS_WRITE(group_handle,itran)
---
> 
> !!! added by zf2
>         if(myrank_toroidal==0) then 
>           ADIOS_WRITE(group_handle,radial)
>           call adios_write(group_handle,"mtheta"//char(0),(mtheta+1), adios_err)
>           ADIOS_WRITE(group_handle,itran)
>         endif
116a139,143
> 
> !!! added by zf2
> !!!     ADIOS_WRITE(group_handle,(mzeta+kp)*nproc_toroidal)
> !!!     ADIOS_WRITE(group_handle,(mzeta+kp)*myrank_toroidal)
> 
117a145
> 
118a147,151
> 
> !!! added by zf2
>         ADIOS_WRITE(group_handle,mgrid*nproc_toroidal)
>         ADIOS_WRITE(group_handle,mgrid*myrank_toroidal)
>     
121a155,159
> 
>     CALL write_end_for_group(group_prof_handle,istep)
> 
>     CALL close_start_for_group(group_prof_handle,istep)
> 
123c161,164
<         call adios_close(group_handle)
---
>         call adios_close(group_handle,adios_err)
> 
>     CALL close_end_for_group(group_prof_handle,istep)
> 
153a195
> !!!  write(output_fname,'("phi_dir/PHI_",i5.5,"_",i5.5,".bp")')mstepall+istep,myrank_toroidal !!(1+(mstepall+istep-ndiag)/isnap)
157,158c199,206
<   write(dirstr,'("/Potential_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
<   dirstr=trim(dirstr)//char(0)
---
> 
> !!! modified by zf2
> !!!  write(dirstr,'("/Potential_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
> !!!  dirstr=trim(dirstr)//char(0)
>   dirstr="Potential"//char(0)
> 
> !!!  write(dirstr,'("/Potential_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
> !!!  dirstr=trim(dirstr)//char(0)
160,161c208,214
<   call adios_get_group (type_id, "output3d.1"//char(0))
<   call adios_set_path (type_id,dirstr)
---
> 
> !!! modified by zf2 for timing
>   CALL MPI_BARRIER(toroidal_comm,merror)
> 
>   CALL open_start_for_group(group_prof_handle, "output3d.1"//char(0),istep)
> 
>   !!!call adios_get_group (type_id, "output3d.1"//char(0))
163c216,222
<   call adios_open (group_handle, type_id, output_fname)
---
>   call adios_open (group_handle, "output3d.1"//char(0), output_fname,"w"//char(0),adios_err)
>   call adios_set_path (group_handle,dirstr,adios_err)
> 
>   CALL open_end_for_group(group_prof_handle,istep)
> 
> !  CALL write_start_for_group(group_prof_handle,istep)
> 
176,177c235,248
<   call adios_write(group_handle,"mtheta"//char(0),(mtheta+1))
<   call adios_write (group_handle, "phi"//char(0), dataout)
---
> 
> !!! added by zf2
>   ADIOS_WRITE(group_handle,mgrid*nproc_toroidal)
>   ADIOS_WRITE(group_handle,mgrid*myrank_toroidal)
> 
>   if(myrank_toroidal==0) then
>     call adios_write(group_handle,"mtheta"//char(0),(mtheta+1),adios_err)
>   endif
>   call adios_write (group_handle, "phi"//char(0), dataout,adios_err)
> 
>   CALL write_end_for_group(group_prof_handle,istep)
> 
>   CALL close_start_for_group(group_prof_handle,istep)
> 
179a251,252
>   CALL close_end_for_group(group_prof_handle,istep)
> 
