/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>   /* _SC_PAGE_SIZE, _SC_AVPHYS_PAGES */

#if defined(__APPLE__)
#    include <mach/mach.h>
#endif

#include "core/buffer.h"
#include "core/adios_logger.h"
#include "public/adios_error.h"

#define BYTE_ALIGN 8

int adios_databuffer_resize (struct adios_file_struct *fd, uint64_t size)
{
    int retval = 0;
    // try to alloc/realloc a buffer to requested size
    // align usable buffer to BYTE_ALIGN bytes
    void * b = realloc (fd->allocated_bufptr, size +  BYTE_ALIGN - 1);
    if (b)
    {
        fd->allocated_bufptr = b;
        uint64_t p = (uint64_t) fd->allocated_bufptr;
        fd->buffer = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
        fd->buffer_size = size;

    }
    else
    {
        retval = 1;
        log_warn ("Cannot allocate %llu bytes for buffered output "
                "of group %s in adios_group_size(). Continue buffering "
                "with buffer size %llu MB\n",
                size, fd->group->name, fd->buffer_size/1048576);
    }

    return retval;
}

int adios_databuffer_free (struct adios_file_struct *fd)
{
    if (fd->allocated_bufptr)
        free (fd->allocated_bufptr);
    fd->allocated_bufptr = 0;
    fd->buffer = 0;
    fd->buffer_size = 0;
    fd->offset = 0;
    fd->bytes_written = 0;
}



// buffer sizing may be problematic.  To get a more accurate picture, check:
// http://chandrashekar.info/vault/linux-system-programs.html
static uint64_t adios_buffer_size_requested = 0;
static uint64_t adios_buffer_size_max = 0;
static uint64_t adios_buffer_size_remaining = 0;
static int adios_buffer_alloc_percentage = 0;  // 1 = yes, 0 = no
static enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when = ADIOS_BUFFER_ALLOC_UNKNOWN;

void      adios_buffer_size_requested_set (uint64_t v)  { adios_buffer_size_requested = v; }
uint64_t  adios_buffer_size_requested_get (void)        { return adios_buffer_size_requested; }
void      adios_buffer_size_max_set (uint64_t v)        { adios_buffer_size_max = v; }
void      adios_buffer_size_remaining_set (uint64_t v)  { adios_buffer_size_remaining = v; }
void      adios_buffer_alloc_percentage_set (int v)     { adios_buffer_alloc_percentage = v; }
void      adios_buffer_alloc_when_set (enum ADIOS_BUFFER_ALLOC_WHEN v)   { adios_buffer_alloc_when = v; }
enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when_get (void)   { return adios_buffer_alloc_when; }

#if defined (__APPLE__) 
// See e.g. http://www.opensource.apple.com/source/system_cmds/system_cmds-496/vm_stat.tproj/vm_stat.c
// for the code for the vm_stat command.
// http://www.opensource.apple.com/source/xnu/xnu-792.6.61/osfmk/man/host_statistics.html?txt
// describes the host_statistics function
// Added by  Dorian Krause <dorian.krause@usi.ch>
static inline size_t adios_get_avphys_pages ()
{
    // Since we are only interested in the number of free pages
    // it is fine to work with the "older" host_statistics()
    // instead of host_statistics64(). The advantage is that the
    // first function is also provided on older (e.g., Mac OS X 10.5)
    // systems
    vm_statistics_data_t   host_info;
    mach_msg_type_number_t host_info_outCnt;

    // See mach/host_info.h
    host_info_outCnt = HOST_VM_INFO_COUNT;
    if (host_statistics(mach_host_self(),
                HOST_VM_INFO,
                (host_info_t)&host_info,
                &host_info_outCnt) != KERN_SUCCESS ) {
        log_error("adios_get_avphys_pages (): host_statistics failed.\n");
        return 0;   // Best we can do
    }

    // on Mac OSX 10.4 (Tiger), there is no speculative page counting
    // VM_PAGE_QUERY_PAGE_SPECULATIVE is defined in 10.5's mach/vm_statistics.h (included in mach.h)
#   if defined (VM_PAGE_QUERY_PAGE_SPECULATIVE)
    return host_info.free_count - host_info.speculative_count;
#   else
    return host_info.free_count;
#   endif
}
#else
// See e.g. http://chandrashekar.info/vault/linux-system-programs.html
static inline size_t adios_get_avphys_pages ()
{
    return sysconf (_SC_AVPHYS_PAGES);
}
#endif

int adios_set_buffer_size ()
{
    //if (!adios_buffer_size_max) // not called before
    if (adios_buffer_size_max < adios_buffer_size_requested) // not called before
    {
        long pagesize;
        long pages;

        pagesize = sysconf (_SC_PAGE_SIZE);
        pages =  adios_get_avphys_pages ();

        if (adios_buffer_alloc_percentage)
        {
            adios_buffer_size_max =   (pages * pagesize / 100.0)
                * adios_buffer_size_requested;
        }
        else
        {
            if (pagesize * pages >= adios_buffer_size_requested)
            {
                // sufficient memory, do nothing
                adios_buffer_size_max = adios_buffer_size_requested;
            }
            else
            {
                adios_error (err_no_memory,
                        "adios_allocate_buffer (): insufficient memory: "
                        "%llu requested, %llu available.  Using "
                        "available.\n",
                        adios_buffer_size_requested,
                        (uint64_t)(((uint64_t) pagesize) * pages));
                adios_buffer_size_max = (uint64_t)((uint64_t) pagesize) * pages;
            }
        }

        adios_buffer_size_remaining = adios_buffer_size_max;

        return 1;
    }
    else
    {
        log_debug ("adios_allocate_buffer already called. No changes made.\n");
        return 1;
    }
}

uint64_t adios_method_buffer_alloc (uint64_t size)
{
    if (adios_buffer_size_remaining >= size)
    {
        adios_buffer_size_remaining -= size;

        return size;
    }
    else
    {
        uint64_t remaining = adios_buffer_size_remaining;

        adios_buffer_size_remaining = 0;

        return remaining;
    }
}

int adios_method_buffer_free (uint64_t size)
{
    if (size + adios_buffer_size_remaining > adios_buffer_size_max)
    {
        adios_error (err_invalid_buffer, 
                "ERROR: attempt to return more bytes to buffer "
                "pool than were originally available\n");

        adios_buffer_size_remaining = adios_buffer_size_max;

        return 0;
    }
    else
    {
        adios_buffer_size_remaining += size;

        return 1;
    }
}

