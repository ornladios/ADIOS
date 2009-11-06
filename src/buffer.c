#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>   /* _SC_PAGE_SIZE, _SC_AVPHYS_PAGES */

#include "buffer.h"

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

int adios_set_buffer_size ()
{
    if (!adios_buffer_size_max) // not called before
    {
        long pagesize;
        long pages;

#if defined (__APPLE__)
# include <sys/sysctl.h>
	int mib[2];
	uint64_t memsize;
	size_t len;
 
	mib[0] = CTL_HW;
	mib[1] = HW_MEMSIZE; /*uint64_t: physical ram size */
	len = sizeof(memsize);
	sysctl(mib, 2, &memsize, &len, NULL, 0);
	printf("- memsize  = %10i k\n", memsize/1024);
	printf("- memsize  = %10i MB\n", memsize/1024/1024);
	printf("- memsize  = %10i GB\n", memsize/1024/1024/1024);
	
	return memsize/1024;
        //return (long)memsize/1024;   /*this type cast didn't work*/	
	

#else
        pagesize = sysconf (_SC_PAGE_SIZE);
        pages = sysconf (_SC_AVPHYS_PAGES);
#endif
	
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
                fprintf (stderr, "adios_allocate_buffer (): insufficient memory: "
                                 "%llu requested, %llu available.  Using "
                                 "available.\n"
                        ,adios_buffer_size_requested
                        ,(uint64_t)(((uint64_t) pagesize) * pages)
                        );
                 adios_buffer_size_max = (uint64_t)((uint64_t) pagesize) * pages;
           }
        }

        adios_buffer_size_remaining = adios_buffer_size_max;

        return 1;
    }
    else
    {
        fprintf (stderr, "adios_allocate_buffer already called. "
                         "No changes made.\n"
                );

        return 0;
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
        fprintf (stderr, "ERROR: attempt to return more bytes to buffer "
                         "pool than were originally available\n"
                );

        adios_buffer_size_remaining = adios_buffer_size_max;

        return 0;
    }
    else
    {
        adios_buffer_size_remaining += size;

        return 1;
    }
}

