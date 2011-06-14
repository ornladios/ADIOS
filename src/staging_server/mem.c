#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>   /* _SC_PAGE_SIZE, _SC_AVPHYS_PAGES */

#define __MEM_C__
#include "mem.h"

// Return the available memory but max the input argument
uint64_t mem_get_available (uint64_t max_allowed)
{
    long pagesize;
    long pages;
    uint64_t avail;

    pagesize = sysconf (_SC_PAGE_SIZE);
    pages    = sysconf (_SC_AVPHYS_PAGES); 

    avail = (uint64_t)pagesize * (uint64_t)pages;

    return ( max_allowed == 0 || avail <= max_allowed ? 
             avail : 
             max_allowed);
}

