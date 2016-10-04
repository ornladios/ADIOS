#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>


#include "config.h"
#include "adios_clock.h"

void adios_nanosleep (int sec, int nanosec)
{
#if HAVE_NANOSLEEP
    struct timespec treq = {.tv_sec=sec, .tv_nsec=nanosec};
    struct timespec trem;
    int r;
    r = nanosleep(&treq, &trem);
    //log_debug("adios_nanosleep: Nanoslept for %d.%9.9d sec, r=%d, errno=%d\n",
    //          treq.tv_sec, treq.tv_nsec, r, errno);
    while (r == -1 && errno == EINTR) {
        treq.tv_sec = trem.tv_sec;
        treq.tv_nsec = trem.tv_nsec;
        r = nanosleep (&treq, &trem);
    }
#else
    if (sec>0) {
        //log_debug("adios_nanosleep: Slept for %d seconds\n");
        sleep(sec);
    } else {
        //log_debug("adios_nanosleep: Slept for 1 second\n");
        sleep(1);
    }

#endif
}   


#include <sys/time.h>
static struct timeval adios_timer_tp;
double adios_gettime_double()
{
    gettimeofday(&adios_timer_tp, NULL);
    return  ((double)adios_timer_tp.tv_sec + ((double)adios_timer_tp.tv_usec)/1000000.0);
}

unsigned long adios_gettime_ms()
{
    gettimeofday(&adios_timer_tp, NULL);
    unsigned long ms = lround(adios_timer_tp.tv_usec/1000) + adios_timer_tp.tv_sec*1000;
    return ms;
}

#if HAVE_CLOCK_GETTIME
    int adios_clock_gettime(clockid_t clk_id, struct timespec *ts)
    {
        return clock_gettime(clk_id, ts);
    }
#else
#   if HAVE_CLOCK_GET_TIME
#       ifdef __MACH__
#          include <mach/clock.h>
#          include <mach/mach.h>
#       else
#           error "Don't know how to compile adios_clock_gettime if this is not OSX and clock_gettime() is not present"
#       endif
        int adios_clock_gettime(clockid_t clk_id, struct timespec *ts)
        {
            // OS X does not have clock_gettime, use clock_get_time
            clock_serv_t cclock;
            mach_timespec_t mts;
            host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
            clock_get_time(cclock, &mts);
            mach_port_deallocate(mach_task_self(), cclock);
            ts->tv_sec = mts.tv_sec;
            ts->tv_nsec = mts.tv_nsec;
            return 0;
        }
#   else
#       error "Don't know how to compile adios_clock_gettime if neither clock_gettime() nor clock_get_time() is available"
#   endif
#endif



