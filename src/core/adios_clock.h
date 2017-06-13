#ifndef ADIOS_CLOCK_H_
#define ADIOS_CLOCK_H_

#include <time.h>

/* sleep for a bit */
void adios_nanosleep (int sec, int nanosec);

/* get current time as double (in seconds) */
double adios_gettime_double();

/* get current time in milliseconds as unsigned long */
unsigned long adios_gettime_ms();

#ifndef HAVE_CLOCKID_T
    typedef int clockid_t;
#endif
#ifndef CLOCK_REALTIME
#   define CLOCK_REALTIME 0
#endif
#ifndef CLOCK_MONOTONIC
#   define CLOCK_MONOTONIC 1
#endif
#ifndef CLOCK_PROCESS_CPUTIME_ID
#   define CLOCK_PROCESS_CPUTIME_ID 2
#endif

/* create a timespec variable and pass the pointer here */
int adios_clock_gettime(clockid_t clk_id, struct timespec *ts);

#endif /* ADIOS_CLOCK_H */
