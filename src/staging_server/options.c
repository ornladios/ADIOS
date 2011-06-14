#include "globals.h"

#include <getopt.h>
#include <string.h>
#include <errno.h>


static struct option options[] = {
    {"help",                 no_argument,          NULL,    'h'},
    {"verbose",              no_argument,          NULL,    'v'},
    {"maxmem",               required_argument,    NULL,    'm'},
    {"maxblock",             required_argument,    NULL,    'b'},
    {"logpath",              required_argument,    NULL,    'l'},
    {"logbyranks",           no_argument,          NULL,    'r'},
    {NULL,                   0,                    NULL,    0}
};


static const char *optstring = "hvm:b:l:r";

// help function
static void display_help(char * prgname) 
{
    printf(
"usage: %s [OPTIONS]\n"
"\nRun ADIOS BP staging server to which an ADIOS application can stage data\n"
"to write to disk in ADIOS BP format\n"
"\n"
"Options:\n"
"  --maxmem    | -m         Maximum buffer size to be used by one process of\n"
"                             the server (in MBs) to pull data from the\n"
"                             client application.\n"
"  --maxblock  | -b         Request maximum size for block size of one write.\n"
"                             (in MBs). To be maximum, it should be less then\n"
"                             half of maxmem and greater than any data block\n"
"                             from any client process.\n"
"  --logpath   | -l         Path for errors and verbose logging\n"
"                             default is stderr\n"
"  --logbyranks | -r        Separate log file per process\n"
"\n"
"Help options\n"
"  --help      | -h         Print this help.\n"
"  --verbose   | -v         Print log about what this program is doing.\n"
"                             Use multiple -v to increase logging level.\n"
       ,prgname);
}

static void init_globals() 
{
    user_max_memory_allowed = 0; // 0 = use all available memory
    verbose = 0; // be silent
}

/* Process arguments 
   return 0: OK
          255: exit with 0
          other: exit with other as errorcode
*/
int options_process_args( int argc, char *argv[] ) 
{
    int i, timearg=false;
    long int tmp;
    
    init_globals();
    
    /* other variables */
    int c, last_c='_';
    int last_opt = -1;
    /* Process the arguments */
    while ((c = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
        switch (c) {
            case 'b':
                errno = 0;
                tmp = strtol(optarg, (char **)NULL, 0);
                if (errno) {
                    fprintf(stderr, "Error: could not convert --maxblock value: %s\n", optarg);
                    return 1;
                }
                user_max_block_size  = (uint64_t)tmp;
                user_max_block_size *= 1048576; // argument was given as MB
                break;
            case 'h':
                display_help(argv[0]);
                return 255;
                break;
            case 'l':
                logpath = strdup(optarg);
                break;
            case 'm':
                errno = 0;
                tmp = strtol(optarg, (char **)NULL, 0);
                if (errno) {
                    fprintf(stderr, "Error: could not convert --maxmem value: %s\n", optarg);
                    return 1;
                }
                user_max_memory_allowed  = (uint64_t)tmp;
                user_max_memory_allowed *= 1048576; // argument was given as MB
                break;
            case 'r':
                logfile_separate_ranks = true;
                break;
            case 'v':
                verbose++;
                break;
            case 1:
                fprintf(stderr, "Error: Unrecognized argument: %s\n", optarg);
                return 1;
                break;
            default:
                printf("Processing default: %c\n", c);
                break;
        } /* end switch */ 
        last_c = c;
    } /* end while */
    return 0;
}


