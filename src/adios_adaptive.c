#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/vfs.h>
#include <sys/ioctl.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>
#include <assert.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

static int adios_adaptive_initialized = 0;

#define PRINT_MESSAGES 1

struct adios_adaptive_data_struct
{
    int f;
    MPI_File fh;
    MPI_Request req;
    MPI_Status status;
    MPI_Comm group_comm;
    int rank;
    int size;

    struct adios_bp_buffer_struct_v1 b;

    struct adios_index_process_group_struct_v1 * old_pg_root;
    struct adios_index_var_struct_v1 * old_vars_root;
    struct adios_index_attribute_struct_v1 * old_attrs_root;

    uint64_t vars_start;
    uint64_t vars_header_size;
    uint64_t biggest_size;     // largest writer's data size (round up)
    uint16_t storage_targets;  // number of storage targets being used
    uint16_t split_groups;     // number of files to split output into
    uint16_t split_group_size; // how big the groups are at max (smaller at end)

    // these correspond to the series of parameters available to this transport
    // and control how it performs
    int16_t files_number;       // # files user is creating as part of output
    int16_t max_storage_targets;// number of OSTs in parallel filesystem
    int16_t max_stripe_count;   // max OSTs per file (system max)
    int16_t min_stripe_count;   // min OSTs per file (user desire)
    int16_t overlap_factor;     // percentage of OSTs to overlap for efficiency
    enum {split_files_unknown
         ,split_files_min
         ,split_files_max
         } split_files_count;   // how many 'split' files to generate
    int16_t split_target_count; // number of OSTs for each 'split' file

    pthread_t coordinator;     // only used if in this proc
    pthread_t sub_coordinator; // only used if in this proc
    pthread_mutex_t mutex;     // need to make sure only one thread does comm

    int groups;          // how many groups we split into
    int group;           // which group are we a member of
    int group_size;      // how large our group is
    int sub_coord_rank;  // what is the rank of our sub coordinator
    int coord_rank;      // what is the rank of the coordinator
    int stripe_size;     // how big each stripe piece is

    struct adios_file_struct * fd; // link to what was passed in
    struct adios_method_struct * method; // link to main struct

    // adaptive support stuff start
    // messaging between threads ([0] is msg, others are parameters)
    #define PARAMETER_COUNT 6
    volatile uint64_t * writer_flag;
    Queue sub_coordinator_flag;
    Queue coordinator_flag;

    pthread_mutex_t sub_coordinator_mutex;
    pthread_mutex_t coordinator_mutex;

    //volatile uint64_t * w_sub_coordinator_flag; // writer->sub_coord
    //volatile uint64_t * c_sub_coordinator_flag; // coord->sub_coord
    //volatile uint64_t * c_coordinator_flag;
    //volatile uint64_t * w_coordinator_flag;
};

#define COPY_ALL_PARAMS(dst,src) \
{ \
int i; \
for (i = 0; i < PARAMETER_COUNT; i++) \
dst [i] = src [i]; \
}

#define INIT_PARAMS(x) \
{ \
int i; \
for (i = 0; i < PARAMETER_COUNT; i++) \
x [i] = NO_FLAG; \
}


// used to message between threads without a mutex or MPI message
enum MESSAGE_FLAGS
{
     NO_FLAG                 = -1  // no msg
    ,SHUTDOWN_FLAG           = -2  // exit threads
    ,DO_WRITE_FLAG           = -3  // start a write
    ,REGISTER_COMPLETE       = -4  // registration suceed, start write
    ,WRITE_COMPLETE          = -5  // write has finished
    ,ADAPTIVE_WRITE_START    = -6  // start an adaptive write to this group
    ,WRITERS_BUSY            = -7  // all writers are either done or busy
    ,OVERALL_WRITE_COMPLETE  = -8  // all writers are done
    ,SEND_INDEX              = -9  // tell the writers to send their index
    ,REGISTER_FLAG           = -10 // register with the parent
    ,INDEX_SIZE              = -11 // size of index from sub to coord
    ,START_WRITES            = -12 // start the writing process
    ,INDEX_BODY              = -13 // the index contents
};

static
const char * message_to_string (enum MESSAGE_FLAGS m)
{
    static char x [20];
    switch (m)
    {
        case NO_FLAG: return "NO_FLAG";
        case SHUTDOWN_FLAG: return "SHUTDOWN_FLAG";
        case DO_WRITE_FLAG: return "DO_WRITE_FLAG";
        case REGISTER_COMPLETE: return "REGISTER_COMPLETE";
        case WRITE_COMPLETE: return "WRITE_COMPLETE";
        case ADAPTIVE_WRITE_START: return "ADAPTIVE_WRITE_START";
        case WRITERS_BUSY: return "WRITERS_BUSY";
        case OVERALL_WRITE_COMPLETE: return "OVERALL_WRITE_COMPLETE";
        case SEND_INDEX: return "SEND_INDEX";
        case REGISTER_FLAG: return "REGISTER_FLAG";
        case INDEX_SIZE: return "INDEX_SIZE";
        case START_WRITES: return "START_WRITES";
        case INDEX_BODY: return "INDEX_BODY";
        default: sprintf (x, "unknown (%d)", m); return x;
    }
}

static
const char * message_to_string_full (uint64_t * msg)
{
    static char x [100];
    switch (msg [0])
    {
        case NO_FLAG: return "NO_FLAG";

        case SHUTDOWN_FLAG: return "SHUTDOWN_FLAG";

        case DO_WRITE_FLAG: sprintf (x, "DO_WRITE_FLAG (tg: %lld tgr: %lld g: %lld gr: %lld offset: %lld)", msg [1], msg [2], msg [3], msg [4], msg [5]); return x;

        case REGISTER_COMPLETE: return "REGISTER_COMPLETE";

        case WRITE_COMPLETE: sprintf (x, "WRITE_COMPLETE (tg: %lld wg: %lld eo: %lld is: %lld)", msg [1], msg [2], msg [3], msg [4]); return x;

        case ADAPTIVE_WRITE_START: sprintf (x, "ADAPTIVE_WRITE_START (tg: %lld tscr: %lld o: %lld)", msg [1], msg [2], msg [3]); return x;

        case WRITERS_BUSY: sprintf (x, "WRITERS_BUSY (tg: %lld g: %lld)", msg [1], msg [2]); return x;

        case OVERALL_WRITE_COMPLETE: sprintf (x, "OVERALL_WRITE_COMPLETE (eo: %lld)", msg [1]); return x;

        case SEND_INDEX: return "SEND_INDEX";

        case REGISTER_FLAG: sprintf (x, "REGISTER_FLAG (1: %lld 2: %lld)", msg [1], msg [2]); return x;

        case INDEX_SIZE: sprintf (x, "INDEX_SIZE (g: %lld is: %lld)", msg [1], msg [2]); return x;

        case START_WRITES: return "START_WRITES";

        case INDEX_BODY: return "INDEX_BODY";

        default: sprintf (x, "unknown (%d)", msg [0]); return x;
    }
}

// used to specify which thread is the target for the MPI messages
enum MPI_TAG
{
     TAG_WRITER          = 0
    ,TAG_SUB_COORDINATOR = 1
    ,TAG_COORDINATOR     = 2
    ,TAG_FILE_OPEN       = 3
};

// just passing the set of calculated and base values rather than recalc
struct proc_struct
{
    int rank;
    int size;
    int groups;
    int group_size;
    int group;
    int coord_rank;
    int sub_coord_rank;
};

static
void calc_groups (int rank, int size, int groups
                 ,int * group, int * group_size, int * sub_coord_rank
                 ,int * coord_rank
                 )
{
    *group_size = size / groups;
    int larger_groups = size % groups;
    if (rank < larger_groups * (*group_size + 1) || !larger_groups)
    {
        if (larger_groups)
            (*group_size)++;
        *group = rank / (*group_size);
        *sub_coord_rank =   *group_size * *group
                          + *group_size - 1;
    }
    else
    {
        *group =     (larger_groups)
                   + (rank - (larger_groups * (*group_size + 1)))
                 / *group_size;
        *sub_coord_rank =   (larger_groups * (*group_size + 1) - 1)
                          + (   (*group - larger_groups + 1)
                              * *group_size
                             );
    }
    if (*group == groups - 1)
        *sub_coord_rank = size - 1;
    *coord_rank = size - 1;
}

static void buffer_write (char ** buffer, uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,const void * data, uint64_t size
                         )
{
    if (*buffer_offset + size > *buffer_size || *buffer == 0)
    {
        char * b = realloc (*buffer, *buffer_offset + size + 1000);
        if (b)
        {
            *buffer = b;
            *buffer_size = (*buffer_offset + size + 1000);
        }
        else
        {
            fprintf (stderr, "Cannot allocate memory in buffer_write.  "
                             "Requested: %llu\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}
// adaptive support stuff end

static void set_stripe_size (struct adios_adaptive_data_struct * md
                            ,struct adios_file_struct * fd
                            ,const char * filename
                            );

static void adios_var_to_comm (const char * comm_name
                              ,enum ADIOS_FLAG host_language_fortran
                              ,void * data
                              ,MPI_Comm * comm
                              )
{
    if (data)
    {
        int t = *(int *) data;

        if (!comm_name)
        {
            if (!t)
            {
                fprintf (stderr, "communicator not provided and none "
                                 "listed in XML.  Defaulting to "
                                 "MPI_COMM_SELF\n"
                        );

                *comm = MPI_COMM_SELF;
            }
            else
            {
                if (host_language_fortran == adios_flag_yes)
                {
                    *comm = MPI_Comm_f2c (t);
                }
                else
                {
                    *comm = *(MPI_Comm *) data;
                }
            }
        }
        else
        {
            if (!strcmp (comm_name, ""))
            {
                if (!t)
                {
                    fprintf (stderr, "communicator not provided and none "
                                     "listed in XML.  Defaulting to "
                                     "MPI_COMM_SELF\n"
                            );

                    *comm = MPI_COMM_SELF;
                }
                else
                {
                    if (host_language_fortran == adios_flag_yes)
                    {
                        *comm = MPI_Comm_f2c (t);
                    }
                    else
                    {
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
            else
            {
                if (!t)
                {
                    fprintf (stderr, "communicator not provided but one "
                                     "listed in XML.  Defaulting to "
                                     "MPI_COMM_WORLD\n"
                            );

                    *comm = MPI_COMM_WORLD;
                }
                else
                {
                    if (host_language_fortran == adios_flag_yes)
                    {
                        *comm = MPI_Comm_f2c (t);
                    }
                    else
                    {
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
        }
    }
    else
    {
        fprintf (stderr, "coordination-communication not provided. "
                         "Using MPI_COMM_WORLD instead\n"
                );

        *comm = MPI_COMM_WORLD;
    }
}

void adios_adaptive_init (const char * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_adaptive_data_struct * md =
                  (struct adios_adaptive_data_struct *) method->method_data;
    if (!adios_adaptive_initialized)
    {
        adios_adaptive_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_adaptive_data_struct));
    md = (struct adios_adaptive_data_struct *) method->method_data;
    md->f = -1;
    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = 0;
    md->size = 0;
    md->group_comm = MPI_COMM_NULL;
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
    md->vars_start = 0;
    md->vars_header_size = 0;
    md->biggest_size = 0;
    md->storage_targets = 0;
    md->split_groups = 1;
    md->split_group_size = 0;

    md->files_number = 1;          // always can be 1 by default
    md->max_storage_targets = -1;
    md->max_stripe_count = -1;
    md->min_stripe_count = 1;      // always can do only 1 by default
    md->overlap_factor = 0;        // always can not overlap by default
    md->split_target_count = -1;
    md->split_files_count = split_files_unknown;

    md->group = -1;
    md->sub_coord_rank = -1;
    md->fd = 0;
    md->method = 0;
    pthread_mutex_init (&md->mutex, NULL);

    md->writer_flag = malloc (8 * PARAMETER_COUNT);
    queue_init (&md->sub_coordinator_flag, free);
    queue_init (&md->coordinator_flag, free);
    pthread_mutex_init (&md->sub_coordinator_mutex, NULL);
    pthread_mutex_init (&md->coordinator_mutex, NULL);

    // parse the parameters into key=value segments for optional settings
    if (parameters)
    {
        int param_len = strlen (parameters);
        if (param_len > 0)
        {
            char * p = strdup (parameters);
            char * token = strtok (p, ";");

            while (token)
            {
                char * equal_pos = strchr (token, '=');
                if (!equal_pos)
                {
                    continue;
                }
                int equal = equal_pos - token + 1;
                int len = strlen (token);
                char * key = malloc (len);
                char * value = malloc (len);
                strncpy (key, token, equal);
                key [equal - 1] = 0;
                strncpy (value, equal_pos + 1, len - equal);
                value [len - equal] = 0;

                if (key && value)
                {
                    if (strcasecmp (key, "max_storage_targets") == 0)
                    {
                        int v = atoi (value);
                        md->max_storage_targets = v;
                    }
                    else if (strcasecmp (key, "max_stripe_count") == 0)
                    {
                        int v = atoi (value);
                        md->max_stripe_count = v;
                    }
                    else if (strcasecmp (key, "min_stripe_count") == 0)
                    {
                        int v = atoi (value);
                        md->min_stripe_count = v;
                    }
                    else if (strcasecmp (key, "files_number") == 0)
                    {
                        int v = atoi (value);
                        if (v < 1)
                        {
                            fprintf (stderr, "ADAPTIVE: files_number %d "
                                             "too small. defaulting to 1.\n"
                                    ,v
                                    );

                            v = 1;
                        }
                        md->files_number = v;
                    }
                    else if (strcasecmp (key, "overlap_factor") == 0)
                    {
                        int v = atoi (value);
                        if (v < 0)
                        {
                            fprintf (stderr, "ADAPTIVE: overlap_factor %d "
                                             "too small. defaulting to 0.\n"
                                    ,v
                                    );

                            v = 0;
                        } else
                        if (v > 99)
                        {
                            fprintf (stderr, "ADAPTIVE: overlap_factor %d "
                                             "too large. defaulting to 99.\n"
                                    ,v
                                    );

                            v = 99;
                        }
                        md->overlap_factor = v;
                    }
                    else if (strcasecmp (key, "split_target_count") == 0)
                    {
                        int v = atoi (value);
                        md->split_target_count = v;
                    }
                    else if (strcasecmp (key, "split_files_count") == 0)
                    {
                        if (strcasecmp (value, "min") == 0)
                        {
                            md->split_files_count = split_files_min;
                        }
                        else if (strcasecmp (value, "max") == 0)
                        {
                            md->split_files_count = split_files_max;
                        }
                        else
                        {
                            fprintf (stderr, "unkown split_files_count: %s "
                                             "(parameter ignored)\n"
                                    ,value
                                    );
                            md->split_files_count = split_files_unknown;
                        }
                    }
                    else
                    {
                        fprintf (stderr, "ADAPTIVE parameter: key: {%s} "
                                         "value: {%s} not recognized. Ignored\n"
                               ,key, value
                               );
                    }
                }

                free (key);
                free (value);

                token = strtok (NULL, ";");
            }

            free (p);
        }
    }

    adios_buffer_struct_init (&md->b);
}

static void * sub_coordinator_main (void * param);
static void * coordinator_main (void * param);

int adios_adaptive_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method
                   )
{
    struct adios_adaptive_data_struct * md =
                   (struct adios_adaptive_data_struct *) method->method_data;

    // we have to wait for the group_size (should_buffer) to get the comm
    // before we can do an open for any of the modes
    md->fd = fd;
    md->method = method;

    return 1;
}

static
void build_offsets (struct adios_bp_buffer_struct_v1 * b
                   ,MPI_Offset * offsets, int size, char * group_name
                   ,struct adios_index_process_group_struct_v1 * pg_root
                   )
{
    while (pg_root)
    {
        if (!strcasecmp (pg_root->group_name, group_name))
        {
            MPI_Offset size = 0;

            if (pg_root->next)
            {
                size = pg_root->next->offset_in_file - pg_root->offset_in_file;
            }
            else
            {
                size = b->pg_index_offset - pg_root->offset_in_file;
            }

            offsets [pg_root->process_id * 3] = pg_root->offset_in_file;
            offsets [pg_root->process_id * 3 + 1] = size;
            offsets [pg_root->process_id * 3 + 2] = b->version;
        }

        pg_root = pg_root->next;
    }
}

static void
adios_build_file_offset (struct adios_adaptive_data_struct *md
                        ,struct adios_file_struct *fd, char *name)
{
#define SCATTER_PARAMS 5
    if (md->group_comm != MPI_COMM_NULL)
    {
        if (md->rank == 0)
        {
            // make one space for offset and one for size
            MPI_Offset * offsets = malloc(sizeof (MPI_Offset)
                                           * md->size * SCATTER_PARAMS);
            int i;

            offsets [0] = fd->write_size_bytes;
            MPI_Gather (offsets, 1, MPI_LONG_LONG
                       ,offsets, 1, MPI_LONG_LONG
                       ,0, md->group_comm);

            // find the largest and use that as a basis for stripe
            // size for each process writing
            uint64_t biggest_size = 0;
            for (i = 0; i < md->size; i++)
            {
                if (offsets [i] > biggest_size)
                    biggest_size = offsets [i];
            }
            // now round up to the next stripe size increment (Lustre: 64 KB)
#define STRIPE_INCREMENT (64 * 1024)
            if (biggest_size % (STRIPE_INCREMENT))
            {
                biggest_size = (  ((biggest_size / STRIPE_INCREMENT) + 1) 
                                * STRIPE_INCREMENT
                               );
            }
            // also round up the base_offset, just in case
            if (fd->base_offset % (STRIPE_INCREMENT))
            {
                fd->base_offset = (  ((fd->base_offset / STRIPE_INCREMENT) + 1) 
                                   * STRIPE_INCREMENT
                                  );
            }
            md->biggest_size = biggest_size;
#undef STRIPE_INCREMENT
            md->groups = (md->size > md->max_storage_targets) ? md->max_storage_targets : md->size;
            md->split_groups = md->groups;
            calc_groups (md->rank, md->size, md->groups
                        ,&md->group, &md->group_size, &md->sub_coord_rank
                        ,&md->coord_rank
                        );
            set_stripe_size (md, fd, name);
            offsets [1] = biggest_size;
            offsets [2] = md->storage_targets;
            offsets [3] = md->split_groups;
            // need to use this calc since last grouping may be incomplete
            // with a previous group having a piece later in the file.
            //md->b.pg_index_offset =   last_offset + biggest_size;

            md->stripe_size = biggest_size;
            MPI_Scatter (offsets, SCATTER_PARAMS, MPI_LONG_LONG
                        ,offsets, SCATTER_PARAMS, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            fd->base_offset = offsets [0];
            fd->pg_start_in_file = fd->base_offset;
            free (offsets);
        }
        else
        {
            MPI_Offset offset [SCATTER_PARAMS];
            offset [0] = fd->write_size_bytes;

            MPI_Gather (offset, 1, MPI_LONG_LONG
                       ,offset, 1, MPI_LONG_LONG
                       ,0, md->group_comm
                       );

            MPI_Scatter (offset, SCATTER_PARAMS, MPI_LONG_LONG
                        ,offset, SCATTER_PARAMS, MPI_LONG_LONG
                        ,0, md->group_comm
                        );

            md->stripe_size = offset [1];
            md->storage_targets = offset [2];
            md->split_groups = offset [3];
            md->groups = offset [3];
            md->groups = (md->size > md->max_storage_targets) ? md->max_storage_targets : md->size;
            calc_groups (md->rank, md->size, md->groups
                        ,&md->group, &md->group_size, &md->sub_coord_rank
                        ,&md->coord_rank
                        );

            fd->pg_start_in_file = fd->base_offset;
        }
    }
    else
    {
        md->b.pg_index_offset = fd->write_size_bytes;
    }
}

// calc_stripe_info figures out how many OSTs to use based on either an
// assumption of one file for all processes or a different count if the
// parameters are supplied in the XML file. These parameters are defined as
// follows:

// files_number - the number of files being written simultaneously. This also
//                implies that all MPI processes are involved in the write.
// max_storage_targets - the number of OSTs in the system.
//                       On ewok, this is 12. On jaguarpf, this is 672.
// max_stripe_count - the maximum number of OSTs available for a single file.
//                    This is 160 on jaguarpf and 12 on ewok.
// min_stripe_count - the fewest OSTs to use per file. This will override the
//                    overlap_factor, if necessary.
// overlap_factor - what percentage of the allocated portion of the OSTs should
//                  be allowed to overlap with the next file set. For example,
//                  a value of 50 means that half of the next set of OSTs will
//                  also be used for this set (e.g., set 0= 0-15, set 1= 10-25,
//                  set 2=20-35, ...).

// The way to put it into the XML, which is currently pretty unforgiving, is
// like this:
// <transport method="ADAPTIVE" group="restart">max_storage_targets=672;max_stripe_count=160;files_number=3;overlap_factor=50;min_stripe_count=10</transport>

// This says to divide the 672 OSTs so that there are 16 groups. Each group
// will consist of 1/3 portion of the whole plus 1/6 as an overlap factor.
// This will not exceed the max_stripe_count. (672/3 = 224 + 50% = 336, but
// max is 160 so limited to 160 [0-159, 112-272, 225-385, etc.]).

// If this were set for 600 files, then it would be 10 OSTs per file at an
// offset of 2 [0-9, 2-11, 4-13, etc.]

// Initial assumption is that based on the rank, we can guess which set of OSTs
// to use. If this won't work, then we need to add a communication to exchange
// information so that the processes can make that decision (and MPI_All_to_all
// of the rank is sufficient so that each main process can determine the
// ordering and then calculate which set to use).

// Once we figure out the stripe info, we may determine that we are not using
// all possible OSTs and therefore limiting parallelism in writing. To address
// this, two other options may be supplied
// split_target_count - the number of OSTs to use if maximizing the number
//                      of files we create.
// split_files_count - {min|max} Either maximize the number of files (using
//                     the split_target_count to guide the number (no shared
//                     OSTs) or minimize while still using all of the OSTs.
static void calc_stripe_info (struct adios_adaptive_data_struct * md
                             ,int * stripe_offset
                             ,int * stripe_count
                             )
{
    // these 2 are the required parameters to decide to do stripe manipulation
    if (   md->max_storage_targets > 0
        && md->max_stripe_count > 0
       )
    {
        int targets_per_file = md->max_storage_targets / md->files_number;
        if (md->max_storage_targets % md->files_number)
            targets_per_file++;

        int overlap_count = (int) (  targets_per_file
                                   * (md->overlap_factor/100.0)
                                  );

        int net_targets_per_file = targets_per_file + overlap_count;

        if (net_targets_per_file < md->min_stripe_count)
            net_targets_per_file = md->min_stripe_count;
        if (net_targets_per_file > md->max_stripe_count)
            net_targets_per_file = md->max_stripe_count;

        int range = 0;
        int number = 0;

        range = md->size / md->files_number;
        if (md->size % md->files_number)
            range++;
        // leave number the same so that we don't falsely overlap
        number = md->rank / range;

        // offset should start based on the split position or number of files
        // if we split, then this is the base and size of this range of files
        *stripe_offset = targets_per_file * number;
        *stripe_count = net_targets_per_file;

        // split, if we must
        if (md->split_files_count != split_files_unknown)
        {
            switch (md->split_files_count)
            {
                case split_files_min:
                {
                     // if we can split to our advantage, do so
                     if (  net_targets_per_file * md->files_number
                         < md->max_storage_targets
                        )
                     {
                         int groups = md->max_storage_targets / *stripe_count;
                         if (md->max_storage_targets % *stripe_count)
                             groups++;
                         if (md->max_storage_targets > md->size)
                         {
                             groups = md->size / *stripe_count;
                             if (md->size % *stripe_count)
                                 groups++;
                         }
                         int group = md->rank / (md->size / groups);
                         int group_size = md->size / groups;
                         int group_rank = md->rank % group_size;

                         md->split_groups = groups;
                         *stripe_count = group_size;
                     }
                     break;
                }

                case split_files_max:
                {
                    // if we can split to our advantage, do so
                    if (net_targets_per_file > md->split_target_count)
                    {
                        int groups =   md->max_storage_targets
                                     / md->split_target_count;
                        if (md->max_storage_targets % md->split_target_count)
                            groups++;
                        if (md->max_storage_targets > md->size)
                        {
                            groups = md->size / md->split_target_count;
                            if (md->size % md->split_target_count)
                                groups++;
                        }
                        int group = md->rank / (md->size / groups);

                        md->split_groups = groups;
                        *stripe_count = md->split_target_count;
                    }
                    break;
                }
            }
        }
    }
    else
    {
        *stripe_count = UINT16_MAX;
    }
}

// LUSTRE Structure
// from /usr/include/lustre/lustre_user.h
#define LUSTRE_SUPER_MAGIC 0x0BD00BD0
#  define LOV_USER_MAGIC 0x0BD10BD0
#  define LL_IOC_LOV_SETSTRIPE  _IOW ('f', 154, long)
#  define LL_IOC_LOV_GETSTRIPE  _IOW ('f', 155, long)
#define O_LOV_DELAY_CREATE 0100000000

struct lov_user_ost_data {           // per-stripe data structure
        uint64_t l_object_id;        // OST object ID
        uint64_t l_object_gr;        // OST object group (creating MDS number)
        uint32_t l_ost_gen;          // generation of this OST index
        uint32_t l_ost_idx;          // OST index in LOV
} __attribute__((packed));
struct lov_user_md {                 // LOV EA user data (host-endian)
        uint32_t lmm_magic;          // magic number = LOV_USER_MAGIC_V1
        uint32_t lmm_pattern;        // LOV_PATTERN_RAID0, LOV_PATTERN_RAID1
        uint64_t lmm_object_id;      // LOV object ID
        uint64_t lmm_object_gr;      // LOV object group
        uint32_t lmm_stripe_size;    // size of stripe in bytes
        uint16_t lmm_stripe_count;   // num stripes in use for this object
        uint16_t lmm_stripe_offset;  // starting stripe offset in lmm_objects
        struct lov_user_ost_data  lmm_objects[0]; // per-stripe data
} __attribute__((packed));

// do the magic ioctl calls to set Lustre's stripe size
// generate the number of groups, group size, and number of OSTs per file.
// These are returned in md->split_groups, md->split_size, & md->storage_targets
static void set_stripe_size (struct adios_adaptive_data_struct * md
                            ,struct adios_file_struct * fd
                            ,const char * filename
                            )
{
    struct statfs fsbuf;
    int err;

    int f;
    int old_mask;
    int perm;

    old_mask = umask (0);
    umask (old_mask);
    perm = old_mask ^ 0666;

    if (fd->mode == adios_mode_write)
        unlink (filename); // cleanup old stuff
    f = open (filename, O_RDONLY | O_CREAT | O_LOV_DELAY_CREATE, perm);
    // Note: Since each file might have different write_buffer,
    // So we will reset write_buffer even buffer_size != 0
    err = statfs (filename, &fsbuf);
    if (!err && fsbuf.f_type == LUSTRE_SUPER_MAGIC)
    {
        if (f != -1)
        {
            int i;
            int stripe_count;
            int stripe_offset;
            struct lov_user_md lum;
            lum.lmm_magic = LOV_USER_MAGIC;
            // get what Lustre assigns by default
            err = ioctl (f, LL_IOC_LOV_GETSTRIPE, (void *) &lum);
            stripe_count = lum.lmm_stripe_count;
            stripe_offset = lum.lmm_stripe_offset;
            calc_stripe_info (md, &stripe_offset, &stripe_count);

            // fixup for our desires
            lum.lmm_magic = LOV_USER_MAGIC;
            lum.lmm_pattern = 0;
            lum.lmm_stripe_size = md->biggest_size;
            lum.lmm_stripe_count = stripe_count; // from calc_stripe_info
            lum.lmm_stripe_offset = stripe_offset;
            err = ioctl (f, LL_IOC_LOV_SETSTRIPE, (void *) &lum);
            lum.lmm_stripe_count = 0;
            err = ioctl (f, LL_IOC_LOV_GETSTRIPE, (void *) &lum);
            // if err != 0, the must not be Lustre
            if (err == 0)
            {
                md->storage_targets = lum.lmm_stripe_count;
            }
            close (f);

#if 0
            // figure out if we need to split to use all of the OSTs (round up)
            int targets_per_file;
            if (md->storage_targets > md->max_stripe_count)
                targets_per_file = md->max_stripe_count;
            else
                targets_per_file = md->storage_targets;
            // assuming our max targets per file, see if we have to split
	    md->split_groups =   (md->max_storage_targets / md->files_number)
                               / targets_per_file;
            if ((md->max_storage_targets / md->files_number) % targets_per_file)
                md->split_groups++;

            // make sure we don't make more splits than we can handle with
            // our min stripe count
            if (md->size / md->split_groups < md->min_stripe_count)
                md->split_groups = md->size / md->min_stripe_count;
            if (md->size % md->min_stripe_count)
                md->split_groups++;

            // if we are splitting files, fixup now (need to do above to get
            // the real max storage targets from the filesystem).
            if (   md->split_groups > 1
                || md->split_files_count != split_files_unknown
               )
            {
                // adjust the split_groups based on how we want to split
                // if not specified (or min), just use all of the OSTs
                // in as few files as possible
                if (md->split_files_count == split_files_max)
                {
                    if (md->size / md->split_groups < md->split_target_count)
                    {
                        md->split_groups = md->size / md->split_target_count;
                        if (md->size % md->split_target_count)
                            md->split_groups++;
                    }
                }

                // adjust the storage targets to the number of groups
                md->storage_targets =   md->max_storage_targets
                                      / md->split_groups;
                if (md->storage_targets < md->min_stripe_count)
                    md->storage_targets = md->min_stripe_count;
                if (md->storage_targets > md->max_stripe_count)
                    md->storage_targets = md->max_stripe_count;
#else
                md->storage_targets = 1;
#endif

                int * f_split = malloc (sizeof (int) * md->split_groups);
                char split_format [7] = ".%d";
                char split_name [7];
                char * new_name = malloc (strlen (filename) + 7 + 1);
                for (i = 0; i < md->split_groups; i++)
                {
                    sprintf (split_name, split_format, i);
                    sprintf (new_name, "%s%s", filename, split_name);
                    if (fd->mode == adios_mode_write)
                        unlink (new_name);  // clean up old stuff
                    f_split [i] = open (new_name
                                       ,O_RDONLY | O_CREAT | O_LOV_DELAY_CREATE
                                       ,perm
                                       );
                }

                for (i = 0; i < md->split_groups; i++)
                {
                    lum.lmm_magic = LOV_USER_MAGIC;
                    lum.lmm_pattern = 0;
                    lum.lmm_stripe_size = md->biggest_size;
                    lum.lmm_stripe_count = md->storage_targets;
                    lum.lmm_stripe_offset =   stripe_offset
                                            + i * md->storage_targets;
                    err = ioctl (f_split [i], LL_IOC_LOV_SETSTRIPE
                                ,(void *) &lum
                                );
                    close (f_split [i]);
                }
                free (f_split);
                unlink (filename);
#if 0
            }
#endif
        }
    }
}

static void setup_threads_and_register (struct adios_adaptive_data_struct * md
                                       ,struct adios_file_struct * fd
                                       )
{
    uint64_t msg [PARAMETER_COUNT];
    MPI_Request req;
    MPI_Status status;
    
    // we can get the threads setup and get everyone registered though.
    int i;
    for (i = 0; i < PARAMETER_COUNT; i++)
    {
        md->writer_flag [i] = NO_FLAG;
    }

    // spawn worker threads for coordination
    if (md->rank == md->coord_rank) 
    {
        pthread_create (&md->coordinator, NULL, coordinator_main, (void *) md);
    }   
    if (md->rank == md->sub_coord_rank)
    {
        pthread_create (&md->sub_coordinator, NULL, sub_coordinator_main
                       ,(void *) md
                       );
    }                  

    // register writer with sub_coordinator
    if (md->rank != md->sub_coord_rank)
    {
        uint64_t msgx [PARAMETER_COUNT];
        msgx [0] = REGISTER_FLAG;
        msgx [1] = md->rank;
        pthread_mutex_lock (&md->mutex);
        MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                  ,md->sub_coord_rank, TAG_SUB_COORDINATOR
                  ,md->group_comm, &req
                  );

        pthread_mutex_unlock (&md->mutex);
    }
    else
    {
        uint64_t * flag = malloc (8 * PARAMETER_COUNT);
        INIT_PARAMS(flag);
        flag [0] = REGISTER_FLAG;
        flag [1] = md->rank;
        pthread_mutex_lock (&md->sub_coordinator_mutex);
        queue_enqueue (&md->sub_coordinator_flag, flag);
        pthread_mutex_unlock (&md->sub_coordinator_mutex);
    }
}

enum ADIOS_FLAG adios_adaptive_should_buffer (struct adios_file_struct * fd
                                        ,struct adios_method_struct * method
                                        ,void * comm
                                        )
{
    int i;
    struct adios_adaptive_data_struct * md =
                   (struct adios_adaptive_data_struct *) method->method_data;
    char * name;
    int err;
    int flag;    // used for coordinating the MPI_File_open

    int previous;
    int current;
    int next;

    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

    adios_var_to_comm (fd->group->group_comm
                      ,fd->group->adios_host_language_fortran
                      ,comm
                      ,&md->group_comm
                      );
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }
    fd->group->process_id = md->rank;

    if (md->rank == md->size - 1)
        next = -1;
    else
        next = md->rank + 1;
    previous = md->rank - 1;
    current = md->rank;

    fd->base_offset = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL, &md->fh
                                    );
                if (err != MPI_SUCCESS)
                {
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    fprintf (stderr, "MPI open read failed for %s: '%s'\n"
                            ,name, e
                            );
                    free (name);

                    return adios_flag_no;
                }

                MPI_Offset file_size;
                MPI_File_get_size (md->fh, &file_size);
                md->b.file_size = file_size;

                adios_init_buffer_read_version (&md->b);
                MPI_File_seek (md->fh, md->b.file_size - md->b.length
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_version (&md->b, &md->b.version);

                adios_init_buffer_read_index_offsets (&md->b);
                // already in the buffer
                adios_parse_index_offsets_v1 (&md->b);

                adios_init_buffer_read_process_group_index (&md->b);
                MPI_File_seek (md->fh, md->b.pg_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.pg_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_process_group_index_v1 (&md->b
                                                   ,&md->old_pg_root
                                                   );

#if 1
                adios_init_buffer_read_vars_index (&md->b);
                MPI_File_seek (md->fh, md->b.vars_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                adios_init_buffer_read_attributes_index (&md->b);
                MPI_File_seek (md->fh, md->b.attrs_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.attrs_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_attributes_index_v1 (&md->b, &md->old_attrs_root);
#endif

                fd->base_offset = md->b.end_of_pgs;
            }

            if (   md->group_comm != MPI_COMM_NULL
                && md->group_comm != MPI_COMM_SELF
               )
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size * 3
                                                  );
                    memset (offsets, 0, sizeof (MPI_Offset) * md->size * 3);

                    // go through the pg index to build the offsets array
                    build_offsets (&md->b, offsets, md->size
                                  ,fd->group->name, md->old_pg_root
                                  );
                    MPI_Scatter (offsets, 3, MPI_LONG_LONG
                                ,offsets, 3, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
                    md->b.read_pg_offset = offsets [0];
                    md->b.read_pg_size = offsets [1];
                    free (offsets);
                }
                else
                {
                    MPI_Offset offset [3];
                    offset [0] = offset [1] = offset [2] = 0;

                    MPI_Scatter (offset, 3, MPI_LONG_LONG
                                ,offset, 3, MPI_LONG_LONG
                                ,0, md->group_comm
                                );

                    md->b.read_pg_offset = offset [0];
                    md->b.read_pg_size = offset [1];
                    md->b.version = offset [2];
                }
            }

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // note rank 0 is already open
                // don't open it again here

                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mutex);
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mutex);
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                        ,name, e
                        );
                free (name);

                return adios_flag_no;
            }

            break;
        }

        case adios_mode_write:
        {
            fd->base_offset = 0;
            fd->pg_start_in_file = 0;

            // this needs to be before our build_file_offsets because
            // that process will attempt to create the file from scratch
            // and set the stripe parameters
            if (previous == -1)
            {
                unlink (name); // make sure clean
            }

            // figure out the offsets and create the file with proper striping
            // before the MPI_File_open is called
            adios_build_file_offset (md, fd, name);
            setup_threads_and_register (md, fd);

            if (md->split_groups != 1)
            {
                // if we need to do a file split, we need to fixup the name
                name = realloc (name, (  strlen (method->base_path)
                                       + strlen (fd->name) + 1 + 6
                                      )
                               ); // 6 extra for '.XXXXX' file number
                char split_format [7] = ".%d";
                char split_name [7];
                sprintf (split_name, split_format, md->group);
                strcat (name, split_name);
            }

#define BLAST_OPENS 1
#if BLAST_OPENS
            md->f = open (name, O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC);
#else
            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                md->f = open (name, O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC);
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, TAG_FILE_OPEN
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mutex);
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INTEGER, previous, TAG_FILE_OPEN
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, TAG_FILE_OPEN
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mutex);
                }
                md->f = open (name, O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC);
            }
            if (next != -1)
                MPI_Wait (&md->req, &md->status);
#endif
            err = MPI_SUCCESS;

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                        ,name, e
                        );
                free (name);

                return adios_flag_no;
            }

            break;
        }

        case adios_mode_append:
        {
            int old_file = 1;
            adios_buffer_struct_clear (&md->b);

            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL, &md->fh
                                    );

                if (err != MPI_SUCCESS)
                {
                    old_file = 0;
                    err = MPI_File_open (MPI_COMM_SELF, name
                                        ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                        ,MPI_INFO_NULL, &md->fh
                                        );
                    if (err != MPI_SUCCESS)
                    {
                        char e [MPI_MAX_ERROR_STRING];
                        int len = 0;
                        memset (e, 0, MPI_MAX_ERROR_STRING);
                        MPI_Error_string (err, e, &len);
                        fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                                ,name, e
                                );
                        free (name);

                        return adios_flag_no;
                    }
                }
                MPI_Bcast (&old_file, 1, MPI_INTEGER, 0, md->group_comm);
            }
            else
            {
                if (md->group_comm != MPI_COMM_NULL)
                    MPI_Bcast (&old_file, 1, MPI_INTEGER, 0, md->group_comm);
            }

            if (old_file)
            {
                if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
                {
                    if (err != MPI_SUCCESS)
                    {
                        md->b.file_size = 0;
                    }
                    else
                    {
                        MPI_Offset file_size;
                        MPI_File_get_size (md->fh, &file_size);
                        md->b.file_size = file_size;
                    }

                    adios_init_buffer_read_version (&md->b);
                    MPI_File_seek (md->fh, md->b.file_size - md->b.length
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_version (&md->b, &md->b.version);

                    adios_init_buffer_read_index_offsets (&md->b);
                    // already in the buffer
                    adios_parse_index_offsets_v1 (&md->b);

                    adios_init_buffer_read_process_group_index (&md->b);
                    MPI_File_seek (md->fh, md->b.pg_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.pg_size, MPI_BYTE
                                  ,&md->status
                                  );

                    adios_parse_process_group_index_v1 (&md->b
                                                       ,&md->old_pg_root
                                                       );

                    // find the largest time index so we can append properly
                    struct adios_index_process_group_struct_v1 * p;
                    uint32_t max_time_index = 0;
                    p = md->old_pg_root;
                    while (p)
                    {
                        if (p->time_index > max_time_index)
                            max_time_index = p->time_index;
                        p = p->next;
                    }
                    fd->group->time_index = ++max_time_index;
                    MPI_Bcast (&fd->group->time_index, 1, MPI_INTEGER, 0
                              ,md->group_comm
                              );

                    adios_init_buffer_read_vars_index (&md->b);
                    MPI_File_seek (md->fh, md->b.vars_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                    adios_init_buffer_read_attributes_index (&md->b);
                    MPI_File_seek (md->fh, md->b.attrs_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.attrs_size
                                  ,MPI_BYTE, &md->status
                                  );
                    adios_parse_attributes_index_v1 (&md->b
                                                    ,&md->old_attrs_root
                                                    );

                    fd->base_offset = md->b.end_of_pgs;
                    fd->pg_start_in_file = fd->base_offset;
                }
                else
                {
                    fd->base_offset = 0;
                    fd->pg_start_in_file = 0;
                    MPI_Bcast (&fd->group->time_index, 1, MPI_INTEGER, 0
                              ,md->group_comm
                              );

                }

                MPI_File_close (&md->fh);
            }
            else
            {
                fd->base_offset = 0;
                fd->pg_start_in_file = 0;
            }

            // figure out the offsets and create the file with proper striping
            // before the MPI_File_open is called
            adios_build_file_offset (md, fd, name);

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // we know it exists, because we created it if it didn't
                // when reading the old file so can just open wronly
                // but adding the create for consistency with write mode
                // so it is easier to merge write/append later
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mutex);
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mutex);
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                        ,name, e
                        );
                free (name);

                return adios_flag_no;
            }

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);

            free (name);

            return adios_flag_no;
        }
    }

    free (name);

    if (fd->shared_buffer == adios_flag_no && fd->mode != adios_mode_read)
    {
        // write the process group header
        adios_write_process_group_header_v1 (fd, fd->write_size_bytes);

        MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
        MPI_File_write (md->fh, fd->buffer, fd->bytes_written, MPI_BYTE
                       ,&md->status
                       );
        int count;
        MPI_Get_count (&md->status, MPI_BYTE, &count);
        if (count != fd->bytes_written)
        {
            fprintf (stderr, "a:MPI method tried to write %llu, "
                             "only wrote %d\n"
                    ,fd->bytes_written
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);

        // setup for writing vars
        adios_write_open_vars_v1 (fd);
        md->vars_start = fd->base_offset;
        md->vars_header_size = fd->offset;
        fd->base_offset += fd->offset;
        MPI_File_seek (md->fh, md->vars_header_size, MPI_SEEK_CUR);
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);
    }

    return fd->shared_buffer;
}

void adios_adaptive_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )
{
    struct adios_adaptive_data_struct * md =
                     (struct adios_adaptive_data_struct *) method->method_data;

    if (v->got_buffer == adios_flag_yes)
    {
        if (data != v->data)  // if the user didn't give back the same thing
        {
            if (v->free_data == adios_flag_yes)
            {
                free (v->data);
                adios_method_buffer_free (v->data_size);
            }
        }
        else
        {
            // we already saved all of the info, so we're ok.
            return;
        }
    }

    if (fd->shared_buffer == adios_flag_no)
    {
        // var payload sent for sizing information
        adios_write_var_header_v1 (fd, v);

        MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                       ,MPI_BYTE, &md->status
                       );
        int count;
        MPI_Get_count (&md->status, MPI_BYTE, &count);
        if (count != fd->bytes_written)
        {
            fprintf (stderr, "b:MPI method tried to write %llu, "
                             "only wrote %d\n"
                    ,fd->bytes_written
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);

        // write payload
        // adios_write_var_payload_v1 (fd, v);
        uint64_t var_size = adios_get_var_size (v, fd->group, v->data);
        MPI_File_write (md->fh, v->data, var_size, MPI_BYTE, &md->status);
        MPI_Get_count (&md->status, MPI_BYTE, &count);
        if (count != var_size)
        {
            fprintf (stderr, "c:MPI method tried to write %llu, "
                             "only wrote %d\n"
                    ,var_size
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);
    }
}

void adios_adaptive_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    uint64_t mem_allowed;

    if (*size == 0)
    {
        *buffer = 0;

        return;
    }

    if (v->data && v->free_data)
    {
        adios_method_buffer_free (v->data_size);
        free (v->data);
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size)
    {
        *buffer = malloc (*size);
        if (!*buffer)
        {
            adios_method_buffer_free (mem_allowed);
            fprintf (stderr, "Out of memory allocating %llu bytes for %s\n"
                    ,*size, v->name
                    );
            v->got_buffer = adios_flag_no;
            v->free_data = adios_flag_no;
            v->data_size = 0;
            v->data = 0;
            *size = 0;
            *buffer = 0;
        }
        else
        {
            v->got_buffer = adios_flag_yes;
            v->free_data = adios_flag_yes;
            v->data_size = mem_allowed;
            v->data = *buffer;
        }
    }
    else
    {
        adios_method_buffer_free (mem_allowed);
        fprintf (stderr, "OVERFLOW: Cannot allocate requested buffer of %llu "
                         "bytes for %s\n"
                ,*size
                ,v->name
                );
        *size = 0;
        *buffer = 0;
    }
}

void adios_adaptive_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
    v->data_size = buffer_size;
}

static void adios_mpi_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
    struct adios_adaptive_data_struct * md =
                    (struct adios_adaptive_data_struct *) method->method_data;
    struct adios_var_struct * v = fd->group->vars;

    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;

    switch (md->b.version)
    {
        case 1:
        {
            // the three section headers
            struct adios_process_group_header_struct_v1 pg_header;
            struct adios_vars_header_struct_v1 vars_header;
            struct adios_attributes_header_struct_v1 attrs_header;

            struct adios_var_header_struct_v1 var_header;
            struct adios_var_payload_struct_v1 var_payload;
            struct adios_attribute_struct_v1 attribute;

            int i;

            adios_init_buffer_read_process_group (&md->b);
            MPI_File_seek (md->fh, md->b.read_pg_offset
                          ,MPI_SEEK_SET
                          );
            MPI_File_read (md->fh, md->b.buff, md->b.read_pg_size, MPI_BYTE
                          ,&md->status
                          );
            adios_parse_process_group_header_v1 (&md->b, &pg_header);

            adios_parse_vars_header_v1 (&md->b, &vars_header);

            for (i = 0; i < vars_header.count; i++)
            {
                memset (&var_payload, 0
                       ,sizeof (struct adios_var_payload_struct_v1)
                       );
                adios_parse_var_data_header_v1 (&md->b, &var_header);

                struct adios_var_struct * v1 = v;
                while (v1)
                {
                    if (   strcasecmp (var_header.name, v1->name)
                        || strcasecmp (var_header.path, v1->path)
                       )
                    {
                        v1 = v1->next;
                    }
                    else
                    {
                        break;
                    }
                }

                if (v1)
                {
                    var_payload.payload = v1->data;
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,&var_payload
                                                    ,v1->data_size
                                                    );
                }
                else
                {
                    printf ("MPI read: skipping name: %s path: %s\n"
                           ,var_header.name, var_header.path
                           );
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,NULL, 0
                                                    );
                }

                adios_clear_var_header_v1 (&var_header);
            }

#if 1
            adios_parse_attributes_header_v1 (&md->b, &attrs_header);

            for (i = 0; i < attrs_header.count; i++)
            {
                adios_parse_attribute_v1 (&md->b, &attribute);
                adios_clear_attribute_v1 (&attribute);
            }
#endif
            adios_clear_process_group_header_v1 (&pg_header);

            break;
        }

        default:
            fprintf (stderr, "MPI read: file version unknown: %u\n"
                    ,md->b.version
                    );
            return;
    }

    adios_buffer_struct_clear (&md->b);
}

void adios_adaptive_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_adaptive_data_struct * md =
                (struct adios_adaptive_data_struct *) method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            // read the index to find the place to start reading
            adios_mpi_do_read (fd, method);
            struct adios_var_struct * v = fd->group->vars;
            while (v)
            {
                v->data = 0;
                v = v->next;
            }

            break;
        }

        case adios_mode_write:
        {
            MPI_Request req;
            MPI_Status status;
            uint64_t msg [PARAMETER_COUNT];
            int err;
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset;

            // setup so we can use split files or not
            int index_rank = md->rank;
            int index_size = md->size;
            MPI_Comm index_comm = md->group_comm;

            if (fd->shared_buffer == adios_flag_no)
            {
                MPI_Offset new_off;
                // set it up so that it will start at 0, but have correct sizes
                MPI_File_get_position (md->fh, &new_off);
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_vars_v1 (fd);
                // fd->vars_start gets updated with the size written
                MPI_File_seek (md->fh, md->vars_start, MPI_SEEK_SET);
                MPI_File_write (md->fh, fd->buffer, md->vars_header_size
                               ,MPI_BYTE, &md->status
                               );
                int count;
                MPI_Get_count (&md->status, MPI_BYTE, &count);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "d:MPI method tried to write %llu, "
                                     "only wrote %d\n"
                            ,md->vars_header_size
                            ,count
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
                adios_shared_buffer_free (&md->b);

                adios_write_open_attributes_v1 (fd);
                md->vars_start = new_off;
                md->vars_header_size = fd->offset;
                MPI_File_seek (md->fh, new_off + md->vars_header_size
                              ,MPI_SEEK_SET
                              ); // go back to end, but after attr header
                fd->base_offset += fd->offset;  // add size of header
                fd->offset = 0;
                fd->bytes_written = 0;

                while (a)
                {
                    adios_write_attribute_v1 (fd, a);
                    MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                                   ,MPI_BYTE, &md->status
                                   );
                    MPI_Get_count (&md->status, MPI_BYTE, &count);
                    if (count != fd->bytes_written)
                    {
                        fprintf (stderr, "e:MPI method tried to write %llu, "
                                         "only wrote %d\n"
                                ,fd->bytes_written
                                ,count
                                );
                    }
                    fd->base_offset += count;
                    fd->offset = 0;
                    fd->bytes_written = 0;
                    adios_shared_buffer_free (&md->b);

                    a = a->next;
                }

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);
                MPI_File_seek (md->fh, md->vars_start, MPI_SEEK_SET);
                // fd->vars_start gets updated with the size written
                MPI_File_write (md->fh, fd->buffer, md->vars_header_size
                               ,MPI_BYTE, &md->status
                               );
                MPI_Get_count (&md->status, MPI_BYTE, &count);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "f:MPI method tried to write %llu, "
                                     "only wrote %d\n"
                            ,md->vars_header_size
                            ,count
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
            }

            char * index_buffer = 0;
            uint64_t index_buffer_size = 0;
            uint64_t index_buffer_offset = 0;

            // registration should be complete, so start writing
            if (md->rank == md->coord_rank)
            {
                // wait for registration to complete
                while (md->writer_flag [0] != REGISTER_COMPLETE)
                    ;
                md->writer_flag [0] = NO_FLAG;

                uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                flag [0] = START_WRITES;
                pthread_mutex_lock (&md->coordinator_mutex);
                queue_enqueue (&md->coordinator_flag, flag);
                pthread_mutex_unlock (&md->coordinator_mutex);
            }

            // wait to be told to start writing
            if (md->rank != md->sub_coord_rank)
            {
                uint64_t msgx [PARAMETER_COUNT];
                pthread_mutex_lock (&md->mutex);
                MPI_Recv (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                          ,md->sub_coord_rank, TAG_WRITER
                          ,md->group_comm, &status
                          );

                pthread_mutex_unlock (&md->mutex);
                COPY_ALL_PARAMS(msg,msgx);
            }
            else
            {
                while (md->writer_flag [0] != DO_WRITE_FLAG)
                    ;
                COPY_ALL_PARAMS(msg,md->writer_flag);
                md->writer_flag [0] = NO_FLAG;
            }

            // set our base offset for building the index
            fd->base_offset = msg [5] * md->stripe_size;
            // build index appending to any existing index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                 ,&md->old_attrs_root
                                 );

            // we need the size of the buffer for responding to the write
            adios_write_index_v1 (&index_buffer, &index_buffer_size
                                 ,&index_buffer_offset
                                 ,0, md->old_pg_root
                                 ,md->old_vars_root
                                 ,md->old_attrs_root
                                 );

            if (msg [1] == msg [3]) // same file
            {
                if (md->f == -1) printf ("we got a bad file handle\n");
                lseek (md->f, md->stripe_size * msg [5], SEEK_SET);
                ssize_t s = write (md->f, fd->buffer, fd->bytes_written);
                if (s != fd->bytes_written)
                {
                    fprintf (stderr, "Need to do multi-write 1 (tried: %llu wrote: %lld) errno %d\n", fd->bytes_written, s, errno);
                }
            }
            else // we are adaptive writing and need to use the other file
            {
                MPI_File fh;
                char * new_name;

                new_name = malloc (  strlen (method->base_path)
                                   + strlen (fd->name) + 1 + 6
                                  ); // 6 extra for '.XXXXX' file number
                char split_format [10] = "%s%s.%d";
                sprintf (new_name, split_format, method->base_path
                        ,fd->name, md->group
                        );

                int f = open (new_name, O_WRONLY | O_LARGEFILE);
                if (f != -1)
                {
                    lseek (f, md->stripe_size * msg [5], SEEK_SET);
                    ssize_t s = write (f, fd->buffer, fd->bytes_written);
                    if (s != fd->bytes_written)
                    {
                        fprintf (stderr, "Need to do multi-write 2\n");
                    }
                    close (f);
                }
                else
                {
                    fprintf (stderr, "ADAPTIVE WRITE FAILURE. File: %s\n", new_name);
                }
            }

            // respond to sub coord(s) we are done
            uint64_t new_offset = msg [5] + 1;
            int source = msg [2]; // who to tell our index to
            if (md->rank != msg [2])
            {
                uint64_t msgx [PARAMETER_COUNT];
                msgx [0] = WRITE_COMPLETE;
                msgx [1] = msg [1];
                msgx [2] = msg [3];
                msgx [3] = new_offset;
                msgx [4] = index_buffer_offset;
                pthread_mutex_lock (&md->mutex);
                MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                          ,msg [2], TAG_SUB_COORDINATOR
                          ,md->group_comm, &req
                          );

                pthread_mutex_unlock (&md->mutex);
            }
            else
            {
                uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                INIT_PARAMS(flag);
                flag [0] = WRITE_COMPLETE;
                flag [1] = msg [1];
                flag [2] = msg [3];
                flag [3] = new_offset;
                flag [4] = index_buffer_offset;
                pthread_mutex_lock (&md->sub_coordinator_mutex);
                queue_enqueue (&md->sub_coordinator_flag, flag);
                pthread_mutex_unlock (&md->sub_coordinator_mutex);
            }

            // tell our sub_coord if we are adaptive writing
            if (msg [1] != msg [3])
            {
                if (md->rank != msg [4])
                {
                    uint64_t msgx [PARAMETER_COUNT];
                    msgx [0] = WRITE_COMPLETE;
                    msgx [1] = msg [1];
                    msgx [2] = msg [3];
                    msgx [3] = new_offset;
                    msgx [4] = index_buffer_offset;
                    pthread_mutex_lock (&md->mutex);
                    MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                              ,msg [4], TAG_SUB_COORDINATOR
                              ,md->group_comm, &req
                              );

                    pthread_mutex_unlock (&md->mutex);
                }
                else
                {
                    uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                    INIT_PARAMS(flag);
                    flag [0] = WRITE_COMPLETE;
                    flag [1] = msg [1];
                    flag [2] = msg [3];
                    flag [3] = new_offset;
                    flag [4] = index_buffer_offset;
                    pthread_mutex_lock (&md->sub_coordinator_mutex);
                    queue_enqueue (&md->sub_coordinator_flag, flag);
                    pthread_mutex_unlock (&md->sub_coordinator_mutex);
                }
            }

            // wait to do index
            if (md->rank != source)
            {
                uint64_t msgx [PARAMETER_COUNT];
                int message_available = 0;
                while (!message_available)
                    MPI_Iprobe (MPI_ANY_SOURCE, TAG_WRITER
                               ,md->group_comm
                               ,&message_available, &status
                               );

                pthread_mutex_lock (&md->mutex);
                MPI_Recv (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                          ,source, TAG_WRITER
                          ,md->group_comm, &status
                          );
                assert (msgx [0] == SEND_INDEX);

                pthread_mutex_unlock (&md->mutex);
                msg [0] = msgx [0];
                // do not get the rest of the parameters because
                // we need the old ones
            }
            else
            {
                while (md->writer_flag [0] == NO_FLAG)
                    ;
                msg [0] = md->writer_flag [0];
                assert (msg [0] == SEND_INDEX);
                // do not get the rest of the parameters because
                // we need the old ones
                md->writer_flag [0] = NO_FLAG;
            }

            // send index
            if (md->rank != msg [2])
            {
                pthread_mutex_lock (&md->mutex);
                MPI_Isend (index_buffer, index_buffer_offset, MPI_BYTE, msg [2]
                          ,TAG_SUB_COORDINATOR, md->group_comm, &req
                          );

                pthread_mutex_unlock (&md->mutex);
            }
            else
            {
                uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                INIT_PARAMS(flag);
                flag [0] = INDEX_BODY;
                flag [1] = (uint64_t) index_buffer;
                pthread_mutex_lock (&md->sub_coordinator_mutex);
                queue_enqueue (&md->sub_coordinator_flag, flag);
                pthread_mutex_unlock (&md->sub_coordinator_mutex);
            }

            // finished with output. Shutdown the system
            if (md->rank == md->sub_coord_rank)
            {
                uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                INIT_PARAMS(flag);
                flag [0] = SHUTDOWN_FLAG;
                pthread_mutex_lock (&md->sub_coordinator_mutex);
                queue_enqueue (&md->sub_coordinator_flag, flag);
                pthread_mutex_unlock (&md->sub_coordinator_mutex);

                err = pthread_join (md->sub_coordinator, NULL);
                if (err != 0) printf ("join sub coord error: %d\n", err);
            }
            if (md->rank == md->coord_rank)
            {
                uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                INIT_PARAMS(flag);
                flag [0] = SHUTDOWN_FLAG;
                pthread_mutex_lock (&md->coordinator_mutex);
                queue_enqueue (&md->coordinator_flag, flag);
                pthread_mutex_unlock (&md->coordinator_mutex);

                err = pthread_join (md->coordinator, NULL);
                if (err != 0) printf ("join sub coord error: %d\n", err);
            }

            if (buffer)
            {
                free (buffer);
                buffer = 0;
                buffer_size = 0;
                buffer_offset = 0;
            }

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
            adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                                 ,md->old_attrs_root
                                 );
            new_pg_root = 0;
            new_vars_root = 0;
            new_attrs_root = 0;
            md->old_pg_root = 0;
            md->old_vars_root = 0;
            md->old_attrs_root = 0;

            if (index_buffer)
                free (index_buffer);

            break;
        }

        case adios_mode_append:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset;

            if (fd->shared_buffer == adios_flag_no)
            {
                MPI_Offset new_off;
                // set it up so that it will start at 0, but have correct sizes
                MPI_File_get_position (md->fh, &new_off);
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_vars_v1 (fd);
                // fd->vars_start gets updated with the size written
                MPI_File_seek (md->fh, md->vars_start, MPI_SEEK_SET);
                MPI_File_write (md->fh, fd->buffer, md->vars_header_size
                               ,MPI_BYTE, &md->status
                               );
                int count;
                MPI_Get_count (&md->status, MPI_BYTE, &count);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "d:MPI method tried to write %llu, "
                                     "only wrote %d\n"
                            ,md->vars_header_size
                            ,count
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
                adios_shared_buffer_free (&md->b);

                adios_write_open_attributes_v1 (fd);
                md->vars_start = new_off;
                md->vars_header_size = fd->offset;
                MPI_File_seek (md->fh, new_off + md->vars_header_size
                              ,MPI_SEEK_SET
                              ); // go back to end, but after attr header
                fd->base_offset += fd->offset;  // add size of header
                fd->offset = 0;
                fd->bytes_written = 0;

                while (a)
                {
                    adios_write_attribute_v1 (fd, a);
                    MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                                   ,MPI_BYTE, &md->status
                                   );
                    MPI_Get_count (&md->status, MPI_BYTE, &count);
                    if (count != fd->bytes_written)
                    {
                        fprintf (stderr, "e:MPI method tried to write %llu, "
                                         "only wrote %d\n"
                                ,fd->bytes_written
                                ,count
                                );
                    }
                    fd->base_offset += count;
                    fd->offset = 0;
                    fd->bytes_written = 0;
                    adios_shared_buffer_free (&md->b);

                    a = a->next;
                }

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);
                MPI_File_seek (md->fh, md->vars_start, MPI_SEEK_SET);
                // fd->vars_start gets updated with the size written
                MPI_File_write (md->fh, fd->buffer, md->vars_header_size
                               ,MPI_BYTE, &md->status
                               );
                MPI_Get_count (&md->status, MPI_BYTE, &count);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "f:MPI method tried to write %llu, "
                                     "only wrote %d\n"
                            ,md->vars_header_size
                            ,count
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
            }

            // build index appending to any existing index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                 ,&md->old_attrs_root
                                 );
            // if collective, gather the indexes from the rest and call
            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    int * index_sizes = malloc (4 * md->size);
                    int * index_offsets = malloc (4 * md->size);
                    char * recv_buffer = 0;
                    uint32_t size = 0;
                    uint32_t total_size = 0;
                    int i;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->group_comm
                               );

                    for (i = 0; i < md->size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->group_comm
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    for (i = 1; i < md->size; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }
                    md->b.buff = buffer_save;
                    md->b.length = buffer_size_save;
                    md->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);
                }
                else
                {
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
                                         );

                    MPI_Gather (&buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->group_comm
                               );
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->group_comm
                                );
                }
            }

            if (fd->shared_buffer == adios_flag_yes)
            {
                // everyone writes their data
                MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
                MPI_File_write (md->fh, fd->buffer, fd->bytes_written, MPI_BYTE
                               ,&md->status
                               );
            }

            if (md->rank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,index_start, md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
                                     );
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                MPI_File_seek (md->fh, md->b.pg_index_offset, MPI_SEEK_SET);
                MPI_File_write (md->fh, buffer, buffer_offset, MPI_BYTE
                               ,&md->status
                               );
            }

            free (buffer);

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
            adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                                 ,md->old_attrs_root
                                 );
            new_pg_root = 0;
            new_vars_root = 0;
            new_attrs_root = 0;
            md->old_pg_root = 0;
            md->old_vars_root = 0;
            md->old_attrs_root = 0;

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);
        }
    }

    if (md && md->f)
    {
        close (md->f);
        md->f = -1;
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {
        md->group_comm = MPI_COMM_NULL;
    }

    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->group_comm = MPI_COMM_NULL;

    adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                         ,md->old_attrs_root
                         );
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
}

void adios_adaptive_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_adaptive_data_struct * md =
                  (struct adios_adaptive_data_struct *) method->method_data;

    if (adios_adaptive_initialized)
    {
        adios_adaptive_initialized = 0;
        pthread_mutex_destroy (&md->mutex);
        // these need a cast on the free because they are volatile
        if (md->writer_flag)
            free ((void *) md->writer_flag);
#if 0
        if (md->w_sub_coordinator_flag)
            free ((void *) md->w_sub_coordinator_flag);
        if (md->c_sub_coordinator_flag)
            free ((void *) md->c_sub_coordinator_flag);
        if (md->c_coordinator_flag)
            free ((void *) md->c_coordinator_flag);
        if (md->w_coordinator_flag)
            free ((void *) md->w_coordinator_flag);
#else
        queue_destroy (&md->sub_coordinator_flag);
        queue_destroy (&md->coordinator_flag);
        pthread_mutex_destroy (&md->sub_coordinator_mutex);
        pthread_mutex_destroy (&md->coordinator_mutex);
#endif
    }
}

void adios_adaptive_end_iteration (struct adios_method_struct * method)
{
}

void adios_adaptive_start_calculation (struct adios_method_struct * method)
{
}

void adios_adaptive_stop_calculation (struct adios_method_struct * method)
{
}

// each sub coordinator will do the following:
// 1. wait for a registration indication from all writers it manages
// 2. send each writer, one at a time, a note to write and wait for
//    a completion notice.
// 3. Possibly process coordinator messages to tell a writer to go to a
//    different file. Wait for response for writer and then respond to
//    coordinator
// 4. When done tell coordinator and wait for shutdown message
// 5. shutdown
static void * sub_coordinator_main (void * param)
{
    struct adios_adaptive_data_struct * md = (struct adios_adaptive_data_struct *) param;

    MPI_Request req;
    MPI_Status status;
    int i;
    uint64_t msg [PARAMETER_COUNT];
    int message_available;
    int source;

    // register the sub_coordinator with the coordinator
    if (md->rank != md->coord_rank)
    {
        uint64_t msgx [PARAMETER_COUNT];
        msgx [0] = REGISTER_FLAG;
        msgx [1] = md->group;
        msgx [2] = md->rank;
        pthread_mutex_lock (&md->mutex);
        MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG, md->coord_rank
                  ,TAG_COORDINATOR, md->group_comm, &req
                  );

        pthread_mutex_unlock (&md->mutex);
    }
    else
    {
        uint64_t * flag = malloc (8 * PARAMETER_COUNT);
        INIT_PARAMS(flag);
        flag [0] = REGISTER_FLAG;
        flag [1] = md->group;
        flag [2] = md->rank;
        pthread_mutex_lock (&md->coordinator_mutex);
        queue_enqueue (&md->coordinator_flag, flag);
        pthread_mutex_unlock (&md->coordinator_mutex);
    }

    // get the registration from the writers
    i = 0;
    while (i < md->group_size)
    {
        MPI_Iprobe (MPI_ANY_SOURCE, TAG_SUB_COORDINATOR, md->group_comm
                   ,&message_available, &status
                   );
        if (message_available)
        {
            uint64_t msgx [PARAMETER_COUNT];
            MPI_Recv (msgx, PARAMETER_COUNT, MPI_LONG_LONG, status.MPI_SOURCE
                     ,TAG_SUB_COORDINATOR
                     ,md->group_comm, &status
                     );

            COPY_ALL_PARAMS(msg,msgx);
            source = status.MPI_SOURCE;
            i++;
        }
        else
        {
            if (queue_size (&md->sub_coordinator_flag) != 0)
            {
                message_available = 1;
                pthread_mutex_lock (&md->sub_coordinator_mutex);
                uint64_t * flag;
                queue_dequeue (&md->sub_coordinator_flag, &flag);
                assert (flag [0] == REGISTER_FLAG);
                pthread_mutex_unlock (&md->sub_coordinator_mutex);
                COPY_ALL_PARAMS(msg,flag);
                free (flag);
                source = md->rank;
                i++;
            }
        }

        if (message_available)
        {
#if PRINT_MESSAGES
           printf ("B:writer registered group: %d writer: %d msg: %s\n"
                  ,md->group, source, message_to_string_full (msg)
                  );
#endif
            msg [0] = NO_FLAG;
            message_available = 0;
        }
    }

    // wait for START_WRITES
    if (md->rank != md->coord_rank)
    {
        uint64_t msgx [PARAMETER_COUNT];
        pthread_mutex_lock (&md->mutex);
        MPI_Recv (msgx, PARAMETER_COUNT, MPI_LONG_LONG, md->coord_rank
                  ,TAG_SUB_COORDINATOR
                  ,md->group_comm, &status
                  );

        pthread_mutex_unlock (&md->mutex);
        assert (msgx [0] == START_WRITES);
    }
    else
    {
        while (queue_size (&md->sub_coordinator_flag) == 0)
            ;
        pthread_mutex_lock (&md->sub_coordinator_mutex);
        uint64_t * flag;
        queue_dequeue (&md->sub_coordinator_flag, &flag);
        pthread_mutex_unlock (&md->sub_coordinator_mutex);
        assert (flag [0] == START_WRITES);
        free (flag);
    }

    // do writing
    int completed_writing = 0;  // track when we have notified coordinator we
                                // have completed (avoid multiple notifies)
    int writers_size = (int) (md->group_size * 1.20); // add 20% for adaptation
    int * writers = (int *) malloc (writers_size * sizeof (int));
    int current_writer = md->sub_coord_rank - md->group_size + 1;  // start at first so that we can always write using next_writer
    int next_writer = current_writer + 1; // track the next one so the adaptive
                                      // writer knows who to start next
    int writers_served = 0;
    int active_writers = 0; // track how many of our procs are writing

    int * adaptive_writers = 0;       // keep track of adaptive ranks pending
    int adaptive_writers_size = 0;    // track size of array
    int adaptive_writers_being_served = 0;  // how many pending

    uint64_t current_offset = 0;

    int * index_sizes = (int *) malloc (sizeof (int) * writers_served);
    int largest_index = 0;

    // for encoding the local index for writing and sending to the coord
    char * buffer = 0;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;

int xxx = 0;
    int currently_writing = 0;
    do
    {
if (!(xxx++ % 10000000)) printf ("AAAA %d\n", md->group);
        MPI_Iprobe (MPI_ANY_SOURCE, TAG_SUB_COORDINATOR, md->group_comm
                   ,&message_available, &status
                   );
        if (message_available)
        {
            uint64_t msgx [PARAMETER_COUNT];
            pthread_mutex_lock (&md->mutex);
            MPI_Recv (msgx, PARAMETER_COUNT, MPI_LONG_LONG, status.MPI_SOURCE
                     ,TAG_SUB_COORDINATOR
                     ,md->group_comm, &status
                     );
            pthread_mutex_unlock (&md->mutex);

            source = status.MPI_SOURCE;
            COPY_ALL_PARAMS(msg,msgx);
        }
        else
        {
            if (queue_size (&md->sub_coordinator_flag) != 0)
            {
                message_available = 1;
                uint64_t * flag;
                pthread_mutex_lock (&md->sub_coordinator_mutex);
                queue_dequeue (&md->sub_coordinator_flag, &flag);
                pthread_mutex_unlock (&md->sub_coordinator_mutex);
                COPY_ALL_PARAMS(msg,flag);
                free (flag);
                source = md->rank;
            }
        }

        if (message_available)
        {
#if PRINT_MESSAGES
            printf ("sc: %d source: %2d %s A\n", md->group, source, message_to_string_full (msg));
#endif
            message_available = 0;
            switch (msg [0])
            {
                case WRITE_COMPLETE:
                {
                    // if this group was writing it, we were tracking it
                    if (msg [2] == md->group)
                    {
                        active_writers--;
printf ("sc: %d active writers remaining: %d\n", md->group, active_writers);
                    }

                    // if the target group is our group (we wrote to our file)
                    if (msg [1] == md->group)
                    {
                        // save the writer in our list for indexing
                        if (writers_size <= writers_served)
                        {
                            writers_size = (int) (writers_size * 1.20);
                            writers = (int *) realloc (writers,   writers_size
                                                                * sizeof (int)
                                                      );
                        }
                        writers [writers_served] = source;
                        index_sizes [writers_served] = msg [4];
                        if (index_sizes [writers_served] > largest_index)
                            largest_index = index_sizes [writers_served];

                        writers_served++;
printf ("sc: %d saved writer for indexing: %lld\n", md->group, source);
                    }

                    // if it is adaptive, tell the coordinator it is done
                    if (msg [2] == md->group && msg [1] != msg [2])
                    {
                         printf ("What do we do here? group: %d msg [1] %lld msg [2] %lld msg [3] %lld\n", md->group, msg [1], msg [2], msg [3]);
                        // tell coordinator done with this one
                        // and if we have more capacity
                        // only if we were writing to our group
                        //if (msg [1] != msg [2])
                        {
                            if (md->rank != md->coord_rank)
                            {
printf ("sending write_complete: 1 g: %d 1: %lld 2: %lld 3: %lld\n", md->group, msg [1], msg [2], msg [3]);
                                uint64_t msgx [PARAMETER_COUNT];
                                msgx [0] = WRITE_COMPLETE;
                                msgx [1] = msg [1];
                                msgx [2] = msg [2];
                                msgx [3] = msg [3];
                                pthread_mutex_lock (&md->mutex);
                                MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                                          ,md->coord_rank, TAG_COORDINATOR
                                          ,md->group_comm, &req
                                          );

                                pthread_mutex_unlock (&md->mutex);
                            }
                            else
                            {
printf ("sending write_complete: 2 g: %d 1: %lld 2: %lld 3: %lld\n", md->group, msg [1], msg [2], msg [3]);
                                uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                                INIT_PARAMS(flag);
                                flag [0] = WRITE_COMPLETE;
                                flag [1] = msg [1];
                                flag [2] = msg [2];
                                flag [3] = msg [3];
                                pthread_mutex_lock (&md->coordinator_mutex);
                                queue_enqueue (&md->coordinator_flag, flag);
                                pthread_mutex_unlock (&md->coordinator_mutex);
                            }
                        }
#if 0
                        else
                        {
                            printf ("we shouldn't ever get here\n");
                        }
#endif
                    }

                    // if it was local, start the next write or send
                    // the write complete for this group (at bottom)
                    if (msg [1] == msg [2])
                    {
                        assert (msg [1] == md->group);
                        currently_writing = 0;
                        current_writer = next_writer++;
                    }
                    break;
                }

                case ADAPTIVE_WRITE_START:
                {
printf ("next_writer %d md->rank %d\n", next_writer, md->rank);
                    if (next_writer > md->rank)
                    {
                        // tell coordinator we are done and can't do it
                        if (md->rank != md->coord_rank)
                        {
                            uint64_t msgx [PARAMETER_COUNT];
                            msgx [0] = WRITERS_BUSY;
                            msgx [1] = msg [1];
                            msgx [2] = md->group;
                            pthread_mutex_lock (&md->mutex);
                            MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                                      ,md->coord_rank, TAG_COORDINATOR
                                      ,md->group_comm, &req
                                      );

                            pthread_mutex_unlock (&md->mutex);
                        }
                        else
                        {
                            uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                            INIT_PARAMS(flag);
                            flag [0] = WRITERS_BUSY;
                            flag [1] = msg [1];
                            flag [2] = md->group;
                            pthread_mutex_lock (&md->coordinator_mutex);
                            queue_enqueue (&md->coordinator_flag, flag);
                            pthread_mutex_unlock (&md->coordinator_mutex);
                        }
                    }
                    else
                    {
#if 0
printf ("AA\n");
                    if (adaptive_writers_size <= adaptive_writers_being_served + 1)
                    {
printf ("BB\n");
                        adaptive_writers_size += 10;
                        adaptive_writers = (int *) realloc (adaptive_writers
                                                           ,adaptive_writers_size
                                                           );
                    }
                    adaptive_writers [adaptive_writers_being_served++] = next_writer;
#endif
                    if (next_writer < md->rank)
                    {
printf ("CC\n");
                        active_writers++;
                        uint64_t msgx [PARAMETER_COUNT];
                        msgx [0] = DO_WRITE_FLAG;
                        msgx [1] = msg [1];
                        msgx [2] = msg [2];
                        msgx [3] = md->group;
                        msgx [4] = md->rank;
                        msgx [5] = msg [3];
                        pthread_mutex_lock (&md->mutex);
                        MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                                  ,next_writer, TAG_WRITER
                                  ,md->group_comm, &req
                                  );

                        pthread_mutex_unlock (&md->mutex);
                        next_writer++;
                    }
                    else
                    {
printf ("DD\n");
                        if (next_writer == md->rank)
                        {
printf ("EE\n");
                            active_writers++;
                            while (md->writer_flag [0] != NO_FLAG)
                                ;
                            md->writer_flag [1] = msg [1];
                            md->writer_flag [2] = msg [2];
                            md->writer_flag [3] = md->group;
                            md->writer_flag [4] = md->rank;
                            md->writer_flag [5] = msg [3];
                            md->writer_flag [0] = DO_WRITE_FLAG;
                            next_writer++;
                        }
                    }
printf ("FF\n");
                    }
                    break;
                }

                case OVERALL_WRITE_COMPLETE:
                {
                    // tell writers to enter send index mode
                    for (i = 0; i < writers_served; i++)
                    {
                        if (writers [i] != md->rank)
                        {
                            uint64_t msgx [PARAMETER_COUNT];
                            msgx [0] = SEND_INDEX;
                            pthread_mutex_lock (&md->mutex);
                            MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                                      ,writers [i], TAG_WRITER
                                      ,md->group_comm, &req
                                      );

                            pthread_mutex_unlock (&md->mutex);
                        }
                        else
                        {
                            while (md->writer_flag [0] != NO_FLAG)
                                ;
                            md->writer_flag [0] = SEND_INDEX;
                        }
                    }

                    char * buf = malloc (largest_index + 1);
                    buf [largest_index] = 0;
                    struct adios_bp_buffer_struct_v1 b;
                    struct adios_index_process_group_struct_v1 * new_pg_root;
                    struct adios_index_var_struct_v1 * new_vars_root;
                    struct adios_index_attribute_struct_v1 * new_attrs_root;
                    new_pg_root = 0;
                    new_vars_root = 0;
                    new_attrs_root = 0;
                    for (i = 0; i < writers_served; i++)
                    {
                        buf [index_sizes [i]] = 0;
                        if (writers [i] != md->rank)
                        {
                            int message_available = 0;
                            while (!message_available)
                                MPI_Iprobe (MPI_ANY_SOURCE, TAG_SUB_COORDINATOR
                                           ,md->group_comm
                                           ,&message_available, &status
                                           );

                            pthread_mutex_lock (&md->mutex);
                            MPI_Recv (buf, index_sizes [i], MPI_BYTE
                                      ,writers [i]
                                      ,TAG_SUB_COORDINATOR, md->group_comm
                                      ,&status
                                      );

                            pthread_mutex_unlock (&md->mutex);
                            b.buff = buf;
                        }
                        else
                        {
                            uint64_t * flag;
                            while (queue_size (&md->sub_coordinator_flag) == 0)
                                ;
                            pthread_mutex_lock (&md->sub_coordinator_mutex);
                            queue_dequeue (&md->sub_coordinator_flag, &flag);
                            b.buff = (char *) (flag [1]);
                            pthread_mutex_unlock (&md->sub_coordinator_mutex);
                            free (flag);
                        }

                        // merge buf into the index
                        b.length = index_sizes [i];
                        b.offset = 0;

                        adios_parse_process_group_index_v1 (&b ,&new_pg_root);
                        adios_parse_vars_index_v1 (&b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&b ,&new_attrs_root);
                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }
                    free (buf);

                    uint64_t only_index_buffer_offset;
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
                                         );
                    only_index_buffer_offset = buffer_offset;
                    adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                    lseek (md->f, md->stripe_size * msg [1], SEEK_SET);
                    ssize_t s = write (md->f, buffer, buffer_offset);
                    if (s != buffer_offset)
                    {
                        fprintf (stderr, "Need to do multi-write 3\n");
                    }

                    uint64_t index_size = buffer_offset;

                    // send index to the coordinator for global use
                    if (md->rank != md->coord_rank)
                    {
                        uint64_t msgx [PARAMETER_COUNT];
                        msgx [0] = INDEX_SIZE;
                        msgx [1] = md->group;
                        msgx [2] = only_index_buffer_offset;
                        pthread_mutex_lock (&md->mutex);
                        MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                                  ,md->coord_rank, TAG_COORDINATOR
                                  ,md->group_comm, &req
                                  );

                        pthread_mutex_unlock (&md->mutex);
                        MPI_Wait (&req, &status); // need to wait so receive
                                                  // is setup
                        pthread_mutex_lock (&md->mutex);
                        MPI_Isend (buffer, only_index_buffer_offset, MPI_BYTE
                                  ,md->coord_rank, TAG_COORDINATOR
                                  ,md->group_comm, &req
                                  );

                        pthread_mutex_unlock (&md->mutex);
                    }
                    else
                    {
                        uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                        INIT_PARAMS(flag);
                        flag [0] = INDEX_SIZE;
                        flag [1] = md->group;
                        flag [2] = only_index_buffer_offset;
                        flag [3] = (uint64_t) buffer;
                        pthread_mutex_lock (&md->coordinator_mutex);
                        queue_enqueue (&md->coordinator_flag, flag);
                        pthread_mutex_unlock (&md->coordinator_mutex);
                    }
                    free (buffer);

                    break;
                }
                case SHUTDOWN_FLAG:
                {
                    break;
                }
            }
#if PRINT_MESSAGES
            printf ("sc: %d source: %2d %s B\n", md->group, source, message_to_string_full (msg));
#endif
        }
        if (!currently_writing)
        {
            if (current_writer < md->rank)
            {
                active_writers++;
                uint64_t msgx [PARAMETER_COUNT];
                msgx [0] = DO_WRITE_FLAG;
                msgx [1] = md->group;
                msgx [2] = md->rank;
                msgx [3] = md->group;
                msgx [4] = md->rank;
                msgx [5] = current_offset++;
                pthread_mutex_lock (&md->mutex);
                MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                          ,current_writer, TAG_WRITER
                          ,md->group_comm, &req
                          );

                pthread_mutex_unlock (&md->mutex);
                currently_writing = 1;
            }
            else
            {
                if (current_writer == md->rank)
                {
                    active_writers++;
                    while (md->writer_flag [0] != NO_FLAG)
                        ;
                    md->writer_flag [1] = md->group; //msg [1];
                    md->writer_flag [2] = md->rank;
                    md->writer_flag [3] = md->group; //msg [1];
                    md->writer_flag [4] = md->rank;
                    md->writer_flag [5] = current_offset++;
                    md->writer_flag [0] = DO_WRITE_FLAG;
                    currently_writing = 1;
                }
                else
                {
                    if (!completed_writing && !active_writers)
                    {
                        completed_writing = 1;
                        if (md->rank != md->coord_rank)
                        {
printf ("sending write_complete: 3 g: %d 1: %lld 2: %lld 3: %lld\n", md->group, msg [1], msg [2], msg [3]);
                            uint64_t msgx [PARAMETER_COUNT];
                            msgx [0] = WRITE_COMPLETE;
                            msgx [1] = md->group;//msg [1];
                            msgx [2] = md->group;//msg [2];
                            msgx [3] = msg [3];
                            pthread_mutex_lock (&md->mutex);
                            MPI_Isend (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                                      ,md->coord_rank, TAG_COORDINATOR
                                      ,md->group_comm, &req
                                      );

                            pthread_mutex_unlock (&md->mutex);
                        }
                        else
                        {
printf ("sending write_complete: 4 g: %d 1: %lld 2: %lld 3: %lld\n", md->group, msg [1], msg [2], msg [3]);
                            uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                            INIT_PARAMS(flag);
                            flag [0] = WRITE_COMPLETE;
                            flag [1] = md->group; //msg [1];
                            flag [2] = md->group; //msg [2];
                            flag [3] = msg [3];
                            pthread_mutex_lock (&md->coordinator_mutex);
                            queue_enqueue (&md->coordinator_flag, flag);
                            pthread_mutex_unlock (&md->coordinator_mutex);
                        }
                    }
                }
            }
#if 0
            if (currently_writing)
            {
                current_writer = next_writer++;
            }
#endif
        }
    } while (msg [0] != SHUTDOWN_FLAG);
    free (writers);
    free (index_sizes);

    pthread_exit (NULL);
    return NULL;
}

// one coordinator overall that only talks with the sub_coordinators
// 1. wait for registration from each of the sub_coordinators.
// 2. wait for messages from sub_coordinators indicating completion
// 3. Whenever a subcoordinator indicates it is done, remove it from the
//    pending set and add it to the adaptation set.
// 4. Go through the pending set telling one to tell it to use on of the
//    adaptation set. Wait for a response from that adaptive write.
// 5. When all sub_coordinators report completed, tell sub_coordinators to
//    build index, close, and send a copy of the index to the coordinator.
// 6. collect the sub indices and make the overall index file.
// 7. wait for shutdown message.
// 8. shutdown.
static void * coordinator_main (void * param)
{
    struct adios_adaptive_data_struct * md = (struct adios_adaptive_data_struct *) param;

    uint64_t msg [PARAMETER_COUNT];
    int source;
    int message_available;
    MPI_Status status;
    MPI_Request req [md->groups]; // use a different request for each sub coord
    int i;
    enum group_state
    {
         STATE_WRITING = 0
        ,STATE_COMPLETE = 1
        ,STATE_WRITERS_ALL_OCCUPIED = 2
    };

    int sub_coord_ranks [md->groups];

    char group_state [md->groups];
    int groups_complete = 0;
    uint64_t group_offset [md->groups];
    //uint64_t index_size [md->groups];  // I don't think we need this anymore

    memset (group_offset, 0, md->groups * sizeof (uint64_t));
    memset (group_state, STATE_WRITING, md->groups);
    // process sub_coordinator registrations
    i = 0;
    while (i < md->groups)
    {
        MPI_Iprobe (MPI_ANY_SOURCE, TAG_COORDINATOR, md->group_comm
                   ,&message_available, &status
                   );
        if (message_available)
        {
            uint64_t msgx [PARAMETER_COUNT];
printf ("X 1\n");
            pthread_mutex_lock (&md->mutex);
            MPI_Recv (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                     ,MPI_ANY_SOURCE, TAG_COORDINATOR
                     ,md->group_comm, &status
                     );
            pthread_mutex_unlock (&md->mutex);
printf ("X 2\n");
            COPY_ALL_PARAMS(msg,msgx);
            assert (msg [0] == REGISTER_FLAG);

            source = msgx [2];
            i++;
        }
        else
        {
            if (queue_size (&md->coordinator_flag) != 0)
            {
printf ("X 3\n");
                message_available = 1;
                uint64_t * flag;
                pthread_mutex_lock (&md->coordinator_mutex);
                queue_dequeue (&md->coordinator_flag, &flag);
                pthread_mutex_unlock (&md->coordinator_mutex);
                COPY_ALL_PARAMS(msg,flag);
                source = flag [2];
                assert (flag [0] == REGISTER_FLAG);
                assert (flag [2] == md->rank);
                free (flag);
                i++;
printf ("X 4\n");
            }
        }

        if (message_available)
        {
#if PRINT_MESSAGES
            printf ("B:sub_coordinator registered rank: %2d group: %d msg: %s\n"
                   ,source, msg [1], message_to_string_full (msg)
                   );
#endif
            sub_coord_ranks [msg [1]] = source;
            msg [0] = 0;
            message_available = 0;
        }
    }

    md->writer_flag [0] = REGISTER_COMPLETE;

    // wait for registration to complete before continuing
    while (queue_size (&md->coordinator_flag) == 0)
        ;
    uint64_t * flag;
    pthread_mutex_lock (&md->coordinator_mutex);
    queue_dequeue (&md->coordinator_flag, &flag);
    pthread_mutex_unlock (&md->coordinator_mutex);
    assert (flag [0] == START_WRITES);
    free (flag);

    for (i = 0; i < md->groups; i++)
    {
        if (sub_coord_ranks [i] != md->rank)
        {
            uint64_t msgx [PARAMETER_COUNT];
            msgx [0] = START_WRITES;
            pthread_mutex_lock (&md->mutex);
            MPI_Isend (msgx, PARAMETER_COUNT
                      ,MPI_LONG_LONG
                      ,sub_coord_ranks [i]
                      ,TAG_SUB_COORDINATOR
                      ,md->group_comm, &req [i]
                      );

            pthread_mutex_unlock (&md->mutex);
        }
        else
        {
            flag = malloc (8 * PARAMETER_COUNT);
            INIT_PARAMS(flag);
            flag [0] = START_WRITES;
            pthread_mutex_lock (&md->sub_coordinator_mutex);
            queue_enqueue (&md->sub_coordinator_flag, flag);
            pthread_mutex_unlock (&md->sub_coordinator_mutex);
        }
    }

    uint64_t largest_index_size = 0;
    int index_sizes_received = 0;
    char * index_buf = 0;     // for receiving from remote
    ssize_t index_size = 0;

    char * buffer = 0;          // for building overall
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;

    struct adios_file_index_format_v2
    {
        uint32_t file_number;
        uint16_t name_len;
        char * name;
        uint64_t offset;
        uint64_t length;
    };

    uint64_t file_index_count = 0;
    char start_index_collection = 0;
    char index_collection_started = 0;
    uint64_t adaptive_writes_outstanding = 0;

    struct adios_file_index_format_v2 * file_index =
       (struct adios_file_index_format_v2 *)
              malloc (sizeof (struct adios_file_index_format_v2) * md->groups);

    // process messages for the coordinator either on or off process
int xxx = 0;
    do
    {
if (!(xxx++ % 10000000)) printf ("BBBB adaptive_writes_outstanding: %lld\n", adaptive_writes_outstanding);
        MPI_Iprobe (MPI_ANY_SOURCE, TAG_COORDINATOR, md->group_comm
                   ,&message_available, &status
                   );
        if (message_available)
        {
            uint64_t msgx [PARAMETER_COUNT];
            pthread_mutex_lock (&md->mutex);
            MPI_Recv (msgx, PARAMETER_COUNT, MPI_LONG_LONG
                     ,status.MPI_SOURCE, TAG_COORDINATOR
                     ,md->group_comm, &status
                     );

            pthread_mutex_unlock (&md->mutex);
            source = status.MPI_SOURCE;
            COPY_ALL_PARAMS(msg,msgx);
        }
        else
        {
            if (queue_size (&md->coordinator_flag) != 0)
            {
                message_available = 1;
                uint64_t * flag;
                pthread_mutex_lock (&md->coordinator_mutex);
                queue_dequeue (&md->coordinator_flag, &flag);
                pthread_mutex_unlock (&md->coordinator_mutex);
                COPY_ALL_PARAMS(msg,flag);
                source = md->rank;
                free (flag);
            }
        }

        if (message_available)
        {
#if PRINT_MESSAGES
            printf ("c: source: %2d msg: %s A\n"
                   ,source, message_to_string_full (msg)
                   );
            assert (msg [0] != ADAPTIVE_WRITE_START);
#endif
            message_available = 0;
            switch (msg [0])
            {
                case WRITE_COMPLETE:
                {
                    // what part of the file finished writing last is unknown
                    // so only save the largest final offset
                    if (msg [3] > group_offset [msg [1]])
                        group_offset [msg [1]] = msg [3];

                    // we just finished this group, so
                    // start an adaptive write for this file
                    if (msg [1] == msg [2])
                    {
                        groups_complete++;
                        group_state [msg [1]] = STATE_COMPLETE;
                        if (groups_complete != md->groups)
                        {
                            int i = (msg [1] + 1) % md->groups;
                            while (i != msg [1])
                            {
                                if (group_state [i] == STATE_WRITING)
                                {
                                    // tell the subcoordinator to write
                                    // to this file
                                    adaptive_writes_outstanding++;
                                    if (sub_coord_ranks [i] != md->rank)
                                    {
printf ("new adaptive writes: %lld A1 ++ to: %d for: %d\n", adaptive_writes_outstanding, i, msg [1]);
                                        uint64_t msgx [PARAMETER_COUNT];
                                        msgx [0] = ADAPTIVE_WRITE_START;
                                        msgx [1] = msg [1];
                                        msgx [2] = sub_coord_ranks [msg [1]];
                                        msgx [3] = group_offset [msg [1]];
                                        pthread_mutex_lock (&md->mutex);
                                        MPI_Isend (msgx, PARAMETER_COUNT
                                                  ,MPI_LONG_LONG
                                                  ,sub_coord_ranks [i]
                                                  ,TAG_SUB_COORDINATOR
                                                  ,md->group_comm, &req [i]
                                                  );

                                        pthread_mutex_unlock (&md->mutex);
                                    }
                                    else
                                    {
printf ("new adaptive writes: %lld A2 ++ to: %d for: %d\n", adaptive_writes_outstanding, i, msg [1]);
                                        uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                                        INIT_PARAMS(flag);
                                        flag [0] = ADAPTIVE_WRITE_START;
                                        flag [1] = msg [1];
                                        flag [2] = sub_coord_ranks [msg [1]];
                                        flag [3] = group_offset [msg [1]];
                                        pthread_mutex_lock (&md->sub_coordinator_mutex);
                                        queue_enqueue (&md->sub_coordinator_flag, flag);
                                        pthread_mutex_unlock (&md->sub_coordinator_mutex);
                                    }
                                    break;
                                }
                                i = (i + 1) % md->groups;
                            }
                            if (i == msg [1]) // we didn't find someplace to write
                                printf ("c: ending adaptive write for group: %lld\n", msg [1]);
                        }
                        else // start index collection if no adaptive left
                        {
                            if (!adaptive_writes_outstanding)
                            {
printf ("START INDEX COLLECTION A\n");
                                start_index_collection = 1;
printf ("START INDEX COLLECTION B\n");
                            }
                        }
                    }
                    else  // move to the next adaptive writer
                    {
                        adaptive_writes_outstanding--;
printf ("new adaptive writes: %lld B -- msg [1] %lld msg [2] %lld\n", adaptive_writes_outstanding, msg [1], msg [2]);
                        int i = (msg [2] + 1) % md->groups;
                        while (i != msg [2])
                        {
                            if (group_state [i] == STATE_WRITING)
                            {
                                // tell the subcoordinator to write
                                // to this file
                                adaptive_writes_outstanding++;
                                if (sub_coord_ranks [i] != md->rank)
                                {
printf ("new adaptive writes: %lld C1 ++ to: %d for: %d\n", adaptive_writes_outstanding, i, msg [1]);
                                    uint64_t msgx [PARAMETER_COUNT];
                                    msgx [0] = ADAPTIVE_WRITE_START;
                                    msgx [1] = msg [1];
                                    msgx [2] = sub_coord_ranks [msg [1]];
                                    msgx [3] = group_offset [msg [1]];
                                    pthread_mutex_lock (&md->mutex);
                                    MPI_Isend (msgx, PARAMETER_COUNT
                                              ,MPI_LONG_LONG
                                              ,sub_coord_ranks [i]
                                              ,TAG_SUB_COORDINATOR
                                              ,md->group_comm, &req [i]
                                              );

                                    pthread_mutex_unlock (&md->mutex);
                                }
                                else
                                {
printf ("new adaptive writes: %lld C2 ++ to: %d for: %d\n", adaptive_writes_outstanding, i, msg [1]);
                                    uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                                    INIT_PARAMS(flag);
                                    flag [0] = ADAPTIVE_WRITE_START;
                                    flag [1] = msg [1];
                                    flag [2] = sub_coord_ranks [msg [1]];
                                    flag [3] = group_offset [msg [1]];
                                    pthread_mutex_lock (&md->sub_coordinator_mutex);
                                    queue_enqueue (&md->sub_coordinator_flag, flag);
                                    pthread_mutex_unlock (&md->sub_coordinator_mutex);
                                }

                                break;
                            }
                            i = (i + 1) % md->groups;
                        }
                        if (i == msg [2]) // we didn't find someplace to write
                            printf ("c: ending adaptive write for group: %lld\n", msg [1]);
printf ("i: %d msg [1] %lld wb adaptive_writes_outstanding: %d\n", i, msg [1], adaptive_writes_outstanding);
                        if (i == msg [1] && !adaptive_writes_outstanding && groups_complete == md->groups)
                        {
printf ("START INDEX COLLECTION Y\n");
                            start_index_collection = 1;
                        }
                    }
                    break;
                }

                case WRITERS_BUSY: // we told it to adaptive write, but
                                   // all are writing already for that group
                {
                    adaptive_writes_outstanding--;
printf ("new adaptive writes: %lld D -- msg [1] %lld msg [2] %lld\n", adaptive_writes_outstanding, msg [1], msg [2]);
                    // the writer that told us this can't take any more
                    // adaptive requests so mark as occupied (only if still
                    // marked as writing)
                    if (group_state [msg [2]] == STATE_WRITING)
                        group_state [msg [2]] = STATE_WRITERS_ALL_OCCUPIED;

                    // look for the next group we can ask to write
                    int i = (msg [2] + 1) % md->groups;
                    while (i != msg [2])  // look at all
                    {
                        if (   i != msg [1]  // don't write to the source group
                            && group_state [i] == STATE_WRITING
                           )
                        {
                            // tell the subcoordinator to write
                            // to this file
                            adaptive_writes_outstanding++;
                            if (sub_coord_ranks [i] != md->rank)
                            {
printf ("new adaptive writes: %lld E1 ++ to: %d for: %d\n", adaptive_writes_outstanding, i, msg [1]);
                                uint64_t msgx [PARAMETER_COUNT];
                                msgx [0] = ADAPTIVE_WRITE_START;
                                msgx [1] = msg [1];
                                msgx [2] = sub_coord_ranks [msg [1]];
                                msgx [3] = group_offset [msg [1]];
                                pthread_mutex_lock (&md->mutex);
                                MPI_Isend (msgx, PARAMETER_COUNT
                                          ,MPI_LONG_LONG
                                          ,sub_coord_ranks [i]
                                          ,TAG_SUB_COORDINATOR
                                          ,md->group_comm, &req [i]
                                          );

                                pthread_mutex_unlock (&md->mutex);
                                break;
                            }
                            else
                            {
printf ("new adaptive writes: %lld E2 ++ to: %d for: %d\n", adaptive_writes_outstanding, i, msg [1]);
                                uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                                INIT_PARAMS(flag);
                                flag [0] = ADAPTIVE_WRITE_START;
                                flag [1] = msg [1];
                                flag [2] = sub_coord_ranks [msg [1]];
                                flag [3] = group_offset [msg [1]];
                                pthread_mutex_lock (&md->sub_coordinator_mutex);
                                queue_enqueue (&md->sub_coordinator_flag, flag);
                                pthread_mutex_unlock (&md->sub_coordinator_mutex);
                                break;
                            }
                        }
                        i = (i + 1) % md->groups;
                    }
printf ("i: %d msg [1] %lld wb adaptive_writes_outstanding: %d\n", i, msg [1], adaptive_writes_outstanding);
                    if (i == msg [2] && !adaptive_writes_outstanding && groups_complete == md->groups)
                    {
printf ("START INDEX COLLECTION wb\n");
                        start_index_collection = 1;
                    }
                    break;
                }

                case INDEX_SIZE:
                {
                    int source_group = msg [1];
                    uint64_t proc_index_size = msg [2];

                    index_sizes_received++;

                    if (largest_index_size < proc_index_size)
                        largest_index_size = proc_index_size;

                    if (index_size < largest_index_size)
                    {
                        index_size = largest_index_size;
                        if (index_buf)
                            free (index_buf);
                        index_buf = malloc (largest_index_size + 1);
                    }

                    char * new_name;
                    int new_name_len;

                    new_name = malloc (  strlen (md->method->base_path)
                                       + strlen (md->fd->name) + 1 + 6
                                      );
                    char split_format [10] = "%s%s.%lld";
                    sprintf (new_name, split_format, md->method->base_path
                            ,md->fd->name, source_group
                            );
                    new_name_len = strlen (new_name);
                    uint64_t buffer_offset_tmp = buffer_offset;
                    buffer_offset += 8 + 4 + 2 + new_name_len;

                    if (source_group != md->group)
                    {
                        pthread_mutex_lock (&md->mutex);
                        MPI_Recv (index_buf, proc_index_size
                                  ,MPI_BYTE
                                  ,sub_coord_ranks [source_group]
                                  ,TAG_COORDINATOR
                                  ,md->group_comm, &status
                                  );
                        pthread_mutex_unlock (&md->mutex);
                        buffer_write (&buffer, &buffer_size, &buffer_offset
                                     ,index_buf, proc_index_size
                                     );
                    }
                    else
                    {
                        buffer_write (&buffer, &buffer_size, &buffer_offset
                                     ,(void *) msg [3]
                                     ,proc_index_size
                                     );
                    }

                    file_index [file_index_count].file_number = file_index_count + 1;
                    file_index [file_index_count].name_len = new_name_len;
                    file_index [file_index_count].name = new_name;
                    file_index [file_index_count].offset = buffer_offset_tmp;
                    file_index [file_index_count].length = proc_index_size
                                                         + new_name_len
                                                         + 2   // name len
                                                         + 4   // number
                                                         + 8;  // length
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].offset, 8
                                 );
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].file_number, 4
                                 );
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].name_len, 2
                                 );
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].name
                                 ,new_name_len
                                 );
                    file_index_count++;

                    // if we have recieved all, write to the file
                    if (index_sizes_received == md->groups)
                    {
                        adios_write_version_v2 (&buffer, &buffer_size, &buffer_offset);

                        char * new_name;
                        int old_mask;
                        int perm;

                        old_mask = umask (0);
                        umask (old_mask);
                        perm = old_mask ^ 0666;

                        new_name = malloc (  strlen (md->method->base_path)
                                           + strlen (md->fd->name) + 1
                                          );
                        char split_format [10] = "%s%s";
                        sprintf (new_name, split_format, md->method->base_path
                                ,md->fd->name
                                );
                        int f = open (new_name, O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC, perm);
                        if (f == -1) printf ("oops! %s\n", new_name);
                        ssize_t s = write (f, buffer, buffer_offset);
                        if (s != buffer_offset)
                        {
                            fprintf (stderr, "Need to do multi-write 4\n");
                        }
                        close (f);
                        free (new_name);

                        for (i = 0; i < md->groups; i++)
                            free (file_index [i].name);
                        free (file_index);
                    }
                    break;
                }
                case SHUTDOWN_FLAG:
                {
                    break;
                }
                default:
                {
                    printf ("Unknown coordinator message: %lld\n", msg [0]);
                    break;
                }
            }
            if (start_index_collection && !index_collection_started && !adaptive_writes_outstanding)
            {
                index_collection_started = 1;
                int i;
                for (i = 0; i < md->groups; i++)
                {
                    if (i != md->group)
                    {
                        uint64_t msgx [PARAMETER_COUNT];
                        msgx [0] = OVERALL_WRITE_COMPLETE;
                        msgx [1] = group_offset [i];
                        pthread_mutex_lock (&md->mutex);
                        MPI_Isend (msgx, PARAMETER_COUNT
                                  ,MPI_LONG_LONG
                                  ,sub_coord_ranks [i]
                                  ,TAG_SUB_COORDINATOR
                                  ,md->group_comm, &req [i]
                                  );

                        pthread_mutex_unlock (&md->mutex);
                    }
                    else
                    {
                        uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                        INIT_PARAMS(flag);
                        flag [0] = OVERALL_WRITE_COMPLETE;
                        flag [1] = group_offset [i];
                        pthread_mutex_lock (&md->sub_coordinator_mutex);
                        queue_enqueue (&md->sub_coordinator_flag, flag);
                        pthread_mutex_unlock (&md->sub_coordinator_mutex);
                    }
                }
            }
#if PRINT_MESSAGES
            printf ("c: source: %2d msg: %s B\n"
                   ,source, message_to_string_full (msg)
                   );
#endif
        }
    } while (msg [0] != SHUTDOWN_FLAG);

    //if (index_buf)
    //    free (index_buf);
    if (buffer)
        free (buffer);

    pthread_exit (NULL);
    return NULL;
}
