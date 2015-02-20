/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 *
 * Author: Qing Liu
 */


/******************************/
/* A read method for BP files */
/******************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>
#include "public/adios_types.h"
#include "public/adios_read.h"
#include "public/adios_error.h"
#include "core/bp_utils.h"
#include "core/bp_types.h"
#include "core/adios_read_hooks.h"
#include "core/futils.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define READ_CLOSE 0
#define MAX_DIMS 32
#define DEBUG 0
#define INVALID_VALUE -1
#define HT_SIZE 100

int adios_read_bp_staged1_init_method (MPI_Comm comm, PairStruct * param)
{
    return 0;
}

int adios_read_bp_staged1_finalize_method ()
{
    return 0;
}

ADIOS_FILE * adios_read_bp_staged1_open (const char * fname, MPI_Comm comm, enum ADIOS_LOCKMODE lock_mode, float timeout_sec)
{
    return 0;
}

ADIOS_FILE * adios_read_bp_staged1_open_file (const char * fname, MPI_Comm comm)
{
    return 0;
}

int adios_read_bp_staged1_close (ADIOS_FILE *fp)
{
    return 0;
}

int adios_read_bp_staged1_advance_step (ADIOS_FILE *fp, int last, float timeout_sec)
{
    return 0;
}

void adios_read_bp_staged1_release_step (ADIOS_FILE *fp)
{

}

ADIOS_VARINFO * adios_read_bp_staged1_inq_var_byid (const ADIOS_FILE *gp, int varid)
{
    return 0;
}

int adios_read_bp_staged1_inq_var_stat (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo, int per_step_stat, int per_block_stat)
{
    return 0;
}

int adios_read_bp_staged1_inq_var_blockinfo (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    return 0;
}

int adios_read_bp_staged1_schedule_read_byid (const ADIOS_FILE * fp, const ADIOS_SELECTION * sel, int varid, int from_steps, int nsteps, void * data)
{
    return 0;
}

int adios_read_bp_staged1_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    return 0;
}

int adios_read_bp_staged1_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    return 0;
}

int adios_read_bp_staged1_get_attr_byid (const ADIOS_FILE * fp, int attrid, enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    return 0;
}

int  adios_read_bp_staged1_get_dimension_order (const ADIOS_FILE *fp)
{
    return 0;
}

void adios_read_bp_staged1_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{

}

void adios_read_bp_staged1_get_groupinfo (const ADIOS_FILE *fp, int *ngroups, char ***group_namelist,  uint32_t **nvars_per_group,  uint32_t **nattrs_per_group)
{

}

int adios_read_bp_staged1_is_var_timed (const ADIOS_FILE *fp, int varid)
{
    int retval = 0;
    /*
        retval = this variable had time dimension at write
    */
    return retval;
}

#if 0
typedef struct read_info
{
    int ndim;
    uint64_t * start_notime;
    uint64_t * count_notime;
    int ndim_notime;
    int file_tdim;
    int start_time;
    int stop_time;
    uint64_t * dims;
} read_info;

typedef struct read_args
{
    int varid;
    int ndims;
    uint64_t * start;
    uint64_t * count;
    void * data;
    uint64_t size;
    int file_idx;
    uint64_t offset;
    void * parent;
} read_args;

typedef struct candidate_reader
{
    int rank;
    read_args * ra;
    struct candidate_reader * next;
} candidate_reader;

typedef struct bucket
{
    candidate_reader ** h;
    int size;
    int offset;
} bucket;

struct proc_struct
{
    ADIOS_GROUP * gp;
    int rank;
    int new_rank;
    int size;
    int groups;
    int group_size;
    int group;
    int n_total_sf;
    int n_my_sf;
    int f1;
    int f2;
    MPI_Comm group_comm;
    MPI_Comm new_comm;  /* MPI communicator for a group */
    MPI_Comm new_comm2; /* MPI communicator for all the aggregators */
    int aggregator_rank;
    int aggregator_new_rank;
    candidate_reader * local_read_request_list;
    candidate_reader * split_read_request_list1;
    candidate_reader * split_read_request_list2;
    bucket * ht;
    void * b;
    uint32_t num_aggregators;
    uint64_t chunk_size;
    int * aggregator_rank_array;
    int read_close_received;
    int group_close_received;
};

// used to specify which thread is the target for the MPI messages
enum MPI_TAG
{
     TAG_CONTROL = 0
    ,TAG_DATA = 1
};

static void swap_order (int n, uint64_t * array, int * tdim);
static void _swap_order (int n, uint64_t * array);
static int isTimeless (int tdim);
static void getReadInfo (ADIOS_GROUP * gp
                        ,int varid
                        ,uint64_t * start
                        ,uint64_t * count
                        ,read_info * ri
                        );
static void freeReadInfo (read_info * ri);
static int adios_read_bp_staged1_get_dimensioncharacteristics (
                           struct adios_index_characteristic_struct_v1 * c
                          ,uint64_t * ldims
                          ,uint64_t * gdims
                          ,uint64_t * offsets
                          );
static void adios_read_bp_staged1_get_dimensions (struct adios_index_var_struct_v1 * v
                                                 ,int ntsteps, int file_is_fortran
                                                 ,int * ndim, uint64_t ** dims
                                                 ,int * tdim
                                                 );
static struct adios_index_var_struct_v1 * adios_find_var_byid (ADIOS_GROUP * gp, int varid);
static int get_var_start_index (struct adios_index_var_struct_v1 * v, int t);
static int get_var_stop_index (struct adios_index_var_struct_v1 * v, int t);
static void getDataAddress (ADIOS_GROUP * gp, int varid
                    ,const uint64_t * start
                    ,const uint64_t * count
                    ,int * file_idx
                    ,uint64_t * offset
                    ,uint64_t * payload_size
                    );
static MPI_File * get_BP_file_handle(struct BP_file_handle * l, uint32_t file_index);
static void add_BP_file_handle (struct BP_file_handle ** l, struct BP_file_handle * n);
int64_t adios_read_bp_staged1_read_var_byid2 (ADIOS_GROUP * gp
                                            ,int varid
                                            ,const uint64_t * start
                                            ,const uint64_t * count
                                            ,void * data
                                            );
static void copy_buffer (ADIOS_GROUP * gp
                        ,candidate_reader * parent
                        ,candidate_reader * child
                        );
static ADIOS_VARINFO * _inq_var_byid (BP_FILE * fh, int varid);
void adios_read_bp_staged1_free_varinfo (ADIOS_VARINFO *vp);
static int get_num_subfiles (BP_FILE * fh);
static candidate_reader * parse_buffer (struct proc_struct * p, void * b, int len);

static void split_read_request (ADIOS_GROUP * gp
                               ,candidate_reader * r
                               ,candidate_reader ** h1
                               ,candidate_reader ** h2
                               );

static void build_hash (struct proc_struct * p);
static void free_hash (struct proc_struct * p);
int is_fortran_file (ADIOS_GROUP * gp);
int has_subfiles (ADIOS_GROUP * gp);
static void free_candidate_reader_list (candidate_reader * h);

static int isAggregator (struct proc_struct * p)
{
    return (p->rank == p->aggregator_rank);
}

static void list_print_readers (struct proc_struct * p, candidate_reader * h)
{
    read_args * ra;

    while (h)
    {
        ra = h->ra;
        printf ("%d:%d [varid: %d, ndims = %d, start[0]: %llu, start[1]: %llu, count[0]: %llu, count[1]: %llu, file_index: %d, offset = %llu]\n", p->rank, h->rank, ra->varid, ra->ndims, ra->start[0], ra->start[1], ra->count[0], ra->count[1], ra->file_idx, ra->offset);
        h = h->next;
    }
    printf ("\n");
}

static void list_print_my_readers (struct proc_struct * p, candidate_reader * h)
{
    read_args * ra;

    while (h)
    {   
        ra = h->ra;
        if (p->f1 <= ra->file_idx && ra->file_idx <= p->f2)
        {
            printf ("%d [varid: %d, ndims = %d, start[0]: %llu, start[1]: %llu, count[0]: %llu, count[1]: %llu, file_index: %d, offset = %llu]\n", p->rank, ra->varid, ra->ndims, ra->start[0], ra->start[1], ra->count[0], ra->count[1], ra->file_idx, ra->offset);
            printf ("%d [varid: %d, ndims = %d, start[0]: %llu, count[0]: %llu, file_index: %d, offset = %llu]\n", p->rank, ra->varid, ra->ndims, ra->start[0], ra->count[0], ra->file_idx, ra->offset);
        }
        h = h->next;
    }
    printf ("\n");
}

static void list_insert_reader (candidate_reader ** h, candidate_reader * q)
{
    candidate_reader * head;
    if (!h || !q)
    {
        printf ("Error: list_insert_reader ()\n");
        return;
    }

    head = * h;
    if (!head)
    {
        * h = q;
        q->next = NULL;

        return;
    }
#if 0
    while (head->next)
    {
        head = head->next;
    }

    head->next = q;
    q->next = NULL;
#else
    q->next = head;
    * h = q;
#endif
    return;
}

/* This routine append q to list h. We pass t, which points to the tail of h
   to speed up list insertion.
*/
static void list_append_reader_list (candidate_reader ** h, candidate_reader ** t, candidate_reader * q)
{
    candidate_reader * tail;
    if (!h || !t || !q)
    {
        printf ("Error: list_append_reader_list: h: %d, q: %d\n", h == 0, q == 0);
        return;
    }

    tail = * t;
    // tail is null also means head is null.
    if (!tail)
    {
        * h = q;
    }
    else
    {
        tail->next = q;
    }

    while (q->next)
    {
        q = q->next;
    }

    * t = q;

    return;
}

static candidate_reader * list_remove_reader (candidate_reader ** h)
{
    candidate_reader * t = NULL;
    if (h && * h)
    {
        t = * h;
        *h = (* h)->next;
    }

    return t;
}

static int list_get_length (candidate_reader * h)
{
    int l = 0;

    while (h)
    {
        h = h->next;
        l++;
    }

    return l;
}

static int calc_data_size (struct proc_struct * p)
{
    int i, size = 0;
    candidate_reader * h = p->local_read_request_list;

    // message type
//    size += 4;

    // count
//    size += 4;

    while (h)
    {
        // rank + varid + ndim + start + count
        size += 4 + 4 + 4 + h->ra->ndims * (8 + 8) + 8 + 4 + 8;
        h = h->next;
    }

    return size;
}

static void buffer_write (void ** buffer, void * data, int size)
{
    memcpy (* buffer, data, size);
    * buffer = * buffer + size;
}

// *****************************************************************************
static void _buffer_write (char ** buffer, uint64_t * buffer_size
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

// *****************************************************************************
static void _buffer_write32 (char ** buffer, uint32_t * buffer_size
                            ,uint32_t * buffer_offset
                            ,const void * data, uint32_t size
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
                             "Requested: %u\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}

// *****************************************************************************
static void _buffer_read (char * buffer, uint64_t * buffer_offset
                         ,void * data, uint64_t size
                         )
{
    memcpy (data, buffer + *buffer_offset, size);
    *buffer_offset += size;
}

static void sort_read_requests (struct proc_struct * p, candidate_reader ** h)
{
    int file_idx;
    uint64_t offset;
    candidate_reader * r = * h;
    candidate_reader * n = 0, * t, * t_prev, * next;
    while (r)
    {
//printf ("[%d]: r->ra->offset = %llu\n", p->rank, r->ra->offset);
        t = n;
        t_prev = 0;
        next = r->next;

        file_idx = r->ra->file_idx;
        offset = r->ra->offset;

        while (t && t->ra->file_idx <= file_idx)
        {
            if (t->ra->file_idx == file_idx && t->ra->offset > offset)
            {
                break;
            }

            t_prev = t;
            t = t->next;
        }

        if (!t_prev)
        {
            r->next = n;
            n = r;
        }
        else
        {
            t_prev->next = r;
            r->next = t;
        }
//if (p->rank == 1) list_print_readers (p, n);
        r = next;
    }

    * h = n;
}

static int is_short_list (candidate_reader * r)
{
    /* 0 element */
    if (!r)
    {
        return 1;
    }

    /* 1 element */
    if (!r->next)
    {
       return 1;
    }

    /* 2 elements */
    if (!r->next->next)
    {
        return 1;
    }

    return 0;
}

/* Divide the list into sub-lists */
static void divide_list (struct proc_struct * p, candidate_reader * h
                        ,candidate_reader ** l1
                        ,candidate_reader ** l2)
{
    candidate_reader * r1 = 0, * r2 = 0, * r, * hn, * rn;
    int i = 0;

    /* split into l1 and l2 in turns. */
    while (h)
    {
        hn = h->next;

        r = (i % 2 == 0) ? r1 : r2;
  
        h->next = r;
        (i % 2 == 0) ? (r1 = h) : (r2 = h);

        h = hn;

        i++;
    }

    * l1 = r1;
    * l2 = r2;
}

static void combine_list (struct proc_struct * p, candidate_reader **h
                         ,candidate_reader * l1, candidate_reader * l2)
{
    candidate_reader * head = 0, * tail = 0, * l1n, * l2n;

    while (l1 && l2)
    {
        if (l1->ra->file_idx < l2->ra->file_idx
        || (l1->ra->file_idx == l2->ra->file_idx && l1->ra->offset <= l2->ra->offset))
        {
            l1n = l1->next; 
            if (!tail)
            {
                head = tail = l1;
            }
            else
            {
                tail->next = l1;
                tail = tail->next;
            }
 
            tail->next = 0;

            l1 = l1n;
        }
        else
        {
            l2n = l2->next;
            if (!tail)
            {
                head = tail = l2;
            }
            else
            {
                tail->next = l2;
                tail = tail->next;
            }

            tail->next = 0;

            l2 = l2n;
        }
    }

    tail->next = l1 ? l1 : l2; 

    * h = head;
}

static void combine_list_by_rank (struct proc_struct * p, candidate_reader **h
                                 ,candidate_reader * l1, candidate_reader * l2)
{
    candidate_reader * head = 0, * tail = 0, * l1n, * l2n;

    while (l1 && l2)
    {
        if (l1->rank <= l2->rank)
        {
            l1n = l1->next;
            if (!tail)
            {
                head = tail = l1;
            }
            else
            {
                tail->next = l1;
                tail = tail->next;
            }

            tail->next = 0;

            l1 = l1n;
        }
        else
        {
            l2n = l2->next;
            if (!tail)
            {
                head = tail = l2;
            }
            else
            {
                tail->next = l2;
                tail = tail->next;
            }

            tail->next = 0;

            l2 = l2n;
        }
    }

    tail->next = l1 ? l1 : l2;

    * h = head;
}

static void merge_sort_read_requests (struct proc_struct * p, candidate_reader **h)
{
    candidate_reader * r, * r0 = 0, * r1 = 0;

    divide_list (p, * h, &r0, &r1);

    if (is_short_list (r0) || is_short_list (r1))
    {
        sort_read_requests (p, &r0);
        sort_read_requests (p, &r1);
    }
    else
    {
        merge_sort_read_requests (p, &r0);
        merge_sort_read_requests (p, &r1);
    }

    combine_list (p, h, r0, r1);
}

/*************************************************
 This routine sorts split_read_request_list1 by
 the r->rank of each entry. This is to make the
 subsequent data shuffling easy.
**************************************************/
static void sort_list_by_rank (struct proc_struct * p, candidate_reader ** list)
{
    candidate_reader * r = * list;
    candidate_reader * next = 0, * h = 0, * l = 0, * m = 0;

    while (r)
    {
        next = r->next;

        if (!h)
        {
            h = r;
            h->next = 0;
        }
        else
        {
            l = h;
            m = 0;

            while (l && l->rank < r->rank)
            {
                m = l;
                l = l->next;
            }

            if (!m)
            {
                r->next = h;
                h = r;
            }
            else
            {
                m->next = r;
                r->next = l;
            }
        }

        r = next;
    }

    * list = h;
}

static void merge_sort_list_by_rank (struct proc_struct * p, candidate_reader ** h)
{
    candidate_reader * r, * r0 = 0, * r1 = 0;

    divide_list (p, * h, &r0, &r1);

    if (is_short_list (r0) || is_short_list (r0)) 
    {
        sort_list_by_rank (p, &r0);
        sort_list_by_rank (p, &r1);
    }
    else
    {
        merge_sort_list_by_rank (p, &r0);
        merge_sort_list_by_rank (p, &r1);
    }

    combine_list_by_rank (p, h, r0, r1);
}

/***************************************************
 This routine sorts list by varid. If there is a tie,
 we break it by start then count. This is to make the
 subsequent data shuffling easy.
**************************************************/
static void sort_list_by_varid (candidate_reader ** list)
{
    candidate_reader * r = * list;
    candidate_reader * next = 0, * h = 0, * l = 0, * m = 0;

    while (r)
    {
        next = r->next;

        if (!h)
        {
            h = r;
            h->next = 0;
        }
        else
        {
            l = h;
            m = 0;

            while (l && l->ra->varid < r->ra->varid)
            {
                m = l;
                l = l->next;
            }

            if (!m)
            {
                r->next = h;
                h = r;
            }
            else
            {
                m->next = r;
                r->next = l;
            }
        }

        r = next;
    }

    * list = h;
}

/* Helper routine that maps reader rank to a group */
static int rank_to_group_mapping (struct proc_struct * p, int rank)
{
    int grp_size = p->size / p->groups;
    int remain = p->size - grp_size * p->groups;
    int g;

    if (remain == 0)
    {
        g = rank / grp_size;
    }
    else
    {
        if (rank < (grp_size + 1) * remain)
        {
            g = rank / (grp_size + 1);
        }
        else
        {
            g = remain + (rank - (grp_size + 1) * remain) / grp_size;
        }
    }

    return g;
}

/* Helper routine that maps file to a group */
static int fidx_to_group_mapping (struct proc_struct * p, int fidx)
{
    int grp_size = p->n_total_sf / p->groups;
    int remain = p->n_total_sf - grp_size * p->groups;;
    int g;

    if (remain == 0)
    {
        g = fidx / grp_size;
    }
    else
    {
        if (fidx < (grp_size + 1) * remain)
        {
            g = fidx / (grp_size + 1);
        }
        else
        {
            g = remain + (fidx - (grp_size + 1) * remain) / grp_size;
        }
    }

    return g;
}
/* Check whether r's read parameters are the same as (varid, ndims, start, count). */
static int is_same_read_request (candidate_reader * r, int varid, int ndims
                                ,uint64_t * start, uint64_t * count)
{
    int i;

    if (r->ra->varid == varid && r->ra->ndims == ndims)
    {
        for (i = 0; i < ndims; i++)
        {
            if (start[i] != r->ra->start[i] || count[i] != r->ra->count[i])
            {
                return 0;
            }
        }
        
        return 1;
    }

    return 0;
}


/* Check whether r's read parameters are the same as (varid, ndims, start, count). */
static candidate_reader * find_local_read_request (struct proc_struct * p, int varid, int ndims
                                                  , uint64_t * start, uint64_t * count)
{
    int i, j, idx;
    candidate_reader * r;

    idx = start[0] % HT_SIZE;

    for (i = 0; i < p->ht[idx].size; i++)
    {
        r = p->ht[idx].h[i];
        if (r->ra->varid == varid && r->ra->ndims == ndims)
        {
            for (j = 0; j < ndims; j++)
            {
                if (start[j] != r->ra->start[j] || count[j] != r->ra->count[j])
                {
                    break;
                }
             }
  
             if (j == ndims)
             { 
                 return r;
             }
        }
    }

    return 0;
}

/* This routine returns a vector that contains the file index
   that a read request will access.
*/
static int * get_subfile_idx_vector (struct proc_struct * p
                                    ,candidate_reader * r
                                    ,int * len
                                    )
{
    candidate_reader * n1 = 0, * n2 = 0, * t;
    int c, * idx_vector, prev;

    split_read_request (p->gp, r, &n1, &n2);

    sort_read_requests (p, &n1);
    sort_read_requests (p, &n2);

    /* First pass to calculate the size of vector */
    c = 0;
    t = n1;
    prev = -1;
    while (t)
    {
        if (prev != t->ra->file_idx)
        {
            c++;
        }

        prev = t->ra->file_idx;
        t = t->next;
    }

    t = n2;
    prev = -1;
    while (t)
    {
        if (prev != t->ra->file_idx)
        {
            c++;
        }

        prev = t->ra->file_idx;
        t = t->next;
    }

    * len = c;
    idx_vector = (int *) malloc (4 * c);

    t = n1;
    prev = -1;
    c = 0;
    while (t)
    {
        if (prev != t->ra->file_idx)
        {
            idx_vector[c] = t->ra->file_idx;
            c++;
        }

        prev = t->ra->file_idx;
        t = t->next;
    }

    t = n2;
    prev = -1;
    while (t)
    {
        if (prev != t->ra->file_idx)
        {
            idx_vector[c] = t->ra->file_idx;
            c++;
        }

        prev = t->ra->file_idx;
        t = t->next;
    }

    return idx_vector;
}

static void free_subfile_idx_vector (int * idx_vector)
{
    free (idx_vector);
}

/* This routine copies a new request with file_idx
   field set to file_idx. It ignores data field and 
   set it to 0.
*/
candidate_reader * copy_request (candidate_reader * r, int file_idx)
{
    candidate_reader * n =  (candidate_reader *) malloc (sizeof (candidate_reader));
    assert (n);

    n->rank = r->rank;

    n->ra = (read_args *) malloc (sizeof (read_args));
    assert (n->ra);

    n->ra->varid = r->ra->varid;
    n->ra->ndims = r->ra->ndims;

    n->ra->start = (uint64_t *) malloc (8 * n->ra->ndims);
    assert (n->ra->start);
    memcpy (n->ra->start, r->ra->start, 8 * n->ra->ndims);

    n->ra->count = (uint64_t *) malloc (8 * n->ra->ndims);
    assert (n->ra->count);
    memcpy (n->ra->count, r->ra->count, 8 * n->ra->ndims);

    n->ra->data = 0;
    n->ra->size = r->ra->size;
    n->ra->file_idx = file_idx;
    n->ra->offset = r->ra->offset;
    n->ra->parent = 0;
    n->next = 0;
 
    return n;
}

/* This routine shuffles read requests by the file id they access.
*/
static void shuffle_request (struct proc_struct * p)
{
    candidate_reader * h1, * h2 = 0, * r, * n, * next, * prev;
    int len, i, n_grps, * fidx_vector, size;
    char * b = 0, * recv_buffer2 = 0;
    uint32_t * recv_sizes = 0;
    uint32_t offset, total_size = 0, buffer_size = 0;
    int * offsets = 0, * sizes = 0, * offsets2 = 0, * sizes2 = 0;
    int g, g_prev, offset_prev, varid, ndims;

    // h1 is local read request list
    // h2 is constructed based upon h1 to make request shuffle easier.
    // Note: if a request needs to be served by two aggregators, we will
    // make two copies of it in h2.
    h1 = p->local_read_request_list;
    while (h1)
    {
        next = h1->next;

        fidx_vector = get_subfile_idx_vector (p, h1, &len);

/*
if (p->rank == 0)
{
    for (i = 0; i < len; i++)
    {
        printf ("fidx_vector[%d] = %d", i, fidx_vector[i]);
    }
}
printf ("\n");
*/
        for (i = 0; i < len; i++)
        {
            g = fidx_to_group_mapping (p, fidx_vector[i]);
            g_prev = (i == 0 ? -1 : fidx_to_group_mapping (p, fidx_vector[i - 1]));

            if (i != 0 && g == g_prev)
            {
                continue;
            }

            n = copy_request (h1, fidx_vector[i]);

            r = h2;
            prev = 0;

            if (!r)
            {
                h2 = n;
            }
            else
            {
                while (r)
                {
//printf ("loop 2: %d\n", r->ra->file_idx);
                    if (fidx_vector[i] < r->ra->file_idx)
                    {
                        if (!prev)
                        {
                            n->next = h2;
                            h2 = n;
                            prev = 0;
                        }
                        else
                        {
                            prev->next = n;
                            n->next = r;
                        }

                        break;
                    }

                    if (!r->next)
                    {
                        r->next = n;
                        break;
                    }

                    prev = r;
                    r = r->next;
                }
            }
        }

        free_subfile_idx_vector (fidx_vector);

        h1 = next;
    }
/*
if (p->rank == 0)
{
printf ("****************\n");
        list_print_readers (p, h2);
printf ("****************\n");
}
*/
    MPI_Comm_size (p->new_comm2, &size);

    offsets = malloc (size * 4);
    sizes = malloc (size * 4);
    assert (offsets);
    assert (sizes);

    /* Init offsets and sizes */
    for (i = 0; i < size; i++)
    {
        offsets[i] = INVALID_VALUE;
        sizes[i] = 0;
    }

    // write h2 to buffer
    g_prev = -1;
    offset_prev = 0;
    offset = 0;

    r = h2;
    while (r)
    {
        next = r->next;

        g = fidx_to_group_mapping (p, r->ra->file_idx);

        /* if a request will access files that are under control
           of other groups/aggregators. 
         */
        if (g != p->group)
        {
            if (offsets[g] == INVALID_VALUE)
            {
                offsets[g] = offset;

                if (g_prev != -1)
                {
                    sizes[g_prev] = offset - offset_prev;
                }

                g_prev = g;
                offset_prev = offset;
            }

            sizes[g] += 4 + 4 + 4 + r->ra->ndims * 8 + r->ra->ndims * 8 + 4 + 4;

            _buffer_write32 (&b, &buffer_size, &offset, &r->rank, 4);
            _buffer_write32 (&b, &buffer_size, &offset, &r->ra->varid, 4);
            _buffer_write32 (&b, &buffer_size, &offset, &r->ra->ndims, 4);
            _buffer_write32 (&b, &buffer_size, &offset, r->ra->start, r->ra->ndims * 8);
            _buffer_write32 (&b, &buffer_size, &offset, r->ra->count, r->ra->ndims * 8);
            _buffer_write32 (&b, &buffer_size, &offset, &r->ra->size, 4);
            _buffer_write32 (&b, &buffer_size, &offset, &r->ra->file_idx, 4);

        }

        free (r->ra->start);
        free (r->ra->count);
        free (r->ra);
        free (r);

        r = next;
    }

    recv_sizes = malloc (size * size * 4);
    assert (recv_sizes);

    MPI_Allgather (sizes, size, MPI_INT
                  ,recv_sizes, size, MPI_INT
                  ,p->new_comm2
                  );

    offsets2 = malloc (size * 4);
    sizes2 = malloc (size * 4);

    assert (offsets2);
    assert (sizes2);

    len = 0;
    offsets2[0] = 0;

    for (i = 0; i < size; i++)
    {
        len += recv_sizes[i * size + p->group];
        sizes2[i] = recv_sizes[i * size + p->group];
        if (i > 0)
        {
            offsets2[i] = offsets2[i - 1] + sizes2[i - 1];
        }
    }

    recv_buffer2 = malloc (len);
    assert (recv_buffer2);

    MPI_Alltoallv (b, sizes, offsets, MPI_BYTE
                  ,recv_buffer2, sizes2, offsets2, MPI_BYTE
                  ,p->new_comm2
                  );

    /* Free sender buffer and let 'be' point to receiver buffer */
    free (b);
    b = recv_buffer2;
    offset = 0;

    while (offset < len)
    {
        n = (candidate_reader *) malloc (sizeof (candidate_reader));
        assert (n);
   
        n->rank = * (uint32_t *) (b + offset);
        offset += 4;

        n->ra = (read_args *) malloc (sizeof (read_args));
        assert (n->ra);

        n->ra->varid = * (uint32_t *) (b + offset);
        offset += 4;

        n->ra->ndims = * (uint32_t *) (b + offset);
        offset += 4;

        n->ra->start = malloc (n->ra->ndims * 8);
        n->ra->count = malloc (n->ra->ndims * 8);
        assert (n->ra->start);
        assert (n->ra->count);

        memcpy (n->ra->start, b + offset, n->ra->ndims * 8);
        offset += n->ra->ndims * 8;

        memcpy (n->ra->count, b + offset, n->ra->ndims * 8);
        offset += n->ra->ndims * 8;

        n->ra->size = * (uint32_t *) (b + offset);
        offset += 4;
        
        n->ra->file_idx = * (uint32_t *) (b + offset);
        offset += 4;

        n->ra->data = 0;
        n->ra->offset = 0;

        if (n->rank != p->rank)
        {
            list_insert_reader (&p->local_read_request_list, n);
        }
        else
        {
            free (n->ra->start);
            free (n->ra->count);
            free (n->ra);
            free (n);
        }

    }

    free (sizes);
    free (sizes2);
    free (offsets);
    free (offsets2);
    free (recv_sizes);
    free (recv_buffer2);
}
 
/* This routine shuffles data among aggregators.
   It is called after do_read(), which populates
   split_read_request_list1, which contains data
   for both the local group and other remotes groups.
   In the end of the routine, it copies read data from
   split_read_request_list1 and split_read_request_list2
   to local_read_request_list.
*/
static void shuffle_data (struct proc_struct * p)
{
    candidate_reader * r, * parent, * recv_list, * n;
    char * b = 0, * recv_buffer2 = 0;
    uint32_t * recv_sizes = 0;
    uint32_t offset, total_size = 0, len, buffer_size = 0;
    int * offsets = 0, * sizes = 0, * offsets2 = 0, * sizes2 = 0;
    int i, size, g, g_prev, offset_prev, varid, ndims;
    uint64_t * start, * count;

    MPI_Comm_size (p->new_comm2, &size);

    offsets = malloc (size * 4);
    sizes = malloc (size * 4);
    assert (offsets);
    assert (sizes);

    /* Init offsets and sizes */
    for (i = 0; i < size; i++)
    {
        offsets[i] = INVALID_VALUE;
        sizes[i] = 0;
    }
    
    /* Sort the split_read_request_list1 by the r->rank.
       This makes the subsequent buffer construction (for the alltoallv call) easier,
       since the data that fall under an aggregator will be consecutive list entries.
     */
    //sort_list_by_rank (p, &p->split_read_request_list1);
    merge_sort_list_by_rank (p, &p->split_read_request_list1);

    g_prev = -1;
    offset_prev = 0;
    offset = 0;

    /* Reminder: any read request in split_read_request_list1 is the one
       that needs to be read from disk. The list can contain
         1. requests from the local group.
         2. requests from other groups.
       These requests are constructed via shuffle_requests call.
       After the list is sorted by rank (sort_list_by_rank), requests that
       are from the same group will be grouped together. The buffer constructed looks like

       -----------------------------------------------------------------------------------------------------------------
       | all the data aggregator 0 needs | all the data aggregator 1 needs | ......| all the data aggregator n-1 needs |
       -----------------------------------------------------------------------------------------------------------------

       sizes[g]  -  the size of each segment above
       offsets[g]-  the offset of each segment above
    */
    r = p->split_read_request_list1;
    while (r)
    {
        /* r->rank - the rank that waits for the data.
           g - the group that we need to send the data to
         */
        g = rank_to_group_mapping (p, r->rank);

        /* To make sure we handle a request that belongs to other aggregators.
           p->group is the calling aggregator's group id.
         */
        if (g != p->group)
        {
            if (offsets[g] == INVALID_VALUE)
            {
                offsets[g] = offset;

                if (g_prev != -1)
                {
                    sizes[g_prev] = offset - offset_prev;
                }

                g_prev = g;
                offset_prev = offset;
            }

            sizes[g] += 4 + 4 + r->ra->ndims * 8 + r->ra->ndims * 8 + 4 + r->ra->size;

            _buffer_write32 (&b, &buffer_size, &offset, &r->ra->varid, 4);
            _buffer_write32 (&b, &buffer_size, &offset, &r->ra->ndims, 4);
            _buffer_write32 (&b, &buffer_size, &offset, r->ra->start, r->ra->ndims * 8);
            _buffer_write32 (&b, &buffer_size, &offset, r->ra->count, r->ra->ndims * 8);
            _buffer_write32 (&b, &buffer_size, &offset, &r->ra->size, 4);
            _buffer_write32 (&b, &buffer_size, &offset, r->ra->data, r->ra->size);

            /* NOTE: r-ra->data in child is not needed anymore.
               Release it right away so that the overall memory consumption 
               is minimized.
            */
            free (r->ra->data);
            r->ra->data = 0;
        }

        r = r->next;
    }
#if 0
if (p->rank = 2)
    for (i = 0; i < size; i++)
    {
        if (sizes[i] != 0)
        {
            printf ("aggregator: %d offsets[%d] = %d, sizes[%d] = %d  ", p->rank, i, offsets[i], i, sizes[i]);
        }
    }
    printf ("\n");
#endif

    recv_sizes = malloc (size * size * 4);
    assert (recv_sizes);

    /* Since each aggregator maintains the 'sizes' array, the entire thing is of size * size */
    MPI_Allgather (sizes, size, MPI_INT
                  ,recv_sizes, size, MPI_INT
                  ,p->new_comm2
                  );

    offsets2 = malloc (size * 4);
    sizes2 = malloc (size * 4);

    assert (offsets2);
    assert (sizes2);

    len = 0;
    offsets2[0] = 0;

    /* len - the size of buffer that contains all segments that is targeted at this group
       sizes2 - the size of each segment that is targeted at this group
       offsets2 - the offset of each segment that is targeted at this group
     */
      
    for (i = 0; i < size; i++)
    {
        len += recv_sizes[i * size + p->group];
        sizes2[i] = recv_sizes[i * size + p->group];
        if (i > 0)
        {
            offsets2[i] = offsets2[i - 1] + sizes2[i - 1];
        }
    }

    recv_buffer2 = malloc (len);
    assert (recv_buffer2);

#if DEBUG
    for (i = 0; i < size; i++)
    {
        printf ("aggregator : %d, offsets2[%d] = %d, sizes2[%d] = %d\n", p->rank, i, offsets2[i], i, sizes2[i]);
    }
#endif

#if DEBUG
    struct timeval t0, t1, t2;
    gettimeofday (&t0, NULL);
#endif

    MPI_Alltoallv (b, sizes, offsets, MPI_BYTE
                  ,recv_buffer2, sizes2, offsets2, MPI_BYTE
                  ,p->new_comm2
                  );
#if DEBUG
    gettimeofday (&t1, NULL);
    printf ("alltoall = %7.5g\n", t1.tv_sec - t0.tv_sec + (double)(t1.tv_usec - t0.tv_usec)/1000000);
#endif
    /* Free sender buffer and let 'be' point to receiver buffer */
    free (b);
    b = recv_buffer2;
    offset = 0;

#if 1
    build_hash (p);

    while (offset < len)
    {
        varid = * (uint32_t *) (b + offset);
        offset += 4;

        ndims = * (uint32_t *) (b + offset);
        offset += 4;

        start = malloc (ndims * 8);
        count = malloc (ndims * 8);
        assert (start);
        assert (count);

        memcpy (start, b + offset, ndims * 8);
        offset += ndims * 8;

        memcpy (count, b + offset, ndims * 8);
        offset += ndims * 8;
/*
        r = p->split_read_request_list2;
        while (r)
        {
            if (!r->ra->data && is_same_read_request (r, varid, ndims, start, count))
            {
                uint32_t l = * (uint32_t *) (b + offset);
                offset += 4;

                r->ra->data = malloc (l);
                memcpy (r->ra->data, b + offset, l);

                offset += l;
                break;

            }

            r = r->next;
        }

        assert (r);
*/

        r = find_local_read_request (p, varid, ndims, start, count);
        assert (r);

        uint32_t l = * (uint32_t *) (b + offset);
        offset += 4;

        r->ra->data = malloc (l);
        memcpy (r->ra->data, b + offset, l);
        offset += l;

        free (start);
        free (count);
    }

    free_hash (p);
#else
    recv_list = 0;
    while (offset < len)
    {
        n = (candidate_reader *) malloc (sizeof (candidate_reader));
        assert (n);

        n->ra = (read_args *) malloc (sizeof (read_args));
        assert (n->ra);

        n->ra->varid = * (uint32_t *) (b + offset);
        offset += 4;

        n->ra->ndims = * (uint32_t *) (b + offset);
        offset += 4;

        n->ra->start = malloc (n->ra->ndims * 8);
        n->ra->count = malloc (n->ra->ndims * 8);
        assert (n->ra->start);
        assert (n->ra->count);

        memcpy (n->ra->start, b + offset, n->ra->ndims * 8);
        offset += n->ra->ndims * 8;

        memcpy (n->ra->count, b + offset, n->ra->ndims * 8);
        offset += n->ra->ndims * 8;

        n->ra->size = * (uint32_t *) (b + offset);
        offset += 4;

        /* NOTE: here we simply point to the data in the receive buffer
           to avoid double memory. So later in freeing the list  
         */
        n->ra->data = b + offset; 
        offset += n->ra->size;

        n->next = 0;

        list_insert_reader (&recv_list, n);
    }

    /* Sort these lists so that one-to-one matching is established */
    sort_list_by_varid (&recv_list);
    sort_list_by_varid (&p->split_read_request_list2);

    r = p->split_read_request_list2;
    n = recv_list;
    while (n)
    {
        memcpy (r->ra->data, n->ra->data, n->ra->size);
        n = n->next;
    }

    /* Note: this following call doesn't free ra->data */
    free_candidate_reader_list (recv_list);
#endif
#if DEBUG
    gettimeofday (&t2, NULL);
    printf (" after alltoall = %7.5g\n", t2.tv_sec - t1.tv_sec + (double)(t2.tv_usec - t1.tv_usec)/1000000);
#endif
    free (sizes);
    free (sizes2);
    free (offsets);
    free (offsets2);
    free (recv_sizes);
    free (recv_buffer2);

    /* Traverse split_read_request_list1 and copy data from child
       to parent. Check whether r->ra->data is null or not. If it 
       is null, it is a request from other aggregators.
    */
/*
if (p->rank == 0)
{
candidate_reader * tr = p->local_read_request_list;

while (tr)
{
    if (tr->ra->data)
        printf ("[%d]: data = %4.6f %4.6f\n", tr->rank, * (double *) tr->ra->data, * ((double *) tr->ra->data + 1));
    else
        printf ("reader : %d\n", tr->rank);
    tr = tr->next;
}

//list_print_readers (p, p->local_read_request_list);
printf ("*********************\n");
}
*/
//while (1);
    r = p->split_read_request_list1;
    while (r)
    {
        if (r->ra->data)
        {
            copy_buffer (p->gp, r->ra->parent, r);

            free (r->ra->data);
            r->ra->data = 0;
        }

        r = r->next;
    }

    /* Traverse split_read_request_list2 and copy data from child
       to parent.
    */
    r = p->split_read_request_list2;
/*
if (p->rank == 0)
{
candidate_reader * tr = p->split_read_request_list2;
while (tr)
{
    printf ("data = %4.6f %4.6f\n", * (double *) tr->ra->data, * ((double *) tr->ra->data + 1));
    tr = tr->next;
}
}
*/

    while (r)
    {
        copy_buffer (p->gp, r->ra->parent, r);

        free (r->ra->data);
        r->ra->data = 0;
     
        r = r->next;
    }
}

/* Send read data from aggregator to member processors */
static void send_read_data1 (struct proc_struct * p)
{
    int g;
    void * b = 0;

    candidate_reader * r = p->local_read_request_list;

    while (r)
    {
        g = rank_to_group_mapping (p, r->rank);

        /*  g == p->group       -> requests are from processors that belong to the group
            p->rank != r->rank  -> requests are NOT from aggregator itself
        */
        if (g == p->group && p->rank != r->rank)
        {
            assert (r->ra->data);

            MPI_Send (r->ra->data, r->ra->size, MPI_BYTE
                     ,r->rank - p->aggregator_rank, TAG_DATA, p->new_comm
                     );  
        }

        r = r->next;
    }
}

/* Send read data from aggregator to member processors */
static void send_read_data (struct proc_struct * p)
{
    int g, i;
    char * b = 0, * recv_buff = 0;
    uint32_t offset = 0, buffer_size = 0;
    int size, * sizes = 0, * offsets = 0;
    candidate_reader * r = p->local_read_request_list;

    MPI_Comm_size (p->new_comm, &size);

    sizes = (int *) malloc (size * 4);
    offsets = (int *) malloc (size * 4);
    assert (sizes);
    assert (offsets);

    for (i = 0; i < size; i++)
    {
        sizes[i] = 0;
        offsets[i] = INVALID_VALUE;
    }

    while (r)
    {
        g = rank_to_group_mapping (p, r->rank);

        /*  g == p->group       -> requests are from processors that belong to the group
            p->rank != r->rank  -> requests are NOT from aggregator itself
        */

        if (g == p->group && p->rank != r->rank)
        {
            assert (r->ra->data);

            if (offsets[r->rank - p->aggregator_rank] == INVALID_VALUE)
            {
                offsets[r->rank - p->aggregator_rank] = offset;
            }

            sizes[r->rank - p->aggregator_rank] += r->ra->size;

            _buffer_write32 (&b, &buffer_size, &offset, r->ra->data, r->ra->size);

            /* Free it immediately to avoid double buffering */
            free (r->ra->data);
            r->ra->data = 0;
        }

        r = r->next;
    }

    MPI_Scatterv (b, sizes, offsets, MPI_BYTE
                 ,recv_buff, 0, MPI_BYTE
                 ,0, p->new_comm
                 );

    free (b);
    free (sizes);
    free (offsets);
}

/* Receive read data from aggregator */
static void get_read_data1 (struct proc_struct * p)
{
    MPI_Status status;
    candidate_reader * r = p->local_read_request_list;

    if (p->rank != p->aggregator_rank)
    {
        while (r) 
        {
            MPI_Recv (r->ra->data, r->ra->size, MPI_BYTE
                     ,MPI_ANY_SOURCE, MPI_ANY_TAG, p->new_comm
                     ,&status
                     );

            r = r->next;
        }
    }
}

/* Receive read data from aggregator */
static void get_read_data (struct proc_struct * p)
{
    void * b = 0, * recv_buff = 0;
    int * sizes = 0, * offsets = 0;
    int size = 0;
    MPI_Status status;
    candidate_reader * r = p->local_read_request_list;

    if (p->rank == p->aggregator_rank)
    {
        return;
    }

    r = p->local_read_request_list;
    while (r)
    {
        size += r->ra->size;
        r = r->next;
    }

    recv_buff = malloc (size);
    if (recv_buff == 0)
    {
        printf ("Warning: the size of data is 0\n");
        return;
    }

    MPI_Scatterv (b, sizes, offsets, MPI_BYTE
                 ,recv_buff, size, MPI_BYTE
                 ,0, p->new_comm
                 );

    b = recv_buff;
    r = p->local_read_request_list;
    while (r)
    {
        memcpy (r->ra->data, b, r->ra->size);
        b += r->ra->size;

        r = r->next;
    }

    free (recv_buff);
}

static candidate_reader * parse_buffer (struct proc_struct * p, void * b, int len)
{
    candidate_reader * h = 0;
    int i, j, type, count, varid, ndims;
    uint32_t l = 0;
    void * b1 = b; 
    void * buf;
    candidate_reader * r;

    while (l < len)
    {
        r = (candidate_reader *) malloc (sizeof (candidate_reader));
        assert (r);

        r->rank = * (uint32_t *) b;
        b += 4;
        l += 4;

        r->ra = (read_args *) malloc (sizeof (read_args));
        assert (r->ra);
   
        r->ra->varid = * (uint32_t *) b;
        b += 4;
        l += 4;

        r->ra->ndims = * (uint32_t *) b;
        b += 4;
        l += 4;

        r->ra->start = (uint64_t *) malloc (r->ra->ndims * 8);
        r->ra->count = (uint64_t *) malloc (r->ra->ndims * 8);
        assert (r->ra->start);
        assert (r->ra->count);

        memcpy (r->ra->start, b, r->ra->ndims * 8);
        b += r->ra->ndims * 8;
        l += r->ra->ndims * 8;

        memcpy (r->ra->count, b, r->ra->ndims * 8);
        b += r->ra->ndims * 8;
        l += r->ra->ndims * 8;

        r->ra->size = * (uint64_t *) b;
        b += 8;
        l += 8;

        r->ra->file_idx = * (uint32_t *) b;
        b += 4;
        l += 4;

        r->ra->offset = * (uint64_t *) b;
        b += 8;
        l += 8;

        /* DO NOT ALLOCATE PARENT BUFFER AT THIS POINT. Delay it to the point when necessary.
           Do it when copying child buffer to parent (and after it, child buffer
           will be released) so that the overall memory consumption won't be tripled.
        */
        r->ra->data = 0;

        r->ra->parent = 0;

        r->next = 0;

        if (r->rank != p->rank)
        {
            list_insert_reader (&h, r);
        }
        else
        {
            free (r->ra->start);
            free (r->ra->count);
            free (r->ra);
            free (r);
        }
    }

    return h;
}

/* This routine split a read requests into smaller
   requests which fits PG boundary. For example, a BP
   file contains two PG's with one PG contains an array [0,99]
   and the other contains [100,199]. A read request of [0,199] will
   be split into two requests, i.e., [0,99] and [100,199].
   h1 points to those sub-requests that access subfiles [f1,f2].
   h2 points to those that don't access [f1,f2], which means data
   will be read by other aggregators and sent over.
*/
static void split_read_request (ADIOS_GROUP * gp
                               ,candidate_reader * r
                               ,candidate_reader ** h1
                               ,candidate_reader ** h2
                               )
{
    struct BP_GROUP * gh;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * v;
    int i, j, k, idx, t, varid;
    int start_time, stop_time, start_idx, stop_idx, f_idx;
    int ndim, ndim_notime, has_subfile, file_is_fortran;
    uint64_t size, * dims = 0;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t count_notime[32], start_notime[32], * start, * count;
    int file_tdim = -1, read_arg_tdim, is_global = 0, flag;
    struct proc_struct * p;
    int f1, f2;

    * h1 = * h2 = 0;
    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;
    p = fh->priv;
    f1 = p->f1;
    f2 = p->f2;
    varid = r->ra->varid;
    start = r->ra->start;
    count = r->ra->count;

    file_is_fortran = is_fortran_file (gp);
    has_subfile = has_subfiles (gp);

    v = adios_find_var_byid (gp, varid);

    /* Get dimensions and flip if caller != writer language */
    adios_read_bp_staged1_get_dimensions (v
                                        ,fh->tidx_stop - fh->tidx_start + 1
                                        ,file_is_fortran
                                        ,&ndim
                                        ,&dims
                                        ,&file_tdim
                                        );

    if (file_is_fortran)
    {
        swap_order (ndim, dims, &file_tdim);
    }

    if (isTimeless (file_tdim))
    {
        start_time = fh->tidx_start;
        stop_time = fh->tidx_stop;
        ndim_notime = ndim;

        memcpy (start_notime, start, ndim * 8);
        memcpy (count_notime, count, ndim * 8);
    }
    else
    {
        read_arg_tdim = futils_is_called_from_fortran () ? ndim - 1 : 0;

        start_time = start[read_arg_tdim] + fh->tidx_start;
        stop_time = start_time + count[read_arg_tdim] - 1;
        ndim_notime = ndim - 1;

        memcpy (start_notime
               ,futils_is_called_from_fortran () ? start : start + 1
               ,ndim_notime * 8
               );
        memcpy (count_notime
               ,futils_is_called_from_fortran () ? count : count + 1
               ,ndim_notime * 8
               );
    }

    if (futils_is_called_from_fortran ())
    {
        _swap_order (ndim_notime, count_notime);
        _swap_order (ndim_notime, start_notime);
    }

    for (t = start_time; t <= stop_time; t++)
    {
        start_idx = get_var_start_index (v, t);
        stop_idx = get_var_stop_index (v, t);

        if (start_idx < 0 || stop_idx < 0)
        {
            fprintf (stderr,"Variable (id=%d) has no data at %d time step\n",
                varid, t);
            continue;
        }

        if (ndim_notime == 0)
        {
            /* THIS IS A SCALAR VARIABLE */
            idx = 0;

            if ( (v->characteristics[start_idx + idx].file_index >= f1
                 && v->characteristics[start_idx + idx].file_index <= f2)
               || rank_to_group_mapping (p, r->rank) == p->group
               )
            {
                /* MALLOC CHILD REQUEST */
                candidate_reader * n = (candidate_reader *) malloc (sizeof (candidate_reader));
                assert (n);

                n->rank = r->rank;
                n->ra = (read_args *) malloc (sizeof (read_args));
                assert (n->ra);

                n->ra->varid = r->ra->varid;
                n->ra->ndims = r->ra->ndims;
                n->ra->file_idx = v->characteristics[start_idx + idx].file_index;
                n->ra->offset = v->characteristics[start_idx + idx].payload_offset;
                n->ra->parent = r;
                n->ra->start = 0;
                n->ra->count = 0;
                n->ra->data = 0;
                n->next = 0;

                if (v->characteristics[start_idx + idx].file_index >= f1
                 && v->characteristics[start_idx + idx].file_index <= f2)
                {
                    list_insert_reader (h1, n);
                }
                else
                {
                    list_insert_reader (h2, n);
                }
            }
    
            if (isTimeless (file_tdim) )
                break;
            else
                continue;
        }

         /* READ AN ARRAY VARIABLE */
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));

        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
        {
            if (((  v->characteristics[start_idx + idx].file_index < f1
                   || v->characteristics[start_idx + idx].file_index > f2
                   )
                  && rank_to_group_mapping (p, r->rank) != p->group
                  )
               )
            {
                continue;
            }

            idx_table[idx] = 1;
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = adios_read_bp_staged1_get_dimensioncharacteristics(&(v->characteristics[start_idx + idx])
                                                                          ,ldims
                                                                          ,gdims
                                                                          ,offsets
                                                                          );
            if (!is_global)
            {
                memcpy (gdims, ldims, ndim * 8);
            }

            if (file_is_fortran)
            {
                _swap_order (ndim, gdims);
                _swap_order (ndim, ldims);
                _swap_order (ndim, offsets);
            }

            if (!isTimeless (file_tdim))
            {
                for (i = file_tdim; i < ndim - 1; i++)
                {
                    ldims[i] = ldims[i + 1];
                    if (file_is_fortran)
                    {
                        gdims[i] = gdims[i + 1];
                        offsets[i] = offsets[i + 1];
                    }
                }
            }

            /*
            printf("ldims   = "); for (j = 0; j<ndim; j++) printf("%d ",ldims[j]); printf("\n");
            printf("gdims   = "); for (j = 0; j<ndim; j++) printf("%d ",gdims[j]); printf("\n");
            printf("offsets = "); for (j = 0; j<ndim; j++) printf("%d ",offsets[j]); printf("\n");
            printf("count_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
            printf("start_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
            */
            for (j = 0; j < ndim_notime; j++)
            {
                if ( (count_notime[j] > gdims[j])
                  || (start_notime[j] > gdims[j])
                  || (start_notime[j] + count_notime[j] > gdims[j]))
                {
                    fprintf (stderr, "Error: Variable (id=%d, %s) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu]), count_notime[1] = %llu, count_notime[2] = %llu\n",
                        varid, v->var_name, j + 1, count_notime[j], start_notime[j], gdims[j] - 1, count_notime[1], count_notime[2]);
                    return;
                }

                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start_notime[j]
                        && offsets[j] < start_notime[j] + count_notime[j])
                    || (offsets[j] < start_notime[j]
                        && offsets[j] + ldims[j] > start_notime[j] + count_notime[j])
                    || (offsets[j] + ldims[j] > start_notime[j]
                        && offsets[j] + ldims[j] <= start_notime[j] + count_notime[j]);

                idx_table[idx] = idx_table[idx] && flag;
            }

            if (!idx_table[idx])
            {
                continue;
            }

            /* determined how many (fastest changing) dimensions can we read in in one read */
            int hole_break;
            for (i = ndim_notime - 1; i > -1; i--)
            {
                if (offsets[i] == start_notime[i] && ldims[i] == count_notime[i])
                {
                }
                else
                {
                    break;
                }
            }

            hole_break = i;

            candidate_reader * n = (candidate_reader *) malloc (sizeof (candidate_reader));
            assert (n);

            n->rank = r->rank;
            n->ra = (read_args *) malloc (sizeof (read_args));
            assert (n->ra);

            n->ra->varid = r->ra->varid;
            n->ra->ndims = r->ra->ndims;

            n->ra->start = (uint64_t *) malloc (ndim_notime * 8);
            assert (n->ra->start);

            n->ra->count = (uint64_t *) malloc (ndim_notime * 8);
            assert (n->ra->count);

            n->ra->file_idx = v->characteristics[start_idx + idx].file_index;
            n->ra->offset = v->characteristics[start_idx + idx].payload_offset;
            n->ra->parent = r;
            n->ra->data = 0;
            n->next = 0;

            memcpy (n->ra->start, start_notime, ndim_notime * 8);
            memcpy (n->ra->count, count_notime, ndim_notime * 8);

            if (hole_break == -1)
            {
            }
            else if (hole_break == 0)
            {
                int isize;
                uint64_t size_in_dset = 0;
                uint64_t offset_in_dset = 0;

                isize = offsets[0] + ldims[0];
                if (start_notime[0] >= offsets[0])
                {
                    // head is in
                    if (start_notime[0] < isize)
                    {
                        if (start_notime[0] + count_notime[0] > isize)
                        {
                            n->ra->count[0] = isize - start_notime[0];
                        //    size_in_dset = isize - start_notime[0];
                        }
                        else
                        {
                            n->ra->count[0] = count_notime[0];
                        //    size_in_dset = count_notime[0];
                        }
                        n->ra->start[0] = start_notime[0];
                        //offset_in_dset = start_notime[0] - offsets[0];
                    }
                }
                else
                {
                    // middle is in
                    if (isize < start_notime[0] + count_notime[0])
                    {
                        n->ra->count[0] = ldims[0];
                    //    size_in_dset = ldims[0];
                    }
                    else
                    {
                    // tail is in
                        n->ra->count[0] = count_notime[0] + start_notime[0] - offsets[0];
                    //    size_in_dset = count_notime[0] + start_notime[0] - offsets[0];
                    }
                    n->ra->start[0] = offsets[0];
                    //offset_in_dset = 0;
                }

            }
            else
            {
                uint64_t stride_offset = 0;
                int isize;
                uint64_t size_in_dset[MAX_DIMS];
                uint64_t offset_in_dset[MAX_DIMS];
                uint64_t offset_in_var[MAX_DIMS];

                memset(size_in_dset, 0 , MAX_DIMS * 8);
                memset(offset_in_dset, 0 , MAX_DIMS * 8);
                memset(offset_in_var, 0 , MAX_DIMS * 8);

                for (i = 0; i < ndim_notime; i++)
                {
                    isize = offsets[i] + ldims[i];
                    if (start_notime[i] >= offsets[i])
                    {
                        // head is in
                        if (start_notime[i] < isize)
                        {
                            if (start_notime[i] + count_notime[i] > isize)
                            {
                                n->ra->count[i] = isize - start_notime[i];
                            //    size_in_dset[i] = isize - start_notime[i];
                            }
                            else
                            {
                                n->ra->count[i] = count_notime[i];
                            //    size_in_dset[i] = count_notime[i];
                            }
                            //n->ra->start[i] = start_notime[i] - offsets[i];
                            n->ra->start[i] = start_notime[i];
                            //offset_in_dset[i] = start_notime[i] - offsets[i];
                            offset_in_var[i] = 0;
                        }
                    }
                    else
                    {
                        // middle is in
                        if (isize < start_notime[i] + count_notime[i])
                        {
                            n->ra->count[i] = ldims[i];
                        //    size_in_dset[i] = ldims[i];
                        }
                        else
                        {
                            // tail is in
                            n->ra->count[i] = count_notime[i] + start_notime[i] - offsets[i];
                        //    size_in_dset[i] = count_notime[i] + start_notime[i] - offsets[i];
                        }
                        n->ra->start[i] = offsets[i];
                        //offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start_notime[i];
                    }
                }

            }
           
            n->ra->size = bp_get_type_size (v->type, v->characteristics[start_idx + idx].value);
            for (i = 0; i < ndim_notime; i++)
            {
                n->ra->size *= n->ra->count[i];
            }
        
            if (v->characteristics[start_idx + idx].file_index >= f1
               && v->characteristics[start_idx + idx].file_index <= f2)
            {
                list_insert_reader (h1, n);
            }
            else
            {
                list_insert_reader (h2, n);
            }
        }

        free (idx_table);

        if (isTimeless (file_tdim))
            break;
    } // end for (timestep ... loop over timesteps

    if (dims)
    {
        free (dims);
    }
}

#if 1
/* Get number of vars in the file */
int get_var_num (ADIOS_GROUP * gp)
{
    struct BP_GROUP * gh = (struct BP_GROUP *) gp->gh;
    struct adios_index_var_struct_v1 * var_root = gh->vars_root;
    int n = 0;

    while (var_root)
    {
        var_root = var_root->next;
        n++;
    }

    return n;
}

/* Get the number of buckets needed for a hash table */
int get_num_bucket (ADIOS_GROUP * gp, int varid)
{
    struct adios_index_var_struct_v1 * v = adios_find_var_byid (gp, varid);

    return v->characteristics_count;
}

/* To build hashtable for each individual var */
void build_var_hash (ADIOS_GROUP * gp, int varid, int ** htable)
{
    int * ht, nbucket;
 
    if (!htable)
    {
        return;
    }

    nbucket = get_num_bucket (gp, varid);  
    * htable = (int *) malloc (4 * nbucket);
    assert (* htable);
}

/* To build a metadata hash table.
   The hashtable maps 'key' to 'index'.
*/
static void build_hash (struct proc_struct * p)
{
    candidate_reader * r;
    int idx, i, offset;

    p->ht = (bucket *) malloc (HT_SIZE * sizeof (bucket));
    assert (p->ht);

    for (i = 0; i < HT_SIZE; i++)
    {
        p->ht[i].h = 0;
        p->ht[i].size = 0;
        p->ht[i].offset = 0;
    }

    /* First pass to figure out the size of each bucket*/
    r = p->split_read_request_list2;
    while (r)
    {
        idx = r->ra->start[0] % HT_SIZE;
        p->ht[idx].size++; 

        r = r->next;
    }

    /* Second pass to construct HT */
    r = p->split_read_request_list2;
    while (r)
    {
        idx = r->ra->start[0] % HT_SIZE;

        if (!p->ht[idx].h)
        {
            p->ht[idx].h = calloc (p->ht[idx].size, sizeof (candidate_reader *));
            assert (p->ht[idx].h);
        }

        p->ht[idx].h[p->ht[idx].offset++] = r;

        r = r->next;
    }
    

}

/* To free hashtable */
void free_hash (struct proc_struct * p)
{
    int i;

    for (i = 0; i < HT_SIZE; i++)
    {
        if (p->ht[i].h)
        {
            free (p->ht[i].h);
            p->ht[i].h = 0;
        }
    }

    if (p->ht)
    {
        free (p->ht);
        p->ht = 0;
    }
} 
#endif

/* Split orginal read requests into split_read_request_list1
   and split_read_request_list2.
*/
static void process_read_requests (struct proc_struct * p)
{
    candidate_reader * h = p->local_read_request_list;
    candidate_reader * n1, * n2, * t1 = 0, * t2 = 0;
int c = 0;
struct timeval t3, t4, t5;

    while (h)
    {
        n1 = n2 = 0;

gettimeofday (&t3, NULL);
        split_read_request (p->gp, h, &n1, &n2);

gettimeofday (&t4, NULL);
        if (n1)
        {
            list_append_reader_list (&p->split_read_request_list1, &t1, n1);
        }

        if (n2)
        {
            list_append_reader_list (&p->split_read_request_list2, &t2, n2);
        }
gettimeofday (&t5, NULL);

//printf ("[%3d] insdie process_read_request = %f,%f\n", p->rank, t4.tv_sec - t3.tv_sec + (double)(t4.tv_usec - t3.tv_usec)/1000000, t5.tv_sec - t4.tv_sec + (double)(t5.tv_usec - t4.tv_usec)/1000000);
        h = h->next;
c++;
    }

//    printf ("[%3d] count = %d\n", p->rank, c);
}

int get_num_subfiles (BP_FILE * fh)
{
    struct adios_bp_buffer_struct_v1 * b = fh->b;
    struct adios_index_var_struct_v1 ** vars_root = &(fh->vars_root);
    struct bp_minifooter * mh = &(fh->mfooter);
    struct adios_index_var_struct_v1 ** root;
    int i, j, n = 0;

    root = vars_root;
    for (i = 0; i < mh->vars_count; i++)
    {
        for (j = 0; j < (*root)->characteristics_count; j++)
        {
            if ((*root)->characteristics [j].file_index > n)
            {
                n = (*root)->characteristics [j].file_index;
            }
        }
    }

    return n + 1;
}


int is_fortran_file (ADIOS_GROUP * gp)
{
    struct BP_GROUP * gh;
    BP_FILE * fh;

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;

    return (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
}

int has_subfiles (ADIOS_GROUP * gp)
{
    struct BP_GROUP * gh;
    BP_FILE * fh;

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;

    return (fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE);
}

static enum ADIOS_DATATYPES get_type (ADIOS_GROUP * gp, int varid)
{
    struct adios_index_var_struct_v1 * v = adios_find_var_byid (gp, varid);

    return v->type;
}

static int get_type_size (ADIOS_GROUP * gp, int varid)
{
    struct adios_index_var_struct_v1 * v;

    v = adios_find_var_byid (gp, varid);

    return bp_get_type_size (v->type, v->characteristics [0].value);
}

/* This routine copies data from child request buffer to parent request buffer */
static void copy_buffer (ADIOS_GROUP * gp, candidate_reader * parent, candidate_reader * child)
{
    uint64_t * p_start, * p_count, * c_start, * c_count; 
    int i, j, k, varid, size_unit, break_dim;
    uint64_t datasize, dset_stride, var_stride, total_size = 0, items_read, size;
    uint64_t * count_notime, * start_notime;
    read_info ri;

    varid = parent->ra->varid;

    // PARENT read request
    p_start = parent->ra->start;
    p_count = parent->ra->count;

    // CHILD request after split
    c_start = child->ra->start;
    c_count = child->ra->count;

    memset (&ri, 0, sizeof (read_info));

    /* parent requests are the original requests and hence
       are needed to remove time dimension.
    */
    getReadInfo (gp, varid, p_start, p_count, &ri);
    start_notime = ri.start_notime;
    count_notime = ri.count_notime;

    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ri.ndim_notime; i++)
    {
        items_read *= count_notime[i];
    }
    
    size_unit = get_type_size (gp, varid);

    if (ri.ndim_notime == 0)
    {
        if (get_type (gp, varid) == adios_string)
        {
            // strings are stored without \0 in file
            // get_type_size() return the length that includes \0
            size_unit--;
        }

        /* ALLOCATE MEMORY IN PARENT REQUEST */
        if (parent->ra->data == 0)
        {
            parent->ra->data = malloc (parent->ra->size);
            assert (parent->ra->data);
        }

        memcpy (parent->ra->data, child->ra->data, size_unit);

        if (get_type (gp, varid) == adios_string)
        {
            ((char *)parent->ra->data + total_size)[size_unit] = '\0';
        }

        total_size += size_unit;

        if (isTimeless (ri.file_tdim))
        {
//          break;
        }
        else
        {
//          continue;
        }
    }

    /* READ AN ARRAY VARIABLE */
    datasize = 1;
    var_stride = 1;
    dset_stride = 1;
    break_dim =  ri.ndim_notime - 1;

    while (break_dim > -1)
    {
        if (start_notime[break_dim] == c_start[break_dim] && c_count[break_dim] == count_notime[break_dim])
        {
            datasize *= c_count[break_dim];
        }
        else
            break;

        break_dim--;
    }

    if (break_dim == -1)
    {
        if (parent->ra->data == 0)
        {
            parent->ra->data = malloc (parent->ra->size);
            assert (parent->ra->data);
        }
   
        memcpy (parent->ra->data, child->ra->data, child->ra->size);
    }
    else if (break_dim == 0) 
    {
        uint64_t write_offset = (c_start[0] - start_notime[0]) * size_unit;

        for (i = 1; i < ri.ndim_notime; i++)
        {
            write_offset *= count_notime[i];
        }

        /* ALLOCATE PARENT BUFFER IF HAVEN'T DONE SO */
        if (parent->ra->data == 0)
        {
            parent->ra->data = malloc (parent->ra->size);
        }

        assert (parent->ra->data);
        assert (child->ra->data);

        memcpy (parent->ra->data + write_offset, child->ra->data, child->ra->size);
    }
    else 
    {
        uint64_t stride_offset = 0, var_offset = 0, dset_offset = 0;
        uint64_t * size_in_dset, * offset_in_dset, * offset_in_var;

        size_in_dset = malloc (ri.ndim_notime * 8);
        offset_in_dset = malloc (ri.ndim_notime * 8);
        offset_in_var = malloc (ri.ndim_notime * 8);

        assert (size_in_dset);
        assert (offset_in_dset);
        assert (offset_in_var);

        datasize = 1;
        var_stride = 1;
        var_offset = 0;

        for (i = 0; i < ri.ndim_notime; i++)
        {
            offset_in_dset[i] = 0;
            size_in_dset[i] = c_count[i];
            offset_in_var[i] = c_start[i] - start_notime[i];
            var_offset = offset_in_var[i] + var_offset * count_notime[i];
        }
  
        for (i = ri.ndim_notime - 1; i >= break_dim; i--)
        {
            datasize *= size_in_dset[i];
            dset_stride *= c_count[i];
            var_stride *= count_notime[i];
        }
    
        if (parent->ra->data == 0)
        {
            parent->ra->data = malloc (parent->ra->size);
            assert (parent->ra->data);
        }

        copy_data (parent->ra->data
                  ,child->ra->data
                  ,0
                  ,break_dim
                  ,size_in_dset
                  ,c_count
                  ,count_notime
                  ,var_stride
                  ,dset_stride
                  ,var_offset
                  ,dset_offset
                  ,datasize
                  ,size_unit 
                  );

        free (size_in_dset);
        free (offset_in_dset);
        free (offset_in_var);
    }

/*
    data = (char *)data + (items_read * size_unit);
*/
    if (isTimeless (ri.file_tdim))
//      break;

    freeReadInfo (&ri);
}
void adios_read_bp_staged1_read_buffer (ADIOS_GROUP * gp
                                      ,uint64_t buffer_offset
                                      ,candidate_reader * r
                                      ,candidate_reader * s
                                      )
{
    struct BP_GROUP * gh;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * v;
    uint64_t * r_start, * r_count, * s_start, * s_count; 
    int i, j, k, idx, t;
    int varid, start_idx, stop_idx, has_subfile, file_is_fortran;
    uint64_t ldims[MAX_DIMS], gdims[MAX_DIMS], offsets[MAX_DIMS];
    uint64_t datasize, dset_stride, var_stride, total_size = 0, items_read, size;
    uint64_t * count_notime, * start_notime;
    int is_global = 0, size_unit, break_dim, idx_check1, idx_check2;
    uint64_t slice_offset, slice_size;
    void * data;
    read_info ri;

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;
    varid = r->ra->varid;

    // orginal read request
    r_start = r->ra->start;
    r_count = r->ra->count;

    // new read request after split
    s_start = s->ra->start;
    s_count = s->ra->count;

    file_is_fortran = is_fortran_file (gp);
    has_subfile = has_subfiles (gp);

    v = adios_find_var_byid (gp, varid);

    memset (&ri, 0, sizeof (read_info));

    getReadInfo (gp, varid, r_start, r_count, &ri);
//    start_notime = ri.start_notime;
//    count_notime = ri.count_notime;
    start_notime = s_start;
    count_notime = s_count;

    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ri.ndim_notime; i++)
    {
        items_read *= count_notime[i];
    }
    
    size_unit = bp_get_type_size (v->type, v->characteristics [0].value);

    /* For each timestep, do reading separately (they are stored in different sets of process groups */
    for (t = ri.start_time; t <= ri.stop_time; t++)
    {
        start_idx = get_var_start_index (v, t);
        stop_idx = get_var_stop_index (v, t);

        if (start_idx < 0 || stop_idx < 0)
        {
            adios_error (err_no_data_at_timestep
                        ,"Variable (id=%d) has no data at %d time step"
                        ,varid, t
                        );
            continue;
        }

        if (ri.ndim_notime == 0)
        {
            /* READ A SCALAR VARIABLE */
            slice_size = 1 * size_unit;
            idx = 0;

            if (v->type == adios_string)
            {
                // strings are stored without \0 in file
                // size_of_type here includes \0 so decrease by one
                size_unit--;
            }

            slice_offset = v->characteristics[start_idx + idx].payload_offset;

            /* ALLOCATE MEMORY IN CHILD REQUEST */
            s->ra->data = malloc (slice_size);
            assert (s->ra->data);

            data = s->ra->data;

            memcpy (data, fh->b->buff + slice_offset - buffer_offset, size_unit);
            if (fh->mfooter.change_endianness == adios_flag_yes)
            {
                change_endianness (data, size_unit, v->type);
            }

            if (v->type == adios_string)
            {
                // add \0 to the end of string
                // size_of_type here is the length of string
                // FIXME: how would this work for strings written over time?
                ((char *)data + total_size)[size_unit] = '\0';
            }

            total_size += size_unit;

            if (isTimeless (ri.file_tdim))
            {
                break;
            }
            else
            {
                continue;
            }
        }

         /* READ AN ARRAY VARIABLE */
        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
        {
            int flag1, flag2;
            uint64_t payload_size = size_unit;

            datasize = 1;
            var_stride = 1;
            dset_stride = 1;
            idx_check1 = 1;
            idx_check2 = 1;
    
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = adios_read_bp_staged1_get_dimensioncharacteristics (&(v->characteristics[start_idx + idx])
                                                                           ,ldims
                                                                           ,gdims
                                                                           ,offsets
                                                                           );

            if (!is_global)
            {
                memcpy (gdims, ldims, ri.ndim * 8);
            }

            if (file_is_fortran)
            {
                _swap_order (ri.ndim, gdims);
                _swap_order (ri.ndim, ldims);
                _swap_order (ri.ndim, offsets);
            }

            if (!isTimeless (ri.file_tdim))
            {
                for (i = ri.file_tdim; i < ri.ndim - 1; i++)
                {
                    ldims[i] = ldims[i + 1];
                    if (file_is_fortran)
                    {
                        gdims[i] = gdims[i + 1];
                        offsets[i] = offsets[i + 1];
                    }
                }
            }

            for (j = 0; j < ri.ndim_notime; j++)
            {
                payload_size *= ldims [j];
    
                if ( (count_notime[j] > gdims[j]) 
                  || (start_notime[j] > gdims[j]) 
                  || (start_notime[j] + count_notime[j] > gdims[j]))
                {
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], gdims[j] - 1);
                    return;
                }
    
                /* check if there is any data in this pg and this dimension to read in */
                flag1 = (offsets[j] >= start_notime[j] 
                        && offsets[j] < start_notime[j] + count_notime[j])
                    || (offsets[j] < start_notime[j]
                        && offsets[j] + ldims[j] > start_notime[j] + count_notime[j]) 
                    || (offsets[j] + ldims[j] > start_notime[j] 
                        && offsets[j] + ldims[j] <= start_notime[j] + count_notime[j]);

                idx_check1 = idx_check1 && flag1;

                flag2 = (offsets[j] >= s_start[j]
                        && offsets[j] < s_start[j] + s_count[j])
                    || (offsets[j] < s_start[j]
                        && offsets[j] + ldims[j] > s_start[j] + s_count[j])
                    || (offsets[j] + ldims[j] > s_start[j]
                        && offsets[j] + ldims[j] <= s_start[j] + s_count[j]);

                idx_check2 = idx_check2 && flag2;
            }

            if (!idx_check1)
            {
                continue;
            }

            break_dim =  ri.ndim_notime - 1;
            while (break_dim > -1)
            {
                if (start_notime[break_dim] == 0 && ldims[break_dim] == count_notime[break_dim])
                {
                    datasize *= ldims[break_dim];
                }
                else
                    break;

                break_dim--;
            }

            slice_offset = 0;
            slice_size = 0;

            if (break_dim == -1)
            {
                slice_size = payload_size;
    
                slice_offset = v->characteristics[start_idx + idx].payload_offset;

                if (idx_check2)
                {
                    s->ra->data = malloc (slice_size);
                    assert (s->ra->data);
   
                    data = s->ra->data;
 
                    memcpy (data, fh->b->buff + slice_offset - buffer_offset, slice_size);
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {
                        change_endianness (data, slice_size, v->type);
                    }
                }
            }
            else if (break_dim == 0) 
            {
                int isize;
                uint64_t size_in_dset = 0, offset_in_dset = 0, offset_in_var = 0, write_offset;
    
                isize = offsets[0] + ldims[0];
                if (start_notime[0] >= offsets[0])
                {
                    // head is in
                    if (start_notime[0] < isize)
                    {
                        if (start_notime[0] + count_notime[0] > isize)
                        {
                            size_in_dset = isize - start_notime[0];
                        }
                        else
                        {
                            size_in_dset = count_notime[0];
                        }
                        offset_in_var = 0;
                        offset_in_dset = start_notime[0] - offsets[0];
                    }
                }
                else
                {
                    // middle is in
                    if (isize < start_notime[0] + count_notime[0])
                    {
                        size_in_dset = ldims[0];
                    }
                    else
                    {
                    // tail is in
                        size_in_dset = count_notime[0] + start_notime[0] - offsets[0];
                    }
                    offset_in_var = offsets[0] - start_notime[0];
                    offset_in_dset = 0;
                }
    
                slice_size = size_in_dset * datasize * size_unit;
                slice_offset = v->characteristics[start_idx + idx].payload_offset 
                                 + offset_in_dset * datasize * size_unit;
/*
                write_offset = offset_in_var * size_unit;
                for (i = 1; i < ri.ndim_notime; i++)
                {
                    write_offset *= count_notime[i];
                }
*/
                if (idx_check2)
                {
                    s->ra->data = malloc (slice_size);
                    assert (s->ra->data);
   
                    data = s->ra->data;
/*
                    memcpy (data + write_offset, fh->b->buff + slice_offset - buffer_offset, slice_size);
*/
                    memcpy (data, fh->b->buff + slice_offset - buffer_offset, slice_size);
/*
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {
                        change_endianness ((char *) data + write_offset, slice_size, v->type);
                    }
*/
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {
                        change_endianness ((char *) data, slice_size, v->type);
                    }
                }
            }
            else 
            {
                uint64_t stride_offset = 0;
                int isize;
                uint64_t size_in_dset[MAX_DIMS];
                uint64_t offset_in_dset[MAX_DIMS];
                uint64_t offset_in_var[MAX_DIMS];

                memset (size_in_dset, 0, MAX_DIMS * 8);
                memset (offset_in_dset, 0, MAX_DIMS * 8);
                memset (offset_in_var, 0, MAX_DIMS * 8);

                for (i = 0; i < ri.ndim_notime ; i++)
                {
                    isize = offsets[i] + ldims[i];
                    if (start_notime[i] >= offsets[i])
                    {
                        // head is in
                        if (start_notime[i] < isize)
                        {
                            if (start_notime[i] + count_notime[i] > isize)
                            {
                                size_in_dset[i] = isize - start_notime[i];
                            }
                            else
                            {
                                size_in_dset[i] = count_notime[i];
                            }
                            offset_in_dset[i] = start_notime[i] - offsets[i];
                            offset_in_var[i] = 0;
                        }
                        else
                        {
                        }
                    }
                    else
                    {
                        // middle is in
                        if (isize < start_notime[i] + count_notime[i])
                        {
                            size_in_dset[i] = ldims[i];
                        }
                        else
                        {
                            // tail is in
                            size_in_dset[i] = count_notime[i] + start_notime[i] - offsets[i];
                        }
                        offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start_notime[i];
                    }
                }
    
                datasize = 1;
                var_stride = 1;
    
                for (i = ri.ndim_notime - 1; i >= break_dim; i--)
                {
                    datasize *= size_in_dset[i];
                    dset_stride *= ldims[i];
                    var_stride *= count_notime[i];
                }
    
                uint64_t start_in_payload = 0, end_in_payload = 0, s_temp = 1;
                uint64_t var_offset = 0, dset_offset = 0;

                for (i = ri.ndim_notime - 1; i > -1; i--)
                {
                    start_in_payload += s_temp * offset_in_dset[i] * size_unit;
                    end_in_payload += s_temp * (offset_in_dset[i] + size_in_dset[i] - 1) * size_unit;
                    s_temp *= ldims[i];
                }

                slice_size = end_in_payload - start_in_payload + 1 * size_unit;
/*
                slice_size = datasize * size_unit; 
*/
                slice_offset =  v->characteristics[start_idx + idx].payload_offset
                                  + start_in_payload;
 
                for (i = 0; i < ri.ndim_notime ; i++)
                {
                    offset_in_dset[i] = 0;
                }
    
                for (i = 0; i < ri.ndim_notime ; i++)
                {
                    var_offset = offset_in_var[i] + var_offset * count_notime[i];
                    dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
                }

                if (idx_check2)
                {
/*
                    s->ra->data = malloc (slice_size);
*/
/*
int rank;
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
printf ("[%d: %d] s->ra->size = %llu, datasize = %llu, slice_size = %llu\n", rank, r->rank, s->ra->size, datasize, slice_size);
*/
                    s->ra->data = malloc (s->ra->size);

                    assert (s->ra->data);

                    data = s->ra->data;
/*
                    memcpy (data, fh->b->buff + slice_offset - buffer_offset, slice_size);


                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {   
                        change_endianness ((char *) data, slice_size, v->type);
                    }
*/

/*
printf ("copy_data: %7.5g %7.5g %7.5g %7.5g %7.5g %7.5g %7.5g %7.5g %7.5g %7.5g\n"
       , * (double *) (fh->b->buff + slice_offset - buffer_offset)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 1)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 2)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 3)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 4)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 5)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 6)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 7)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 8)
       , * ((double *) (fh->b->buff + slice_offset - buffer_offset) + 9)
       );
*/

/*
printf ("var_stride = %lu, dset_stride = %lu, var_offset = %lu, dset_offset = %lu, datasize = %lu, size_unit = %d\n", var_stride, dset_stride, var_offset, dset_offset, datasize, size_unit); 
*/
                    copy_data (data
                              ,fh->b->buff + slice_offset - buffer_offset
                              ,0
                              ,break_dim
                              ,size_in_dset
                              ,ldims
                              ,count_notime
                              ,var_stride
                              ,dset_stride
                              ,var_offset
                              ,dset_offset
                              ,datasize
                              ,size_unit 
                              );

//printf ("after copy_data = %7.5f %7.5f %7.5f %7.5f\n", * (double *) data, * ((double *) data + 1), * ((double *) data + 2), * ((double *) data + 3));
/*
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
if (rank != 0)
while (1);
*/
/*
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {   
                        change_endianness ((char *) data, slice_size, v->type);
                    }
*/

                }
            }
        }  // for idx ... loop over pgs

/*
        // shift target pointer for next read in
        data = (char *)data + (items_read * size_unit);
*/
        if (isTimeless (ri.file_tdim))
            break;
    } // for t

    freeReadInfo (&ri);
}

void adios_read_bp_staged1_read_chunk (ADIOS_GROUP * gp, int file_idx, uint64_t chunk_offset, uint64_t size)
{
    struct BP_GROUP * gh;
    BP_FILE * fh;
    MPI_File * sfh;
    MPI_Status status;
    int has_subfile;

struct timeval t0, t1;
gettimeofday (&t0, NULL);

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;

    has_subfile = fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE;

    bp_realloc_aligned(fh->b, size);
    fh->b->offset = 0;

    if (has_subfile)
    {
        sfh = get_BP_file_handle (fh->sfh, file_idx);

        if (!sfh)                                                                         
        {                                                                                   
            int err;                                                                        
            char * ch, * name_no_path, * name;                                              
            struct BP_file_handle * new_h =                                                
                  (struct BP_file_handle *) malloc (sizeof (struct BP_file_handle));

            new_h->file_index = file_idx;      
            new_h->next = 0;                                                                
            if (ch = strrchr (fh->fname, '/'))                                              
            {                                                                               
                name_no_path = malloc (strlen (ch + 1) + 1);                                
                strcpy (name_no_path, ch + 1);                                              
            }                                                                               
            else                                                                            
            {                                                                               
                name_no_path = malloc (strlen (fh->fname) + 1);                             
                strcpy (name_no_path, fh->fname);                                           
            }

            name = malloc (strlen (fh->fname) + 5 + strlen (name_no_path) + 1 + 10 + 1);
            sprintf (name, "%s.dir/%s.%d", fh->fname, name_no_path, new_h->file_index);

            err = MPI_File_open (MPI_COMM_SELF                                              
                                ,name                                                       
                                ,MPI_MODE_RDONLY                                            
                                ,(MPI_Info)MPI_INFO_NULL                                    
                                ,&new_h->fh                                                 
                                );                                                          
            if (err != MPI_SUCCESS)                                                          
            {
                fprintf (stderr, "can not open file %S\n", name);
                return;
            }
                                                                                
            add_BP_file_handle (&fh->sfh, new_h);
            sfh = &new_h->fh;

            free (name_no_path);
            free (name);
        }

        MPI_File_seek (*sfh
                      ,(MPI_Offset)chunk_offset
                      ,MPI_SEEK_SET
                      );
        MPI_File_read (*sfh
                      ,fh->b->buff
                      ,size
                      ,MPI_BYTE
                      ,&status
                      );
    }
    else
    {
        MPI_File_seek (fh->mpi_fh
                      ,(MPI_Offset)chunk_offset
                      ,MPI_SEEK_SET
                      );
        MPI_File_read (fh->mpi_fh
                      ,fh->b->buff
                      ,size
                      ,MPI_BYTE
                      ,&status
                      );
    }

    fh->b->offset = 0;

    gettimeofday (&t1, NULL);

//    printf ("read chunk time = %f \n", t1.tv_sec - t0.tv_sec + (double)(t1.tv_usec - t0.tv_usec)/1000000 );

}

/* The read request link list needs to be sorted beforehand
*/
static void do_read (struct proc_struct * p)
{
    void * data = 0;
    int i, counter = 0;
    int file_idx;
    uint64_t offset, payload_size;
#if DEBUG
    struct timeval t0;
    struct timeval t1;
    double t2, t3, t4, t5;
#endif
    candidate_reader * s = p->split_read_request_list1, * f_start = s, * f_end = s;
    candidate_reader * o_start = s, * o_end = s, * o_prev_end = 0, * parent = 0;

#if DEBUG
    gettimeofday (&t0, NULL);
    t2 = MPI_Wtime();
#endif

    while (f_start)
    {
        f_end = f_start;

        // Find a set of reqeusts that fall into the same file, i.e., [f_start, f_end)
        while (f_end && f_end->ra->file_idx == f_start->ra->file_idx)
        {
            f_end = f_end->next;
        }

        o_start = f_start;
        o_end = f_start;
        o_prev_end = 0;

        if (o_start->ra->file_idx >= p->f1 && o_start->ra->file_idx <= p->f2)
        {
            while (o_start != f_end)         
            {
                // Find a set of requests that fall into chunk size, i.e., [o_start, o_end)
                while (o_end && o_end != f_end && o_end->ra->offset - o_start->ra->offset <= p->chunk_size)
                {
                    o_prev_end = o_end;
                    o_end = o_end->next;
                }

                // Calculate the var payload size of the last request
                getDataAddress (p->gp, o_prev_end->ra->varid, o_prev_end->ra->start, o_prev_end->ra->count, &file_idx, &offset, &payload_size);
#if DEBUG
                //printf ("o_start.offset = %llu\n", o_start->ra->offset);
                //printf ("o_prev_end.offset = %llu\n", o_prev_end->ra->offset);

                t4 = MPI_Wtime ();
#endif
                // read a chunk from file into internal buffer
                adios_read_bp_staged1_read_chunk (p->gp, o_start->ra->file_idx, o_start->ra->offset, o_prev_end->ra->offset - o_start->ra->offset + payload_size); 

#if DEBUG
                t5 = MPI_Wtime ();
                //    printf ("read chunk = %f \n", t5 - t4);
#endif
                s = o_start;
                do
                {
                    parent = s->ra->parent;
                    // copy data from internal buffer to user buffer
                    adios_read_bp_staged1_read_buffer (p->gp, o_start->ra->offset, parent, s);

                    s = s->next;
                } while (s != o_end);

                o_start = o_end;
                o_prev_end = 0;
            }
        }

        f_start = f_end;
    }

#if DEBUG
    gettimeofday (&t1, NULL);
    t3 = MPI_Wtime ();

//    printf ("while time = %f \n", t1.tv_sec - t0.tv_sec + (double)(t1.tv_usec - t0.tv_usec)/1000000 );
//    printf ("while time = %f \n", t3 - t2);
#endif
}

static void free_candidate_reader_list (candidate_reader * h)
{
    candidate_reader * n;
    while (h)
    {
        n = h->next;

        if (h->ra && h->ra->start)
        {
            free (h->ra->start);
        }

        if (h->ra && h->ra->count)
        {
            free (h->ra->count);
        }

        if (h->ra)
        {
            free (h->ra);
        }

        free (h);
        h = n;
    }
}

static void free_proc_struct (struct proc_struct * p)
{
    candidate_reader * h = p->local_read_request_list, * n;

    while (h)
    {
        n = h->next;

        if (h->ra && h->ra->start)
        {
            free (h->ra->start);
        }

        if (h->ra && h->ra->count)
        {
            free (h->ra->count);
        }

        /* DO NOT free if the request is from aggregator itself */
        if (h->ra && h->ra->data && h->rank != p->rank)
        {
            free (h->ra->data);
        }

        if (h->ra)
        {
            free (h->ra);
        }

        free (h);
        h = n;
    }

    /* Free both split_read_request_list1 and split_read_request_list2 */
    h = p->split_read_request_list1;
    while (h)
    {
        n = h->next;
    
        if (h->ra && h->ra->start)
        {
            free (h->ra->start);
        }
    
        if (h->ra && h->ra->count)
        {
            free (h->ra->count); 
        }

        /* h->ra->data should have already been freed after copy_buffer() */
        if (h->ra && h->ra->data)
        {
            free (h->ra->data);
        }

        if (h->ra)
        {
            free (h->ra);
        }

        free (h);
        h = n;
    }

    h = p->split_read_request_list2;
    while (h)
    {
        n = h->next;

        if (h->ra && h->ra->start)
        {
            free (h->ra->start);
        }

        if (h->ra && h->ra->count)
        {
            free (h->ra->count);
        }

        /* h->ra->data should have already been freed after copy_buffer() */
        if (h->ra && h->ra->data)
        {
            free (h->ra->data);
        }

        if (h->ra)
        {
            free (h->ra);
        }

        free (h);
        h = n;
    }

    if (p->aggregator_rank_array)
    {
        free (p->aggregator_rank_array);
    }

    if (p->b)
    {
        free (p->b);
    }
}

static void init_read (BP_FILE * fh)
{
    int thread_level, i, remain;
    int color1, color2;
    char * env_str;

    struct proc_struct * p = (struct proc_struct *) malloc (sizeof (struct proc_struct));
    assert (p);

    fh->priv = p;
    MPI_Comm_rank (fh->comm, &p->rank);
    MPI_Comm_size (fh->comm, &p->size);

    env_str = getenv ("num_aggregators");
    if (!env_str)
    {
        fprintf (stderr, "Environment variable \"num_aggregators\" hasn't been set.\n");
        exit(0);
    }

    p->num_aggregators = atoi (env_str);

    if (p->rank == 0)
    {
        printf ("%d aggregators are used.\n", p->num_aggregators);
    }

    env_str = getenv ("chunk_size");
    if (!env_str)
    {
        fprintf (stderr, "Environment variable \"chunk_size\" (in MB) hasn't been set.\n");
        exit(0);
    }

    p->chunk_size = 1024 * 1024 * atoi (env_str);
    p->gp = 0; // p->gp is set in gopen call
    p->groups = (p->num_aggregators > p->size || p->num_aggregators <= 0) ? p->size : p->num_aggregators;
    p->group_size = p->size / p->groups;
    remain = p->size - p->group_size * p->groups;

    p->aggregator_rank_array = (int *) malloc (p->groups * sizeof (int));
    for (i = 0; i < p->groups; i++)
    {
        if (remain == 0)
        {
            p->aggregator_rank_array[i] = p->group_size * i;
        }
        else
        {
            if (i < remain)
            {
                p->aggregator_rank_array[i] = (p->group_size + 1) * i;
            }
            else
            {
                p->aggregator_rank_array[i] = remain * (p->group_size + 1) + (i - remain) * p->group_size;
            }
        }
    }

    if (remain == 0)
    {
        color1 = p->rank / p->group_size;
        color2 = p->rank % p->group_size;
 
        p->aggregator_rank = color1 * p->group_size;
    }
    else
    {
        if (p->rank < (p->group_size + 1) * remain)
        {
            color1 = p->rank / (p->group_size + 1);
            color2 = p->rank % (p->group_size + 1);
        
            p->aggregator_rank = color1 * (p->group_size + 1);
            p->group_size++;
        }
        else
        {
            color1 = remain + (p->rank - (p->group_size + 1) * remain) / p->group_size;
            color2 = (p->rank - (p->group_size + 1) * remain) % p->group_size;
        
            p->aggregator_rank = remain * (p->group_size + 1) + (color1 - remain) * p->group_size;
        }
    }

    p->group = color1;

    MPI_Comm_split (fh->comm, color1, p->rank, &p->new_comm);
    MPI_Comm_split (fh->comm, color2, p->rank, &p->new_comm2);
    MPI_Comm_rank (p->new_comm, &p->new_rank);

    p->aggregator_new_rank = 0;
    p->group_comm = fh->comm;

    p->local_read_request_list = 0;
    p->split_read_request_list1 = 0;
    p->split_read_request_list2 = 0;
    p->ht = 0;
    p->b = 0;
    p->read_close_received = 0;
    p->group_close_received = 0;

    return;
}

static MPI_File * get_BP_file_handle(struct BP_file_handle * l, uint32_t file_index)
{
    if (!l)
        return 0;

    while (l)
    {
        if (l->file_index == file_index)
            return &l->fh;

        l = l->next;
    }

    return 0;
}

static void add_BP_file_handle (struct BP_file_handle ** l, struct BP_file_handle * n)
{
    if (!n)
        return;

    n->next = *l;
    *l = n;
}


static void close_all_BP_files (struct BP_file_handle * l)
{
    struct BP_file_handle * n;

    while (l)
    {
        n = l->next;

        MPI_File_close (&l->fh);
        free (l);

        l = n;
    }
}

/* Return 0: if file is little endian, 1 if file is big endian 
 * We know if it is different from the current system, so here
 * we determine the current endianness and report accordingly.
 */
static int adios_read_bp_staged1_get_endianness( uint32_t change_endianness )
{
   int LE = 0;
   int BE = !LE;
   int i = 1;
   char *p = (char *) &i;
   int current_endianness;
   if (p[0] == 1) // Lowest address contains the least significant byte
       current_endianness = LE;
   else
       current_endianness = BE;
    if (change_endianness == adios_flag_yes)
        return !current_endianness;
    else
        return current_endianness;
}

static int getNumSubfiles (const char * fname)
{
   char * dirname;
   int n = 0;
#include <sys/types.h>
#include <dirent.h>
    dirname = malloc (strlen (fname) + 5);
    sprintf (dirname, "%s%s", fname, ".dir");

    DIR * mydir = opendir(dirname);
    struct dirent * entry;

    n = 0;
    while((entry = (struct dirent *)readdir (mydir)))
    {
        n++;
    }

    // exclude . and ..
    n -= 2;
    closedir(mydir);
    free (dirname);

    return n;
}

int adios_read_bp_staged1_init (MPI_Comm comm) { return 0; }
int adios_read_bp_staged1_finalize () { return 0; }

static void broadcast_fh_buffer (BP_FILE * fh)
{
    struct bp_index_pg_struct_v1 * pgs_root = fh->pgs_root, * pg;
    struct adios_index_var_struct_v1 * vars_root = fh->vars_root, * v;
    struct adios_index_attribute_struct_v1 * attrs_root = fh->attrs_root;
    char * buffer;
    uint64_t buffer_size, buffer_offset = 0;
    int i, j, timedim;
    uint16_t len;
    uint8_t flag;
    struct proc_struct * p = (struct proc_struct *) fh->priv;


    bp_realloc_aligned (fh->b, 0);
/*
    buffer = fh->b->buff;
*/
    buffer = 0;
    buffer_size = 0;
    buffer_offset = 0;

    if (isAggregator (p))
    {
        _buffer_write (&buffer, &buffer_size, &buffer_offset, &p->num_aggregators, 4); // n_sf
        _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->gvar_h->group_count, 2); //group_count 
        _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->mfooter.pgs_count, 8); //vars_count 
        _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->mfooter.vars_count, 4); //vars_count 
        _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->mfooter.attrs_count, 4); //attrs_count 

        for (i = 0; i < fh->gvar_h->group_count; i++)
        {
            len = strlen (fh->gvar_h->namelist[i]);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset, fh->gvar_h->namelist[i], len); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->gvar_h->var_counts_per_group[i], 2); // var_counts_per_group
        }

        _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->gattr_h->group_count, 2); //group_count 
        for (i = 0; i < fh->gattr_h->group_count; i++)
        {
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->gattr_h->attr_counts_per_group[i], 2); // attr_counts_per_group
        }

        for (i = 0; i < fh->mfooter.vars_count; i++)
        {
            len = strlen (fh->gvar_h->var_namelist[i]);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset, fh->gvar_h->var_namelist[i], len); // namelist
        }

        for (i = 0; i < fh->mfooter.attrs_count; i++)
        {
            len = strlen (fh->gattr_h->attr_namelist[i]);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset, fh->gattr_h->attr_namelist[i], len); // namelist
        }   

        _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->tidx_start, 4); // tidx_start
        _buffer_write (&buffer, &buffer_size, &buffer_offset, &fh->tidx_stop, 4); // tidx_start

        pgs_root = fh->pgs_root;
        while (pgs_root)
        {
            len = strlen (pgs_root->group_name);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2); // group_name len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, pgs_root->group_name, len); // group_name

            flag = (pgs_root->adios_host_language_fortran == adios_flag_yes) ? 'y' : 'n';
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &flag, 1);

            _buffer_write (&buffer, &buffer_size, &buffer_offset, &pgs_root->process_id, 4);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &pgs_root->time_index, 4);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &pgs_root->offset_in_file, 8);

            pgs_root = pgs_root->next;
        }

        vars_root = fh->vars_root;
        while (vars_root)
        {
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &vars_root->id, 4); // id

            len = strlen (vars_root->group_name);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2); // group_name len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, vars_root->group_name, len); // group_name

            len = strlen (vars_root->var_name);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2); // var_name len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, vars_root->var_name, len); // var_name

            len = strlen (vars_root->var_path);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2); // var_path len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, vars_root->var_path, len); // var_path
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &vars_root->type, 4); // type

            ADIOS_VARINFO * vi = _inq_var_byid (fh, vars_root->id);
            assert (vi);

            _buffer_write (&buffer, &buffer_size, &buffer_offset, &vi->ndim, 4); // ndim
            if (vi->ndim)
            {
                _buffer_write (&buffer, &buffer_size, &buffer_offset, vi->dims, vi->ndim * 8); // dims
            }

            _buffer_write (&buffer, &buffer_size, &buffer_offset, &vi->timedim, 4); // ndim
         
            len = (vars_root->characteristics[0].value == 0) ? 0 : bp_get_type_size (vars_root->type, vars_root->characteristics[0].value);
            if (vars_root->type == adios_string)
            {
                len--;
            }
/*
int rank;
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
if (rank == 0)
fprintf (stderr, "bc %s bo 1 = %llu, bo 2 = %llu, len = %d\n", vars_root->var_name, bo, buffer_offset, len);
*/
            _buffer_write (&buffer, &buffer_size, &buffer_offset, &len, 2);

            if (len)
            {
//                _buffer_write (&buffer, &buffer_size, &buffer_offset, vi->value, len); // ndim
                _buffer_write (&buffer, &buffer_size, &buffer_offset, vars_root->characteristics[0].value, len); // ndim

            }

            adios_read_bp_staged1_free_varinfo (vi);

            vars_root = vars_root->next;
        }
    }

    MPI_Bcast (&buffer_offset, 8, MPI_BYTE, 0, p->new_comm);
    if (!isAggregator (p))
    {
/*
        bp_realloc_aligned (fh->b, buffer_offset);
        assert (fh->b->buff);

        buffer = fh->b->buff;
*/
        buffer = malloc (buffer_offset);
        assert (buffer);
    }

    MPI_Bcast (buffer, buffer_offset, MPI_BYTE, 0, p->new_comm);

    if (!isAggregator (p))
    {
        uint16_t len, group_count, var_counts_per_group;

        buffer_offset = 0;

        _buffer_read (buffer, &buffer_offset, &p->num_aggregators, 4);

        fh->gvar_h = (struct BP_GROUP_VAR *) malloc (sizeof (struct BP_GROUP_VAR));
        assert (fh->gvar_h);
        fh->gvar_h->time_index = 0;
        fh->gvar_h->var_offsets = 0;
        fh->gvar_h->pg_offsets = 0;

        _buffer_read (buffer, &buffer_offset, &fh->gvar_h->group_count, 2); //group_count 
        _buffer_read (buffer, &buffer_offset, &fh->mfooter.pgs_count, 8); //pgs_count 
        _buffer_read (buffer, &buffer_offset, &fh->mfooter.vars_count, 4); //vars_count 
        _buffer_read (buffer, &buffer_offset, &fh->mfooter.attrs_count, 4); //attrs_count 

        fh->gvar_h->namelist = (char **) malloc (fh->gvar_h->group_count * sizeof (char *));
        fh->gvar_h->var_counts_per_group = (uint16_t *) malloc (fh->gvar_h->group_count * 2);

        for (i = 0; i < fh->gvar_h->group_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &len, 2); // len
            fh->gvar_h->namelist[i] = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, fh->gvar_h->namelist[i], len); // namelist
            fh->gvar_h->namelist[i][len] = '\0';

            _buffer_read (buffer, &buffer_offset, &fh->gvar_h->var_counts_per_group[i], 2); // var_counts_per_group
        }

        fh->gattr_h = (struct BP_GROUP_ATTR *) malloc (sizeof (struct BP_GROUP_ATTR));
        assert (fh->gattr_h);
        fh->gattr_h->attr_offsets = 0;

        _buffer_read (buffer, &buffer_offset, &fh->gattr_h->group_count, 2); //group_count 

        fh->gattr_h->attr_counts_per_group = (uint16_t *) malloc (fh->gattr_h->group_count * 2);
        fh->gattr_h->namelist = fh->gvar_h->namelist;

        for (i = 0; i < fh->gattr_h->group_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &fh->gattr_h->attr_counts_per_group[i], 2); // attr_counts_per_group
        }

        fh->gvar_h->var_namelist = (char **) malloc (fh->mfooter.vars_count * sizeof (char *));
        for (i = 0; i < fh->mfooter.vars_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &len, 2); // len
            fh->gvar_h->var_namelist[i] = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, fh->gvar_h->var_namelist[i], len); // namelist
            fh->gvar_h->var_namelist[i][len] = '\0';
        }

        fh->gattr_h->attr_namelist = (char **) malloc (fh->mfooter.attrs_count * sizeof (char *));
        for (i = 0; i < fh->mfooter.attrs_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &len, 2); // len
            fh->gattr_h->attr_namelist[i] = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, fh->gattr_h->attr_namelist[i], len); // namelist
            fh->gattr_h->attr_namelist[i][len] = '\0';
        }

        _buffer_read (buffer, &buffer_offset, &fh->tidx_start, 4);
        _buffer_read (buffer, &buffer_offset, &fh->tidx_stop, 4);

        pgs_root = 0;
        int pgs_count = fh->mfooter.pgs_count;

        for (i = 0; i < pgs_count; i++)
        {
            pg = (struct bp_index_pg_struct_v1 *) malloc (sizeof (struct bp_index_pg_struct_v1));
            assert (pg);

            _buffer_read (buffer, &buffer_offset, &len, 2);

            pg->group_name = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, pg->group_name, len);
            pg->group_name[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &flag, 1);
            pg->adios_host_language_fortran = (flag == 'y') ? adios_flag_yes : adios_flag_no;

            _buffer_read (buffer, &buffer_offset, &pg->process_id, 4);
            _buffer_read (buffer, &buffer_offset, &pg->time_index, 4);
            _buffer_read (buffer, &buffer_offset, &pg->offset_in_file, 8);

            pg->time_index_name = 0;
            pg->next = 0;

            if (!pgs_root)
            {
                pgs_root = pg;
                fh->pgs_root = pg;
            }
            else
            {
                pgs_root->next = pg;
                pgs_root = pg;
            }
        }

        vars_root = 0;
        for (i = 0; i < fh->mfooter.vars_count; i++)
        {
            v = (struct adios_index_var_struct_v1 *) malloc (sizeof (struct adios_index_var_struct_v1));
            assert (v);
uint64_t bo = buffer_offset;
            _buffer_read (buffer, &buffer_offset, &v->id, 4);

            _buffer_read (buffer, &buffer_offset, &len, 2);
            v->group_name = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, v->group_name, len);
            v->group_name[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &len, 2);
            v->var_name = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, v->var_name, len);
            v->var_name[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &len, 2);
            v->var_path = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, v->var_path, len);
            v->var_path[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &v->type, 4);

            v->characteristics_count = 1; 
            v->characteristics = (struct adios_index_characteristic_struct_v1 *)
                                     malloc (sizeof (struct adios_index_characteristic_struct_v1));
            assert (v->characteristics);
 
            v->characteristics->stats = 0;
            _buffer_read (buffer, &buffer_offset, &(v->characteristics->dims.count), 4);

            v->characteristics->dims.dims = 0;

int rank;
MPI_Comm_rank (MPI_COMM_WORLD, &rank);

            if (v->characteristics->dims.count)
            {
                v->characteristics->dims.dims = (uint64_t *) malloc (v->characteristics->dims.count * 3 * 8);
                assert (v->characteristics->dims.dims);

                uint64_t * gdims = (uint64_t *) malloc (v->characteristics->dims.count * 8);
                assert (gdims);

                _buffer_read (buffer, &buffer_offset, gdims, v->characteristics->dims.count * 8);
                for (j = 0; j < v->characteristics->dims.count; j++)
                {
                    v->characteristics->dims.dims[j * 3] = gdims[j];
                    v->characteristics->dims.dims[j * 3 + 1] = gdims[j];
                    v->characteristics->dims.dims[j * 3 + 2] = 0;
                }
                free (gdims);
            }
            _buffer_read (buffer, &buffer_offset, &timedim, 4);
/*
if (rank == 1)
fprintf (stderr, "bc 1 v->id = %d, bo 1 = %llu, bo 2 = %llu, v->var_name = %s\n", v->id, bo, buffer_offset, v->var_name);
*/
            _buffer_read (buffer, &buffer_offset, &len, 2);
    
            if (len)
            {
                if (v->type == adios_string)
                {
                    v->characteristics->value = malloc (len + 1);
                }
                else
                {
                    v->characteristics->value = malloc (len);
                }

                _buffer_read (buffer, &buffer_offset, v->characteristics->value, len);

                if (v->type == adios_string)
                {
                    ((char * )(v->characteristics->value))[len] = '\0';
                }
            }
            else
            {
                v->characteristics->value = 0;
            }
            v->next = 0;

            if (!vars_root)
            {
                vars_root = v;
                fh->vars_root = v;
            }
            else
            {
                vars_root->next = v;
                vars_root = v;
            }
        }
    }

    if (buffer)
    {
        free (buffer);
    }
}

ADIOS_FILE * adios_read_bp_staged1_fopen (const char * fname, MPI_Comm comm)
{
    int i, rank, remain;
    BP_FILE * fh;
    ADIOS_FILE * fp;
    uint64_t header_size;
    struct proc_struct * p;
    char * env_str;

    adios_errno = 0;

    fh = (BP_FILE *) malloc (sizeof (BP_FILE));
    assert (fh);

    fh->fname = (fname ? strdup (fname) : 0L);
    //FIXME
    fh->mpi_fh = 0;
    fh->sfh = 0;
    fh->comm = comm;
    fh->gvar_h = 0;
    fh->pgs_root = 0;
    fh->vars_root = 0;
    fh->attrs_root = 0;
    fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    assert (fh->b);

    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    assert (fp);

    adios_buffer_struct_init (fh->b);

    init_read (fh);

    p = (struct proc_struct *) fh->priv;

    if (isAggregator (p))
    {
        if (bp_read_open (fname, p->new_comm2, fh))
            return NULL;

        if (p->rank == 0)
        {
            if (bp_read_minifooter (fh))
            {
                return NULL;
            }
        }

        MPI_Bcast (&fh->mfooter, sizeof (struct bp_minifooter), MPI_BYTE, 0, p->new_comm2);

        header_size = fh->mfooter.file_size - fh->mfooter.pgs_index_offset;

        if (p->rank != 0)
        {
            if (!fh->b->buff)
            {
                bp_alloc_aligned (fh->b, header_size);
                if (!fh->b->buff)
                   return NULL;
                memset (fh->b->buff, 0, header_size);
                fh->b->offset = 0;
            }
        }
    
        MPI_Bcast (fh->b->buff, fh->mfooter.file_size - fh->mfooter.pgs_index_offset, MPI_BYTE, 0, p->new_comm2);

        /* Everyone parses the index on its own */
        bp_parse_pgs (fh);
        bp_parse_vars (fh);
        bp_parse_attrs (fh);

        /* find out how many subfiles */
        p->n_total_sf = get_num_subfiles (fh);
        p->n_my_sf = p->n_total_sf / p->groups;
        remain = p->n_total_sf - p->n_my_sf * p->groups;

        if (remain == 0)
        {   
            p->f1 = p->group * p->n_my_sf;
            p->f2 = (p->group + 1) * p->n_my_sf - 1;
        }  
        else
        {
            if (p->group < remain)
            {   
                p->f1 = p->group * (p->n_my_sf + 1);
                p->f2 = p->f1 + p->n_my_sf;

                p->n_my_sf++;
            }
            else
            {   
                p->f1 = remain * (p->n_my_sf + 1) + (p->group - remain) * p->n_my_sf;
                p->f2 = p->f1 + p->n_my_sf - 1;
            }
        }

//printf ("rank: %d, f1 = %d, f2 = %d\n", p->rank, p->f1, p->f2);

        /* fill out ADIOS_FILE struct */
        fp->vars_count = fh->mfooter.vars_count;
        fp->attrs_count = fh->mfooter.attrs_count;
        fp->tidx_start = fh->tidx_start;
        fp->ntimesteps = fh->tidx_stop - fh->tidx_start + 1;
        fp->file_size = fh->mfooter.file_size;
        fp->version = fh->mfooter.version & ADIOS_VERSION_NUM_MASK;
        fp->endianness = adios_read_bp_staged1_get_endianness (fh->mfooter.change_endianness);
    }

    broadcast_fh_buffer (fh);

    fp->fh = (uint64_t) fh;
    fp->groups_count = fh->gvar_h->group_count;
    alloc_namelist (&fp->group_namelist, fp->groups_count);

    for (i = 0;i < fp->groups_count; i++)
    {
        if (!fp->group_namelist[i])
        {
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_fopen()", fp->groups_count);
            adios_read_bp_fclose (fp);
            return NULL;
        }
        else
        {
            strcpy (fp->group_namelist[i], fh->gvar_h->namelist[i]);
        }
    }

    return fp;
}

/* This function can be called if user places 
   the wrong sequences of dims for a var 
*/   
void adios_read_bp_staged1_reset_dimension_order (ADIOS_FILE *fp, int is_fortran)
{
    BP_FILE * fh = (BP_FILE *)(fp->fh);
    struct bp_index_pg_struct_v1 ** root = &(fh->pgs_root);
    struct bp_minifooter * mh = &(fh->mfooter);
    uint64_t i;

    for (i = 0; i < mh->pgs_count; i++) {
        is_fortran ? ((*root)->adios_host_language_fortran = adios_flag_yes) 
               : ((*root)->adios_host_language_fortran = adios_flag_no);
        root = &(*root)->next;
    }
}

int adios_read_bp_staged1_fclose (ADIOS_FILE *fp) 
{
    BP_FILE * fh = (BP_FILE *) fp->fh;
    struct BP_GROUP_VAR * gh = fh->gvar_h;
    struct BP_GROUP_ATTR * ah = fh->gattr_h;
    struct adios_index_var_struct_v1 * vars_root = fh->vars_root, *vr;
    struct adios_index_attribute_struct_v1 * attrs_root = fh->attrs_root, *ar;
    struct bp_index_pg_struct_v1 * pgs_root = fh->pgs_root, *pr;
    int i,j;

    MPI_File mpi_fh = fh->mpi_fh;

    adios_errno = 0;
    if (fh->mpi_fh) 
        MPI_File_close (&mpi_fh);

    if (fh->sfh)
        close_all_BP_files (fh->sfh);

    if (fh->b) {
        adios_posix_close_internal (fh->b);
        free(fh->b);
    }

    /* Free variable structures */
    /* alloc in bp_utils.c: bp_parse_vars() */
    while (vars_root) {
        vr = vars_root;
        vars_root = vars_root->next;
        for (j = 0; j < vr->characteristics_count; j++) {
            // alloc in bp_utils.c:bp_parse_characteristics() <- bp_get_characteristics_data()

            if (vr->characteristics[j].dims.dims)
                free (vr->characteristics[j].dims.dims);

            if (vr->characteristics[j].value)
                free (vr->characteristics[j].value);
            // NCSU - Clearing up statistics
            if (vr->characteristics[j].stats)
            {
                uint8_t k = 0, idx = 0;
                uint8_t i = 0, count = adios_get_stat_set_count(vr->type);

                while (vr->characteristics[j].bitmap >> k)
                {
                    if ((vr->characteristics[j].bitmap >> k) & 1)
                    {
                        for (i = 0; i < count; i ++)
                        {
                            if (k == adios_statistic_hist)
                            {
                                struct adios_index_characteristics_hist_struct * hist = (struct adios_index_characteristics_hist_struct *) (vr->characteristics [j].stats[i][idx].data);
                                free (hist->breaks);
                                free (hist->frequencies);
                                free (hist);
                            }
                            else
                            free (vr->characteristics[j].stats [i][idx].data);
                        }
                        idx ++;
                    }
                    k ++;
                }

                for (i = 0; i < count; i ++)
                    free (vr->characteristics[j].stats [i]);

                free (vr->characteristics[j].stats);
                vr->characteristics[j].stats = 0;
            }

        }

        // NCSU ALACRITY-ADIOS - Clear transform metadata
        adios_transform_clear_transform_characteristic(&vr->characteristics[j]);

        if (vr->characteristics) 
            free (vr->characteristics);
        if (vr->group_name) 
            free (vr->group_name);
        if (vr->var_name) 
            free (vr->var_name);
        if (vr->var_path) 
            free (vr->var_path);
        free(vr);
    }

    /* Free attributes structures */
    /* alloc in bp_utils.c bp_parse_attrs() */

    while (attrs_root) {
        ar = attrs_root;
        attrs_root = attrs_root->next;
        for (j = 0; j < ar->characteristics_count; j++) {
            if (ar->characteristics[j].value)
                free (ar->characteristics[j].value);
        }
        if (ar->characteristics) 
            free (ar->characteristics);
        if (ar->group_name) 
            free (ar->group_name);
        if (ar->attr_name) 
            free (ar->attr_name);
        if (ar->attr_path) 
            free (ar->attr_path);
        free(ar);
    }


    /* Free process group structures */
    /* alloc in bp_utils.c bp_parse_pgs() first loop */
    //printf ("pgs: %d\n", fh->mfooter.pgs_count);
    while (pgs_root) {
        pr = pgs_root;
        pgs_root = pgs_root->next;
        //printf("%d\tpg pid=%d addr=%x next=%x\n",i, pr->process_id, pr, pr->next);
        if (pr->group_name)
            free(pr->group_name);
        if (pr->time_index_name)
            free(pr->time_index_name);
        free(pr);
    }

    /* Free variable structures in BP_GROUP_VAR */
    if (gh) {
        for (j=0;j<2;j++) { 
            for (i=0;i<gh->group_count;i++) {
                if (gh->time_index && gh->time_index[j] && gh->time_index[j][i])
                    free(gh->time_index[j][i]);
            }
            if (gh->time_index && gh->time_index[j])
                free(gh->time_index[j]);
        }
        
        if (gh->time_index)
            free (gh->time_index);
    
        for (i=0;i<gh->group_count;i++) { 
            if (gh->namelist[i])
                free(gh->namelist[i]);
        }
        if (gh->namelist)
            free (gh->namelist);

        for (i=0;i<fh->mfooter.vars_count;i++) {
            if (gh->var_namelist[i])
                free(gh->var_namelist[i]);
            if (gh->var_offsets && gh->var_offsets[i]) 
                free(gh->var_offsets[i]);
        }
        if (gh->var_namelist)
            free (gh->var_namelist);

        if (gh->var_offsets) 
            free(gh->var_offsets);

        if (gh->var_counts_per_group)
            free(gh->var_counts_per_group);

        if (gh->pg_offsets)
            free (gh->pg_offsets);

        free (gh);
    }

    /* Free attribute structures in BP_GROUP_ATTR */
    if (ah) {
        for (i = 0; i < fh->mfooter.attrs_count; i++) {
            if (ah->attr_offsets && ah->attr_offsets[i]) 
                free(ah->attr_offsets[i]);
            if (ah->attr_namelist && ah->attr_namelist[i]) 
                free(ah->attr_namelist[i]);
        }
        if (ah->attr_offsets)
            free(ah->attr_offsets);

        if (ah->attr_namelist)
            free(ah->attr_namelist);
        if (ah->attr_counts_per_group) 
            free(ah->attr_counts_per_group);

        free(ah);
    }

    if (fh->fname)
        free (fh->fname);
    if (fh->priv)
        free (fh->priv);    
    if (fh)
        free (fh);    

    free_namelist ((fp->group_namelist),fp->groups_count);
    free(fp);
    return 0;
}


ADIOS_GROUP * adios_read_bp_staged1_gopen (ADIOS_FILE *fp, const char * grpname)
{
    BP_FILE * fh = (BP_FILE *) fp->fh;
    int grpid, rank, nproc; 
    ADIOS_GROUP * gp;
    struct proc_struct * p = (struct proc_struct *) fh->priv;

    adios_errno = 0;
    for (grpid=0;grpid<(fh->gvar_h->group_count);grpid++) {
        if (!strcmp(fh->gvar_h->namelist[grpid], grpname))
            break; 
    }

    if (grpid >= fh->gvar_h->group_count) {
        adios_error ( err_invalid_group, "Invalid group name %s", grpname);
        return NULL;
    }

    gp = adios_read_bp_staged1_gopen_byid (fp, grpid);

    p->gp = gp; 

    return gp;
}

ADIOS_GROUP * adios_read_bp_staged1_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    BP_FILE * fh = (BP_FILE *) fp->fh;
    struct BP_GROUP * gh;
    ADIOS_GROUP * gp;
    int i, offset;

    adios_errno = 0;
    if (grpid < 0 || grpid >= fh->gvar_h->group_count) {
        adios_error ( err_invalid_group, "Invalid group index %d", grpid);
        return NULL;
    }

    gh = (struct BP_GROUP *) malloc(sizeof(struct BP_GROUP));
    if (!gh) {
        adios_error ( err_no_memory, "Could not allocate memory for group info");
        return NULL;
    }

    gp = (ADIOS_GROUP *) malloc(sizeof(ADIOS_GROUP));
    if (!gp) {
        adios_error ( err_no_memory, "Could not allocate memory for group info");
        free(gh);
        return NULL;
    }

    /* set offset index of variables (which is a long list of all vars in all groups) in this group */
    offset = 0;
    for (i=0; i<grpid; i++)
        offset += fh->gvar_h->var_counts_per_group[i];


    /* gh->vars_root will point to the list of vars in this group */
    gh->vars_root = fh->vars_root; 
    for (i=0; i<offset; i++)
        gh->vars_root = gh->vars_root->next;

    gh->group_id = grpid;
    gh->vars_offset = offset;
    gh->vars_count = fh->gvar_h->var_counts_per_group[grpid];

    /* set offset of attributes in this group */
    offset = 0;
    for(i=0;i<grpid;i++)
        offset += fh->gattr_h->attr_counts_per_group[i];

    /* gh->attrs_root will point to the list of vars in this group */
    gh->attrs_root = fh->attrs_root; 
    for (i=0; i<offset; i++)
        gh->attrs_root = gh->attrs_root->next;

    gh->attrs_offset = offset;
    gh->attrs_count = fh->gattr_h->attr_counts_per_group[grpid];

    gh->fh = fh; 

    /* fill out ADIOS_GROUP struct */
    gp->grpid = grpid;
    gp->gh = (uint64_t) gh;
    gp->fp = fp;
    gp->vars_count = gh->vars_count;
    gp->attrs_count = gh->attrs_count;

    offset = gh->vars_offset;
    alloc_namelist (&(gp->var_namelist), gp->vars_count);
    for (i=0;i<gp->vars_count;i++) {
        if (!gp->var_namelist[i]) { 
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_gopen()", gp->vars_count);
            adios_read_bp_staged1_gclose(gp);
            return NULL;
        }
        else
            strcpy(gp->var_namelist[i], gh->fh->gvar_h->var_namelist[i+offset]);
    }

    offset = gh->attrs_offset;
    alloc_namelist (&(gp->attr_namelist), gp->attrs_count);
    for (i=0;i<gp->attrs_count;i++) {
        if (!gp->attr_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_gopen()", gp->vars_count);
            adios_read_bp_staged1_gclose(gp);
            return NULL;
        }
        else {
            strcpy(gp->attr_namelist[i], gh->fh->gattr_h->attr_namelist[i+offset]);
        }
    }

    return gp;
}
                   
int adios_read_bp_staged1_gclose (ADIOS_GROUP * gp)
{
    struct BP_GROUP * gh = (struct BP_GROUP *) gp->gh;
    BP_FILE * fh = gh->fh;
    struct proc_struct * p = (struct proc_struct *) fh->priv;
    candidate_reader * h = p->local_read_request_list, * t = 0;
    int i, type, count, varid, ndims, total_size, size = calc_data_size (p);
    void * buf;

    p->b = malloc (size);
    assert (p->b);

    buf = p->b;

    /*  1. All rank packs their read requests.
        2. An aggregator gets packed requests from all the processors in the group
        3. Subsequently, all aggregators share with each other the requests it
           received previously.
        Note: the two-round communication is designed to avoid scalability issue so
        that MPI_Allgatherv() call is done on a smaller scale.
     */    
    while (h)
    {
        buffer_write (&buf, &p->rank, 4);

        varid = h->ra->varid;
        ndims = h->ra->ndims;

        buffer_write (&buf, &varid, 4);
        buffer_write (&buf, &ndims, 4);
        buffer_write (&buf, h->ra->start, ndims * 8);
        buffer_write (&buf, h->ra->count, ndims * 8);
        buffer_write (&buf, &h->ra->size, 8);
        buffer_write (&buf, &h->ra->file_idx, 4);
        buffer_write (&buf, &h->ra->offset, 8);

        h = h->next;
    }

    int * sizes = malloc (p->group_size * 4);
    int * offsets = malloc (p->group_size * 4);
    void * recv_buffer;

    MPI_Gather (&size, 1, MPI_INT
               ,sizes, 1, MPI_INT
               ,p->aggregator_new_rank, p->new_comm
               );

    if (isAggregator (p))
    {
        total_size = 0;
        offsets[0] = 0;

        for (i = 0; i < p->group_size; i++)
        {
            total_size += sizes[i];
            if (i > 0)
            {
                offsets[i] = offsets[i - 1] + sizes[i - 1];
            }
        }

        recv_buffer = malloc (total_size);
        assert (recv_buffer);
    }

    MPI_Gatherv (p->b, size, MPI_BYTE
                ,recv_buffer, sizes, offsets
                ,MPI_BYTE, p->aggregator_new_rank, p->new_comm
                );


    if (isAggregator (p))
    {
#if 1
        struct timeval t3, t31, t32, t33, t4;
        gettimeofday (&t3, NULL);
#endif

        t = p->local_read_request_list;
        if (!t)
        {
            p->local_read_request_list = parse_buffer (p, recv_buffer, total_size);
        }
        else
        {
            while (t->next)
            {
                t = t->next;
            }

            t->next = parse_buffer (p, recv_buffer, total_size);
        }

#if 1
        gettimeofday (&t31, NULL);
#endif
        shuffle_request (p);

#if 1
        gettimeofday (&t32, NULL);
#endif

        free (recv_buffer);

#if DEBUG
        gettimeofday (&t33, NULL);
#endif

        process_read_requests (p);
#if DEBUG
        gettimeofday (&t4, NULL);

        printf ("[%3d] process read request = %f (%f,%f,%f,%f)\n", p->rank, t4.tv_sec - t3.tv_sec + (double)(t4.tv_usec - t3.tv_usec)/1000000
                                                                    , t31.tv_sec - t3.tv_sec + (double)(t31.tv_usec - t3.tv_usec)/1000000
                                                                    , t32.tv_sec - t31.tv_sec + (double)(t32.tv_usec - t31.tv_usec)/1000000
                                                                    , t33.tv_sec - t32.tv_sec + (double)(t33.tv_usec - t32.tv_usec)/1000000
                                                                    , t4.tv_sec - t33.tv_sec + (double)(t4.tv_usec - t33.tv_usec)/1000000);

#endif

    }

    free (sizes);
    free (offsets);

    if (isAggregator (p))
    {
#if DEBUG
        struct timeval t00;
        gettimeofday (&t00, NULL);
#endif

        //sort_read_requests (p, &p->split_read_request_list1);
        merge_sort_read_requests (p, &p->split_read_request_list1);
#if DEBUG
        struct timeval t0, t01, t1, t2;
        gettimeofday (&t0, NULL);
#endif
        do_read (p);

#if DEBUG
        gettimeofday (&t01, NULL);
#endif
        shuffle_data (p);

//while (1);
#if DEBUG
        gettimeofday (&t1, NULL);
        printf ("[%3d] do_read time = (%f,%f,%f)\n", p->rank, t0.tv_sec - t00.tv_sec + (double)(t0.tv_usec - t00.tv_usec)/1000000, t01.tv_sec - t0.tv_sec + (double)(t01.tv_usec - t0.tv_usec)/1000000, t1.tv_sec - t01.tv_sec + (double)(t1.tv_usec - t01.tv_usec)/1000000);
#endif
        send_read_data (p);

#if DEBUG
        gettimeofday (&t2, NULL);
        printf ("[%3d] send_read_data time = %f\n", p->rank, t2.tv_sec - t1.tv_sec + (double)(t2.tv_usec - t1.tv_usec)/1000000);
#endif
    }
    else
    {
        get_read_data (p); 
    }

    free_proc_struct (p);
    free (gh);
    free_namelist ((gp->var_namelist), gp->vars_count);
    free_namelist ((gp->attr_namelist), gp->attrs_count);
    free (gp);

    return 0;
}

int adios_read_bp_staged1_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    // Find the attribute: full path is stored with a starting / 
    // Like in HDF5, we need to match names given with or without the starting /
    // startpos is 0 or 1 to indicate if the argument has starting / or not
    int attrid;
    int vstartpos = 0, fstartpos = 0; 
    struct BP_GROUP * gh = (struct BP_GROUP *)gp->gh;
    int offset;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_get_attr()");
        return adios_errno;
    }
    if (!attrname) {
        adios_error (err_invalid_attrname, "Null pointer passed as attribute name to adios_get_attr()!");
        return adios_errno;
    }

    offset = gh->attrs_offset;

    if (attrname[0] == '/') 
        vstartpos = 1;
    for (attrid=0; attrid<(gp->attrs_count);attrid++) {
        //if (gp->attr_namelist[attrid][0] == '/') 
        if (gh->fh->gattr_h->attr_namelist[attrid+offset][0] == '/') 
            fstartpos = 1;
        //if (!strcmp(gp->attr_namelist[attrid]+fstartpos, attrname+vstartpos))
        if (!strcmp(gh->fh->gattr_h->attr_namelist[attrid+offset]+fstartpos, attrname+vstartpos))
            break; 
    }
    if (attrid >= gp->attrs_count) {
        adios_error ( err_invalid_attrname, "Invalid attribute name %s", attrname);
        return adios_errno;
    }

    return adios_read_bp_staged1_get_attr_byid(gp, attrid, type, size, data);
}

int adios_read_bp_staged1_get_attr_byid (ADIOS_GROUP * gp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    int    i, offset, count;
    struct BP_GROUP * gh;
    BP_FILE * fh;
    struct adios_index_attribute_struct_v1 * attr_root;
    struct adios_index_var_struct_v1 * var_root;
    int    file_is_fortran;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_get_attr()");
        return adios_errno;
    }
    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return adios_errno;
    }
    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return adios_errno;
    }
    if (attrid < 0 || attrid >= gh->attrs_count) {
        adios_error (err_invalid_attrid, "Invalid attribute id %d (allowed 0..%d)", attrid, gh->attrs_count);
        return adios_errno;
    }

    attr_root = gh->attrs_root; /* need to traverse the attribute list of the group */
    for (i = 0; i < attrid && attr_root; i++)
        attr_root = attr_root->next;
    if (i != attrid) {
        adios_error (err_corrupted_attribute, "Attribute id=%d is valid but was not found in internal data structures!",attrid);
        return adios_errno; 
    }

    file_is_fortran = is_fortran_file (gp);

    if (attr_root->characteristics[0].value) {
        /* Attribute has its own value */
        *size = bp_get_type_size (attr_root->type, attr_root->characteristics[0].value);
        *type = attr_root->type;
        *data = (void *) malloc (*size);  
        if (*data)
            memcpy(*data, attr_root->characteristics[0].value, *size);
    }
    else if (attr_root->characteristics[0].var_id) {
        /* Attribute is a reference to a variable */
        /* FIXME: var ids are not unique in BP. If a group of variables are written several
           times under different path using adios_set_path(), the id of a variable is always
           the same (should be different). As a temporary fix, we look first for a matching
           id plus path between an attribute and a variable. If not found, then we look for
           a match on the ids only.*/
        var_root = gh->vars_root; 
        while (var_root) {
            if (var_root->id == attr_root->characteristics[0].var_id && 
                !strcmp(var_root->var_path, attr_root->attr_path))
                break;
            var_root = var_root->next;
        }
        if (!var_root) {
            var_root = gh->vars_root; 
            while (var_root) {
                if (var_root->id == attr_root->characteristics[0].var_id)
                    break;
                var_root = var_root->next;
            }
        }

        if (!var_root) {
            adios_error (err_invalid_attribute_reference, 
                   "Attribute %s/%s in group %s is a reference to variable ID %d, which is not found", 
                   attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                   attr_root->characteristics[0].var_id);
            return adios_errno;
        }

        /* default values in case of error */
        *data = NULL;
        *size = 0;
        *type = attr_root->type;

        /* FIXME: variable and attribute type may not match, then a conversion is needed. */
        /* Cases:
                1. attr has no type, var is byte array     ==> string
                2. attr has no type, var is not byte array ==> var type
                3. attr is string, var is byte array       ==> string
                4. attr type == var type                   ==> var type 
                5. attr type != var type                   ==> attr type and conversion needed 
        */
        /* Error check: attr cannot reference an array in general */
        if (var_root->characteristics[0].dims.count > 0) {
            if ( (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
                 (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
                 (var_root->characteristics[0].dims.count == 1)) {
                 ; // this conversions are allowed
            } else {
                adios_error (err_invalid_attribute_reference, 
                    "Attribute %s/%s in group %s, typeid=%d is a reference to an %d-dimensional array variable "
                    "%s/%s of type %s, which is not supported in ADIOS",
                    attr_root->attr_path, attr_root->attr_name, attr_root->group_name, attr_root->type,
                    var_root->characteristics[0].dims.count,
                    var_root->var_path, var_root->var_name, common_read_type_to_string(var_root->type));
                return adios_errno;
            }
        }

        if ( (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
             (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
             (var_root->characteristics[0].dims.count == 1) ) {
            /* 1D byte arrays are converted to string */
            /* 1. read in variable */
            char varname[512];
            char *tmpdata;
            uint64_t start, count;
            int status;
            start = 0; 
            count = var_root->characteristics[0].dims.dims[0];
            snprintf(varname, 512, "%s/%s", var_root->var_path, var_root->var_name);
            tmpdata = (char *) malloc (count+1);
            if (tmpdata == NULL) {
                adios_error (err_no_memory, 
                      "Cannot allocate memory of %lld bytes for reading in data for attribute %s/%s of group %s.",
                      count, attr_root->attr_path, attr_root->attr_name, attr_root->group_name);
                return adios_errno;
            }

            status = adios_read_bp_staged1_read_var (gp, varname, &start, &count, tmpdata);
            
            if (status < 0) {
                char *msg = strdup(adios_get_last_errmsg());
                adios_error ((enum ADIOS_ERRCODES) status, 
                      "Cannot read data of variable %s/%s for attribute %s/%s of group %s: %s",
                      var_root->var_path, var_root->var_name, 
                      attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                      msg);
                free(tmpdata);
                free(msg);
                return status;
            }

            *type = adios_string;
            if (file_is_fortran) {
                /* Fortran byte array to C string */
                *data = futils_fstr_to_cstr( tmpdata, (int)count); /* FIXME: supports only 2GB strings... */
                *size = strlen( (char *)data );
                free(tmpdata);
            } else {
                /* C array to C string */
                tmpdata[count] = '\0';
                *size = count+1;
                *data = tmpdata;
            }
        } else {
            /* other types are inherited */
            *type = var_root->type;
            *size = bp_get_type_size (var_root->type, var_root->characteristics[0].value);
            *data = (void *) malloc (*size);  
            if (*data)
                memcpy(*data, var_root->characteristics[0].value, *size);
        }
    }

    return 0;
}

static void _swap_order (int n, uint64_t * array)
{
    int i;
    uint64_t tmp;

    for (i = 0; i< n / 2; i++)
    {
        tmp = array[i];
        array[i] = array[n - 1 - i];
        array[n - 1 - i] = tmp;
    }
}

/* Reverse the order in an array in place.
   use swapping from Fortran/column-major order to ADIOS-read-api/C/row-major order and back
*/
static void swap_order (int n, uint64_t * array, int * tdim)
{
    int i;
    uint64_t tmp;

    _swap_order (n, array);

    if (* tdim > -1)
    {
        * tdim = (n - 1) - * tdim;
    }
}

/* Look up variable id based on variable name.
   Return index 0..gp->vars_count-1 if found, -1 otherwise
*/
static int adios_read_bp_staged1_find_var(ADIOS_GROUP *gp, const char *varname)
{
    // Find the variable: full path is stored with a starting / 
    // Like in HDF5, we need to match names given with or without the starting /
    // startpos is 0 or 1 to indicate if the argument has starting / or not
    int varid;
    int vstartpos = 0, fstartpos; 
    struct BP_GROUP * gh = (struct BP_GROUP *)gp->gh;
    int offset;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group");
        return -1;
    }
    if (!varname) {
        adios_error (err_invalid_varname, "Null pointer passed as variable name!");
        return -1;
    }

    /* Search in gp->fh->gvar_h->var_namelist instead of gp->var_namelist, because that
       one comes back from the user. One idiot (who writes this comment) sorted the
       list in place after gopen and before inq_var.
    */
    offset = gh->vars_offset;

    if (varname[0] == '/') 
        vstartpos = 1;
    for (varid=0; varid<(gp->vars_count);varid++) {
        fstartpos = 0;
        /* if (gp->var_namelist[varid][0] == '/') */
        fstartpos = 0;
        if (gh->fh->gvar_h->var_namelist[varid+offset][0] == '/')
            fstartpos = 1;
        /*if (!strcmp(gp->var_namelist[varid]+fstartpos, varname+vstartpos))*/
        if (!strcmp(gh->fh->gvar_h->var_namelist[varid+offset]+fstartpos, varname+vstartpos))
            break; 
    }
    if (varid >= gp->vars_count) {
        adios_error (err_invalid_varname, "Invalid variable name %s", varname);
        return -1;
    }
    return varid;
}

// NCSU - For custom memory allocation 
#define CALLOC(var, num, sz, comment)\
{\
    var = calloc (num, sz); \
    if (!var)    {\
        adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for ", comment, " in common_read_get_characteristics"); \
        return; \
    }\
}

#define MALLOC(var,sz,comment)\
{\
    var = malloc (sz); \
    if (!var)    {\
        adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for ", comment, " in common_read_get_characteristics"); \
        return; \
    }\
}\

// NCSU - Reading the statistics
/** Get value and statistics, allocate space for them too */
static void adios_read_bp_staged1_get_characteristics (struct adios_index_var_struct_v1 * var_root, ADIOS_VARINFO *vi)
{
    int i, j, c, count = 1;
    int size, sum_size, sum_type;

    vi->value = NULL;

    vi->gmin = vi->gmax = NULL;
    vi->gavg = NULL;
    vi->mins = vi->maxs = NULL;
    vi->avgs = NULL;
    vi->gstd_dev = NULL;
    vi->std_devs = NULL;
    vi->hist = NULL;

    // set value for scalars
    if (var_root->characteristics [0].value) {
        size = bp_get_type_size(var_root->type, var_root->characteristics [0].value);
        vi->value = (void *) malloc (size);
        
        if (vi->value)
           memcpy(vi->value, var_root->characteristics [0].value, size);
        else {
            adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for value in common_read_get_characteristics");
            return;
        }
    } else {
        vi->value = NULL;
    }
/*
    int npgs = var_root->characteristics_count, timestep, ntimes = -1;
    uint64_t gcnt = 0, * cnts;

    double *gsum = NULL, *gsum_square = NULL;
    double **sums = NULL, **sum_squares = NULL;

    int16_t map[32];
    memset (map, -1, sizeof(map));

    // Bitmap shows which statistical information has been calculated
    i = j = 0;
    while (var_root->characteristics[0].bitmap >> j)
    {
        if ((var_root->characteristics[0].bitmap >> j) & 1)
            map [j] = i ++;
        j ++;
     }

    if(vi->timedim >= 0)
    {    
        ntimes = vi->dims[0];

        if (map[adios_statistic_min] != -1)
        {
            MALLOC(vi->mins, ntimes * sizeof(void *), "minimum per timestep");
            for (i = 0; i < ntimes; i++)
                vi->mins[i] = 0;
        }

        if (map[adios_statistic_max] != -1)
        {
            MALLOC(vi->maxs, ntimes * sizeof(void *), "maximum per timestep");
            for (i = 0; i < ntimes; i++)
                vi->maxs[i] = 0;
        }

        if (map[adios_statistic_sum] != -1)
        {
            MALLOC(sums, ntimes * sizeof(double *), "summation per timestep");
            MALLOC(vi->avgs, ntimes * sizeof(double *), "average per timestep");

            for (i = 0; i < ntimes; i++)
                sums[i] = vi->avgs[i] = 0;

            CALLOC(cnts, ntimes, sizeof(uint64_t), "count of elements per timestep");
        }

        if (map[adios_statistic_sum_square] != -1)
        {
            MALLOC(sum_squares, ntimes * sizeof(double *), "summation per timestep");
            MALLOC(vi->std_devs, ntimes * sizeof(double *), "standard deviation per timestep");

            for (i = 0; i < ntimes; i++)
                vi->std_devs[i] = sum_squares[i] = 0;
        }
    }

    if (map[adios_statistic_hist] != -1 && (var_root->characteristics[0].stats[0][map[adios_statistic_hist]].data))
    {
        struct adios_index_characteristics_stat_struct * stats = var_root->characteristics[0].stats[0];
        struct adios_index_characteristics_hist_struct * hist = stats[map[adios_statistic_hist]].data;
        int num_breaks = hist->num_breaks;

        MALLOC(vi->hist, sizeof(struct ADIOS_HIST), "histogram");
        MALLOC(vi->hist->breaks, num_breaks * sizeof(double), "break points of histogram");    
        MALLOC(vi->hist->gfrequencies, (num_breaks + 1) * sizeof(uint32_t), "global frequencies of histogram");

        vi->hist->num_breaks = hist->num_breaks;
        vi->hist->min = hist->min;
        vi->hist->max = hist->max;

        memcpy(vi->hist->breaks, hist->breaks, num_breaks * sizeof(double));
        CALLOC(vi->hist->gfrequencies, (num_breaks + 1), bp_get_type_size(adios_unsigned_integer, ""), "global frequency");

        if (ntimes > 0)
        {
            MALLOC(vi->hist->frequenciess, (ntimes * sizeof(int32_t *)), "frequencies for timesteps");
            for(i = 0; i < ntimes; i++)
                CALLOC(vi->hist->frequenciess[i], (num_breaks + 1), bp_get_type_size(adios_unsigned_integer, ""), "frequency at timestep");
        }
    }

    size = bp_get_type_size (var_root->type, "");    
    sum_size = bp_get_type_size (adios_double, "");

    if (var_root->type == adios_complex || var_root->type == adios_double_complex)
    {
        int type;
        count = 3;
        timestep = 0;

        if (var_root->type == adios_complex)
            type = adios_double;
        else
            type = adios_long_double;

        // Only a double precision returned for all complex values
        size = bp_get_type_size (adios_double, "");    

           for (i=0; i<var_root->characteristics_count; i++)
        {
            if (ntimes > 0)
                timestep = var_root->characteristics[i].time_index - 1;

            if (!var_root->characteristics[i].stats)
                continue;

            struct adios_index_characteristics_stat_struct ** stats = var_root->characteristics[i].stats;

            if ((map[adios_statistic_finite] != -1) && (* ((uint8_t *) stats[0][map[adios_statistic_finite]].data) == 0))
                continue;

            if (map[adios_statistic_min] != -1 && stats[0][map[adios_statistic_min]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_min]].data);

                if(!vi->gmin) {
                    MALLOC (vi->gmin, count * size, "global minimum")
                    for (c = 0; c < count; c ++)
                           ((double * ) vi->gmin)[c] = data[c]; 

                } else {
                    for (c = 0; c < count; c ++)
                        if (data[c] < ((double *) vi->gmin)[c])
                               ((double * ) vi->gmin)[c] = data[c]; 
                }

                if (ntimes > 0) {
                    if(!vi->mins[timestep]) {
                        MALLOC (vi->mins[timestep], count * size, "minimum per timestep")
                        for (c = 0; c < count; c ++)
                               ((double **) vi->mins)[timestep][c] = data[c]; 

                    } else {
                        for (c = 0; c < count; c ++)
                            if (data[c] < ((double **) vi->mins)[timestep][c])
                                   ((double **) vi->mins)[timestep][c] = data[c]; 
                    }
                }
            }

            if (map[adios_statistic_max] != -1 && stats[0][map[adios_statistic_max]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_max]].data);

                if(!vi->gmax) {
                    MALLOC (vi->gmax, count * size, "global minimum")
                    for (c = 0; c < count; c ++)
                        ((double * ) vi->gmax)[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        if (data[c] > ((double *) vi->gmax)[c])
                            ((double * ) vi->gmax)[c] = data[c];
                }

                if (ntimes > 0) {
                    if(!vi->maxs[timestep]) {
                        MALLOC (vi->maxs[timestep], count * size, "minimum per timestep")
                        for (c = 0; c < count; c ++)
                            ((double **) vi->maxs)[timestep][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            if (data[c] > ((double **) vi->maxs)[timestep][c])
                                ((double **) vi->maxs)[timestep][c] = data[c];
                    }
                }
            }

            if (map[adios_statistic_sum] != -1 && stats[0][map[adios_statistic_sum]].data)
            {    
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_sum]].data);

                if(!gsum) {
                    MALLOC(gsum, count * sum_size, "global summation")
                    for (c = 0; c < count; c ++)
                           gsum[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        gsum[c] = gsum[c] + data[c];
                }

                if (ntimes > 0) {
                    if(!sums[timestep]) {
                        MALLOC(sums[timestep], count * sum_size, "summation per timestep")
                        for (c = 0; c < count; c ++)
                            sums[timestep][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            sums[timestep][c] = sums[timestep][c] + data[c];
                    }
                }
            }

            if (map[adios_statistic_sum_square] != -1 && stats[0][map[adios_statistic_sum_square]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_sum_square]].data);

                if(!gsum_square) {
                    MALLOC(gsum_square, count * sum_size, "global summation of squares")
                    for (c = 0; c < count; c ++)
                        gsum_square[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                           gsum_square[c] = gsum_square[c] + data[c]; 
                }

                if (ntimes > 0) {
                    if(!sum_squares[timestep]) {
                        MALLOC(sum_squares[timestep], count * sum_size, "summation of square per timestep")
                        for (c = 0; c < count; c ++)
                            sum_squares[timestep][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            sum_squares[timestep][c] = sum_squares[timestep][c] + data[c]; 
                    }
                }
            }

            if (map[adios_statistic_cnt] != -1 && stats[0][map[adios_statistic_cnt]].data)
            {
                if (ntimes > 0)
                    cnts[timestep] += * ((uint32_t *) stats[0][map[adios_statistic_cnt]].data);
                gcnt += * (uint32_t *) stats[0][map[adios_statistic_cnt]].data;
            }
        }

        if(ntimes > 0 && vi->gmin && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            // min, max, summation exists only for arrays
            // Calculate average / timestep

            for(timestep = 0; timestep < ntimes; timestep ++) {
                MALLOC(vi->avgs[timestep], count * sum_size, "average per timestep")
                for (c = 0; c < count; c ++)
                    vi->avgs[timestep][c] = sums[timestep][c] / cnts[timestep];

                MALLOC(vi->std_devs[timestep], count * sum_size, "standard deviation per timestep")
                for (c = 0; c < count; c ++)
                    vi->std_devs[timestep][c] = sqrt((sum_squares[timestep][c] / cnts[timestep]) - (vi->avgs[timestep][c] * vi->avgs[timestep][c]));

                free (sums[timestep]);
                free (sum_squares[timestep]);
            }
        }

        // Calculate global average
        if(vi->gmin && gsum && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            MALLOC(vi->gavg, count * sum_size, "global average")

            if(gcnt > 0)
                for (c = 0; c < count; c ++)
                    vi->gavg[c] = gsum[c] / gcnt;
            else
                for (c = 0; c < count; c ++)
                    vi->gavg[c] = gsum[c];

            MALLOC(vi->gstd_dev, count * sum_size, "global average")
            if(vi->gavg && gcnt > 0)
                for (c = 0; c < count; c ++)
                    vi->gstd_dev[c] = sqrt(gsum_square[c] / gcnt - (vi->gavg[c] * vi->gavg[c]));
            else
                for (c = 0; c < count; c ++)
                    vi->gstd_dev[c] = 0;
        }
    }
    else
    {
        timestep = 0;
           for (i=0; i<var_root->characteristics_count; i++)
        {
            if (ntimes > 0)
                timestep = var_root->characteristics[i].time_index - 1;
                //timestep = i / (npgs / ntimes);

            if (!var_root->characteristics[i].stats)
                continue;

            struct adios_index_characteristics_stat_struct * stats = var_root->characteristics[i].stats[0];
            struct adios_index_characteristics_hist_struct * hist = stats[map[adios_statistic_hist]].data;

            if (map[adios_statistic_finite] != -1 && (* ((uint8_t *) stats[map[adios_statistic_finite]].data) == 0))
                continue;

            if (map[adios_statistic_min] != -1 && stats[map[adios_statistic_min]].data)
            {
                if(!vi->gmin) {
                    MALLOC (vi->gmin, size, "global minimum")
                       memcpy(vi->gmin, stats[map[adios_statistic_min]].data, size);

                } else if (adios_lt(var_root->type, stats[map[adios_statistic_min]].data, vi->gmin)){
                       memcpy(vi->gmin, stats[map[adios_statistic_min]].data, size);
                }

                if (ntimes > 0) {
                    if(!vi->mins[timestep]) {
                        MALLOC (vi->mins[timestep], size, "minimum per timestep")
                        memcpy(vi->mins[timestep], stats[map[adios_statistic_min]].data, size);

                    } else if (adios_lt(var_root->type, stats[map[adios_statistic_min]].data, vi->mins[timestep])) {
                        memcpy(vi->mins[timestep], stats[map[adios_statistic_min]].data, size);
                    }
                }
            }

            if (map[adios_statistic_max] != -1 && stats[map[adios_statistic_max]].data)
            {
                if(!vi->gmax) {
                    MALLOC (vi->gmax, size, "global maximum")
                    memcpy(vi->gmax, stats[map[adios_statistic_max]].data, size);

                } else if (adios_lt(var_root->type, vi->gmax, stats[map[adios_statistic_max]].data))
                       memcpy(vi->gmax, stats[map[adios_statistic_max]].data, size);
                
                if (ntimes > 0) {
                    if(!vi->maxs[timestep]) {
                        MALLOC (vi->maxs[timestep], size, "maximum per timestep")
                        memcpy(vi->maxs[timestep], stats[map[adios_statistic_max]].data, size);

                    } else if (adios_lt(var_root->type, vi->maxs[timestep], stats[map[adios_statistic_max]].data)) {
                        memcpy(vi->maxs[timestep], stats[map[adios_statistic_max]].data, size);    
                    }
                }
            }
    
            if (map[adios_statistic_sum] != -1 && stats[map[adios_statistic_sum]].data)
            {    
                if(!gsum) {
                    MALLOC(gsum, sum_size, "global summation")
                       memcpy(gsum, stats[map[adios_statistic_sum]].data, sum_size);

                } else {
                    *gsum = *gsum + * ((double *) stats[map[adios_statistic_sum]].data);
                }

                if (ntimes > 0) {
                    if(!sums[timestep]) {
                        MALLOC(sums[timestep], sum_size, "summation per timestep")
                        memcpy(sums[timestep], stats[map[adios_statistic_sum]].data, sum_size);

                    } else {
                        *sums[timestep] = *sums[timestep] + * ((double *) stats[map[adios_statistic_sum]].data);
                    }
                }
            }

            if (map[adios_statistic_sum_square] != -1 && stats[map[adios_statistic_sum_square]].data)
            {
                if(!gsum_square) {
                    MALLOC(gsum_square, sum_size, "global summation of squares")
                    memcpy(gsum_square, stats[map[adios_statistic_sum_square]].data, sum_size);

                } else {
                       *gsum_square = *gsum_square + * ((double *) stats[map[adios_statistic_sum_square]].data);
                }

                if (ntimes > 0) {
                    if(!sum_squares[timestep]) {
                        MALLOC(sum_squares[timestep], sum_size, "summation of square per timestep")
                        memcpy(sum_squares[timestep], stats[map[adios_statistic_sum_square]].data, sum_size);

                    } else {
                        *sum_squares[timestep] = *sum_squares[timestep] + * ((double *) stats[map[adios_statistic_sum_square]].data);
                    }
                }
            }

            if(map[adios_statistic_hist] != -1 && stats[map[adios_statistic_hist]].data)
            {
                for(j = 0; j <= vi->hist->num_breaks; j++)
                {    
                    uint32_t freq = hist->frequencies[j];
                    vi->hist->gfrequencies[j] += freq;
                    if (ntimes > 0)
                        vi->hist->frequenciess[timestep][j] += freq;
                }
            }

            if (map[adios_statistic_cnt] != -1 && stats[map[adios_statistic_cnt]].data)
            {
                if (ntimes > 0)
                    cnts[timestep] += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
                gcnt += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
            }
        }

        if(ntimes > 0 && vi->gmin && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            // min, max, summation exists only for arrays
            // Calculate average / timestep

            for(timestep = 0; timestep < ntimes; timestep ++) {
                MALLOC(vi->avgs[timestep], sum_size, "average per timestep")
                *(vi->avgs[timestep]) = *(sums[timestep]) / cnts[timestep];

                MALLOC(vi->std_devs[timestep], sum_size, "standard deviation per timestep")
                *(vi->std_devs[timestep]) = sqrt(*(sum_squares[timestep]) / cnts[timestep] - ((*(vi->avgs[timestep]) * (*(vi->avgs[timestep])))));

                free (sums[timestep]);
                free (sum_squares[timestep]);
            }
        }

        // Calculate global average
        if(vi->gmin && gsum && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            MALLOC(vi->gavg, sum_size, "global average")
            if(gcnt > 0)
                *vi->gavg = *gsum / gcnt;
            else
                vi->gavg = gsum;

            MALLOC(vi->gstd_dev, sum_size, "global average")
            if(vi->gavg && gcnt > 0)
                *vi->gstd_dev = sqrt(*gsum_square / gcnt - ((*(vi->gavg)) * (*(vi->gavg))));
            else
                *vi->gstd_dev = 0;
        }
    }

    if (!vi->value && vi->gmin) {
        vi->value = vi->gmin; // arrays have no value but we assign here the minimum
    }

    if(!vi->gmin) {
        vi->gmin = vi->value; // scalars have value but not min
    }
    if(!vi->gmax) {
        vi->gmax = vi->value; // scalars have value but not max
    }

    if (sums && gsum) {
        free (sums);
        free (gsum);
    }

    if (sum_squares && gsum_square) {
        free (sum_squares);
        free (gsum_square);
    } 
*/
}

/* get local and global dimensions and offsets from a variable characteristics 
   return: 1 = it is a global array, 0 = local array
*/
static int adios_read_bp_staged1_get_dimensioncharacteristics (
                           struct adios_index_characteristic_struct_v1 * c
                          ,uint64_t * ldims
                          ,uint64_t * gdims
                          ,uint64_t * offsets
                          )
{
    int is_global = 0;
    int ndim = c->dims.count;
    int i;

    for (i = 0; i < ndim; i++)
    {
        ldims[i] = c->dims.dims[i * 3];
        gdims[i] = c->dims.dims[i * 3 + 1];
        offsets[i] = c->dims.dims[i * 3 + 2];

        if (gdims[i])
        {
            is_global = 1;
        }
    }

    return is_global;
}

static void adios_read_bp_staged1_get_dimensions (struct adios_index_var_struct_v1 * v
                                                 ,int ntsteps, int file_is_fortran
                                                 ,int * ndim, uint64_t ** dims
                                                 ,int * tdim
                                                 )
{
    int i, k;
    int is_global; // global array or just an array written by one process?
    uint64_t ldims[32], gdims[32], offsets[32];

    /* Get dimension information */    
    * ndim = v->characteristics [0].dims.count; 
/*
int rank;
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
fprintf (stderr, "[%3d] ndim 2 = %d\n", rank, * ndim);
*/
    * dims = NULL;
    * tdim = -1;

    if (* ndim == 0)
    {
        return;
    }

    * dims = (uint64_t *) malloc (sizeof (uint64_t) * (* ndim));
    memset (* dims, 0, sizeof (uint64_t) * (* ndim));

    is_global = adios_read_bp_staged1_get_dimensioncharacteristics (&(v->characteristics[0])
                                                                   ,ldims
                                                                   ,gdims
                                                                   ,offsets
                                                                   );
    if (!is_global)
    {
        for (i = 0; i < * ndim; i++)
        { 
            if (ldims[i] == 1 && v->characteristics_count > 1)
            {
                * tdim = i;
                (* dims)[i] = ntsteps;
            }
            else
            {
                (* dims)[i] = ldims[i];
            }

            gdims[i] = ldims[i];
        }
    }         
    else
    {
        if (gdims[* ndim - 1] == 0)
        {
            * tdim = (!file_is_fortran) ? 0 : * ndim - 1;
            (* dims)[* tdim] = ntsteps;

            if (* ndim > 1 && ldims[* tdim] != 1)
            {
                fprintf (stderr, "ndim = %d\n", * ndim);
                fprintf (stderr, "name = %s\n", v->var_name);
                fprintf (stderr,"ADIOS Error: this is a BP file with C ordering but we didn't find "
                        "an array to have time dimension in the first dimension. l:g:o = (");
                for (i = 0; i < * ndim; i++)
                {
                    fprintf (stderr,"%llu:%llu:%llu%s", ldims[i], gdims[i], offsets[i]
                                                      , (i < * ndim - 1 ? ", " : "") );
                }

                fprintf (stderr, ")\n");
            }

            if (!file_is_fortran)
            {
                for (i = 1; i < * ndim; i++)
                {
                    (* dims)[i] = gdims[i - 1];
                }
            }
            else
            {
                for (i = 0; i < * ndim - 1; i++)
                {
                    (* dims)[i] = gdims[i];
                }
            }
        }
        else
        {
            for (i = 0; i < * ndim; i++)
            {
                (* dims)[i] = gdims[i];
            }
        }
    }
}

ADIOS_VARINFO * adios_read_bp_staged1_inq_var (ADIOS_GROUP *gp, const char * varname) 
{
    int varid = adios_read_bp_staged1_find_var(gp, varname);
    if (varid < 0)
        return NULL;
    return adios_read_bp_staged1_inq_var_byid(gp, varid);
}

static ADIOS_VARINFO * _inq_var_byid (BP_FILE * fh, int varid)
{
    ADIOS_VARINFO * vi;
    int file_is_fortran;
    struct adios_index_var_struct_v1 * var_root = fh->vars_root;
    int i,k;

    assert (fh);

    vi = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    assert (vi);

    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    
    for (i = 1; i < varid && var_root; i++)
    {
        var_root = var_root->next;
    }

    if (i!=varid) {
        adios_error (err_corrupted_variable, "Variable id=%d is valid but was not found in internal data structures!",varid);
        return NULL; 
    }

    vi->varid = varid;
    vi->type = var_root->type;
    vi->value = 0;
    vi->gmin = 0;
    vi->gmax = 0;

    if (!var_root->characteristics_count)
    {
        free(vi);
        return NULL;
    }
    vi->characteristics_count = var_root->characteristics_count;

    /* Get value or min/max */
    adios_read_bp_staged1_get_dimensions (var_root, fh->tidx_stop - fh->tidx_start + 1, file_is_fortran, 
                          &(vi->ndim), &(vi->dims), &(vi->timedim));
#if 0
    if (file_is_fortran != futils_is_called_from_fortran()) {
        /* If this is a Fortran written file and this is called from C code,  
           or this is a C written file and this is called from Fortran code ==>
           We need to reverse the order of the dimensions */
        swap_order(vi->ndim, vi->dims, &(vi->timedim));
        /*printf("File was written from %s and read now from %s, so we swap order of dimensions\n",
                (file_is_fortran ? "Fortran" : "C"), (futils_is_called_from_fortran() ? "Fortran" : "C"));*/
    }
#endif
    adios_read_bp_staged1_get_characteristics (var_root, vi);

    return vi;
}

ADIOS_VARINFO * adios_read_bp_staged1_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    struct BP_GROUP      * gh;
    BP_FILE       * fh;
    ADIOS_VARINFO * vi;
    int file_is_fortran;
    struct adios_index_var_struct_v1 * var_root;
    int i,k;

    adios_errno = 0;
    assert (gp);

    gh = (struct BP_GROUP *) gp->gh;
    assert (gh);

    fh = gh->fh;
    assert (fh);

    if (varid < 0 || varid >= gh->vars_count) {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return NULL;
    }
    vi = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    assert (vi);

    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    
    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i=0; i<varid && var_root; i++) {
        var_root = var_root->next;
    }

    if (i!=varid) {
        adios_error (err_corrupted_variable, "Variable id=%d is valid but was not found in internal data structures!",varid);
        return NULL; 
    }

    vi->varid = varid;
    vi->type = var_root->type;
    vi->value = 0;
    vi->gmin = 0;
    vi->gmax = 0;

    if (!var_root->characteristics_count) {
        adios_error (err_corrupted_variable, "Variable %s does not have information on dimensions", 
              gp->var_namelist[varid]);
        free(vi);
        return NULL;
    }
    vi->characteristics_count = var_root->characteristics_count;

    /* Get value or min/max */

    adios_read_bp_staged1_get_dimensions (var_root, fh->tidx_stop - fh->tidx_start + 1, file_is_fortran, 
                          &(vi->ndim), &(vi->dims), &(vi->timedim));

    if (file_is_fortran != futils_is_called_from_fortran()) {
        /* If this is a Fortran written file and this is called from C code,  
           or this is a C written file and this is called from Fortran code ==>
           We need to reverse the order of the dimensions */
        swap_order(vi->ndim, vi->dims, &(vi->timedim));
        /*printf("File was written from %s and read now from %s, so we swap order of dimensions\n",
                (file_is_fortran ? "Fortran" : "C"), (futils_is_called_from_fortran() ? "Fortran" : "C"));*/
    }

    adios_read_bp_staged1_get_characteristics (var_root, vi);

    return vi;
}

void adios_read_bp_staged1_free_varinfo (ADIOS_VARINFO *vp)
{
    if (vp) {
        if (vp->dims)   free(vp->dims);
        if (vp->value)  free(vp->value);
        if (vp->gmin && vp->gmin != vp->value)   free(vp->gmin);
        if (vp->gmax && vp->gmax != vp->value)   free(vp->gmax);
        //if (vp->mins)   free(vp->mins);
        //if (vp->maxs)   free(vp->maxs);
        free(vp);
    }
}

#define MPI_FILE_READ_OPS                           \
        bp_realloc_aligned(fh->b, slice_size);      \
        fh->b->offset = 0;                          \
                                                    \
        MPI_File_seek (fh->mpi_fh                   \
                      ,(MPI_Offset)slice_offset     \
                      ,MPI_SEEK_SET                 \
                      );                            \
                                                    \
        MPI_File_read (fh->mpi_fh                   \
                      ,fh->b->buff                  \
                      ,slice_size                   \
                      ,MPI_BYTE                     \
                      ,&status                      \
                      );                            \
        fh->b->offset = 0;                          \

//We also need to be able to read old .bp which doesn't have 'payload_offset'
#define MPI_FILE_READ_OPS1                                                                  \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) var_root->characteristics[start_idx + idx].offset       \
                      ,MPI_SEEK_SET);                                                       \
        MPI_File_read (fh->mpi_fh, fh->b->buff, 8, MPI_BYTE, &status);                      \
        tmpcount= *((uint64_t*)fh->b->buff);                                                \
                                                                                            \
        bp_realloc_aligned(fh->b, tmpcount + 8);                                            \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) (var_root->characteristics[start_idx + idx].offset)     \
                      ,MPI_SEEK_SET);                                                       \
        MPI_File_read (fh->mpi_fh, fh->b->buff, tmpcount + 8, MPI_BYTE, &status);           \
        fh->b->offset = 0;                                                                  \
        adios_parse_var_data_header_v1 (fh->b, &var_header);                                \


// To read subfiles
#define MPI_FILE_READ_OPS2                                                                  \
        bp_realloc_aligned(fh->b, slice_size);                                              \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File * sfh;                                                                     \
        sfh = get_BP_file_handle (fh->sfh                                                   \
                                 ,var_root->characteristics[start_idx + idx].file_index     \
                                 );                                                         \
        if (!sfh)                                                                           \
        {                                                                                   \
            int err;                                                                        \
            char * ch, * name_no_path, * name;                                              \
            struct BP_file_handle * new_h =                                                 \
                  (struct BP_file_handle *) malloc (sizeof (struct BP_file_handle));        \
            new_h->file_index = var_root->characteristics[start_idx + idx].file_index;      \
            new_h->next = 0;                                                                \
            if (ch = strrchr (fh->fname, '/'))                                              \
            {                                                                               \
                name_no_path = malloc (strlen (ch + 1) + 1);                                \
                strcpy (name_no_path, ch + 1);                                              \
            }                                                                               \
            else                                                                            \
            {                                                                               \
                name_no_path = malloc (strlen (fh->fname) + 1);                             \
                strcpy (name_no_path, fh->fname);                                           \
            }                                                                               \
                                                                                            \
            name = malloc (strlen (fh->fname) + 5 + strlen (name_no_path) + 1 + 10 + 1);    \
            sprintf (name, "%s.dir/%s.%d", fh->fname, name_no_path, new_h->file_index);     \
                                                                                            \
            err = MPI_File_open (MPI_COMM_SELF                                              \
                                ,name                                                       \
                                ,MPI_MODE_RDONLY                                            \
                                ,(MPI_Info)MPI_INFO_NULL                                    \
                                ,&new_h->fh                                                 \
                                );                                                          \
           if (err != MPI_SUCCESS)                                                          \
           {                                                                                \
               fprintf (stderr, "can not open file %S\n", name);                            \
               return -1;                                                                   \
           }                                                                                \
                                                                                            \
           add_BP_file_handle (&fh->sfh                                                     \
                              ,new_h                                                        \
                              );                                                            \
           sfh = &new_h->fh;                                                                \
                                                                                            \
           free (name_no_path);                                                             \
           free (name);                                                                     \
        }                                                                                   \
                                                                                            \
        MPI_File_seek (*sfh                                                                 \
                      ,(MPI_Offset)slice_offset                                             \
                      ,MPI_SEEK_SET                                                         \
                      );                                                                    \
        MPI_File_read (*sfh                                                                 \
                      ,fh->b->buff                                                          \
                      ,slice_size                                                           \
                      ,MPI_BYTE                                                             \
                      ,&status                                                              \
                      );                                                                    \
        fh->b->offset = 0;                                                                  \

// Search for the start var index.
static int get_var_start_index (struct adios_index_var_struct_v1 * v, int t)
{
    int i = 0;

    while (i < v->characteristics_count) {
        if (v->characteristics[i].time_index == t) {
            return i;
        }

        i++;
    }

    return -1;
}

// Search for the stop var index
static int get_var_stop_index (struct adios_index_var_struct_v1 * v, int t)
{
    int i = v->characteristics_count - 1;

    while (i > -1) {
        if (v->characteristics[i].time_index == t) {
            return i;
        }

        i--;
    }

    return -1;
}

/****************************************************
  Find the var associated with the given variable id 
*****************************************************/
static struct adios_index_var_struct_v1 * adios_find_var_byid (ADIOS_GROUP * gp, int varid)
{
    struct BP_GROUP * gh;
    struct adios_index_var_struct_v1 * var_root;
    int i;

    gh = (struct BP_GROUP *) gp->gh;
    var_root = gh->vars_root;

    for (i = 0; i < varid && var_root; i++)
    {
        var_root = var_root->next;
    }

    if (i != varid)
    {
        adios_error (err_corrupted_variable,
               "Variable id=%d is valid but was not found in internal data structures!",
               varid);
        return NULL;
    }

    return var_root;
}

/* Check whether it has time */
static int isTimeless (int tdim)
{
    return tdim <= -1;
}

/* This routine is used to Populate read_info data structure  */
static void getReadInfo (ADIOS_GROUP * gp
                        ,int varid
                        ,uint64_t * start
                        ,uint64_t * count
                        ,read_info * ri
                        )
{
    struct BP_GROUP * gh;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * v;
    int i, j, k, idx, t;
    int start_idx, stop_idx, f_idx;
    int has_subfile, file_is_fortran;
    uint64_t size, * dims;
    int file_tdim = -1, read_arg_tdim, is_global = 0, flag;

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;
    v = adios_find_var_byid (gp, varid);

    file_is_fortran = is_fortran_file (gp);
    has_subfile = has_subfiles (gp);

    /* Get dimensions and flip if caller != writer language */
    adios_read_bp_staged1_get_dimensions (v
                                        ,fh->tidx_stop - fh->tidx_start + 1
                                        ,file_is_fortran
                                        ,&ri->ndim
                                        ,&ri->dims
                                        ,&ri->file_tdim
                                        );

    if (file_is_fortran)
    {
        swap_order (ri->ndim, ri->dims, &ri->file_tdim);
    }

    if (isTimeless (ri->file_tdim))
    {
        ri->start_time = fh->tidx_start;
        ri->stop_time = fh->tidx_stop;
        ri->ndim_notime = ri->ndim;

        ri->start_notime = malloc (ri->ndim * 8);
        ri->count_notime = malloc (ri->ndim * 8);
        assert (ri->start_notime);
        assert (ri->count_notime);

        memcpy (ri->start_notime, start, ri->ndim * 8);
        memcpy (ri->count_notime, count, ri->ndim * 8);
    }
    else
    {
        read_arg_tdim = futils_is_called_from_fortran () ? ri->ndim - 1 : 0;

        ri->start_time = start[read_arg_tdim] + fh->tidx_start;
        ri->stop_time = ri->start_time + count[read_arg_tdim] - 1;
        ri->ndim_notime = ri->ndim - 1;

        ri->start_notime = malloc (ri->ndim_notime * 8);
        ri->count_notime = malloc (ri->ndim_notime * 8);
        assert (ri->start_notime);
        assert (ri->count_notime);

        memcpy (ri->start_notime
               ,futils_is_called_from_fortran () ? start : start + 1
               , ri->ndim_notime * 8
               );
        memcpy (ri->count_notime
               ,futils_is_called_from_fortran () ? count : count + 1
               ,ri->ndim_notime * 8
               );
    }

    if (futils_is_called_from_fortran ())
    {
        _swap_order (ri->ndim_notime, ri->count_notime);
        _swap_order (ri->ndim_notime, ri->start_notime);
    }


}

static void freeReadInfo (read_info * ri)
{
    if (ri->start_notime)
    {
        free (ri->start_notime);
    }

    if (ri->count_notime)
    {
        free (ri->count_notime);
    }

    if (ri->dims)
    {
        free (ri->dims);
    }
}

/**************************************************
* Get the subfile index and associate data offset *
***************************************************/
static void getDataAddress (ADIOS_GROUP * gp, int varid
                           ,const uint64_t * start
                           ,const uint64_t * count
                           ,int * file_idx
                           ,uint64_t * offset
                           ,uint64_t * payload_size
                           )
{
    struct BP_GROUP * gh;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * v;
    int i, j, k, idx, t;
    int start_time, stop_time, start_idx, stop_idx, f_idx;
    int ndim, ndim_notime, has_subfile, file_is_fortran;
    uint64_t size, * dims = 0;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t count_notime[32], start_notime[32];
    int file_tdim = -1, read_arg_tdim, is_global = 0, flag;

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;

    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    has_subfile = fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE;

    v = adios_find_var_byid (gp, varid);

    /* Get dimensions and flip if caller != writer language */
    adios_read_bp_staged1_get_dimensions (v
                                         ,fh->tidx_stop - fh->tidx_start + 1
                                         ,file_is_fortran
                                         ,&ndim
                                         ,&dims
                                         ,&file_tdim
                                         );

    if (file_is_fortran)
    {
        swap_order (ndim, dims, &file_tdim);
    }

    if (isTimeless (file_tdim))
    {
        start_time = fh->tidx_start;
        stop_time = fh->tidx_stop;
        ndim_notime = ndim;

        memcpy (start_notime, start, ndim * 8);
        memcpy (count_notime, count, ndim * 8);
    }
    else
    {
        read_arg_tdim = futils_is_called_from_fortran () ? ndim - 1 : 0;

        start_time = start[read_arg_tdim] + fh->tidx_start;
        stop_time = start_time + count[read_arg_tdim] - 1;
        ndim_notime = ndim - 1;

        memcpy (start_notime
               ,futils_is_called_from_fortran () ? start : start + 1
               ,ndim_notime * 8
               );
        memcpy (count_notime
               ,futils_is_called_from_fortran () ? count : count + 1
               ,ndim_notime * 8
               );
    }

    for (t = start_time; t <= stop_time; t++)
    {
        start_idx = get_var_start_index (v, t);
        stop_idx = get_var_stop_index (v, t);

        if (start_idx < 0 || stop_idx < 0)
        {
            adios_error (err_no_data_at_timestep,"Variable (id=%d) has no data at %d time step",
                varid, t);
            continue;
        }

        if (ndim_notime == 0)
        {
            /* THIS IS A SCALAR VARIABLE */
            idx = 0;

            if (isTimeless (file_tdim) )
                break;
            else
                continue;
        }

         /* READ AN ARRAY VARIABLE */
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));

        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
        {
            idx_table[idx] = 1;
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = adios_read_bp_staged1_get_dimensioncharacteristics(&(v->characteristics[start_idx + idx])
                                                                          ,ldims
                                                                          ,gdims
                                                                          ,offsets
                                                                          );
            if (!is_global)
            {
                memcpy (gdims, ldims, ndim * 8);
            }

            if (file_is_fortran)
            {
                _swap_order (ndim, gdims);
                _swap_order (ndim, ldims);
                _swap_order (ndim, offsets);
            }

            if (!isTimeless (file_tdim))
            {
                for (i = file_tdim; i < ndim - 1; i++)
                {
                    ldims[i] = ldims[i + 1];
                    if (file_is_fortran)
                    {
                        gdims[i] = gdims[i + 1];
                        offsets[i] = offsets[i + 1];
                    }
                }
            }

            /*
            printf("ldims   = "); for (j = 0; j<ndim; j++) printf("%d ",ldims[j]); printf("\n");
            printf("gdims   = "); for (j = 0; j<ndim; j++) printf("%d ",gdims[j]); printf("\n");
            printf("offsets = "); for (j = 0; j<ndim; j++) printf("%d ",offsets[j]); printf("\n");
            printf("count_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
            printf("start_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
            */
            for (j = 0; j < ndim_notime; j++)
            {
                if ( (count_notime[j] > gdims[j])
                  || (start_notime[j] > gdims[j])
                  || (start_notime[j] + count_notime[j] > gdims[j]))
                {
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j + 1, count_notime[j], start_notime[j], gdims[j] - 1);
                    return;
                }

                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start_notime[j]
                        && offsets[j] < start_notime[j] + count_notime[j])
                    || (offsets[j] < start_notime[j]
                        && offsets[j] + ldims[j] > start_notime[j] + count_notime[j])
                    || (offsets[j] + ldims[j] > start_notime[j]
                        && offsets[j] + ldims[j] <= start_notime[j] + count_notime[j]);

                idx_table[idx] = idx_table[idx] && flag;
            }

            //FIXME
            if (idx_table[idx])
            {
                free (idx_table);
                if (dims)
                {
                    free (dims);
                }

                * file_idx = v->characteristics[start_idx + idx].file_index;
                * offset = v->characteristics[start_idx + idx].payload_offset;
                * payload_size = bp_get_type_size (v->type, v->characteristics[start_idx + idx].value);
                for (j = 0; j < ndim_notime; j++)
                {
                    * payload_size *= ldims[j];
                }
                return;
            }
        }

        free (idx_table);

        if (isTimeless (file_tdim))
            break;
    } // end for (timestep ... loop over timesteps

    if (dims)
    {
        free (dims);
    }
}

int64_t adios_read_bp_staged1_read_var (ADIOS_GROUP * gp
                                       ,const char * varname
                                       ,const uint64_t * start
                                       ,const uint64_t * count
                                       ,void * data
                                       )
{
    struct BP_GROUP * gh;
    BP_FILE * fh;
    int i, varid, has_subfile, rank, nproc;
    uint64_t ds, payload_size;
    struct proc_struct * p;
    candidate_reader * r;

    adios_errno = 0;

    assert (gp);

    gh = (struct BP_GROUP *) gp->gh;
    assert (gh);

    fh = gh->fh;
    assert (fh);

    varid = adios_read_bp_staged1_find_var (gp, varname);
    if (varid < 0 || varid >= gh->vars_count)
    {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return adios_errno;
    }

    p = (struct proc_struct *) fh->priv;
    assert (p);

    read_args * ra = (read_args *) malloc (sizeof (read_args));
    assert (ra);

    ADIOS_VARINFO * vi = adios_read_bp_staged1_inq_var_byid (gp, varid);
    ra->varid = varid;
    ra->ndims = vi->ndim;

    ra->start = (uint64_t *) malloc (ra->ndims * 8);
    memcpy (ra->start, start, ra->ndims * 8);

    ra->count = (uint64_t *) malloc (ra->ndims * 8);
    memcpy (ra->count, count, ra->ndims * 8);

    ra->data = data;

    ra->size = bp_get_type_size (vi->type, 0);

    for (i = 0; i < vi->ndim; i++)
    {
        ra->size *= ra->count[i];
    }

    ra->parent = 0;
    ra->file_idx = -1;
    ra->offset = 0;
//    getDataAddress (gp, varid, start, count, &ra->file_idx, &ra->offset, &payload_size);

    adios_read_bp_staged1_free_varinfo (vi);

    r = (candidate_reader *) malloc (sizeof (candidate_reader));
    assert (r);

    r->rank = p->rank;
    r->ra = ra;
    r->next = NULL;

    list_insert_reader (&p->local_read_request_list, r);

    return 0;
}

/***********************************************
 * This routine is to read in data in a 'local *
 * array fashion (as opposed to global array)  *
 *     Q. Liu, 11/2010                         *
 ***********************************************/
int64_t adios_read_bp_staged1_read_local_var (ADIOS_GROUP * gp, const char * varname,
                                      int vidx, const uint64_t * start,
                                      const uint64_t * count, void * data)
{
    struct BP_GROUP      * gh;
    BP_FILE       * fh;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    int    i,j,k, t, varid, start_idx, idx;
    int    ndim, ndim_notime, has_subfile, file_is_fortran;
    uint64_t size, * dims;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=0, items_read;
    uint64_t count_notime[32], start_notime[32];
    int timedim = -1, temp_timedim, is_global = 0, size_of_type;
    uint64_t slice_offset, slice_size, tmpcount = 0;
    uint64_t datatimeoffset = 0; // offset in data to write a given timestep
    MPI_Status status;

    adios_errno = 0;
    if (!gp)
    {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return adios_errno;
    }

    gh = (struct BP_GROUP *) gp->gh;
    if (!gh)
    {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return adios_errno;
    }

    fh = gh->fh;
    if (!fh)
    {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return adios_errno;
    }

    varid = adios_read_bp_staged1_find_var(gp, varname);
    if (varid < 0 || varid >= gh->vars_count)
    {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return adios_errno;
    }

    /* Check if file is written out by Fortran or C */
    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);

    /* Check whether we need to handle subfiles */
    has_subfile = fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE;

    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i = 0; i< varid && var_root; i++)
    {
        var_root = var_root->next;
    }

    if (i != varid)
    {
        adios_error (err_corrupted_variable, 
                     "Variable id=%d is valid but was not found in internal data structures!",
                     varid);
        return adios_errno; 
    }

    if (vidx < 0 || vidx >= var_root->characteristics_count)
    {
        adios_error (err_out_of_bound, "idx=%d is out of bound", vidx);
    }

    ndim = var_root->characteristics [vidx].dims.count;

    /* count_notime/start_notime are working copies of count/start */
    for (i = 0; i < ndim; i++)
    {
        count_notime[i] = count[i];
        start_notime[i] = start[i];
    }

    ndim_notime = ndim;

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if (futils_is_called_from_fortran())
    {
        timedim = -1;
        swap_order(ndim_notime, count_notime, &timedim);
        swap_order(ndim_notime, start_notime, &timedim);
    }
    
    /* items_read = how many data elements are we going to read in total */
    items_read = 1;
    for (i = 0; i < ndim_notime; i++)
        items_read *= count_notime[i];

    size_of_type = bp_get_type_size (var_root->type, var_root->characteristics [vidx].value);

    /* READ A SCALAR VARIABLE */
    if (ndim_notime == 0)
    {
        slice_size = size_of_type;
        start_idx = 0; // OPS macros below need it
        idx = vidx; // OPS macros below need it

        if (var_root->type == adios_string)
        {
            // strings are stored without \0 in file
            // size_of_type here includes \0 so decrease by one
            size_of_type--;
        }

        /* Old BP files don't have payload_offset characteristic */
        if (var_root->characteristics[vidx].payload_offset > 0)
        {
            slice_offset = var_root->characteristics[vidx].payload_offset;

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
        }
        else
        {
            slice_offset = 0;
            MPI_FILE_READ_OPS1
        }

        memcpy((char *)data + total_size, fh->b->buff + fh->b->offset, size_of_type);

        if (fh->mfooter.change_endianness == adios_flag_yes)
            change_endianness((char *)data + total_size
                             ,size_of_type
                             ,var_root->type
                             );

        if (var_root->type == adios_string)
        {
            // add \0 to the end of string
            // size_of_type here is the length of string
            // FIXME: how would this work for strings written over time?
            ((char*)data + total_size)[size_of_type] = '\0';
        }

        total_size += size_of_type;

        return total_size;
    } /* READ A SCALAR VARIABLE END HERE */

    /* READ AN ARRAY VARIABLE */
    uint64_t write_offset = 0;
    int npg = 0;
    tmpcount = 0;
    int flag;
    datasize = 1;
    nloop = 1;
    var_stride = 1;
    dset_stride = 1;
    uint64_t payload_size = size_of_type;

    /* To get ldims for the index vidx */
    adios_read_bp_staged1_get_dimensioncharacteristics( &(var_root->characteristics[vidx]),
                                                ldims, gdims, offsets);

    /* Again, a Fortran written file has the dimensions in Fortran order we need to swap here */
    /* Only local dims are needed for reading local vars */ 
    if (file_is_fortran)
    {
        i=-1;
        swap_order(ndim, ldims, &(i));
    }

    /*
    printf("ldims   = "); for (j = 0; j < ndim; j++) printf("%d ",ldims[j]); printf("\n");
    printf("count_notime   = "); for (j = 0; j < ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
    printf("start_notime   = "); for (j = 0; j < ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
    */        

    for (j = 0; j < ndim_notime; j++)
    {
        payload_size *= ldims [j];
    
        if ( (start_notime[j] > ldims[j]) 
            || (start_notime[j] + count_notime[j] > ldims[j]))
        {
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], ldims[j] - 1);
                    return adios_errno;
        }
    }

    /* determined how many (fastest changing) dimensions can we read in in one read */
    int break_dim =  ndim_notime - 1;
    while (break_dim > -1)
    {
        if (start_notime[break_dim] == 0 && ldims[break_dim] == count_notime[break_dim])
        {
            datasize *= ldims[break_dim];
        }
        else
            break;
        
        break_dim--;
    }
    
    slice_offset = 0;
    slice_size = 0;
    /* Note: MPI_FILE_READ_OPS  - for reading single BP file.
     *       MPI_FILE_READ_OPS2 - for reading those with subfiles.
     *       MPI_FILE_READ_OPS1 - for reading old version of BP files
     *                            which don't contain "payload_offset"
     * Whenever to use OPS macro, start_idx and idx variable needs to be
     * properly set.
     */
    
    start_idx = 0;
    idx = vidx;

    if (break_dim <= 0) 
    {
        /* The slowest changing dimensions should not be read completely but
           we still need to read only one block */
   
        uint64_t size_in_dset = count_notime[0];
        uint64_t offset_in_dset = start_notime[0];

        slice_size = (break_dim == -1 ? datasize * size_of_type : size_in_dset * datasize * size_of_type);
    
        if (var_root->characteristics[start_idx + idx].payload_offset > 0)
        {
            slice_offset = var_root->characteristics[start_idx + idx].payload_offset 
                         + offset_in_dset * datasize * size_of_type;

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
        }
        else
        {
            slice_offset = 0;
            MPI_FILE_READ_OPS1
        }

        memcpy ((char *)data, fh->b->buff + fh->b->offset, slice_size);
        if (fh->mfooter.change_endianness == adios_flag_yes)
        {
            change_endianness((char *)data + write_offset, slice_size, var_root->type);
        }
    }
    else 
    {
        uint64_t stride_offset = 0;
        uint64_t * size_in_dset, * offset_in_dset, * offset_in_var;
        uint64_t start_in_payload, end_in_payload, s;
        uint64_t var_offset;
        uint64_t dset_offset;

        size_in_dset = (uint64_t *) malloc (8 * ndim_notime);
        offset_in_dset = (uint64_t *) malloc (8 * ndim_notime);
        offset_in_var = (uint64_t *) malloc (8 * ndim_notime);
 
        if (size_in_dset == 0 || offset_in_dset == 0 || offset_in_var == 0)
        {
             adios_error (err_no_memory, "Malloc failed in %s at %d\n"
                         , __FILE__, __LINE__
                         );
             return adios_errno;
        }

        for (i = 0; i < ndim_notime ; i++)
        {
            size_in_dset[i] = count_notime[i];
            offset_in_dset[i] = start_notime[i];
            offset_in_var[i] = 0;
        }
 
        datasize = 1;
        var_stride = 1;
        for (i = ndim_notime - 1; i >= break_dim; i--)
        {
            datasize *= size_in_dset[i];
            dset_stride *= ldims[i];
            var_stride *= count_notime[i];
        }

        /* Calculate the size of the chunk we are trying to read in */
        start_in_payload = 0;
        end_in_payload = 0;
        s = 1;
        for (i = ndim_notime - 1; i >= 0; i--)
        {
            start_in_payload += s * offset_in_dset[i] * size_of_type;
            end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
            s *= ldims[i];
        }
        slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
 
        if (var_root->characteristics[start_idx + idx].payload_offset > 0)
        {
            slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                          + start_in_payload;
            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
 
            for ( i = 0; i < ndim_notime ; i++)
            {
                offset_in_dset[i] = 0;
            }
        }
        else
        {
            slice_offset =  start_in_payload;
            MPI_FILE_READ_OPS1
        }

        var_offset = 0;
        dset_offset = 0;
        for (i = 0; i < ndim_notime ; i++)
        {
            var_offset = offset_in_var[i] + var_offset * count_notime[i];
            dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
        }

        copy_data (data
                  ,fh->b->buff + fh->b->offset
                  ,0
                  ,break_dim
                  ,size_in_dset
                  ,ldims
                  ,count_notime
                  ,var_stride
                  ,dset_stride
                  ,var_offset
                  ,dset_offset
                  ,datasize
                  ,size_of_type 
                  );

        free (size_in_dset);
        free (offset_in_dset);
        free (offset_in_var);
    }
    
    total_size += items_read * size_of_type;

    return total_size;
}

// The purpose of keeping this function is to be able
// to read in old BP files. Can be deleted later on.
int64_t adios_read_bp_staged1_read_var_byid1 (ADIOS_GROUP    * gp,
                             int              varid,
                             const uint64_t  * start,
                             const uint64_t  * count,
                             void           * data)
{
    struct BP_GROUP      * gh;
    BP_FILE       * fh;
    int file_is_fortran;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    int    i,j,k, idx, timestep;
    int    start_time, stop_time;
    int    pgoffset, pgcount, next_pgoffset,start_idx, stop_idx;
    int    ndim, ndim_notime;  
    uint64_t size;
    uint64_t *dims;
    uint64_t ldims[32];
    uint64_t gdims[32];
    uint64_t offsets[32];
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=0, items_read;
    uint64_t count_notime[32], start_notime[32];
    MPI_Status status;
    int timedim = -1, temp_timedim, timedim_c;
    int rank;
    int is_global = 0;
    int size_of_type;
    uint64_t slice_offset;
    uint64_t slice_size;
    uint64_t tmpcount = 0;
    uint64_t datatimeoffset = 0; // offset in data to write a given timestep

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return adios_errno;
    }
    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return adios_errno;
    }
    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return adios_errno;
    }
    if (varid < 0 || varid >= gh->vars_count) {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return adios_errno;
    }
    
    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    
    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i=0; i<varid && var_root; i++) {
        var_root = var_root->next;
    }

    if (i!=varid) {
        adios_error (err_corrupted_variable, "Variable id=%d is valid but was not found in internal data structures!",varid);
        return adios_errno; 
    }

    /* Get dimensions and flip if caller != writer language */
    adios_read_bp_staged1_get_dimensions (var_root, fh->tidx_stop - fh->tidx_start + 1, file_is_fortran, 
                          &ndim, &dims, &timedim);

    /* Here the cases in which .bp written from Fortran and C are considered separately.
       1) bp written from Fortran */
    if (file_is_fortran)
    {
        /* Get the timesteps we need to read */
        if (timedim > -1) 
        {
            if (timedim != ndim - 1)
            {
                adios_error (err_no_data_at_timestep,"Variable (id=%d) has wrong time dimension index",
                      varid);
                return adios_errno;
            }
            if (futils_is_called_from_fortran())
            {
                start_time = start[timedim] + fh->tidx_start;
                stop_time = start_time + count[timedim] - 1;
            }
            else
            {
                start_time = start[0] + fh->tidx_start;
                stop_time = start_time + count[0] - 1;
            }
        }
        else 
        {
            /* timeless variable, but we still need to handle the case that
               var is not written in the first few timesteps. 
               This happens in Pixie3D.  */
            for (i = 0; i < fh->mfooter.time_steps; i++)
            {
                pgoffset = fh->gvar_h->time_index[0][gh->group_id][i];
                if (i < fh->mfooter.time_steps - 1)
                    next_pgoffset = fh->gvar_h->time_index[0][gh->group_id][i + 1];
                else
                    next_pgoffset = -1;

                if (fh->gvar_h->pg_offsets[pgoffset] < var_root->characteristics[0].offset
                && (i == fh->mfooter.time_steps - 1 
                   ||fh->gvar_h->pg_offsets[next_pgoffset] > var_root->characteristics[0].offset)
                )
                {
                    start_time = fh->tidx_start + i;
                    stop_time = start_time;
                    break;
                }
            }
        }

        /* flip dims and timedim to C order */
        swap_order(ndim, dims, &timedim);

        /* Take out the time dimension from start[] and count[] */
        /* if we have time dimension */
        if (timedim > -1)
        {
            j = 0;
            if (futils_is_called_from_fortran())
                temp_timedim = ndim - 1;
            else
                temp_timedim = 0;

            for (i = 0; i < temp_timedim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            i++; // skip timedim
            for (; i < ndim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            ndim_notime = ndim-1;
        }
        else
        /* if we don't have time dimension */
        {
            for (i = 0; i < ndim; i++)
            {
                count_notime[i] = count[i];
                start_notime[i] = start[i];
            }
            ndim_notime = ndim;
        }
    }
    /* 2) .bp written by C */
    else
    {
        /* Get the timesteps we need to read */
        if (timedim > -1) 
        {
            /* timedim has to be the 1st dimension. To be extended to handle 
               the cases timedim at any dimension */
            if (timedim != 0)
            {
                adios_error (err_no_data_at_timestep,"Variable (id=%d) has wrong time dimension",
                      varid);
                return adios_errno;
            }

            if (futils_is_called_from_fortran())
            {
                start_time = start[ndim - 1] + fh->tidx_start;
                stop_time = start_time + count[ndim -1] - 1;
            }
            else
            {
                start_time = start[0] + fh->tidx_start;
                stop_time = start_time + count[0] - 1;
            }

            start_time = start[timedim] + fh->tidx_start;
            stop_time = start_time + count[timedim] - 1;
        }
        else 
        {
            /* timeless variable */
            start_time = fh->tidx_start;
            stop_time = fh->tidx_start;
        }

        /* No need to flip dims, timedim as they are already in C order. */
        //swap_order(ndim, dims, &timedim);

        /* Take out the time dimension from start[] and count[] */
        if (timedim == -1) /* timeless variable */ 
        {
            for (i = 0; i < ndim; i++) 
            {
                count_notime[i] = count[i];
                start_notime[i] = start[i];
            }
            ndim_notime = ndim;
        }
        /* if we have time dimension */
        else
        {
            j = 0;
            if (futils_is_called_from_fortran())
                temp_timedim = ndim - 1;
            else
                temp_timedim = 0;

            for (i = 0; i < temp_timedim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            i++; // skip timedim
            for (; i < ndim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            ndim_notime = ndim - 1;
        }
    }

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if ( futils_is_called_from_fortran()) {
        swap_order(ndim_notime, count_notime, &timedim);
        swap_order(ndim_notime, start_notime, &timedim);
    }
    
    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ndim_notime; i++)
        items_read *= count_notime[i];
    
    MPI_Comm_rank(gh->fh->comm, &rank);

    size_of_type = bp_get_type_size (var_root->type, var_root->characteristics [0].value);

    /* For each timestep, do reading separately (they are stored in different sets of process groups */
    for (timestep = start_time; timestep <= stop_time; timestep++) {

        // pgoffset = the starting offset for the given time step
        // pgcount  = number of process groups of that time step
        pgoffset = fh->gvar_h->time_index[0][gh->group_id][timestep - fh->tidx_start];
        pgcount = fh->gvar_h->time_index[1][gh->group_id][timestep - fh->tidx_start];

        start_idx = -1;
        for (i=0;i<var_root->characteristics_count;i++) {
            if (   (  var_root->characteristics[i].offset > fh->gvar_h->pg_offsets[pgoffset])
                && (  (pgoffset + pgcount == fh->mfooter.pgs_count) 
                    ||(  var_root->characteristics[i].offset < fh->gvar_h->pg_offsets[pgoffset + 1]))
               ) 
            {
                start_idx = i;
                break;
            }
        }
/*
printf ("var_root->characteristics_count = %d\n", var_root->characteristics_count);
printf ("pg offset 0 = %lld\n", fh->gvar_h->pg_offsets[pgoffset]);
printf ("pg offset 1 = %lld\n", fh->gvar_h->pg_offsets[pgoffset + 1]);
printf ("var offset 3 = %lld\n", var_root->characteristics[3].offset);
printf ("var offset 2 = %lld\n", var_root->characteristics[2].offset);
printf ("var offset 1 = %lld\n", var_root->characteristics[1].offset);
printf ("pgcount = %lld\n", pgcount);
*/
        for (i=var_root->characteristics_count-1;i>-1;i--) {
            if (   (  var_root->characteristics[i].offset > fh->gvar_h->pg_offsets[pgoffset])
                && (  (pgoffset + pgcount == fh->mfooter.pgs_count)
                    ||(  var_root->characteristics[i].offset < fh->gvar_h->pg_offsets[pgoffset + pgcount]))
               )
            {
                stop_idx = i;
                break;
            }
        }

        if (start_idx<0) {
            adios_error (err_no_data_at_timestep,"Variable (id=%d) has no data at %d time step",
                varid, timestep);
            return adios_errno;
        }

        if (ndim_notime == 0) {
            /* READ A SCALAR VARIABLE */

            slice_size = size_of_type;
            idx = 0; // macros below need it

            if (var_root->type == adios_string) {
                // strings are stored without \0 in file
                // size_of_type here includes \0 so decrease by one
                size_of_type--;
            }

            if (var_root->characteristics[start_idx+idx].payload_offset > 0) {
                slice_offset = var_root->characteristics[start_idx+idx].payload_offset;
                MPI_FILE_READ_OPS
            } else {
                slice_offset = 0;
                MPI_FILE_READ_OPS1
            }

            memcpy((char *)data+total_size, fh->b->buff + fh->b->offset, size_of_type);
            if (fh->mfooter.change_endianness == adios_flag_yes) {
                change_endianness((char *)data+total_size, size_of_type, var_root->type);
            }

            if (var_root->type == adios_string) {
                // add \0 to the end of string
                // size_of_type here is the length of string
                // FIXME: how would this work for strings written over time?
                ((char*)data+total_size)[size_of_type] = '\0';
            }

            total_size += size_of_type;
            continue;
        }

        /* READ AN ARRAY VARIABLE */
        //int * idx_table = (int *) malloc (sizeof(int) * pgcount);
        //int * idx_table = (int *) malloc (sizeof(int) * (var_root->characteristics_count - start_idx));
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));

        uint64_t write_offset = 0;
        int npg = 0;
        tmpcount = 0;
        if (pgcount > var_root->characteristics_count)
            pgcount = var_root->characteristics_count;

        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++) {
            int flag;
            datasize = 1;
            nloop = 1;
            var_stride = 1;
            dset_stride = 1;
            idx_table[idx] = 1;
            uint64_t payload_size = size_of_type;
    
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = adios_read_bp_staged1_get_dimensioncharacteristics( &(var_root->characteristics[start_idx + idx]),
                                                            ldims, gdims, offsets);
            if (!is_global) {
                // we use gdims below, which is 0 for a local array; set to ldims here
                for (j = 0; j< ndim; j++) {
                    gdims[j]=ldims[j];
                }
            }

            /* Again, a Fortran written file has the dimensions in Fortran order we need to swap here */
            //if (file_is_fortran != futils_is_called_from_fortran()) {
            if (file_is_fortran ) {
                i=-1;
                swap_order(ndim, gdims,   &(i)); // i is dummy 
                swap_order(ndim, ldims,   &(i));
                swap_order(ndim, offsets, &(i));
            }
            
            /* take out the time dimension */
            /* For C, gdims and offset are one size shorter because the timedim part is missing,
               so we take it out only for fortran files
            */
            if (timedim > -1) {
                for (i = timedim; i < ndim-1; i++) {
                    ldims[i] = ldims[i+1];
                    if (file_is_fortran) {
                        gdims[i] = gdims[i+1];
                        offsets[i] = offsets[i+1];
                    }
                }
            }
            /*
            printf("ldims   = "); for (j = 0; j<ndim; j++) printf("%d ",ldims[j]); printf("\n");
            printf("gdims   = "); for (j = 0; j<ndim; j++) printf("%d ",gdims[j]); printf("\n");
            printf("offsets = "); for (j = 0; j<ndim; j++) printf("%d ",offsets[j]); printf("\n");
            printf("count_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
            printf("start_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
            */
            for (j = 0; j < ndim_notime; j++) {
    
                payload_size *= ldims [j];
    
                if ( (count_notime[j] > gdims[j]) 
                  || (start_notime[j] > gdims[j]) 
                  || (start_notime[j] + count_notime[j] > gdims[j])){
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], gdims[j] - 1);
                    return adios_errno;
                }
    
                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start_notime[j] 
                        && offsets[j] < start_notime[j] + count_notime[j])
                    || (offsets[j] < start_notime[j]
                        && offsets[j] + ldims[j] > start_notime[j] + count_notime[j]) 
                    || (offsets[j] + ldims[j] > start_notime[j] 
                        && offsets[j] + ldims[j] <= start_notime[j] + count_notime[j]);
                idx_table [idx] = idx_table[idx] && flag;
            }
            
            if ( !idx_table[idx] ) {
                continue;
            }
            ++npg;

            /* determined how many (fastest changing) dimensions can we read in in one read */
            int hole_break; 
            for (i = ndim_notime - 1; i > -1; i--) {
                if (offsets[i] == start_notime[i] && ldims[i] == count_notime[i]) {
                    datasize *= ldims[i];
                }
                else
                    break;
            }
    
            hole_break = i;
            slice_offset = 0;
            slice_size = 0;

            if (hole_break == -1) {
                /* The complete read happens to be exactly one pg, and the entire pg */
                /* This means we enter this only once, and npg=1 at the end */
                /* This is a rare case. FIXME: cannot eliminate this? */
                slice_size = payload_size;

                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset;
                    MPI_FILE_READ_OPS
                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
    
                memcpy( (char *)data, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness(data, slice_size, var_root->type);
                }
            }
            else if (hole_break == 0) 
            {
                /* The slowest changing dimensions should not be read completely but
                   we still need to read only one block */
                int isize;
                uint64_t size_in_dset = 0;
                uint64_t offset_in_dset = 0;
                uint64_t offset_in_var = 0;
    
                isize = offsets[0] + ldims[0];
                if (start_notime[0] >= offsets[0]) {
                    // head is in
                    if (start_notime[0]<isize) {
                        if (start_notime[0] + count_notime[0] > isize)
                            size_in_dset = isize - start_notime[0];
                        else
                            size_in_dset = count_notime[0];
                        offset_in_dset = start_notime[0] - offsets[0];
                        offset_in_var = 0;
                    }
                }
                else {
                    // middle is in
                    if (isize < start_notime[0] + count_notime[0])
                        size_in_dset = ldims[0];
                    else
                    // tail is in
                        size_in_dset = count_notime[0] + start_notime[0] - offsets[0];
                    offset_in_dset = 0;
                    offset_in_var = offsets[0] - start_notime[0];
                }
    
                slice_size = size_in_dset * datasize * size_of_type;
                write_offset = offset_in_var * datasize * size_of_type;

                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset 
                                 + offset_in_dset * datasize * size_of_type;
                    MPI_FILE_READ_OPS
                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
    
                memcpy ((char *)data + write_offset, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness((char *)data + write_offset, slice_size, var_root->type);
                }
    
                //write_offset +=  slice_size;
            }
            else 
            {

                uint64_t stride_offset = 0;
                int isize;
                uint64_t size_in_dset[10];
                uint64_t offset_in_dset[10];
                uint64_t offset_in_var[10];
                memset(size_in_dset, 0 , 10 * 8);
                memset(offset_in_dset, 0 , 10 * 8);
                memset(offset_in_var, 0 , 10 * 8);
                int hit = 0;
                for ( i = 0; i < ndim_notime ; i++) {
                    isize = offsets[i] + ldims[i];
                    if (start_notime[i] >= offsets[i]) {
                        // head is in
                        if (start_notime[i]<isize) {
                            if (start_notime[i] + count_notime[i] > isize)
                                size_in_dset[i] = isize - start_notime[i];
                            else
                                size_in_dset[i] = count_notime[i];
                            offset_in_dset[i] = start_notime[i] - offsets[i];
                            offset_in_var[i] = 0;
                            hit = 1 + hit * 10;
                        }
                        else
                            hit = -1;
                    }
                    else {
                        // middle is in
                        if (isize < start_notime[i] + count_notime[i]) {
                            size_in_dset[i] = ldims[i];
                            hit = 2 + hit * 10;
                        }
                        else {
                            // tail is in
                            size_in_dset[i] = count_notime[i] + start_notime[i] - offsets[i];
                            hit = 3 + hit * 10;
                        }
                        offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start_notime[i];
                    }
                }
    
                datasize = 1;
                var_stride = 1;
    
                for ( i = ndim_notime-1; i >= hole_break; i--) {
                    datasize *= size_in_dset[i];
                    dset_stride *= ldims[i];
                    var_stride *= count_notime[i];
                }
    
                uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
                for (i = ndim_notime - 1; i > -1; i--) {
                    start_in_payload += s * offset_in_dset[i] * size_of_type;
                    end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
                    s *= ldims[i];
                }
    
                slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                                  + start_in_payload;
                    MPI_FILE_READ_OPS
    
                    for ( i = 0; i < ndim_notime ; i++) {
                        offset_in_dset[i] = 0;
                    }
                } else {
                    slice_offset =  start_in_payload;
                    MPI_FILE_READ_OPS1
                }

                uint64_t var_offset = 0;
                uint64_t dset_offset = 0;
                for ( i = 0; i < hole_break; i++) {
                    nloop *= size_in_dset[i];
                }
    
                for ( i = 0; i < ndim_notime ; i++) {
                    var_offset = offset_in_var[i] + var_offset * count_notime[i];
                    dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
                }

                copy_data (data
                          ,fh->b->buff + fh->b->offset
                          ,0
                          ,hole_break
                          ,size_in_dset
                          ,ldims
                          ,count_notime
                          ,var_stride
                          ,dset_stride
                          ,var_offset
                          ,dset_offset
                          ,datasize
                          ,size_of_type 
                          );

            }
        }  // end for (idx ... loop over pgs
    
        free (idx_table);

        total_size += items_read * size_of_type;
        // shift target pointer for next read in
        data = (char *)data + (items_read * size_of_type);

    } // end for (timestep ... loop over timesteps

    free (dims);

    return total_size;
}

int64_t adios_read_bp_staged1_read_var_byid2 (ADIOS_GROUP    * gp,
                                      int            varid,
                                      const uint64_t * start,
                                      const uint64_t * count,
                                      void           * data)
{
    struct BP_GROUP      * gh;
    BP_FILE       * fh;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    int    i,j,k, idx, t;
    int    start_time, stop_time;
    int    start_idx, stop_idx;
    int    ndim, ndim_notime, has_subfile, file_is_fortran;
    uint64_t size, * dims;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=0, items_read;
    uint64_t count_notime[32], start_notime[32];
    int timedim = -1, temp_timedim, is_global = 0, size_of_type;
    uint64_t slice_offset, slice_size, tmpcount = 0;
    uint64_t datatimeoffset = 0; // offset in data to write a given timestep
    MPI_Status status;

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;

    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    has_subfile = fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE;

    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i = 0; i< varid && var_root; i++) {
        var_root = var_root->next;
    }

    if (i!=varid) {
        adios_error (err_corrupted_variable, 
               "Variable id=%d is valid but was not found in internal data structures!",
               varid);
        return adios_errno; 
    }

    /* Get dimensions and flip if caller != writer language */
    adios_read_bp_staged1_get_dimensions (var_root, fh->tidx_stop - fh->tidx_start + 1, file_is_fortran, 
                          &ndim, &dims, &timedim);

    /* In a Fortran written files, dimensions are in reversed order for C */
    //if ( file_is_fortran != futils_is_called_from_fortran() ) 
    if (file_is_fortran) 
        swap_order(ndim, dims, &timedim);

    /* Take out the time dimension from start[] and count[] */
    if (timedim == -1) {
        /* For timeless var, we still search from fh->tidx_start to fh->tidx_stop
           to handle the situation that some variables are dumped out in selected timesteps
        */
        start_time = fh->tidx_start;
        stop_time = fh->tidx_stop;

        for (i = 0; i < ndim; i++) {
             count_notime[i] = count[i];
             start_notime[i] = start[i];
        }
        ndim_notime = ndim;
    } else {
        j = 0;
        if (futils_is_called_from_fortran())
            temp_timedim = ndim - 1;
        else
            temp_timedim = 0;

        start_time = start[temp_timedim] + fh->tidx_start;
        stop_time = start_time + count[temp_timedim] - 1;

        for (i = 0; i < temp_timedim; i++) {
             count_notime[j] = count[i];
             start_notime[j] = start[i];
             j++;
        }
        i++; // skip timedim
        for (; i < ndim; i++) {
             count_notime[j] = count[i];
             start_notime[j] = start[i];
             j++;
        }
        ndim_notime = ndim-1;
    }

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if ( futils_is_called_from_fortran()) {
        swap_order(ndim_notime, count_notime, &timedim);
        swap_order(ndim_notime, start_notime, &timedim);
    }
    
    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ndim_notime; i++)
        items_read *= count_notime[i];
    
    size_of_type = bp_get_type_size (var_root->type, var_root->characteristics [0].value);

    /* For each timestep, do reading separately (they are stored in different sets of process groups */
    for (t = start_time; t <= stop_time; t++) {
        start_idx = get_var_start_index(var_root, t);
        stop_idx = get_var_stop_index(var_root, t);

        if (start_idx < 0 || stop_idx < 0) {
            adios_error (err_no_data_at_timestep,"Variable (id=%d) has no data at %d time step",
                varid, t);
//            return adios_errno;
            continue;
        }

        if (ndim_notime == 0) {
            /* READ A SCALAR VARIABLE */
            slice_size = size_of_type;
            idx = 0; // macros below need it

            if (var_root->type == adios_string) {
                // strings are stored without \0 in file
                // size_of_type here includes \0 so decrease by one
                size_of_type--;
            }

            if (var_root->characteristics[start_idx+idx].payload_offset > 0) {
                slice_offset = var_root->characteristics[start_idx+idx].payload_offset;
                if (!has_subfile) {
                    MPI_FILE_READ_OPS
                } else {
                    MPI_FILE_READ_OPS2
                }
            } else {
                slice_offset = 0;
                MPI_FILE_READ_OPS1
            }

            memcpy((char *)data+total_size, fh->b->buff + fh->b->offset, size_of_type);

            memcpy((char *)data+total_size, var_root->characteristics[start_idx+idx].value, size_of_type);
            if (fh->mfooter.change_endianness == adios_flag_yes) {
                change_endianness((char *)data+total_size, size_of_type, var_root->type);
            }

            if (var_root->type == adios_string) {
                // add \0 to the end of string
                // size_of_type here is the length of string
                // FIXME: how would this work for strings written over time?
                ((char*)data+total_size)[size_of_type] = '\0';
            }

            total_size += size_of_type;
            
            if (timedim == -1)
                break;
            else
                continue;
        }

         /* READ AN ARRAY VARIABLE */
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));

        uint64_t write_offset = 0;
        int npg = 0;
        tmpcount = 0;
        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++) {
            int flag;
            datasize = 1;
            nloop = 1;
            var_stride = 1;
            dset_stride = 1;
            idx_table[idx] = 1;
            uint64_t payload_size = size_of_type;
    
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = adios_read_bp_staged1_get_dimensioncharacteristics( &(var_root->characteristics[start_idx + idx]),
                                                            ldims, gdims, offsets);
            if (!is_global) {
                // we use gdims below, which is 0 for a local array; set to ldims here
                for (j = 0; j< ndim; j++) {
                    gdims[j]=ldims[j];
                }
            }

            /* Again, a Fortran written file has the dimensions in Fortran order we need to swap here */
            //if (file_is_fortran != futils_is_called_from_fortran()) {
            if (file_is_fortran) {
                i=-1;
                swap_order(ndim, gdims,   &(i)); // i is dummy 
                swap_order(ndim, ldims,   &(i));
                swap_order(ndim, offsets, &(i));
            }
            
            /* take out the time dimension */
            /* For C, gdims and offset are one size shorter because the timedim part is missing,
               so we take it out only for fortran files
            */
            if (timedim > -1) {
                for (i = timedim; i < ndim-1; i++) {
                    ldims[i] = ldims[i+1];
                    if (file_is_fortran) {
                        gdims[i] = gdims[i+1];
                        offsets[i] = offsets[i+1];
                    }
                }
            }

            /*
            printf("ldims   = "); for (j = 0; j<ndim; j++) printf("%d ",ldims[j]); printf("\n");
            printf("gdims   = "); for (j = 0; j<ndim; j++) printf("%d ",gdims[j]); printf("\n");
            printf("offsets = "); for (j = 0; j<ndim; j++) printf("%d ",offsets[j]); printf("\n");
            printf("count_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
            printf("start_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
            */
                
            for (j = 0; j < ndim_notime; j++) {
    
                payload_size *= ldims [j];
    
                if ( (count_notime[j] > gdims[j]) 
                  || (start_notime[j] > gdims[j]) 
                  || (start_notime[j] + count_notime[j] > gdims[j])){
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], gdims[j] - 1);
                    return adios_errno;
                }
    
                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start_notime[j] 
                        && offsets[j] < start_notime[j] + count_notime[j])
                    || (offsets[j] < start_notime[j]
                        && offsets[j] + ldims[j] > start_notime[j] + count_notime[j]) 
                    || (offsets[j] + ldims[j] > start_notime[j] 
                        && offsets[j] + ldims[j] <= start_notime[j] + count_notime[j]);
                idx_table [idx] = idx_table[idx] && flag;
            }
            
            if ( !idx_table[idx] ) {
                continue;
            }
            ++npg;

            /* determined how many (fastest changing) dimensions can we read in in one read */
            int hole_break; 
            for (i = ndim_notime - 1; i > -1; i--) {
                if (offsets[i] == start_notime[i] && ldims[i] == count_notime[i]) {
                    datasize *= ldims[i];
                }
                else
                    break;
            }
    
            hole_break = i;
            slice_offset = 0;
            slice_size = 0;

            if (hole_break == -1) {
                /* The complete read happens to be exactly one pg, and the entire pg */
                /* This means we enter this only once, and npg=1 at the end */
                /* This is a rare case. FIXME: cannot eliminate this? */
                slice_size = payload_size;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset;
                    if (!has_subfile) {
                        MPI_FILE_READ_OPS
                    } else {
                        MPI_FILE_READ_OPS2
                    }
                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
 
                memcpy( (char *)data, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness(data, slice_size, var_root->type);
                }
            }
            else if (hole_break == 0) 
            {
                /* The slowest changing dimensions should not be read completely but
                   we still need to read only one block */
                int isize;
                uint64_t size_in_dset = 0;
                uint64_t offset_in_dset = 0;
    
                isize = offsets[0] + ldims[0];
                if (start_notime[0] >= offsets[0]) {
                    // head is in
                    if (start_notime[0]<isize) {
                        if (start_notime[0] + count_notime[0] > isize)
                            size_in_dset = isize - start_notime[0];
                        else
                            size_in_dset = count_notime[0];
                        offset_in_dset = start_notime[0] - offsets[0];
                    }
                }
                else {
                    // middle is in
                    if (isize < start_notime[0] + count_notime[0])
                        size_in_dset = ldims[0];
                    else
                    // tail is in
                        size_in_dset = count_notime[0] + start_notime[0] - offsets[0];
                    offset_in_dset = 0;
                }
    
                slice_size = size_in_dset * datasize * size_of_type;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset 
                                 + offset_in_dset * datasize * size_of_type;
                    if (!has_subfile) {
                        MPI_FILE_READ_OPS
                    } else {
                        MPI_FILE_READ_OPS2
                    }

                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
    
                memcpy ((char *)data + write_offset, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness((char *)data + write_offset, slice_size, var_root->type);
                }
    
                write_offset +=  slice_size;
            }
            else 
            {

                uint64_t stride_offset = 0;
                int isize;
                uint64_t size_in_dset[10];
                uint64_t offset_in_dset[10];
                uint64_t offset_in_var[10];
                memset(size_in_dset, 0 , 10 * 8);
                memset(offset_in_dset, 0 , 10 * 8);
                memset(offset_in_var, 0 , 10 * 8);
                int hit = 0;
                for ( i = 0; i < ndim_notime ; i++) {
                    isize = offsets[i] + ldims[i];
                    if (start_notime[i] >= offsets[i]) {
                        // head is in
                        if (start_notime[i]<isize) {
                            if (start_notime[i] + count_notime[i] > isize)
                                size_in_dset[i] = isize - start_notime[i];
                            else
                                size_in_dset[i] = count_notime[i];
                            offset_in_dset[i] = start_notime[i] - offsets[i];
                            offset_in_var[i] = 0;
                            hit = 1 + hit * 10;
                        }
                        else
                            hit = -1;
                    }
                    else {
                        // middle is in
                        if (isize < start_notime[i] + count_notime[i]) {
                            size_in_dset[i] = ldims[i];
                            hit = 2 + hit * 10;
                        }
                        else {
                            // tail is in
                            size_in_dset[i] = count_notime[i] + start_notime[i] - offsets[i];
                            hit = 3 + hit * 10;
                        }
                        offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start_notime[i];
                    }
                }
    
                datasize = 1;
                var_stride = 1;
    
                for ( i = ndim_notime-1; i >= hole_break; i--) {
                    datasize *= size_in_dset[i];
                    dset_stride *= ldims[i];
                    var_stride *= count_notime[i];
                }
    
                uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
                for (i = ndim_notime - 1; i > -1; i--) {
                    start_in_payload += s * offset_in_dset[i] * size_of_type;
                    end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
                    s *= ldims[i];
                }
    
                slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                                  + start_in_payload;
                    if (!has_subfile) {
                        MPI_FILE_READ_OPS
                    } else {
                        MPI_FILE_READ_OPS2
                    }
 
                    for ( i = 0; i < ndim_notime ; i++) {
                        offset_in_dset[i] = 0;
                    }
                } else {
                    slice_offset =  start_in_payload;
                    MPI_FILE_READ_OPS1
                }
    
                uint64_t var_offset = 0;
                uint64_t dset_offset = 0;
                for ( i = 0; i < hole_break; i++) {
                    nloop *= size_in_dset[i];
                }
    
                for ( i = 0; i < ndim_notime ; i++) {
                    var_offset = offset_in_var[i] + var_offset * count_notime[i];
                    dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
                }
    
                copy_data (data
                          ,fh->b->buff + fh->b->offset
                          ,0
                          ,hole_break
                          ,size_in_dset
                          ,ldims
                          ,count_notime
                          ,var_stride
                          ,dset_stride
                          ,var_offset
                          ,dset_offset
                          ,datasize
                          ,size_of_type 
                          );
            }
        }  // end for (idx ... loop over pgs
    
        free (idx_table);
    
        total_size += items_read * size_of_type;
        // shift target pointer for next read in
        data = (char *)data + (items_read * size_of_type);

        if (timedim == -1)
            break;
    } // end for (timestep ... loop over timesteps

    free (dims);

    return total_size;
}

int64_t adios_read_bp_staged1_read_var_byid (ADIOS_GROUP    * gp,
                                     int            varid,
                                     const uint64_t  * start,
                                     const uint64_t  * count,
                                     void            * data)
{
    struct BP_GROUP      * gh;
    BP_FILE       * fh;
    int has_time_index_characteristic;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return adios_errno;
    }

    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return adios_errno;
    }

    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return adios_errno;
    }

    has_time_index_characteristic = fh->mfooter.version & ADIOS_VERSION_HAVE_TIME_INDEX_CHARACTERISTIC;
    if (!has_time_index_characteristic) {
        // read older file format. Can be deleted later on.
        return adios_read_bp_staged1_read_var_byid1(gp, varid, start, count, data);
    } else {
        return adios_read_bp_staged1_read_var_byid2(gp, varid, start, count, data);
    }
}
#endif
