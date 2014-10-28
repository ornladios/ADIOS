/* Auto-generated on Sat Sep 13 16:09:17 EDT 2014 */

#ifndef ADIOS_ICEE_H
#define ADIOS_ICEE_H

typedef struct icee_varinfo_rec {
    char* varname;
    int varid;
    int type;
    int typesize;
    int ndims;
    uint64_t* gdims;
    uint64_t* ldims;
    uint64_t* offsets;
    uint64_t varlen;
    char* data;
    struct icee_varinfo_rec * next;
} icee_varinfo_rec_t, *icee_varinfo_rec_ptr_t;

static FMField icee_varinfo_field_list[] =
{
    {"varname", "string", sizeof(char*), FMOffset(icee_varinfo_rec_ptr_t, varname)},
    {"varid", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, varid)},
    {"type", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, type)},
    {"typesize", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, typesize)},
    {"ndims", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, ndims)},
    {"gdims", "integer[ndims]", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, gdims)},
    {"ldims", "integer[ndims]", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, ldims)},
    {"offsets", "integer[ndims]", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, offsets)},
    {"varlen", "integer", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, varlen)},
    {"data", "char[varlen]", sizeof(char), FMOffset(icee_varinfo_rec_ptr_t, data)},
    {"next", "*icee_varinfo", sizeof(struct icee_varinfo_rec ), FMOffset(icee_varinfo_rec_ptr_t, next)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec icee_varinfo_format_list[] =
{
    {"icee_varinfo", icee_varinfo_field_list, sizeof(icee_varinfo_rec_t), NULL},
    {NULL, NULL}
};

/* Auto-generated on Sat Sep 13 16:09:17 EDT 2014 */

typedef struct icee_fileinfo_rec {
    char* fname;
    int nvars;
    int nchunks;
    int comm_size;
    int comm_rank;
    int merge_count;
    int timestep;
    struct icee_varinfo_rec * varinfo;
    struct icee_fileinfo_rec * next;
} icee_fileinfo_rec_t, *icee_fileinfo_rec_ptr_t;

static FMField icee_fileinfo_field_list[] =
{
    {"fname", "string", sizeof(char*), FMOffset(icee_fileinfo_rec_ptr_t, fname)},
    {"nvars", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, nvars)},
    {"nchunks", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, nchunks)},
    {"comm_size", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, comm_size)},
    {"comm_rank", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, comm_rank)},
    {"merge_count", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, merge_count)},
    {"timestep", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, timestep)},
    {"varinfo", "*icee_varinfo", sizeof(struct icee_varinfo_rec ), FMOffset(icee_fileinfo_rec_ptr_t, varinfo)},
    {"next", "*icee_fileinfo", sizeof(struct icee_fileinfo_rec ), FMOffset(icee_fileinfo_rec_ptr_t, next)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec icee_fileinfo_format_list[] =
{
    {"icee_fileinfo", icee_fileinfo_field_list, sizeof(icee_fileinfo_rec_t), NULL},
    {"icee_varinfo", icee_varinfo_field_list, sizeof(icee_varinfo_rec_t), NULL},
    {NULL, NULL}
};

/* Auto-generated on Sat Sep 13 16:09:17 EDT 2014 */

typedef struct icee_clientinfo_rec {
    char* client_host;
    int client_port;
    int stone_id;
} icee_clientinfo_rec_t, *icee_clientinfo_rec_ptr_t;

static FMField icee_clientinfo_field_list[] =
{
    {"client_host", "string", sizeof(char*), FMOffset(icee_clientinfo_rec_ptr_t, client_host)},
    {"client_port", "integer", sizeof(int), FMOffset(icee_clientinfo_rec_ptr_t, client_port)},
    {"stone_id", "integer", sizeof(int), FMOffset(icee_clientinfo_rec_ptr_t, stone_id)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec icee_clientinfo_format_list[] =
{
    {"icee_clientinfo", icee_clientinfo_field_list, sizeof(icee_clientinfo_rec_t), NULL},
    {NULL, NULL}
};

/*
 * Thread pool implementation
 * Credit: Multithreaded Programming Guide by Oracle
 * http://docs.oracle.com/cd/E19253-01/816-5137/6mba5vqn3/index.html
 */

/*
 * Declarations for the clients of a thread pool.
 */

#include <pthread.h>
#include <sys/time.h>

typedef unsigned int uint_t;

/*
 * The thr_pool_t type is opaque to the client.
 * It is created by thr_pool_create() and must be passed
 * unmodified to the remainder of the interfaces.
 */
typedef	struct thr_pool	thr_pool_t;

/*
 * Create a thread pool.
 *	min_threads:	the minimum number of threads kept in the pool,
 *			always available to perform work requests.
 *	max_threads:	the maximum number of threads that can be
 *			in the pool, performing work requests.
 *	linger:		the number of seconds excess idle worker threads
 *			(greater than min_threads) linger before exiting.
 *	attr:		attributes of all worker threads (can be NULL);
 *			can be destroyed after calling thr_pool_create().
 * On error, thr_pool_create() returns NULL with errno set to the error code.
 */
extern	thr_pool_t	*thr_pool_create(uint_t min_threads, uint_t max_threads,
				uint_t linger, pthread_attr_t *attr);

/*
 * Enqueue a work request to the thread pool job queue.
 * If there are idle worker threads, awaken one to perform the job.
 * Else if the maximum number of workers has not been reached,
 * create a new worker thread to perform the job.
 * Else just return after adding the job to the queue;
 * an existing worker thread will perform the job when
 * it finishes the job it is currently performing.
 *
 * The job is performed as if a new detached thread were created for it:
 *	pthread_create(NULL, attr, void *(*func)(void *), void *arg);
 *
 * On error, thr_pool_queue() returns -1 with errno set to the error code.
 */
extern	int	thr_pool_queue(thr_pool_t *pool,
			void *(*func)(void *), void *arg);

/*
 * Wait for all queued jobs to complete.
 */
extern	void	thr_pool_wait(thr_pool_t *pool);

/*
 * Cancel all queued jobs and destroy the pool.
 */
extern	void	thr_pool_destroy(thr_pool_t *pool);

#endif
