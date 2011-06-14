#define _LARGEFILE64_SOURCE
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>

#include "config.h"

#include "bpwriter_thread.h"
#include "globals.h"
#include "precedence.h"


/* bpwriter thread is started from bpworker thread
   In a non-threaded environment, bpwriter_main() is called from bpworker.
*/

static int fd; /* output file */
static int fd_idx; /* rank 0 writes extra index file */
static off64_t current_offset; /* remember latest offset to avoid unnecessary seeks */
static char *fpath;    /* Remember file path for error messages */
static char *basepath; /* Remember base file path for rank 0 writing the index */
static bool writer_inited = false;
queue_t *woq; // Worker to Writer queue, extern to bpworker */
queue_t *wrq; // Writer to Worker queue, extern to bpworker */
//struct writer_info_struct writer_info; /* extern to bpworker */
static bool should_exit = false; // Worker has requested exit -> true

#if HAVE_PTHREAD
void bpwriter_init (void);
void bpwriter_finalize(); 
int bpwriter_main(void);
void * bpwriter_thread_main (void *arg)
{   
    int num_event = 1;  //1:ensure entering loop below and avoid first sleep
    struct timespec treq = {.tv_sec=0, .tv_nsec=2000000}; // 2ms sleeptime
    struct timespec trem;
    bpwriter_init();
    while (!gd.terminate && !should_exit) {
        if (num_event == 0)
            nanosleep(&treq, &trem);
        num_event = 1;
        while (num_event) { // complete all outstanding writing requests
            num_event = bpwriter_main();
        }
    }
    if (gd.terminate) {
        log_debug("rank %d: writer thread: exit due to failure.\n", gd.mpi_rank);
    } else {
        log_debug("rank %d: writer thread: exit due to request.\n", gd.mpi_rank);
    }

    bpwriter_finalize();
    pthread_exit(NULL);
    return NULL; // just to avoid compiler warning
}
#endif

void bpwriter_init (void)
{
    if (writer_inited)
        return;
    /*
    int x = pthread_mutex_init(&writer_info.mutex, NULL);
    if (x != 0) {   
        log_error("rank %0: ERROR: writer info mutex initialization failed.", gd.mpi_rank);
    }
    pthread_mutex_lock(&writer_info.mutex);
    writer_info.n_outstanding_requests = 0;
    pthread_mutex_unlock(&writer_info.mutex);
    */
    writer_inited = true;
    fpath = NULL;
    basepath = NULL;
    fd = -1;
    fd_idx = -1;
    should_exit = false;
}

void bpwriter_finalize() {
    if (fpath) free(fpath);
    if (basepath) free(basepath);
}

void bpwriter_open (struct writer_request_open * o);
int bpwriter_write (int fd, uint64_t addr, uint64_t size, uint64_t offset);

int bpwriter_main(void)
{
    int num_event = 0;
    if (!writer_inited) // non-threaded env: transport calls this function only
        bpwriter_init();


    struct writer_request *wreq = dequeue(woq);
    if (wreq != NULL) {
        switch (wreq->type) {

            case WRITER_REQUEST_OPEN:
                bpwriter_open (&wreq->request.open);
                /*
                pthread_mutex_lock(&writer_info.mutex);
                writer_info.n_outstanding_requests--;
                pthread_mutex_unlock(&writer_info.mutex);
                */
                break;

            case WRITER_REQUEST_WRITE:
                if (fd > 0) {
                    log_info("rank %d: write block size=%lld offs=%lld\n", 
                            gd.mpi_rank,
                            wreq->request.write.size,
                            wreq->request.write.offset);

                    bpwriter_write (fd, wreq->request.write.addr, 
                                    wreq->request.write.size,
                                    wreq->request.write.offset);
                }
                if (wreq->request.write.pso != NULL) {
                    // 'signal' completion to sender
                    log_debug("rank %d: writer set pso=%lld to weo=%lld\n", gd.mpi_rank,
                              *(wreq->request.write.pso), wreq->request.write.weo);
                    *(wreq->request.write.pso) = wreq->request.write.weo;
                }
                /*
                pthread_mutex_lock(&writer_info.mutex);
                writer_info.n_outstanding_requests--;
                pthread_mutex_unlock(&writer_info.mutex);
                */
                break;

            case WRITER_REQUEST_WRITEIDX:
                /* append local index to end of file */
                if (fd > 0) {
                    log_info("rank %d: write index size=%lld offs=%lld\n", 
                            gd.mpi_rank, 
                            wreq->request.writeidx.lsize,
                            wreq->request.writeidx.loffset);

                    bpwriter_write (fd, wreq->request.writeidx.laddr, 
                                    wreq->request.writeidx.lsize,
                                    wreq->request.writeidx.loffset);
                }
                // free buffer allocated in indexing thread (instead of sending notifications around)
                free ((char *)wreq->request.writeidx.laddr);

                /* on rank 0, write global metadata file */
                if (gd.mpi_rank == 0) {
                    if (fd_idx > 0) {
                        log_info("rank %d: write global index size=%lld to %s\n", 
                                gd.mpi_rank, wreq->request.writeidx.lsize, basepath);
                        ssize_t n = write (fd_idx, (void *)wreq->request.writeidx.gaddr, 
                                           (size_t)wreq->request.writeidx.gsize);
                        if (n < 0) {
                            log_error("rank %d: I/O ERROR: Failed to write block to file %s, "
                                    "size=%lld, addr=%lld, file offset=0, errno=%d: %s\n",
                                    gd.mpi_rank, basepath, wreq->request.writeidx.gsize, 
                                    wreq->request.writeidx.gaddr, errno, strerror(errno));
                        }
                        close(fd_idx);
                    }
                    free ((char *)wreq->request.writeidx.gaddr);
                    /*
                    pthread_mutex_lock(&writer_info.mutex);
                    writer_info.n_outstanding_requests--;
                    pthread_mutex_unlock(&writer_info.mutex);
                    */
                } 
                break;

            case WRITER_REQUEST_CLOSE:
                if (fd > 0) close(fd);
                if (fd_idx > 0) close(fd_idx);
                /*
                pthread_mutex_lock(&writer_info.mutex);
                writer_info.n_outstanding_requests--;
                pthread_mutex_unlock(&writer_info.mutex);
                */
                break;

            case WRITER_REQUEST_EXIT:
                should_exit = true;
                break;

            default:
                log_error ("rank %d: Writer ERROR unexpected request from worker, event type=%d\n",
                        gd.mpi_rank, wreq->type);
                break;
        }
        free(wreq);
        num_event++;
    }
    return num_event;
}


void bpwriter_open (struct writer_request_open * o)
{
    log_debug("rank %d: %s: o=%x, path=%s\n", gd.mpi_rank, __func__, o, o->path);
    char path[MAXPATHLEN];
    char *dirname, *basename;
    getbasename (o->path, &dirname, &basename); // need the file name only
    
    // create path.dir/filename.<rank>
    strncpy (path, o->path, MAXPATHLEN);
    strncat (path, ".dir", MAXPATHLEN);
    int res;
    if (gd.mpi_rank == 0) {
        log_info("rank %d: create dir %s\n", gd.mpi_rank, path);
        res = mkdir(path, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
        if (res) {
            if (errno != EEXIST) {
                log_error("rank %d: mkdir %s failed: %s\n", path, strerror(errno));
                res = errno;
            } else {
                res = 0; // directory already exist, not a problem
            }
        }
    }
    // broadcast res(=errno) from rank 0 to all 
    precedence_get(gd.prec_writer_indexing, 0);
    MPI_Comm newcomm;
    MPI_Comm_dup(gd.mpi_comm, &newcomm);
    MPI_Bcast (&res, 1, MPI_INT, 0, newcomm);
    precedence_release(gd.prec_writer_indexing);

    if (!res) {
        // add filename to path
        char fname[MAXPATHLEN];
        snprintf (fname, MAXPATHLEN, "%s/%s.%d", path, basename, gd.mpi_rank);
        log_info("rank %d: create file %s\n", gd.mpi_rank, fname);
        fd  = open (fname, 
                O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH
                );
        current_offset = 0;

        // rank 0: create the metadata file with o->path as name
        if (gd.mpi_rank == 0) {
            log_info("rank %d: create metadata file %s\n", gd.mpi_rank, o->path);
            fd_idx  = open (o->path, 
                    O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                    S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH
                    );
        }

    } else {
        fd = -1;
        fd_idx = -1;
    }

    /* remember file names */
    if (fpath) free(fpath);
    if (basepath) free(basepath);
    fpath = strdup (path);
    basepath = strdup (o->path);
}


int bpwriter_write (int fd, uint64_t addr, uint64_t size, uint64_t offset)
{
    int ret;
    if (addr > 0 && size > 0) {
        ret=0;
        if (current_offset != (off64_t) offset) {
            lseek64 (fd, offset, SEEK_SET);
            log_warn("rank %d: Writer Warning: had to seek from current "
                    "offset %lld to offset %lld. The implementation should "
                    "ensure there is not seeks, so this indicates some error\n",
                    gd.mpi_rank, current_offset, offset);
        }
        ssize_t n = write (fd, (void *)addr, (size_t)size);
        if (n < 0) {
            ret = errno;
            log_error("rank %d: I/O ERROR: Failed to write block to file %s, "
                      "size=%lld, addr=%lld, file offset=%lld, errno=%d: %s\n",
                      gd.mpi_rank, fpath, size, addr, offset, errno, strerror(errno));
        }
        current_offset = (off64_t)(offset + size);
    } else {
        ret=1;
    }
    return ret;
}

