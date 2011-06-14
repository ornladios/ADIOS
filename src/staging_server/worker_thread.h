#ifndef __WORKER_THREAD_H__
#define __WORKER_THREAD_H__

/* Worker thread should tell if
   1. index is needed or not
   2. size of rdma buffer used for client data
   This is called in the initial main thread, 
   NOT from the worker thread!
*/
void worker_global_init (void);

#if HAVE_PTHREAD
/* Worker thread entry point, if we have threads */
void * worker_thread_main (void *arg);
#else
/* Non-thread env: transport loop calls this function regularly 
   Returns: how many events the worker 'thread' processed
*/
int worker_main (void *arg);
void worker_finalize (void);
#endif

#endif
