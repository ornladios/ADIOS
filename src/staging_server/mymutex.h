#ifndef __MYMUTEX_H__
#define __MYMUTEX_H__

#include "config.h"

#if HAVE_PTHREAD
#   include "pthread.h"
#else
#   warning DEFINE FAKE PHTREAD MUTEX AND COND_VARS
#   define pthread_mutex_t int
#   define pthread_mutexattr_t int
    inline int pthread_mutex_init(pthread_mutex_t *mutex, pthread_mutexattr_t *attr) { return 0;}
    inline int pthread_mutex_destroy(pthread_mutex_t *mutex) { return 0;}
    inline int pthread_mutex_lock(pthread_mutex_t *mutex) { return 0;}
    inline int pthread_mutex_trylock(pthread_mutex_t *mutex) { return 0;}
    inline int pthread_mutex_unlock(pthread_mutex_t *mutex) { return 0;}
#   define pthread_cond_t int
#   define pthread_condattr_t int
    inline int pthread_cond_init(pthread_cond_t *cond, pthread_condattr_t *attr) {return 0;}
    inline int pthread_cond_signal(pthread_cond_t *cond) {return 0;}
    inline int pthread_cond_broadcast(pthread_cond_t *cond) {return 0;}
    inline int pthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex) {return 0;}
#endif


#endif
