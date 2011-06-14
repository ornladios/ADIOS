#ifndef __PRECEDENCE_H__
#define __PRECEDENCE_H__

#include <stdlib.h>  // malloc
#include "mymutex.h"
//#include "ardma_logger.h"

struct precedence_struct {
    pthread_mutex_t mutex;
    pthread_cond_t  cond_var;
    int max_precedence;
    int current_precedence;
};

inline struct precedence_struct * precedence_create (int max_prec)
{
    struct precedence_struct *p = (struct precedence_struct *) 
                            malloc (sizeof (struct precedence_struct));
    if (p) { 
        pthread_mutex_init(&p->mutex, NULL);
        pthread_cond_init (&p->cond_var, NULL);
        p->max_precedence = max_prec;
        p->current_precedence = 0;
    }
    return p;
}

inline void precedence_get (struct precedence_struct *p, int my_prec)
{
    if (!p) return;
    pthread_mutex_lock (&p->mutex);
    while (p->current_precedence != my_prec) {
        // wait until thread with lower prec makes a get/release
        pthread_cond_wait (&p->cond_var, &p->mutex);
    } 
    /* 
    if (p->current_precedence > my_prec) {
        // logical error, this guy should have called before
        log_error ("ERROR in precedence: thread called with precedence %d "
                   "when precedence level is already %d. "
                   "Behavior is undefined from this point.\n",
                   my_prec, p->current_precedence);
    };
    */
    // this thread has now control until it calls release
    p->current_precedence = my_prec;
    return;
}

inline void precedence_release (struct precedence_struct *p)
{
    if (!p) return;
    if (p->current_precedence < p->max_precedence) {
        p->current_precedence++;
        // wake up those who still wait (at next unlock)
        pthread_cond_broadcast (&p->cond_var);
    } else {
        p->current_precedence = 0;
    }
    pthread_mutex_unlock (&p->mutex);
}


#endif
