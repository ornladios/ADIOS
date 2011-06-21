#include <stdlib.h>  // malloc
#include "mymutex.h"
#include "precedence.h"

//static const char tabs[] = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";

struct precedence_struct {
    pthread_mutex_t mutex;
    pthread_cond_t  cond_var;
    int max_precedence;
    int current_precedence;
};

struct precedence_struct * precedence_create (int max_prec)
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

void precedence_get (struct precedence_struct *p, int my_prec)
{
    if (!p) return;
    //printf("%.*s%d:   lock...\n", 2*my_prec, tabs, my_prec);
    pthread_mutex_lock (&p->mutex);
    //printf("%.*s%d:   locked.\n", 2*my_prec, tabs, my_prec);
    while (p->current_precedence != my_prec) {
        // wait until thread with lower prec makes a get/release
        //printf("%.*s%d:   wait...\n", 2*my_prec, tabs, my_prec);
        pthread_cond_wait (&p->cond_var, &p->mutex);
        //printf("%.*s%d:   woke up.\n", 2*my_prec, tabs, my_prec);
    } 
    // this thread has now control until it calls release
    //pthread_mutex_unlock (&p->mutex);
    p->current_precedence = my_prec;
    return;
}

void precedence_release (struct precedence_struct *p)
{
    if (!p) return;
    //pthread_mutex_lock (&p->mutex);
    int my_prec = p->current_precedence;
    if (p->current_precedence < p->max_precedence) {
        p->current_precedence++;
    } else {
        p->current_precedence = 0;
    }
    // wake up those who still wait (at next unlock)
    //printf("%.*s%d:   signal to prec=%d...\n", 2*my_prec, tabs, my_prec, p->current_precedence);
    pthread_cond_broadcast (&p->cond_var);
    //printf("%.*s%d:   unlock...\n", 2*my_prec, tabs, my_prec);
    pthread_mutex_unlock (&p->mutex);
    //printf("%.*s%d:   unlocked.\n", 2*my_prec, tabs, my_prec);
}


