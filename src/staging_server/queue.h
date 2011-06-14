/*************************************************************************
* FILE:     queue.h
* DESCRIPTION: Header file for a basic queue implementation.
*************************************************************************/
#ifndef __QUEUE_H__
#define __QUEUE_H__

/*******************************************************************************
* Includes
*******************************************************************************/
#include "config.h"
#include <stdio.h>
#include <stdlib.h>

#include "mymutex.h"
/*
#if HAVE_PTHREAD
#  include "pthread.h"
#else
# define pthread_mutex_t int
# define pthread_mutexattr_t int
inline int pthread_mutex_init(pthread_mutex_t *mutex, pthread_mutexattr_t *attr) { return 0;}
inline int pthread_mutex_destroy(pthread_mutex_t *mutex) { return 0;}
inline int pthread_mutex_lock(pthread_mutex_t *mutex) { return 0;}
inline int pthread_mutex_trylock(pthread_mutex_t *mutex) { return 0;}
inline int pthread_mutex_unlock(pthread_mutex_t *mutex) { return 0;}
#endif
*/

/*******************************************************************************
* Queue Constants
*******************************************************************************/
extern const int QUEUE_UNBOUNDED;


/*************************************************************************
* Queue Datastructures
*************************************************************************/
struct queue_element
{
    void                       *data;
    struct queue_element *next;
};
typedef struct queue_element queue_element_t;

struct queue
{
    int                    depth;
    int                    max_depth;
    queue_element_t *head;
    queue_element_t *tail;
    pthread_mutex_t        mutex;
};
typedef struct queue queue_t;


/*************************************************************************
* Function Prototypes
*************************************************************************/

/*====================================================================
* queue_init() - Create a thread-safe queue. Setting the max_depth=-1
*    means the queue size is unbounded.
*    WARNING: gen_pthread_init() must have been previously called.
**==================================================================*/
queue_t *queue_init(int max_depth);


/*====================================================================
* enqueue() - thread safe enqueue of data.
*    RETURNS: 0=success, 1=queue_full 
**==================================================================*/
int enqueue(queue_t * q, void *data);


/*====================================================================
* dequeue() - thread safe data dequeue.
*    RETURNS: the queued data or NULL on an empty queue.
**==================================================================*/
void *dequeue(queue_t * q);


/*====================================================================
* queue_destroy() - frees the specified queue datastructure. Elements
*    within the queue are NOT explicitly freed. Empty the queue if
*    these are your only pointers to the queued data!
**==================================================================*/
void queue_destroy(queue_t *);

int  is_queue_empty(queue_t *);

void *peakqueue(queue_t * q);

#endif                                                           /* QUEUE_H */
/*******************************************************************************
* End of File
*******************************************************************************/
