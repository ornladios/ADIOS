/******************************************************************************
* FILE:        queue.c
* DESCRIPTION: A basic queue implementation.  The operations are thread safe.
******************************************************************************/
#include "queue.h"

/*******************************************************************************
* Constants
*******************************************************************************/
const int  QUEUE_UNBOUNDED = -1;

/*******************************************************************************
* Implementation
*******************************************************************************/

/*====================================================================
* queue_init()
**==================================================================*/
queue_t *
queue_init(int max_depth)
{
    queue_t *q;
    
    if ((max_depth < 1) && (max_depth != QUEUE_UNBOUNDED))
    {
        fprintf(stderr, "ERROR: queue_init called with [%d] as argument.\n", max_depth);
        exit(1);
    }

    q = (queue_t *) malloc(sizeof(queue_t));
    if (q == NULL)
    {
        fprintf(stderr, "ERROR: unable to malloc in init_queue.\n");
        exit(1);
    }

    q->depth = 0;
    q->max_depth = max_depth;
    q->head = NULL;
    q->tail = NULL;
    int x = pthread_mutex_init(&q->mutex, NULL);
    if (x != 0)
    {
        fprintf(stderr, "ERROR: queue mutex initialization failed.");
        exit(1);
    }

    return (q);
}


/*====================================================================
* enqueue()
**==================================================================*/
int
enqueue(queue_t * q, void *data)
{
    queue_element_t *qe;
    int             retval = 0;

    if (q == NULL)
    {
        fprintf(stderr, "ERROR: attempted enqueue to NULL queue.\n");
        exit(1);
    }
    if (data == NULL)
    {
        fprintf(stderr, "ERROR: attempted to enqueue NULL data.\n");
        exit(1);
    }

    qe = (queue_element_t *) malloc(sizeof(queue_element_t));
    if (qe == NULL)
    {
        fprintf(stderr, "ERROR: unable to malloc for enqueue.\n");
        exit(1);
    }
    qe->data = data;
    qe->next = NULL;

    pthread_mutex_lock(&q->mutex);
    if ((q->max_depth == QUEUE_UNBOUNDED) || (q->depth < q->max_depth))
    {
        if (q->head == NULL)
        {
            q->head = qe;
        }
        else
        {
            q->tail->next = qe;
        }
        q->tail = qe;
        q->depth++;
    }
    else
    {
        // queue full
        retval = 1;
    }
    pthread_mutex_unlock(&q->mutex);

    return retval;
}


/*====================================================================
* dequeue()
**==================================================================*/
void *
dequeue(queue_t * q)
{
    queue_element_t *qe;
    void *data;

    if (q == NULL)
    {
        fprintf(stderr, "ERROR: dequeue called on NULL queue.\n");
        exit(1);
    }

    qe = NULL;

    pthread_mutex_lock(&q->mutex);
    if (q->head == NULL)
    {
        // empty queue
        data = NULL;
    }
    else
    {
        // queue not empty
        qe = q->head;
        data = qe->data;

        q->head = q->head->next;
        if (q->head == NULL)
        {
            q->tail = NULL;
        }
        q->depth--;
    }
    pthread_mutex_unlock(&q->mutex);

    if (qe != NULL)
        free(qe);
    return data;
}


/*====================================================================
* queue_destroy()
**==================================================================*/
void
queue_destroy(queue_t * q)
{
    if (q == NULL)
    {
        fprintf(stderr, "WARNING: Called queue_destroy() against a NULL pointer.\n");
        return;
    }

    pthread_mutex_destroy(&q->mutex);
    if (q->head != NULL || q->tail != NULL)
        fprintf(stderr, "WARNING: Destroying a non-empty queue.\n");
    free(q);
}

int
is_queue_empty(queue_t * q)
{
    pthread_mutex_lock(&q->mutex);
    if (q->depth != 0)
    {
        pthread_mutex_unlock(&q->mutex);
        return q->depth;
    }
    pthread_mutex_unlock(&q->mutex);
    return 0;
}

void *
peakqueue(queue_t * q)
{
    void *data;

    if (q == NULL)
    {
        fprintf(stderr, "ERROR: dequeue called on NULL queue.\n");
        exit(1);
    }

    pthread_mutex_lock(&q->mutex);
    if (q->head == NULL)
    {
        // empty queue
        data = NULL;
    }
    else
    {
        data = q->head->data;
    }
    pthread_mutex_unlock(&q->mutex);

    return data;

}


/*******************************************************************************
* END OF FILE
*******************************************************************************/
