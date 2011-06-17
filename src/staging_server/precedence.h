#ifndef __PRECEDENCE_H__
#define __PRECEDENCE_H__

struct precedence_struct;

struct precedence_struct * precedence_create (int max_prec);
void precedence_get (struct precedence_struct *p, int my_prec);
void precedence_release (struct precedence_struct *p);

#endif
