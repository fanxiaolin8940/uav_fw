#ifndef __TIMEUTILS_H__
#define __TIMEUTILS_H__

#include <time.h>

void timespec_add(struct timespec *ta, struct timespec *tb);
void timespec_add_us(struct timespec *t, long us);
int timespec_cmp(struct timespec *a, struct timespec *b);
int timespec_sub(struct timespec *d, struct timespec *a, struct timespec *b);

#endif 
