#include "time_utils.h"

void timespec_add(struct timespec *t_a, struct timespec *t_b)
{
   struct timespec temp;
   temp.tv_nsec = t_a->tv_nsec + t_b->tv_nsec;
   temp.tv_sec = t_a->tv_sec + t_b->tv_sec;

   while(temp.tv_nsec > 1e9)
   {
       temp.tv_sec++;
       temp.tv_nsec = temp.tv_nsec - 1e9;
   }

   t_a->tv_nsec = temp.tv_nsec;
   t_a->tv_sec = temp.tv_sec;
}

void timespec_add_us(struct timespec *t, long us)
{
    t->tv_nsec += us*1000;
    if (t->tv_nsec > 1000000000) {
        t->tv_nsec = t->tv_nsec - 1000000000;// + ms*1000000;
        t->tv_sec += 1;
    }
}

int timespec_cmp(struct timespec *a, struct timespec *b)
{
    if (a->tv_sec > b->tv_sec) return 1;
    else if (a->tv_sec < b->tv_sec) return -1;
    else if (a->tv_sec == b->tv_sec) {
        if (a->tv_nsec > b->tv_nsec) return 1;
        else if (a->tv_nsec == b->tv_nsec) return 0;
        else return -1;
    }
}

int timespec_sub(struct timespec *d, struct timespec *a, struct timespec *b)
{
    d->tv_nsec = a->tv_nsec - b->tv_nsec;
    d->tv_sec =  a->tv_sec - b->tv_sec;
    if (a->tv_nsec < b->tv_nsec) {
        d->tv_nsec += 1000000000;
        d->tv_sec -= 1;
    }
    return 1;
}
