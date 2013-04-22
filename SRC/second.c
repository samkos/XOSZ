#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <time.h>
 
double second()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

double second_()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

double second_cpu() {
  clock_t ticks;
  ticks = clock();
 
  
  return (double)ticks/CLOCKS_PER_SEC;
}
