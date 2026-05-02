#include <stdio.h>
#include <time.h>
//usage
//clock_t start = clock() ;
//do_some_work() ;
//clock_t end = clock() ;
//double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
void C_wraptime(double *res) {

   *res=(double) clock()/(double)CLOCKS_PER_SEC;
    printf("time is %d\n:",*res);
}

