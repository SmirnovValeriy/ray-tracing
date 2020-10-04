#include "Time.hpp"


Timer::Timer(): bench_t_end(0), bench_t_start(0) {}

double get_time() {
    struct timeval Tp;
    /* structure including fields:
       long tv_sec - seconds since Jan. 1, 1970
       long tv_usec - microseconds */
    int stat;
    stat = gettimeofday(&Tp, nullptr); /* get current time and data in seconds and microseconds (fill fields of struct Tp) */
    if (stat != 0) /* 0 - success, (-1) - failure */
        printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void Timer::start() {
    bench_t_start = get_time(); /* remember time of start */
}

void Timer::stop() {
    bench_t_end = get_time(); /* remember time of end */
}

void Timer::print_time() const {
    printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start); /* print runtime */
}