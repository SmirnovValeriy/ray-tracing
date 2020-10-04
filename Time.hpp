#pragma once
#include <ctime>
#include <sys/time.h>
#include <cstdio>


double get_time();

class Timer {
    double bench_t_start, bench_t_end;

public:
    Timer();
    void start();
    void stop();
    void print_time() const;
};