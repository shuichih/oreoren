#include <cstdio>
#include "Timer.h"
#include <ctime>


Timer::Timer()
{
    Restart();
}

void Timer::Restart()
{
#ifdef TIMER_USE_CLOCK
    start_ = clock();
#else
    time(&start_);
#endif
    ticking_ = true;
}

long Timer::Stop()
{
#ifdef TIMER_USE_CLOCK
    end_ = clock();
#else
    time(&end_);
#endif
    ticking_ = false;
    
    return Elapsed();
}

long Timer::Elapsed()
{
    long msec;
    if (ticking_) {
#ifdef TIMER_USE_CLOCK
        msec = (clock() - start_) / (CLOCKS_PER_SEC / 1000);
#else
        time_t now;
        time(&now);
        msec = difftime(now, start_) * 1000;
#endif
    } else {
#ifdef TIMER_USE_CLOCK
        msec = (end_ - start_) / (CLOCKS_PER_SEC / 1000);
#else
        msec = difftime(end_, start_) * 1000;
#endif
    }
    
    return msec;
}

void Timer::PrintElapsed(const char* str)
{
    long elapsed = Elapsed();
    long min = elapsed / 1000 / 60;
    long sec = (elapsed / 1000) % 60;
    long ms = elapsed % 1000;
    printf("%s%ldm %lds %ldms\n", str, min, sec, ms);
}
