#ifndef Timer_h
#define Timer_h

#include <ctime>
#include <string>

//#define TIMER_USE_CLOCK

/**
 * Timerクラス
 */
class Timer
{
public:
    Timer();
    Timer(const char* pMsg);
    ~Timer();
    
    void Restart();
    long Elapsed();
    long Stop();
    void PrintElapsed(const char* str=NULL);
    
private:
    std::string msg_;
    bool ticking_;
#ifdef TIMER_USE_CLOCK
    clock_t start_;
    clock_t end_;
#else
    time_t start_;
    time_t end_;
#endif
};

#endif
