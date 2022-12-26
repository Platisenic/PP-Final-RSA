#ifndef TIMER_H
#define TIMER_H

#define ENABLE_TIMER
#ifdef ENABLE_TIMER
    /*
     * If ENABLE_TIMER is defined we define a clock and two macros START and
     * STOP. Use START() to start a timer and STOP("some message") to stop
     * it and print the time elapsed since START was called in ms to stdout.
     */
    clock_t timer;
    #define START() timer = std::clock();
    #define STOP(msg) \
        std::cerr << msg << " in ";\
        std::cerr << (1000.0 * (std::clock() - timer)/CLOCKS_PER_SEC);\
        std::cerr << " ms" << std::endl;
#else
    // Else the macros are defined as no-ops.
    #define START(x)
    #define STOP(x)
#endif // ENABLE_TIMER

#define ENABLE_TIMER_SUM
#ifdef ENABLE_TIMER_SUM
    clock_t timer_sum;
    int i = 0;
    double timeSum[10] ={0.0};
    #define START_SUM() timer_sum = std::clock();
    #define STOP_SUM(index) \
        timeSum[index] += (1000.0 * (std::clock() - timer_sum)/CLOCKS_PER_SEC);
    #define SHOW() \
        for(i=0; i<10; i++){if(timeSum[i]!=0.0) std::cerr << "timeSum[" << i << "]: " << timeSum[i] <<" ms" << std::endl;}
    #define CLEAN_SUM() \
        for(i=0; i<10; i++){timeSum[i]=0.0;}
#else
    // Else the macros are defined as no-ops.
    #define START_SUM(x)
    #define STOP_SUM(x)
    #define SHOW()
    #define CLEAN_SUM()
#endif // ENABLE_TIMER_SUM

#endif // TIMER_H
