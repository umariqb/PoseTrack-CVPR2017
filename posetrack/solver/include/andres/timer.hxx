#pragma once
#ifndef ANDRES_TIMER_HXX
#define ANDRES_TIMER_HXX

#include <chrono>

namespace andres {

template<class T = double>
class Timer {
public:
    typedef T value_type;

    Timer();
    void reset();
    void start();
    void stop();
    value_type elapsedSeconds() const;

private:
    value_type seconds_;

    decltype(std::chrono::high_resolution_clock::now()) time_;
};

template<class T>
inline
Timer<T>::Timer()
:   seconds_()
{}

template<class T>
inline void
Timer<T>::reset() {
    seconds_ = .0;
}

template<class T>
inline void
Timer<T>::start() {
    time_ = std::chrono::high_resolution_clock::now();
}

template<class T>
inline void
Timer<T>::stop() {
    seconds_ += std::chrono::duration_cast<std::chrono::duration<value_type>>(
        std::chrono::high_resolution_clock::now() - time_
    ).count();
}

template<class T>
inline typename Timer<T>::value_type
Timer<T>::elapsedSeconds() const {
    return seconds_;
}

} // namespace andres

#endif
