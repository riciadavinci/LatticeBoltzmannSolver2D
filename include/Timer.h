#pragma once

#include <chrono>
#include <string>

// TODO: Add custom exceptions

class Timer
{
public:
    Timer(const std::string& name);
    void start();
    void stop();
    void reset();
    double elapsed_seconds();

private:
    std::string m_name;
    std::chrono::steady_clock::time_point m_start_time;
    std::chrono::steady_clock::time_point m_stop_time;
    bool m_is_running {false};
    bool m_is_reset {true};
};
