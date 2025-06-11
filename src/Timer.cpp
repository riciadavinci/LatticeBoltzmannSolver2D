#include "Timer.h"

Timer::Timer(const std::string& name) 
{/* Timer Constructor */}

void Timer::start() {
    if ((!m_is_running) && (m_is_reset)) {
        m_is_running = true;
        m_is_reset = false;
        m_start_time = std::chrono::system_clock::now();
    }
}

void Timer::stop() {
    if (!m_is_running) {
        return;
    }
    m_stop_time = std::chrono::system_clock::now();
    m_is_running = false;
}

void Timer::reset() {
    m_is_reset = true;
}


std::chrono::seconds Timer::elapsed_seconds() {
    if (m_is_running) {
        auto current_time = std::chrono::system_clock::now();
        return std::chrono::duration_cast<std::chrono::seconds>(current_time - m_start_time);
    } 
    else {
        return std::chrono::duration_cast<std::chrono::seconds>(m_stop_time - m_start_time);
    }
}
