#include "Timer.h"

Timer::Timer(const std::string& name) 
    : m_name {name}
{/* Timer Constructor */}

void Timer::start() {
    if ((!m_is_running) && (m_is_reset)) {
        m_is_running = true;
        m_is_reset = false;
        m_start_time = std::chrono::high_resolution_clock::now();
    }
}

void Timer::stop() {
    if (!m_is_running) {
        return;
    }
    m_stop_time = std::chrono::high_resolution_clock::now();
    m_is_running = false;
    auto a = std::chrono::high_resolution_clock::now();

}

void Timer::reset() {
    m_is_reset = true;
}

double Timer::elapsed_seconds() {
    if (m_is_running) {
        auto current_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(current_time - m_start_time).count();
    }
    else {
        return std::chrono::duration<double>(m_stop_time - m_start_time).count();
    }
}
