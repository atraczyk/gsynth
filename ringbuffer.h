#pragma once

#include <thread>
#include <atomic>
#include <cinttypes>

using namespace std;

static constexpr size_t SIZEBUF = 16000;

template<typename T>
class ring_buffer
{
public:

    ring_buffer() {
        write.store(0);
        read.store(0);
    }

    bool tryPush(T val) {
        const auto current_tail = m_write.load();
        const auto next_tail = increment(current_tail);

        if (current_tail - m_read.load() <= m_capacity - 1) {
            m_buffer.get()[current_tail % m_capacity] = val;
            m_write.store(next_tail);
            return true;
        }

        return false;
    }

    void push(T val) {
        while (!try_push(val));
    }

    bool tryPop(T* element) {
        const auto currentHead = m_read.load();

        if (currentHead != m_write.load()) {
            *element = m_buffer.get()[currentHead % m_capacity];
            m_read.store(increment(currentHead));
            return true;
        }

        return false;
    }

    T pop() {
        T ret;
        while (!try_pop(&ret));
        return ret;
    }

private:
    std::atomic<T> write;
    std::atomic<T> read;

    size_t size = SIZEBUF;
    T buffer[SIZEBUF];

    size_t increment(int n) {
        return (n + 1);
    }
};

using RingBuffer = ring_buffer<AudioSample>;