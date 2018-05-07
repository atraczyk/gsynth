#pragma once

#include <thread>
#include <atomic>
#include <cinttypes>

using AudioSample = int16_t;
static constexpr size_t maxBufLength = 4096;

template<typename T>
class lock_free_queue
{
public:
    lock_free_queue() {
        write_.store(0);
        read_.store(0);
    };

    size_t getBufSize() {
        return write_ - read_;
    };

    bool try_push(T val);
    void push(T val);

    bool try_pop(T& element);
    T pop();

private:
    std::atomic<size_t> write_;
    std::atomic<size_t> read_;

    size_t capacity_ = maxBufLength;
    T buffer_[maxBufLength];
};

template<typename T>
bool lock_free_queue<T>::try_push(T val)
{
    const auto current_tail = write_.load();
    const auto next_tail = 1 + current_tail;

    if (current_tail - read_.load() <= capacity_ - 1) {
        buffer_[current_tail % capacity_] = val;
        write_.store(next_tail);
        return true;
    }

    return false;
}

template<typename T>
void lock_free_queue<T>::push(T val)
{
    while (!tryPush(val));
}

template<typename T>
bool lock_free_queue<T>::try_pop(T& element)
{
    const auto currentHead = read_.load();

    if (currentHead != write_.load()) {
        element = buffer_[currentHead % capacity_];
        read_.store(1 + currentHead);
        return true;
    }

    return false;
}

template<typename T>
T lock_free_queue<T>::pop()
{
    T ret;
    while (!tryPop(&ret));
    return ret;
}

using RingBuffer = lock_free_queue<AudioSample>;