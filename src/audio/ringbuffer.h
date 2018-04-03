#pragma once

#include <thread>
#include <atomic>
#include <cinttypes>

using namespace std;

using AudioSample = int16_t;
static constexpr size_t SIZEBUF = 16000;

class RingBuffer
{
public:
    RingBuffer() {
        write_.store(0);
        read_.store(0);
    }

    size_t getBufSize() {
        return write_ - read_;
    }

    bool tryPush(AudioSample val);
    void push(AudioSample val);

    bool tryPop(AudioSample* element);
    AudioSample pop();

private:
    std::atomic<size_t> write_;
    std::atomic<size_t> read_;

    size_t capacity_ = SIZEBUF;
    AudioSample buffer_[SIZEBUF];
};