#include "ringbuffer.h"

bool RingBuffer::tryPush(AudioSample val)
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

void RingBuffer::push(AudioSample val)
{
    while (!tryPush(val));
}

bool RingBuffer::tryPop(AudioSample* element)
{
    const auto currentHead = read_.load();

    if (currentHead != write_.load()) {
        *element = buffer_[currentHead % capacity_];
        read_.store(1 + currentHead);
        return true;
    }

    return false;
}

AudioSample RingBuffer::pop()
{
    AudioSample ret;
    while (!tryPop(&ret));
    return ret;
}