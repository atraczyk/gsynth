#pragma once

#include <sstream>
#include <vector>
#include <string>
#include <cstddef>

#include <ciso646>

typedef int16_t AudioSample;
static const constexpr unsigned DEFAULT_SAMPLE_RATE = 44100;

struct AudioFormat {
    unsigned sample_rate;
    unsigned nb_channels;

    constexpr AudioFormat(unsigned sr, unsigned c) : sample_rate(sr), nb_channels(c) {}

    inline bool operator == (const AudioFormat &b) const {
        return ((b.sample_rate == sample_rate) && (b.nb_channels == nb_channels));
    }

    inline bool operator != (const AudioFormat &b) const {
        return !(*this == b);
    }

    static const AudioFormat NONE() { return AudioFormat{ 0, 0 }; }
    static const AudioFormat MONO() { return AudioFormat{ DEFAULT_SAMPLE_RATE, 1 }; }
    static const AudioFormat STEREO() { return AudioFormat{ DEFAULT_SAMPLE_RATE, 2 }; }
};

class AudioBuffer {
public:
    AudioBuffer() : AudioBuffer{ 0, AudioFormat::NONE() } {}

    AudioBuffer(size_t sample_num, AudioFormat format);

    /**
    * Copy constructor
    */
    AudioBuffer(const AudioBuffer& other, bool copy_content = false);

    /**
    * Move constructor
    */
    AudioBuffer(AudioBuffer&& other) : sampleRate_(other.sampleRate_), samples_(std::move(other.samples_)) {};

    /**
    * Copy operator
    */
    AudioBuffer& operator=(const AudioBuffer& other);

    /**
    * Move operator
    */
    AudioBuffer& operator=(AudioBuffer&& other);

    inline size_t size() const
    {
        return frames() * channels() * sizeof(AudioSample);
    }

    /**
    * Returns the sample rate (in samples/sec) associated to this buffer.
    */
    int getSampleRate() const;

    /**
    * Set the sample rate (in samples/sec) associated to this buffer.
    */
    void setSampleRate(int sr);

    /**
    * Returns the number of channels in this buffer.
    */
    inline unsigned channels() const
    {
        return static_cast<unsigned>(samples_.size());
    }

    /**
    * Set the buffer format (channels and sample rate).
    * No data conversion is performed.
    */
    void setFormat(AudioFormat format);

    inline AudioFormat getFormat() const
    {
        return AudioFormat(sampleRate_, channels());
    }

    /**
    * Returns the number of (multichannel) frames in this buffer.
    */
    inline size_t frames() const
    {
        if (not samples_.empty())
            return samples_[0].size();
        else
            return 0;
    }

    inline size_t capacity() const
    {
        return frames() * channels();
    }

    void resize(size_t sample_num);

    void clear()
    {
        for (auto& c : samples_)
            c.clear();
    }

    /**
    * Set all samples in this buffer to 0. Buffer size is not changed.
    */
    void reset()
    {
        for (auto& c : samples_)
            std::fill(c.begin(), c.end(), 0);
    }

    inline std::vector<std::vector<AudioSample> > &getData() {
        return samples_;
    }

    /**
    * Returns pointers to non-interleaved raw data.
    * Caller should not store result because pointer validity is
    * limited in time.
    */
    inline const std::vector<AudioSample*> getDataRaw()
    {
        const unsigned chans = channels();
        std::vector<AudioSample*> raw_data(chans, nullptr);
        for (unsigned i = 0; i<chans; i++)
            raw_data[i] = samples_[i].data();
        return raw_data;
    }

    /**
    * Convert fixed-point channel to float and write in the out buffer (Float 32-bits).
    * The out buffer must be at least of size capacity()*sizeof(float) bytes.
    *
    * @returns Number of samples writen.
    */
    size_t channelToFloat(float* out, const int& channel) const;

    /**
    * Write null data (silence) to the out buffer (fixed-point 16-bits).
    * The out buffer must be at least of size capacity()*sizeof(AudioSample) bytes.
    *
    * @returns Number of samples writen.
    */
    size_t fillWithZero(AudioSample* out) const;

    /**
    * convert float planar data to signed 16
    */
    void convertFloatPlanarToSigned16(uint8_t** extended_data, size_t frame_num, unsigned nb_channels = 1);


    /**
    * In-place gain transformation.
    *
    * @param gain: 0.0 -> 1.0 scale
    */
    void applyGain(double gain);

    /**
    * Mix samples from the other buffer within this buffer (in-place simple addition).
    * If the other buffer has more channels than this one, only the first this.channels() channels are imported.
    * If the other buffer has less channels than this one, behavior depends on upmix.
    * Sample rate is not considered by this function.
    *
    * TODO: some kind of check for overflow/saturation.
    *
    * @param other: the other buffer to mix in this one.
    * @param upmix: if true, upmixing occurs when other.channels() < this.channels().
    *              If false, only the first other.channels() channels are edited in this buffer.
    *
    * @returns Number of samples modified.
    */
    size_t mix(const AudioBuffer& other, bool upmix = true);

    /**
    * Copy sample_num samples from in (from sample sample pos_in) to this buffer (at sample sample pos_out).
    * If sample_num is -1 (the default), the entire in buffer is copied.
    *
    * Buffer sample number is increased if required to hold the new requested samples.
    */
    size_t copy(AudioBuffer& in, int sample_num = -1, size_t pos_in = 0, size_t pos_out = 0, bool upmix = true);

    /**
    * Copy sample_num samples from in to this buffer (at sample pos_out).
    * Input data is treated as mono and samples are duplicated in the case of a multichannel buffer.
    *
    * Buffer sample number is increased if required to hold the new requested samples.
    */
    size_t copy(AudioSample* in, size_t sample_num, size_t pos_out = 0);

private:
    int sampleRate_;

    // buffers holding data for each channels
    std::vector<std::vector<AudioSample> > samples_;
};