#pragma once

class LowPassFilter
{
public:
    LowPassFilter(unsigned sampleRate, double Fc, double Q, double peakGainDB);
    ~LowPassFilter();

    void setFc(double Fc);
    void setQ(double Q);
    double processSample(double in);

private:
    double a0_, a1_, a2_, b1_, b2_;
    double Fc_, Q_, peakGain_;
    double z1_, z2_;
    unsigned sampleRate_;

    void configureCoefficients();
};
