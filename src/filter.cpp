#include "filter.h"

#include "common.h"

LowPassFilter::LowPassFilter(unsigned sampleRate, double Fc, double Q, double peakGainDB)
{
    z1_ = z2_ = 0.0;
    sampleRate_ = sampleRate;
    Q_ = Q / sampleRate_;
    Fc_ = Fc;
    peakGain_ = peakGainDB;
    configureCoefficients();
}

LowPassFilter::~LowPassFilter()
{
}

void LowPassFilter::setQ(double Q)
{
    Q_ = Q / sampleRate_;
    configureCoefficients();
}

void LowPassFilter::setFc(double Fc)
{
    Fc_ = Fc;
    configureCoefficients();
}

double LowPassFilter::processSample(double in)
{
    double out = in * a0_ + z1_;
    z1_ = in * a1_ + z2_ - b1_ * out;
    z2_ = in * a2_ - b2_ * out;
    return out;
}

void LowPassFilter::configureCoefficients()
{
    double norm;
    double V = pow(10, fabs(peakGain_) / 20.0);
    double K = tan(M_PI * Fc_);

    norm = 1 / (1 + K / Q_ + K * K);
    a0_ = K * K * norm;
    a1_ = 2 * a0_;
    a2_ = a0_;
    b1_ = 2 * (K * K - 1) * norm;
    b2_ = (1 - K / Q_ + K * K) * norm;
}