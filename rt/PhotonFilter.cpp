#include <cmath>
#include "PhotonFilter.h"

//

ConeFilter::ConeFilter(float k)
{
    k_ = k;
    normalizer_ = 1.f / (1.f - (2.f / 3.f*k_));
}

float ConeFilter::Weight(float dist, float maxDist) const
{
    return 1.f - dist / (k_ * maxDist);
}

float ConeFilter::Normalizer() const
{
    return normalizer_;
}
