#ifndef Sampler_H
#define Sampler_H

#include "Common.h"
#include "Random.h"

/**
 *
 */
class Sampler
{
public:
    Sampler();
    Vec3 SampleHemisphere(float e, const Vec3& w);
    Vec3 SampleHemisphere(float e);
    Vec3 SampleHemisphere();
    
private:
    Random rand;
};

#endif

