#include <cmath>
#include <cstdlib>
#include "LightSource.h"

LightSource::LightSource(Vec position, float intensity)
: position_(position)
, intensity_(intensity)
{
    intensity_ = intensity;
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position.x + position.y + position.z + intensity);
}

Ray LightSource::GenerateRay()
{
	Vec dir;
    float theta = M_PI * erand48(xi_);
    float phi = 2.0f*M_PI * erand48(xi_);
    dir.x = sinf(theta) * cosf(phi);
    dir.y = sinf(theta) * sinf(phi);
    dir.z = cosf(theta);
    
	return Ray(position_, dir);
}
