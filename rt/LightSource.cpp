#include <cmath>
#include <cstdlib>
#include "LightSource.h"

LightSource::LightSource()
: position_()
, intensity_(10000, 10000, 10000)
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position_.x + position_.y + position_.z
                              + intensity_.x + intensity_.y + intensity_.z);
}

LightSource::LightSource(const Vec3& position, const Vec3& intensity)
: position_(position)
, intensity_(intensity)
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position.x + position.y + position.z
                              + intensity.x + intensity.y + intensity.z);
}

Ray LightSource::GenerateRay() const
{
	Vec3 dir;
    real theta = (real)M_PI * (real)erand48(xi_);
    real phi = 2.0f*(real)M_PI * (real)erand48(xi_);
    dir.x = sinf(theta) * cosf(phi);
    dir.y = sinf(theta) * sinf(phi);
    dir.z = cosf(theta);
    
	return Ray(position_, dir);
}
