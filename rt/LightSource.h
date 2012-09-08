#ifndef LightSource_h
#define LightSource_h

#include "Common.h"

/**
 * point light source
 */
class LightSource
{
public:
    LightSource();
    LightSource(const Vec3& position, const Vec3& intensity);
    
    Ray GenerateRay() const;
    
	Vec3 position_;
    Vec3 intensity_;
private:
	mutable unsigned short xi_[3];
};

#endif
