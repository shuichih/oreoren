#ifndef LightSource_h
#define LightSource_h

#include "Common.h"

/**
 * point light source
 */
class LightSource
{
public:
    LightSource(Vec position, float intensity);
    
    Ray GenerateRay();
    
	Vec position_;
    float intensity_;
	unsigned short xi_[3];
};

#endif
