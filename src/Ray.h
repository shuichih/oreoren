#ifndef Ray_H
#define Ray_H

#include "Common.h"
#include "simd.h"

class Random;

class Ray
{
public:
    Vec3 o;
    Vec3 d;
    Vec3 idir;
    i32 sign[3]; // sign of the direction, 0:positive 1:negative
    i32 pad;
#ifdef RAY_USE_SIMD
    __m128 xmo;
    __m128 xmInvDir;
#endif
    
#ifdef BBOX_USE_SIMD
    union {
        __m128 xmDirSignMask; // sign of the direction 0:positive FF:negative
        __m128i xmDirSignMask_i;
    };
#endif
    
    Ray(Vec3 o_, Vec3 d_);
    Ray(const Ray& rhs);
    void SetDirection(const Vec3& dir);
    Vec3 PointAtParameter(float t) const;
    // generate cos distribution ray
    static Vec3 CosRay(const Vec3& w, Random& rand);
};

#endif // Ray

