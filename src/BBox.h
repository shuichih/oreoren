#ifndef BBox_h
#define BBox_h

#include "Common.h"
#include "simd.h"

class Ray;

/**
 * Axis-aligned bounding box.
 */
class BBox
{
public:
    static BBox Surround(const BBox& b1, const BBox& b2);
    
public:
    BBox() {};
    BBox(const Vec3& min, const Vec3& max);
    BBox(float minX, float minY, float minZ, float maxX, float maxY, float maxZ);
    
    Vec3& Min() { return pp[0]; }
    Vec3& Max() { return pp[1]; }
    const Vec3& Min() const { return pp[0]; }
    const Vec3& Max() const { return pp[1]; }
    Vec3 Center() const { return (pp[0] + pp[1]) * 0.5f; };
    Vec3 Size() const { return pp[1] - pp[0]; };
    virtual bool RayIntersect(const Ray &r, float tmin, float tmax) const;
    
    Vec3 pp[2]; // [0]:min [1]:max
#ifdef BBOX_USE_SIMD
    __m128 xmpp0;
    __m128 xmpp1;
#endif
};

#endif
