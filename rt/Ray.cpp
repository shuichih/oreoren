#include "Ray.h"
#include "stdlib.h"
#include "Random.h"

Ray::Ray(Vec3 o_, Vec3 d_) : o(o_), d(d_)
{
#ifdef RAY_USE_SIMD
    xmo = _mm_set_ps(o.x, o.y, o.z, 0);
#endif
    SetDirection(d_);
}

void Ray::SetDirection(const Vec3& dir)
{
    d = dir;
    idir = Vec3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
#ifdef RAY_USE_SIMD
    xmInvDir = _mm_set_ps(idir.x, idir.y, idir.z, 0);
#endif
    
#ifdef BBOX_USE_SIMD
    xmDirSignMask = _mm_setzero_ps();
#endif
    if (dir.x > 0) {
        sign[0] = 0;
    } else {
        sign[0] = 1;
#ifdef BBOX_USE_SIMD
        xmDirSignMask_i = _mm_add_epi32(xmDirSignMask_i, _mm_set_epi32(0xffffffff, 0, 0, 0));
#endif
    }
    if (dir.y > 0) {
        sign[1] = 0;
    } else {
        sign[1] = 1;
#ifdef BBOX_USE_SIMD
        xmDirSignMask_i = _mm_add_epi32(xmDirSignMask_i, _mm_set_epi32(0, 0xffffffff, 0, 0));
#endif
    }
    if (dir.z > 0) {
        sign[2] = 0;
    } else {
        sign[2] = 1;
#ifdef BBOX_USE_SIMD
        xmDirSignMask_i = _mm_add_epi32(xmDirSignMask_i, _mm_set_epi32(0, 0, 0xffffffff, 0));
#endif
    }
}

Vec3 Ray::PointAtParameter(float t) const
{
    return o + t * d;
}
    
// generate cos distribution ray
Vec3 Ray::CosRay(const Vec3& w, Random& rand)
{
    real r1 = 2.f*(PI*rand.F32());
    real r2 = rand.F32(); // => 1-cos^2θ = 1-sqrt(1-r_2)^2 = r_2
    real r2s = sqrtf(r2);    // => sinθ = sqrt(1-cos^2θ) = sqrt(r_2)
    Vec3 u = ((fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) % w).normalize(); // binormal
    Vec3 v = w % u; // tangent
    
    // ucosφsinθ + vsinφsinθ + wcosθ
    return (u*cosf(r1)*r2s + v*sinf(r1)*r2s + w*sqrtf(1-r2)).normalize();
}

