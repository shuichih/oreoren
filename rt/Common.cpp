#include "Common.h"
#include <cstdlib>

Vec3 operator*(const float f, const Vec3& v)
{
    return Vec3(v.x * f, v.y * f, v.z * f);
}

// generate cos distribution ray
Vec3 Ray::CosRay(const Vec3& w, unsigned short seed[3])
{
    real r1 = 2.f*(real)(M_PI*erand48(seed));
    real r2 = (real)erand48(seed); // => 1-cos^2θ = 1-sqrt(1-r_2)^2 = r_2
    real r2s = sqrtf(r2);    // => sinθ = sqrt(1-cos^2θ) = sqrt(r_2)
    Vec3 u = ((fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) % w).normalize(); // binormal
    Vec3 v = w % u; // tangent
    
    // ucosφsinθ + vsinφsinθ + wcosθ
    return (u*cosf(r1)*r2s + v*sinf(r1)*r2s + w*sqrtf(1-r2)).normalize();
}
