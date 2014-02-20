#include "Sampler.h"
#include "Scene.h"


//--------------------------------------------------------------------------------

// Ray::CosRay, GlossyRay()をSamplerを使って置き換え
Sampler::Sampler()
{
}

Vec3 Sampler::SampleHemisphere(float e, const Vec3& w)
{
    Vec3 u = (fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) ^ w;
    u.normalize(); // binormal
    Vec3 v = w ^ u; // tangent
    
    Vec3 s = SampleHemisphere(e);
    return u * s.x + v * s.y + w * s.z;
}

/**
 * 指数eに従ったcosine分布で半球上の点をサンプリング
 * e==0だと班球上で一様な分布になる
 * UVW座標系でW軸が半球の中央
 *
 * Shirley, Realistic Ray Tracing, N
 * Suffern, Ray Tracing from the Ground Up, 7.3
 */
Vec3 Sampler::SampleHemisphere(const float e)
{
    real r1 = 2.f * PI * rand.F32(); // phi(azimuth angle)
    real r2 = rand.F32(); // theta(polar angle)
    
    // theta = cos^-1(powf(1 - r2, 1 / (e + 1))) らしい
    // uniform分布のとき theta = cos^-1(1 - r2) で、
    // 半球の外側にサンプルが多く集まるようにcos^1が入ってる感じ
    real cosTheta = powf(1.f - r2, 1.0f / (e + 1.0f));
    real sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
    real rx = cosf(r1) * sinTheta;
    real ry = sinf(r1) * sinTheta;
    real rz = cosTheta;
    return Vec3(rx, ry, rz);
}

/**
 * 指数eに従った一様分布で半球上の点をサンプリング
 * UVW座標系でW軸が半球の中央
 *
 * Shirley, Realistic Ray Tracing, N
 * Suffern, Ray Tracing from the Ground Up, 7.3
 */
Vec3 Sampler::SampleHemisphere()
{
    real r1 = 2.f * PI * rand.F32(); // phi
    real r2 = rand.F32(); // theta
    real cosTheta = 1.f - r2;
    real sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
    real rx = cosf(r1) * sinTheta;
    real ry = sinf(r1) * sinTheta;
    real rz = cosTheta;
    return Vec3(rx, ry, rz);
}
