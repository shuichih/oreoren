// from Kevin Suffern, Ray Tracing from the Ground Up, 13.9, 26.2

#ifndef BRDF_H
#define BRDF_H

#include "Common.h"

struct HitRecord;
class Sampler;

//--------------------------------------------------------------------------------
/**
 * BRDF interface
 */
class IBRDF
{
public:
    virtual ~IBRDF();
    
    /**
     * A BRDF itself.
     * This will not be defined if the BRDF contains a delta function. (e.g. perfect specular)
     * @param hr a HitRecord that has a shading(hit) point.
     * @param wi the direction of incoming light / outgoing ray.
     * @param wo the direction of refrected light / incoming ray.
     */
    virtual RGB F(const HitRecord& hr, const Vec3& wi, const Vec3& wo) const = 0;
    
    /**
     * Sample A BRDF to compute light transport
     *  and compute reflected ray direction.
     * @param hr a HitRecord that has a shading(hit) point.
     * @param wi sampled direction
     * @param wo light direction
     */
    virtual RGB SampleF(const HitRecord& hr, Vec3& wi, const Vec3& wo, float& pdf) const = 0;
    
    /**
     * Bihemispherical reflectance
     *
     * perfectly diffuse surfaceのambient illuminationをモデル化するためには
     * bihemispherical reflectance: ρhhを用いる。これは入射放射輝度が等方で位置非依存なときの
     * 半球全域に反射するfluxと、半球全域からの合計入射fluxの比である。
     * (Ray Tracing from the Ground Up, 13.8)
     */
    virtual RGB Rho(const HitRecord& hr, const Vec3& wo) const = 0;
};

//--------------------------------------------------------------------------------
/**
 * Lambertian BRDF
 */
class Lambertian : public IBRDF
{
public:
    Lambertian();
    Lambertian(float kd, const RGB& cd);
    void SetKd(float k) { kd_ = k; }
    void SetCr(const RGB& c) { cd_ = c; }
    virtual RGB F(const HitRecord& hr, const Vec3& wi, const Vec3& wo) const;
    virtual RGB SampleF(const HitRecord& hr, Vec3& wi, const Vec3& wo, float& pdf) const;
    virtual RGB Rho(const HitRecord& hr, const Vec3& wo) const;
    
private:
    float kd_; // diffuse reflection coefficient [0, 1]
    RGB cd_;   // diffuse color
    Sampler* pSampler_;
};

//--------------------------------------------------------------------------------
/**
 * Perfect Specular BRDF
 */
class PerfectSpecular : public IBRDF
{
public:
    PerfectSpecular();
    PerfectSpecular(float kr, const RGB& cr);
    void SetKr(float k) { kr_ = k; }
    void SetCr(const RGB& c) { cr_ = c; }
    virtual RGB F(const HitRecord& hr, const Vec3& wi, const Vec3& wo) const;
    virtual RGB SampleF(const HitRecord& hr, Vec3& wi, const Vec3& wo, float& pdf) const;
    virtual RGB Rho(const HitRecord& hr, const Vec3& wo) const;
    
private:
    float kr_; // specular reflection coefficient [0, 1]
    RGB cr_;   // specular color
};

//--------------------------------------------------------------------------------
/**
 * Glossy Specular BRDF
 */
class GlossySpecular : public IBRDF
{
public:
    GlossySpecular();
    GlossySpecular(float kr, const RGB& cr, float exp);
    virtual RGB F(const HitRecord& hr, const Vec3& wi, const Vec3& wo) const;
    virtual RGB SampleF(const HitRecord& hr, Vec3& wi, const Vec3& wo, float& pdf) const;
    virtual RGB Rho(const HitRecord& hr, const Vec3& wo) const;

private:
    float ks_; // specular reflection coefficient [0, 1]
    RGB cs_;   // specular color
    float exp_;
    Sampler* pSampler_;
};

#endif

