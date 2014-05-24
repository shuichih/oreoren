#ifndef Light_h
#define Light_h

#include "Common.h"
#include "Scene.h"
#include "BVH.h"
#include "Sampler.h"

class Random;

//--------------------------------------------------------------------------------
enum LightType
{
    Lit_Invalid = -1,
    Lit_Ambient = 0,
    Lit_Point = 1,
    Lit_Area = 2,
    Lit_Sphere = 3
};

//--------------------------------------------------------------------------------
/**
 * Light Interface.
 */
class Light
{
public:
    
public:
    
    virtual ~Light() {}
    
    LightType GetType() const { return type_; };
    const Vec3& GetFlux() const { return flux_; };
    virtual Ray GenerateRay(Random& rand) const = 0;
    virtual Vec3 DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra, Random& rand) const = 0;
    
    
protected:
    Light()
        : type_(Lit_Invalid)
        , flux_(10000, 10000, 10000)
    {}
    
    Light(LightType type, const Vec3& flux)
        : type_(type)
        , flux_(flux)
    {}
    
    LightType type_;
    Vec3 flux_; // flux
};

//--------------------------------------------------------------------------------
/**
 * Ambient Light.
 */
class AmbientLight : public Light
{
public:
    AmbientLight(const Vec3& intensity);
    
    virtual Ray GenerateRay(Random& rand) const;
    virtual Vec3 DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra, Random& rand) const;
};

//--------------------------------------------------------------------------------
/**
 * Point Light.
 */
class PointLight : public Light
{
public:
    PointLight(const Vec3& position, const Vec3& intensity);
    
    virtual Ray GenerateRay(Random& rand) const;
    virtual Vec3 DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra, Random& rand) const;
    
	Vec3 position_;
};

//--------------------------------------------------------------------------------
/**
 * Area Light.
 */
class AreaLight : public Light
{
public:
    //AreaLight();
    AreaLight(
        const Vec3& p0,
        const Vec3& p1,
        const Vec3& p2,
        const Vec3& p3,
        const Vec3& intensity,
        int nSamples
    );
    
    virtual Ray GenerateRay(Random& rand) const;
    virtual Vec3 DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra, Random& rand) const;
    const Vec3& GetIrradiance() const { return irradiance_; };
    
	Vec3 p_[4];
    Vec3 normal_;
    float area_;
    float pdf_;
    Vec3 irradiance_;
    int nSamples_;
    Sampler sampler_;
    
private:
    float CalcTriangleArea(const Vec3& p0, const Vec3& p1, const Vec3& p2);
};

// AreaLightShape
class AreaLightShape : public Triangle
{
public:
    AreaLightShape(AreaLight* pLitSrc, const Vec3 points[3], const RGB& color, OldMaterial* pMtl);
    virtual ~AreaLightShape() {};
    
    virtual Vec3 SelfIrradiance() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    
private:
    AreaLight* pLight_;
};

//--------------------------------------------------------------------------------
class SphereLightShape;

/**
 * Sphere Light.
 */
class SphereLight : public Light
{
public:
    //SphereLight();
    SphereLight(const Vec3& position, float radius, const Vec3& intensity, int nSamples);
    
    void SetShape(SphereLightShape* pShape) { pShape_ = pShape; }
    virtual Ray GenerateRay(Random& rand) const;
    virtual Vec3 DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra, Random& rand) const;
    const Vec3& GetIrradiance() const { return irradiance_; };
    
	Vec3 position_;
    float radius_;
    Vec3 irradiance_;
    SphereLightShape* pShape_;
    int nSamples_;
};

// SphereLightShape
class SphereLightShape : public Sphere
{
public:
    SphereLightShape(SphereLight* pLitSrc, float radius, const Vec3& pos, const RGB& color, OldMaterial* pMtl);
    virtual ~SphereLightShape() {};
    
    virtual Vec3 SelfIrradiance();
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    
    SphereLight* pLight_;
};

#endif
