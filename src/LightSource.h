#ifndef LightSource_h
#define LightSource_h

#include "Common.h"
#include "Scene.h"
#include "BVH.h"
#include "Sampler.h"

class Random;

//--------------------------------------------------------------------------------
enum LightSourceType
{
    Lit_Invalid = -1,
    Lit_Point = 0,
    Lit_Area = 1,
    Lit_Sphere = 2
};

//--------------------------------------------------------------------------------
/**
 * Light Source Interface.
 */
class LightSource
{
public:
    
public:
    
    virtual ~LightSource() {}
    
    LightSourceType GetType() const { return type_; };
    const Vec3& GetFlux() const { return flux_; };
    virtual Ray GenerateRay(Random& rand) const = 0;
    virtual Vec3 DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra, Random& rand) const = 0;
    
    
protected:
    LightSource()
        : type_(Lit_Invalid)
        , flux_(10000, 10000, 10000)
    {}
    
    LightSource(LightSourceType type, const Vec3& flux)
        : type_(type)
        , flux_(flux)
    {}
    
    LightSourceType type_;
    Vec3 flux_; // flux
};

/**
 * Point Light Source.
 */
class PointLightSource : public LightSource
{
public:
    PointLightSource();
    PointLightSource(const Vec3& position, const Vec3& intensity);
    
    virtual Ray GenerateRay(Random& rand) const;
    virtual Vec3 DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra, Random& rand) const;
    
	Vec3 position_;
};

//--------------------------------------------------------------------------------
/**
 * Area Light Source.
 */
class AreaLightSource : public LightSource
{
public:
    //AreaLightSource();
    AreaLightSource(
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
    AreaLightShape(AreaLightSource* pLitSrc, const Vec3 points[3], const RGB& color, OldMaterial* pMtl);
    virtual ~AreaLightShape() {};
    
    virtual Vec3 SelfIrradiance() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    
private:
    AreaLightSource* pLightSource_;
};

//--------------------------------------------------------------------------------
class SphereLightShape;

/**
 * Sphere Light Source.
 */
class SphereLightSource : public LightSource
{
public:
    //SphereLightSource();
    SphereLightSource(const Vec3& position, float radius, const Vec3& intensity, int nSamples);
    
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
    SphereLightShape(SphereLightSource* pLitSrc, float radius, const Vec3& pos, const RGB& color, OldMaterial* pMtl);
    virtual ~SphereLightShape() {};
    
    virtual Vec3 SelfIrradiance();
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    
    SphereLightSource* pLightSource_;
};

#endif
