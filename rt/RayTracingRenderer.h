#ifndef RayTracingRenderer_h
#define RayTracingRenderer_h

#include "Common.h"
#include "IRenderer.h"
#include "Config.h"

struct HitRecord;
class Scene;
class BVH;

class RayTracingRenderer : public IRenderer
{
public:
    
    RayTracingRenderer();
    ~RayTracingRenderer();
    virtual void SetConfig(const Config& config);
    virtual void Run(Vec3* pColorBuf, const Scene& scene, BVH* pBVH);
    
private:
    
    bool Intersect(const Ray& r, HitRecord& rec);
    void RayTracing(Vec3* pColorBuf);
    Vec3 Irradiance(const Ray &r, int depth);
    
    const Config* pConfig_;
    const RayTracingConfig* pRtConfig_;
    unsigned short xi_[3];
    const Scene* pScene_;
    BVH* pBVH_;
};

#endif
