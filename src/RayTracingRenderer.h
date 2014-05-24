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
    virtual void Run(Image& image, const Scene& scene);
    
private:
    
    void RayTracing(Image& image);
    Vec3 Irradiance(const Ray &r, int depth, Random& rand);
    
    const Config* pConfig_;
    const RayTracingConfig* pRtConfig_;
    const Scene* pScene_;
    BVH* pBVH_;
};

#endif
