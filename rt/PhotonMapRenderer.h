#ifndef _PhotonMapRenderer_
#define _PhotonMapRenderer_

#include "Common.h"
#include "IRenderer.h"
#include "Config.h"

class Photon_map;
class PhotonFilter;
struct HitRecord;
class Scene;

class PhotonMapRenderer : public IRenderer
{
public:
    
    //
    
    PhotonMapRenderer();
    ~PhotonMapRenderer();
    virtual void SetConfig(const Config& config);
    virtual void Run(Vec3* pColorBuf, const Scene& scene);

private:
    
    struct PathInfo
    {
        int depth;
        int diffuseDepth;
        int glossyDepth;
        int specularDepth;
        int lightNo;
        PathInfo()
        : depth(0)
        , diffuseDepth(0)
        , glossyDepth(0)
        , specularDepth(0)
        , lightNo(0)
        {}
    };
    void PhotonTracing();
    void TracePhoton(const Ray& r, const Vec3& power, PathInfo& pathInfo);
    bool Intersect(const Ray& r, HitRecord& rec);
    Vec3 CosImportanceSamplingRay(const Vec3& n);
    Vec3 GlossyRay(const Vec3& w, float exponent);
    void RayTracing(Vec3* pColorBuf);
    Vec3 Irradiance(const Ray &r, PathInfo& pathInfo);
    
    const Config* pConfig_;
    const PhotonMapConfig* pPmConfig_;
    Photon_map* pPhotonMap_;
    unsigned short xi_[3];
    PhotonFilter* pFilter_;
    const Scene* pScene_;
};

#endif
