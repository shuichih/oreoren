#ifndef _PhotonMapRenderer_
#define _PhotonMapRenderer_

#include "Common.h"
#include "IRenderer.h"
#include "Config.h"

class Photon_map;
class PhotonFilter;
class Scene;
class Random;
struct HitRecord;

class PhotonMapRenderer : public IRenderer
{
public:
    
    //
    
    PhotonMapRenderer();
    ~PhotonMapRenderer();
    virtual void SetConfig(const Config& config);
    virtual void Run(Image& image, const Scene& scene);

private:
    
    struct PathInfo
    {
        u8 depth;
        u8 diffuseDepth;
        u8 glossyDepth;
        u8 specularDepth;
        u8 lightNo;
        PathInfo()
        : depth(0)
        , diffuseDepth(0)
        , glossyDepth(0)
        , specularDepth(0)
        , lightNo(0)
        {}
    };
    void PhotonTracing();
    void TracePhoton(const Ray& r, const Vec3& power, PathInfo pathInfo, Random& rand);
    bool Intersect(const Ray& r, HitRecord& rec);
    Vec3 CosImportanceSamplingRay(const Vec3& n);
    Vec3 GlossyRay(const Vec3& w, float exponent, Random& rand);
    void RayTracing(Image& image);
    Vec3 Irradiance(const Ray &r, PathInfo pathInfo, Random& rand);
    
    const Config* pConf_;
    const PhotonMapRendererConfig* pPmRenConf_;
    const PhotonMapConfig* pPmConf_;
    Photon_map* pPhotonMap_;
    PhotonFilter* pFilter_;
    const Scene* pScene_;
};

#endif
