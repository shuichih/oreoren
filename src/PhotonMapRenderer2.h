#ifndef _PhotonMapRenderer2_
#define _PhotonMapRenderer2_

#include "Common.h"
#include "IRenderer.h"
#include "Config.h"

class Photon_map;
class PhotonFilter;
struct HitRecord;
class Scene;
class Random;
class Sampler;
class Image;

class PhotonMapRenderer2 : public IRenderer
{
public:
    
    //
    
    PhotonMapRenderer2();
    ~PhotonMapRenderer2();
    virtual void SetConfig(const Config& config);
    virtual void Run(Image& image, const Scene& scene);

private:
    
    enum TraceFlag {
        Trace_Indirect = 0,
        Trace_Caustic = 1,
        Trace_Shadow   = 2,
    };
    
    struct PathInfo
    {
        u8 depth;
        u8 diffuseDepth;
        u8 glossyDepth;
        u8 specularDepth;
        u8 refractionDepth;
        
        PathInfo()
        : depth(0)
        , diffuseDepth(0)
        , glossyDepth(0)
        , specularDepth(0)
        , refractionDepth(0)
        {}
    };
    
    void InitializePhotonMap(Photon_map** ppPm, const PhotonMapConfig& pmConf, PhotonFilter** ppFilter);
    void PhotonTracing();
    void PhotonTracing_(Photon_map& photonMap, const PhotonMapConfig& pmConfig,
                        double sumIntensity, TraceFlag traceFlag);
    void TracePhoton(const Ray& r, const Vec3& power, PathInfo pathInfo, Random& rand);
    Vec3 CosImportanceSamplingRay(const Vec3& n);
    Vec3 GlossyRay(const Vec3& w, float exponent);
    void RayTracing(Image& image);
    Vec3 Irradiance(const Ray &r, PathInfo pathInfo, Random& rand);
    
    const Config* pConf_;
    const PhotonMapRendererConfig* pPmRenConf_;
    const PhotonMapConfig* pPmConf_; // indirectPmConfig
    const PhotonMapConfig* pCausticPmConf_;
    const PhotonMapConfig* pShadowPmConf_;
    const Scene* pScene_;
    Photon_map* pCurrPm_;
    Photon_map* pPhotonMap_; // indirectPhotonMap
    Photon_map* pCausticPhotonMap_;
    Photon_map* pShadowPhotonMap_;
    PhotonFilter* pFilter_;
    PhotonFilter* pCausticFilter_;
    
    TraceFlag traceFlag_;
    int lightNo_;
    Sampler* pSampler_;
};

#endif
