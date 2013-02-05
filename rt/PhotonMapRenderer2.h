#ifndef _PhotonMapRenderer2_
#define _PhotonMapRenderer2_

#include "Common.h"
#include "IRenderer.h"
#include "Config.h"

class Photon_map;
class PhotonFilter;
struct HitRecord;
class Scene;

class PhotonMapRenderer2 : public IRenderer
{
public:
    
    //
    
    PhotonMapRenderer2();
    ~PhotonMapRenderer2();
    virtual void SetConfig(const Config& config);
    virtual void Run(Vec3* pColorBuf, const Scene& scene);

private:
    
    enum TraceFlag {
        Trace_Indirect = 0,
        Trace_Coarstic = 1,
        Trace_Shadow   = 2,
    };
    
    struct PathInfo
    {
        int depth;
        int diffuseDepth;
        int glossyDepth;
        int specularDepth;
        int refractionDepth;
        
        PathInfo()
        : depth(0)
        , diffuseDepth(0)
        , glossyDepth(0)
        , specularDepth(0)
        , refractionDepth(0)
        {}
    };
    
    //typedef void (PhotonMapRenderer2::*TracePhotonFunc)(const Ray& r, const Vec3& power, PathInfo& pathInfo);
    
    void InitializePhotonMap(Photon_map** ppPm, const PhotonMapConfig& pmConf, PhotonFilter** ppFilter);
    void PhotonTracing();
    void PhotonTracing_(Photon_map& photonMap, const PhotonMapConfig& pmConfig, int nLit,
                        double sumIntensity, TraceFlag traceFlag);
    void TracePhoton(const Ray& r, const Vec3& power, PathInfo& pathInfo);
    bool Intersect(const Ray& r, HitRecord& rec);
    Vec3 CosImportanceSamplingRay(const Vec3& n);
    Vec3 GlossyRay(const Vec3& w, float exponent);
    void RayTracing(Vec3* pColorBuf);
    Vec3 Irradiance(const Ray &r, PathInfo& pathInfo);
    
    unsigned short xi_[3];
    const Config* pConfig_;
    const Scene* pScene_;
    Photon_map* pCurrPm_;
    Photon_map* pPhotonMap_; // indirectPhotonMap
    Photon_map* pCoarsticPhotonMap_;
    Photon_map* pShadowPhotonMap_;
    const PhotonMapConfig* pPmConfig_; // indirectPmConfig
    const PhotonMapConfig* pCoarsticPmConfig_;
    const PhotonMapConfig* pShadowPmConfig_;
    PhotonFilter* pFilter_;
    PhotonFilter* pCoarsticFilter_;
    std::vector<HitRecord> shadyHits_;
    
    TraceFlag traceFlag_;
    int lightNo_;
};

#endif
