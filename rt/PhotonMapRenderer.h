#ifndef _PhotonMapRenderer_
#define _PhotonMapRenderer_

#include "Common.h"
#include "IRenderer.h"
#include "Config.h"

class Photon_map;
class PhotonFilter;
struct HitRecord;
class Scene;
class BVH;

class PhotonMapRenderer : public IRenderer
{
public:
    
    //
    
    PhotonMapRenderer();
    ~PhotonMapRenderer();
    virtual void SetConfig(const Config& config);
    virtual void Run(Vec3* pColorBuf, const Scene& scene, BVH* pBVH);

private:
    
    void PhotonTracing(const Ray& r, float power[3], int depth);
    bool Intersect(const Ray& r, HitRecord& rec);
    void RayTracing(Vec3* pColorBuf);
    Vec3 Irradiance(const Ray &r, int depth);
    
    const Config* pConfig_;
    const PhotonMapConfig* pPmConfig_;
    Photon_map* pPhotonMap_;
    unsigned short xi_[3];
    PhotonFilter* pFilter_;
    const Scene* pScene_;
    BVH* pBVH_;
};

#endif
