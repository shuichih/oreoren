#ifndef _PhotonMapRenderer_
#define _PhotonMapRenderer_

#include "Common.h"

class Photon_map;
class PhotonFilter;
struct HitRecord;
class Scene;

class PhotonMapRenderer
{
public:
    
    struct Config
    {
        u32 screenWidth;
        u32 screenHeight;
        u32 nSubPixelsSqrt;
        u32 nPhotons;
        u32 nEstimatePhotons;
        float estimateDist;
        float coneFilterK;
        u32 maxPhotonBounce;
        u32 maxRayBounce;
    };
    
    //
    
    PhotonMapRenderer();
    ~PhotonMapRenderer();
    void SetConfig(const Config& config);
    Config GetDefaultConfig();
    void SetCamera(const Ray& camera, real fovY);
    Vec3* Run(const Scene& scene);

private:
    
    void PhotonTracing(const Ray& r, float power[3], int depth);
    bool Intersect(const Ray& r, HitRecord& out);
    Vec3* RayTracing();
    Vec3 Irradiance(const Ray &r, int depth);
    
    struct Config config_;
    struct Config defaultConfig_;
    Photon_map* pPhotonMap_;
    unsigned short xi_[3];
    PhotonFilter* pFilter_;
    Ray camera_;
    real fovY_;
    const Scene* pScene_;
};

#endif
