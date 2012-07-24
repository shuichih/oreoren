#ifndef _PhotonMapRenderer_
#define _PhotonMapRenderer_

#include "Common.h"

class Photon_map;
struct HitRecord;

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
        const PhotonFilter* pFilter;
    };
    
    //
    
    PhotonMapRenderer();
    ~PhotonMapRenderer();
    void SetConfig(const Config& config);
    Config GetDefaultConfig();
    Vec* Run();

private:
    
    void PhotonTracing(const Ray& r, float power[3], int depth);
    bool Intersect(const Ray& r, HitRecord& out);
    Vec* RayTracing();
    Vec Irradiance(const Ray &r, int depth);
    
    struct Config config_;
    struct Config defaultConfig_;
    Photon_map* pPhotonMap_;
    unsigned short xi_[3];
};

#endif
