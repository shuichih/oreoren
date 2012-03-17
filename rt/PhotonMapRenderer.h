#ifndef _PhotonMapRenderer_
#define _PhotonMapRenderer_

#include "Common.h"

class Photon_map;

class PhotonMapRenderer
{
public:
    
    struct Config
    {
        unsigned int screenWidth;
        unsigned int screenHeight;
        unsigned int nSamplePerPixel;
        unsigned int nPhotons;
        unsigned int nEstimatePhotons;
        float estimateDist;
        const PhotonFilter* pFilter;
    };
    
    //
    
    PhotonMapRenderer();
    ~PhotonMapRenderer();
    void SetConfig(const Config& config);
    Config GetDefaultConfig();
    unsigned char* Run();

private:
    
    void PhotonTracing(const Ray& r, float power[3], int depth);
    bool Intersect(const Ray& r, double& t, int& id);
    unsigned char* RayTracing();
    Vec Irradiance(const Ray &r, int depth);
    
    struct Config config_;
    struct Config defaultConfig_;
    Photon_map* pPhotonMap_;
    unsigned short xi_[3];
};

#endif
