#ifndef _PhotonMapRenderer_
#define _PhotonMapRenderer_

#include "Common.h"

class Photon_map;

class PhotonMapRenderer
{
public:
    PhotonMapRenderer();
    ~PhotonMapRenderer();
    unsigned char* Run(Photon_map* pPhotonMap, unsigned int nPhotons, int w, int h);

private:
    
    void PhotonTracing(const Ray& r, float power[3], int depth);
    bool Intersect(const Ray& r, double& t, int& id);
    unsigned char* RayTracing(int w, int h, int samps);
    Vec Irradiance(const Ray &r, int depth);
    
    unsigned short xi_[3];
    Photon_map* pPhotonMap_;
    unsigned int nPhotons_;
    unsigned int nEstimatePhotons_;
    float estimateDist_;
    PhotonFilter* pFilter_;
};

#endif
