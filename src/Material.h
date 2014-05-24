#ifndef Material_H
#define Material_H

#include "Common.h"

struct HitRecord;
class Lambertian;

//--------------------------------------------------------------------------------
class Material
{
public:
    virtual RGB PhotonMapShade(HitRecord& hr);
};

//--------------------------------------------------------------------------------
class Matte : public Material
{
public:
    Matte();
    void SetKa(float k);
    void SetKd(float k);
    void SetCd(const RGB& c);
    virtual RGB PhotonMapShade(HitRecord& hr);
    
private:
    Lambertian* pAmbientBrdf_;
    Lambertian* pDiffuseBrdf_;
};

#endif


