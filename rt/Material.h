#ifndef Material_H
#define Material_H

#include "Common.h"
#include <string>

enum Refl_t
{
    DIFF,
    SPEC,
    REFR,
    PHONGMETAL,
    LIGHT
};

class Material
{
public:
    Material()
    {
        Reset();
    }
    void Reset()
    {
        name = "noname";
        refl = DIFF;
        color = Vec3(.5f, .5f, .5f);
        refractiveIndex = 1.f;
    }
    std::string name;
    Refl_t refl;
    RGB color;
    float refractiveIndex;
};

#endif

