#ifndef OldMaterial_H
#define OldMaterial_H

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

class OldMaterial
{
public:
    OldMaterial(std::string name, Refl_t refl, RGB color, float refIdx)
    {
        this->name = name;
        this->refl = refl;
        this->color = color;
        this->refractiveIndex = refIdx;
    }
    OldMaterial()
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

