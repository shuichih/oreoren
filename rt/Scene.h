

#ifndef _Scene_H_
#define _Scene_H_

#include "Common.h"

enum Refl_t
{
    DIFF, SPEC, REFR, PHONGMETAL
}; // material types, used in radiance()

struct HitRecord
{
    real t;
    Vec normal;
    Vec color;
    Refl_t refl;
};

class Shape
{
public:
    virtual bool intersect(const Ray& r, HitRecord& rec) const = 0;
};

class Sphere : public Shape
{
public:
    Sphere(real rad_, Vec p_, Vec c_, Refl_t refl_);
    virtual bool intersect(const Ray &r, HitRecord& rec) const;

    real rad;   // radius
    Vec p;      // position
    Vec c;      // color
    Refl_t refl;    // reflection type (DIFFuse, SPECular, REFRactive)
};

class Triangle : public Shape
{
public:
    Triangle(const Vec& _p0, const Vec& _p1, const Vec& _p2, const RGB& _color, Refl_t _refl);
    bool intersect(const Ray& r, HitRecord& rec) const;

    Vec p0;
    Vec p1;
    Vec p2;
    Vec normal;
    RGB color;
    Refl_t refl;
};

#endif
