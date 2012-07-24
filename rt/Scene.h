

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
    virtual bool intersect(const Ray& r, HitRecord& rec) const;

    Vec p0;
    Vec p1;
    Vec p2;
    Vec normal;
    RGB color;
    Refl_t refl;
};

struct Vertex
{
    Vec pos;
    Vec normal;
};

class Mesh;

class MeshTriangle : public Shape
{
public:
    MeshTriangle();
    ~MeshTriangle();

    virtual bool intersect(const Ray& r, HitRecord& rec) const;
    
    u32 indices[3];
    Vec normal;
    Mesh* pMesh;
};

class Mesh : public Shape
{
public:
    Mesh(u32 nVertices, u32 nFaces);
    ~Mesh();

    virtual bool intersect(const Ray& r, HitRecord& rec) const;
    
    void scale(real x, real y, real z);
    void translate(real x, real y, real z);

    Vertex*         pVertices;
    MeshTriangle*   pFaces;
    u32             nVertices;
    u32             nFaces;
};

#endif
