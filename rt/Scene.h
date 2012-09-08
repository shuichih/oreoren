

#ifndef _Scene_H_
#define _Scene_H_

#include "Common.h"
#include "LightSource.h"
#include <vector>

enum Refl_t
{
    DIFF, SPEC, REFR, PHONGMETAL
}; // material types, used in radiance()

struct HitRecord
{
    real t;
    Vec3 normal;
    Vec3 color;
    Refl_t refl;
};

class Shape
{
public:
    virtual ~Shape();
    virtual bool intersect(const Ray& r, HitRecord& rec) const = 0;
};

class Sphere : public Shape
{
public:
    Sphere();
    Sphere(const Sphere& rhs);
    Sphere(real rad_, Vec3 p_, Vec3 c_, Refl_t refl_);
    virtual ~Sphere();
    virtual bool intersect(const Ray &r, HitRecord& rec) const;

    real rad;   // radius
    Vec3 p;      // position
    Vec3 c;      // color
    Refl_t refl;    // reflection type (DIFFuse, SPECular, REFRactive)
};

class Triangle : public Shape
{
public:
    Triangle(const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, const RGB& _color, Refl_t _refl);
    Triangle();
    virtual ~Triangle();
    virtual bool intersect(const Ray& r, HitRecord& rec) const;
    void CalcNormal();
    
    Vec3 p0;
    Vec3 p1;
    Vec3 p2;
    Vec3 normal;
    RGB color;
    Refl_t refl;
};

struct Vertex
{
    Vec3 pos;
    Vec3 normal;
};

class Mesh;

class MeshTriangle : public Shape
{
public:
    MeshTriangle();
    virtual ~MeshTriangle();

    virtual bool intersect(const Ray& r, HitRecord& rec) const;
    
    u32 indices[3];
    Vec3 normal;
    Mesh* pMesh;
};

class Mesh : public Shape
{
public:
    Mesh(u32 nVertices, u32 nFaces);
    virtual ~Mesh();

    virtual bool intersect(const Ray& r, HitRecord& rec) const;
    
    void scale(Vec3 scl);
    void scale(real x, real y, real z);
    void translate(Vec3 transl);
    void translate(real x, real y, real z);

    Vertex*         pVertices;
    MeshTriangle*   pFaces;
    u32             nVertices;
    u32             nFaces;
};

class Scene
{
public:
    Scene();
    ~Scene();
    
    void AddLightSource(const LightSource* pLitSrc);
    void AddShape(const Shape* pShape);
    
    std::vector<const LightSource*> litSrcs_;
    std::vector<const Shape*> shapes_;
};

#endif
