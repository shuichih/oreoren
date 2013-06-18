#ifndef _Scene_H_
#define _Scene_H_

#include "Common.h"
#include "BBox.h"
#include <vector>
#include "BVHType.h"


class IShape;
class LightSource;
class BVH;

// material type
// @todo introduce Material class
enum Refl_t
{
    DIFF,
    SPEC,
    REFR,
    PHONGMETAL,
    LIGHT
};

enum ColorUnit
{
    CU_Mesh,
    CU_Face,
    CU_Vertex
};

enum ShapeType
{
    ST_SPHERE,
    ST_TRIANGLE,
    ST_MESH_TRIANGLE,
    ST_MESH,
    ST_BBVH,
    ST_QBVH_SISD,
    ST_QBVH_SIMD,
};

struct HitRecord
{
    real t;
    Vec3 normal;
    Vec3 color;
    Refl_t refl;
    bool hitLit;
    const IShape* pShape;
    HitRecord()
    : hitLit(false)
    {}
};

class IShape
{
public:
    virtual ~IShape();
    
    virtual BBox BoundingBox() const = 0;
    virtual int RayCast(std::vector<HitRecord>& shapes, int nHits, const Ray& r, float tmin, float tmax) const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const = 0;
    virtual bool IsBVH() const { return false; };
    virtual int GetChildNum() const { return 0; }
    virtual const IShape** GetChildren() const { return NULL; }
    virtual ShapeType GetType() const = 0;
};

class Sphere : public IShape
{
public:
    Sphere();
    Sphere(const Sphere& rhs);
    Sphere(real radius_, Vec3 p_, Vec3 c_, Refl_t refl_);
    virtual ~Sphere();
    
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    
    real rad;   // radius
    Vec3 p;      // position
    Vec3 c;      // color
    Refl_t refl;    // reflection type (DIFFuse, SPECular, REFRactive)
};

class Triangle : public IShape
{
public:
    Triangle();
    Triangle(const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, const RGB& _color, Refl_t _refl);
    virtual ~Triangle();
    
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
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
    RGB color;
};

class Mesh;

class MeshTriangle : public IShape
{
public:
    MeshTriangle();
    virtual ~MeshTriangle();

    void CalcFaceNormal();
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    //void Reverse();
    
    u32 indices[3];
    Vec3 normal;
    RGB color_;
    Mesh* pMesh;
};

class Mesh : public IShape
{
public:
    Mesh(u32 nVertices, u32 nFaces);
    virtual ~Mesh();

    void SetUseFaceNormal(bool useFaceNormal);
    bool GetUseFaceNormal();
    void CalcFaceNormals();
    void CalcVertexNormals();
    void CalcBoundingBox();
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    virtual int GetChildNum() const;
    virtual const IShape** GetChildren() const;  // for BVH
    
    void scale(Vec3& scl);
    void scale(real x, real y, real z);
    void translate(Vec3& transl);
    void translate(real x, real y, real z);
    void rotateXYZ(Vec3& rot);
    //void ReverseFaces();
    
    Vertex*         pVertices;
    MeshTriangle*   pFaces;
    const IShape**  ppFaces; // for BVH
    u32             nVertices;
    u32             nFaces;
    BBox            bbox_;
    Refl_t          material_;
    RGB             color_;
    ColorUnit       colorUnit_;
    bool            useFaceNormal_;
};

class Scene
{
public:
    Scene();
    ~Scene();
    
    void AddLightSource(const LightSource* pLitSrc);
    void AddShape(const IShape* pShape);
    bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    int RayCast(std::vector<HitRecord>& shapes, int nHits, const Ray& r, float tmin, float tmax) const;
    void BuildBVH(BVHType bvhType);
    
    std::vector<const LightSource*> litSrcs_;
    std::vector<const IShape*> shapes_;
    IShape* pBVH_;
};

#endif
