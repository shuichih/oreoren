#ifndef _Scene_H_
#define _Scene_H_

#include "Common.h"
#include "BBox.h"
#include <vector>
#include <string>
#include <map>
#include "BVHType.h"
#include "Material.h"


class IShape;
class LightSource;
class BVH;

// material type
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
    RGB color;
    Material* pMaterial;
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
    virtual void SetMaterial(Material* pMtl) = 0;
    virtual Material* GetMaterial() const = 0;
};

class ShapeBase : public IShape
{
public:
    virtual ~ShapeBase();
    
    ShapeBase(Material* pMtl);
    virtual void SetMaterial(Material* pMtl);
    virtual Material* GetMaterial() const;
    
    Material* pMaterial;
};

class Sphere : public ShapeBase
{
public:
    Sphere(real radius, Vec3 position, Material* pMtl);
    virtual ~Sphere();
    
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    
    real radius;
    Vec3 position;
};

class Triangle : public ShapeBase
{
public:
    Triangle();
    Triangle(const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, Material* pMaterial);
    virtual ~Triangle();
    
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    void CalcNormal();
    
    Vec3 p0;
    Vec3 p1;
    Vec3 p2;
    Vec3 normal;
};

struct Vertex
{
    Vec3 pos;
    Vec3 normal;
    RGB color;
};

class Mesh;

class MeshTriangle : public ShapeBase
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
    Mesh* pMesh;
};

class Mesh : public ShapeBase
{
public:
    Mesh(u32 nVertices, u32 nFaces, Material* pMtl);
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
    ColorUnit       colorUnit_;
    bool            useFaceNormal_;
};

class Scene
{
public:
    typedef std::map<std::string, Material*> MaterialMap;
    
    Scene();
    ~Scene();
    
    void AddLightSource(const LightSource* pLitSrc);
    void AddShape(IShape* pShape);
    bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    int RayCast(std::vector<HitRecord>& shapes, int nHits, const Ray& r, float tmin, float tmax) const;
    void BuildBVH(BVHType bvhType);
    inline u32 GetShapeNum() const { return (u32)shapes_.size(); }
    inline const IShape* GetShape(u32 index) const { return shapes_[index]; }
    inline u32 GetLightNum() const { return (u32)litSrcs_.size(); }
    inline const LightSource* GetLight(u32 index) const { return litSrcs_[index]; }
    Material* GetMaterial(std::string name);
    Material* GetDefaultMaterial();
    Material* GetLightMaterial();
    bool AddMaterial(std::string name, Material* pMtl);
    
private:
    Material defaultMaterial_;
    Material lightMaterial_;
    MaterialMap materialMap_;
    std::vector<const LightSource*> litSrcs_;
    std::vector<const IShape*> shapes_;
    IShape* pBVH_;
};

#endif
