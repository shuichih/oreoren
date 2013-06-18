#ifndef QBVHBase_h
#define QBVHBase_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"

/**
 * Quad Bounding Volume Hierarchy Base.
 */
class QBVHBase : public IShape
{
public:
    static int QSplit(const IShape** pShapes, int nShapes, float pivot, int axis);
    
public:
    QBVHBase();
    QBVHBase(const IShape* s1, const IShape* s2);
    QBVHBase(const IShape* s1, const IShape* s2, const BBox& bbox);
    virtual ~QBVHBase();
    
    virtual BBox BoundingBox() const;
    
private:
    
    struct Triangle
    {
        Vec3 p[3];
        const MeshTriangle *pMeshTriangle;
    };
    
    struct Leaf
    {
        int iTriangles;
        int iOtherPrims;
        u8 nTriangles;
        u8 nOtherPrims;
    };
    
    //
    
    void Reset();
    bool IsLeafType(ShapeType shapeType);
    bool IsOtherPrimType(ShapeType shapeType);
    int CountLeafShapes(const IShape** pShapes, int nShapes);
    const IShape** FlattenLeafShapes(const IShape** ppFlatten, const IShape** ppShapes, int nShapes);
    BBox SurroundBBox(const IShape** pShapes, int nShapes);
    int BuildLeaf(u8& nMeshTris, u8& nOtherPrims, const IShape** ppShapes, int nShapes);
    int BuildOtherPrimitive(const IShape** pShapes, int nShapes);
    Leaf& GetNewLeaf();
    bool IntersectLeaf(Leaf& leaf, const Ray& r, float tmin, HitRecord& rec) const;
    bool IntersectTriangle(Triangle& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    int RayCastLeaf(std::vector<HitRecord>& hits, int nHits, Leaf& leaf, const Ray &r, float tmin, float tmax) const;
    int LargestAxis(const BBox& bbox);
    
    //

    Triangle* pTriangles_;
    const IShape** ppOtherPrims_;
    Leaf* pLeaves_;
    int nNodeCapacity_;
    int nLeafCapacity_;
    int nTriangles_;
    int nOtherPrims_;
    int iNodes_;
    int iLeaves_;
    int iTriangles_;
    int iOtherPrims_;
    BBox bbox_;
    mutable int depth_;
};

#endif
