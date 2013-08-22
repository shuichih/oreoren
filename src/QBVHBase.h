#ifndef QBVHBase_h
#define QBVHBase_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"

/**
 * SISD version QBVH node.
 */
struct SISD_QBVH_NODE
{
    BBox bboxes[4]; // 24*4=96 bytes
    int children[4]; // 16bytes
    int axis0; // 4bytes
    int axis1; // 4bytes
    int axis2; // 4bytes
    int reserved;
    
    SISD_QBVH_NODE()
    {
        Reset();
    }
    
    void Reset()
    {
        for (int i=0; i<4; i++) {
            children[i] = INT_MIN;
        }
        axis0 = 0;
        axis1 = 0;
        axis2 = 0;
        reserved = 0;
    }
};

/**
    SIMD version.
    total size == 128bytes
*/
struct SIMD_QBVH_NODE
{
    __m128 bboxes[2][3]; // min-max xyz 4 children
    int children[4]; // 16bytes
    int axis0; // 4bytes
    int axis1; // 4bytes
    int axis2; // 4bytes
    int reserved; // 4bytes
    
    SIMD_QBVH_NODE()
    {
        Reset();
    }
    
    void Reset()
    {
        for (int i=0; i<2; i++) {
            for (int j=0; j<3; j++) {
                bboxes[i][j] = _mm_set1_ps(0);
            }
        }
        for (int i=0; i<4; i++) {
            children[i] = INT_MIN;
        }
        axis0 = 0;
        axis1 = 0;
        axis2 = 0;
        reserved = 0;
    }
};

/**
* Quad Bounding Volume Hierarchy Base.
 */
template <typename NODE_T>
class QBVHBase : public ShapeBase
{
public:
    static int QSplit(const IShape** pShapes, int nShapes, float pivot, int axis);
    
public:
    QBVHBase();
    virtual ~QBVHBase();
    
    virtual BBox BoundingBox() const;
    void Build(const IShape** pShapes, int nShapes);
    
protected:
    
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
    
    // child index access
    void SetChildEmpty(int& child);
    void SetChildNode(int& child, int index);
    void SetChildLeaf(int& child, int index);
    bool IsChildEmpty(int child) const;
    int GetChildIndex(int child) const;
    bool IsChildNode(int child) const;
    
    // for build
    void Reset();
    bool IsLeafType(ShapeType shapeType);
    bool IsOtherPrimType(ShapeType shapeType);
    bool IsMeshTriangleType(ShapeType shapeType);
    int CountLeafShapes(const IShape** pShapes, int nShapes);
    const IShape** FlattenLeafShapes(const IShape** ppFlatten, const IShape** ppShapes, int nShapes);
    BBox SurroundBBox(const IShape** pShapes, int nShapes);
    int LargestAxis(const BBox& bbox);
    int BuildLeaf(u8& nMeshTris, u8& nOtherPrims, const IShape** ppShapes, int nShapes);
    void ShrinkBuffersToFit();
    int AddNewNode();
    Leaf& AddNewLeaf();
    
    // for Intersection
    bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    bool IntersectLeaf(Leaf& leaf, const Ray& r, float tmin, HitRecord& rec) const;
    bool IntersectTriangle(Triangle& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    
    // for RayCast
    int RayCast(std::vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const;
    int RayCastLeaf(std::vector<HitRecord>& hits, int nHits, Leaf& leaf, const Ray &r, float tmin, float tmax) const;
    
    // pure virtual
    virtual void BuildBranch(int iCurrNode, const IShape** pShapes, int nShapes) = 0;
    virtual void MakeWholeBBox() = 0;
    virtual bool IntersectBranch(NODE_T& rNode, const Ray &r, float tmin, HitRecord& rec) const = 0;
    virtual int RayCastBranch(std::vector<HitRecord>& hits, int nHits, NODE_T& rNode, const Ray &r, float tmin, float tmax) const = 0;
    
    //

    NODE_T* pNodes_;
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
