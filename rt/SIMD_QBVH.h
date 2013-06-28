#ifndef SIMD_QBVH_h
#define SIMD_QBVH_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"

/**
 * Quad Bounding Volume Hierarchy using SIMD instructions.
 */
class SIMD_QBVH : public IShape
{
public:
    static int QSplit(const IShape** pShapes, int nShapes, float pivot, int axis);
    
public:
    SIMD_QBVH();
    SIMD_QBVH(const IShape* s1, const IShape* s2);
    SIMD_QBVH(const IShape* s1, const IShape* s2, const BBox& bbox);
    virtual ~SIMD_QBVH();
    
    void Build(const IShape** pShapes, int nShapes);
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    virtual int RayCast(std::vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const;
    //void LimitMinScale(float minScale);
    
private:
    
    // SIMD version. total size == 128bytes
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
    
    struct SIMD_TRIANGLE
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
    BBox SurroundBBox(__m128 bboxes[2][3]);
    BBox SurroundBBox(const IShape** ppShapes, int nShapes);
    void BuildBranch(int iNode, const IShape** pShapes, int nShapes);
    int BuildLeaf(u8& nMeshTris, u8& nOtherPrims, const IShape** ppShapes, int nShapes);
    int BuildOtherPrimitive(const IShape** pShapes, int nShapes);
    int AddNewNode();
    Leaf& AddNewLeaf();
    bool IntersectBranch(SIMD_QBVH_NODE& rNode, const Ray &r, float tmin, HitRecord& rec) const;
    bool IntersectLeaf(Leaf& leaf, const Ray& r, float tmin, HitRecord& rec) const;
    bool IntersectTriangle(SIMD_TRIANGLE& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    int IntersectSIMD(
        const __m128 bboxes[2][3],  // min-max[2] of xyz[3] of boxes[4]
        const __m128 org[3],        // ray origin
        const __m128 idir[3],       // ray inverse direction
        const int sign[3],          // ray xyz direction -> +:0,-:1
        __m128 tmin, __m128 tmax    // ray range tmin-tmax
    ) const;
    int RayCastBranch(std::vector<HitRecord>& hits, int nHits, SIMD_QBVH_NODE& rNode, const Ray &r, float tmin, float tmax) const;
    int RayCastLeaf(std::vector<HitRecord>& hits, int nHits, Leaf& leaf, const Ray &r, float tmin, float tmax) const;
    int LargestAxis(const BBox& bbox);
    
    //

    SIMD_QBVH_NODE* pNodes_;
    SIMD_TRIANGLE* pTriangles_;
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
