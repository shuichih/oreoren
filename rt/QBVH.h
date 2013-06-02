#ifndef QBVH_h
#define QBVH_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"

/**
 * Bounding Volume Hierarchy
 */
class QBVH : public IShape
{
public:
    static int QSplit(const IShape** pShapes, int nShapes, float pivot, int axis);
    
public:
    QBVH();
    QBVH(const IShape* s1, const IShape* s2);
    QBVH(const IShape* s1, const IShape* s2, const BBox& bbox);
    virtual ~QBVH();
    
    void Build(const IShape** pShapes, int nShapes);
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    virtual int RayCast(std::vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const;
    virtual bool IsQBVH() const { return true; };
    //void LimitMinScale(float minScale);
    
private:
    
    // SISD version
    struct SISD_QBVH_NODE
    {
        //float bboxes_[2][3][4]; // min-max xyz 4 children
        BBox bboxes[4]; // 24*4=96 bytes
        int children[4]; // 16bytes
        int axes; // pPrimsSISDを含めて128bytesに収めるためにまとめた
        int reserved0;

        // 末端ノードを作成するときに、SIMD処理できないMeshTriangle以外の
        // プリミティブがあれば、それらを格納する配列を作成する。
        int nOtherPrims;
        int iOtherPrims;
        //const IShape* pOtherPrims; // 8bytes(x64)
        int reserved1;
        
        SISD_QBVH_NODE()
        {
            Reset();
        }
        
        void Reset()
        {
            for (int i=0; i<4; i++) {
                children[i] = INT_MIN;
            }
            axes = 0;
            nOtherPrims = 0;
            iOtherPrims = 0;
        }
    };
    
    struct SISD_TRIANGLE
    {
        Vec3 p[3];
        const MeshTriangle *pMeshTriangle;
    };
    
    /*
    struct Leaf
    {
        SISD_TRIANGLE* pTriangles;
        IShape* pOtherPrims;
        int nTriangles;
        int nOtherPrims;
    };
    */
    
#if 0
    // SIMD version. total size == 128bytes
    __m128 bboxes_[2][3]; // min-max xyz 4 children
    int children_[4];
    int axes_; // pPrimsSISDを含めて128bytesに収めるためにまとめた
    int reserved_;
    
    // 末端ノードを作成するときに、SIMD処理できないMeshTriangle以外の
    // プリミティブがあれば、それらを格納する配列を作成する。
    IShape* pPrimsSISD_;
#endif

    //
    
    void Reset();
    bool IsLeafType(ShapeType shapeType);
    int CountLeafShapes(const IShape** pShapes, int nShapes);
    const IShape** FlattenLeafShapes(const IShape** ppFlatten, const IShape** ppShapes, int nShapes);
    BBox SurroundBBox(const IShape** pShapes, int nShapes);
    void BuildBranch(SISD_QBVH_NODE& rNode, const IShape** pShapes, int nShapes);
    void BuildLeaf(int& nMeshTri, int& nOtherPrims, const IShape** pShapes, int nShapes);
    bool IntersectBranch(SISD_QBVH_NODE& rNode, const Ray &r, float tmin, HitRecord& rec) const;
    bool IntersectLeaf(SISD_TRIANGLE* pTri, int nTri, const Ray& r, float tmin, HitRecord& rec) const;
    bool IntersectTriangle(SISD_TRIANGLE& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    int RayCastBranch(std::vector<HitRecord>& hits, int nHits, SISD_QBVH_NODE& rNode, const Ray &r, float tmin, float tmax) const;
    int RayCastLeaf(std::vector<HitRecord>& hits, int nHits, SISD_TRIANGLE* pTri, int nTri, const Ray &r, float tmin, float tmax) const;
    
    //

    SISD_QBVH_NODE* pNodes_;
    SISD_TRIANGLE* pTriangles_;
    const IShape** ppOtherPrims_;
    int nNodes_;
    int nTriangles_;
    int nOtherPrims_;
    int iNodes_;
    int iTriangles_;
    int iOtherPrims_;
    BBox bbox_;
};

#endif
