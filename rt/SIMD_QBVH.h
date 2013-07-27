#ifndef SIMD_QBVH_h
#define SIMD_QBVH_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"
#include "QBVHBase.h"

/**
 * Quad Bounding Volume Hierarchy using SIMD instructions.
 */
class SIMD_QBVH : public QBVHBase<SIMD_QBVH_NODE>
{
public:
    SIMD_QBVH();
    SIMD_QBVH(const IShape* s1, const IShape* s2);
    SIMD_QBVH(const IShape* s1, const IShape* s2, const BBox& bbox);
    virtual ~SIMD_QBVH();
    
    virtual ShapeType GetType() const;
    
protected:
    virtual void BuildBranch(int iCurrNode, const IShape** pShapes, int nShapes);
    virtual void MakeWholeBBox();
    
private:
    bool IntersectBranch(SIMD_QBVH_NODE& rNode, const Ray &r, float tmin, HitRecord& rec) const;
    int RayCastBranch(std::vector<HitRecord>& hits, int nHits, SIMD_QBVH_NODE& rNode, const Ray &r, float tmin, float tmax) const;
    inline int IntersectSIMD(
        const __m128 bboxes[2][3],  // min-max[2] of xyz[3] of boxes[4]
        const __m128 orig[3],       // ray origin, xyz[3]
        const __m128 idir[3],       // ray inverse direction, xyz[3]
        const int sign[3],          // ray xyz direction -> +:0,-:1
        __m128 tmin, __m128 tmax    // ray range tmin-tmax
    ) const;
};

#endif
