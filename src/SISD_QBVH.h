#ifndef SISD_QBVH_h
#define SISD_QBVH_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"
#include "QBVHBase.h"

    
/**
 * Quad Bounding Volume Hierarchy.
 */
class SISD_QBVH : public QBVHBase<SISD_QBVH_NODE>
{
public:
    
    SISD_QBVH();
    SISD_QBVH(const IShape* s1, const IShape* s2);
    SISD_QBVH(const IShape* s1, const IShape* s2, const BBox& bbox);
    virtual ~SISD_QBVH();
    
    virtual ShapeType GetType() const;
    
protected:
    virtual void BuildBranch(int iCurrNode, const IShape** pShapes, int nShapes);
    virtual void MakeWholeBBox();
    
private:
    bool IntersectBranch(SISD_QBVH_NODE& rNode, const Ray &r, float tmin, HitRecord& rec) const;
    int RayCastBranch(std::vector<HitRecord>& hits, int nHits, SISD_QBVH_NODE& rNode, const Ray &r, float tmin, float tmax) const;
};

#endif
