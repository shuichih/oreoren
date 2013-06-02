#ifndef BVH_h
#define BVH_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"

/**
 * Bounding Volume Hierarchy
 */
class BVH : public IShape
{
public:
    static int QSplit(const IShape** pShapes, int nShapes, float pivot, int axis);
    
public:
    BVH();
    BVH(const IShape** pShapes, int nShapes);
    BVH(const IShape* s1, const IShape* s2);
    BVH(const IShape* s1, const IShape* s2, const BBox& bbox);
    virtual ~BVH();
    
    virtual ShapeType GetType() const;
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    virtual int RayCast(std::vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const;
    virtual bool IsBVH() const { return true; };
    void LimitMinScale(float minScale);
    // virtual bool ShadowHit(const Ray& r, float tmin, float tmax) const;
    
private:
    const IShape* BuildBranch(const IShape** pShapes, int nShapes, int axis = 0);

public:
    const IShape* pLeft_;
    const IShape* pRight_;
    BBox bbox_;
};

#endif
