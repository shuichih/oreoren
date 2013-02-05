#ifndef BVH_h
#define BVH_h

#include "Common.h"
#include "Scene.h"
#include "BBox.h"

/**
 * Bounding Volume Hierarchy
 */
class BVH : public Shape
{
public:
    static int QSplit(const Shape** pShapes, int nShapes, float pivot, int axis);
    
public:
    BVH();
    BVH(const Shape** pShapes, int nShapes);
    BVH(const Shape* s1, const Shape* s2);
    BVH(const Shape* s1, const Shape* s2, const BBox& bbox);
    virtual ~BVH();
    
    virtual BBox BoundingBox() const;
    virtual bool Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const;
    virtual int RayCast(std::vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const;
    virtual bool IsBVH() const { return true; };
    void LimitMinScale(float minScale);
    // virtual bool ShadowHit(const Ray& r, float tmin, float tmax) const;
    
private:
    const Shape* BuildBranch(const Shape** pShapes, int nShapes, int axis = 0);

public:
    const Shape* pLeft_;
    const Shape* pRight_;
    BBox bbox_;
};

#endif
