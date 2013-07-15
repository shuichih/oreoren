#include "Scene.h"
#include "BVH.h"
#include <cassert>

using namespace std;

BVH::BVH()
: pLeft_(NULL)
, pRight_(NULL)
{
}

BVH::~BVH()
{
    if (pLeft_->IsBVH()) delete pLeft_;
    if (pRight_->IsBVH()) delete pRight_;
}

inline BVH::BVH(const IShape* pLeft, const IShape* pRight, const BBox& bbox)
: bbox_(bbox)
, pLeft_(pLeft)
, pRight_(pRight)
{
}

inline BVH::BVH(const IShape* pLeft, const IShape* pRight)
{
    pLeft_ = pLeft;
    pRight_ = pRight;
    // LeftとRightを包含するBBoxを作る
    bbox_ = BBox::Surround(pLeft->BoundingBox(), pRight->BoundingBox());
}

BVH::BVH(const IShape** pShapes, int nShapes)
{
    if (nShapes == 0) {
        pLeft_ = NULL;
        pRight_ = NULL;
        return;
    }
    if (nShapes == 1) {
        *this = BVH(pShapes[0], pShapes[0]);
        return;
    }
    else if (nShapes == 2) {
        *this = BVH(pShapes[0], pShapes[1]);
        return;
    }
    
    // qsplitのピボットとして使用するためにBBoxの中間点を求める
    bbox_ = pShapes[0]->BoundingBox();
    for (int i=1; i<nShapes; i++) {
        bbox_ = BBox::Surround(bbox_, pShapes[i]->BoundingBox());
    }
    
    // 空間的な中心をpivotとする
    Vec3 pivot = (bbox_.Max() + bbox_.Min()) * 0.5f;
    
    int mid_point = QSplit(pShapes, nShapes, pivot.x, 0); // まずX方向で
    
    // 子ノードをビルド
    pLeft_ = BuildBranch(pShapes, mid_point, 1); // Y方向で
    pRight_ = BuildBranch(&pShapes[mid_point], nShapes - mid_point, 1); // Y方向で
}

const IShape* BVH::BuildBranch(const IShape** pShapes, int nShapes, int axis)
{
    if (nShapes == 1) {
        int nChild = pShapes[0]->GetChildNum();
        return (nChild) ? BuildBranch(pShapes[0]->GetChildren(), nChild) : pShapes[0];
    }
    if (nShapes == 2) {
        int nChild = pShapes[0]->GetChildNum();
        const IShape* pLeft = (nChild) ? BuildBranch(pShapes[0]->GetChildren(), nChild) : pShapes[0];
        nChild = pShapes[1]->GetChildNum();
        const IShape* pRight = (nChild) ? BuildBranch(pShapes[1]->GetChildren(), nChild) : pShapes[1];
        
        return new BVH(pLeft, pRight);
    }
    
    BBox bbox = pShapes[0]->BoundingBox();
    for (int i=1; i<nShapes; i++) {
        bbox = BBox::Surround(bbox, pShapes[i]->BoundingBox());
    }
    
    // 空間的な中心をpivotとする
    Vec3 pivot = (bbox.Max() + bbox.Min()) * 0.5f;
    
    // axis方向で分割
    int mid_point = QSplit(pShapes, nShapes, pivot.x, 0); // まずX方向で
    
    // 子ノードをビルド
    const IShape* pLeft = BuildBranch(pShapes, mid_point, (axis+1)%3);
    const IShape* pRight = BuildBranch(&pShapes[mid_point], nShapes - mid_point, (axis+1)%3);
    
    return new BVH(pLeft, pRight, bbox);
}

ShapeType BVH::GetType() const
{
    return ST_BBVH;
}

BBox BVH::BoundingBox() const
{
    return bbox_;
}

bool BVH::Intersect(const Ray &r, float tmin, float tmax, HitRecord& rec) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return false;
    
    rec.t = tmax;
    
    bool isahit1 = pLeft_->Intersect(r, tmin, tmax, rec);
    bool isahit2 = pRight_->Intersect(r, tmin, rec.t, rec); // rec.t Leftが当たってたらそれより手前のみ判定
    
    return (isahit1 || isahit2);
}

int BVH::RayCast(vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return nHits;
    
    nHits = pLeft_->RayCast(hits, nHits, r, tmin, tmax);
    nHits = pRight_->RayCast(hits, nHits, r, tmin, tmax);
    
    return nHits;
}


// 直接照明実装するとき使用?
//bool BVH::ShadowHit(const Ray &r, float tmin, float tmax)
//{
//    if (!(bbox.RayIntersect(const Ray &r, float tmin, float tmax)))
//        return false;
//
//    if (pRight_->shadowHit(r, tmin, tmax)) return true;
//    return pLeft_->ShadowHit(r, tmin, tmax);
//}

void BVH::LimitMinScale(float minScale)
{
    if ((bbox_.Max().x - bbox_.Min().x) < minScale) {
        bbox_.Min().x -= minScale * 0.5f;
        bbox_.Max().x += minScale * 0.5f;
    }
    if ((bbox_.Max().y - bbox_.Min().y) < minScale) {
        bbox_.Min().y -= minScale * 0.5f;
        bbox_.Max().y += minScale * 0.5f;
    }
    if ((bbox_.Max().z - bbox_.Min().z) < minScale) {
        bbox_.Min().z -= minScale * 0.5f;
        bbox_.Max().z += minScale * 0.5f;
    }
    
    if (pLeft_->IsBVH()) ((BVH*)pLeft_)->LimitMinScale(minScale);
    if (pRight_->IsBVH()) ((BVH*)pRight_)->LimitMinScale(minScale);
}

int BVH::QSplit(const IShape** pShapes, int nShapes, float pivot, int axis)
{
    BBox bbox;
    double centroid;
    int mid_idx = 0;
    for (int i=0; i<nShapes; i++) {
        bbox = pShapes[i]->BoundingBox();
        centroid = ((bbox.Min()).e[axis] + (bbox.Max()).e[axis]) * 0.5f;
        if (centroid < pivot)
        {
            // pivotを基準にmid_idxの左右に集める
            const IShape* pTemp = pShapes[i];
            pShapes[i] = pShapes[mid_idx];
            pShapes[mid_idx] = pTemp;
            mid_idx++;
        }
    }
    
    // 全てのshapeがpivotの右または左にあったら、mid_idxを配列の真ん中にする
    if (mid_idx == 0 || mid_idx == nShapes) {
        mid_idx = nShapes / 2;
    }
    
    return mid_idx;
}

void BVH::SetMaterial(Material* pMtl)
{
}

Material* BVH::GetMaterial() const
{
    return NULL;
}
