#include "Scene.h"
#include "QBVHBase.h"
#include <cassert>
#include "simd.h"
#include "Ray.h"

using namespace std;

#define INITIAL_NODE_NUM 10000 // about 1.2MB
#define INITIAL_LEAF_NUM 30000 // about 0.72MB

QBVHBase::QBVHBase()
: pTriangles_(NULL)
, ppOtherPrims_(NULL)
, nNodeCapacity_(0)
, nTriangles_(0)
, nOtherPrims_(0)
, iNodes_(0)
, iTriangles_(0)
, iLeaves_(0)
, iOtherPrims_(0)
, depth_(0)
{
    nNodeCapacity_ = INITIAL_NODE_NUM;
    nLeafCapacity_ = INITIAL_LEAF_NUM;
    pLeaves_ = new Leaf[nLeafCapacity_];
}

QBVHBase::~QBVHBase()
{
}

// Build時に使うバッファを初期化
void QBVHBase::Reset()
{
    delete [] pTriangles_;
    delete [] ppOtherPrims_;
    pTriangles_ = NULL;
    ppOtherPrims_ = NULL;
}

// シェイプが葉になるタイプならtrueを返す
bool QBVHBase::IsLeafType(ShapeType shapeType)
{
    return (shapeType == ST_MESH_TRIANGLE
         || shapeType == ST_SPHERE
         || shapeType == ST_TRIANGLE);
}

// シェイプがMeshTriangle以外の葉になるタイプならtrueを返す
bool QBVHBase::IsOtherPrimType(ShapeType shapeType)
{
    return (shapeType == ST_SPHERE
         || shapeType == ST_TRIANGLE);
}

// 葉ノード総数をカウント
// MESH, QBVHは葉にならない
int QBVHBase::CountLeafShapes(const IShape** pShapes, int nShapes)
{
    int ret = 0;
    for (int i=0; i<nShapes; i++) {
        int nChildren = pShapes[i]->GetChildNum();
        if (nChildren == 0) {
            if (IsLeafType(pShapes[i]->GetType())) {
                ret++;
            }
        } else {
            ret += CountLeafShapes(pShapes[i]->GetChildren(), nChildren);
        }
    }
    return ret;
}

// 葉のシェイプだけをフラットな配列にする
const IShape** QBVHBase::FlattenLeafShapes(const IShape** ppFlatten, const IShape** ppShapes, int nShapes)
{
    for (int i=0; i<nShapes; i++) {
        int nChildren = ppShapes[i]->GetChildNum();
        if (nChildren == 0) {
            *ppFlatten = ppShapes[i];
            ppFlatten++;
        } else {
            ppFlatten = FlattenLeafShapes(ppFlatten, ppShapes[i]->GetChildren(), nChildren);
        }
    }
    return ppFlatten;
}

// 葉ノード構築
int QBVHBase::BuildLeaf(u8& nMeshTris, u8& nOtherPrims, const IShape** ppShapes, int nShapes)
{
    nMeshTris = 0;
    nOtherPrims = 0;
    for (int i=0; i<nShapes; i++) {
        if (ppShapes[i]->GetType() == ST_MESH_TRIANGLE) {
            const MeshTriangle* pMeshTri = reinterpret_cast<const MeshTriangle*>(ppShapes[i]);
            assert(iTriangles_ < nTriangles_);
            Triangle& rTri = pTriangles_[iTriangles_];
            rTri.pMeshTriangle = pMeshTri;
            for (int j=0; j<3; j++) {
                Vertex& v = pMeshTri->pMesh->pVertices[pMeshTri->indices[j]];
                rTri.p[j].x = v.pos.x;
                rTri.p[j].y = v.pos.y;
                rTri.p[j].z = v.pos.z;
            }
            nMeshTris++;
            iTriangles_++;
        } else {
            assert(IsOtherPrimType(ppShapes[i]->GetType()));
            assert(iOtherPrims_ < nOtherPrims_);
            ppOtherPrims_[iOtherPrims_] = ppShapes[i];
            nOtherPrims++;
            iOtherPrims_++;
        }
    }
    return nMeshTris;
}

// 全Shapeを含むBBoxを返す
BBox QBVHBase::SurroundBBox(const IShape** ppShapes, int nShapes)
{
    if (nShapes == 0) {
        return BBox();
    }
    
    BBox bbox = ppShapes[0]->BoundingBox();
    for (int i=1; i<nShapes; i++) {
        bbox = BBox::Surround(bbox, ppShapes[i]->BoundingBox());
    }
    return bbox;
}

// 新しい葉を取得
QBVHBase::Leaf& QBVHBase::GetNewLeaf()
{
    // ensure capacity
    if (iLeaves_ >= nLeafCapacity_) {
        nLeafCapacity_ *= 2;
        Leaf* pNewLeaves = new Leaf[nLeafCapacity_];
        for (int i=0; i<iLeaves_; i++) {
            pNewLeaves[i] = pLeaves_[i];
        }
        delete pLeaves_;
        pLeaves_ = pNewLeaves;
    }
    Leaf& rLeaf = pLeaves_[iLeaves_];
    iLeaves_++;
    return rLeaf;
}

// ツリー全体のBBを返す
BBox QBVHBase::BoundingBox() const
{
    return bbox_;
}

//　葉の交差判定
bool QBVHBase::IntersectLeaf(Leaf& leaf, const Ray& r, float tmin, HitRecord& rec) const
{
    bool ret = false;
    Triangle* pTri = &pTriangles_[leaf.iTriangles];
    int nTri = leaf.nTriangles;
    for (int i=0; i<nTri; i++) {
        bool hit = IntersectTriangle(pTri[i], r, tmin, rec.t, rec);
        ret = ret || hit;
    }
    
    const IShape** ppOP = &ppOtherPrims_[leaf.iOtherPrims];
    int nOP = leaf.nOtherPrims;
    for (int i=0; i<nOP; i++) {
        bool hit = ppOP[i]->Intersect(r, tmin, rec.t, rec);
        ret = ret || hit;
    }
    
    return ret;
}

// 三角形交差判定
bool QBVHBase::IntersectTriangle(Triangle& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const
{
    // based PHISICALLY BASED RENDERING 2ND EDITION, 3.6.2
    Vec3 e1 = tri.p[1] - tri.p[0];
    Vec3 e2 = tri.p[2] - tri.p[0];
    Vec3 s1 = r.d % e2;
    real divisor = s1.dot(e1);
    
    if (divisor == 0.0f)
        return false;
    
    real invDivisor = 1.0f / divisor;
    
    // 重心座標を算出して範囲内にあるかチェック
    Vec3 d = r.o - tri.p[0];
    real b1 = d.dot(s1) * invDivisor;
    if (b1 < 0.0f || b1 > 1.0f)
        return false;
    
    Vec3 s2 = d % e1;
    real b2 = r.d.dot(s2) * invDivisor;
    if (b2 < 0.0f || b1 + b2 > 1.0f)
        return false;
    
    // t算出
    real t = e2.dot(s2) * invDivisor;
    
    // 光線始点より後ろならヒットしない
    if (t < tmin || t > tmax)
        return false;
    
    rec.t = t;
    real b0 = 1.f - b1 - b2;
    const MeshTriangle* pMt = tri.pMeshTriangle;
    Mesh* pMesh = pMt->pMesh;
    if (pMesh->GetUseFaceNormal()) {
        rec.normal = pMt->normal; // face normal
    } else {
        rec.normal = ((pMesh->pVertices[pMt->indices[0]].normal * b0)
                   +  (pMesh->pVertices[pMt->indices[1]].normal * b1)
                   +  (pMesh->pVertices[pMt->indices[2]].normal * b2)).normalize();
    }
    
    if (pMesh->colorUnit_ == CU_Mesh) {
        rec.color = pMesh->color_;
    } else if (pMesh->colorUnit_ == CU_Face) {
        rec.color = pMt->color_;
    } else {
        rec.color = ((pMesh->pVertices[pMt->indices[0]].color * b0)
                  +  (pMesh->pVertices[pMt->indices[1]].color * b1)
                  +  (pMesh->pVertices[pMt->indices[2]].color * b2)) / 3;
    }
    rec.refl = pMesh->material_;
    return true;
}

// 葉に対するレイキャスト
int QBVHBase::RayCastLeaf(vector<HitRecord>& hits, int nHits, Leaf& leaf, const Ray &r, float tmin, float tmax) const
{
    Triangle* pTri = &pTriangles_[leaf.iTriangles];
    int nTri = leaf.nTriangles;
    for (int i=0; i<nTri; i++) {
        if (hits.size() == nHits)
        {
            hits.resize((nHits+1) * 2);
        }
        HitRecord& rec = hits.at(nHits);
        if (IntersectTriangle(pTri[i], r, tmin, tmax, rec)) {
            nHits++;
        }
    }
 
    const IShape** ppOP = &ppOtherPrims_[leaf.iOtherPrims];
    int nOP = leaf.nOtherPrims;
    for (int i=0; i<nOP; i++) {
        nHits = ppOP[i]->RayCast(hits, nHits, r, tmin, tmax);
    }
    
    return nHits;
}

// シェイプをpivotの左右に分ける
int QBVHBase::QSplit(const IShape** pShapes, int nShapes, float pivot, int axis)
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

// BBoxの一番大きい軸
int QBVHBase::LargestAxis(const BBox& bbox)
{
    Vec3 bboxSize = bbox.pp[1] - bbox.pp[0];
    int axis;
    if (bboxSize.x > bboxSize.y) {
        axis = (bboxSize.x > bboxSize.z) ? 0 : 2;
    } else {
        axis = (bboxSize.y > bboxSize.z) ? 1 : 2;
    }
    return axis;
}
