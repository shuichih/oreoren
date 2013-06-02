#include "Scene.h"
#include "QBVH.h"
#include <cassert>
#include "simd.h"
#include "Ray.h"

using namespace std;

#define NODE_MAX 100000 // about 12MB

QBVH::QBVH()
: pTriangles_(NULL)
, ppOtherPrims_(NULL)
, nNodes_(0)
, nTriangles_(0)
, nOtherPrims_(0)
, iNodes_(0)
, iTriangles_(0)
, iOtherPrims_(0)
{
    /*
    for (int mm=0; mm<2; mm++) {
        for (int axis=0; axis<3; axis++) {
            for (int child=0; child<4; child++) {
                bboxes[mm][axis][child] = 0;
            }
        }
    }
    */
    
    nNodes_ = NODE_MAX;
    pNodes_ = new SISD_QBVH_NODE[nNodes_]; // about 12MB
}

QBVH::~QBVH()
{
    delete pNodes_;
}

ShapeType QBVH::GetType() const
{
    return ST_QBVH_SISD;
}

void QBVH::Reset()
{
    for (int i=0; i<nNodes_; i++) {
        pNodes_[i].Reset();
    }
    delete [] pTriangles_;
    delete [] ppOtherPrims_;
    pTriangles_ = NULL;
    ppOtherPrims_ = NULL;
}

bool QBVH::IsLeafType(ShapeType shapeType)
{
    return (shapeType == ST_MESH_TRIANGLE
         || shapeType == ST_SPHERE
         || shapeType == ST_TRIANGLE);
}

// 葉ノード総数をカウント
// MESH, QBVHは葉にならない
int QBVH::CountLeafShapes(const IShape** pShapes, int nShapes)
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
const IShape** QBVH::FlattenLeafShapes(const IShape** ppFlatten, const IShape** ppShapes, int nShapes)
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

// 構築
void QBVH::Build(const IShape** pShapes, int nShapes)
{
    Reset();
    
    // 全ての葉のシェイプをフラットな配列に
    int nLeafShapes = CountLeafShapes(pShapes, nShapes);
    const IShape** ppLeafShapes = new const IShape*[nLeafShapes];
    FlattenLeafShapes(ppLeafShapes, pShapes, nShapes);
    
    // MeshTriangleとその他のシェイプの数をカウント
    nTriangles_ = 0;
    for (int i=0; i<nLeafShapes; i++) {
        ShapeType st = ppLeafShapes[i]->GetType();
        if (st == ST_MESH_TRIANGLE) {
            nTriangles_++;
        }
        else if (st == ST_TRIANGLE || st == ST_SPHERE) {
            nOtherPrims_++;
        }
        else
        {
            assert(false);
        }
    }
    //nTriangles_*=2;
    
    // MeshTriangleとその他のシェイプの領域を確保
    pTriangles_ = new SISD_TRIANGLE[nTriangles_];
    ppOtherPrims_ = new const IShape*[nOtherPrims_];
    printf("LeafShape=%d MeshTriangle=%d OtherPrim=%d\n", nLeafShapes, nTriangles_, nOtherPrims_);
    
    // 枝構築
    SISD_QBVH_NODE& rNode = pNodes_[0];
    iNodes_++;
    BuildBranch(rNode, ppLeafShapes, nLeafShapes);
    
    // 全体のBBoxを作成
    BBox bbox0 = BBox::Surround(rNode.bboxes[0], rNode.bboxes[1]);
    BBox bbox1 = BBox::Surround(rNode.bboxes[2], rNode.bboxes[3]);
    bbox_ = BBox::Surround(bbox0, bbox1);
}

// 枝ノード構築
void QBVH::BuildBranch(SISD_QBVH_NODE& rNode, const IShape** pShapes, int nShapes)
{
    // qsplitのピボットとして使用するためにBBoxの中間点を求める
    BBox bbox = SurroundBBox(pShapes, nShapes);
    
    // BBoxの幅が大きい順に軸をソート
    int axes[] = { 0, 1, 2 };
    Vec3 bboxSize = bbox.pp[1] - bbox.pp[0];
    if (bboxSize.z > bboxSize.y) {
        Swap(bboxSize.z, bboxSize.y);
        Swap(axes[2], axes[1]);
        if (bboxSize.y > bboxSize.x) {
            Swap(bboxSize.y, bboxSize.x);
            Swap(axes[1], axes[0]);
        }
        if (bboxSize.z > bboxSize.y) {
            Swap(bboxSize.z, bboxSize.y);
            Swap(axes[2], axes[1]);
        }
    }
    
    // 分割軸。交差判定で使ってないけど一応入れとく
    rNode.axes = (axes[2] << 16) | (axes[1] << 8) | axes[0];
    
    //printf("axes=%d %d %d\n", axes[0], axes[1], axes[2]);
    
    // 空間的な中心をpivotとする
    Vec3 pivot = (bbox.Max() + bbox.Min()) * 0.5f;
    
    // pShapesの要素をpivotの左右に振り分ける
    int mid_point = QSplit(pShapes, nShapes, pivot.e[axes[0]], axes[0]);
    
    // leftとrightそれぞれのbboxを求める
    BBox bbox_l = SurroundBBox(pShapes, mid_point);
    BBox bbox_r = SurroundBBox(&pShapes[mid_point], nShapes-mid_point);
    
    // leftとrightそれぞれpivot決める
    // @todo axisは子のBBで一番大きい軸を使う
    float pivot_l = ((bbox_l.Max() + bbox_l.Min()) * 0.5f).e[axes[1]];
    float pivot_r = ((bbox_r.Max() + bbox_r.Min()) * 0.5f).e[axes[1]];
    
    // leftとrightをそれぞれさらにQSplit
    int mid_point_l = QSplit(pShapes, mid_point, pivot_l, axes[1]);
    int mid_point_r = QSplit(&pShapes[mid_point], nShapes-mid_point, pivot_r, axes[1]);
    
    // 各子ノードを構築
    const IShape** ppChildShapes[4] = {
        &pShapes[0],
        &pShapes[mid_point_l],
        &pShapes[mid_point],
        &pShapes[mid_point+mid_point_r]
    };
    int nChildren[4] = {
        mid_point_l,
        mid_point-mid_point_l,
        mid_point_r,
        nShapes-(mid_point+mid_point_r)
    };
    
    //printf("%d %d %d %d  %d %d %d %d  %d\n", 0, mid_point_l, mid_point, mid_point+mid_point_r,
    //       mid_point_l, mid_point-mid_point_l, mid_point_r, nShapes-(mid_point + mid_point_r), nShapes);
    
    rNode.iOtherPrims = iOtherPrims_;
    //printf("rNode=%p iOP=%d\n", &rNode, iOtherPrims_);
    
    for (int i=0; i<4; i++) {
        
        // bbox求める
        rNode.bboxes[i] = SurroundBBox(ppChildShapes[i], nChildren[i]);
        
        // 子が空なら特別なindexをセット
        if (nChildren[i] == 0) {
            rNode.children[i] = INT_MIN;
        }
        // 子が4より多ければさらに枝に分割
        else if (nChildren[i] > 4) {
            assert(iNodes_ < nNodes_);
            
            SISD_QBVH_NODE& rChildNode = pNodes_[iNodes_];
            rNode.children[i] = iNodes_;
            iNodes_++;
            BuildBranch(rChildNode, ppChildShapes[i], nChildren[i]);
        }
        // 子が4以下だったらtriangle listを作る
        else {
            int nMeshTris = 0;
            int nOtherPrims = 0;
            rNode.children[i] = INT_MIN; // 符号で枝と葉を区別
            assert(iTriangles_ < (1<<27));
            rNode.children[i] |= iTriangles_; // 三角形配列へのインデックス: 27bit
            
            BuildLeaf(nMeshTris, nOtherPrims, ppChildShapes[i], nChildren[i]);
            
            if (nMeshTris > 0) {
                rNode.children[i] |= (nMeshTris << 27); // 子の数: 4bit
            } else {
                rNode.children[i] = INT_MIN;
            }
            rNode.nOtherPrims += nOtherPrims;
            //printf("b %x\n", rNode.children[i]);
            //printf("nMT=%d\n", nMeshTris);
        }
    }
    //printf("rNode=%p iTri=%d nOP=%d\n", &rNode, iTriangles_, rNode.nOtherPrims);
}

// 葉ノード構築
void QBVH::BuildLeaf(int& nMeshTris, int& nOtherPrims, const IShape** ppShapes, int nShapes)
{
    for (int i=0; i<nShapes; i++) {
        if (ppShapes[i]->GetType() == ST_MESH_TRIANGLE) {
            const MeshTriangle* pMeshTri = reinterpret_cast<const MeshTriangle*>(ppShapes[i]);
            assert(iTriangles_ < nTriangles_);
            SISD_TRIANGLE& rTri = pTriangles_[iTriangles_];
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
            //printf("single prim %p\n", ppShapes[i]);
            assert(ppShapes[i]->GetType() == ST_TRIANGLE || ppShapes[i]->GetType() == ST_SPHERE);
            assert(iOtherPrims_ < nOtherPrims_);
            ppOtherPrims_[iOtherPrims_] = ppShapes[i];
            printf("ppOP[%d]=%p\n", iOtherPrims_, ppOtherPrims_[iOtherPrims_]);
            nOtherPrims++;
            iOtherPrims_++;
        }
    }
}

// 全Shapeを含むBBoxを返す
BBox QBVH::SurroundBBox(const IShape** ppShapes, int nShapes)
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

BBox QBVH::BoundingBox() const
{
    return bbox_;
}

#if 0
inline int QBVH::IntersectSIMD(
    const __m128 bboxes[2][3],  // 4 boxes : min-max[2] of xyz[3] of boxes[4]
    const __m128 org[3],        // ray origin
    const __m128 idir[3],       // ray inveresed direction
    const int sign[3],          // ray xyz direction -> +:0,-:1
    __m128 tmin, __m128 tmax    // ray range tmin-tmax
)
{
    // x coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[0]][0],org[0]), idir[0])
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[0]][0], org[0]), idir[0])
    );

    // y coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[1]][1],org[1]), idir[1])
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[1]][1], org[1]), idir[1])
    );

    // z coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[2]][2],org[2]), idir[2])
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[2]][2], org[2]), idir[2])
    );

    return _mm_movemask_ps(_mm_cmpge_ps(tmax, tmin));//tmin<tmaxとなれば交差
}

int QBVH::IntersectSIMD(
    const __m128 bboxes[2][3],  // 4 boxes : min-max[2] of xyz[3] of boxes[4]
    const __m128 org[3],        // ray origin
    const __m128 idir[3],       // ray inveresed direction
    const int sign[3],          // ray xyz direction -> +:0,-:1
    __m128 tmin, __m128 tmax    // ray range tmin-tmax
)
{
}

#endif

bool QBVH::IntersectBranch(SISD_QBVH_NODE& rNode, const Ray &r, float tmin, HitRecord& rec) const
{
    bool ret = false;
    for (int i=0; i<4; i++) {
        if (rNode.children[i] == INT_MIN && rNode.nOtherPrims == 0) {
            // 空ノード
            continue;
        }
        
        if (rNode.bboxes[i].RayIntersect(r, tmin, rec.t)) {
            // 枝
            if ((rNode.children[i] & (0x01 << 31)) == 0) {
                bool hit = IntersectBranch(pNodes_[rNode.children[i]], r, tmin, rec);
                ret = ret || hit;
            }
            // 葉
            else {
                int nTri = (rNode.children[i] & (0x0F << 27)) >> 27;
                //printf("i %x\n", rNode.children[i]);
                //printf("a %x\n", a);
                int idx = rNode.children[i] & 0x07FFFFFF;
                //printf("nTri=%d idx=%d\n", nTri, idx);
                bool hit = IntersectLeaf(&pTriangles_[idx], nTri, r, tmin, rec);
                ret = ret || hit;
                
                //printf("I rNode=%p\n", &rNode);
                for (int isp=0; isp<rNode.nOtherPrims; isp++) {
                    hit = ppOtherPrims_[rNode.iOtherPrims+isp]->Intersect(r, tmin, rec.t, rec);
                    ret = ret || hit;
                }
            }
        }
    }
    
    return ret;
}

// 三角形交差判定
// @todo Scene::Triangleと共通化
bool QBVH::IntersectTriangle(SISD_TRIANGLE& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const
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
    
//if (tri.p[0].z < 0) return false;
    
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

//　葉の交差判定
bool QBVH::IntersectLeaf(SISD_TRIANGLE* pTri, int nTri, const Ray& r, float tmin, HitRecord& rec) const
{
    bool ret = false;
    for (int i=0; i<nTri; i++) {
        bool hit = IntersectTriangle(pTri[i], r, tmin, rec.t, rec);
        ret = ret || hit;
    }
    //if (ret) {
    //    printf("I\n");
    //}
    return ret;
}

// 交差判定
bool QBVH::Intersect(const Ray &r, float tmin, float tmax, HitRecord& rec) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return false;
    
    rec.t = tmax;
    
    return IntersectBranch(pNodes_[0], r, tmin, rec);
}

// 引数整理
int QBVH::RayCastBranch(vector<HitRecord>& hits, int nHits, SISD_QBVH_NODE& rNode, const Ray &r, float tmin, float tmax) const
{
    for (int i=0; i<4; i++) {
        if (rNode.children[i] == INT_MIN) {
            // 空ノード
            continue;
        }
        
        if (rNode.bboxes[i].RayIntersect(r, tmin, tmax)) {
            // 枝
            if ((rNode.children[i] & (0x01 << 31)) == 0) {
                return RayCastBranch(hits, nHits, pNodes_[rNode.children[i]], r, tmin, tmax);
            }
            // 葉
            else {
                int nTri = (rNode.children[i] & (0x0F << 27)) >> 27;
                return RayCastLeaf(hits, nHits, &pTriangles_[i], nTri, r, tmin, tmax);
            }
        }
    }
    
    return nHits;
}

// 葉に対するレイキャスト
int QBVH::RayCastLeaf(vector<HitRecord>& hits, int nHits, SISD_TRIANGLE* pTri, int nTri, const Ray &r, float tmin, float tmax) const
{
    if (hits.size() == nHits)
    {
        hits.resize((nHits+1) * 2);
    }
    
    HitRecord& rec = hits.at(nHits);
    for (int i=0; i<nTri; i++) {
        if (IntersectTriangle(pTri[i], r, tmin, tmax, rec)) {
            nHits++;
        }
    }
    
    return nHits;
}

// レイキャスト
int QBVH::RayCast(vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return nHits;
    
    return RayCastBranch(hits, nHits, pNodes_[0], r, tmin, tmax);
}

#if 0
void QBVH::LimitMinScale(float minScale)
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
    
    /*
    if (pLeft_->IsBVH()) ((QBVH*)pLeft_)->LimitMinScale(minScale);
    if (pRight_->IsBVH()) ((QBVH*)pRight_)->LimitMinScale(minScale);
    */
    // @todo implement
}
#endif

int QBVH::QSplit(const IShape** pShapes, int nShapes, float pivot, int axis)
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


