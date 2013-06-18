#include "Scene.h"
#include "SIMD_QBVH.h"
#include <cassert>
#include "simd.h"
#include "Ray.h"

// @todo Branch -> Node に rename?

using namespace std;

#define INITIAL_NODE_NUM 10000 // about 1.2MB
#define INITIAL_LEAF_NUM 30000 // about 0.72MB

SIMD_QBVH::SIMD_QBVH()
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
    pNodes_ = new SIMD_QBVH_NODE[nNodeCapacity_];
    nLeafCapacity_ = INITIAL_LEAF_NUM;
    pLeaves_ = new Leaf[nLeafCapacity_];
}

SIMD_QBVH::~SIMD_QBVH()
{
    delete pNodes_;
}

ShapeType SIMD_QBVH::GetType() const
{
    return ST_QBVH_SIMD;
}

// Build時に使うバッファを初期化
void SIMD_QBVH::Reset()
{
    for (int i=0; i<nNodeCapacity_; i++) {
        pNodes_[i].Reset();
    }
    delete [] pTriangles_;
    delete [] ppOtherPrims_;
    pTriangles_ = NULL;
    ppOtherPrims_ = NULL;
}

// シェイプが葉になるタイプならtrueを返す
bool SIMD_QBVH::IsLeafType(ShapeType shapeType)
{
    return (shapeType == ST_MESH_TRIANGLE
         || shapeType == ST_SPHERE
         || shapeType == ST_TRIANGLE);
}

// シェイプがMeshTriangle以外の葉になるタイプならtrueを返す
bool SIMD_QBVH::IsOtherPrimType(ShapeType shapeType)
{
    return (shapeType == ST_SPHERE
         || shapeType == ST_TRIANGLE);
}

// 葉ノード総数をカウント
// MESH, QBVHは葉にならない
int SIMD_QBVH::CountLeafShapes(const IShape** pShapes, int nShapes)
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
const IShape** SIMD_QBVH::FlattenLeafShapes(const IShape** ppFlatten, const IShape** ppShapes, int nShapes)
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
void SIMD_QBVH::Build(const IShape** pShapes, int nShapes)
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
        else if (IsOtherPrimType(st)) {
            nOtherPrims_++;
        }
        else
        {
            assert(false);
        }
    }
    
    // MeshTriangleとその他のシェイプの領域を確保
    pTriangles_ = new SIMD_TRIANGLE[nTriangles_];
    ppOtherPrims_ = new const IShape*[nOtherPrims_];
    printf("LeafShape=%d MeshTriangle=%d OtherPrim=%d\n", nLeafShapes, nTriangles_, nOtherPrims_);
    
    // 構築
    SIMD_QBVH_NODE& rNode = pNodes_[0];
    iNodes_++;
    BuildBranch(rNode, ppLeafShapes, nLeafShapes);
    
    // 全体のBBoxを作成
    bbox_ = SurroundBBox(rNode.bboxes);
    
    // realloc
    printf("Node=%d Leaf=%d\n", iNodes_, iLeaves_ );
    SIMD_QBVH_NODE* pNewNodes = new SIMD_QBVH_NODE[iNodes_];
    for (int i=0; i<iNodes_; i++) {
        pNewNodes[i] = pNodes_[i];
    }
    delete pNodes_;
    Leaf* pNewLeaves = new Leaf[iLeaves_];
    for (int i=0; i<iLeaves_; i++) {
        pNewLeaves[i] = pLeaves_[i];
    }
    delete pLeaves_;
    pNodes_ = pNewNodes;
    pLeaves_ = pNewLeaves;
}

// 枝ノード構築
void SIMD_QBVH::BuildBranch(SIMD_QBVH_NODE& rNode, const IShape** pShapes, int nShapes)
{
    depth_++;
    
    // qsplitのピボットとして使用するためにBBoxの中間点を求める
    BBox bbox = SurroundBBox(pShapes, nShapes); // @todo 上から渡した方がいい
    
    // BBoxの一番大きい軸
    int axis = LargestAxis(bbox);
    
    //printf("axes=%d %d %d\n", axes[0], axes[1], axes[2]);
    
    // 空間的な中心をpivotとする
    Vec3 pivot = (bbox.Max() + bbox.Min()) * 0.5f;
    
    // pShapesの要素をpivotの左右に振り分ける
    int mid_point = QSplit(pShapes, nShapes, pivot.e[axis], axis);
    
    // leftとrightそれぞれのbboxを求める
    BBox bbox_l = SurroundBBox(pShapes, mid_point);
    BBox bbox_r = SurroundBBox(&pShapes[mid_point], nShapes-mid_point);
    
    // 左右の分割軸
    int axis_l = LargestAxis(bbox_l);
    int axis_r = LargestAxis(bbox_r);
    
    // 交差判定では使ってないけど一応保存しておく
    rNode.axis0 = axis;
    rNode.axis1 = axis_l;
    rNode.axis2 = axis_r;
    
    // leftとrightそれぞれpivot決める
    float pivot_l = ((bbox_l.Max() + bbox_l.Min()) * 0.5f).e[axis_l];
    float pivot_r = ((bbox_r.Max() + bbox_r.Min()) * 0.5f).e[axis_r];
    
    // leftとrightをそれぞれさらに振り分ける
    int mid_point_l = QSplit(pShapes, mid_point, pivot_l, axis_l);
    int mid_point_r = QSplit(&pShapes[mid_point], nShapes-mid_point, pivot_r, axis_r);
    
    // 子ノードを構築
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
    
    BBox childBBoxes[4];
    for (int i=0; i<4; i++) {
        
        // bbox求める
        childBBoxes[i] = SurroundBBox(ppChildShapes[i], nChildren[i]);
        
        // 子が空なら特別なindexをセット
        if (nChildren[i] == 0) {
            //printf("-\n");
            rNode.children[i] = INT_MIN;
        }
        // 子が4より多ければさらに枝に分割
        else if (nChildren[i] > 4) {
            //printf("=\n");
            //assert(iNodes_ < nNodeCapacity_);
            rNode.children[i] = iNodes_;
            SIMD_QBVH_NODE& rChildNode = GetNewNode();
            BuildBranch(rChildNode, ppChildShapes[i], nChildren[i]);
        }
        // 子が4以下だったら葉を作る
        else {
            //printf("*\n");
            rNode.children[i] = INT_MIN; // 符号で枝と葉を区別
            assert(iLeaves_ < (1<<30));
            rNode.children[i] |= iLeaves_;
            Leaf& rLeaf = GetNewLeaf();
            
            rLeaf.iTriangles = iTriangles_;
            rLeaf.iOtherPrims = iOtherPrims_;
            BuildLeaf(rLeaf.nTriangles, rLeaf.nOtherPrims, ppChildShapes[i], nChildren[i]);
        }
    }
    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            rNode.bboxes[i][j] = _mm_set_ps(
                childBBoxes[3].pp[i].e[j],
                childBBoxes[2].pp[i].e[j],
                childBBoxes[1].pp[i].e[j],
                childBBoxes[0].pp[i].e[j]
                                            );
        }
    }
    
    depth_--;
    //printf("rNode=%p iTri=%d iOP=%d\n", &rNode, iTriangles_, iOtherPrims_);
}

// 葉ノード構築
int SIMD_QBVH::BuildLeaf(u8& nMeshTris, u8& nOtherPrims, const IShape** ppShapes, int nShapes)
{
    nMeshTris = 0;
    nOtherPrims = 0;
    for (int i=0; i<nShapes; i++) {
        if (ppShapes[i]->GetType() == ST_MESH_TRIANGLE) {
            const MeshTriangle* pMeshTri = reinterpret_cast<const MeshTriangle*>(ppShapes[i]);
            assert(iTriangles_ < nTriangles_);
            SIMD_TRIANGLE& rTri = pTriangles_[iTriangles_];
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

// 全Shapeを内包するBBoxを返す
BBox SIMD_QBVH::SurroundBBox(const IShape** ppShapes, int nShapes)
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

// 4つのBBoxを内包するBBoxを返す
BBox SIMD_QBVH::SurroundBBox(__m128 bboxes[2][3])
{
    //Um128 bbox_min = bboxes[0][0];
    //Um128 bbox_max = bboxes[1][0];
    Vec3 minv;
    Vec3 maxv;
    for (int xyz=0; xyz<3; xyz++) {
        Um128 bbox_min = bboxes[0][xyz];
        Um128 bbox_max = bboxes[1][xyz];
        minv.e[xyz] = bbox_min.e[0]; // 0:最初のBox
        maxv.e[xyz] = bbox_max.e[0];
        // 残り3つのBoxをループ
        for (int box=1; box<4; box++) {
            if (bbox_min.e[box] < minv.e[xyz]) {
                minv.e[xyz] = bbox_min.e[box];
            }
            if (bbox_max.e[box] > maxv.e[xyz]) {
                maxv.e[xyz] = bbox_max.e[box];
            }
        }
    }
    return BBox(minv, maxv);
}

// 新しい枝を取得
SIMD_QBVH::SIMD_QBVH_NODE& SIMD_QBVH::GetNewNode()
{
    // ensure capacity
    if (iNodes_ >= nNodeCapacity_) {
        nNodeCapacity_ *= 2;
        SIMD_QBVH_NODE* pNewNodes = new SIMD_QBVH_NODE[nNodeCapacity_];
        for (int i=0; i<iNodes_; i++) {
            pNewNodes[i] = pNodes_[i];
        }
        delete pNodes_;
        pNodes_ = pNewNodes;
    }
    SIMD_QBVH_NODE& rNode = pNodes_[iNodes_];
    iNodes_++;
    return rNode;
}

// 新しい葉を取得
SIMD_QBVH::Leaf& SIMD_QBVH::GetNewLeaf()
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
BBox SIMD_QBVH::BoundingBox() const
{
    return bbox_;
}

// SIMD命令を使って4つのBBoxを一度に交差判定
// reference:
//   - Kaze Renderer (@TODO URL here)
//   - Shirley, Realistic Ray Tracing
//   - Implemented sisd version in BBox::RayIntersect()
inline int SIMD_QBVH::IntersectSIMD(
    const __m128 bboxes[2][3],  // min-max[2] of xyz[3] of boxes[4]
    const __m128 orig[3],       // ray origin, xyz[3]
    const __m128 idir[3],       // ray inverse direction, xyz[3]
    const int sign[3],          // ray xyz direction -> +:0,-:1
    __m128 tmin, __m128 tmax    // ray range tmin-tmax
) const
{
    // x coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[0]][0], orig[0]), idir[0])
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[0]][0], orig[0]), idir[0])
    );

    // y coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[1]][1], orig[1]), idir[1])
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[1]][1], orig[1]), idir[1])
    );

    // z coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[2]][2], orig[2]), idir[2])
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[2]][2], orig[2]), idir[2])
    );

    return _mm_movemask_ps(_mm_cmpge_ps(tmax, tmin));//tmin<tmaxとなれば交差
}


// 枝との交差判定
bool SIMD_QBVH::IntersectBranch(SIMD_QBVH_NODE& rNode, const Ray &r, float tmin, HitRecord& rec) const
{
    depth_++;
    
    bool ret = false;
    __m128 xmtmin = _mm_set1_ps(tmin);
    __m128 xmtmax = _mm_set1_ps(rec.t);
    __m128 orig[3] = { _mm_set1_ps(r.o.x), _mm_set1_ps(r.o.y), _mm_set1_ps(r.o.z) };
    __m128 idir[3] = { _mm_set1_ps(r.idir.x), _mm_set1_ps(r.idir.y), _mm_set1_ps(r.idir.z) };
    
    int hit = IntersectSIMD(rNode.bboxes, orig, idir, r.sign, xmtmin, xmtmax);
    
    /*
    Um128 dbg_minx = rNode.bboxes[0][0];
    Um128 dbg_miny = rNode.bboxes[0][1];
    Um128 dbg_minz = rNode.bboxes[0][2];
    Um128 dbg_maxx = rNode.bboxes[1][0];
    Um128 dbg_maxy = rNode.bboxes[1][1];
    Um128 dbg_maxz = rNode.bboxes[1][2];
    */
    for (int i=0; i<4; i++) {
        /*
        printf("%d  %f %f %f  %f %f %f\n", depth_,
            dbg_minx.e[i], dbg_miny.e[i], dbg_minz.e[i],
            dbg_maxx.e[i], dbg_maxy.e[i], dbg_maxz.e[i]);
        */
        
        if (((1 << i) & hit) && (rNode.children[i] != INT_MIN)) {
            //printf("hit!\n");
            // 枝
            if ((rNode.children[i] & (0x01 << 31)) == 0) {
                bool hit = IntersectBranch(pNodes_[rNode.children[i]], r, tmin, rec);
                ret = ret || hit;
            }
            // 葉
            else {
                //printf("Leaf!\n");
                int iLeaves = (rNode.children[i] & 0x7FFFFFFF);
                Leaf& rLeaf = pLeaves_[iLeaves];
                assert((rLeaf.nTriangles + rLeaf.nOtherPrims) > 0);
                bool hit = IntersectLeaf(rLeaf, r, tmin, rec);
                ret = ret || hit;
            }
        }
    }

    depth_--;
    return ret;
}

//　葉の交差判定
bool SIMD_QBVH::IntersectLeaf(Leaf& leaf, const Ray& r, float tmin, HitRecord& rec) const
{
    bool ret = false;
    SIMD_TRIANGLE* pTri = &pTriangles_[leaf.iTriangles];
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
bool SIMD_QBVH::IntersectTriangle(SIMD_TRIANGLE& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const
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

// 交差判定
bool SIMD_QBVH::Intersect(const Ray &r, float tmin, float tmax, HitRecord& rec) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return false;
    
    rec.t = tmax;
    
    return IntersectBranch(pNodes_[0], r, tmin, rec);
}

// 枝とのレイキャスト
int SIMD_QBVH::RayCastBranch(vector<HitRecord>& hits, int nHits, SIMD_QBVH_NODE& rNode, const Ray &r, float tmin, float tmax) const
{
    __m128 xmtmin = _mm_set1_ps(tmin);
    __m128 xmtmax = _mm_set1_ps(tmax);
    __m128 orig[3] = { _mm_set1_ps(r.o.x), _mm_set1_ps(r.o.y), _mm_set1_ps(r.o.z) };
    __m128 idir[3] = { _mm_set1_ps(r.idir.x), _mm_set1_ps(r.idir.y), _mm_set1_ps(r.idir.z) };
    
    int hit = IntersectSIMD(rNode.bboxes, orig, idir, r.sign, xmtmin, xmtmax);
    
    for (int i=0; i<4; i++) {
        if ((1 << i) & hit) {
            // 枝
            if ((rNode.children[i] & (0x01 << 31)) == 0) {
                nHits = RayCastBranch(hits, nHits, pNodes_[rNode.children[i]], r, tmin, tmax);
            }
            // 葉
            else {
                int iLeaves = (rNode.children[i] & 0x7FFFFFFF);
                Leaf& rLeaf = pLeaves_[iLeaves];
                assert((rLeaf.nTriangles + rLeaf.nOtherPrims) > 0);
                nHits = RayCastLeaf(hits, nHits, rLeaf, r, tmin, tmax);
            }
        }
    }
    
    return nHits;
}

// 葉に対するレイキャスト
int SIMD_QBVH::RayCastLeaf(vector<HitRecord>& hits, int nHits, Leaf& leaf, const Ray &r, float tmin, float tmax) const
{
    SIMD_TRIANGLE* pTri = &pTriangles_[leaf.iTriangles];
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

// レイキャスト
int SIMD_QBVH::RayCast(vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return nHits;
    
    return RayCastBranch(hits, nHits, pNodes_[0], r, tmin, tmax);
}

// シェイプをpivotの左右に分ける
int SIMD_QBVH::QSplit(const IShape** pShapes, int nShapes, float pivot, int axis)
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
int SIMD_QBVH::LargestAxis(const BBox& bbox)
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