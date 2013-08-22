#include "Scene.h"
#include "QBVHBase.h"
#include <cassert>
#include "simd.h"
#include "Ray.h"
#include "Material.h"

using namespace std;

#define INITIAL_NODE_NUM 10000 // about 1.2MB
#define INITIAL_LEAF_NUM 30000 // about 0.72MB

template <typename NODE_T>
QBVHBase<NODE_T>::QBVHBase()
: ShapeBase(NULL)
, pTriangles_(NULL)
, ppOtherPrims_(NULL)
, pNodes_(NULL)
, pLeaves_(NULL)
{
    Reset();
    
    nNodeCapacity_ = INITIAL_NODE_NUM;
    nLeafCapacity_ = INITIAL_LEAF_NUM;
    pLeaves_ = new Leaf[nLeafCapacity_];
}

template <typename NODE_T>
QBVHBase<NODE_T>::~QBVHBase()
{
    Reset();
}

// Build時に使うバッファを初期化
template <typename NODE_T>
void QBVHBase<NODE_T>::Reset()
{
    delete [] pTriangles_;
    delete [] ppOtherPrims_;
    pTriangles_ = NULL;
    ppOtherPrims_ = NULL;
    nTriangles_ = 0;
    nOtherPrims_ = 0;
    iNodes_ = 0;
    iTriangles_ = 0;
    iLeaves_ = 0;
    iOtherPrims_ = 0;
    depth_ = 0;
    
    delete pNodes_;
    delete pLeaves_;
    nNodeCapacity_ = INITIAL_NODE_NUM;
    pNodes_ = new NODE_T[nNodeCapacity_];
    nLeafCapacity_ = INITIAL_LEAF_NUM;
    pLeaves_ = new Leaf[nLeafCapacity_];
}

// シェイプが葉になるタイプならtrueを返す
template <typename NODE_T>
bool QBVHBase<NODE_T>::IsLeafType(ShapeType shapeType)
{
    return (shapeType == ST_MESH_TRIANGLE
         || shapeType == ST_SPHERE
         || shapeType == ST_TRIANGLE);
}

// シェイプがMeshTriangle以外の葉になるタイプならtrueを返す
template <typename NODE_T>
bool QBVHBase<NODE_T>::IsOtherPrimType(ShapeType shapeType)
{
    return (shapeType == ST_SPHERE
         || shapeType == ST_TRIANGLE);
}

template <typename NODE_T>
bool QBVHBase<NODE_T>::IsMeshTriangleType(ShapeType shapeType)
{
    return shapeType == ST_MESH_TRIANGLE;
}

// BBoxの一番大きい軸
template <typename NODE_T>
int QBVHBase<NODE_T>::LargestAxis(const BBox& bbox)
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

//--------------------------------------------------------------------------------------------------
// Build
//--------------------------------------------------------------------------------------------------

// 葉ノード総数をカウント
// MESH, QBVHは葉にならない
template <typename NODE_T>
int QBVHBase<NODE_T>::CountLeafShapes(const IShape** pShapes, int nShapes)
{
    int ret = 0;
    for (int i=0; i<nShapes; i++) {
        int nChildren = pShapes[i]->GetChildNum();
        if (nChildren == 0) {
            ShapeType st = pShapes[i]->GetType();
            if (IsLeafType(st)) {
                ret++;
            }
            if (IsOtherPrimType(st)) {
                nOtherPrims_++;
            }
            if (IsMeshTriangleType(st)) {
                nTriangles_++;
            }
        } else {
            ret += CountLeafShapes(pShapes[i]->GetChildren(), nChildren);
        }
    }
    return ret;
}

// 葉のシェイプだけをフラットな配列にする
template <typename NODE_T>
const IShape** QBVHBase<NODE_T>::FlattenLeafShapes(const IShape** ppFlatten, const IShape** ppShapes, int nShapes)
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

// 子を空にセット
template <typename NODE_T>
void QBVHBase<NODE_T>::SetChildEmpty(int& child)
{
    child = INT_MIN;
}

// 子に枝インデックスをセット
template <typename NODE_T>
void QBVHBase<NODE_T>::SetChildNode(int& child, int index)
{
    assert(index < 0x80000000);
    
    // index == 0を扱えるようにbit not
    child = ~index & 0x7FFFFFFF;
}

// 子に葉インデックスをセット
template <typename NODE_T>
void QBVHBase<NODE_T>::SetChildLeaf(int& child, int index)
{
    assert(index < 0x80000000);
    
    // index == 0を扱えるようにbit not
    // 最上位ビットで枝と葉を区別
    child = 0x80000000 | (~index & 0x7FFFFFFF);
}

// 子が空ならtrue
template <typename NODE_T>
bool QBVHBase<NODE_T>::IsChildEmpty(int child) const
{
    return child == INT_MIN;
}


// 子のインデックスを取得
template <typename NODE_T>
int QBVHBase<NODE_T>::GetChildIndex(int child) const
{
    // index == 0を扱えるようにbit not
    return ~child & 0x7FFFFFFF;
}

// 子が枝ならtrue
template <typename NODE_T>
bool QBVHBase<NODE_T>::IsChildNode(int child) const
{
    // 最上位ビットで枝と葉を区別
    return (0x80000000 & child) == 0;
}

// 構築
template <typename NODE_T>
void QBVHBase<NODE_T>::Build(const IShape** pShapes, int nShapes)
{
    Reset();
    
    // 全ての葉のシェイプをフラットな配列に
    int nLeafShapes = CountLeafShapes(pShapes, nShapes);
    const IShape** ppLeafShapes = new const IShape*[nLeafShapes];
    FlattenLeafShapes(ppLeafShapes, pShapes, nShapes);
    
    // MeshTriangleとその他のシェイプの領域を確保
    pTriangles_ = new Triangle[nTriangles_];
    ppOtherPrims_ = new const IShape*[nOtherPrims_];
    printf("LeafShape=%d MeshTriangle=%d OtherPrim=%d\n", nLeafShapes, nTriangles_, nOtherPrims_);
    
    // 構築
    iNodes_++;
    BuildBranch(0, ppLeafShapes, nLeafShapes);
    
    // 全体のBBoxを作成
    MakeWholeBBox();
    
    ShrinkBuffersToFit();
}

// 葉ノード構築
template <typename NODE_T>
int QBVHBase<NODE_T>::BuildLeaf(u8& nMeshTris, u8& nOtherPrims, const IShape** ppShapes, int nShapes)
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
template <typename NODE_T>
BBox QBVHBase<NODE_T>::SurroundBBox(const IShape** ppShapes, int nShapes)
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

// 新しい枝を取得
template <typename NODE_T>
int QBVHBase<NODE_T>::AddNewNode()
{
    // ensure capacity
    if (iNodes_ >= nNodeCapacity_) {
        nNodeCapacity_ *= 2;
        NODE_T* pNewNodes = new NODE_T[nNodeCapacity_];
        for (int i=0; i<iNodes_; i++) {
            pNewNodes[i] = pNodes_[i];
        }
        delete pNodes_;
        pNodes_ = pNewNodes;
    }
    return ++iNodes_;
}

// 新しい葉を取得
template <typename NODE_T>
typename QBVHBase<NODE_T>::Leaf& QBVHBase<NODE_T>::AddNewLeaf()
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

// バッファを使用サイズに合わせてシュリンクする
template <typename NODE_T>
void QBVHBase<NODE_T>::ShrinkBuffersToFit()
{
    // realloc
    printf("Node=%d Leaf=%d\n", iNodes_, iLeaves_);
    NODE_T* pNewNodes = new NODE_T[iNodes_];
    for (int i=0; i<iNodes_; i++) {
        pNewNodes[i] = pNodes_[i];
    }
    delete pNodes_;
    pNodes_ = pNewNodes;
    
    Leaf* pNewLeaves = new Leaf[iLeaves_];
    for (int i=0; i<iLeaves_; i++) {
        pNewLeaves[i] = pLeaves_[i];
    }
    delete pLeaves_;
    
    pLeaves_ = pNewLeaves;
}

// ツリー全体のBBを返す
template <typename NODE_T>
BBox QBVHBase<NODE_T>::BoundingBox() const
{
    return bbox_;
}

//--------------------------------------------------------------------------------------------------
// Intersection
//--------------------------------------------------------------------------------------------------

// 交差判定
template <typename NODE_T>
bool QBVHBase<NODE_T>::Intersect(const Ray &r, float tmin, float tmax, HitRecord& rec) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return false;
    
    rec.t = tmax;
    
    return IntersectBranch(pNodes_[0], r, tmin, rec);
}

//　葉の交差判定
template <typename NODE_T>
bool QBVHBase<NODE_T>::IntersectLeaf(Leaf& leaf, const Ray& r, float tmin, HitRecord& rec) const
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
template <typename NODE_T>
bool QBVHBase<NODE_T>::IntersectTriangle(Triangle& tri, const Ray& r, float tmin, float tmax, HitRecord& rec) const
{
    // @toto Scene　Triangleと共通化
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
        rec.pMaterial = pMesh->GetMaterial();
        rec.color = pMesh->GetMaterial()->color;
    } else if (pMesh->colorUnit_ == CU_Face) {
        rec.pMaterial = pMt->pMaterial;
        rec.color = pMt->pMaterial->color;
    } else {
        rec.pMaterial = pMt->pMaterial;
        rec.color = ((pMesh->pVertices[pMt->indices[0]].color * b0)
                  +  (pMesh->pVertices[pMt->indices[1]].color * b1)
                  +  (pMesh->pVertices[pMt->indices[2]].color * b2)) / 3;
    }
    return true;
}

// レイキャスト
template <typename NODE_T>
int QBVHBase<NODE_T>::RayCast(vector<HitRecord>& hits, int nHits, const Ray &r, float tmin, float tmax) const
{
    if (!bbox_.RayIntersect(r, tmin, tmax)) return nHits;
    
    return RayCastBranch(hits, nHits, pNodes_[0], r, tmin, tmax);
}

// 葉に対するレイキャスト
template <typename NODE_T>
int QBVHBase<NODE_T>::RayCastLeaf(vector<HitRecord>& hits, int nHits, Leaf& leaf, const Ray &r, float tmin, float tmax) const
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
template <typename NODE_T>
int QBVHBase<NODE_T>::QSplit(const IShape** pShapes, int nShapes, float pivot, int axis)
{
    int mid_idx = 0;
    for (int i=0; i<nShapes; i++) {
        BBox bbox = pShapes[i]->BoundingBox();
        // 1頂点でもpivotより小さければ、小さい方に振り分ける
        // そうしないとレイの方向符号を使った子ノードの処理順ソートした場合に、
        //double centroid = ((bbox.Min()).e[axis] + (bbox.Max()).e[axis]) * 0.5f;
        double centroid = bbox.Min().e[axis];
        
        if (centroid <= pivot)
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

//--------------------------------------------------------------------------------------------------
// 明示的インスタンス化
//--------------------------------------------------------------------------------------------------
template class QBVHBase<SISD_QBVH_NODE>;
template class QBVHBase<SIMD_QBVH_NODE>;
