#include "Scene.h"
#include "SIMD_QBVH.h"
#include <cassert>
#include "simd.h"
#include "Ray.h"

using namespace std;

SIMD_QBVH::SIMD_QBVH()
{
}

SIMD_QBVH::~SIMD_QBVH()
{
}

ShapeType SIMD_QBVH::GetType() const
{
    return ST_QBVH_SIMD;
}

// 枝ノード構築
void SIMD_QBVH::BuildBranch(int iCurrNode, const IShape** pShapes, int nShapes)
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
    {
        SIMD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
        pNode->axis0 = axis;
        pNode->axis1 = axis_l;
        pNode->axis2 = axis_r;
    }
    
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
            SIMD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
            SetChildEmpty(pNode->children[i]);
        }
        // 子が4より多ければさらに枝に分割
        else if (nChildren[i] > 4) {
            //printf("=\n");
            //assert(iNodes_ < nNodeCapacity_);
            SIMD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
            SetChildNode(pNode->children[i], iNodes_);
            AddNewNode(); // increments iNodes_
            BuildBranch(iNodes_ - 1, ppChildShapes[i], nChildren[i]);
        }
        // 子が4以下だったら葉を作る
        else {
            //printf("*\n");
            SIMD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
            SetChildLeaf(pNode->children[i], iLeaves_);
            Leaf& rLeaf = AddNewLeaf();
            rLeaf.iTriangles = iTriangles_;
            rLeaf.iOtherPrims = iOtherPrims_;
            BuildLeaf(rLeaf.nTriangles, rLeaf.nOtherPrims, ppChildShapes[i], nChildren[i]);
        }
    }
    
    // 子のBBoxを内包するBBoxをこのノードのBBoxとする
    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            pNodes_[iCurrNode].bboxes[i][j] = _mm_set_ps(
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

// 全体のBBoxを作成
void SIMD_QBVH::MakeWholeBBox()
{
    Vec3 minv;
    Vec3 maxv;
    for (int xyz=0; xyz<3; xyz++) {
        Um128 bbox_min = pNodes_[0].bboxes[0][xyz];
        Um128 bbox_max = pNodes_[0].bboxes[1][xyz];
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
    bbox_ = BBox(minv, maxv);
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
        tmin, _mm_mul_ps(_mm_sub_ps(bboxes[sign[1]][1], orig[1]), idir[1])
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
        
        if (((1 << i) & hit) && !IsChildEmpty(rNode.children[i])) {
            //printf("hit!\n");
            int iChild = GetChildIndex(rNode.children[i]);
            
            // 枝
            if (IsChildNode(rNode.children[i])) {
                bool hit = IntersectBranch(pNodes_[iChild], r, tmin, rec);
                ret = ret || hit;
            }
            // 葉
            else {
                //printf("Leaf!\n");
                Leaf& rLeaf = pLeaves_[iChild];
                assert((rLeaf.nTriangles + rLeaf.nOtherPrims) > 0);
                bool hit = IntersectLeaf(rLeaf, r, tmin, rec);
                ret = ret || hit;
            }
        }
    }

    depth_--;
    return ret;
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
            int iChild = GetChildIndex(rNode.children[i]);
            
            // 枝
            if (IsChildNode(rNode.children[i])) {
                nHits = RayCastBranch(hits, nHits, pNodes_[iChild], r, tmin, tmax);
            }
            // 葉
            else {
                Leaf& rLeaf = pLeaves_[iChild];
                assert((rLeaf.nTriangles + rLeaf.nOtherPrims) > 0);
                nHits = RayCastLeaf(hits, nHits, rLeaf, r, tmin, tmax);
            }
        }
    }
    
    return nHits;
}
