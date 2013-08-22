#include "Scene.h"
#include "SISD_QBVH.h"
#include <cassert>
#include "simd.h"
#include "Ray.h"

using namespace std;


SISD_QBVH::SISD_QBVH()
{
}

SISD_QBVH::~SISD_QBVH()
{
}

ShapeType SISD_QBVH::GetType() const
{
    return ST_QBVH_SISD;
}

// 枝ノード構築
void SISD_QBVH::BuildBranch(int iCurrNode, const IShape** pShapes, int nShapes)
{
    depth_++;
    
    // qsplitのピボットとして使用するためにBBoxの中間点を求める
    BBox bbox = SurroundBBox(pShapes, nShapes);
    
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
        SISD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
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
    
    for (int i=0; i<4; i++) {
        
        // bbox求める
        pNodes_[iCurrNode].bboxes[i] = SurroundBBox(ppChildShapes[i], nChildren[i]);
        
        // 子が空なら特別なindexをセット
        if (nChildren[i] == 0) {
            //printf("-\n");
            SISD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
            SetChildEmpty(pNode->children[i]);
        }
        // 子が4より多ければさらに枝に分割
        else if (nChildren[i] > 4) {
            //printf("=\n");
            SISD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
            SetChildNode(pNode->children[i], iNodes_);
            AddNewNode(); // increments iNodes_
            BuildBranch(iNodes_ - 1, ppChildShapes[i], nChildren[i]);
        }
        // 子が4以下だったら葉を作る
        else {
            //printf("*\n");
            SISD_QBVH_NODE* pNode = &pNodes_[iCurrNode];
            SetChildLeaf(pNode->children[i], iLeaves_);
            Leaf& rLeaf = AddNewLeaf();
            rLeaf.iTriangles = iTriangles_;
            rLeaf.iOtherPrims = iOtherPrims_;
            BuildLeaf(rLeaf.nTriangles, rLeaf.nOtherPrims, ppChildShapes[i], nChildren[i]);
        }
    }
    
    depth_--;
    //printf("rNode=%p iTri=%d iOP=%d\n", &rNode, iTriangles_, iOtherPrims_);
}

// 全体のBBoxを作成
void SISD_QBVH::MakeWholeBBox()
{
    BBox bbox0 = BBox::Surround(pNodes_[0].bboxes[0], pNodes_[0].bboxes[1]);
    BBox bbox1 = BBox::Surround(pNodes_[0].bboxes[2], pNodes_[0].bboxes[3]);
    bbox_ = BBox::Surround(bbox0, bbox1);
}

// 枝との交差判定
bool SISD_QBVH::IntersectBranch(SISD_QBVH_NODE& rNode, const Ray &r, float tmin, HitRecord& rec) const
{
    depth_++;
    
    bool ret = false;
    for (int i=0; i<4; i++) {
        if (IsChildEmpty(rNode.children[i])) {
            // 空ノード
            continue;
        }
        
        /*
        printf("%d  %f %f %f  %f %f %f\n", depth_,
            rNode.bboxes[i].pp[0].e[0], rNode.bboxes[i].pp[0].e[1], rNode.bboxes[i].pp[0].e[2],
            rNode.bboxes[i].pp[1].e[0], rNode.bboxes[i].pp[1].e[1], rNode.bboxes[i].pp[1].e[2]
        );
        */
        if (rNode.bboxes[i].RayIntersect(r, tmin, rec.t)) {
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
int SISD_QBVH::RayCastBranch(vector<HitRecord>& hits, int nHits, SISD_QBVH_NODE& rNode, const Ray &r, float tmin, float tmax) const
{
    for (int i=0; i<4; i++) {
        if (IsChildEmpty(rNode.children[i])) {
            // 空ノード
            continue;
        }
        
        if (rNode.bboxes[i].RayIntersect(r, tmin, tmax)) {
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

