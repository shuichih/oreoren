#include "BBox.h"
#include "simd.h"
#include "Ray.h"
#include <cassert>

BBox::BBox(const Vec3& min, const Vec3& max)
{
    pp[0] = min;
    pp[1] = max;
#ifdef BBOX_USE_SIMD
    xmpp0 = _mm_set_ps(pp[0].x, pp[0].y, pp[0].z, 0);
    xmpp1 = _mm_set_ps(pp[1].x, pp[1].y, pp[1].z, 0);
#endif
}

BBox::BBox(float minX, float minY, float minZ, float maxX, float maxY, float maxZ)
{
    pp[0].x = minX;
    pp[0].y = minY;
    pp[0].z = minZ;
    pp[1].x = maxX;
    pp[1].y = maxY;
    pp[1].z = maxZ;
#ifdef BBOX_USE_SIMD
    xmpp0 = _mm_set_ps(pp[0].x, pp[0].y, pp[0].z, 0);
    xmpp1 = _mm_set_ps(pp[1].x, pp[1].y, pp[1].z, 0);
#endif
}

// 2つのBBoxを内包する最小のBBoxを作成
BBox BBox::Surround(const BBox& b1, const BBox& b2)
{
    return BBox(
        Vec3(
            b1.Min().x < b2.Min().x ? b1.Min().x : b2.Min().x,
            b1.Min().y < b2.Min().y ? b1.Min().y : b2.Min().y,
            b1.Min().z < b2.Min().z ? b1.Min().z : b2.Min().z
        ),
       Vec3(
            b1.Max().x > b2.Max().x ? b1.Max().x : b2.Max().x,
            b1.Max().y > b2.Max().y ? b1.Max().y : b2.Max().y,
            b1.Max().z > b2.Max().z ? b1.Max().z : b2.Max().z
        )
    );
}

bool BBox::RayIntersect(const Ray &r, float tmin, float tmax) const
{
#ifdef BBOX_USE_SIMD
    // interval_min, interval_maxを徐々に限定して行くからxyzまとめてやるのは無理
    // Horizontal operations aren't supported in SSE...
    __m128 xmtmin = _mm_set1_ps(tmin);
    __m128 xmtmax = _mm_set1_ps(tmax);
    __m128 xmppmin = _mm_blendv_ps(xmpp0, xmpp1, r.xmDirSignMask);
    __m128 xmppmax = _mm_blendv_ps(xmpp1, xmpp0, r.xmDirSignMask);
    xmtmin = _mm_max_ps(
        xmtmin,
        _mm_mul_ps(_mm_sub_ps(xmppmin, r.xmo), r.xmInvDir)
    );
    xmtmax = _mm_min_ps(
        xmtmax,
        _mm_mul_ps(_mm_sub_ps(xmppmax, r.xmo), r.xmInvDir)
    );
    //__m128 tmp1 = _mm_cmpge_ps(xmtmax, xmtmin);
    //int tmp2 = _mm_movemask_ps(tmp1);
    //int tmp3 = tmp2 & 0x0e;
    //bool tmp = (_mm_movemask_ps(_mm_cmpge_ps(xmtmax, xmtmin)) & 0x0e) == 0x0e;
    //bool tmp4 = tmp3 == 0x0e;
    return (_mm_movemask_ps(_mm_cmpge_ps(xmtmax, xmtmin)) & 0x0e) == 0x0e;
#else
    
    float interval_min = tmin;
    float interval_max = tmax;
    
    // X軸方向のBB内区間をRayのパラメータ空間(Ray方向のt)で計算
    int sign = r.sign[0];
    float t0 = (pp[sign].x - r.o.x) * r.idir.x;
    float t1 = (pp[1-sign].x - r.o.x) * r.idir.x;
    // BB内に入っている区間にtの範囲を限定していく
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    // BB範囲外ならヒットしない
    if (interval_min > interval_max) return false;
    
    // Y軸、Z軸についても同様に処理
    sign = r.sign[1];
    t0 = (pp[sign].y - r.o.y) * r.idir.y;
    t1 = (pp[1-sign].y - r.o.y) * r.idir.y;
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    if (interval_min > interval_max) return false;
    
    sign = r.sign[2];
    t0 = (pp[sign].z - r.o.z) * r.idir.z;
    t1 = (pp[1-sign].z - r.o.z) * r.idir.z;
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    if (interval_min > interval_max) return false;
    
    return true;
#endif
}


