#include "BBox.h"


BBox::BBox(const Vec3& min, const Vec3& max)
{
    pp[0] = min;
    pp[1] = max;
}

BBox::BBox(float minX, float minY, float minZ, float maxX, float maxY, float maxZ)
{
    pp[0].x = minX;
    pp[0].y = minY;
    pp[0].z = minZ;
    pp[1].x = maxX;
    pp[1].y = maxY;
    pp[1].z = maxZ;
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
    float interval_min = tmin;
    float interval_max = tmax;
    
    // X軸方向のBB内区間をRayのパラメータ空間(Ray方向のt)で計算
    int dirSign = r.dirSign[0];
    float t0 = (pp[dirSign].x - r.o.x) * r.invDir.x;
    float t1 = (pp[1-dirSign].x - r.o.x) * r.invDir.x;
    // BB内に入っている区間にtの範囲を限定していく
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    // BB範囲外ならヒットしない
    if (interval_min > interval_max) return false;
    
    // Y軸、Z軸についても同様に処理
    dirSign = r.dirSign[1];
    t0 = (pp[dirSign].y - r.o.y) * r.invDir.y;
    t1 = (pp[1-dirSign].y - r.o.y) * r.invDir.y;
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    if (interval_min > interval_max) return false;
    
    dirSign = r.dirSign[2];
    t0 = (pp[dirSign].z - r.o.z) * r.invDir.z;
    t1 = (pp[1-dirSign].z - r.o.z) * r.invDir.z;
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    if (interval_min > interval_max) return false;
    
    return true;
}


