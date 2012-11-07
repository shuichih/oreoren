#include <cmath>
#include <cstdlib>
#include "LightSource.h"
#include <cassert>

//--------------------------------------------------------------------------------

PointLightSource::PointLightSource()
: position_()
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position_.x + position_.y + position_.z
                              + intensity_.x + intensity_.y + intensity_.z);
}

PointLightSource::PointLightSource(const Vec3& position, const Vec3& intensity)
: LightSource(Lit_Point, intensity)
, position_(position)
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position.x + position.y + position.z
                              + intensity.x + intensity.y + intensity.z);
}

Ray PointLightSource::GenerateRay() const
{
	Vec3 dir;
    real theta = (real)M_PI * (real)erand48(xi_);
    real phi = 2.0f*(real)M_PI * (real)erand48(xi_);
    dir.x = sinf(theta) * cosf(phi);
    dir.y = sinf(theta) * sinf(phi);
    dir.z = cosf(theta);
    
	return Ray(position_, dir);
}


//--------------------------------------------------------------------------------

/*
AreaLightSource::AreaLightSource()
: LightSource(Lit_Area, Vec3(10000, 10000, 10000))
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(intensity_.x + intensity_.y + intensity_.z);
}
*/

AreaLightSource::AreaLightSource(
    const Vec3& p0,
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const Vec3& intensity
)
: LightSource(Lit_Area, intensity)
{
    p_[0] = p0;
    p_[1] = p1;
    p_[2] = p2;
    p_[3] = p3;
    normal_ = ((p1 - p0) % (p2 - p0));
    if (normal_.length() == 0) {
        normal_ = ((p0 - p2) % (p0 - p3));
        if (normal_.length() == 0) {
            assert(false);
            normal_ = Vec3(0, -1, 0);
        }
    }
    normal_.normalize();
    
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(intensity.x + intensity.y + intensity.z);
}

Ray AreaLightSource::GenerateRay() const
{
    real beta = (real)erand48(xi_);
    real gamma = 1.0f - beta;
//    Vec3 pos = ((real)erand48(xi_) < 0.5f) ?
//        (p_[0] + p_[1] * beta + p_[2] * gamma) * 0.5f :
//        (p_[0] + p_[2] * beta + p_[3] * gamma) * 0.5f;
    //Vec3 pos = (p_[0] + p_[1] * beta + p_[2] * gamma) * 0.5f;
    Vec3 pos = (p_[0] + p_[2] * beta + p_[3] * gamma) * 0.5f;
    
	Vec3 dir = Ray::CosRay(normal_, xi_);
    
    // @todo ライトのポリゴン化、レイトレしてライトポリゴンに当たったら明るさは直接計算して返す。
    
	return Ray(pos, dir);
}

//--------------------------------------------------------------------------------

AreaLightShape::AreaLightShape(AreaLightSource* pLitSrc, const Vec3 ps[3], const RGB& color, Refl_t refl)
: Triangle(ps[0], ps[1], ps[2], color, refl)
{
    Vec3 edge1 = p1 - p0;
    Vec3 edge2 = p2 - p0;
    float area = (edge2 % edge1).length() * 0.5f;
    irradiance_ = pLitSrc->GetIntensity() / area;
}

Vec3 AreaLightShape::SelfIrradiance()
{
    return irradiance_;
}

bool AreaLightShape::Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const
{
    if (rec.hitLit)
    {
        rec.pShape = this;
        return Triangle::Intersect(r, tmin, tmax, rec);
    }
    return false;
}

