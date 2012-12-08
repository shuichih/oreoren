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

Vec3 PointLightSource::DirectLight(const Vec3& pos, const Vec3& normal) const
{
    Vec3 ldir = position_ - pos;
    float r2 = ldir.square_length();
    ldir.normalize();
    return normal.dot(ldir) * intensity_ / (4 * PI * r2);
}


//--------------------------------------------------------------------------------

AreaLightSource::AreaLightSource(
    const Vec3& p0,
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const Vec3& intensity,
    int nSamples
)
: LightSource(Lit_Area, intensity)
, nSamples_(nSamples)
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
    
    // Heron's formula
    float s0 = CalcTriangleArea(p0, p1, p2);
    float s1 = CalcTriangleArea(p0, p2, p3);
    area_ = s0 + s1;
    
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(intensity.x + intensity.y + intensity.z);
}

float AreaLightSource::CalcTriangleArea(const Vec3& p0, const Vec3& p1, const Vec3& p2)
{
    // Heron's formula
    float e0 = (p1 - p0).length();
    float e1 = (p2 - p1).length();
    float e2 = (p0 - p2).length();
    float s = 0.5f * (e0 + e1 + e2);
    return sqrtf(s * (s-e0) * (s-e1) * (s-e2));
}

Ray AreaLightSource::GenerateRay() const
{
    real b0 = (real)erand48(xi_);
    real b1 = (1.0f - b0) * (real)erand48(xi_);
    real b2 = 1.0f - b0 - b1;
    Vec3 pos = ((real)erand48(xi_) < 0.5f) ?
        (p_[0] * b0 + p_[1] * b1 + p_[2] * b2):
        (p_[0] * b0 + p_[2] * b1 + p_[3] * b2);
    
	Vec3 dir = Ray::CosRay(normal_, xi_);
    
	return Ray(pos, dir);
}

Vec3 AreaLightSource::DirectLight(const Vec3& pos, const Vec3& normal) const
{
    // 遮蔽なし

    if (normal.dot(this->normal_) > 0.0f) {
        return Vec3();
    }

    Vec3 irrad;
    const float nSampleInv = 1.0f / nSamples_;
    for (int i=0; i<nSamples_; i++) {
        real b0 = (real)erand48(xi_);
        real b1 = (1.0f - b0) * (real)erand48(xi_);
        real b2 = 1.0f - b0 - b1;
        Vec3 lpos = ((real)erand48(xi_) < 0.5f) ?
            (p_[0] * b0 + p_[1] * b1 + p_[2] * b2):
            (p_[0] * b0 + p_[2] * b1 + p_[3] * b2);
        
        Vec3 ldir = lpos - pos;
        float r2 = ldir.square_length();
        ldir.normalize();
        float cosX = this->normal_.dot(-1 * ldir);
        //if (cosX < 0.0f) continue;
        float cosY = normal.dot(ldir);
        //if (cosY < 0.0f) continue;
        irrad += (intensity_ * nSampleInv) * cosX * cosY / r2;
    }
    return irrad;
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

