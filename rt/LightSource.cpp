#include <cmath>
#include <cstdlib>
#include <cassert>
#include "LightSource.h"
#include "ONB.h"
#include "vector3.h"
#include "Ray.h"

//--------------------------------------------------------------------------------

PointLightSource::PointLightSource()
: position_()
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position_.x + position_.y + position_.z
                              + flux_.x + flux_.y + flux_.z);
}

PointLightSource::PointLightSource(const Vec3& position, const Vec3& flux)
: LightSource(Lit_Point, flux)
, position_(position)
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position.x + position.y + position.z
                              + flux.x + flux.y + flux.z);
}

Ray PointLightSource::GenerateRay() const
{
	Vec3 dir;
    real theta = acosf(2 * (real)erand48(xi_) - 1);
    real phi = 2.0f*(real)M_PI * (real)erand48(xi_);
    dir.x = sinf(theta) * cosf(phi);
    dir.y = sinf(theta) * sinf(phi);
    dir.z = cosf(theta);
    
	return Ray(position_, dir);
}

Vec3 PointLightSource::DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra) const
{
    // penumbra
    // == 1.0, 完全に光が当たっている領域: サンプリング、シャドウレイなし
    // == 0.0, penumbra完全に影の領域　: 黒
    // else  , 半影領域              : サンプリング、シャドウレイあり
   
    if (penumbra == 0.f) {
        return Vec3();
    }
    
    HitRecord rec;
    Vec3 ldir = position_ - pos;
    float r2 = ldir.lengthSquared();
    float r = sqrtf(ldir.lengthSquared());
    ldir.normalize();
    Ray ray(pos, ldir);
    float pdf = 1.f / (4 * PI * r2); // pdfではないような?
    if (penumbra == 1.f || !scene.Intersect(ray, EPSILON, r, rec)) {
        return normal.dot(ldir) * flux_ * pdf;
    } else {
        return Vec3();
    }
}


//--------------------------------------------------------------------------------

AreaLightSource::AreaLightSource(
    const Vec3& p0,
    const Vec3& p1,
    const Vec3& p2,
    const Vec3& p3,
    const Vec3& flux,
    int nSamples
)
: LightSource(Lit_Area, flux)
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
    
    float s0 = CalcTriangleArea(p0, p1, p2);
    float s1 = CalcTriangleArea(p0, p2, p3);
    area_ = s0 + s1;
    pdf_ = 1.f / area_;
    irradiance_ = flux_ / area_;
    
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(flux.x + flux.y + flux.z);
}

float AreaLightSource::CalcTriangleArea(const Vec3& p0, const Vec3& p1, const Vec3& p2)
{
    // Heron's formula
    float e0 = (p1 - p0).length();
    float e1 = (p2 - p1).length();
    float e2 = (p0 - p2).length();
    float s = 0.5f * (e0 + e1 + e2);
    return sqrtf(s * (s-e0) * (s-e1) * (s-e2));
    
    // こっちでもいいような
//    Vec3 edge1 = p1 - p0;
//    Vec3 edge2 = p2 - p0;
//    float area = (edge2 % edge1).length() * 0.5f; // 三角形の面積だから*0.5f
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

Vec3 AreaLightSource::DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra) const
{
    // penumbra
    // == 1.0, 完全に光が当たっている領域: サンプリング、シャドウレイなし
    // == 0.0, penumbra完全に影の領域　: 黒
    // else  , 半影領域              : サンプリング、シャドウレイあり
    
    if (penumbra == 0.f) {
        return Vec3();
    }
    
    if (normal.dot(this->normal_) > 0.0f) {
        return Vec3();
    }

    Vec3 irrad;
    const float nSampleInv = 1.0f / nSamples_;
    HitRecord rec;
    
    for (int i=0; i<nSamples_; i++) {
        real b0 = (real)erand48(xi_);
        real b1 = (1.0f - b0) * (real)erand48(xi_);
        real b2 = 1.0f - b0 - b1;
        Vec3 lpos = ((real)erand48(xi_) < 0.5f) ?
            (p_[0] * b0 + p_[1] * b1 + p_[2] * b2):
            (p_[0] * b0 + p_[2] * b1 + p_[3] * b2);
        
        Vec3 ldir = lpos - pos;
        float r2 = ldir.lengthSquared();
        float r = sqrtf(r2);
        ldir.normalize();
        
        Ray ray(pos, ldir);
        if (penumbra == 1.f || !scene.Intersect(ray, EPSILON, r, rec)) {
            float cosX = this->normal_.dot(-1.f * ldir);
            float cosY = normal.dot(ldir);
            irrad += (irradiance_ * nSampleInv) * cosX * cosY / (pdf_ * r2);
        }
    }
    return irrad;
}


//--------------------------------------------------------------------------------

AreaLightShape::AreaLightShape(AreaLightSource* pLitSrc, const Vec3 ps[3], const RGB& color, Refl_t refl)
: Triangle(ps[0], ps[1], ps[2], color, refl)
{
    pLightSource_ = pLitSrc;
}

Vec3 AreaLightShape::SelfIrradiance() const
{
    return pLightSource_->GetIrradiance();
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


//--------------------------------------------------------------------------------

/*
SphereLightSource::SphereLightSource()
: position_()
, radius_(1)
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position_.x + position_.y + position_.z
                              + flux_.x + flux_.y + flux_.z);
}
*/

SphereLightSource::SphereLightSource(const Vec3& position, float radius, const Vec3& flux, int nSamples)
: LightSource(Lit_Sphere, flux)
, position_(position)
, radius_(radius)
, nSamples_(nSamples)
{
	xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)(position.x + position.y + position.z
                              + flux.x + flux.y + flux.z);
    irradiance_ = flux_ / (4 * PI * radius_ * radius_);
}

Ray SphereLightSource::GenerateRay() const
{
	Vec3 dir;
    real theta = acosf(2 * (real)erand48(xi_) - 1);
    real phi = 2.0f*(real)M_PI * (real)erand48(xi_);
    dir.x = sinf(theta) * cosf(phi);
    dir.y = sinf(theta) * sinf(phi);
    dir.z = cosf(theta);
    
	return Ray(position_ + dir*radius_, dir);
}

Vec3 SphereLightSource::DirectLight(const Vec3& pos, const Vec3& normal, const Scene& scene, float penumbra) const
{
    // penumbra
    // == 1.0, 完全に光が当たっている領域: サンプリング、シャドウレイなし
    // == 0.0, penumbra完全に影の領域　: 黒
    // else  , 半影領域              : サンプリング、シャドウレイあり
    
    if (penumbra == 0.f) {
        return Vec3();
    }
    
    // from Realistic Ray Tracing Second Edition, Chapter 13
    
    // シェーディング点から光源中心までのベクトル
    Vec3 xc = position_ - pos;
    float xc_len = xc.length();
    Vec3 w = xc.normalize();
    
    // シェーディング点から見た球光源の立体角
    float rxc2 = radius_ / xc_len;
    rxc2 = rxc2 * rxc2;
    float cos_a_max = sqrtf(1 - rxc2);
    //float a_max = acosf(cos_a_max);
    
    // 立体角でのdensity function
    float q = 1.f / (2.f*PI * (1.f - cos_a_max));
    
    Vec3 irrad;
    HitRecord rec;
    for (int i=0; i<nSamples_; i++) {
        
        // 球光源への立体角の範囲内でdensity function qに沿ったランダム方向を決める
        float xi1 = (float)erand48(xi_);
        float xi2 = (float)erand48(xi_);
        float cos_a = 1 + xi1 * (cos_a_max - 1);
        float phi = 2 * PI * xi2;
        float sin_a = sqrtf(1 - cos_a*cos_a);
        ONB uvw;
        uvw.InitFromW(w);
        Vec3 a = uvw.u_ * cosf(phi) * sin_a
               + uvw.v_ * sinf(phi) * sin_a
               + uvw.w_ * cos_a;
        
        // ランダム方向と球光源の交点と、交点の法線
        HitRecord rec;
        bool ret = pShape_->Intersect(Ray(pos, a), 0, REAL_MAX, rec);
        if (!ret) {
            rec.t = xc_len - radius_; // 誤差等で万一当たってない場合
        }
        Vec3 lx = pos + a * rec.t;
        Vec3 ln = (lx - position_).normalize();
        
        // 領域でのdensity function
        float xlx2 = (pos-lx).lengthSquared();
        float pdf = q * (-1.f*w).dot(ln) / xlx2;
    
        //
        rec = HitRecord();
        Ray ray(pos, a);
        if (penumbra == 1.f || !scene.Intersect(ray, EPSILON, xc_len, rec)) {
            float cosX = ln.dot(-1.f * a);
            float cosY = normal.dot(a);
            Vec3 lit_irrad = flux_ / (4 * PI * radius_ * radius_ * xlx2);
            irrad += cosX * cosY * lit_irrad / (pdf * nSamples_);
        }
    }
    
    return irrad;
}

//--------------------------------------------------------------------------------

SphereLightShape::SphereLightShape(SphereLightSource* pLitSrc, float radius, const Vec3& pos, const RGB& color, Refl_t refl)
: Sphere(radius, pos, color, refl)
{
    pLightSource_ = pLitSrc;
    pLightSource_->SetShape(this);
}

Vec3 SphereLightShape::SelfIrradiance()
{
    return pLightSource_->GetIrradiance();
}

bool SphereLightShape::Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const
{
    if (rec.hitLit)
    {
        rec.pShape = this;
        return Sphere::Intersect(r, tmin, tmax, rec);
    }
    return false;
}
