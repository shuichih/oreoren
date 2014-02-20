#include "BRDF.h"
#include "Scene.h"
#include "Sampler.h"
#include "ONB.h"


IBRDF::~IBRDF()
{
}

//--------------------------------------------------------------------------------

Lambertian::Lambertian(float kd, const RGB& cd)
: kd_(kd)
, cd_(cd)
, pSampler_(new Sampler())
{
}

RGB Lambertian::F(const HitRecord& hr, const Vec3& wi, const Vec3& wo) const
{
    return kd_ * cd_ * PI_INV;
}

/**
 * from Kevin Suffern, Ray Tracing from the Ground Up, 26.2.2
 * サンプルする関数
 */
RGB Lambertian::SampleF(const HitRecord& hr, Vec3& wi, const Vec3& wo, float& pdf) const
{
    ONB onb = ONB::InitFromW(hr.normal);
    Vec3 sp = pSampler_->SampleHemisphere();
    wi = onb.Transform(sp);
    
    pdf = hr.normal.dot(wi) * PI_INV;
    
    return kd_ * cd_ * PI_INV;
}


RGB Lambertian::Rho(const HitRecord& hr, const Vec3& wo) const
{
    return kd_ * cd_;
}

//--------------------------------------------------------------------------------

PerfectSpecular::PerfectSpecular(float kr, const RGB& cr)
: kr_(kr)
, cr_(cr)
{
}

RGB PerfectSpecular::F(const HitRecord& hr, const Vec3& wi, const Vec3& wo) const
{
    // undefined, because of the BRDF contains delta function.
    return RGB();
}

RGB PerfectSpecular::SampleF(const HitRecord& hr, Vec3& wi, const Vec3& wo, float& pdf) const
{
    float ndotwo = hr.normal.dot(wo);
    wi = -wo + 2.0f * hr.normal * ndotwo;
    
    pdf = hr.normal.dot(wi);
    
    return kr_ * cr_;
}

RGB PerfectSpecular::Rho(const HitRecord& hr, const Vec3& wo) const
{
    // undefined, because of the BRDF contains delta function.
    return RGB();
}

//--------------------------------------------------------------------------------
GlossySpecular::GlossySpecular(float ks, const RGB& cs, float exp)
: ks_(ks)
, cs_(cs)
, exp_(exp)
{
}

RGB GlossySpecular::F(const HitRecord& hr, const Vec3& wi, const Vec3& wo) const
{
    // undefined, because of the BRDF contains delta function.
    return RGB();
}

RGB GlossySpecular::SampleF(const HitRecord& hr, Vec3& wi, const Vec3& wo, float& pdf) const
{
    float ndotwo = hr.normal.dot(wo);
    Vec3 r = -wo + 2.0f * hr.normal * ndotwo;
    
    ONB onb = ONB::InitFromW(hr.normal);
    Vec3 sp = pSampler_->SampleHemisphere();
    wi = onb.Transform(sp);
    
    if (hr.normal.dot(wi) < 0.f) { // reflected ray is below surface
        wi = onb.Transform(-sp);
    }
    
    float phongLobe = powf(r.dot(wi), exp_);
    pdf = phongLobe * hr.normal.dot(wi);
    
    return ks_ * cs_ * phongLobe;
    
    // Glossyのexpに応じた反射方向にするのではなく、diffuseと同じe=1の半球サンプリングにして
    // 方向の代わりに色にphongLobeを掛けている。
    // expに応じた反射方向にして色はks*csとして、だと、鋭いスペキュラの場合ノイズが出やすいんだろうか?
    // 両方試して比べるといいかもしれない。
    // @todo Phong materialの所を読んでサンプリング方法を確認
}

RGB GlossySpecular::Rho(const HitRecord& hr, const Vec3& wo) const
{
    // undefined, because of the BRDF contains delta function.
    return RGB();
}

// @todo Glossy Reflection
/*
 Glossyではpdfを使うので計算
 Glossyで直接照明をPhongで計算するか?
 GlossyReflectorがPhongを継承してる
 Phongの直接照明+Glossyの照明としている
 Light方向への反射をサンプルすると遅いので、ライト以外の物体への(からの)反射(==間接光)とLight
 からの直接照明を分けている。
 よってGlossyのtraceではライトに直接当たらないようにする必要がある。
 直接光計算は本当に早いのか。エリアライトの場合多数サンプルが必要だが、それでも?
 Phongの実装を見れば早さが分かりそう
 PerfectSpecularのTraceでは当てないとライトのシルエット(ハイライト)は出なくなる
 しかし明るすぎでライトの縁にジャギーが出る
 Mental Rayとかでは置いたライトが映るんだろうか？
 
 GlossyReflectorはambient, diffuse, specularも含んでるが, phongを継承したから..
 複数BRDFの複合マテリアルという事?
 
 @todo GlossySpecular BRDFを実装
 @todo ２５章序盤およびPhongの章を読み込む
 @todo 14.7 Materials読む
 
 GlossyReflector, Matte, Reflectiveが同じ役割っぽいが共通の基底クラスがある? Materialとか。
 area_light_shade()とpath_shade()で
 
 global_shadeってなんだ
 */
// @todo Diffuse+Glossy reflection
// @todo various materials
// キューブテクスチャでさっさとIBL, サンプルするだけ、スケール適当に掛けていい感じに出す
