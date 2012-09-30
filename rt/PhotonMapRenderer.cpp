#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "Common.h"
#include "Scene.h"
#include "PhotonMap.h"
#include "LightSource.h"
#include "PhotonMapRenderer.h"
#include "PhotonFilter.h"
#include "BVH.h"
#include "Config.h"

using namespace std;

//----------------------------------------------------------------
PhotonMapRenderer::PhotonMapRenderer()
: pBVH_(NULL)
, pPhotonMap_(NULL)
, pFilter_(NULL)
{
}

PhotonMapRenderer::~PhotonMapRenderer()
{
    delete pPhotonMap_;
    delete pFilter_;
}

void PhotonMapRenderer::SetConfig(const Config& config)
{
    pConfig_ = &config;
    pPmConfig_ = &pConfig_->photonMapConf;
    
    if (pPhotonMap_ != NULL) delete pPhotonMap_;
    pPhotonMap_ = new Photon_map(pPmConfig_->nPhotons * pPmConfig_->maxPhotonBounce); // ちと多め
    
    if (pFilter_ != NULL) delete pFilter_;
    if (pPmConfig_->enableConeFilter) {
        pFilter_ = new ConeFilter(pPmConfig_->coneFilterK);
        ((ConeFilter*)pFilter_)->SetK(pPmConfig_->coneFilterK);
        pPhotonMap_->SetFilter(pFilter_);
    }
    else {
        pPhotonMap_->SetFilter(NULL);
    }
    
    pPhotonMap_->SetEstimateEllipseScale(pPmConfig_->estimateEllipseScale);
}

bool PhotonMapRenderer::Intersect(const Ray& r, HitRecord& rec)
{
    if (pPmConfig_->useBVH && pBVH_) {
        return pBVH_->Intersect(r, EPSILON, REAL_MAX, rec);
    }
    else {
        rec.t = REAL_MAX;
        size_t nShapes = pScene_->shapes_.size();
        const vector<const Shape*>& shapes = pScene_->shapes_;
        for (size_t i=nShapes; i--;) {
            shapes[i]->Intersect(r, EPSILON, rec.t, rec);
        }
        
        return rec.t < REAL_MAX;
    }
    
    assert(true);
    return false;
}

void PhotonMapRenderer::PhotonTracing(const Ray& r, float power[3], int depth)
{
    // max refl
    if (++depth > pPmConfig_->maxPhotonBounce)
    {
        return;
    }

    HitRecord rec;
//    if (!Intersect(r, t, id))
    if (!Intersect(r, rec))
        return;
    
    //const Shape &obj = *g_shapes[id];    // the hit object
    Vec3 x = r.o + r.d * rec.t;
    Vec3 n = rec.normal;
    Vec3 nl = n.dot(r.d) < 0.f ? n : n*-1.f; // 交点の法線, 裏面ヒット考慮
    Vec3 color = rec.color;
    Refl_t refl = rec.refl;
    
    // Ideal DIFFUSE reflection
    if (refl == DIFF) {
        pPhotonMap_->store(power, x.e, r.d.e);
        
        float ave_refl = (color.x + color.y + color.z) / 3.f;
        if ((real)erand48(xi_) < ave_refl)
        {
            // @todo importance samplingになってる？
            real r1 = 2.f*(real)(M_PI*erand48(xi_));
            real r2 = (real)erand48(xi_); // => 1-cos^2θ = 1-sqrt(1-r_2)^2 = r_2
            real r2s = sqrtf(r2);    // => sinθ = sqrt(1-cos^2θ) = sqrt(r_2)
            Vec3 w = nl; // normal
            Vec3 u = ((fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) % w).normalize(); // binormal
            Vec3 v = w % u; // tangent
            
            // ucosφsinθ + vsinφsinθ + wcosθ
            Vec3 d = (u*cosf(r1)*r2s + v*sinf(r1)*r2s + w*sqrtf(1-r2)).normalize();

            real ave_refl_inv = 1.0f / ave_refl;
            power[0] *= color.x * ave_refl_inv;
            power[1] *= color.y * ave_refl_inv;
            power[2] *= color.z * ave_refl_inv;
            PhotonTracing(Ray(x,d), power, depth);
        }
        //fprintf(stderr, "p (%f %f %f) (%f %f %f) (%f %f %f)\n", power[0], power[1], power[2], pos[0], pos[1], pos[2], dir[0], dir[1], dir[2]);
        return;
    }
    // Ideal SPECULAR reflection
    else if (refl == SPEC) {
        Ray refl(x, r.d - n * 2.f * n.dot(r.d));
        PhotonTracing(refl, power, depth);
        return;
    }
    else if (refl == PHONGMETAL) {
        real r1 = 2.f*(real)(M_PI*erand48(xi_));
        real r2 = (real)erand48(xi_);
        real exponent = 10.0f;
        real cosTheta = powf(1.f-r2, 1.0f/(exponent+1.0f));
        real sinTheta = sqrtf(1.0f - cosTheta*cosTheta);
        real rx = cosf(r1) * sinTheta;
        real ry = sinf(r1) * sinTheta;
        real rz = cosTheta;
        Vec3 w = r.d - n * 2.f * n.dot(r.d); // reflected ray
        Vec3 u = ((fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) % w).normalize(); // binormal
        Vec3 v = w % u; // tangent
        
        // ucosφsinθ + vsinφsinθ + wcosθ
        Vec3 d = (u*rx + v*ry + w*rz).normalize();
        
        Vec3 mx = x + n*1e-4f;
        PhotonTracing(Ray(mx, d), power, depth);
        return;
    }
    
    // Ideal dielectric REFRACTION
    Ray reflRay(x, r.d - n * 2.f * n.dot(r.d));
    bool into = n.dot(nl) > 0.f;  // Ray from outside going in?
    real airRefrIdx = 1.f;
    real grassRefrIdx = 1.5f;
    real nnt = into ? airRefrIdx/grassRefrIdx : grassRefrIdx/airRefrIdx;
    real ddn = r.d.dot(nl);   // レイと法線のcos
    real cos2t;
    
    // Total internal reflection
    if ((cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn)) < 0.f) {
        PhotonTracing(reflRay, power, depth);
        return;
    }
    
    // 屈折方向
    Vec3 tdir = (r.d * nnt - n * ((into ? 1.f : -1.f) * (ddn * nnt + sqrtf(cos2t)))).normalize();
    
    real a = grassRefrIdx - airRefrIdx;
    real b = grassRefrIdx + airRefrIdx;
    real R0 = a * a / (b * b); // 垂直反射率 0.25 / 2.5 = 0.1 @todo 根拠調査
    real c = 1.f - (into ? -ddn : tdir.dot(n));
    real fresnel = R0 + (1.f - R0)*c*c*c*c*c;
    // 0.25 と 0.5は結果を綺麗にするためのヒューリスティックな調整値。
    // 屈折する確率も反射する確率も最低限25%にするということ。例えば.1 + (.8 * fresnel)でもよい。
    // この調整が無ければRP, TPは1になって、下でRP, TP掛ける必要はなくなる。
    real P = fresnel;
    if ((real)erand48(xi_) < P)
        PhotonTracing(reflRay, power, depth);          // 反射
    else
        PhotonTracing(Ray(x, tdir), power, depth);  // 屈折

    return;
}


Vec3 PhotonMapRenderer::Irradiance(const Ray &r, int depth)
{
    // max refl
    if (++depth > pPmConfig_->maxRayBounce)
    {
        return Vec3();
    }

    HitRecord rec;
    if (!Intersect(r, rec)) return Vec3();
    //const Shape& obj = *g_shapes[id];       // the hit object
    Vec3 x = r.o + r.d * rec.t;
    Vec3 n = rec.normal;
    Vec3 nl = n.dot(r.d) < 0.f ? n : n * -1.f;   // 交点の法線
    Vec3 f = rec.color;
    Refl_t refl = rec.refl;
    
    
    // 0.5にしたらカラーが反射率になってるから暗くなるだけ。IDEALでない反射は扱えない。カラーと混ぜるとかもない。
    // Ideal DIFFUSE reflection
    if (refl == DIFF){
        float pos[3] = { x.x, x.y, x.z };
        float irrad[3] = { 0.f, 0.f, 0.f };
        pPhotonMap_->irradiance_estimate(irrad, pos, nl, pPmConfig_->estimateDist, pPmConfig_->nEstimatePhotons);
        //fprintf(stderr, "irrad %f %f %f\r", irrad[0], irrad[1], irrad[2]);
        return Vec3(f.x * irrad[0], f.y * irrad[1], f.z * irrad[2]);
    }
    else if (refl == SPEC) {
        // Ideal SPECULAR reflection
        return Irradiance(Ray(x,r.d-n*2.f*n.dot(r.d)), depth);
    }
    else if (refl == PHONGMETAL) {
        Vec3 irrad;
        for (int i = 0; i < 16; i++) {
            // Imperfect SPECULAR reflection
            real r1 = 2.f*(real)(M_PI*erand48(xi_));
            real r2 = (real)erand48(xi_);
            real exponent = 10.0f;
            real cosTheta = powf(1.f-r2, 1.0f/(exponent+1.0f));
            real sinTheta = sqrtf(1.0f - cosTheta*cosTheta);
            real rx = cosf(r1) * sinTheta;
            real ry = sinf(r1) * sinTheta;
            real rz = cosTheta;
            Vec3 w = r.d - n * 2.f * n.dot(r.d); // reflected ray
            Vec3 u = ((fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) % w).normalize(); // binormal
            Vec3 v = w % u; // tangent
            
            // ucosφsinθ + vsinφsinθ + wcosθ
            Vec3 rd = (u*rx + v*ry + w*rz).normalize();
            Vec3 mx = x + n*1e-4f; // 自己ヒットしないようにちょっと浮かす
            irrad += Irradiance(Ray(mx, rd), depth);
        }
        return irrad / 16;
    }
    
    // Ideal dielectric REFRACTION
    Ray reflRay(x, r.d - n * 2.f * n.dot(r.d));
    bool into = n.dot(nl) > 0.f; // Ray from outside going in?
    real airRefrIdx = 1.f;
    real grassRefrIdx = 1.5f;
    real nnt = into ? airRefrIdx/grassRefrIdx : grassRefrIdx/airRefrIdx;
    real ddn = r.d.dot(nl); // レイと法線のcos
    real cos2t;
    
    // Total internal reflection
    if ((cos2t = 1.f-nnt * nnt * (1.f - ddn * ddn)) < 0.f)
        return Irradiance(reflRay, depth);
    
    // 屈折方向
    Vec3 tdir = (r.d * nnt - n * ((into ? 1.f : -1.f) * (ddn * nnt + sqrtf(cos2t)))).normalize();
    
    real a = grassRefrIdx - airRefrIdx;
    real b = grassRefrIdx + airRefrIdx;
    real R0 = a * a / (b * b); // 垂直反射率 0.25 / 2.5 = 0.1 @todo 根拠調査
    real c = 1.f - (into ? -ddn : tdir.dot(n));
    real fresnel = R0 + (1.f - R0)*c*c*c*c*c;
    real Tr = 1.f - fresnel;
    // 反射屈折両方トレース
    return Irradiance(reflRay, depth) * fresnel + Irradiance(Ray(x,tdir), depth) * Tr;
}

void PhotonMapRenderer::RayTracing(Vec3* pColorBuf)
{
    const u32 w = pConfig_->windowWidth;
    const u32 h = pConfig_->windowHeight;
    const u32 nSub = pPmConfig_->nSubPixelsSqrt;
    const real subPixelFactor = 1.0f / (real)(nSub*nSub);
    
    // バッファクリア
    for (int y=0; y<h; y++) {
        for (int x=0; x<w; x++) {
            pColorBuf[w*y + x] = Vec3();
        }
    }
    
    Ray camRay = Ray(pConfig_->camera.position, pConfig_->camera.direction);
    real fovY = pConfig_->camera.fovY;
    
    // 投影面のXY軸
    const Vec3 proj_plane_axis_x = Vec3(w * fovY / h, 0.f, 0.f);
    const Vec3 proj_plane_axis_y = (proj_plane_axis_x % camRay.d).normalize() * fovY;
    
    // Loop over image rows
    #pragma omp parallel for num_threads(4) schedule(dynamic, 1)
    for (int y=0; y<h; y++) {
        #pragma omp critical
        {
            fprintf(stderr, "RayTracing (%d spp) %5.2f%%\n", nSub*nSub, 100.f * y / (h-1));
        }

        xi_[0] = 0;
        xi_[1] = 0;
        xi_[2] = y*y*y;
        
        // Loop cols
        for (unsigned short x=0; x<w; x++) {
            
            int i = (h-y-1) * w + x; // カラーバッファのインデックス
            for (int sy=0; sy<nSub; sy++) {         // subpixel rows
                for (int sx=0; sx<nSub; sx++) {     // subpixel cols
                    
                    // r1, r2 = 0 to 2
                    // dx, dy = -1 to 1  中心に集まったサンプリング --> tent filter
#if USE_TENT_FILTER
                    real r1 = 2*erand48(xi_), dx = (r1 < 1) ? sqrtf(r1)-1 : 1-sqrtf(2-r1);
                    real r2 = 2*erand48(xi_), dy = (r2 < 1) ? sqrtf(r2)-1 : 1-sqrtf(2-r2);
#else
                    real dx = 0;
                    real dy = 0;
#endif
                    // sx = 0 or 1  (...nSub == 2の場合)
                    // dx = -1 to 1
                    // (sx+.5 + dx)/2 --> .5でサブピクセルの中心に。dxでフィルタの揺らぎ。
                    // sx+.5 = 0.5 or 1.5
                    // sx+.5 + dx = -0.5 to 1.5 or 0.5 to 2.5
                    // (sx+.5 + dx)/2 = -0.25 to 0.75 or 0.25 to 1.25、前者0.25中心に+-0.5範囲の揺らぎ、後者0.75中心に+-0.5範囲のゆらぎ
                    // +x / w ピクセルの位置へ、0to1へ。 -.5 0to1から-0.5to0.5へ。
                    // cx* 投影面中心座標から+-0.5の範囲を走査するため
                    // / 2.fはsx, dxが倍の
                    //Vec3 d = cx*( ( (sx+.5f + dx)/2.f + x)/w - .5f) +
                    //        cy*( ( (sy+.5f + dy)/2.f + y)/h - .5f) + camRay.d;
                    const float toCenter = 0.5f;
                    const float sx2 = (sx+toCenter + dx) / (float)nSub;
                    const float sy2 = (sy+toCenter + dy) / (float)nSub;
                    Vec3 d = proj_plane_axis_x * ((x + sx2) / w - .5f)
                           + proj_plane_axis_y * ((y + sy2) / h - .5f)
                           + camRay.d;
                    
                    // 140は視点から投影面までの距離
                    Vec3 r = Irradiance(Ray(camRay.o + d * 140, d.normalize()), 0);

                    // Camera rays are pushed ^^^^^ forward to start in interior
                    // トーンマップとか特にやってない。クランプしてるだけ。
                    pColorBuf[i] += r * subPixelFactor;
                }
            }
            //fprintf(stderr, "\r%f %f %f", c[i].x, c[i].y, c[i].z);
        }
    }
    
}

void PhotonMapRenderer::Run(Vec3* pColorBuf, const Scene& scene, BVH* pBVH)
{
    xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)pPmConfig_->nPhotons;
    
    pScene_ = &scene;
    pBVH_ = pBVH;

    // すべてのライトの合計の明るさを求める
    u32 nLit = (u32)pScene_->litSrcs_.size();
    double sumIntensity = 0;
    for (int i=0; i<nLit; i++) {
        const Vec3& intensity = pScene_->litSrcs_[i]->intensity_;
        sumIntensity += intensity.sum(); // sumでなく輝度を使った方が精度が上がる
    }
    
    // 各ライトからライトの明るさに応じてフォトンをばらまく
    u32 iPhoton = 0;
    for (int i=0; i<nLit; i++) {
        const LightSource* pLit = pScene_->litSrcs_[i];
        float nPhotonRatio = (float)(pLit->intensity_.sum() / sumIntensity);
        u32 nPhotons = (u32)(pPmConfig_->nPhotons * nPhotonRatio);
        #pragma omp parallel for num_threads(4) schedule(dynamic, 1)
        for (int j=0; j<nPhotons; j++, iPhoton++) {
            if (iPhoton % (pPmConfig_->nPhotons / 100) == 0) {
                #pragma omp critical
                {
                    fprintf(stderr, "PhotonTracing %5.2f%%\n", 100. * iPhoton / pPmConfig_->nPhotons);
                }
            }
            
            Ray ray = pLit->GenerateRay();
            Vec3 power = pLit->intensity_;
            PhotonTracing(ray, power.e, 0);
        }
        pPhotonMap_->scale_photon_power(1.0f / nPhotons); // 前回スケールした範囲は除外される
    }
    pPhotonMap_->balance();
    
    RayTracing(pColorBuf);
}
