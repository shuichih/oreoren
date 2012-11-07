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
#include "Timer.h"


// Direct
// Coarstic  LSD+E
// Indirect  LD(S|D)*DE
// Shadow
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
    
    return false;
}

Vec3 PhotonMapRenderer::GlossyRay(const Vec3& w, float exponent)
{
    real r1 = 2.f*(real)(M_PI*erand48(xi_));
    real r2 = (real)erand48(xi_);
    real cosTheta = powf(1.f-r2, 1.0f/(exponent+1.0f));
    real sinTheta = sqrtf(1.0f - cosTheta*cosTheta);
    real rx = cosf(r1) * sinTheta;
    real ry = sinf(r1) * sinTheta;
    real rz = cosTheta;
    Vec3 u = ((fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) % w).normalize(); // binormal
    Vec3 v = w % u; // tangent
    
    // ucosφsinθ + vsinφsinθ + wcosθ
    return (u*rx + v*ry + w*rz).normalize();
}
        

void PhotonMapRenderer::TracePhoton(const Ray& r, const Vec3& power, PathInfo& pathInfo)
{
    // max refl
    if (++pathInfo.depth > pPmConfig_->maxPhotonBounce)
    {
        return;
    }

    HitRecord rec;
    rec.hitLit = false;
    if (!Intersect(r, rec))
        return;
    
    Vec3 x = r.o + r.d * rec.t;
    Vec3 n = rec.normal;
    Vec3 nl = n.dot(r.d) < 0.f ? n : n*-1.f; // 交点の法線, 裏面ヒット考慮
    Vec3 color = rec.color;
    Refl_t refl = rec.refl;
    
    // Ideal DIFFUSE reflection
    if (refl == DIFF || refl == LIGHT) { // 光源表面はLambert面という事にしておく
        pathInfo.diffuseDepth++;
        // 直接光のフォトンはストアしない(集光模様フォトン除く)
//        if (pathInfo.diffuseDepth > 1 || (pathInfo.specularDepth + pathInfo.glossyDepth) > 0) {
//        if (pathInfo.diffuseDepth > 1 && (pathInfo.specularDepth + pathInfo.glossyDepth) == 0) {
        pPhotonMap_->store(power.e, x.e, r.d.e);
//        }
        
        float ave_refl = (color.x + color.y + color.z) / 3.f;
        if ((real)erand48(xi_) < ave_refl)
        {
            // 法線となす角のcosに比例する分布で反射方向を決める
            Vec3 d = Ray::CosRay(nl, xi_);
            
            real ave_refl_inv = 1.0f / ave_refl;
            Vec3 refPower = power.mult(color) * PI_INV * ave_refl_inv;
            TracePhoton(Ray(x,d), refPower, pathInfo);
        }
        //fprintf(stderr, "p (%f %f %f) (%f %f %f) (%f %f %f)\n", power[0], power[1], power[2], pos[0], pos[1], pos[2], dir[0], dir[1], dir[2]);
        return;
    }
    // Ideal SPECULAR reflection
    else if (refl == SPEC) {
        pathInfo.specularDepth++;
        Ray refl(x, r.d - n * 2.f * n.dot(r.d));
        TracePhoton(refl, power, pathInfo);
        return;
    }
    else if (refl == PHONGMETAL) {
        pathInfo.glossyDepth++;
        
        Vec3 rdir = r.d - n * 2.f * n.dot(r.d); // reflected ray
        Vec3 d = GlossyRay(rdir, 10);
        Vec3 mx = x + n*1e-4f;
        TracePhoton(Ray(mx, d), power, pathInfo);
        return;
    }
    
    pathInfo.specularDepth++;
    
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
        TracePhoton(reflRay, power, pathInfo);
        return;
    }
    
    // 屈折方向
    Vec3 tdir = (r.d * nnt - n * ((into ? 1.f : -1.f) * (ddn * nnt + sqrtf(cos2t)))).normalize();
    real a = grassRefrIdx - airRefrIdx;
    real b = grassRefrIdx + airRefrIdx;
    real R0 = a * a / (b * b); // 垂直反射率 0.25 / 2.5 = 0.1 @todo 根拠調査
    real c = 1.f - (into ? -ddn : tdir.dot(n));
    real fresnel = R0 + (1.f - R0)*c*c*c*c*c;
    real P = fresnel;
    if ((real)erand48(xi_) < P)
        TracePhoton(reflRay, power, pathInfo);       // 反射
    else
        TracePhoton(Ray(x, tdir), power.mult(color), pathInfo);  // 屈折

    return;
}


Vec3 PhotonMapRenderer::Irradiance(const Ray &r, PathInfo& pathInfo)
{
    // max refl
    if (++pathInfo.depth > pPmConfig_->maxRayBounce)
    {
        return Vec3();
    }

    HitRecord rec;
    rec.hitLit = true;
    if (!Intersect(r, rec)) return Vec3();
    //const Shape& obj = *g_shapes[id];       // the hit object
    Vec3 x = r.o + r.d * rec.t;
    Vec3 n = rec.normal;
    Vec3 nl = n.dot(r.d) < 0.f ? n : n * -1.f;   // 交点の法線
    Vec3 color = rec.color;
    Refl_t refl = rec.refl;
    
    // 0.5にしたらカラーが反射率になってるから暗くなるだけ。IDEALでない反射は扱えない。カラーと混ぜるとかもない。
    // Ideal DIFFUSE reflection
    if (refl == DIFF) {
        pathInfo.diffuseDepth++;
        Vec3 irrad;
        //fprintf(stderr, "irrad %f %f %f\r", irrad[0], irrad[1], irrad[2]);
        if (pPmConfig_->finalGethering && pathInfo.diffuseDepth <= 1) {
            for (int i=0; i<pPmConfig_->nFinalGetheringRays; i++) {
                real r1 = 2.f*(real)(M_PI*erand48(xi_));
                real r2 = (real)erand48(xi_); // => 1-cos^2θ = 1-sqrt(1-r_2)^2 = r_2
                real r2s = sqrtf(r2);    // => sinθ = sqrt(1-cos^2θ) = sqrt(r_2)
                Vec3 w = nl; // normal
                Vec3 u = ((fabs(w.x) > .1f ? Vec3(0.f, 1.f, 0.f) : Vec3(1.f, 0.f, 0.f)) % w).normalize(); // binormal
                Vec3 v = w % u; // tangent
                
                // ucosφsinθ + vsinφsinθ + wcosθ
                Vec3 d = (u*cosf(r1)*r2s + v*sinf(r1)*r2s + w*sqrtf(1-r2)).normalize();
                PathInfo pathInfo2(pathInfo);
                irrad += Irradiance(Ray(x, d), pathInfo2);// * nl.dot(d) / nl.dot(d);
                // 入射方向と法線のcosθ掛けるのとimportance samplingしたので確率密度関数cosθで割るのとで
                // 相殺するような。
            }
            irrad /= pPmConfig_->nFinalGetheringRays;
        } else {
#if 0
            //if (pathInfo.diffuseDepth == 1) {
            // @todo 複数ライト対応
            Vec3 ldir = pScene_->litSrcs_[0]->position_ - x;
            float r2 = ldir.square_length();
            ldir.normalize();
            irrad = nl.dot(ldir) * (pScene_->litSrcs_[0]->GetIntensity() / (4 * PI * r2));
            
            //} else {
                //pPhotonMap_->irradiance_estimate(irrad.e, x.e, nl, pPmConfig_->estimateDist, pPmConfig_->nEstimatePhotons);
            //}
#else
            pPhotonMap_->irradiance_estimate(irrad.e, x.e, nl, pPmConfig_->estimateDist, pPmConfig_->nEstimatePhotons);
#endif
            
        }
        return Vec3(irrad.x * color.x, irrad.y * color.y, irrad.z * color.z);
    }
    else if (refl == SPEC) {
        // Ideal SPECULAR reflection
        return Irradiance(Ray(x,r.d-n*2.f*n.dot(r.d)), pathInfo);
    }
    else if (refl == PHONGMETAL) {
        // Imperfect SPECULAR reflection
        
        // 指数的に追跡回数が増えるのを防ぐ
        if (pathInfo.glossyDepth >= pPmConfig_->nMaxGlossyBounce) {
            // Ideal SPECULAR reflectionで近似
            return Irradiance(Ray(x,r.d-n*2.f*n.dot(r.d)), pathInfo);
        }
        pathInfo.glossyDepth++;
        
        Vec3 irrad;
        const u32 nSamples = pPmConfig_->nGlossyRays;
        for (int i = 0; i < nSamples; i++) {
            
            Vec3 rdir = r.d - n * 2.f * n.dot(r.d); // reflected ray
            rdir = GlossyRay(rdir, 10);
            Vec3 mx = x + n*1e-4f; // 自己ヒットしないようにちょっと浮かす
            irrad += Irradiance(Ray(mx, rdir), pathInfo);
        }
        return irrad / nSamples;
    }
    else if (refl == REFR) {
    
        // Ideal dielectric REFRACTION
        Vec3 rdir = r.d - n * 2.f * n.dot(r.d);
        //Vec3 mx = x + rdir * EPSILON; // 自己ヒット抑止にレイ始点をちょっと進めてみる
        Ray reflRay(x, rdir);
        bool into = n.dot(nl) > 0.f; // Ray from outside going in?
        real airRefrIdx = 1.f;
        real grassRefrIdx = 1.5f;
        real nnt = into ? airRefrIdx/grassRefrIdx : grassRefrIdx/airRefrIdx;
        real ddn = r.d.dot(nl); // レイと法線のcos
        real cos2t;
        
        // Total internal reflection
        if ((cos2t = 1.f-nnt * nnt * (1.f - ddn * ddn)) < 0.f)
            return Irradiance(reflRay, pathInfo);
        
        // 屈折方向
        Vec3 tdir = (r.d * nnt - n * ((into ? 1.f : -1.f) * (ddn * nnt + sqrtf(cos2t)))).normalize();
        //mx = x + rdir * EPSILON; // 自己ヒット抑止にレイ始点をちょっと進めてみる
        
        real a = grassRefrIdx - airRefrIdx;
        real b = grassRefrIdx + airRefrIdx;
        real R0 = a * a / (b * b); // 垂直反射率 0.25 / 2.5 = 0.1 @todo 根拠調査
        real c = 1.f - (into ? -ddn : tdir.dot(n));
        real fresnel = R0 + (1.f - R0)*c*c*c*c*c;
        real Tr = 1.f - fresnel;
        // 反射屈折両方トレース
        PathInfo pathInfo2(pathInfo);
        return Irradiance(reflRay, pathInfo) * fresnel
             + Irradiance(Ray(x,tdir), pathInfo2).mult(color) * Tr;
    }
    
    // refl == LIGHT
    return ((AreaLightShape*)rec.pShape)->SelfIrradiance();
}

void PhotonMapRenderer::PhotonTracing()
{
    // すべてのライトの合計の明るさを求める
    u32 nLit = (u32)pScene_->litSrcs_.size();
    double sumIntensity = 0;
    for (int i=0; i<nLit; i++) {
        const Vec3& intensity = pScene_->litSrcs_[i]->GetIntensity();
        sumIntensity += intensity.sum(); // sumでなく輝度を使った方が精度が上がる
    }
    
    // 各ライトからライトの明るさに応じてフォトンをばらまく
    u32 iPhoton = 0;
    for (int i=0; i<nLit; i++) {
        const LightSource* pLit = pScene_->litSrcs_[i];
        float nPhotonRatio = (float)(pLit->GetIntensity().sum() / sumIntensity);
        u32 nPhotons = (u32)(pPmConfig_->nPhotons * nPhotonRatio);
        
        const int nPhotonsPerThread =
            pPmConfig_->nTracePhotonsPerThread > 0 ? pPmConfig_->nTracePhotonsPerThread : nPhotons;
        int nThread = ceilf(nPhotons / (float)nPhotonsPerThread);
        
        #pragma omp parallel for num_threads(4) schedule(dynamic, 1)
        for (int t=0; t<nThread; t++) {
            int nPhotonThisThread = (t == nThread-1) ?
                (nPhotons-((nThread-1)*nPhotonsPerThread)) : nPhotonsPerThread;
                 
            for (int j=0; j<nPhotonThisThread; j++) {
                #pragma omp atomic
                iPhoton++;
                
                if (iPhoton % (pPmConfig_->nPhotons / 100) == 0) {
                    #pragma omp critical
                    {
                        fprintf(stderr, "PhotonTracing %5.2f%%\n", 100. * iPhoton / pPmConfig_->nPhotons);
                    }
                }
                
                
                Ray ray = pLit->GenerateRay();
                Vec3 power = pLit->GetIntensity();
                PathInfo pathInfo;
                TracePhoton(ray, power, pathInfo);
            }
        }
        
        printf("%d photons traced.\n", iPhoton);
        
        pPhotonMap_->scale_photon_power(1.0f / nPhotons); // 前回スケールした範囲は除外される
    }
    pPhotonMap_->balance();
    
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
                    real dx = 0;
                    real dy = 0;
                    if (pPmConfig_->useTentFilter) {
                        real r1 = (float)(2*erand48(xi_));
                        real r2 = (float)(2*erand48(xi_));
                        dx = (r1 < 1) ? sqrtf(r1)-1 : 1-sqrtf(2-r1);
                        dy = (r2 < 1) ? sqrtf(r2)-1 : 1-sqrtf(2-r2);
                    } else {
                        dx = 0;
                        dy = 0;
                    }
                    
                    const float toCenter = 0.5f;
                    const float sx2 = (sx+toCenter + dx) / (float)nSub;
                    const float sy2 = (sy+toCenter + dy) / (float)nSub;
                    Vec3 d = proj_plane_axis_x * ((x + sx2) / w - .5f)
                           + proj_plane_axis_y * ((y + sy2) / h - .5f)
                           + camRay.d;
                    
                    // 140は視点から投影面までの距離
                    PathInfo pathInfo;
                    Vec3 r = Irradiance(Ray(camRay.o + d * 140, d.normalize()), pathInfo);

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
    
    {
        Timer timer;
        PhotonTracing();
        timer.PrintElapsed("PhotonTracing time: ");
    }
    {
        Timer timer;
        RayTracing(pColorBuf);
        timer.PrintElapsed("RayTracing time: ");
    }
}
