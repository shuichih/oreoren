//siren

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "Common.h"
#include "Scene.h"
#include "PhotonMap.h"
#include "LightSource.h"
#include "PhotonMapRenderer2.h"
#include "PhotonFilter.h"
#include "Config.h"
#include "Timer.h"
#include "Ray.h"


// Direct
// Caustic  LSD+E
// Indirect  LD(S|D)*DE
// Shadow
using namespace std;

//----------------------------------------------------------------
PhotonMapRenderer2::PhotonMapRenderer2()
: pCurrPm_(NULL)
, pPhotonMap_(NULL)
, pCausticPhotonMap_(NULL)
, pShadowPhotonMap_(NULL)
, pFilter_(NULL)
, pCausticFilter_(NULL)
{
}

PhotonMapRenderer2::~PhotonMapRenderer2()
{
    delete pPhotonMap_;
    delete pFilter_;
}

void PhotonMapRenderer2::SetConfig(const Config& config)
{
    pConf_ = &config;
    pPmRenConf_ = &pConf_->pmRendererConf;
    pPmConf_ = &pConf_->photonMapConf;
    pCausticPmConf_ = &pConf_->causticPmConf;
    pShadowPmConf_ = &pConf_->shadowPmConf;
    
    InitializePhotonMap(&pPhotonMap_, *pPmConf_, &pFilter_);
    InitializePhotonMap(&pCausticPhotonMap_, *pCausticPmConf_, &pCausticFilter_);
    InitializePhotonMap(&pShadowPhotonMap_, *pShadowPmConf_, NULL);
}

void PhotonMapRenderer2::InitializePhotonMap(Photon_map** ppPm, const PhotonMapConfig& pmConf, PhotonFilter** ppFilter)
{
    if (*ppPm != NULL) delete *ppPm;
    *ppPm = new Photon_map(pmConf.nMaxStorePhotons);
    
    if (ppFilter) {
        if (*ppFilter) delete *ppFilter;
        
        if (pmConf.enableConeFilter) {
            *ppFilter = new ConeFilter(pmConf.coneFilterK);
            ((ConeFilter*)(*ppFilter))->SetK(pmConf.coneFilterK);
            (*ppPm)->SetFilter(*ppFilter);
        }
    }
    else {
        (*ppPm)->SetFilter(NULL);
    }
    (*ppPm)->SetEstimateEllipseScale(pmConf.estimateEllipseScale);
}
    
Vec3 PhotonMapRenderer2::GlossyRay(const Vec3& w, float exponent)
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

void PhotonMapRenderer2::TracePhoton(const Ray& r, const Vec3& power, PathInfo& pathInfo)
{
    if (++pathInfo.depth > pPmConf_->maxPhotonBounce)
    {
        return;
    }

    HitRecord rec;
    rec.hitLit = false;
    if (!pScene_->Intersect(r, EPSILON, REAL_MAX, rec))
        return;
    
    Vec3 x = r.o + r.d * rec.t;
    Vec3 n = rec.normal;
    Vec3 nl = n.dot(r.d) < 0.f ? n : n*-1.f; // 交点の法線, 裏面ヒット考慮
    Vec3 color = rec.color;
    Refl_t refl = rec.refl;
    
    // 影フォトン
    if (traceFlag_ == Trace_Shadow) {
        
        // ここを通るのは必ず直接光
        assert(pathInfo.diffuseDepth == 0);
        
        // 直接光をストア
        pCurrPm_->store(power.e, x.e, r.d.e, true, lightNo_);
        
        // 影フォトンをストア
        int nHits = pScene_->RayCast(shadyHits_, 0, r, rec.t+EPSILON, REAL_MAX);
        //printf("%d\n", nHits);
        for (int i=0; i<nHits; i++) {
            HitRecord& sh_rec = shadyHits_.at(i);
            if (sh_rec.refl == DIFF && sh_rec.normal.dot(r.d) < 0) {
                //printf("hit\n");
                Vec3 sh_x = r.o + r.d * sh_rec.t;
                pCurrPm_->store((-1.f*power).e, sh_x.e, r.d.e, false, lightNo_);
            }
        }
        
        // 再帰せずここで終了
        return;
    }
    
    // Ideal DIFFUSE reflection
    if (refl == DIFF || refl == LIGHT) { // 光源表面はLambert面という事にしておく
        pathInfo.diffuseDepth++;
        bool direct = (pathInfo.diffuseDepth == 1);
        
        if (traceFlag_ == Trace_Caustic) {
            
            // 屈折を経たPhotonのみをストア
            if (pathInfo.refractionDepth > 0)
            {
                pCurrPm_->store(power.e, x.e, r.d.e, true, lightNo_);
            }
            
            // 間接光の集光模様は影響が低そうなのでストアしない
            
        } else {
            if (!direct) {
                // 間接光のフォトンをストア
                pCurrPm_->store(power.e, x.e, r.d.e, false, lightNo_);
            }

            // 拡散反射
            float ave_refl = (color.x + color.y + color.z) / 3.f;
            if ((real)erand48(xi_) < ave_refl)
            {
                Vec3 d = Ray::CosRay(nl, xi_);
                real ave_refl_inv = 1.0f / ave_refl;
                Vec3 refPower = power.mult(color) * ave_refl_inv;
                TracePhoton(Ray(x,d), refPower, pathInfo);
            }
        }
        
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
    
    pathInfo.refractionDepth++;
    
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
    real R0 = a * a / (b * b);
    real c = 1.f - (into ? -ddn : tdir.dot(n));
    real Re = R0 + (1.f - R0)*c*c*c*c*c;
    real P = Re;
    // @todo ロシアンルーレット使う.青い本P79
    if ((real)erand48(xi_) < P)
        TracePhoton(reflRay, power, pathInfo);       // 反射
    else
        TracePhoton(Ray(x, tdir), power.mult(color), pathInfo);  // 屈折

    return;
}

Vec3 PhotonMapRenderer2::Irradiance(const Ray &r, PathInfo& pathInfo)
{
    // max refl
    if (++pathInfo.depth > pPmRenConf_->maxRayBounce)
    {
        return Vec3();
    }

    HitRecord rec;
    rec.hitLit = true;
    if (!pScene_->Intersect(r, EPSILON, REAL_MAX, rec)) return Vec3();
    Vec3 x = r.o + r.d * rec.t;
    Vec3 n = rec.normal;
    Vec3 nl = n.dot(r.d) < 0.f ? n : n * -1.f;   // 交点の法線
    Vec3 color = rec.color;
    Refl_t refl = rec.refl;
    
    // Ideal DIFFUSE reflection
    if (refl == DIFF) {
        pathInfo.diffuseDepth++;
        Vec3 irrad;
        const float BRDF = PI_INV;
        
        const float peDist = pShadowPmConf_->estimateDist;
        const int peNum = pShadowPmConf_->nEstimatePhotons;
        float peRatio = pShadowPhotonMap_->penumbra_estimate(lightNo_, x.e, nl, peDist, peNum);
        
        //printf("pe=%f\n", peRatio);
        if (pPmRenConf_->drawShadowEstimate)
        {
            // 影度合い別に可視化
            if (peRatio == 1.f) {
                // direct
                return Vec3(0, .2f, .5f);
            }
            else if (peRatio == 0) {
                // umbra
                return Vec3(0, 0, 0);
            }
            else {
                // penumbra
                return Vec3(.5f, .5f, .5f);
            }
        }
        
        // Direct Light
        if (pPmRenConf_->directLight) {
            for (int i=0; i<pScene_->litSrcs_.size(); i++) {
                irrad += BRDF * pScene_->litSrcs_[i]->DirectLight(x, nl, *pScene_, peRatio);
            }
        }
        
        // Indirect Light
        if (pPmRenConf_->indirectLight) {
            Vec3 tmpIrrad;
            pPhotonMap_->irradiance_estimate(tmpIrrad.e, x.e, nl, pPmConf_->estimateDist, pPmConf_->nEstimatePhotons);
            irrad += tmpIrrad;
        }
        
        // Caustics
        if (pPmRenConf_->caustics) {
            Vec3 tmpIrrad;
            pCausticPhotonMap_->irradiance_estimate(tmpIrrad.e, x.e, nl, pCausticPmConf_->estimateDist, pCausticPmConf_->nEstimatePhotons);
            irrad += tmpIrrad;
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
        if (pathInfo.glossyDepth >= pPmRenConf_->nMaxGlossyBounce) {
            // Ideal SPECULAR reflectionで近似
            return Irradiance(Ray(x,r.d-n*2.f*n.dot(r.d)), pathInfo);
        }
        pathInfo.glossyDepth++;
        
        Vec3 irrad;
        const u32 nSamples = pPmRenConf_->nGlossyRays;
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
        real R0 = a * a / (b * b);
        real c = 1.f - (into ? -ddn : tdir.dot(n));
        real Re = R0 + (1.f - R0)*c*c*c*c*c;
        // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の2乗の分だけ変化する
        // 屈折率が単位立体角あたりの値だから
        // と思ったけど、この関数ではレイが放射輝度を運んでいるわけではないので、いらないっぽい
        //real nnt2 = nnt * nnt;
        real Tr = (1.f - Re);// * nnt2;
        // 反射屈折両方トレース
        PathInfo pathInfo2(pathInfo);
        return Irradiance(reflRay, pathInfo) * Re
             + Irradiance(Ray(x,tdir), pathInfo2).mult(color) * Tr;
    }
    
    // refl == LIGHT
    return ((AreaLightShape*)rec.pShape)->SelfIrradiance();
}

void PhotonMapRenderer2::PhotonTracing()
{
    // すべてのライトの合計の明るさを求める
    u32 nLit = (u32)pScene_->litSrcs_.size();
    double sumFlux = 0;
    for (int i=0; i<nLit; i++) {
        const Vec3& flux = pScene_->litSrcs_[i]->GetFlux();
        sumFlux += flux.sum(); // sumでなく輝度を使った方が精度が上がる
    }
    
    // 各ライトからライトの明るさに応じてフォトンをばらまく
    if (pPmConf_->enable) {
        printf("<Indirect Photon Map>\n");
        PhotonTracing_(*pPhotonMap_, *pPmConf_, nLit, sumFlux, Trace_Indirect);
    }
    if (pCausticPmConf_->enable) {
        printf("<Caustic Photon Map>\n");
        PhotonTracing_(*pCausticPhotonMap_, *pCausticPmConf_, nLit, sumFlux, Trace_Caustic);
    }
    if (pShadowPmConf_->enable) {
        printf("<Shadow Photon Map>\n");
        PhotonTracing_(*pShadowPhotonMap_, *pShadowPmConf_, nLit, sumFlux, Trace_Shadow);
    }
}

// 各ライトからライトの明るさに応じてフォトンをばらまく
void PhotonMapRenderer2::PhotonTracing_(
    Photon_map& photonMap,
    const PhotonMapConfig& pmConfig,
    int nLit,
    double sumFlux,
    PhotonMapRenderer2::TraceFlag traceFlag)
{
    const int c_nPhotonsPerThread = pPmRenConf_->nTracePhotonsPerThread;
    const int c_nPhotons = pmConfig.nPhotons;
    
    pCurrPm_ = &photonMap;
    traceFlag_ = traceFlag;
    
    // 各ライトからライトの明るさに応じてフォトンをばらまく
    u32 iPhoton = 0;
    for (int i=0; i<nLit; i++) {
        const LightSource* pLit = pScene_->litSrcs_[i];
        float nPhotonRatio = (float)(pLit->GetFlux().sum() / sumFlux);
        u32 nPhotons = (u32)(c_nPhotons * nPhotonRatio);
        
        int nPhotonsPerThread = c_nPhotonsPerThread > 0 ? c_nPhotonsPerThread : nPhotons;
        int nThread = ceilf(nPhotons / (float)nPhotonsPerThread);
        
        #pragma omp flush
        #pragma omp parallel for num_threads(4) schedule(dynamic, 1)
        for (int t=0; t<nThread; t++) {
            int nPhotonThisThread = (t == nThread-1) ?
                (nPhotons-((nThread-1)*nPhotonsPerThread)) : nPhotonsPerThread;
                 
            for (int j=0; j<nPhotonThisThread; j++) {
                
                #pragma omp atomic
                iPhoton++;
                
                if (iPhoton % (c_nPhotons / 100) == 0) {
                    #pragma omp critical
                    {
                        fprintf(stderr, "PhotonTracing %5.2f%%\n", 100. * iPhoton / c_nPhotons);
                    }
                }
                
                Ray ray = pLit->GenerateRay();
                Vec3 power = pLit->GetFlux();
                PathInfo pathInfo;
                lightNo_ = i;
                TracePhoton(ray, power, pathInfo);
            }
        }
        #pragma omp flush
        
        printf("%d photons traced.\n", iPhoton);
        
        photonMap.scale_photon_power(1.0f / nPhotons); // 前回スケールした範囲は除外される
    }
    photonMap.balance();
}



void PhotonMapRenderer2::RayTracing(Vec3* pColorBuf)
{
    const u32 w = pConf_->windowWidth;
    const u32 h = pConf_->windowHeight;
    const u32 nSub = pPmConf_->nSubPixelsSqrt;
    const real subPixelFactor = 1.0f / (real)(nSub*nSub);
    
    // バッファクリア
    for (int y=0; y<h; y++) {
        for (int x=0; x<w; x++) {
            pColorBuf[w*y + x] = Vec3();
        }
    }
    
    Ray camRay = Ray(pConf_->camera.position, pConf_->camera.direction);
    real fovY = pConf_->camera.fovY;
    
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
                    if (pPmRenConf_->useTentFilter) {
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

void PhotonMapRenderer2::Run(Vec3* pColorBuf, const Scene& scene)
{
    xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)pPmConf_->nPhotons;
    
    pScene_ = &scene;
    
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
