﻿//siren

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
#include "Material.h"
#include "Random.h"
#include "BmpFileView.h"
#include "StringUtils.h"
#ifdef _WIN32
#include "Windows.h"
#endif


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

Vec3 PhotonMapRenderer2::GlossyRay(const Vec3& w, float exponent, Random& rand)
{
    real r1 = 2.f*(real)(PI*rand.F32());
    real r2 = rand.F32();
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

void PhotonMapRenderer2::TracePhoton(const Ray& r, const Vec3& power, PathInfo pathInfo, Random& rand)
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
    Refl_t refl = rec.pMaterial->refl;
    
    // 影フォトン
    if (traceFlag_ == Trace_Shadow) {
        
        // ここを通るのは必ず直接光
        assert(pathInfo.diffuseDepth == 0);
        
        // 直接光をストア
        pCurrPm_->store(power.e, x.e, r.d.e, true, lightNo_);
        
        // 影フォトンをストア
        std::vector<HitRecord> shadyHits(8);
        int nHits = pScene_->RayCast(shadyHits, 0, r, rec.t+EPSILON, REAL_MAX);
        //printf("%d\n", nHits);
        for (int i=0; i<nHits; i++) {
            HitRecord& sh_rec = shadyHits.at(i);
            if (sh_rec.pMaterial->refl == DIFF && sh_rec.normal.dot(r.d) < 0) {
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
            if (rand.F32() < ave_refl)
            {
                Vec3 d = Ray::CosRay(nl, rand);
                real ave_refl_inv = 1.0f / ave_refl;
                Vec3 refPower = power.mult(color) * ave_refl_inv;
                TracePhoton(Ray(x,d), refPower, pathInfo, rand);
            }
        }
        
        return;
    }
    // Ideal SPECULAR reflection
    else if (refl == SPEC) {
        pathInfo.specularDepth++;
        Ray refl(x, r.d - n * 2.f * n.dot(r.d));
        TracePhoton(refl, power, pathInfo, rand);
        return;
    }
    else if (refl == PHONGMETAL) {
        pathInfo.glossyDepth++;
        
        Vec3 rdir = r.d - n * 2.f * n.dot(r.d); // reflected ray
        Vec3 d = GlossyRay(rdir, 100, rand);
        Vec3 mx = x + n*1e-4f;
        Vec3 refPower = power.mult(color);
        TracePhoton(Ray(mx, d), refPower, pathInfo, rand);
        return;
    }
    
    pathInfo.refractionDepth++;
    
    // Ideal dielectric REFRACTION
    Ray reflRay(x, r.d - n * 2.f * n.dot(r.d));
    bool into = n.dot(nl) > 0.f;  // Ray from outside going in?
    // @todo airRefrIdxを現在レイがある物体のrefrIdxにしてpathInfoに含める
    real airRefrIdx = 1.f;
    real refrIdx = rec.pMaterial->refractiveIndex;
    real nnt = into ? airRefrIdx/refrIdx : refrIdx/airRefrIdx;
    real ddn = r.d.dot(nl);   // レイと法線のcos
    real cos2t;
    
    // Total internal reflection
    if ((cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn)) < 0.f) {
        TracePhoton(reflRay, power, pathInfo, rand);
        return;
    }
    
    // 屈折方向
    Vec3 tdir = (r.d * nnt - n * ((into ? 1.f : -1.f) * (ddn * nnt + sqrtf(cos2t)))).normalize();
    real a = refrIdx - airRefrIdx;
    real b = refrIdx + airRefrIdx;
    real R0 = a * a / (b * b);
    real c = 1.f - (into ? -ddn : tdir.dot(n));
    real Re = R0 + (1.f - R0)*c*c*c*c*c;
    real P = Re;
    // @todo ロシアンルーレット使う.青い本P79
    if (rand.F32() < P)
        TracePhoton(reflRay, power, pathInfo, rand);       // 反射
    else
        TracePhoton(Ray(x, tdir), power.mult(color), pathInfo, rand);  // 屈折

    return;
}

Vec3 PhotonMapRenderer2::Irradiance(const Ray &r, PathInfo pathInfo, Random& rand)
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
    Refl_t refl = rec.pMaterial->refl;
    
    // Ideal DIFFUSE reflection
    if (refl == DIFF) {
        pathInfo.diffuseDepth++;
        Vec3 irrad;
        const float BRDF = PI_INV;
        
        float peRatio = .5f; // means penumbra area
        if (pPmRenConf_->shadowEstimate) {
            const float peDist = pShadowPmConf_->estimateDist;
            const int peNum = pShadowPmConf_->nEstimatePhotons;
            peRatio = pShadowPhotonMap_->penumbra_estimate(lightNo_, x.e, nl, peDist, peNum);
        }
        
        //printf("pe=%f\n", peRatio);
        if (pPmRenConf_->drawShadowEstimate)
        {
            // 影度合い別に可視化
            if (peRatio == 1.f) {
                // direct
                return Vec3(1.f, 1.f, 1.f);
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
            for (u32 i=0; i<pScene_->GetLightNum(); i++) {
                irrad += BRDF * pScene_->GetLight(i)->DirectLight(x, nl, *pScene_, peRatio, rand);
                assert(irrad.e[0] < 0);
            }
        }
        
        // Indirect Light
        if (pPmRenConf_->indirectLight) {
            Vec3 tmpIrrad;
            pPhotonMap_->irradiance_estimate(tmpIrrad.e, x.e, nl, pPmConf_->estimateDist, pPmConf_->nEstimatePhotons);
            irrad += tmpIrrad;
        }
        
        // Caustics
        if (pPmRenConf_->caustic) {
            Vec3 tmpIrrad;
            pCausticPhotonMap_->irradiance_estimate(tmpIrrad.e, x.e, nl, pCausticPmConf_->estimateDist, pCausticPmConf_->nEstimatePhotons);
            irrad += tmpIrrad;
        }
        
        return Vec3(irrad.x * color.x, irrad.y * color.y, irrad.z * color.z);
    }
    else if (refl == SPEC) {
        // Ideal SPECULAR reflection

        return Irradiance(Ray(x,r.d-n*2.f*n.dot(r.d)), pathInfo, rand);
    }
    else if (refl == PHONGMETAL) {
        // Imperfect SPECULAR reflection
        
        // 指数的に追跡回数が増えるのを防ぐ
        if (pathInfo.glossyDepth >= pPmRenConf_->nMaxGlossyBounce) {
            // Ideal SPECULAR reflectionで近似
            return Irradiance(Ray(x,r.d-n*2.f*n.dot(r.d)), pathInfo, rand).mult(color);
        }
        pathInfo.glossyDepth++;
        
        Vec3 irrad;
        const u32 nSamples = pPmRenConf_->nGlossyRays;
        for (u32 i = 0; i < nSamples; i++) {
            
            Vec3 rdir = r.d - n * 2.f * n.dot(r.d); // reflected ray
            rdir = GlossyRay(rdir, 100, rand); // @todo exponentをMaterialに
            Vec3 mx = x + n*1e-4f; // 自己ヒットしないようにちょっと浮かす
            Vec3 tmp = Irradiance(Ray(mx, rdir), pathInfo, rand);
            irrad += tmp;
        }
        return irrad.mult(color) / (float)nSamples;
    }
    else if (refl == REFR) {
    
        // Ideal dielectric REFRACTION
        Vec3 rdir = r.d - n * 2.f * n.dot(r.d);
        //Vec3 mx = x + rdir * EPSILON; // 自己ヒット抑止にレイ始点をちょっと進めてみる
        Ray reflRay(x, rdir);
        bool into = n.dot(nl) > 0.f; // Ray from outside going in?
        real airRefrIdx = 1.f;
        real refrIdx = rec.pMaterial->refractiveIndex;
        real nnt = into ? airRefrIdx/refrIdx : refrIdx/airRefrIdx;
        real ddn = r.d.dot(nl); // レイと法線のcos
        real cos2t;
        
        // Total internal reflection
        if ((cos2t = 1.f-nnt * nnt * (1.f - ddn * ddn)) < 0.f)
            return Irradiance(reflRay, pathInfo, rand);
        
        // 屈折方向
        Vec3 tdir = (r.d * nnt - n * ((into ? 1.f : -1.f) * (ddn * nnt + sqrtf(cos2t)))).normalize();
        //mx = x + rdir * EPSILON; // 自己ヒット抑止にレイ始点をちょっと進めてみる
        
        real a = refrIdx - airRefrIdx;
        real b = refrIdx + airRefrIdx;
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
        return Irradiance(reflRay, pathInfo, rand) * Re
             + Irradiance(Ray(x,tdir), pathInfo2, rand).mult(color) * Tr;
    }
    
    // refl == LIGHT
    return ((AreaLightShape*)rec.pShape)->SelfIrradiance();
}

void PhotonMapRenderer2::PhotonTracing()
{
    // すべてのライトの合計の明るさを求める
    u32 nLit = (u32)pScene_->GetLightNum();
    double sumFlux = 0;
    for (u32 i=0; i<nLit; i++) {
        const Vec3& flux = pScene_->GetLight(i)->GetFlux();
        sumFlux += flux.sum(); // sumでなく輝度を使った方が精度が上がる
    }
    
    // 各ライトからライトの明るさに応じてフォトンをばらまく
    if (pPmConf_->enable && pPmRenConf_->indirectLight) {
        printf("<Indirect Photon Map>\n");
        PhotonTracing_(*pPhotonMap_, *pPmConf_, nLit, sumFlux, Trace_Indirect);
    }
#ifdef _WIN32
    Sleep(1000);
#endif
    if (pCausticPmConf_->enable && pPmRenConf_->caustic) {
        printf("<Caustic Photon Map>\n");
        PhotonTracing_(*pCausticPhotonMap_, *pCausticPmConf_, nLit, sumFlux, Trace_Caustic);
    }
#ifdef _WIN32
    Sleep(1000);
#endif
    if (pShadowPmConf_->enable && pPmRenConf_->shadowEstimate) {
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
    for (int i=0; i<nLit; i++) {
        u32 iPhoton = 0;

        const LightSource* pLit = pScene_->GetLight(i);
        float nPhotonRatio = (float)(pLit->GetFlux().sum() / sumFlux);
        u32 nPhotons = (u32)(c_nPhotons * nPhotonRatio);
        
        int nPhotonsPerThread = c_nPhotonsPerThread > 0 ? c_nPhotonsPerThread : nPhotons;
        int nThread = (int)ceilf(nPhotons / (float)nPhotonsPerThread);
        Random* pRands = new Random[nThread];

        #pragma omp flush
        #pragma omp parallel for num_threads(8) schedule(dynamic, 1) shared(iPhoton)
        //#pragma omp parallel for num_threads(4) schedule(dynamic, 1) shared(iPhoton)
        for (int t=0; t<nThread; t++) {
            pRands[i].SetSeedW(i);
            
            int nPhotonThisThread = (t == nThread-1) ?
                (nPhotons-((nThread-1)*nPhotonsPerThread)) : nPhotonsPerThread;
                 
            for (int j=0; j<nPhotonThisThread; j++) {
                
                #pragma omp atomic
                iPhoton++;
                
                /*
                 if (iPhoton % (c_nPhotons / 100) == 0) {
                    #pragma omp critical
                    {
                        fprintf(stderr, "PhotonTracing %5.2f%%\n", 100. * iPhoton / c_nPhotons);
                    }
                }
                 */
                
                Ray ray = pLit->GenerateRay(pRands[i]);
                Vec3 power = pLit->GetFlux();
                PathInfo pathInfo;
                lightNo_ = i;
                TracePhoton(ray, power, pathInfo, pRands[t]);
            }
        }
        #pragma omp flush
        
        printf("%d photons traced.\n", iPhoton);
        
        photonMap.scale_photon_power(1.0f / nPhotons); // 前回スケールした範囲は除外される
        
        delete[] pRands;
    }
    photonMap.balance();
}

// @todo Appのと共通化, ImageBufferクラスとか導入してそこに移すか？
inline real clamp(real x)
{
    return x < 0.f ? 0.f : x > 1.f ? 1.f : x;
}

// @todo Appのと共通化
inline int toInt(real x)
{
#ifdef USE_FLOAT
    return int(powf(clamp(x), 1/2.2f) * 255.f + .5f);
#else
    return int(pow(clamp(x), 1/2.2) * 255 + .5);
#endif
}

// @todo Appのと共通化
void ConvertToUint(u8* pColorBuf, const Vec3* pRealColorBuf, int w, int h)
{
    // Convert to u8 format
    for (int i = 0, j=0; i < (w*h); ++i, j+=4)
    {
        pColorBuf[j+0] = (u8)(toInt(pRealColorBuf[i].x)); // including Gamma Correction
        pColorBuf[j+1] = (u8)(toInt(pRealColorBuf[i].y));
        pColorBuf[j+2] = (u8)(toInt(pRealColorBuf[i].z));
        pColorBuf[j+3] = 255;
    }
}

void OutputInterimImage(Vec3* pRealColorBuf, int w, int h, const char* filePath)
{
    u8* pColorBuf = new u8[w * h * sizeof(char)*4]; 
    ConvertToUint(pColorBuf, pRealColorBuf, w, h);
    BmpFileView bmpFileView;
    bmpFileView.SetFilePath(filePath);
    bmpFileView.Init(w, h);
    bmpFileView.Present(pColorBuf);

    delete pColorBuf;
}

void PhotonMapRenderer2::RayTracing(Vec3* pColorBuf)
{
    const u32 w = pConf_->windowWidth;
    const u32 h = pConf_->windowHeight;
    const u32 nSub = pPmRenConf_->nSubPixelsSqrt;
    //const real subPixelFactor = 1.0f / (real)(nSub*nSub);
    
    // バッファクリア
    for (u32 y=0; y<h; y++) {
        for (u32 x=0; x<w; x++) {
            pColorBuf[w*y + x] = Vec3();
        }
    }
    
    const Ray camRay = Ray(pConf_->camera.position, pConf_->camera.direction);
    const real fovY = pConf_->camera.fovY;
    
    // 投影面のXY軸
    const Vec3 proj_plane_axis_x = Vec3(w * fovY / h, 0.f, 0.f);
    const Vec3 proj_plane_axis_y = (proj_plane_axis_x % camRay.d).normalize() * fovY;
    
    Random* pRands = new Random[h];
    
    Timer timer;

    int nImage = 0;
    const int nSub2 = nSub*nSub;
    
    for (u32 sy=0; sy<nSub; sy++) {         // subpixel rows
        for (u32 sx=0; sx<nSub; sx++) {     // subpixel cols
                
            const int iSpp = (sy*nSub+sx);

            //#pragma omp parallel for num_threads(4) schedule(dynamic, 1)
            #pragma omp parallel for num_threads(8) schedule(dynamic, 1)
            for (int y=0; y<(int)h; y++) {
                Random& rand = pRands[y];
                rand.SetSeedW(y);
                
                #pragma omp critical
                if (y % nSub*nSub == 0)
                {
                    float progress = 100.f * ((iSpp / (float)nSub2) + (float)y / (h-1) / nSub2);
                    fprintf(stderr, "RayTracing (%d/%d spp) %5.2f%%\n", iSpp+1, nSub2, progress);
                    if (timer.Elapsed() > 60000) {
                        timer.Restart();
                        char pPath[64];
                        StringUtils::Sprintf(pPath, 64, "./image%d.bmp", nImage);
                        OutputInterimImage(pColorBuf, w, h, pPath);
                        nImage++;
                   }
                }
                
                // Loop cols
                for (unsigned short x=0; x<w; x++) {
                    
                    int i = (h-y-1) * w + x; // カラーバッファのインデックス

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
                        real r1 = rand.F32();
                        real r2 = rand.F32();
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
                    Vec3 r = Irradiance(Ray(camRay.o + d * 140, d.normalize()), pathInfo, rand);

                    // Camera rays are pushed ^^^^^ forward to start in interior
                    // トーンマップとか特にやってない。クランプしてるだけ。
                    float sppRatio = 1.f / (iSpp+1);
                    pColorBuf[i] = r * sppRatio + pColorBuf[i] * (1-sppRatio);
                }
            }
            //fprintf(stderr, "\r%f %f %f", c[i].x, c[i].y, c[i].z);
        }
    }

    delete[] pRands;
}

void PhotonMapRenderer2::Run(Vec3* pColorBuf, const Scene& scene)
{
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
