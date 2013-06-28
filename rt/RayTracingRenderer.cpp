#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "Common.h"
#include "Scene.h"
#include "LightSource.h"
#include "RayTracingRenderer.h"
#include "BVH.h"
#include "Config.h"
#include "Ray.h"

using namespace std;

//----------------------------------------------------------------
RayTracingRenderer::RayTracingRenderer()
: pBVH_(NULL)
{
}

RayTracingRenderer::~RayTracingRenderer()
{
}

void RayTracingRenderer::SetConfig(const Config& config)
{
    pConfig_ = &config;
    pRtConfig_ = &config.rayTracingConf;
}

Vec3 RayTracingRenderer::Irradiance(const Ray &r, int depth)
{
    // max refl
    if (++depth > pRtConfig_->maxRayBounce)
    {
        return Vec3();
    }
    
    HitRecord rec;
    if (!pScene_->Intersect(r, EPSILON, REAL_MAX, rec)) return Vec3(0, 0, 1);
    //const IShape& obj = *g_shapes[id];       // the hit object
    Vec3 x = r.o + r.d * rec.t;
    Vec3 n = rec.normal;
    Vec3 nl = n.dot(r.d) < 0.f ? n : n * -1.f;   // 交点の法線
    Vec3 f = rec.color;
    Refl_t refl = rec.refl;
    
    
    // 0.5にしたらカラーが反射率になってるから暗くなるだけ。IDEALでない反射は扱えない。カラーと混ぜるとかもない。
    // Ideal DIFFUSE reflection
    if (refl == DIFF){
        Vec3 litPos;
        LightSourceType litType = pScene_->litSrcs_[0]->GetType();
        if (litType == Lit_Point) {
            litPos = ((PointLightSource*)pScene_->litSrcs_[0])->position_;
        } else if (litType == Lit_Area) {
            AreaLightSource* pLitSrc = (AreaLightSource*)pScene_->litSrcs_[0];
            litPos = (pLitSrc->p_[0] + pLitSrc->p_[1] + pLitSrc->p_[2] + pLitSrc->p_[3]) * 0.25f;
        } else if (litType == Lit_Sphere) {
            litPos = ((SphereLightSource*)pScene_->litSrcs_[0])->position_;
        } else {
            litPos = ((SphereLightSource*)pScene_->litSrcs_[0])->position_;
            assert(false);
        }
        Vec3 L = (litPos - x).normalize();
        float cosLN = L.dot(nl);
        //fprintf(stderr, "irrad %f %f %f\r", irrad[0], irrad[1], irrad[2]);
        return f * cosLN;
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
    real R0 = a * a / (b * b);
    real c = 1.f - (into ? -ddn : tdir.dot(n));
    real Re = R0 + (1.f - R0)*c*c*c*c*c;
    // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の2乗の分だけ変化する
    // 屈折率が単位立体角あたりの値だから
    real nnt2 = nnt * nnt;
    real Tr = (1.f - Re) * nnt2;
    // 反射屈折両方トレース
    return Irradiance(reflRay, depth) * Re + Irradiance(Ray(x,tdir), depth).mult(f) * Tr;
}

void RayTracingRenderer::RayTracing(Vec3* pColorBuf)
{
    const u32 w = pConfig_->windowWidth;
    const u32 h = pConfig_->windowHeight;
    const u32 nSub = pRtConfig_->nSubPixelsSqrt;
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
    // 355, 225, 187, 167
    #pragma omp parallel for num_threads(4) schedule(dynamic, 1)   // OpenMP
    for (int y=0; y<h; y++) {
    //int y = h / 2; {
        fprintf(stderr, "RayTracing (%d spp) %5.2f%%\n", nSub*nSub, 100.f * y / (h-1));
        
        xi_[0] = 0;
        xi_[1] = 0;
        xi_[2] = y*y*y;
        
        // Loop cols
        for (unsigned short x=0; x<w; x++) {
        //unsigned short x = w / 2; {
    
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
                    if (pRtConfig_->useTentFilter) {
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
                    
                    Vec3 r = Irradiance(Ray(camRay.o + d * pRtConfig_->distanceToProjPlane, d.normalize()), 0);
                    
                    // Camera rays are pushed ^^^^^ forward to start in interior
                    // トーンマップとか特にやってない。クランプしてるだけ。
                    pColorBuf[i] += r * subPixelFactor;
                }
            }
            //fprintf(stderr, "\r%f %f %f", c[i].x, c[i].y, c[i].z);
         }
    }
    
}

void RayTracingRenderer::Run(Vec3* pColorBuf, const Scene& scene)
{
    xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = pConfig_->windowWidth * pConfig_->windowHeight; // テキトウ
    
    pScene_ = &scene;
    
    RayTracing(pColorBuf);
}
