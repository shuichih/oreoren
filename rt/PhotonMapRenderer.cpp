#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Common.h"
#include "Scene.h"
#include "SceneData.h"
#include "PhotonMap.h"
#include "LightSource.h"
#include "PhotonMapRenderer.h"
#include "PhotonFilter.h"


//----------------------------------------------------------------
PhotonMapRenderer::PhotonMapRenderer()
{
    defaultConfig_.screenWidth = 256;
    defaultConfig_.screenHeight = 256;
    defaultConfig_.nSubPixelsSqrt = 2;
    defaultConfig_.nPhotons = 100000;
    defaultConfig_.nEstimatePhotons = 100;
    defaultConfig_.estimateDist = 10.f;
    defaultConfig_.pFilter = new ConeFilter(1.1f);
    
    config_ = defaultConfig_;
    
    pPhotonMap_ = new Photon_map(config_.nPhotons);
}

PhotonMapRenderer::~PhotonMapRenderer()
{
    delete defaultConfig_.pFilter;
    delete pPhotonMap_;
}

PhotonMapRenderer::Config PhotonMapRenderer::GetDefaultConfig()
{
    return defaultConfig_;
}

void PhotonMapRenderer::SetConfig(const PhotonMapRenderer::Config &config)
{
    config_ = config;
    delete pPhotonMap_;
    pPhotonMap_ = new Photon_map(config_.nPhotons);
}

bool PhotonMapRenderer::Intersect(const Ray& r, HitRecord& out)
{
    real inf = 1e20;
    out.t = inf;
    
    HitRecord rec;
    for (int i=g_nShapes; i--;) {
        if(g_shapes[i]->intersect(r, rec) && (rec.t < out.t)) {
            out = rec;
        }
    }
    
    return out.t < inf;
}

void PhotonMapRenderer::PhotonTracing(const Ray& r, float power[3], int depth)
{
    // max refl
    if (++depth > 5)
    {
        return;
    }

    HitRecord rec;
//    if (!Intersect(r, t, id))
    if (!Intersect(r, rec))
        return;
    
    //const Shape &obj = *g_shapes[id];    // the hit object
    Vec x = r.o + r.d * rec.t;
    Vec n = rec.normal;
    Vec nl = n.dot(r.d) < 0 ? n : n*-1; // 交点の法線, 裏面ヒット考慮
    Vec color = rec.color;
    Refl_t refl = rec.refl;
    
    // Ideal DIFFUSE reflection
    if (refl == DIFF) {
        float pos[3] = { x.x, x.y, x.z };
        float dir[3] = { r.d.x, r.d.y, r.d.z };
        
        pPhotonMap_->store(power, pos, dir);
        
        float ave_refl = (color.x + color.y + color.z) / 3.f;
        if (erand48(xi_) < ave_refl)
        {
            real r1 = 2*M_PI*erand48(xi_);
            real r2 = erand48(xi_); // => 1-cos^2θ = 1-sqrt(1-r_2)^2 = r_2
            real r2s = sqrt(r2);    // => sinθ = sqrt(1-cos^2θ) = sqrt(r_2)
            Vec w = nl; // normal
            Vec u = ((fabs(w.x) > .1 ? Vec(0,1) : Vec(1)) % w).normalize(); // binormal
            Vec v = w % u; // tangent
            
            // ucosφsinθ + vsinφsinθ + wcosθ
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).normalize();

            real ave_refl_inv = 1.0 / ave_refl;
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
        Ray refl(x, r.d - n * 2 * n.dot(r.d));
        PhotonTracing(refl, power, depth);
        return;
    }
    else if (refl == PHONGMETAL) {
        real r1 = 2*M_PI*erand48(xi_);
        real r2 = erand48(xi_);
        real exponent = 10.0;
        real cosTheta = pow(1-r2, 1.0/(exponent+1.0));
        real sinTheta = sqrt(1.0 - cosTheta*cosTheta);
        real rx = cos(r1) * sinTheta;
        real ry = sin(r1) * sinTheta;
        real rz = cosTheta;
        Vec w = r.d - n * 2 * n.dot(r.d); // reflected ray
        Vec u = ((fabs(w.x) > .1 ? Vec(0,1) : Vec(1)) % w).normalize(); // binormal
        Vec v = w % u; // tangent
        
        // ucosφsinθ + vsinφsinθ + wcosθ
        Vec d = (u*rx + v*ry + w*rz).normalize();
        
        Vec mx = x + n*1e-4;
        PhotonTracing(Ray(mx, d), power, depth);
        return;
    }
    
    // Ideal dielectric REFRACTION
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d));
    bool into = n.dot(nl) > 0;  // Ray from outside going in?
    real airRefrIdx = 1;
    real grassRefrIdx = 1.5;
    real nnt = into ? airRefrIdx/grassRefrIdx : grassRefrIdx/airRefrIdx;
    real ddn = r.d.dot(nl);   // レイと法線のcos
    real cos2t;
    
    // Total internal reflection
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) {
        PhotonTracing(reflRay, power, depth);
        return;
    }
    
    // 屈折方向
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
    
    real a = grassRefrIdx - airRefrIdx;
    real b = grassRefrIdx + airRefrIdx;
    real R0 = a * a / (b * b); // 垂直反射率 0.25 / 2.5 = 0.1 @todo 根拠調査
    real c = 1 - (into ? -ddn : tdir.dot(n));
    real fresnel = R0 + (1 - R0)*c*c*c*c*c;
    // 0.25 と 0.5は結果を綺麗にするためのヒューリスティックな調整値。
    // 屈折する確率も反射する確率も最低限25%にするということ。例えば.1 + (.8 * fresnel)でもよい。
    // この調整が無ければRP, TPは1になって、下でRP, TP掛ける必要はなくなる。
    real P = fresnel;
    if (erand48(xi_) < P)
        PhotonTracing(reflRay, power, depth);          // 反射
    else
        PhotonTracing(Ray(x, tdir), power, depth);  // 屈折

    return;
}


Vec PhotonMapRenderer::Irradiance(const Ray &r, int depth)
{
    // max refl
    if (++depth > 5)
    {
        return Vec();
    }

    HitRecord rec;
    if (!Intersect(r, rec)) return Vec();
    //const Shape& obj = *g_shapes[id];       // the hit object
    Vec x = r.o + r.d * rec.t;
    Vec n = rec.normal;
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;   // 交点の法線
    Vec f = rec.color;
    Refl_t refl = rec.refl;
    
    
    // 0.5にしたらカラーが反射率になってるから暗くなるだけ。IDEALでない反射は扱えない。カラーと混ぜるとかもない。
    // Ideal DIFFUSE reflection
    if (refl == DIFF){
        float pos[3] = { x.x, x.y, x.z };
        float nrm[3] = { nl.x, nl.y, nl.z };
        float irrad[3] = { 0, 0, 0 };
        pPhotonMap_->irradiance_estimate(irrad, pos, nrm, config_.estimateDist, config_.nEstimatePhotons);
        //fprintf(stderr, "irrad %f %f %f\r", irrad[0], irrad[1], irrad[2]);
        return Vec(f.x * irrad[0], f.y * irrad[1], f.z * irrad[2]);
    }
    else if (refl == SPEC) {
        // Ideal SPECULAR reflection
        return Irradiance(Ray(x,r.d-n*2*n.dot(r.d)), depth);
    }
    else if (refl == PHONGMETAL) {
        Vec irrad;
        for (int i = 0; i < 16; i++) {
            // Imperfect SPECULAR reflection
            real r1 = 2*M_PI*erand48(xi_);
            real r2 = erand48(xi_);
            real exponent = 10.0;
            real cosTheta = pow(1-r2, 1.0/(exponent+1.0));
            real sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            real rx = cos(r1) * sinTheta;
            real ry = sin(r1) * sinTheta;
            real rz = cosTheta;
            Vec w = r.d - n * 2 * n.dot(r.d); // reflected ray
            Vec u = ((fabs(w.x) > .1 ? Vec(0,1) : Vec(1)) % w).normalize(); // binormal
            Vec v = w % u; // tangent
            
            // ucosφsinθ + vsinφsinθ + wcosθ
            Vec rd = (u*rx + v*ry + w*rz).normalize();
            Vec mx = x + n*1e-4; // 自己ヒットしないようにちょっと浮かす
            irrad += Irradiance(Ray(mx, rd), depth);
        }
        return irrad / 16;
    }
    
    // Ideal dielectric REFRACTION
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d));
    bool into = n.dot(nl) > 0; // Ray from outside going in?
    real airRefrIdx = 1;
    real grassRefrIdx = 1.5;
    real nnt = into ? airRefrIdx/grassRefrIdx : grassRefrIdx/airRefrIdx;
    real ddn = r.d.dot(nl); // レイと法線のcos
    real cos2t;
    
    // Total internal reflection
    if ((cos2t = 1-nnt * nnt * (1 - ddn * ddn)) < 0)
        return Irradiance(reflRay, depth);
    
    // 屈折方向
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
    
    real a = grassRefrIdx - airRefrIdx;
    real b = grassRefrIdx + airRefrIdx;
    real R0 = a * a / (b * b); // 垂直反射率 0.25 / 2.5 = 0.1 @todo 根拠調査
    real c = 1 - (into ? -ddn : tdir.dot(n));
    real fresnel = R0 + (1 - R0)*c*c*c*c*c;
    real Tr = 1 - fresnel;
    // 反射屈折両方トレース
    return Irradiance(reflRay, depth) * fresnel + Irradiance(Ray(x,tdir), depth) * Tr;
}

Vec* PhotonMapRenderer::RayTracing()
{
    const u32 w = config_.screenWidth;
    const u32 h = config_.screenHeight;
    const u32 nSub = config_.nSubPixelsSqrt;
    const real subPixelFactor = 1.0 / (real)(nSub*nSub);
    
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());      // cam pos, dir
    Vec cx = Vec(w * .5135 / h);
    Vec cy = (cx % cam.d).normalize() * .5135; // .5135は視野角
    Vec* c = new Vec[w*h];

    #pragma omp parallel for schedule(dynamic, 1) private(r)    // OpenMP

    // Loop over image rows
    for (int y=0; y<h; y++) {
        fprintf(stderr, "RayTracing (%d spp) %5.2f%%\n", nSub*nSub, 100.f * y / (h-1));

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
                    real r1 = 2*erand48(xi_), dx = (r1 < 1) ? sqrt(r1)-1 : 1-sqrt(2-r1);
                    real r2 = 2*erand48(xi_), dy = (r2 < 1) ? sqrt(r2)-1 : 1-sqrt(2-r2);
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
                    Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                            cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;

                    // 140は多分投影面までの距離
                    Vec r = Irradiance(Ray(cam.o + d * 140, d.normalize()), 0);

                    // Camera rays are pushed ^^^^^ forward to start in interior
                    // トーンマップとか特にやってない。クランプしてるだけ。
                    c[i] += r * subPixelFactor;
                }
            }
            //fprintf(stderr, "\r%f %f %f", c[i].x, c[i].y, c[i].z);
        }
    }
    
    return c;
}


#include "MeshLoader.h"

Vec* PhotonMapRenderer::Run()
{
    xi_[0] = 0;
	xi_[1] = 0;
	xi_[2] = (unsigned short)config_.nPhotons;
    
    g_nShapes = g_nSpheres + g_nTriangles;
    g_shapes = new Shape*[g_nShapes];
    for (int i = 0; i < g_nSpheres; i++) {
        g_shapes[i] = &g_spheres[i];
    }
    for (int i = 0; i < g_nTriangles; i++) {
        g_shapes[i + g_nSpheres] = &g_triangles[i];
    }
    
    {
        //const char* OBJ_FILE = "/Users/shuichih/Dev/rt/venusm.obj";
        const char* OBJ_FILE = "/Users/shuichih/Dev/rt/cornell_box.obj";
        ObjLoader loader;
        Mesh* pMesh = loader.Load(OBJ_FILE);
        if (pMesh == NULL)
        {
            printf("File Load Error: %s\n", OBJ_FILE);
            return NULL;
        }
        pMesh->scale(-0.10, 0.10, -0.10);
        pMesh->translate(80, 10, 70);
        g_shapes[g_nShapes-1] = pMesh;
    }
    

	LightSource litSrc(Vec(50.0f, 81.0f, 111.6f), 10000);
	//LightSource litSrc(Vec(50.0f, 81.0f, 81.6f), 10000);
    u32 nLitPhotons = config_.nPhotons; // ライトが複数ならnPhotonsをintensityの比率等で割り振る
	for (int iPhoton = 0; iPhoton < config_.nPhotons; iPhoton++) {
        if (iPhoton % (config_.nPhotons / 100) == 0)
            fprintf(stderr, "PhotonTracing %5.2f%%\n", 100. * iPhoton / config_.nPhotons);
        
		Ray ray = litSrc.GenerateRay();
        float power[3] = { litSrc.intensity_, litSrc.intensity_, litSrc.intensity_ };
		PhotonTracing(ray, power, 0);
	}
    pPhotonMap_->scale_photon_power(1.0f / nLitPhotons);
    pPhotonMap_->balance();
    
    Vec* pColorBuf = RayTracing();
    delete [] g_shapes;
    return pColorBuf;
}
