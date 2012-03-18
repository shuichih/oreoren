#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  // Remove "-fopenmp" for g++ version <4.2
#include "Common.h"
#include "smallpt_fmt.h"
#include "Scene.h"

Sphere g_spheres2[] = {
    //Scene: radius, position, emission, color, material
    Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),    //Left
    Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),    //Rght
    Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),    //Back
    Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),    //Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),    //Botm
    Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),    //Top
    Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),    //Mirr
    Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),    //Glas
    Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12), Vec(),  DIFF)     //Lite
};

inline real clamp(real x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(real x)
{
    return int(pow(clamp(x),1/2.2)*255+.5);
}

inline bool intersect(const Ray& r, real& t, int& id)
{
    int n = sizeof(g_spheres2) / sizeof(g_spheres2[0]);
    real d;
    real inf = 1e20;
    t = inf;

    for(int i=n; i--;) {
        if((d=g_spheres2[i].intersect(r)) && d < t) {
            t=d;
            id=i;
        }
    }

    return t < inf;
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi);
Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    real t;   // distance to intersection
    int id=0;   // id of intersected object
    if (!intersect(r, t, id)) return Vec(); // if miss, return black
    const Sphere &obj = g_spheres2[id];        // the hit object
    Vec x=r.o+r.d*t; // 交点
    Vec n=(x-obj.p).norm(); // 交点の法線
    Vec nl=n.dot(r.d) < 0 ? n : n*-1; // 交点の法線
    Vec f=obj.c;

    real p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
    // max refl
    if (++depth > 5)
    {
        // pは最大のカラー要素
        // 数値の組み合わせ的に、ライトなら必ず下、ミラーとガラスなら上、それ以外ならカラーの最大要素を
        // 1にスケール(拡大)した色か、obj.e(黒)を返す。スケールしてるのは、そうしないとobj.eに入る
        // 確率と元のカラーの暗さ(1に満たなさ)でダブルで暗くする効果が出てしまうから。
        // コードを短くするためのトリックだと思うが分かりにくい。マテリアルにフラグ設けるかemission+colorを
        // 返すでもいい気がする。ミラーとガラスの場合1が返されてるが、カラーを1にしとくじゃ以下のコードのどこかで
        // 都合が悪いんだろうか。
        if (erand48(Xi) < p) f=f*(1/p);
        else return obj.e;
    }
    // 0.5にしたらカラーが反射率になってるから暗くなるだけ。IDEALでない反射は扱えない。カラーと混ぜるとかもない。
    //R.R.
    if (obj.refl == DIFF){
        // Ideal DIFFUSE reflection
        // cosθ = sqrt(1-r_2)
        //    φ = 2πr_1
        // (r1とr2がRealistic Ray Tracingとは逆)
        real r1 = 2*M_PI*erand48(Xi);
        real r2 = erand48(Xi); // => 1-cos^2θ = 1-sqrt(1-r_2)^2 = r_2
        real r2s = sqrt(r2);   // => sinθ = sqrt(1-cos^2θ) = sqrt(r_2)
        Vec w = nl; // normal
        Vec u = ((fabs(w.x) > .1 ? Vec(0,1) : Vec(1)) % w).norm(); // binormal
        Vec v = w % u; // tangent

        // ucosφsinθ + vsinφsinθ + wcosθ
        Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();

        return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
    }
    else if (obj.refl == SPEC) {
        // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
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
        return obj.e + f.mult(radiance(reflRay, depth, Xi));
    
    // 屈折方向
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    
    real a = grassRefrIdx - airRefrIdx;
    real b = grassRefrIdx + airRefrIdx;
    real R0 = a * a / (b * b); // 垂直反射率 0.25 / 2.5 = 0.1 @todo 根拠調査
    real c = 1 - (into ? -ddn : tdir.dot(n));
    real fresnel = R0 + (1 - R0)*c*c*c*c*c;
    real Tr = 1 - fresnel;
    // 0.25 と 0.5はどっから出てきた？ 結果を綺麗にするためのヒューリスティックな調整？
    // そうぽい。この調整が無ければRP, TPは1になって、下のでRP, TP掛ける必要はなくなる。
    // 屈折する確率も反射する確率も最低限25%にするということ。例えば.1 + (.8 * fresnel)でもよい。
    real P = .25 + (.5 * fresnel);
    real RP = fresnel / P;
    real TP = Tr / (1-P);
    Vec retRadiance;
    if (depth > 2) {
        if (erand48(Xi) < P)
            retRadiance = obj.e + radiance(reflRay, depth, Xi) * RP;      // 反射
        else
            retRadiance = obj.e + radiance(Ray(x,tdir), depth, Xi) * TP; // 屈折
    } else {
        // 浅いときは品質のために反射屈折両方トレース
        retRadiance = obj.e + radiance(reflRay, depth, Xi) * fresnel + radiance(Ray(x,tdir), depth, Xi) * Tr;
    }
    return retRadiance;
    /*
    Ray reflRay(x, r.d-n*2*n.dot(r.d)); // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;
    // Ray from outside going in?
    real nc=1, nt=1.5, nnt=into ? nc/nt : nt/nc, ddn=r.d.dot(nl), cos2t;

    if ((cos2t=1-nnt*nnt*(1-ddn*ddn)) < 0) // Total internal reflection
    return obj.e + f.mult(radiance(reflRay,depth,Xi));

    Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();

    real a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));

    real Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);

    return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ?  // Russian roulette
    radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
    radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
     */
}


u8* rt(int w, int h, int samps){
    
    //int w=1024, h=768;
    //int samps = (argc == 2) ? atoi(argv[1])/4 : 1;   // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());      // cam pos, dir
    Vec cx = Vec(w * .5135 / h);
    Vec cy = (cx % cam.d).norm() * .5135; // .5135は視野角っぽい
    Vec r;
    Vec* c = new Vec[w*h];

    #pragma omp parallel for schedule(dynamic, 1) private(r)    // OpenMP

    // Loop over image rows
    for (int y=0; y<h; y++) {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps*4, 100. * y / (h-1));

        unsigned short Xi[3]={ 0,0,y*y*y }; // 乱数seed, 行ごとの関連が(少なくとも見た目上)なくなるように設定
        
        // Loop cols
        for (unsigned short x=0; x<w; x++) {
            
            int i = (h-y-1) * w + x; // カラーバッファのインデックス
            for (int sy=0; sy<2; sy++) {     // 2x2 subpixel rows
                for (int sx=0; sx<2; sx++) {       // 2x2 subpixel cols
                    r = Vec();
                    
                    for (int s=0; s<samps; s++) {
                        // r1, r2 = 0 to 2
                        // dx, dy = -1 to 1  中心に集まったサンプリング --> tent filter
                        real r1 = 2*erand48(Xi), dx = (r1 < 1) ? sqrt(r1)-1 : 1-sqrt(2-r1);
                        real r2 = 2*erand48(Xi), dy = (r2 < 1) ? sqrt(r2)-1 : 1-sqrt(2-r2);
                        // (sx+.5 + dx)/2 --> .5でサブピクセルの中心に。dxでフィルタの揺らぎ。
                        // sx+.5 = 0.5 or 1.5
                        // sx+.5 + dx = -0.5 to 1.5 or 0.5 to 2.5
                        // (sx+.5 + dx)/2 = -0.25 to 0.75 or 0.25 to 1.25、前者0.25中心に+-0.5範囲の揺らぎ、後者0.75中心に+-0.5範囲のゆらぎ
                        // +x / w ピクセルの位置へ、0to1へ。 -.5 0to1から-0.5to0.5へ。
                        // cx* 投影面中心座標から+-0.5の範囲を走査するため
                        Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                                cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;

                        // 140は多分投影面までの距離
                        r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
                    }
                    // Camera rays are pushed ^^^^^ forward to start in interior
                    // トーンマップとか特にやってない。クランプしてるだけ。
                    c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
            }
            //fprintf(stderr, "\r%f %f %f", c[i].x, c[i].y, c[i].z);
        }
    }
    
    static u8* pColorBuf = new u8[w * h * 4];
    for (int i = 0, j=0; i < (w*h); ++i, j+=4)
    {
        //pColorBuf[j+0] = (u8)(toInt(c[i].x));
#if 0
        pColorBuf[j+0] = 255;//(i%4 == 0) ? 0 : 255;
        pColorBuf[j+1] = 0;//(u8)(toInt(c[i].y));
        pColorBuf[j+2] = 0;//(u8)(toInt(c[i].z));
        pColorBuf[j+3] = 255;
#else
        pColorBuf[j+0] = (u8)(toInt(c[i].x));
        pColorBuf[j+1] = (u8)(toInt(c[i].y));
        pColorBuf[j+2] = (u8)(toInt(c[i].z));
        pColorBuf[j+3] = 255;
#endif
        //fprintf(stderr, "\r%f %f %f", c[i].x, c[i].y, c[i].z);
        //fprintf(stderr, "\r%d %d %d", pColorBuf[j+0], pColorBuf[j+1], pColorBuf[j+2]);
    }
    delete [] c;
    
    return pColorBuf;

#if 0
    FILE *f = fopen("/Users/shuichih/image.ppm", "w");  // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);

    for (int i=0; i<w*h; i++)
        fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); //toIntでガンマ補正してる...

    fclose(f);
#endif
}
