#include "PerlinNoise.h"
#include "Common.h"
#include "Random.h"
#include "stdio.h"

static int p[512];
static const int permutation[] = { 151,160,137,91,90,15,
    131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
    190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
    88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
    77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
    102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
    135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
    5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
    223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
    129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
    251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
    49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};

static double fade(double t) { return t * t * t * (t * (t * 6 - 15) + 10); }
static double lerp(double t, double a, double b) { return a + t * (b - a); }
static double grad(int hash, double x, double y, double z) {
    int h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE
    double u = h<8 ? x : y,                 // INTO 12 GRADIENT DIRECTIONS.
    v = h<4 ? y : h==12||h==14 ? x : z;
    return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
}

static double noise(double x, double y, double z) {
    int X = int(x) & 255;                              // FIND UNIT CUBE THAT
    int Y = int(y) & 255;                              // CONTAINS POINT.
    int Z = int(z) & 255;
    x -= int(x);                                       // FIND RELATIVE X,Y,Z
    y -= int(y);                                       // OF POINT IN CUBE.
    z -= int(z);
    double u = fade(x);                                // COMPUTE FADE CURVES
    double v = fade(y);                                // FOR EACH OF X,Y,Z.
    double w = fade(z);
    int A = p[X  ]+Y, AA = p[A]+Z, AB = p[A+1]+Z;      // HASH COORDINATES OF
    int B = p[X+1]+Y, BA = p[B]+Z, BB = p[B+1]+Z;      // THE 8 CUBE CORNERS,
    //double d = grad(p[AA  ], x  , y  , z);
    //printf("%f %f %f  %f\n", x, y, z, d);
    return lerp(w, lerp(v, lerp(u, grad(p[AA  ], x  , y  , z   ),  // AND ADD
                                   grad(p[BA  ], x-1, y  , z   )), // BLENDED
                           lerp(u, grad(p[AB  ], x  , y-1, z   ),  // RESULTS
                                   grad(p[BB  ], x-1, y-1, z   ))),// FROM  8
                   lerp(v, lerp(u, grad(p[AA+1], x  , y  , z-1 ),  // CORNERS
                                   grad(p[BA+1], x-1, y  , z-1 )), // OF CUBE
                           lerp(u, grad(p[AB+1], x  , y-1, z-1 ),
                                   grad(p[BB+1], x-1, y-1, z-1 ))));
}
// 線形フェード
static float LinearFade(float t)
{
    return t;
}

// cosフェード
static float CosineFade(float t)
{
    return (1.f - cosf(t * PI)) * .5f;
}

// 3次エルミートフェード
static float Hermite3dFade(float t)
{
    return t * t * (3 - 2 * t);
}

// 5次エルミートフェード
static float Hermite5dFade(float t)
{
    return t * t * t * (t * (t * 6 - 15) + 10);
}

// 線形補間
static float Lerp(float a, float b, float t)
{
    return (1-t) * a + b * t;
}

// Cubic補間
static float CubicInterporate(float v0, float v1, float v2, float v3, float t)
{
    float p = (v3 - v2) - (v0 - v1);
    float q = (v0 - v1) - p;
    float r = v2 - v0;
    return t * (t * (t * p + q) + r) + v1;
}

// 補間関数ポインタ配列
static float (*pFadeFuncs[])(float t) =
{
    LinearFade,
    CosineFade,
    Hermite3dFade,
    Hermite5dFade
};

// ctor
PerlinNoise2D::PerlinNoise2D()
: interpType_(Interp_Cubic)
{
    Random rand;
    permutation_ = new u8[512];
    for (int i=0; i<256; i++) {
        float to0_255 = 0.9999999f / (0xFFFFFFFF / 0xFF);
        permutation_[i] = permutation_[256+i] = u8(rand.Generate() * to0_255);
    }
    
    for (int i=0; i < 256 ; i++) {
        p[256+i] = p[i] = permutation[i];
    }
}

// 補間方法設定
void PerlinNoise2D::SetInterporatorType(InterporatorType interpType)
{
    interpType_ = interpType;
}

// 2つのintから -1.0 to 1.0 の値を返す
float PerlinNoise2D::IntNoise(int x, int y)
{
    // 57, 15731, 789221, 1376312589は適当な素数で変えてもいい
    // 最後の除算は最終結果を-1.0〜1.0に収めるため
    const double div = 1.0 / 1073741824.0;
    int n = x + y * 57;
    n = (n << 13) ^ n;
    return (float)(1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) * div);
}

// パーリンノイズ算出
float PerlinNoise2D::Noise(float x, float y)
{
#if 0
    int ix = int(x) & 255;
    int iy = int(y) & 255;
    float fx = x - int(x);
    float fy = y - int(y);
    float ret;
    
    if (interpType_ == Interp_Cubic)
    {
        float w[4];
        for (int i = 0; i < 4; i++) {
            float v0 = IntNoise(ix-1, iy-1 + i);
            float v1 = IntNoise(ix,   iy-1 + i);
            float v2 = IntNoise(ix+1, iy-1 + i);
            float v3 = IntNoise(ix+2, iy-1 + i);
            w[i] = CubicInterporate(v0, v1, v2, v3, fx);
        }
        ret = CubicInterporate(w[0], w[1], w[2], w[3], fy);
    }
    else
    {
        printf("%d %d %f %f\n", ix, iy, fx, fy);
        const float to0_1 = 1.f / 255.f;
        int A = permutation_[ix  ] + iy;
        int B = permutation_[ix+1] + iy;
        float v0 = permutation_[A] * to0_1;
        float v1 = permutation_[B] * to0_1;
        float v2 = permutation_[A+1] * to0_1;
        float v3 = permutation_[B+1] * to0_1;
        float (*fade)(float) = pFadeFuncs[interpType_];
        float f = fade(fx);
        v0 = Lerp(v0, v1, f);
        v1 = Lerp(v2, v3, f);
        ret = Lerp(v0, v1, fade(fy));
    }
    return ret;
#else
    float ret = (float)noise(x, y, 0);
    //printf("%f\n", ret);
    return ret;
#endif
}
