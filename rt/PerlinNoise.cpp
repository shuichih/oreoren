#include "PerlinNoise.h"
#include "Common.h"
#include "Random.h"

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
        /*
        int A = permutation_[ix  ] + y;
        int B = permutation_[ix+1] + y;
        int AA = permutation_[A];
        int BA = permutation_[B];
        */
        float v0 = IntNoise(ix,   iy);
        float v1 = IntNoise(ix+1, iy);
        float v2 = IntNoise(ix,   iy+1);
        float v3 = IntNoise(ix+1, iy+1);
        float (*fade)(float) = pFadeFuncs[interpType_];
        float f = fade(fx);
        v0 = Lerp(v0, v1, f);
        v1 = Lerp(v2, v3, f);
        ret = Lerp(v0, v1, fade(fy));
    }
    return ret;
}

