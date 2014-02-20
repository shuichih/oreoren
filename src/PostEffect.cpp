#include "PostEffect.h"

ToneMap::ToneMap()
: keyValue_(0.18f)
, delta_(0.0001f)
, smallestWhiteLum_(-1.f)
{
}


ToneMap::~ToneMap()
{
}

void ToneMap::SetKeyValue(real keyValue)
{
    keyValue_ = keyValue;
}

void ToneMap::SetDelta(real delta)
{
    delta_ = delta;
}

void ToneMap::SetSmallestWhiteLuminance(real l)
{
    smallestWhiteLum_ = l;
}

void ToneMap::Apply(Vec3* pBuffer, i32 bufferWidth, i32 bufferHeight)
{
    const Vec3 RGB2Y  (+0.29900f, +0.58700f, +0.11400f);
    const Vec3 RGB2Cb (-0.16874f, -0.33126f, +0.50000f);
    const Vec3 RGB2Cr (+0.50000f, -0.41869f, -0.08131f);
    const Vec3 YCbCr2R(+1.00000f, +0.00000f, +1.40200f);
    const Vec3 YCbCr2G(+1.00000f, -0.34414f, -0.71414f);
    const Vec3 YCbCr2B(+1.00000f, +1.77200f, +0.00000f);
    
    real lmSum = 0.0f;
    real delta = delta_;
    i32 num = bufferWidth * bufferHeight;
    real lmMax = 0;
    for (i32 i = 0; i < num; i++)
    {
        real lm = pBuffer[i].dot(RGB2Y);
        if (lm <= 0) {
            // @todo debug 負の値入ってる
            //int a = 1;
        }
        lmSum += logf(delta + lm);
        if (lm > lmMax) {
            lmMax = lm;
        }
    }
    
    real lwm = expf(lmSum / num);
    real scale = keyValue_ / lwm;
    real lmWhite = (smallestWhiteLum_ == -1) ? scale*lmMax : smallestWhiteLum_;
    real wL2Inv = 1.0f / (lmWhite * lmWhite);
    //real wL2Inv = 0;
    for (i32 i = 0; i < num; i++)
    {
        Vec3& pixel = pBuffer[i];
        Vec3 YCbCr;
        
        // 色成分はそのまま
        YCbCr.y = RGB2Cb.dot(pixel);
        YCbCr.z = RGB2Cr.dot(pixel);
        
        // 輝度成分を調整
        real lm = RGB2Y.dot(pixel);
        real L = scale * lm;
        YCbCr.x = (L * (1 + L * wL2Inv) / (1 + L));
        
        // RGBにする
        pixel.x = YCbCr2R.dot(YCbCr);
        pixel.y = YCbCr2G.dot(YCbCr);
        pixel.z = YCbCr2B.dot(YCbCr);
    }
}

