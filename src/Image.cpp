#include "Image.h"
#include "vector3.h"
#include <cassert>
#include <cstring>

using namespace std;



namespace {
    inline real clamp(real x)
    {
        return x < 0.f ? 0.f : x > 1.f ? 1.f : x;
    }
}


Image::Image(int width, int height, Image::Format format)
: width_(width)
, height_(height)
, format_(format)
{
    if (format == RGBA_8) {
        pBuf_ = new u8[width * height * 4];
        memset(pBuf_, 0, sizeof(width * height * 4));
    } else if (format == RGB_F32) {
        pBuf_ = (u8*)(new Vec3[width * height]);
    }
}

Image::~Image()
{
    if (pBuf_) {
        delete pBuf_;
        pBuf_ = NULL;
    }
}

void Image::Clear()
{
    int sz = 0;
    if (format_ == RGBA_8) sz = sizeof(u8) * 4;
    else if (format_ == RGB_F32) sz = sizeof(Vec3);
    else {
        assert(false);
        return;
    }
    
    memset(pBuf_, 0, width_ * height_ * sz);
}

bool Image::ConvertFrom(const Image& image)
{
    if (width_ != image.width()) return false;
    if (height_ != image.height()) return false;
    
    if (format_ == RGBA_8 && image.format() == RGB_F32) {
        RGBA8* pDest = (RGBA8*)pBuf_;
        const Vec3* pSrc = (const Vec3*)image.buffer();
        for (int i=0; i<(width_*height_); ++i)
        {
            pDest[i].r = (u8)(clamp(pSrc[i].x) * 255.999f);
            pDest[i].g = (u8)(clamp(pSrc[i].y) * 255.999f);
            pDest[i].b = (u8)(clamp(pSrc[i].z) * 255.999f);
            pDest[i].a = 255;
        }
    } else {
        assert(false);
        return false;
    }
    return true;
}

bool Image::GammaCorrection()
{
    float invGamma = 1/2.2f;
    if (format_ == RGB_F32) {
        Vec3* pDest = (Vec3*)pBuf_;
        for (int i=0, j=0; i<(width_*height_); ++i, j+=4) {
            pDest[i].x = powf(pDest[i].x, invGamma);
            pDest[i].y = powf(pDest[i].y, invGamma);
            pDest[i].z = powf(pDest[i].z, invGamma);
        }
    }
    else {
        return false;
    }
    return true;
}


