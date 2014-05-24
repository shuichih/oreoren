#ifndef Rayzer_Image_H
#define Rayzer_Image_H

#include "Common.h"

/**
 * 画像クラス
 */
class Image
{
public:
    enum Format
    {
        RGB_F32,
        RGBA_8
    };
    
    Image(int width, int height, Format format);
    ~Image();
    
    void Clear();
    bool ConvertFrom(const Image& image);
    bool GammaCorrection();
    
    int width() const { return width_; }
    int height() const { return height_; }
    Format format() const { return format_; }
    u8* buffer() { return pBuf_; }
    const u8* buffer() const { return pBuf_; }
    
private:
    struct RGBA8 {
        u8 r;
        u8 g;
        u8 b;
        u8 a;
    };

    int width_;
    int height_;
    Format format_;
    u8* pBuf_;
};

#endif // Rayzer_Image_H
