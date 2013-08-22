#ifndef _PerlinNoise_H_
#define _PerlinNoise_H_

#include "Common.h"

enum InterporatorType
{
    Interp_Linear,
    Interp_Cosine,
    Interp_Hermite3d,
    Interp_Hermite5d,
    Interp_Cubic
};

// 2次元パーリンノイズ
class PerlinNoise2D
{
public:
    PerlinNoise2D();
    
    float IntNoise(int x, int y);
    float Noise(float x, float y);
    void SetInterporatorType(InterporatorType interpType);
    
private:
    
    //float IntNoise(int x, int y);
    InterporatorType interpType_;
    u8* permutation_;
};


#endif
