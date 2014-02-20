#ifndef ONB_h
#define ONB_h

#include "Common.h"

//--------------------------------------------------------------------------------
/**
 * Orthonormal Basis
 */
class ONB
{
public:
    ONB();
    ONB(const Vec3& u, const Vec3& v, const Vec3& w);
    
    static ONB InitFromW(const Vec3& w);
    
    /**
     * 基準座標系からONBの座標系に変換する
     */
    inline Vec3 Transform(const Vec3& vec)
    {
        return vec.x * u + vec.y * v + vec.z * w;
    }
    
    Vec3 u;
    Vec3 v;
    Vec3 w;
};

#endif
