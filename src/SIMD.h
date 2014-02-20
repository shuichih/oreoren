#ifndef RT_SIMD_H
#define RT_SIMD_H

//

#define RAY_USE_SIMD
#define VECMATH_DATA_FOR_SIMD
//#define VECMATH_USE_SIMD
//#define MESHTRIANGLE_USE_SIMD
//#define BBOX_USE_SIMD
#define QBVH_USE_SIMD

//

#include "xmmintrin.h" // SSE
#include "emmintrin.h" // SSE2
#include "pmmintrin.h" // SSE3
#include "smmintrin.h" // SSE4.1
//#include "immintrin.h" // for _mm_pow_ps, etc. not exist on xcode
//#include "ia32intrin.h" // SVML but not found in xcode 4.5.2

struct Um128
{
    union
    {
        struct
        {
            float x, y, z, w;
        };
        float e[4];
        __m128 m;
    };
    Um128(const __m128& rhs)
    {
        m = rhs;
    }
    Um128& operator=(const __m128& rhs)
    {
        m = rhs;
        return *this;
    }
};

#endif // RT_SIMD_H
