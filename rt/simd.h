#ifndef RT_SIMD_H
#define RT_SIMD_H

//

#define RAY_USE_SIMD
//#define VECMATH_DATA_FOR_SIMD
//#define VECMATH_USE_SIMD
//#define MESHTRIANGLE_USE_SIMD
//#define BBOX_USE_SIMD
#define QBVH_USE_SIMD

//

//#include "xmmintrin.h" // SSE
//#include "emmintrin.h" // SSE2
#include "pmmintrin.h" // SSE3
#include "smmintrin.h" // SSE4.1
//#include "immintrin.h"
//#include "ia32intrin.h" // SVML but not found


#endif // RT_SIMD_H
