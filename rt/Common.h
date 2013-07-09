#ifndef Common_H
#define Common_H

#include <cmath>
#include <cfloat>
#include "vector3.h"
#include "tuple2.h"

#define USE_FLOAT

//----------------------------------------------------------------
#define ARRAY_SZ(a) sizeof(a) / sizeof(a[0])
#ifdef USE_FLOAT
#define Deg2Rad(deg) ((deg) * ((float)M_PI) / 180)
#define Rad2Deg(rad) ((rad) * 180 / ((float)M_PI))
#else
#define Deg2Rad(deg) ((deg) * M_PI / 180)
#define Rad2Deg(rad) ((rad) * 180 / M_PI)
#endif

//----------------------------------------------------------------

typedef char                i8;
typedef short               i16;
typedef int                 i32;
typedef long long           i64;
typedef unsigned char       u8;
typedef unsigned short      u16;
typedef unsigned int        u32;
typedef unsigned long long  u64;
#ifdef USE_FLOAT
typedef float               real;
const float REAL_MAX = FLT_MAX;
const float REAL_MIN = FLT_MIN;
const float EPSILON = 4e-3f;
const float PI = (float)M_PI;
const float PI_INV = 1.0f / (float)M_PI;
#else
typedef double              real;
const double REAL_MAX = DBL_MAX;
const double REAL_MIN = DBL_MIN;
const double EPSILON = 2e-4;
const double PI = M_PI;
const double PI_INV = 1.0 / M_PI;
#endif

//----------------------------------------------------------------
//using namespace severe3d;

typedef Vector3f Vec3;
typedef Tuple2f Vec2;
typedef Vec3 RGB;

template <typename T>
inline void Swap(T& lhs, T& rhs) {
    T tmp = lhs;
    lhs = rhs;
    rhs = tmp;
}

#endif

