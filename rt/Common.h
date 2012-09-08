#ifndef Common_H
#define Common_H

#include <cmath>
#include <cfloat>

#define USE_FLOAT

//----------------------------------------------------------------
#define ARRAY_SZ(a) sizeof(a) / sizeof(a[0])

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
#else
typedef double              real;
const double REAL_MAX = DBL_MAX;
const double REAL_MIN = DBL_MIN;
#endif

//----------------------------------------------------------------
struct Vec3 {
    
    union
    {
        struct
        {
            real x, y, z;
        };
        real e[3];
    };

    inline Vec3(real x_=0, real y_=0, real z_=0)
        : x(x_), y(y_), z(z_)
    {
        x=x_;
        y=y_;
        z=z_;
    }
    
    inline Vec3 operator+(const Vec3 &b) const
    {
        return Vec3(x+b.x,y+b.y,z+b.z);
    }
   
    inline Vec3 operator-(const Vec3 &b) const
    {
        return Vec3(x-b.x,y-b.y,z-b.z);
    }
    
    inline Vec3 operator*(real b) const
    {
        return Vec3(x*b,y*b,z*b);
    }
    
    inline Vec3 operator/(real b) const
    {
        real bi = 1.0f / b;
        return Vec3(x*bi, y*bi, z*bi);
    }
    
    inline Vec3 mult(const Vec3 &b) const
    {
        return Vec3(x*b.x,y*b.y,z*b.z);
    }
    
    inline real sum() const
    {
        return x + y + z;
    }
    
    inline Vec3& operator+=(const Vec3 &b)
    {
        x += b.x;
        y += b.y;
        z += b.z;
        return *this;
    }
    
    inline Vec3& operator-=(const Vec3& b)
    {
        x -= b.x;
        y -= b.y;
        z -= b.z;
        return *this;
    } 
    
    inline Vec3& operator*=(real b)
    {
        x *= b;
        y *= b;
        z *= b;
        return *this;
    }
    
    inline Vec3& operator/=(real b)
    {
        real bInv = 1.0f / b;
        x *= bInv;
        y *= bInv;
        z *= bInv;
        return *this;
    }
    
    inline Vec3& normalize()
    {
        return *this = *this * (1.f/sqrtf(x*x+y*y+z*z));
    }
    
    inline real dot(const Vec3 &b) const
    {
        return x*b.x+y*b.y+z*b.z;
    }
    
    // cross:
    Vec3 operator%(const Vec3&b) const
    {
        return Vec3(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);
    }
};

struct Ray
{
    Vec3 o, d;
    Ray(Vec3 o_, Vec3 d_) : o(o_), d(d_) {}
};

typedef Vec3 RGB;

#endif

