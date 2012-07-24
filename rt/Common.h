#ifndef Common_H
#define Common_H

#include <cmath>

//#define USE_FLOAT

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
#else
typedef double              real;
#endif

//----------------------------------------------------------------
struct Vec {
    
    union
    {
        struct
        {
            real x, y, z;
        };
        real e[3];
    };

    inline Vec(real x_=0, real y_=0, real z_=0)
        : x(x_), y(y_), z(z_)
    {
        x=x_;
        y=y_;
        z=z_;
    }
    
    inline Vec operator+(const Vec &b) const
    {
        return Vec(x+b.x,y+b.y,z+b.z);
    }
   
    inline Vec operator-(const Vec &b) const
    {
        return Vec(x-b.x,y-b.y,z-b.z);
    }
    
    inline Vec operator*(real b) const
    {
        return Vec(x*b,y*b,z*b);
    }
    
    inline Vec operator/(real b) const
    {
        real bi = 1.0 / b;
        return Vec(x*bi, y*bi, z*bi);
    }
    
    inline Vec mult(const Vec &b) const
    {
        return Vec(x*b.x,y*b.y,z*b.z);
    }
    
    inline Vec& operator+=(const Vec &b)
    {
        x += b.x;
        y += b.y;
        z += b.z;
        return *this;
    }
    
    inline Vec& operator-=(const Vec& b)
    {
        x -= b.x;
        y -= b.y;
        z -= b.z;
        return *this;
    } 
    
    inline Vec& operator*=(real b)
    {
        x *= b;
        y *= b;
        z *= b;
        return *this;
    }
    
    inline Vec& operator/=(real b)
    {
        real bInv = 1.0 / b;
        x *= bInv;
        y *= bInv;
        z *= bInv;
        return *this;
    }
    
    inline Vec& normalize()
    {
        return *this = *this * (1/sqrt(x*x+y*y+z*z));
    }
    
    inline real dot(const Vec &b) const
    {
        return x*b.x+y*b.y+z*b.z;
    }
    
    // cross:
    Vec operator%(const Vec&b) const
    {
        return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);
    }
};

struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

typedef Vec RGB;

#endif

