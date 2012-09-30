#ifndef Common_H
#define Common_H

#include <cmath>
#include <cfloat>

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
#else
typedef double              real;
const double REAL_MAX = DBL_MAX;
const double REAL_MIN = DBL_MIN;
const double EPSILON = 2e-4;
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

    inline Vec3()
        : x(0), y(0), z(0)
    {}
    
    inline Vec3(real x_, real y_, real z_)
        : x(x_), y(y_), z(z_)
    {
        x=x_;
        y=y_;
        z=z_;
    }
    
    inline Vec3 operator+(const Vec3 &b) const
    {
        return Vec3(x+b.x, y+b.y, z+b.z);
    }
   
    inline Vec3 operator-(const Vec3 &b) const
    {
        return Vec3(x-b.x,y-b.y,z-b.z);
    }
    
    inline Vec3 operator+(const real b) const
    {
        return Vec3(x+b, y+b, z+b);
    }
    
    inline Vec3 operator-(const real b) const
    {
        return Vec3(x-b, y-b, z-b);
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
    
    inline real square_length()
    {
        return x * x + y * y + z * z;
    }
    
    inline real length()
    {
        return sqrtf(x * x + y * y + z * z);
    }
    
    inline Vec3& normalize()
    {
        return *this = *this / length();
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

Vec3 operator*(const float f, const Vec3& v);

struct Ray
{
    Vec3 o;
    Vec3 d;
    Vec3 invDir;
    i32 dirSign[3]; // dirの符号 0:positive 1:negative
    
    Ray(Vec3 o_, Vec3 d_) : o(o_), d(d_)
    {
        SetDirection(d_);
    }
    
    void SetDirection(const Vec3& dir)
    {
        d = dir;
        invDir = Vec3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
        dirSign[0] = (dir.x > 0) ? 0 : 1;
        dirSign[1] = (dir.y > 0) ? 0 : 1;
        dirSign[2] = (dir.z > 0) ? 0 : 1;
    }

    Vec3 PointAtParameter(float t) const
    {
        return o + t * d;
    }
};


typedef Vec3 RGB;

#endif

