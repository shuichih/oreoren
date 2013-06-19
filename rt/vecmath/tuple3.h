#ifndef SEVERE3D_TUPLE3_H
#define SEVERE3D_TUPLE3_H

#include "../simd.h"

/**
 * x, y, z の３つの要素を持つクラスです。
 */
template <class T>
class Tuple3
{
public:

    union
    {
        struct
        {
            T x, y, z;
        };
        T e[3];
#ifdef VECMATH_DATA_FOR_SIMD
        __m128 m;
#endif
    };

	/**
	 * デフォルトコンストラクタです。
	 */
#ifdef VECMATH_USE_SIMD
	Tuple3()
    {
        m = _mm_setzero_ps();
    }
#else
	Tuple3(): x(0), y(0), z(0) {}
#endif

	/**
	 * x, y, zをとるコンストラクタです。
	 *
	 * @param x xの値
	 * @param y yの値
	 * @param x zの値
	 */
#ifdef VECMATH_USE_SIMD
	Tuple3(T _x, T _y, T _z)
    {
        //m = _mm_set_ps(_x, _y, _z, 0);
        m = _mm_set_ps(0, _z, _y, _x);
    }
#else
	Tuple3(T _x, T _y, T _z): x(_x), y(_y), z(_z) {}
#endif
    
#ifdef VECMATH_USE_SIMD
    /**
     * __m128をとるコンストラクタです。
     */
    Tuple3(__m128 _m) : m(_m) {}
#endif
	
	/**
	 * x, y, zをメンバに設定します。
	 * 
	 * @param x xの値
	 * @param y yの値
	 * @param x zの値
	 */
	void set(T x, T y, T z) {
#ifdef VECMATH_USE_SIMD
        m = _mm_set_ps(0, z, y, x);
#else
		this->x = x;
		this->y = y;
		this->z = z;
#endif
	}

    void sub(const Tuple3& t1, const Tuple3& t2) {
#ifdef VECMATH_USE_SIMD
        this->m = _mm_sub_ps(t1.m, t2.m);
#else
        x = t1.x - t2.x;
        y = t1.y - t2.y;
        z = t1.z - t2.z;
#endif
    }

	/**
	 * このタプルにv1を加算します。
	 */
    void add(const Tuple3& t1) {
#ifdef VECMATH_USE_SIMD
        this->m = _mm_add_ps(this->m, t1.m);
#else
        x += t1.x;
        y += t1.y;
        z += t1.z;
#endif
    }

	/**
	 * このタプルからt1を減じます。
	 */
    void sub(const Tuple3& t1) {
#ifdef VECMATH_USE_SIMD
        this->m = _mm_sub_ps(this->m, t1.m);
#else
        x -= t1.x;
        y -= t1.y;
        z -= t1.z;
#endif
    }

	/**
	 * s倍します。
	 */
    void scale(T s) {
#ifdef VECMATH_USE_SIMD
        __m128 sm = _mm_set1_ps(s);
        this->m = _mm_mul_ps(this->m, sm);
#else
        x *= s;
        y *= s;
        z *= s;
#endif
    }

	/**
	 * 反転します。
	 */
	Tuple3& negate() {
#ifdef VECMATH_USE_SIMD
        __m128 sm = _mm_set1_ps(-1.f);
        this->m = _mm_mul_ps(this->m, sm);
#else
        x = -x;
        y = -y;
        z = -z;
#endif
		return *this;
    }
    
    Tuple3 mult(const Tuple3 &b) const
    {
#ifdef VECMATH_USE_SIMD
        return _mm_mul_ps(this->m, b.m);
#else
        return Tuple3(x*b.x,y*b.y,z*b.z);
#endif
    }
    
    T sum() const
    {
        return x + y + z;
    }
    
	/**
	 * t1 と t2 を alpha:(1-alpha)で線形補間します。
	 */
	void interpolate(const Tuple3<T>& t1, const Tuple3<T>& t2, T alpha)
	{
#ifdef VECMATH_USE_SIMD___
        __m128 xm0 = _mm_set1_ps(alpha);
        __m128 xm1 = _mm_sub_ps(t2.m, t1.m);
        xm0 = _mm_mul_ps(xm0, xm1);
        m = _mm_add_ps(t1.m, xm0);
#else
        *this = t1 + alpha * (t2 - t1);
#endif
	}

    Tuple3& operator+=(const Tuple3& t1)
    {
        add(t1);
        return *this;
    }
    Tuple3& operator-=(const Tuple3& t1)
    {
        sub(t1);
        return *this;
    }
    Tuple3& operator*=(T s)
    {
        scale(s);
        return *this;
    }
    Tuple3& operator/=(T s)
    {
        scale(1.0f/s);
        return *this;
    }
    Tuple3 operator+(const Tuple3& t1) const
    {
        return (Tuple3(*this)).operator+=(t1);
    }
    Tuple3 operator-(const Tuple3& t1) const
    {
        return (Tuple3(*this)).operator-=(t1);
    }
    Tuple3 operator*(T s) const
    {
        return (Tuple3(*this)).operator*=(s);
    }
    Tuple3 operator/(T s) const
    {
        return (Tuple3(*this)).operator/=(s);
    }

	// NOTE:コピーコンストラクタおよび = 演算子はコンパイラによって生成されます。

};

typedef Tuple3<float> Tuple3f;
typedef Tuple3<double> Tuple3d; // CAUTION: isn't supported when VECMATH_USE_SIMD

// グローバル演算子
template <class T>
inline Tuple3<T> operator*(T s, const Tuple3<T>& t1)
{
	return t1 * s;
}


#endif // SEVERE3D_TUPLE3_H

