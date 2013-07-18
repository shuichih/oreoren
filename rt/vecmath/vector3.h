#ifndef SEVERE3D_VECTOR3_H
#define SEVERE3D_VECTOR3_H

#include <cmath>
#include "tuple3.h"
#include "simd.h"

#define VECMATH_NORMALIZE_ACCURACY 1

/**
 * ３次元ベクトルクラスです。
 */
template <class T>
class Vector3 : public Tuple3<T>
{
public:

	/**
	 * デフォルトコンストラクタです。
	 */
	Vector3(): Tuple3<T>() {}

	/**
	 * コンストラクタです。
	 *
	 * @param x Xの値
	 * @param y Yの値
	 * @param z Zの値
	 */
    Vector3(T x, T y, T z): Tuple3<T>(x, y, z) { }

	Vector3(const Tuple3<T>& t): Tuple3<T>(t) { }

	Vector3& operator=(const Tuple3<T>& t) {
#ifdef VECMATH_USE_SIMD
        this->m = t.m;
		return *this;
#else
		this->x = t.x;
		this->y = t.y;
		this->z = t.z;
		return *this;
#endif
	}

	/**
     * このベクトルの長さの２乗を返します。
     * @return the squared length of this vector
     */
    T lengthSquared() const {
#ifdef VECMATH_USE_SIMD
        // 0x71(0111 0001) means using only xyz, and use lowest float as result
        return _mm_cvtss_f32(_mm_dp_ps(this->m, this->m, 0x71));
#else
        return this->x * this->x + this->y * this->y + this->z * this->z;
#endif
    }

    /**
     * このベクトルの長さを返します。
     */
    T length() const {
#ifdef VECMATH_USE_SIMD
        return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(this->m, this->m, 0x71)));
#else
        // CAUTION: float精度
        return (T)sqrtf((float)lengthSquared());
#endif
    }

	/**
	 * このベクトルを正規化します。
	 */
    Vector3& normalize() {
#ifdef VECMATH_USE_SIMD
        // use div, 23bit (7 digits in decimal) precision, much slow
        __m128 len = _mm_sqrt_ps(_mm_dp_ps(this->m, this->m, 0x77));
        this->m = _mm_div_ps(this->m, len);
        return *this;
#else
        T d = length();
        // zero-div may occur.
        this->x /= d;
        this->y /= d;
        this->z /= d;
        return *this;
#endif
    }

	/**
	 * このベクトルを正規化します。
	 */
    Vector3& normalizeFast() {
#ifdef VECMATH_USE_SIMD
        // use rcp with Newton-Raphson method, 22bit (6 digits in decimal) precision.
        // This method avoids _mm_div_ps and is faster than normalize().
        // (from http://www.kaede-software.com/x86_simd_techni/#)
        
        //static const __declspec(align(16)) float valm05[4] = {-0.5f,-0.5f,-0.5f,-0.5f};
        //static const __declspec(align(16)) float val15[4] = {1.5f,1.5f,1.5f,1.5f};
        static const float valm05[4] __attribute__ ((aligned (16))) = {-0.5f,-0.5f,-0.5f,-0.5f};
        static const float val15[4] __attribute__ ((aligned (16))) = {1.5f,1.5f,1.5f,1.5f};
        static const __m128 Q_m0p5 = *(__m128*)valm05;
        static const __m128 Q_1p5 = *(__m128*)val15;
        __m128 xm0 = _mm_dp_ps(this->m, this->m, 0x77);
        __m128 xm1 = _mm_rsqrt_ps(xm0);
        __m128 xm2 = _mm_mul_ps(xm0, xm1);
        xm2 = _mm_mul_ps(xm2, xm1);
        xm2 = _mm_mul_ps(xm2, xm1);
        xm2 = _mm_mul_ps(xm2, Q_m0p5);
        xm0 = _mm_mul_ps(xm1, Q_1p5);
        xm0 = _mm_add_ps(xm0, xm2);
        this->m = _mm_mul_ps(this->m, xm0);
        return *this;
#else
        return normalize();
#endif
    }

                                 
	/**
	 * このベクトルとv1の内積を求めます。
	 * 内積 = |a||b|cosθ です。２つのベクトルが平行に近いほど、大きな値となります。
	 * 
	 */
	T dot(const Vector3& v1) const {
#ifdef VECMATH_USE_SIMD
        return _mm_cvtss_f32(_mm_dp_ps(this->m, v1.m, 0x71));
#else
		return this->x*v1.x + this->y*v1.y + this->z*v1.z;
#endif
	}

	/**
	 * このv1とv2の外積を求め、このベクトルに設定します。
	 * 外積 = |a||b|sinθ です。
	 */
	void cross(const Vector3& v1, const Vector3& v2) {
		// alias safe.
#ifdef VECMATH_USE_SIMD
        this->m = _mm_sub_ps(
            _mm_mul_ps(
                _mm_shuffle_ps(v1.m, v1.m, _MM_SHUFFLE(3, 0, 2, 1)),
                _mm_shuffle_ps(v2.m, v2.m, _MM_SHUFFLE(3, 1, 0, 2))
            ),
            _mm_mul_ps(
                _mm_shuffle_ps(v1.m, v1.m, _MM_SHUFFLE(3, 1, 0, 2)),
                _mm_shuffle_ps(v2.m, v2.m, _MM_SHUFFLE(3, 0, 2, 1))
            )
        );
#else
		this->x = v1.y * v2.z - v1.z * v2.y;
		this->y = v1.z * v2.x - v1.x * v2.z;
		this->z = v1.x * v2.y - v1.y * v2.x;
#endif
	}

	/**
	 * このベクトルtとv1との外積を求めて返します。
	 * 外積 = |a||b|sinθ です。
	 */
	Vector3 cross(const Vector3& v1) {
#ifdef VECMATH_USE_SIMD
        Vector3 v;
        v.cross(*this, v1);
        return v;
#else
		// alias safe.
		return Vector3(
			this->y * v1.z - this->z * v1.y,
			this->z * v1.x - this->x * v1.z,
			this->x * v1.y - this->y * v1.x
		);
#endif
	}

	/**
	 * このベクトルとvec1の角度を求めます。
	 * 求まる角度の単位はラジアンで、範囲は [0 <= θ <= π] となります。
	 *
     * @param v1  the other vector
     * @return 角度 [0,PI]
     */
    T angle(const Vector3& v1) const {
#ifndef SEVERE3D_HIGHLY_ACCTRATE_ANGLE
		// 内積 = |a||b|cosθ なので、内積 / |a| / |b| でcosθが出る。
		// acos(cosθ)して、θを求めている。
        return (T)acosf(dot(v1) / v1.length()); // v2.length());
#else
        // 前記の方法では０とπ付近でacosの誤差が大きいが、以下の実装では誤差が少ない。
		// tanθ = sinθ / cosθ なので、atan(sinθ / cosθ) すれば、θが求まる。
		// atan2はy, xを引数にとり、tanθ = y/xになるようなθを求める。
		// y、xにはそれぞれ外積(|a||b|sinθ)、内積(|a||b|cosθ)を使用することができる。
		// ただし、３次元では外積はベクトルとなるので、その長さを使用(つまりスカラ化)しなければ
		// ならないらしい。
		// ３次元では、最後にatan2の絶対値をとるようにする。幾何的に考えれば、これで
		// ２つのベクトルで決定される平面で見たときの角度が求まることになる。
        Vector3 c;
        c.cross(*this, v1);
        T sin = c.length();

        return abs(atan2(sin, dot(v1)));
#endif
    }

	// 演算結果をTupleでなくVector3にしたいのでオーバーライド
    Vector3& operator+=(const Tuple3<T>& t1) {
        add(t1);
        return *this;
    }
    Vector3& operator-=(const Tuple3<T>& t1) {
        sub(t1);
        return *this;
    }
    Vector3& operator*=(const Tuple3<T>& t1) {
        this->x *= t1.x;
        this->y *= t1.y;
        this->z *= t1.z;
        return *this;
    }
    Vector3& operator/=(const Tuple3<T>& t1) {
        this->x /= t1.x;
        this->y /= t1.y;
        this->z /= t1.z;
        return *this;
    }
    Vector3& operator*=(T s) {
        scale(s);
        return *this;
    }
    Vector3& operator/=(T s) {
        scale(1.0f/s);
        return *this;
    }

    Vector3 operator+(const Tuple3<T>& t1) const {
        return (Vector3(*this)).operator+=(t1);
    }
    Vector3 operator-(const Tuple3<T>& t1) const {
        return (Vector3(*this)).operator-=(t1);
    }
    Vector3 operator*(T s) const {
        return (Vector3(*this)).operator*=(s);
    }
    Vector3 operator/(T s) const {
        return (Vector3(*this)).operator/=(s);
    }

    Vector3 operator-() const {
        return (Vector3(*this)).operator*=(-1);
    }
    
    Vector3 operator%(const Vector3& v1) const {
        return Vector3(*this).cross(v1);
    }

	// NOTE:コピーコンストラクタおよび = 演算子はコンパイラによって生成されます。

};

typedef Vector3<float> Vector3f;
typedef Vector3<double> Vector3d;


// グローバル演算子
template <typename T>
inline
Vector3<T> operator*(T s, const Vector3<T>& v1) {
    return v1 * s;
}

#endif // SEVERE3D_VECTOR3_H

