#ifndef SEVERE3D_TUPLE2_H
#define SEVERE3D_TUPLE2_H


/**
 * x, yの２つの要素を持つクラスです。
 */
template <class T>
class Tuple2
{
public:

	/** 値 */
    union {
        struct {
            T x;
			T y;
        };
        T e[2];
    };


	/**
	 * デフォルトコンストラクタです。
	 */
	Tuple2(): x(0), y(0) {};

	/**
	 * x, y, zをとるコンストラクタです。
	 *
	 * @param x xの値
	 * @param y yの値
	 */
	Tuple2(T _x, T _y): x(_x), y(_y)  {};
	
	/**
	 * x, yをメンバに設定します。
	 * 
	 * @param x xの値
	 * @param y yの値
	 */
	void set(T x, T y) {
		this->x = x;
		this->y = y;
	}

    void sub(const Tuple2& t1, const Tuple2& t2) {
        x = t1.x - t2.x;
        y = t1.y - t2.y;
    }

	/**
	 * このタプルにt1を加算します。
	 */
    void add(const Tuple2& t1) {
        x += t1.x;
        y += t1.y;
    }

	/**
	 * このタプルからt1を減じます。
	 */
    void sub(const Tuple2& t1) {
        x -= t1.x;
        y -= t1.y;
    }

	/**
	 * s倍します。
	 */
    void scale(T s) {
        x *= s;
        y *= s;
    }

	/**
	 * 反転します。
	 */
	Tuple2& negate() {
        x = -x;
        y = -y;
		return *this;
    }

	/**
	 * t1 と t2 を alpha:(1-alpha)で線形補間します。
	 */
	void interpolate(const Tuple2<T>& t1, const Tuple2<T>& t2, T alpha)
	{
		T beta = 1.0f - alpha;
		x = t1.x*alpha + t2.x*beta;
		y = t1.y*alpha + t2.y*beta;
	}

    Tuple2& operator+=(const Tuple2& t1) {
        add(t1);
        return *this;
    }
    Tuple2& operator-=(const Tuple2& t1) {
        sub(t1);
        return *this;
    }
    Tuple2& operator*=(T s) {
        scale(s);
        return *this;
    }
    Tuple2& operator/=(T s) {
        scale(1.0f/s);
        return *this;
    }
    Tuple2 operator+(const Tuple2& t1) const {
        return (Tuple2(*this)).operator+=(t1);
    }
    Tuple2 operator-(const Tuple2& t1) const {
        return (Tuple2(*this)).operator-=(t1);
    }
    Tuple2 operator*(T s) const {
        return (Tuple2(*this)).operator*=(s);
    }
    Tuple2 operator/(T s) const {
        return (Tuple2(*this)).operator/=(s);
    }

	// NOTE:コピーコンストラクタおよび = 演算子はコンパイラによって生成されます。

};

typedef Tuple2<float> Tuple2f;
typedef Tuple2<double> Tuple2d;


#endif // SEVERE3D_TUPLE2_H

