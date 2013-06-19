#ifndef SEVERE3D_QUAT4_H
#define SEVERE3D_QUAT4_H

// 回転軸と回転角度を持つデータ構造
template<class T>
struct AxisAngle4
{
	AxisAngle4(T _x, T _y, T _z, T _angle):
		x(_x), y(_y), z(_z), angle(_angle) {}
	T x;
	T y;
	T z;
	T angle;
};

typedef AxisAngle4<float> AxisAngle4f;
typedef AxisAngle4<double> AxisAngle4d;

/**
 * クォータニオンクラスです。
 */
template<class T>
class Quat4
{
public:
	T w;
	T x;
	T y;
	T z;

public:
    Quat4(): w(0), x(0), y(0), z(0) {}

	Quat4(const AxisAngle4<T>& a)
	{
		set(a);
	}

	Quat4(T _w, T _x, T _y, T _z): w(_w), x(_x), y(_y), z(_z) {}

	~Quat4() {}


	void set(T _w, T _x, T _y, T _z)
	{
		w = _w;
		x = _x;
		y = _y;
		z = _z;
	}

	// 回転軸と角度から、回転クォータニオンを計算してセットします。
	void set(const AxisAngle4<T>& a)
	{
		scaleSet(a, 1);
	}

	// 回転軸と角度から、回転クォータニオンを計算してセットします。
	void scaleSet(const AxisAngle4<T>& a, T scale)
	{
		// 回転クォータニオンは、以下の式で表される
		// cos(θ/2) + R*sin(θ/2); ここで、Rは回転軸ベクトル

		x = a.x;
		y = a.y;
		z = a.z;
		T n = sqrt(x*x + y*y + z*z);
//		T n = 1;
		T s = sin(0.5f * a.angle * scale) / n;	// 正規化のためにnで割る
		x *= s;
		y *= s;
		z *= s;
		w = cos(0.5f * a.angle * scale);
	}

	// 
	void interpolate(const AxisAngle4<T>& a1, const AxisAngle4<T>& a2, T alpha)
	{
		scaleSet(a1, alpha);
		Quat4<T> q;
		q.scaleSet(a2, 1.0f - alpha);
		mul(*this, q);
	}

	void mul(const Quat4<T>& q1, const Quat4<T>& q2)
	{
		// store stack for aliasing-safty
		set(
			q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z,
			q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y,
			q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x,
			q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w
		);
	}

	void mul(const Quat4<T>& q1)
	{
		mul(*this, q1);
	}
	
	void mul(T n)
	{
		mul(*this, n);
	}

	void mul(const Quat4<T>& q, T n)
	{
		w = q.w * n;
		x = q.x * n;
		y = q.y * n;
		z = q.z * n;
	}

    T norm() const
	{
        return x*x + y*y + z*z + w*w;
    }

    void normalize()
	{
        T in = 1 / sqrt(norm());
        // zero-div may occur.
        w = w * in;
        x = x * in;
        y = y * in;
        z = z * in;
    }

	void identity()
	{
		w = 1;
		x = 0;
		y = 0;
		z = 0;
	}



};

typedef Quat4<float>  Quat4f;
typedef Quat4<double> Quat4d;


#endif // SEVERE3D_QUAT4_H
