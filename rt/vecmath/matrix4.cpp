#include "tuple3.h"
#include "matrix4.h"



// デフォルトコンストラクタ
template<class T>
Matrix4<T>::Matrix4():
    m00(0), m01(0), m02(0), m03(0),
    m10(0), m11(0), m12(0), m13(0),
    m20(0), m21(0), m22(0), m23(0),
    m30(0), m31(0), m32(0), m33(0) { }

// コンストラクタ
template<class T>
Matrix4<T>::Matrix4(T n00, T n01, T n02, T n03, 
                    T n10, T n11, T n12, T n13,
                    T n20, T n21, T n22, T n23,
                    T n30, T n31, T n32, T n33):
    m00(n00), m01(n01), m02(n02), m03(n03),
    m10(n10), m11(n11), m12(n12), m13(n13),
    m20(n20), m21(n21), m22(n22), m23(n23),
    m30(n30), m31(n31), m32(n32), m33(n33) { }

// 全要素設定
template<class T>
void Matrix4<T>::set(T n00, T n01, T n02, T n03, 
                     T n10, T n11, T n12, T n13,
                     T n20, T n21, T n22, T n23,
                     T n30, T n31, T n32, T n33)
{
    m00 = n00; m01 = n01; m02 = n02; m03 = n03;
    m10 = n10; m11 = n11; m12 = n12; m13 = n13;
    m20 = n20; m21 = n21; m22 = n22; m23 = n23;
    m30 = n30; m31 = n31; m32 = n32; m33 = n33;
}

// 単位行列
template<class T>
Matrix4<T>& Matrix4<T>::identity()
{
	m00 = 1; m01 = 0; m02 = 0; m03 = 0;
	m10 = 0; m11 = 1; m12 = 0; m13 = 0;
	m20 = 0; m21 = 0; m22 = 1; m23 = 0;
	m30 = 0; m31 = 0; m32 = 0; m33 = 1;

	return *this;
}

// 行列の定数倍
template<class T>
void Matrix4<T>::mul(T scalar)
{
	m00 *= scalar; m01 *= scalar;  m02 *= scalar; m03 *= scalar;
	m10 *= scalar; m11 *= scalar;  m12 *= scalar; m13 *= scalar;
	m20 *= scalar; m21 *= scalar;  m22 *= scalar; m23 *= scalar;
	m30 *= scalar; m31 *= scalar;  m32 *= scalar; m33 *= scalar;
}

// 行列乗算
template<class T>
Matrix4<T>& Matrix4<T>::mul(const Matrix4<T>& m1, const Matrix4<T>& m2)
{
	// alias safe way
	Matrix4		m3;

	m3.m00 = m1.m00*m2.m00 + m1.m01*m2.m10 + m1.m02*m2.m20 + m1.m03*m2.m30;
	m3.m01 = m1.m00*m2.m01 + m1.m01*m2.m11 + m1.m02*m2.m21 + m1.m03*m2.m31;
	m3.m02 = m1.m00*m2.m02 + m1.m01*m2.m12 + m1.m02*m2.m22 + m1.m03*m2.m32;
	m3.m03 = m1.m00*m2.m03 + m1.m01*m2.m13 + m1.m02*m2.m23 + m1.m03*m2.m33;

	m3.m10 = m1.m10*m2.m00 + m1.m11*m2.m10 + m1.m12*m2.m20 + m1.m13*m2.m30;
	m3.m11 = m1.m10*m2.m01 + m1.m11*m2.m11 + m1.m12*m2.m21 + m1.m13*m2.m31;
	m3.m12 = m1.m10*m2.m02 + m1.m11*m2.m12 + m1.m12*m2.m22 + m1.m13*m2.m32;
	m3.m13 = m1.m10*m2.m03 + m1.m11*m2.m13 + m1.m12*m2.m23 + m1.m13*m2.m33;

	m3.m20 = m1.m20*m2.m00 + m1.m21*m2.m10 + m1.m22*m2.m20 + m1.m23*m2.m30;
	m3.m21 = m1.m20*m2.m01 + m1.m21*m2.m11 + m1.m22*m2.m21 + m1.m23*m2.m31;
	m3.m22 = m1.m20*m2.m02 + m1.m21*m2.m12 + m1.m22*m2.m22 + m1.m23*m2.m32;
	m3.m23 = m1.m20*m2.m03 + m1.m21*m2.m13 + m1.m22*m2.m23 + m1.m23*m2.m33;

	m3.m30 = m1.m30*m2.m00 + m1.m31*m2.m10 + m1.m32*m2.m20 + m1.m33*m2.m30;
	m3.m31 = m1.m30*m2.m01 + m1.m31*m2.m11 + m1.m32*m2.m21 + m1.m33*m2.m31;
	m3.m32 = m1.m30*m2.m02 + m1.m31*m2.m12 + m1.m32*m2.m22 + m1.m33*m2.m32;
	m3.m33 = m1.m30*m2.m03 + m1.m31*m2.m13 + m1.m32*m2.m23 + m1.m33*m2.m33;

    m00 = m3.m00; m01 = m3.m01; m02 = m3.m02; m03 = m3.m03;
    m10 = m3.m10; m11 = m3.m11; m12 = m3.m12; m13 = m3.m13;
    m20 = m3.m20; m21 = m3.m21; m22 = m3.m22; m23 = m3.m23;
    m30 = m3.m30; m31 = m3.m31; m32 = m3.m32; m33 = m3.m33;

	return *this;
}

// Ｘ軸周りの回転行列作成
template<class T>
Matrix4<T>& Matrix4<T>::rotX(T angle)
{
    // caution: float精度
	T s = (T)sinf((float)angle);
	T c = (T)cosf((float)angle);

	m00 = 1; m01 = 0;  m02 = 0; m03 = 0;
	m10 = 0; m11 = c;  m12 = s; m13 = 0;
	m20 = 0; m21 = -s; m22 = c; m23 = 0;
	m30 = 0; m31 = 0;  m32 = 0; m33 = 1;
	return *this;
}


// Ｙ軸周りの回転行列作成
template<class T>
Matrix4<T>& Matrix4<T>::rotY(T angle)
{
    // caution: float精度
	T s = (T)(sinf((float)angle));
	T c = (T)(cosf((float)angle));

	m00 = c; m01 = 0; m02 = -s; m03 = 0;
	m10 = 0; m11 = 1; m12 = 0;  m13 = 0;
	m20 = s; m21 = 0; m22 = c;  m23 = 0;
	m30 = 0; m31 = 0; m32 = 0;  m33 = 1;

	return *this;
}


// Ｚ軸周りの回転行列作成
template<class T>
Matrix4<T>& Matrix4<T>::rotZ(T angle)
{
	T s = (T)(sinf((float)angle));
	T c = (T)(cosf((float)angle));

	m00 = c;  m01 = s; m02 = 0; m03 = 0;
	m10 = -s; m11 = c; m12 = 0; m13 = 0;
	m20 = 0;  m21 = 0; m22 = 1; m23 = 0;
	m30 = 0;  m31 = 0; m32 = 0; m33 = 1;

	return *this;
}


// スケーリング
template<class T>
Matrix4<T>& Matrix4<T>::scale(T scale)
{
	m00 = scale; m01 = 0;     m02 = 0;     m03 = 0;
	m10 = 0;     m11 = scale; m12 = 0;     m13 = 0;
	m20 = 0;     m21 = 0;     m22 = scale; m23 = 0;
	m30 = 0;     m31 = 0;     m32 = 0;     m33 = 1;

	return *this;
}


// スケーリング
template<class T>
Matrix4<T>& Matrix4<T>::scale(const Tuple3<T>& scale)
{
	m00 = scale.x; m01 = 0;       m02 = 0;       m03 = 0;
	m10 = 0;       m11 = scale.y; m12 = 0;       m13 = 0;
	m20 = 0;       m21 = 0;       m22 = scale.z; m23 = 0;
	m30 = 0;       m31 = 0;       m32 = 0;       m33 = 1;

	return *this;
}


// 移動
template<class T>
Matrix4<T>& Matrix4<T>::translation(const Vector3<T>& trans)
{
	m00 = 1;       m01 = 0;       m02 = 0;       m03 = 0;
	m10 = 0;       m11 = 1;       m12 = 0;       m13 = 0;
	m20 = 0;       m21 = 0;       m22 = 1;       m23 = 0;
	m30 = trans.x; m31 = trans.y; m32 = trans.z; m33 = 1;

	return *this;
}


// ベクトルを行列で変換。４列目の平行移動成分は使用されない。
// CAUTION: alias unsafe.
template<class T>
void Matrix4<T>::transform(const Vector3<T>& vec, Tuple3<T>* vecOut) const
{
	vecOut->x = vec.x * m00 + vec.y * m10 + vec.z * m20;
	vecOut->y = vec.x * m01 + vec.y * m11 + vec.z * m21;
	vecOut->z = vec.x * m02 + vec.y * m12 + vec.z * m22;
}


// 点を行列で変換。４列目の平行移動成分が使用（プラス）される。
// CAUTION: alias unsafe.
template<class T>
void Matrix4<T>::transform(const Point3<T>& point, Point3<T>* pointOut) const
{
	pointOut->x = point.x * m00 + point.y * m10 + point.z * m20 + m30;
	pointOut->y = point.x * m01 + point.y * m11 + point.z * m21 + m31;
	pointOut->z = point.x * m02 + point.y * m12 + point.z * m22 + m32;
}

// ３次元の点を４次元の同次座標に拡張して行列で変換
// CAUTION: alias unsafe.
template<class T>
void Matrix4<T>::transform(const Point3<T>& point, Point4<T>* pointOut) const
{
	pointOut->x = point.x * m00 + point.y * m10 + point.z * m20 + m30;
	pointOut->y = point.x * m01 + point.y * m11 + point.z * m21 + m31;
	pointOut->z = point.x * m02 + point.y * m12 + point.z * m22 + m32;
	pointOut->w = point.x * m03 + point.y * m13 + point.z * m23 + m33;
}

// ４次元の点を行列で変換
// CAUTION: alias unsafe.
template<class T>
void Matrix4<T>::transform(const Point4<T>& point, Point4<T>* pointOut) const
{
	pointOut->x = point.x * m00 + point.y * m10 + point.z * m20 + point.w * m30;
	pointOut->y = point.x * m01 + point.y * m11 + point.z * m21 + point.w * m31;
	pointOut->z = point.x * m02 + point.y * m12 + point.z * m22 + point.w * m32;
	pointOut->w = point.x * m03 + point.y * m13 + point.z * m23 + point.w * m33;
}

// クォータニオンを行列に変換
template<class T>
void Matrix4<T>::set(const Quat4<T>& q)
{
	T w = q.w, x = q.x, y = q.y, z = q.z;
	T n = x*x + y*y + z*z + w*w;
	T s = 2.0f / n;		// 正規化の為にノルムで割る

	T xs = x*s,  ys = y*s,  zs = z*s;
	T wx = w*xs, wy = w*ys, wz = w*zs;
	T xx = x*xs, xy = x*ys, xz = x*zs;
	T yy = y*ys, yz = y*zs, zz = z*zs;

	m00 = 1.0f - (yy + zz);	m01 = xy + wz;          m02 = xz - wy;          m03 = 0;
	m10 = xy - wz;          m11 = 1.0f - (xx + zz); m12 = yz + wx;          m13 = 0;
	m20 = xz + wy;          m21 = yz - wx;          m22 = 1.0f - (xx + yy); m23 = 0;
	m30 = 0;                m31 = 0;                m32 = 0;                m33 = 1;
}

// 逆行列を求める（正規直交行列限定）
template<class T>
void Matrix4<T>::invertQuick()
{
	// 転置する
	Matrix4<T> mt;
	for(int i=0; i<4 ; i++) {
		for(int j=0; j<4 ; j++) {
			mt.m[i][j] = m[j][i];
		}
	}
	
	// 移動成分を転置行列の回転成分で変換して、-1を掛けると
	mt.m30 = -(m30 * mt.m00 + m31 * mt.m10 + m32 * mt.m20);
	mt.m31 = -(m30 * mt.m01 + m31 * mt.m11 + m32 * mt.m21);
	mt.m32 = -(m30 * mt.m02 + m31 * mt.m12 + m32 * mt.m22);
	mt.m33 = (T)1;
	mt.m03 = mt.m13 = mt.m[2][3] = (T)0;
	*this  = mt;
}

// 逆行列を求める（汎用低速版）
template<class T>
void Matrix4<T>::invert() {
	T s = determinant();
	if (s == 0.0)
	    return;
	s = 1/s;
	// alias-safe way.
	// less *,+,- calculation than expanded expression.
	set(
	    m11*(m22*m33 - m23*m32) + m12*(m23*m31 - m21*m33) + m13*(m21*m32 - m22*m31),
	    m21*(m02*m33 - m03*m32) + m22*(m03*m31 - m01*m33) + m23*(m01*m32 - m02*m31),
	    m31*(m02*m13 - m03*m12) + m32*(m03*m11 - m01*m13) + m33*(m01*m12 - m02*m11),
	    m01*(m13*m22 - m12*m23) + m02*(m11*m23 - m13*m21) + m03*(m12*m21 - m11*m22),

	    m12*(m20*m33 - m23*m30) + m13*(m22*m30 - m20*m32) + m10*(m23*m32 - m22*m33),
	    m22*(m00*m33 - m03*m30) + m23*(m02*m30 - m00*m32) + m20*(m03*m32 - m02*m33),
	    m32*(m00*m13 - m03*m10) + m33*(m02*m10 - m00*m12) + m30*(m03*m12 - m02*m13),
	    m02*(m13*m20 - m10*m23) + m03*(m10*m22 - m12*m20) + m00*(m12*m23 - m13*m22),

	    m13*(m20*m31 - m21*m30) + m10*(m21*m33 - m23*m31) + m11*(m23*m30 - m20*m33),
	    m23*(m00*m31 - m01*m30) + m20*(m01*m33 - m03*m31) + m21*(m03*m30 - m00*m33),
	    m33*(m00*m11 - m01*m10) + m30*(m01*m13 - m03*m11) + m31*(m03*m10 - m00*m13),
	    m03*(m11*m20 - m10*m21) + m00*(m13*m21 - m11*m23) + m01*(m10*m23 - m13*m20),

	    m10*(m22*m31 - m21*m32) + m11*(m20*m32 - m22*m30) + m12*(m21*m30 - m20*m31),
	    m20*(m02*m31 - m01*m32) + m21*(m00*m32 - m02*m30) + m22*(m01*m30 - m00*m31),
	    m30*(m02*m11 - m01*m12) + m31*(m00*m12 - m02*m10) + m32*(m01*m10 - m00*m11),
	    m00*(m11*m22 - m12*m21) + m01*(m12*m20 - m10*m22) + m02*(m10*m21 - m11*m20)
	    );

	mul(s);
}

// 行列式を求める
template<class T>
T Matrix4<T>::determinant() const {
	// less *,+,- calculation than expanded expression.
	return
	    (m00*m11 - m01*m10)*(m22*m33 - m23*m32)
	   -(m00*m12 - m02*m10)*(m21*m33 - m23*m31)
	   +(m00*m13 - m03*m10)*(m21*m32 - m22*m31)
	   +(m01*m12 - m02*m11)*(m20*m33 - m23*m30)
	   -(m01*m13 - m03*m11)*(m20*m32 - m22*m30)
	   +(m02*m13 - m03*m12)*(m20*m31 - m21*m30);
}



// 明示的インスタンス生成
template class Matrix4<float>;
template class Matrix4<double>;


