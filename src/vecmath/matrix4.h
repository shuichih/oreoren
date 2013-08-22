#ifndef SEVERE3D_MATRIX4_H
#define SEVERE3D_MATRIX4_H

#include "tuple3.h"
#include "vector3.h"
#include "point3.h"
#include "point4.h"
#include "quat4.h"



/**
 * ４×４行列クラスです。
 */
template <class T>
class Matrix4
{
public:

	/** 行列の要素 */
    union {
        struct {
            T        m00, m01, m02, m03;
            T        m10, m11, m12, m13;
            T        m20, m21, m22, m23;
            T        m30, m31, m32, m33;
        };
        T m[4][4];
    };

	/**
	 * デフォルトコンストラクタです。
	 */
	Matrix4();
	
	/**
	 * 全ての行列要素を引数にとるコンストラクタです。
	 */
    Matrix4(T m00, T m01, T m02, T m03, 
            T m10, T m11, T m12, T m13,
            T m20, T m21, T m22, T m23,
            T m30, T m31, T m32, T m33);

	/**
	 * 全要素を設定します。
	 */
	void set(T n00, T n01, T n02, T n03, 
             T n10, T n11, T n12, T n13,
             T n20, T n21, T n22, T n23,
             T n30, T n31, T n32, T n33);

	/**
	 * 単位行列を作成します。
	 */
	Matrix4<T>& identity();

	/**
	 * 行列の定数倍を計算し、この行列に設定します。
	 */
	void mul(T scalar);

    /**
     * この行列とm1を乗じた結果を、この行列に設定します。
	 *
	 * @param m1 乗じる行列
     */
    void mul(const Matrix4<T>& m1) {
        mul(*this, m1);
    }

	/**
     * この行列とm1を乗じた結果を、この行列に設定します。
	 *
	 * @param m1 乗じる行列
	 */
	Matrix4<T>& operator*=(const Matrix4<T>& m1) {
        mul(m1);
		return *this;
    }

	/**
     * この行列とm1を乗じた結果を、新たな行列として返します。
	 *
	 * @param m1 乗じる行列
	 */
    Matrix4<T> operator*(const Matrix4<T>& m1) const {
        return (Matrix4(*this)).operator*=(m1);
    }

	/**
     * m1とm2を乗じた結果を、この行列に設定します。
     * 
     * @param m1 １つめの行列
     * @param m2 ２つめの行列
     */
    Matrix4& mul(const Matrix4<T>& m1, const Matrix4<T>& m2);

	/**
	 * Ｘ軸周りの回転行列を作成します。
	 */
	Matrix4& rotX(T angle);

	/**
	 * Ｙ軸周りの回転行列を作成します。
	 */
	Matrix4& rotY(T angle);

	/**
	 * Ｚ軸周りの回転行列を作成します。
	 */
	Matrix4& rotZ(T angle);

	/**
	 * スケーリング行列を作成します。
	 */
	Matrix4& scale(T scale);

	/**
	 * スケーリング行列を作成します。
	 */
	Matrix4& scale(const Tuple3<T>& scale);

	/**
	 * 平行移動行列を作成します。
	 */
	Matrix4& translation(const Vector3<T>& trans);

	/**
     * vecをMatrix4を使用して変換します。結果はvecOutに格納されます。
	 * ベクトルの変換では、定義上、行列の移動成分は無視されます。
     * 
     * @param vec 変換されるベクトル
     * @param vecOut 変換結果が格納されるベクトル
     */
    void transform(const Vector3<T>& vec, Tuple3<T>* vecOut) const;

	/**
     * pointをMatrix4を使用して変換します。結果はpointOutに格納されます。
     * 
     * @param point 変換される点
     * @param pointOut 変換結果が格納される点
     */
    void transform(const Point3<T>& point, Point3<T>* pointOut) const;

	/**
	 * ３次元の点を４次元の同次点[x, y, z, 1]に拡張し、この行列で変換します。
	 * 結果はpointOutに格納されます。
	 */
    void transform(const Point3<T>& point, Point4<T>* pointOut) const;

	/**
	 * ４次元の点を行列で変換
	 */
	void transform(const Point4<T>& point, Point4<T>* pointOut) const;

	/**
	 * クォータニオンを設定
	 */
	void set(const Quat4<T>& q);

	/**
	 * 行列式を計算します。
	 */
	T determinant() const;

	/**
	 * 逆行列にします。（正規直交基底行列限定）
	 */
	void invertQuick();

	/**
	 * 逆行列にします。（汎用低速版）
	 */
	void invert();

	// NOTE:コピーコンストラクタおよび = 演算子はコンパイラによって生成されます。
};

typedef Matrix4<float> Matrix4f;
typedef Matrix4<double> Matrix4d;




#endif // SEVERE3D_MATRIX4_H

