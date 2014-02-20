#ifndef SEVERE3D_POINT3_H
#define SEVERE3D_POINT3_H

#include <cmath>
#include "tuple3.h"


template <typename T> class Vector3;

/**
 * ３次元の点を表現するクラスです。
 */
template <class T>
class Point3 : public Tuple3<T>
{
public:

    /**
     * デフォルトコンストラクタです。
     */
    Point3(): Tuple3<T>() { }

	/**
	 * コンストラクタです。
	 */
    Point3(T x, T y, T z): Tuple3<T>(x, y, z) { }
    
	Point3(const Tuple3<T>& t): Tuple3<T>(t) { }

	Point3& operator=(const Tuple3<T>& t) {
		this->x = t.x;
		this->y = t.y;
		this->z = t.z;
		return *this;
	}

    /**
      * Computes the square of the distance between this point and point p1.
      * @param  p1 the other point
      * @return the square of distance between these two points as a float
      */
    T distanceSquared(const Point3& p1) const {
        T dx = this->x - p1.x;
        T dy = this->y - p1.y;
        T dz = this->z - p1.z;
        return dx*dx + dy*dy + dz*dz;
    }

    /**
      * Returns the distance between this point and point p1.
      * @param p1 the other point
      * @return the distance between these two points as a float
      */
    T distance(const Point3& p1) const {
        return sqrt(distanceSquared(p1));
    }

    Vector3<T> operator-(const Tuple3<T>& t1) const {
		return Vector3<T>(this->x - t1.x, this->y - t1.y, this->z - t1.z);
    }

	// NOTE:コピーコンストラクタおよび = 演算子はコンパイラによって生成されます。

};


typedef Point3<float> Point3f;
typedef Point3<double> Point3d;



#endif // SEVERE3D_POINT3_H
