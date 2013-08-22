#ifndef SEVERE3D_POINT2_H
#define SEVERE3D_POINT2_H

#include <cmath>
#include "vecmath_def.h"
#include "severe3d/base.h"
#include "severe3d/vecmath/tuple2.h"

/**
 * ３次元の点を表現するクラスです。
 */
template <class T>
class Point2 : public Tuple2<T>
{
public:

    /**
     * デフォルトコンストラクタです。
     */
    Point2(): Tuple2<T>() { }

	/**
	 * コンストラクタです。
	 */
    Point2(T x, T y): Tuple2<T>(x, y) { }
    
	Point2(const Tuple2<T>& t): Tuple2<T>(t) { }

	Point2& operator=(const Tuple2& t) {
		this->x = t.x;
		this->y = t.y;
		return *this;
	}

    /**
      * Computes the square of the distance between this point and point p1.
      * @param  p1 the other point
      * @return the square of distance between these two points as a float
      */
    T distanceSquared(const Point2& p1) const {
        T dx = x - p1.x;
        T dy = y - p1.y;
        return dx*dx + dy*dy;
    }

    /**
      * Returns the distance between this point and point p1.
      * @param p1 the other point
      * @return the distance between these two points as a float
      */
    T distance(const Point2& p1) const {
        return sqrt(distanceSquared(p1));
    }

	// NOTE:コピーコンストラクタおよび = 演算子はコンパイラによって生成されます。

};


typedef Point2<float> Point2f;
typedef Point2<double> Point2d;



#endif // SEVERE3D_POINT2_H
