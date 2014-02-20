#ifndef SEVERE3D_POINT4_H
#define SEVERE3D_POINT4_H

#include <cmath>


/**
 * ４次元の点を表現するクラスです。
 */
template <class T>
class Point4
{
public:

	T x;
	T y;
	T z;
	T w;

    Point4(): x(0), y(0), z(0), w(0) {}

    Point4(T _x, T _y, T _z, T _w): x(_x), y(_y), z(_z), w(_w) {}

};


typedef Point4<float> Point4f;
typedef Point4<double> Point4d;


#endif // SEVERE3D_POINT4_H
